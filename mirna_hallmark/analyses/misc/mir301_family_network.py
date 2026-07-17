"""miR-301/130/454 family: regulator-resolution network, family pressures, and
orphan validation prioritization.

This module sits on top of :mod:`mirna_hallmark.analyses.mir301.mir301_family_depth` (which it
reads for the per-arm partial correlations) and answers four linked questions
that the depth module exposed but did not synthesize:

1. **Regulator resolution** — for each focus gene, *among all of its miRNA
   regulators*, how much realized pressure does each family arm carry
   (``global_abs_share``, no-z attribution), how does that rank against its
   abundance rank, and is the share *reflected* in observed CN+proliferation
   adjusted anti-correlation? Each (gene, arm) edge is classified into an
   interpretable quadrant (driver / structural-only / minor-but-real / weak).

2. **Specificity (global vs specific)** — what fraction of an arm's *total*
   realized pressure mass (across every Hallmark target it regulates) lands on
   this one gene. Low = the arm is a promiscuous family hub spreading thin; high
   = the gene is a concentrated target. This is the "is it a global or specific
   target of this miRNA" axis.

3. **Family-as-unit pressure at PAM50 resolution** — collapse the seed family
   into one regulatory unit and report its realized pressure share per gene
   *within each PAM50 subtype*, under both the canonical z-spine
   (``softmax_z_logrpm``) and the new no-z attribution formula
   (``softmax_logrpm``) so the formula change is auditable.

4. **Orphan validation prioritization** — a transparent composite that ranks
   (arm, gene) edges for wet-lab functional validation, combining sequence prior
   (TargetScan), in-tumor coupling (FDR + proliferation surviving), specificity,
   PAM50 concentration, and literature status, with a suggested PAM50-matched
   cell-line context.

Outputs under ``mirna_hallmark/output/tissue_reference/mir301_family_network/``.

Run (depth module must have been run first):
  .venv/bin/python3 -m mirna_hallmark.analyses.mir301.mir301_family_depth
  .venv/bin/python3 -m mirna_hallmark.mir301_family_network
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.mir301.mir301_focus_genes import (
    PRIOR_ORPHAN_HITS,
    PUBLISHED_BRCA,
    evidence_class as _evidence_class,
    he_pairs as _he_pairs,
)
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.analyses.mir301.mir301_family_depth import (
    FAMILY_ARMS,
    OUT_DIR as FAM_DEPTH_DIR,
    _focus_genes,
)
from mirna_hallmark.pressure_build import (
    compute_gene_pressure_contributions,
    load_mirtar_edges,
)
from mirna_hallmark.robustness_checks import _load_targetscan_weights, _pam50_scope_iter
from mirna_hallmark.eval.targetscan_orphan_coupling import _lit_spot_check


OUT_DIR = C.TISSUE_REFERENCE_DIR / "mir301_family_network"

# No-z attribution mode (see config.PRESSURE_ATTRIBUTION_EXPR_MODE). Under this
# mode every contribution is non-negative, so abs == signed share and "who carries
# the pressure" is unambiguous and not gated by the z-term.
ATTR_MODE = C.PRESSURE_ATTRIBUTION_EXPR_MODE
SPINE_MODE = C.PRESSURE_EXPR_MODE  # softmax_z_logrpm, for the formula contrast
# Structural (edge_w·sm, no logrpm) mode for the de-abundance-double-counted "preferred
# regulator" resolution axis. The quadrant's "carries pressure mass" test keys off this,
# not the abundance-driven realized share; abundance stays a separate reported axis.
STRUCT_MODE = C.PRESSURE_STRUCT_EXPR_MODE

# Thresholds for interpretable classification (documented, not magic).
HIGH_SHARE_RANK_PCT = 0.25   # top quartile of a gene's regulators by realized share
SPECIFIC_MIN = 0.05          # gene captures >=5% of the arm's total target mass
PROMISCUOUS_MAX = 0.01       # gene captures <1% -> arm is spread thin over many targets
PAM50_SCOPES = ("LumA", "LumB", "Her2", "Basal")

# PAM50-matched BRCA cell lines for suggested validation context (field-standard
# models; not exhaustive). Used only as a recommendation annotation.
PAM50_CELL_LINES: Dict[str, str] = {
    "LumA": "MCF-7; T-47D",
    "LumB": "BT-474; ZR-75-1",
    "Her2": "SK-BR-3; HCC1954; AU565",
    "Basal": "MDA-MB-231; MDA-MB-468; HCC1937; BT-549",
    "cohort": "MCF-7 (Lum) + MDA-MB-231 (Basal) panel",
}


# --------------------------------------------------------------------------- #
# Inputs from the depth module
# --------------------------------------------------------------------------- #
def _read_depth_partials(depth_dir: Path) -> pd.DataFrame:
    path = depth_dir / "family_all_partials.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}. Run `python -m mirna_hallmark.analyses.mir301.mir301_family_depth` first."
        )
    return pd.read_csv(path, sep="\t")


def _coupling_by_gene_arm(partials: pd.DataFrame) -> pd.DataFrame:
    """Cohort + best-PAM50 coupling per (gene, arm) from single-arm partials."""
    sa = partials.loc[partials["predictor_type"].eq("single_arm")].copy()
    rho_cn = "rho_CPE_HRD_CN" if "rho_CPE_HRD_CN" in sa.columns else "rho_CPE_HRD"
    rho_pro = "rho_e2f_g2m_CN" if "rho_e2f_g2m_CN" in sa.columns else rho_cn
    rows: List[dict] = []
    for (gene, arm), sub in sa.groupby(["gene", "predictor"]):
        coh = sub.loc[sub["scope"].eq("cohort")]
        coh_rho_cn = float(coh[rho_cn].iloc[0]) if not coh.empty else np.nan
        coh_rho_pro = float(coh[rho_pro].iloc[0]) if not coh.empty and rho_pro in coh else np.nan
        coh_negfdr = bool(coh["neg_sig_cn_fdr"].iloc[0]) if not coh.empty and "neg_sig_cn_fdr" in coh else False
        coh_negpro = bool(coh["neg_sig_prolif_cn_fdr"].iloc[0]) if not coh.empty and "neg_sig_prolif_cn_fdr" in coh else False
        # best PAM50 subtype = most negative CN-adjusted rho among FDR-sig subtype scopes
        sub_pam = sub.loc[sub["scope"].isin(PAM50_SCOPES)].copy()
        best_subtype, best_rho, best_negfdr = "", np.nan, False
        if not sub_pam.empty:
            sig_pam = sub_pam.loc[sub_pam.get("neg_sig_cn_fdr", False)]
            pick = sig_pam if not sig_pam.empty else sub_pam
            j = pick[rho_cn].idxmin()
            best_subtype = str(pick.loc[j, "scope"])
            best_rho = float(pick.loc[j, rho_cn])
            best_negfdr = bool(pick.loc[j, "neg_sig_cn_fdr"]) if "neg_sig_cn_fdr" in pick else False
        rows.append(
            {
                "gene": gene,
                "miRNA": arm,
                "cohort_rho_cn": coh_rho_cn,
                "cohort_rho_prolif_cn": coh_rho_pro,
                "cohort_neg_sig_cn_fdr": coh_negfdr,
                "cohort_neg_sig_prolif_cn_fdr": coh_negpro,
                "best_subtype": best_subtype,
                "best_subtype_rho_cn": best_rho,
                "best_subtype_neg_sig_cn_fdr": best_negfdr,
                "subtype_concentration": (best_rho - coh_rho_cn)
                if pd.notna(best_rho) and pd.notna(coh_rho_cn) else np.nan,
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Specificity: global vs specific target of each family arm
# --------------------------------------------------------------------------- #
def _ts_specificity_table(ts_universe: pd.DataFrame, focus: Sequence[str]) -> pd.DataFrame:
    """TargetScan-based specificity for orphans (no miRTar/realized-pressure edge).

    ``ts_specificity`` = this gene's |context++| weight / the arm's *total* TS
    weight summed over every Hallmark gene it is predicted to target. High = a
    sequence-concentrated target; low = the gene is one of many TS sites.
    """
    if ts_universe.empty:
        return pd.DataFrame()
    fam = ts_universe.loc[ts_universe["miRNA"].isin(FAMILY_ARMS)].copy()
    arm_total = fam.groupby("miRNA")["ts_weight"].sum()
    arm_degree_ts = fam.groupby("miRNA")["gene"].nunique()
    f = fam.loc[fam["gene"].isin(set(focus))].copy()
    f["ts_arm_total"] = f["miRNA"].map(arm_total)
    f["ts_arm_degree"] = f["miRNA"].map(arm_degree_ts)
    f["ts_specificity"] = f["ts_weight"] / f["ts_arm_total"].replace(0, np.nan)
    return f[["gene", "miRNA", "ts_weight", "ts_arm_total", "ts_arm_degree", "ts_specificity"]]


def _family_specificity(mirna: pd.DataFrame) -> pd.DataFrame:
    """Per (arm, gene) fraction of the arm's *total* realized pressure mass.

    Computed over every Hallmark gene the family arms regulate (miRTar edges),
    under the no-z attribution mode so the mass is a clean magnitude.
    """
    hs = HallmarkSets.load()
    universe = sorted(hs.universe)
    edges = load_mirtar_edges(universe, resolve_arms=True)
    fam_edges = edges.loc[edges["miRNA"].isin(FAMILY_ARMS)]
    fam_target_genes = sorted(set(fam_edges["gene"]))
    if not fam_target_genes:
        return pd.DataFrame()
    contrib = compute_gene_pressure_contributions(
        fam_target_genes, expr_mode=ATTR_MODE, mirna=mirna,  # type: ignore[arg-type]
    )
    fam = contrib.loc[contrib["miRNA"].isin(FAMILY_ARMS)].copy()
    if fam.empty:
        return fam
    totals = fam.groupby("miRNA")["mean_abs_contribution"].sum()
    degree = fam.groupby("miRNA")["gene"].nunique()
    fam["arm_total_mass"] = fam["miRNA"].map(totals)
    fam["arm_degree"] = fam["miRNA"].map(degree)
    fam["specificity"] = fam["mean_abs_contribution"] / fam["arm_total_mass"].replace(0, np.nan)
    return fam[["gene", "miRNA", "mean_abs_contribution", "arm_total_mass",
                "arm_degree", "specificity"]]


# --------------------------------------------------------------------------- #
# Regulator-resolution edge table
# --------------------------------------------------------------------------- #
def _interpret_edge(row: pd.Series) -> str:
    coupled = bool(row.get("cohort_neg_sig_cn_fdr")) or bool(row.get("cohort_neg_sig_prolif_cn_fdr"))
    # "carries pressure mass" keys off the STRUCTURAL rank (edge_w·sm, abundance spent
    # once) so a merely-abundant arm is not auto-labelled the gene's structural driver;
    # falls back to the realized share rank only if the structural rank is unavailable.
    rank_pct = row.get("struct_share_rank_pct")
    if pd.isna(rank_pct):
        rank_pct = row.get("share_rank_pct")
    high_share = pd.notna(rank_pct) and rank_pct <= HIGH_SHARE_RANK_PCT
    spec = row.get("specificity", np.nan)
    if pd.notna(spec) and spec >= SPECIFIC_MIN:
        spec_tag = "specific"
    elif pd.notna(spec) and spec <= PROMISCUOUS_MAX:
        spec_tag = "promiscuous"
    else:
        spec_tag = "moderate"
    if coupled and high_share:
        base = "driver_candidate"
    elif coupled and not high_share:
        base = "real_but_minor"          # observed repression but crowded out on pressure
    elif (not coupled) and high_share:
        base = "structural_only"         # carries pressure mass but no observed coupling
    else:
        base = "weak_unsupported"
    return f"{base}|{spec_tag}"


def _regulator_resolution(
    focus: Sequence[str], mirna: pd.DataFrame, coupling: pd.DataFrame, spec: pd.DataFrame,
    ts_long: pd.DataFrame, he_pairs,
) -> pd.DataFrame:
    """For each focus gene, family-arm position among ALL its regulators."""
    contrib = compute_gene_pressure_contributions(
        list(focus), expr_mode=ATTR_MODE, mirna=mirna,  # type: ignore[arg-type]
    )
    if contrib.empty:
        return contrib
    # rank every regulator within each gene by realized share and by abundance
    contrib = contrib.sort_values(["gene", "global_abs_share"], ascending=[True, False])
    contrib["share_rank"] = contrib.groupby("gene")["global_abs_share"].rank(
        ascending=False, method="min"
    )
    contrib["n_regulators"] = contrib.groupby("gene")["miRNA"].transform("count")
    contrib["share_rank_pct"] = contrib["share_rank"] / contrib["n_regulators"]
    contrib["abund_rank"] = contrib.groupby("gene")["median_log2rpm"].rank(
        ascending=False, method="min"
    )
    contrib["abund_rank_pct"] = contrib["abund_rank"] / contrib["n_regulators"]

    # STRUCTURAL share rank (edge_w·sm, no logrpm double-count) — the de-abundance axis
    # the quadrant classifies on. Computed on the same edges/mirna under the softmax mode.
    struct = compute_gene_pressure_contributions(
        list(focus), expr_mode=STRUCT_MODE, mirna=mirna,  # type: ignore[arg-type]
    )
    if not struct.empty:
        struct["struct_share_rank"] = struct.groupby("gene")["global_abs_share"].rank(
            ascending=False, method="min"
        )
        struct["struct_n_reg"] = struct.groupby("gene")["miRNA"].transform("count")
        struct["struct_share_rank_pct"] = struct["struct_share_rank"] / struct["struct_n_reg"]
        struct = struct.rename(columns={"global_abs_share": "struct_global_abs_share"})
        contrib = contrib.merge(
            struct[["gene", "miRNA", "struct_global_abs_share",
                    "struct_share_rank", "struct_share_rank_pct"]],
            on=["gene", "miRNA"], how="left",
        )

    fam = contrib.loc[contrib["miRNA"].isin(FAMILY_ARMS)].copy()
    ts_lookup = {(r["miRNA"], r["gene"]): r["ts_weight"] for _, r in ts_long.iterrows()} \
        if not ts_long.empty else {}
    fam["ts_weight"] = fam.apply(lambda r: ts_lookup.get((r["miRNA"], r["gene"]), np.nan), axis=1)
    fam = fam.merge(coupling, on=["gene", "miRNA"], how="left")
    fam = fam.merge(spec[["gene", "miRNA", "arm_degree", "specificity"]],
                    on=["gene", "miRNA"], how="left")
    fam["evidence_class"] = fam.apply(
        lambda r: _evidence_class(r["gene"], r["miRNA"], he_pairs), axis=1
    )
    fam["literature_note"] = fam.apply(
        lambda r: _lit_spot_check(r["gene"], r["miRNA"]), axis=1
    )
    fam["edge_class"] = fam.apply(_interpret_edge, axis=1)
    keep = [
        "gene", "miRNA", "evidence_class", "edge_class",
        "global_abs_share", "share_rank", "n_regulators", "share_rank_pct",
        "struct_global_abs_share", "struct_share_rank", "struct_share_rank_pct",
        "median_log2rpm", "abund_rank", "abund_rank_pct",
        "specificity", "arm_degree", "ts_weight",
        "cohort_rho_cn", "cohort_rho_prolif_cn",
        "cohort_neg_sig_cn_fdr", "cohort_neg_sig_prolif_cn_fdr",
        "best_subtype", "best_subtype_rho_cn", "best_subtype_neg_sig_cn_fdr",
        "subtype_concentration", "literature_note",
    ]
    keep = [c for c in keep if c in fam.columns]
    return fam[keep].sort_values(["gene", "global_abs_share"], ascending=[True, False])


# --------------------------------------------------------------------------- #
# Family-as-unit pressure at PAM50 resolution (+ formula contrast)
# --------------------------------------------------------------------------- #
def _family_pressure_by_subtype(focus: Sequence[str], mirna: pd.DataFrame) -> pd.DataFrame:
    clinical = D.load_clinical_strata()
    rows: List[dict] = []
    for mode in (SPINE_MODE, ATTR_MODE):
        for scope, samples in _pam50_scope_iter(clinical):
            cols = mirna.columns if samples is None else mirna.columns.intersection(list(samples))
            if len(cols) < 20:
                continue
            sub_mirna = mirna[cols]
            contrib = compute_gene_pressure_contributions(
                list(focus), expr_mode=mode, mirna=sub_mirna,  # type: ignore[arg-type]
            )
            if contrib.empty:
                continue
            fam = contrib.loc[contrib["miRNA"].isin(FAMILY_ARMS)]
            for gene, g in fam.groupby("gene"):
                rows.append(
                    {
                        "gene": gene,
                        "expr_mode": mode,
                        "scope": scope,
                        "n_samples": int(len(cols)),
                        "family_abs_share": float(g["global_abs_share"].sum(skipna=True)),
                        "family_signed_share": float(g["global_signed_share"].sum(skipna=True)),
                        "family_pressure_mass": float(g["mean_abs_contribution"].sum(skipna=True)),
                        "n_family_arms_regulating": int(g["miRNA"].nunique()),
                        "top_family_arm": str(g.sort_values("global_abs_share", ascending=False)["miRNA"].iloc[0]),
                    }
                )
    return pd.DataFrame(rows).sort_values(["gene", "expr_mode", "scope"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Orphan validation prioritization
# --------------------------------------------------------------------------- #
def _zscore(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    sd = s.std()
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(0.0, index=s.index)
    return (s - s.mean()) / sd


def _validation_priority(edge_master: pd.DataFrame) -> pd.DataFrame:
    """Transparent composite ranking of (arm, gene) edges for functional validation.

    Operates on the *partials-derived* edge master (every family arm x focus gene,
    TS-orphans included) rather than the realized-pressure table — orphans have no
    miRTar edge so they carry zero realized pressure and would otherwise be dropped.

    Components (z-scored across candidate edges, equal-weight mean):
      - coupling_strength    : -rho (proliferation+CN adjusted if it survives, else CN)
      - ts_weight            : TargetScan sequence prior (site model)
      - ts_specificity       : concentrated TS target -> cleaner single-gene readout
      - subtype_concentration: more negative in a subtype than cohort -> context signal
    novelty_weight upweights orphan/low edges (1.0) over already-published ones (0.6).
    Only edges with any FDR-significant negative coupling (cohort or best subtype) qualify.
    """
    if edge_master.empty:
        return edge_master
    df = edge_master.copy()
    qualifies = (
        df["cohort_neg_sig_cn_fdr"].fillna(False)
        | df["cohort_neg_sig_prolif_cn_fdr"].fillna(False)
        | df["best_subtype_neg_sig_cn_fdr"].fillna(False)
    )
    df = df.loc[qualifies].copy()
    if df.empty:
        return df
    # coupling strength: prefer proliferation-surviving rho, else CN rho, else best subtype
    def _coupling(r):
        if bool(r.get("cohort_neg_sig_prolif_cn_fdr")) and pd.notna(r.get("cohort_rho_prolif_cn")):
            return -float(r["cohort_rho_prolif_cn"])
        if bool(r.get("cohort_neg_sig_cn_fdr")) and pd.notna(r.get("cohort_rho_cn")):
            return -float(r["cohort_rho_cn"])
        if pd.notna(r.get("best_subtype_rho_cn")):
            return -float(r["best_subtype_rho_cn"])
        return 0.0
    df["coupling_strength"] = df.apply(_coupling, axis=1)
    df["subtype_concentration_pos"] = (-df["subtype_concentration"]).clip(lower=0).fillna(0)
    df["novelty_weight"] = np.where(
        df["evidence_class"].isin(["ts_orphan", "mirtar_low"]), 1.0, 0.6
    )
    comp = (
        _zscore(df["coupling_strength"])
        + _zscore(df["ts_weight"])
        + _zscore(df["ts_specificity"])
        + 0.5 * _zscore(df["subtype_concentration_pos"])
    ) / 3.5
    df["validation_priority"] = comp * df["novelty_weight"]
    # recommend a PAM50 context: subtype if it is FDR-sig and concentrated, else cohort
    def _context(r):
        if bool(r.get("best_subtype_neg_sig_cn_fdr")) and pd.notna(r.get("subtype_concentration")) \
                and r["subtype_concentration"] < 0:
            return r["best_subtype"]
        return "cohort"
    df["recommended_context"] = df.apply(_context, axis=1)
    df["recommended_cell_lines"] = df["recommended_context"].map(PAM50_CELL_LINES).fillna(
        PAM50_CELL_LINES["cohort"]
    )
    cols = [
        "gene", "miRNA", "evidence_class", "validation_priority",
        "coupling_strength", "ts_weight", "ts_specificity", "ts_arm_degree",
        "cohort_rho_cn", "cohort_neg_sig_cn_fdr", "cohort_neg_sig_prolif_cn_fdr",
        "best_subtype", "best_subtype_rho_cn", "subtype_concentration", "novelty_weight",
        "recommended_context", "recommended_cell_lines", "literature_note",
    ]
    cols = [c for c in cols if c in df.columns]
    return df[cols].sort_values("validation_priority", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(*, out_dir: Path = OUT_DIR, depth_dir: Path = FAM_DEPTH_DIR,
        top_orphans: int = 25) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    focus = _focus_genes(top_orphans=top_orphans)
    mirna = D.load_mirna_arms()
    he_pairs = _he_pairs()

    partials = _read_depth_partials(depth_dir)
    coupling = _coupling_by_gene_arm(partials)

    hs = HallmarkSets.load()
    ts_universe = _load_targetscan_weights(sorted(hs.universe))
    ts_spec = _ts_specificity_table(ts_universe, focus)
    ts_long = ts_spec.rename(columns={"ts_weight": "ts_weight"})[["gene", "miRNA", "ts_weight"]] \
        if not ts_spec.empty else pd.DataFrame()

    # edge master = every family arm x focus gene with coupling + TS + evidence class
    edge_master = coupling.copy()
    if not ts_spec.empty:
        edge_master = edge_master.merge(ts_spec, on=["gene", "miRNA"], how="left")
    edge_master["evidence_class"] = edge_master.apply(
        lambda r: _evidence_class(r["gene"], r["miRNA"], he_pairs), axis=1
    )
    edge_master["literature_note"] = edge_master.apply(
        lambda r: _lit_spot_check(r["gene"], r["miRNA"]), axis=1
    )

    spec = _family_specificity(mirna)
    reg = _regulator_resolution(focus, mirna, coupling, spec, ts_long, he_pairs)
    fam_press = _family_pressure_by_subtype(focus, mirna)
    priority = _validation_priority(edge_master)

    spec.to_csv(out_dir / "family_arm_specificity.tsv", sep="\t", index=False)
    edge_master.to_csv(out_dir / "family_edge_master.tsv", sep="\t", index=False)
    reg.to_csv(out_dir / "regulator_resolution_edges.tsv", sep="\t", index=False)
    fam_press.to_csv(out_dir / "family_pressure_by_subtype.tsv", sep="\t", index=False)
    priority.to_csv(out_dir / "orphan_validation_priority.tsv", sep="\t", index=False)

    # arm-level network summary (promiscuity / hub structure)
    if not spec.empty:
        arm_summary = (
            spec.groupby("miRNA")
            .agg(arm_degree=("gene", "nunique"),
                 arm_total_mass=("mean_abs_contribution", "sum"),
                 median_specificity=("specificity", "median"))
            .reset_index()
            .sort_values("arm_degree", ascending=False)
        )
        arm_summary.to_csv(out_dir / "family_arm_network_summary.tsv", sep="\t", index=False)
    else:
        arm_summary = pd.DataFrame()

    edge_class_counts = (
        reg.groupby(["edge_class"]).size().reset_index(name="n_edges")
        if not reg.empty else pd.DataFrame()
    )
    edge_class_counts.to_csv(out_dir / "edge_class_counts.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.mir301_family_network",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": str(depth_dir / "family_all_partials.tsv"),
        "attribution_expr_mode": ATTR_MODE,
        "spine_expr_mode": SPINE_MODE,
        "target_norm": C.PRESSURE_TARGET_NORM,
        "share_metric_decision": {
            "primary": "global_abs_share",
            "coherence_diag": "global_signed_share",
            "attribution_mode": ATTR_MODE,
            "rationale": "no-z makes abs==signed; share unambiguous + not z-gated",
        },
        "edge_class_rules": {
            "high_share": f"share_rank_pct <= {HIGH_SHARE_RANK_PCT}",
            "specific": f"specificity >= {SPECIFIC_MIN}",
            "promiscuous": f"specificity <= {PROMISCUOUS_MAX}",
            "coupled": "cohort neg_sig_cn_fdr or neg_sig_prolif_cn_fdr",
        },
        "pam50_cell_lines": PAM50_CELL_LINES,
        "n_focus_genes": len(focus),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[mir301_family_network] focus genes={len(focus)}; edges={len(reg)}")
    if not edge_class_counts.empty:
        print("[mir301_family_network] Edge classes:\n", edge_class_counts.to_string(index=False))
    if not arm_summary.empty:
        print("\n[mir301_family_network] Arm network summary:\n", arm_summary.to_string(index=False))
    if not priority.empty:
        print("\n[mir301_family_network] Top validation candidates:\n",
              priority.head(15).to_string(index=False))
    return {
        "specificity": spec,
        "regulator_resolution": reg,
        "family_pressure_by_subtype": fam_press,
        "validation_priority": priority,
        "arm_summary": arm_summary,
        "edge_class_counts": edge_class_counts,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--depth-dir", type=Path, default=FAM_DEPTH_DIR)
    ap.add_argument("--top-orphans", type=int, default=25)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, depth_dir=args.depth_dir, top_orphans=args.top_orphans)


if __name__ == "__main__":
    main()
