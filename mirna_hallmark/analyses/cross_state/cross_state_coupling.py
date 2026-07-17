"""Cross-state coupling trajectories with composition + proliferation adjustment.

The unique angle of this platform: a pressure model that nominates *which* edges
carry regulatory mass, three tissue states (true-healthy GTEx, field-effect NAT,
tumor), and -- added here -- the confounders that actually apply to normal tissue.

Tumor *purity* (CPE/ABSOLUTE) is undefined for NAT/GTEx, but the real confounders
in normal tissue are **cell composition** (epithelial / immune / stromal fraction)
and **proliferation**, all of which are computable from the mRNA we have in every
state as metagenes. We therefore compute, per edge (miRNA arm -> target), in each
state:

  - raw Spearman rho
  - composition+proliferation-adjusted partial Spearman rho (covariates =
    epithelial, immune, stromal, proliferation metagenes; SAME covariates in all
    three states so the trajectory is comparable; tumor CN/CPE/HRD partials live
    in the dedicated tumor coupling modules)

and classify each edge into a regulatory-trajectory archetype:

  constitutive_repressor       neg & sig in tumor AND in >=1 normal
  tumor_acquired_repressor     neg & sig in tumor, no normal coupling
  repression_replaces_coexpr   positive in normal -> negative in tumor
  lost_repressor               neg in normal, gone in tumor
  normal_coexpression_only     positive in normal, neutral in tumor (lineage program)
  coexpression_all             positive everywhere
  no_consistent_signal         nothing passes

Plus a ``composition_sensitive`` flag: a raw coupling that is killed or sign-flipped
by metagene adjustment (a lineage / composition artifact rather than regulation).

This directly resolves the two motivating cases:
  - miR-9-5p -> PTEN  (NAT +0.43 raw, odd-one-out -> immune-composition test)
  - miR-200c-3p -> ETV5/ITGB4 (strong + in both normals -> epithelial-lineage test)

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_coupling                 # curated focus
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_coupling --gene PTEN
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_coupling --mirna hsa-miR-200c-3p
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_deep_dive import _state_bundles
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges
from mirna_hallmark.robustness_checks import _metagene
from mirna_hallmark.stats import bh_fdr

OUT_DIR = C.TISSUE_REFERENCE_DIR / "cross_state_coupling"
MIN_N = 25
RHO_FLOOR = 0.10
Q_ALPHA = 0.05

# Composition marker metagenes (computable from mRNA in any tissue state).
EPI_MARKERS = ["CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4",
               "CLDN7", "ESRP1", "ESRP2", "GRHL2", "OVOL2"]
IMMUNE_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD2", "CD8A", "CD4", "CD19",
                  "MS4A1", "NKG7", "GZMB", "CD68", "LYZ"]
STROMA_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRB", "ACTA2",
                  "THY1", "LUM", "DCN"]

FOCUS_ARMS = [
    "hsa-miR-210-3p", "hsa-miR-21-5p", "hsa-miR-9-5p", "hsa-miR-200c-3p",
    "hsa-miR-200b-3p", "hsa-miR-429", "hsa-miR-155-5p", "hsa-miR-301a-3p",
    "hsa-miR-301b-3p", "hsa-miR-454-3p", "hsa-miR-130a-3p", "hsa-miR-375",
    "hsa-miR-196a-5p", "hsa-miR-182-5p", "hsa-miR-205-5p",
]
TARGETS_PER_ARM = 20

_METAGENE_GENES = sorted(set(EPI_MARKERS + IMMUNE_MARKERS + STROMA_MARKERS))


def _prolif_metagene(rna: pd.DataFrame, hs: HallmarkSets) -> pd.Series:
    prolif_genes = sorted(set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
                          | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", [])))
    return _metagene(rna, prolif_genes)


def _with_state_batch(state: str, cov: pd.DataFrame) -> pd.DataFrame:
    """Append per-state sequencing-batch covariates to a cross-state covariate frame:
    tumor/NAT analyte plate (from the barcode) and GTEx SMNABTCH/SMGEBTCH. No-op when
    config.TCGA_BATCH_KIND == 'none'."""
    if not getattr(C, "TCGA_BATCH_KIND", "none") or C.TCGA_BATCH_KIND == "none":
        return cov
    from mirna_hallmark.tcga_batch import cross_state_batch_dummies

    b = cross_state_batch_dummies(state, cov.index)
    if b is None or b.empty or b.shape[1] == 0:
        return cov
    return cov.join(b, how="left")


def _state_covariates(state: str, rna: pd.DataFrame, hs: HallmarkSets, source: str) -> pd.DataFrame:
    """Composition + proliferation covariates for one state.

    source='estimate' : actual ESTIMATE ImmuneScore + StromalScore (Yoshihara 2013)
                        + proliferation metagene. The principled composition control;
                        ESTIMATE purity ~ epithelial/tumor fraction, so adjusting for
                        immune+stromal also controls the epithelial-dilution axis.
    source='metagene' : hand-curated epithelial/immune/stromal marker metagenes
                        + proliferation (the earlier robustness proxy).
    """
    prolif = _prolif_metagene(rna, hs).rename("prolif")
    if source in ("estimate", "estimate_epi"):
        from mirna_hallmark.estimate_scores import load_estimate

        est = load_estimate(state)
        cov = est[["ImmuneScore", "StromalScore"]].copy()
        cov["prolif"] = prolif.reindex(cov.index)
        if source == "estimate_epi":
            # ESTIMATE scores immune + stromal only (no epithelial axis); add an
            # explicit epithelial signature so the epithelial-lineage co-expression
            # confound (e.g. miR-200c family) is controlled too.
            cov["epi"] = _metagene(rna, EPI_MARKERS).reindex(cov.index)
        return _with_state_batch(state, cov)
    cov = pd.DataFrame({
        "prolif": prolif,
        "epi": _metagene(rna, EPI_MARKERS),
        "immune": _metagene(rna, IMMUNE_MARKERS),
        "stroma": _metagene(rna, STROMA_MARKERS),
    })
    return _with_state_batch(state, cov)


def _focus_edges(edges: pd.DataFrame, *, gene: Optional[str], mirna: Optional[str]) -> pd.DataFrame:
    if gene:
        sub = edges.loc[edges["gene"].eq(gene)]
    elif mirna:
        sub = edges.loc[edges["miRNA"].eq(mirna)]
    else:
        sub = edges.loc[edges["miRNA"].isin(FOCUS_ARMS)]
        sub = (sub.sort_values("evidence_score", ascending=False)
               .groupby("miRNA", group_keys=False).head(TARGETS_PER_ARM))
    return sub[["miRNA", "gene", "evidence_score"]].drop_duplicates().reset_index(drop=True)


def _couple_vectors(
    y: pd.Series,
    x: pd.Series,
    cov: pd.DataFrame,
) -> Tuple[float, float, float, float, int]:
    """(raw_rho, raw_p, adj_rho, adj_p, n) for paired continuous vectors."""
    from analysis.utils.common.loaders import partial_spearman

    yv = pd.to_numeric(y, errors="coerce")
    xv = pd.to_numeric(x, errors="coerce")
    shared = sorted(set(yv.dropna().index) & set(xv.dropna().index))
    if len(shared) < MIN_N or yv.loc[shared].nunique() < 2 or xv.loc[shared].nunique() < 2:
        return np.nan, np.nan, np.nan, np.nan, len(shared)
    raw_rho, raw_p = spearmanr(xv.loc[shared], yv.loc[shared])
    adj_rho, adj_p, n = partial_spearman(yv.loc[shared], xv.loc[shared], cov.reindex(shared))
    return float(raw_rho), float(raw_p), adj_rho, adj_p, int(n)


def _couple_state(
    mir: pd.DataFrame, rna: pd.DataFrame, cov: pd.DataFrame, arm: str, gene: str,
) -> Tuple[float, float, float, float, int]:
    """(raw_rho, raw_p, adj_rho, adj_p, n) for miRNA abundance vs target expression."""
    if arm not in mir.index or gene not in rna.index:
        return np.nan, np.nan, np.nan, np.nan, 0
    x = pd.to_numeric(mir.loc[arm], errors="coerce")
    y = pd.to_numeric(rna.loc[gene], errors="coerce")
    return _couple_vectors(y, x, cov)


def _couple_state_pressure(
    rna: pd.DataFrame,
    cov: pd.DataFrame,
    gene: str,
    pressure: pd.Series,
) -> Tuple[float, float, float, float, int]:
    """Partial Spearman between per-sample edge pressure and target expression."""
    if gene not in rna.index or pressure is None or pressure.empty:
        return np.nan, np.nan, np.nan, np.nan, 0
    y = pd.to_numeric(rna.loc[gene], errors="coerce")
    return _couple_vectors(y, pressure, cov)


def _sig_row(row: pd.Series, state: str, sign: int, *, rho_col: str, q_col: str) -> bool:
    rho = row.get(rho_col, np.nan)
    q = row.get(q_col, 1.0)
    if pd.isna(rho) or pd.isna(q):
        return False
    if sign < 0:
        return rho < -RHO_FLOOR and q < Q_ALPHA
    return rho > RHO_FLOOR and q < Q_ALPHA


def _edge_archetype_row(
    row: pd.Series,
    *,
    rho_prefix: str = "adj_rho",
    q_prefix: str = "adj_q",
) -> str:
    def sig(st: str, sign: int) -> bool:
        return _sig_row(row, st, sign, rho_col=f"{rho_prefix}_{st}", q_col=f"{q_prefix}_{st}")

    t_neg, t_pos = sig("tumor", -1), sig("tumor", 1)
    g_neg, g_pos = sig("gtex", -1), sig("gtex", 1)
    n_neg, n_pos = sig("nat", -1), sig("nat", 1)
    h_neg = g_neg or n_neg
    h_pos = g_pos or n_pos
    if pd.isna(row.get(f"{rho_prefix}_tumor", np.nan)):
        return "unscored"
    if t_neg and h_neg:
        return "constitutive_repressor"
    if t_neg and h_pos:
        return "repression_replaces_coexpr"
    if t_neg and not h_neg and not h_pos:
        return "tumor_acquired_repressor"
    if not t_neg and h_neg:
        return "lost_repressor"
    if not t_neg and not t_pos and h_pos:
        return "normal_coexpression_only"
    if t_pos and h_pos:
        return "coexpression_all"
    return "no_consistent_signal"


def _nat_decoupled_row(
    row: pd.Series,
    *,
    rho_prefix: str = "adj_rho",
    q_prefix: str = "adj_q",
) -> bool:
    rn = row.get(f"{rho_prefix}_nat", np.nan)
    qn = row.get(f"{q_prefix}_nat", 1.0)
    rt = row.get(f"{rho_prefix}_tumor", np.nan)
    rg = row.get(f"{rho_prefix}_gtex", np.nan)
    if pd.isna(rn) or pd.isna(rt) or pd.isna(rg):
        return False
    if abs(rn) < RHO_FLOOR or qn >= Q_ALPHA:
        return False
    return (np.sign(rn) != np.sign(rt)) and (np.sign(rn) != np.sign(rg))


def _classify(row: pd.Series) -> str:
    return _edge_archetype_row(row, rho_prefix="adj_rho", q_prefix="adj_q")


def _classify_pressure(row: pd.Series) -> str:
    return _edge_archetype_row(row, rho_prefix="adj_rho_pressure", q_prefix="adj_q_pressure")


def _states_mirna_for_pressure(impute_gtex: bool = True) -> Dict[str, pd.DataFrame]:
    """Tumor/NAT/GTEx miRNA matrices for edge pressure (GTEx collapse-repair optional)."""
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
    from mirna_hallmark.mirna_state_class import _impute_gtex_state

    states = _state_matrices()
    if impute_gtex and "gtex" in states:
        states["gtex"], _ = _impute_gtex_state(states["gtex"])
    return states


def _mirna_for_edge_pressure(
    bundles: Dict[str, Tuple[pd.DataFrame, pd.DataFrame]],
    states_mirna_imputed: Dict[str, pd.DataFrame],
) -> Dict[str, pd.DataFrame]:
    """Participant/donor-aligned miRNA matrices for per-edge pressure (GTEx imputed only)."""
    out: Dict[str, pd.DataFrame] = {}
    for s, (mir, rna) in bundles.items():
        if s == "gtex" and s in states_mirna_imputed:
            g = states_mirna_imputed[s].reindex(index=mir.index).reindex(columns=rna.columns)
            out[s] = g
        else:
            out[s] = mir
    return out


def _composition_sensitive(
    row: pd.Series,
    states: Sequence[str],
    *,
    rho_prefix: str = "adj_rho",
) -> bool:
    for s in states:
        raw = row.get(f"raw_rho_{s}", np.nan)
        adj = row.get(f"{rho_prefix}_{s}", np.nan)
        if pd.isna(raw) or pd.isna(adj):
            continue
        if abs(raw) >= 0.2 and (np.sign(raw) != np.sign(adj) or abs(adj) < RHO_FLOOR):
            return True
    return False


def run(*, out_dir: Path = OUT_DIR, gene: Optional[str] = None, mirna: Optional[str] = None,
        cov_source: str = "estimate_epi", quiet: bool = False) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    edges = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    foc = _focus_edges(edges, gene=gene, mirna=mirna)
    genes_needed = sorted(set(foc["gene"]) | set(_METAGENE_GENES)
                          | set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
                          | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", [])))
    bundles = _state_bundles(genes_needed)
    states = list(bundles)
    cov = {s: _state_covariates(s, rna, hs, cov_source) for s, (mir, rna) in bundles.items()}
    cov_desc = {"estimate": "ESTIMATE immune+stromal + prolif",
                "estimate_epi": "ESTIMATE immune+stromal + epithelial sig + prolif",
                "metagene": "prolif/epi/immune/stroma metagenes"}[cov_source]
    print(f"[coupling] {len(foc)} focus edges; states={states}; covariates={cov_desc}")

    from mirna_hallmark.pressure_engine import compute_edge_pressure_map

    pair_keys = set(zip(foc["miRNA"], foc["gene"]))
    edges_sub = edges.loc[
        edges.apply(lambda r: (r["miRNA"], r["gene"]) in pair_keys, axis=1)
    ].copy()
    states_mir = _states_mirna_for_pressure(impute_gtex=True)
    mir_for_pressure = _mirna_for_edge_pressure(bundles, states_mir)
    press_cache: Dict[str, Dict[tuple[str, str], pd.Series]] = {}
    for s in states:
        mir_mat = mir_for_pressure.get(s)
        if mir_mat is None or mir_mat.empty:
            continue
        press_cache[s] = compute_edge_pressure_map(
            edges_sub,
            mir_mat,
            genes=genes_needed,
            expr_mode=C.PRESSURE_ATTRIBUTION_EXPR_MODE,  # type: ignore[arg-type]
            target_norm=C.PRESSURE_TARGET_NORM,  # type: ignore[arg-type]
        )

    rows: List[dict] = []
    for _, e in foc.iterrows():
        arm, g = e["miRNA"], e["gene"]
        row = {"miRNA": arm, "gene": g, "evidence_score": float(e["evidence_score"])}
        key = (str(arm), str(g))
        for s in states:
            mir, rna = bundles[s]
            raw_rho, raw_p, adj_rho, adj_p, n = _couple_state(mir, rna, cov[s], arm, g)
            row[f"raw_rho_{s}"] = raw_rho
            row[f"adj_rho_{s}"] = adj_rho
            row[f"adj_p_{s}"] = adj_p
            row[f"n_{s}"] = n
            pser = press_cache.get(s, {}).get(key)
            _rp, _pp, adj_pr, adj_pp, npres = _couple_state_pressure(rna, cov[s], g, pser)
            row[f"adj_rho_pressure_{s}"] = adj_pr
            row[f"adj_p_pressure_{s}"] = adj_pp
            row[f"n_pressure_{s}"] = npres
        rows.append(row)
    df = pd.DataFrame(rows)
    for s in states:
        if f"adj_p_{s}" in df.columns:
            df[f"adj_q_{s}"] = bh_fdr(df[f"adj_p_{s}"].fillna(1.0).values)
        if f"adj_p_pressure_{s}" in df.columns:
            df[f"adj_q_pressure_{s}"] = bh_fdr(df[f"adj_p_pressure_{s}"].fillna(1.0).values)
    df["archetype"] = df.apply(_classify, axis=1)
    df["archetype_pressure"] = df.apply(_classify_pressure, axis=1)
    df["composition_sensitive"] = df.apply(lambda r: _composition_sensitive(r, states), axis=1)
    df["composition_sensitive_pressure"] = df.apply(
        lambda r: _composition_sensitive(r, states, rho_prefix="adj_rho_pressure"), axis=1
    )
    if {"adj_rho_pressure_tumor", "adj_rho_pressure_gtex"}.issubset(df.columns):
        df["delta_tumor_gtex_pressure"] = (
            df["adj_rho_pressure_tumor"] - df["adj_rho_pressure_gtex"]
        )
    sort_cols = [c for c in ["adj_rho_tumor"] if c in df.columns]
    df = df.sort_values(sort_cols or ["evidence_score"]).reset_index(drop=True)

    tag = f"{(gene or mirna or 'focus').replace('/', '_')}_{cov_source}"
    path = out_dir / f"coupling_trajectories_{tag}.tsv"
    df.to_csv(path, sep="\t", index=False)
    counts = df["archetype"].value_counts().to_dict()
    cov_cols = {
        "estimate": ["ImmuneScore", "StromalScore", "prolif"],
        "estimate_epi": ["ImmuneScore", "StromalScore", "epi", "prolif"],
        "metagene": ["prolif", "epi", "immune", "stroma"],
    }[cov_source]
    manifest = {"module": "mirna_hallmark.analyses.cross_state.cross_state_coupling",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "focus": tag, "cov_source": cov_source, "n_edges": len(df), "states": states,
                "covariates": cov_cols,
                "abundance_coupling": "partial Spearman(miRNA abundance, target expr)",
                "pressure_coupling": (
                    f"partial Spearman(edge pressure c(m,g), target expr); "
                    f"expr_mode={C.PRESSURE_ATTRIBUTION_EXPR_MODE}; target_norm={C.PRESSURE_TARGET_NORM}; "
                    "GTEx miRNA imputed for pressure (collapse-repair)"
                ),
                "rho_floor": RHO_FLOOR, "q_alpha": Q_ALPHA,
                "archetype_counts": counts,
                "archetype_pressure_counts": df["archetype_pressure"].value_counts().to_dict(),
            }
    (out_dir / f"method_manifest_{tag}.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[coupling] ({cov_source}) archetype counts: {counts}")
    print(f"[coupling] ({cov_source}) composition-sensitive edges: {int(df['composition_sensitive'].sum())}/{len(df)}")
    if not quiet:
        _print_stories(df, states)
    return df


def compare_schemes(
    *, out_dir: Path = OUT_DIR, gene: Optional[str] = None, mirna: Optional[str] = None,
    sources: Sequence[str] = ("estimate", "estimate_epi", "metagene"),
) -> pd.DataFrame:
    """Run the focus set under each covariate scheme; quantify archetype concordance."""
    dfs = {s: run(out_dir=out_dir, gene=gene, mirna=mirna, cov_source=s, quiet=True) for s in sources}
    base = None
    for s in sources:
        keep = dfs[s][["miRNA", "gene", "archetype", "composition_sensitive"]].rename(
            columns={"archetype": f"arch_{s}", "composition_sensitive": f"compsens_{s}"})
        base = keep if base is None else base.merge(keep, on=["miRNA", "gene"], how="outer")
    arch_cols = [f"arch_{s}" for s in sources]
    base["n_distinct_calls"] = base[arch_cols].nunique(axis=1)
    base["all_agree"] = base["n_distinct_calls"].eq(1)

    n = len(base)
    print(f"\n=== archetype concordance across schemes ({n} focus edges) ===")
    print(f"all three agree: {int(base['all_agree'].sum())}/{n} "
          f"({100*base['all_agree'].mean():.0f}%)")
    for i in range(len(sources)):
        for j in range(i + 1, len(sources)):
            a, b = sources[i], sources[j]
            same = (base[f"arch_{a}"] == base[f"arch_{b}"]).sum()
            print(f"  {a} vs {b}: {same}/{n} same ({100*same/n:.0f}%)")
    # how many edges leave/enter the "repressor" calls (the biologically actionable set)
    rep = {"constitutive_repressor", "tumor_acquired_repressor", "repression_replaces_coexpr"}
    for s in sources:
        base[f"is_rep_{s}"] = base[f"arch_{s}"].isin(rep)
    print("\nrepressor-call counts (constitutive+tumor_acquired+repl_coexpr):")
    for s in sources:
        print(f"  {s}: {int(base[f'is_rep_{s}'].sum())}")

    # edges whose call changes specifically between estimate and estimate_epi (epi axis)
    if {"arch_estimate", "arch_estimate_epi"}.issubset(base.columns):
        chg = base.loc[base["arch_estimate"] != base["arch_estimate_epi"]]
        print(f"\nestimate -> estimate_epi changed {len(chg)} edges (epithelial-axis sensitivity):")
        if len(chg):
            print(chg.head(15)[["miRNA", "gene", "arch_estimate", "arch_estimate_epi"]].to_string(index=False))

    tag = (gene or mirna or "focus").replace("/", "_")
    path = out_dir / f"archetype_scheme_comparison_{tag}.tsv"
    base.to_csv(path, sep="\t", index=False)
    print(f"\n[compare] wrote {path}")
    return base


def _print_stories(df: pd.DataFrame, states: Sequence[str]) -> None:
    rcols = ["miRNA", "gene", "archetype", "composition_sensitive"]
    rcols += [f"raw_rho_{s}" for s in states] + [f"adj_rho_{s}" for s in states]
    rcols = [c for c in rcols if c in df.columns]
    for arch, title in [
        ("tumor_acquired_repressor", "TUMOR-ACQUIRED repressors (silent in normal, repressive in tumor)"),
        ("constitutive_repressor", "CONSTITUTIVE repressors (repressive in normal + tumor)"),
        ("normal_coexpression_only", "NORMAL CO-EXPRESSION only (lineage program, not regulation)"),
        ("repression_replaces_coexpr", "REPRESSION REPLACES CO-EXPRESSION (+ in normal -> - in tumor)"),
    ]:
        sub = df.loc[df["archetype"].eq(arch)]
        if sub.empty:
            continue
        key = "adj_rho_tumor" if "adj_rho_tumor" in sub.columns else rcols[-1]
        sub = sub.reindex(sub[key].abs().sort_values(ascending=False).index)
        print(f"\n[{arch}] {title}  (n={len(sub)}):")
        print(sub.head(10)[rcols].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene", type=str, default=None)
    ap.add_argument("--mirna", type=str, default=None)
    ap.add_argument("--covariates", choices=["estimate", "estimate_epi", "metagene"],
                    default="estimate_epi",
                    help="composition covariate source (default: ESTIMATE immune+stromal + epithelial sig + prolif). "
                         "Use estimate or metagene for sensitivity audits.")
    ap.add_argument("--compare", action="store_true",
                    help="run all three covariate schemes and quantify archetype concordance")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.compare:
        compare_schemes(out_dir=args.out_dir, gene=args.gene, mirna=args.mirna)
    else:
        run(out_dir=args.out_dir, gene=args.gene, mirna=args.mirna, cov_source=args.covariates)


if __name__ == "__main__":
    main()
