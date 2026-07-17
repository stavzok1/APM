"""Cross-tissue-state miRNA + pressure landscape: tumor / NAT / GTEx-healthy.

A *global* (not family-specific) first-pass landscape over the three references we
now have:

  - TCGA-BRCA primary tumor   (sample type 01, RPM)
  - TCGA-BRCA adjacent normal (NAT, sample type 11, RPM) -- same platform as tumor
  - GTEx v10 healthy breast   (small RNA-seq miRNA TPM)  -- orthogonal true-healthy

Three landscape products:

1. ``mirna_arm_de_landscape.tsv`` -- every miRNA arm: tumor vs NAT median + log2FC +
   Mann-Whitney BH-FDR (same-platform DE), GTEx healthy median + detection, and an
   induced / suppressed / stable class. The cohort-wide map of deregulated miRNAs.

2. ``gene_pressure_by_state.tsv`` + ``hallmark_pressure_by_state.tsv`` -- gene-level
   miRNA pressure recomputed in each state under the NO-Z attribution mode
   (``softmax_logrpm``). No-z is essential here: the z-term is standardized *within*
   each state, so z-based pressure means are ~0 in every state and NOT comparable
   across states; the no-z magnitude is. Aggregated to Hallmark programs with the
   tumor-vs-NAT and tumor-vs-GTEx pressure shift.

3. ``mirna_pressure_share_shift.tsv`` -- each arm's share of the TOTAL realized
   pressure budget per state (Σ|contribution| over the Hallmark universe), and how
   that share moves tumor vs NAT / GTEx. Identifies which miRNAs *become* dominant
   regulators in tumor regardless of any single gene.

Cross-platform note: miRNA GTEx (TPM/URS) vs TCGA (RPM/MIMAT) differ in pipeline +
annotation -> GTEx is used for rank/share/sign, NOT absolute pressure magnitude
against TCGA. Tumor-vs-NAT is the fair same-platform contrast.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.cross_state.cross_state_landscape
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from mirna_hallmark import config as C
from mirna_hallmark.family_normal_reference import (
    FAMILY_ARMS,
    GTEX_MIRNA_TPM,
    _gtex_available,
    _gtex_breast_samples,
    _gtex_donor,
    _gtex_urs_to_arm,
    _load_full_mirna,
    _participant,
    _split_types,
)
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import (
    compute_gene_pressure,
    compute_gene_pressure_contributions,
    load_mirtar_edges,
)

OUT_DIR = C.TISSUE_REFERENCE_DIR / "cross_state_landscape"
ATTR_MODE = C.PRESSURE_ATTRIBUTION_EXPR_MODE  # no-z, cross-state comparable
DE_FC = 0.5
DE_Q = 0.05


def _gtex_all_arm_matrix() -> pd.DataFrame:
    """All miRNA arms x GTEx-breast-donor, log2(TPM+1), indexed on **TCGA arm names**.

    Uses MIMAT-primary join (``gtex_mirna_matrix``) so coupling matches tumor/NAT keys
    and ``healthy_anchor`` abundance paths.
    """
    from mirna_hallmark.gtex_mirna_matrix import gtex_tcga_arm_matrix

    return gtex_tcga_arm_matrix()


def _state_matrices() -> Dict[str, pd.DataFrame]:
    full = _load_full_mirna()
    split = _split_types(full)
    states = {"tumor": split["tumor"], "nat": split["nat"]}
    if _gtex_available():
        g = _gtex_all_arm_matrix()
        if not g.empty:
            states["gtex"] = g
    return states


# --------------------------------------------------------------------------- #
# 1. miRNA arm DE landscape
# --------------------------------------------------------------------------- #
def _mirna_de_landscape(states: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    tumor, nat = states["tumor"], states["nat"]
    gtex = states.get("gtex")
    arms = sorted(set(tumor.index) & set(nat.index))
    rows: List[dict] = []
    for arm in arms:
        t = pd.to_numeric(tumor.loc[arm], errors="coerce").dropna()
        n = pd.to_numeric(nat.loc[arm], errors="coerce").dropna()
        if len(t) < 10 or len(n) < 10:
            continue
        try:
            _, p = mannwhitneyu(t, n, alternative="two-sided")
        except ValueError:
            p = np.nan
        g_med = (
            float(pd.to_numeric(gtex.loc[arm], errors="coerce").median())
            if gtex is not None and arm in gtex.index else np.nan
        )
        rows.append({
            "arm": arm,
            "tumor_median": float(t.median()),
            "nat_median": float(n.median()),
            "log2fc_tumor_vs_nat": float(t.median() - n.median()),
            "mwu_p": float(p),
            "gtex_median_tpm": g_med,
            "gtex_detected": bool(pd.notna(g_med) and g_med > 0),
            "is_family_arm": arm in set(FAMILY_ARMS),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    from mirna_hallmark.stats import bh_fdr
    out["mwu_q"] = bh_fdr(out["mwu_p"].fillna(1.0).values)
    out["de_class"] = np.select(
        [
            (out["log2fc_tumor_vs_nat"] >= DE_FC) & (out["mwu_q"] < DE_Q),
            (out["log2fc_tumor_vs_nat"] <= -DE_FC) & (out["mwu_q"] < DE_Q),
        ],
        ["tumor_induced", "tumor_suppressed"],
        default="stable",
    )
    # --- cross-platform-safe trajectory across the 3 states (rank percentiles) ---
    # GTEx is cross-platform vs TCGA, so the GTEx contrasts use within-state
    # percentile RANK (robust) rather than the cross-scale median delta.
    out["tumor_pct"] = out["tumor_median"].rank(pct=True)
    out["nat_pct"] = out["nat_median"].rank(pct=True)
    g_detect = out["gtex_detected"]
    out["gtex_pct"] = np.where(g_detect, out.loc[g_detect, "gtex_median_tpm"].rank(pct=True).reindex(out.index), np.nan)
    out["gtex_to_nat_pct"] = out["nat_pct"] - out["gtex_pct"]
    out["gtex_to_tumor_pct"] = out["tumor_pct"] - out["gtex_pct"]
    s1 = out["nat_pct"] - out["gtex_pct"]   # healthy -> NAT
    s2 = out["tumor_pct"] - out["nat_pct"]  # NAT -> tumor
    net = out["gtex_to_tumor_pct"]
    out["trajectory"] = np.select(
        [
            net.isna(),
            (net >= 0.15) & (s1 >= 0) & (s2 >= 0),
            (net <= -0.15) & (s1 <= 0) & (s2 <= 0),
            net >= 0.15,
            net <= -0.15,
        ],
        ["gtex_absent_or_undetected", "progressive_up", "progressive_down", "net_up", "net_down"],
        default="flat",
    )
    # NAT field-effect: already shifted >=half-way toward tumor from healthy
    out["nat_field_effect"] = (
        (np.sign(s1) == np.sign(net)) & (net.abs() >= 0.1) & (s1.abs() >= 0.5 * net.abs())
    )
    return out.sort_values("log2fc_tumor_vs_nat", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 2. Cross-state pressure (no-z) + Hallmark aggregation
# --------------------------------------------------------------------------- #
def _gene_pressure_by_state(
    states: Dict[str, pd.DataFrame], genes: Sequence[str], edges: pd.DataFrame,
) -> pd.DataFrame:
    out = {}
    for state, mat in states.items():
        pres = compute_gene_pressure(
            genes, edges=edges, mirna=mat, expr_mode=ATTR_MODE,  # type: ignore[arg-type]
            resolve_arms=False,
        )
        out[state] = pres.mean(axis=1) if not pres.empty else pd.Series(dtype=float)
    df = pd.DataFrame(out)
    df.index.name = "gene"
    if {"tumor", "nat"}.issubset(df.columns):
        df["delta_tumor_nat"] = df["tumor"] - df["nat"]
    if "gtex" in df.columns and "tumor" in df.columns:
        df["delta_tumor_gtex"] = df["tumor"] - df["gtex"]
    return df.reset_index()


def _hallmark_pressure(gene_press: pd.DataFrame, hs: HallmarkSets) -> pd.DataFrame:
    gp = gene_press.set_index("gene")
    state_cols = [c for c in ("tumor", "nat", "gtex") if c in gp.columns]
    rows: List[dict] = []
    for hset, members in hs.sets.items():
        present = [g for g in members if g in gp.index]
        if not present:
            continue
        sub = gp.loc[present, state_cols]
        row = {"hallmark": hset, "n_genes": len(present)}
        for c in state_cols:
            row[f"mean_pressure_{c}"] = float(sub[c].mean())
        if {"tumor", "nat"}.issubset(state_cols):
            row["delta_tumor_nat"] = row["mean_pressure_tumor"] - row["mean_pressure_nat"]
        if "gtex" in state_cols:
            row["delta_tumor_gtex"] = row["mean_pressure_tumor"] - row["mean_pressure_gtex"]
        rows.append(row)
    out = pd.DataFrame(rows)
    sort_col = "delta_tumor_nat" if "delta_tumor_nat" in out.columns else "mean_pressure_tumor"
    return out.sort_values(sort_col, ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 3. miRNA pressure-budget share shift
# --------------------------------------------------------------------------- #
def _pressure_share_shift(
    states: Dict[str, pd.DataFrame], genes: Sequence[str], edges: pd.DataFrame,
) -> pd.DataFrame:
    masses = {}
    for state, mat in states.items():
        contrib = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=mat, expr_mode=ATTR_MODE,  # type: ignore[arg-type]
            resolve_arms=False,
        )
        if contrib.empty:
            continue
        mass = contrib.groupby("miRNA")["mean_abs_contribution"].sum()
        masses[state] = mass / mass.sum()  # share of total pressure budget
    if not masses:
        return pd.DataFrame()
    df = pd.DataFrame(masses)
    df.columns = [f"share_{c}" for c in df.columns]
    df.index.name = "miRNA"
    if {"share_tumor", "share_nat"}.issubset(df.columns):
        df["share_shift_tumor_nat"] = df["share_tumor"].fillna(0) - df["share_nat"].fillna(0)
    df["is_family_arm"] = df.index.isin(set(FAMILY_ARMS))
    sort_col = "share_shift_tumor_nat" if "share_shift_tumor_nat" in df.columns else "share_tumor"
    return df.reset_index().sort_values(sort_col, ascending=False).reset_index(drop=True)


_PAM50 = ["LumA", "LumB", "Her2", "Basal"]


def _pressure_share_by_pam50(
    tumor_mat: pd.DataFrame, genes: Sequence[str], edges: pd.DataFrame,
) -> pd.DataFrame:
    """Per-arm share of the tumor pressure budget WITHIN each PAM50 subtype."""
    from mirna_hallmark import data_loaders as D

    clin = D.load_clinical_strata()
    if "PAM50_final" not in clin.columns:
        return pd.DataFrame()
    part_pam = clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
    col_pam = {c: part_pam.get(_participant(c)) for c in tumor_mat.columns}
    shares = {}
    for sub in _PAM50:
        cols = [c for c in tumor_mat.columns if col_pam.get(c) == sub]
        if len(cols) < 20:
            continue
        contrib = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=tumor_mat[cols], expr_mode=ATTR_MODE,  # type: ignore[arg-type]
            resolve_arms=False,
        )
        if contrib.empty:
            continue
        mass = contrib.groupby("miRNA")["mean_abs_contribution"].sum()
        shares[sub] = mass / mass.sum()
    if not shares:
        return pd.DataFrame()
    df = pd.DataFrame(shares)
    df.columns = [f"share_{c}" for c in df.columns]
    df.index.name = "miRNA"
    sub_cols = list(df.columns)
    df["dominant_subtype"] = df[sub_cols].idxmax(axis=1).str.replace("share_", "", regex=False)
    df["share_spread"] = df[sub_cols].max(axis=1) - df[sub_cols].min(axis=1)
    df["is_family_arm"] = df.index.isin(set(FAMILY_ARMS))
    return df.reset_index().sort_values("share_spread", ascending=False).reset_index(drop=True)


def _acquired_pressure_share(
    tumor_mat: pd.DataFrame, genes: Sequence[str], edges: pd.DataFrame, baseline: pd.Series,
) -> pd.DataFrame:
    """Tumor pressure-budget share under the HEALTHY-ANCHORED mode (acquired-vs-healthy)
    next to the within-tumor z-spine share, so the delta = pressure the tumour *gained*."""
    def _share(mode: str, **kw) -> pd.Series:
        c = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=tumor_mat, expr_mode=mode,  # type: ignore[arg-type]
            resolve_arms=False, **kw)
        if c.empty:
            return pd.Series(dtype=float)
        m = c.groupby("miRNA")["mean_abs_contribution"].sum()
        return m / m.sum()
    z = _share("softmax_z_logrpm").rename("share_zspine")
    h = _share(C.PRESSURE_HEALTHY_ANCHOR_MODE, healthy_baseline=baseline).rename("share_acquired")
    df = pd.concat([z, h], axis=1).fillna(0.0)
    if df.empty:
        return df
    df["delta_acquired_vs_zspine"] = df["share_acquired"] - df["share_zspine"]
    df.index.name = "miRNA"
    df["is_family_arm"] = df.index.isin(set(FAMILY_ARMS))
    return df.reset_index().sort_values("delta_acquired_vs_zspine", ascending=False).reset_index(drop=True)


def _acquired_share_by_pam50(
    tumor_mat: pd.DataFrame, genes: Sequence[str], edges: pd.DataFrame, baseline: pd.Series,
) -> pd.DataFrame:
    """Per-arm share of the HEALTHY-ANCHORED (acquired) pressure budget within each PAM50."""
    from mirna_hallmark import data_loaders as D

    clin = D.load_clinical_strata()
    if "PAM50_final" not in clin.columns:
        return pd.DataFrame()
    part_pam = clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
    col_pam = {c: part_pam.get(_participant(c)) for c in tumor_mat.columns}
    shares = {}
    for sub in _PAM50:
        cols = [c for c in tumor_mat.columns if col_pam.get(c) == sub]
        if len(cols) < 20:
            continue
        contrib = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=tumor_mat[cols], expr_mode=C.PRESSURE_HEALTHY_ANCHOR_MODE,  # type: ignore[arg-type]
            resolve_arms=False, healthy_baseline=baseline,
        )
        if contrib.empty:
            continue
        mass = contrib.groupby("miRNA")["mean_abs_contribution"].sum()
        shares[sub] = mass / mass.sum()
    if not shares:
        return pd.DataFrame()
    df = pd.DataFrame(shares)
    df.columns = [f"share_{c}" for c in df.columns]
    df.index.name = "miRNA"
    sub_cols = list(df.columns)
    df["dominant_subtype"] = df[sub_cols].idxmax(axis=1).str.replace("share_", "", regex=False)
    df["share_spread"] = df[sub_cols].max(axis=1) - df[sub_cols].min(axis=1)
    df["is_family_arm"] = df.index.isin(set(FAMILY_ARMS))
    return df.reset_index().sort_values("share_spread", ascending=False).reset_index(drop=True)


def run(*, out_dir: Path = OUT_DIR, gene_limit: Optional[int] = None,
        healthy_anchor: bool = False) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    genes = sorted(hs.universe)
    if gene_limit:
        genes = genes[:gene_limit]

    states = _state_matrices()
    print(f"[cross_state_landscape] states={list(states)}; "
          f"sizes={{ {', '.join(f'{k}:{v.shape[1]}' for k,v in states.items())} }}")

    de = _mirna_de_landscape(states)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    gene_press = _gene_pressure_by_state(states, genes, edges)
    hallmark_press = _hallmark_pressure(gene_press, hs)
    share_shift = _pressure_share_shift(states, genes, edges)
    share_pam50 = _pressure_share_by_pam50(states["tumor"], genes, edges)

    # Healthy-anchored (acquired-vs-healthy) shares are gated OFF until a reliable
    # true-healthy reference (DIANA-miTED breast) is ingested — GTEx-Portal zeros and
    # field-effect NAT both give wrong baselines. Enable with --healthy-anchor.
    acquired_share = pd.DataFrame()
    acquired_pam50 = pd.DataFrame()
    if healthy_anchor:
        from mirna_hallmark.healthy_anchor import load_baseline
        baseline = load_baseline()
        acquired_share = _acquired_pressure_share(states["tumor"], genes, edges, baseline)
        acquired_pam50 = _acquired_share_by_pam50(states["tumor"], genes, edges, baseline)

    de.to_csv(out_dir / "mirna_arm_de_landscape.tsv", sep="\t", index=False)
    gene_press.to_csv(out_dir / "gene_pressure_by_state.tsv", sep="\t", index=False)
    hallmark_press.to_csv(out_dir / "hallmark_pressure_by_state.tsv", sep="\t", index=False)
    share_shift.to_csv(out_dir / "mirna_pressure_share_shift.tsv", sep="\t", index=False)
    if not share_pam50.empty:
        share_pam50.to_csv(out_dir / "mirna_pressure_share_by_pam50.tsv", sep="\t", index=False)
    if not acquired_share.empty:
        acquired_share.to_csv(out_dir / "mirna_acquired_pressure_share.tsv", sep="\t", index=False)
    if not acquired_pam50.empty:
        acquired_pam50.to_csv(out_dir / "mirna_acquired_share_by_pam50.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.analyses.cross_state.cross_state_landscape",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "states": {k: int(v.shape[1]) for k, v in states.items()},
        "attribution_expr_mode": ATTR_MODE,
        "target_norm": C.PRESSURE_TARGET_NORM,
        "n_genes": len(genes),
        "de_thresholds": {"abs_log2fc": DE_FC, "q": DE_Q, "platform": "tumor vs NAT (same RPM)"},
        "cross_platform_note": "miRNA GTEx TPM/URS vs TCGA RPM/MIMAT -> rank/share/sign; mRNA both TPM",
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    if not de.empty:
        cls = de["de_class"].value_counts().to_dict()
        print(f"[landscape] miRNA DE (tumor vs NAT): {cls}")
        print("  top induced:", ";".join(de.head(8)["arm"]))
        print("  top suppressed:", ";".join(de.tail(8)["arm"]))
        if "trajectory" in de.columns:
            traj = de["trajectory"].value_counts().to_dict()
            fe = int(de["nat_field_effect"].sum())
            print(f"[landscape] 3-state trajectory (rank-based): {traj}; NAT field-effect arms: {fe}")
            pu = de.loc[de["trajectory"].eq("progressive_up")].sort_values("gtex_to_tumor_pct", ascending=False)
            print("  progressive_up (healthy->NAT->tumor):", ";".join(pu.head(8)["arm"]))
    if not hallmark_press.empty and "delta_tumor_nat" in hallmark_press.columns:
        print("\n[landscape] Hallmark pressure gain in tumor vs NAT (top 8):")
        print(hallmark_press.head(8)[["hallmark", "mean_pressure_tumor", "mean_pressure_nat", "delta_tumor_nat"]].to_string(index=False))
        print("[landscape] Hallmark pressure loss (bottom 5):")
        print(hallmark_press.tail(5)[["hallmark", "mean_pressure_tumor", "mean_pressure_nat", "delta_tumor_nat"]].to_string(index=False))
    if not share_shift.empty and "share_shift_tumor_nat" in share_shift.columns:
        print("\n[landscape] miRNA pressure-budget share GAIN in tumor (top 10):")
        print(share_shift.head(10)[["miRNA", "share_tumor", "share_nat", "share_shift_tumor_nat"]].to_string(index=False))
    if not share_pam50.empty:
        print("\n[landscape] most PAM50-DIVERGENT pressure carriers (top 10 by share_spread):")
        cols = ["miRNA"] + [c for c in share_pam50.columns if c.startswith("share_Lum") or c in ("share_Her2", "share_Basal")] + ["dominant_subtype", "share_spread"]
        print(share_pam50.head(10)[cols].to_string(index=False))
    if not acquired_share.empty:
        print("\n[landscape] HEALTHY-ANCHORED acquired-pressure share — top GAINERS vs z-spine:")
        print(acquired_share.head(10)[["miRNA", "share_zspine", "share_acquired", "delta_acquired_vs_zspine"]].to_string(index=False))
        print("[landscape] top LOSERS (constitutively high in healthy):")
        print(acquired_share.tail(8)[["miRNA", "share_zspine", "share_acquired", "delta_acquired_vs_zspine"]].to_string(index=False))
    return {"de": de, "gene_pressure": gene_press, "hallmark_pressure": hallmark_press,
            "share_shift": share_shift, "share_pam50": share_pam50,
            "acquired_share": acquired_share, "acquired_pam50": acquired_pam50}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--gene-limit", type=int, default=None, help="cap genes for a fast smoke run")
    ap.add_argument("--healthy-anchor", action="store_true",
                    help="also emit healthy-anchored (acquired-vs-healthy) shares (needs miTED reference)")
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, gene_limit=args.gene_limit, healthy_anchor=args.healthy_anchor)


if __name__ == "__main__":
    main()
