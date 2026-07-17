"""Copy-number landscape of the *miRNA universe* (locus / arm / family).

This is a first-class deliverable: the entire miRNA universe (not just Hallmark
targets) is characterized for ASCAT3 copy number and stratified by the
subproject strata (PAM50, TNBC/Lehmann, Thornsson immune, stage). We then test
whether locus dosage propagates to mature-arm expression (CNV->expression
concordance), which qualifies the downstream miRNA-pressure interpretation.

We reuse the heavy machinery from
``analysis.cohort_landscapes.cnv.dosage_landscape_cnv`` (entity catalog,
ASCAT3 segment overlap by default, per-stratum summary) rather
than duplicating it, with ``panel_arm_ids=None`` so the full miRNA universe is
covered. The long sample x entity CNV table is cached for reuse.

Outputs (``output/mirna_locus_cnv/``):
- ``tables/sample_entity_cnv.tsv.gz``      -- cached long sample x entity CNV
- ``mirna_cnv_by_stratum.tsv``             -- locus/arm/family CN summary per stratum
- ``mirna_cnv_expr_concordance.tsv``       -- per-arm Spearman/Pearson(CN, arm expression)
- ``mirna_cnv_expr_concordance_scope_audit.tsv`` -- full cohort vs CN>2 / gain|amp sensitivity
- ``figures/mirna_cnv_expr_concordance_rho_boxplot.png`` -- marginal vs partial ρ distribution
- ``concordance_top_arm_target_join.tsv``  -- miRTarBase target-resolved joins for top CN↔expr concordant arms (repressor prior: negative miRNA↔target RNA)
- ``maps/``  -- locus + MIMAT paralog genome maps (Δ vs diploid / weighted healthy Σ2·w_locus)
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S


def _build_long_cnv(out_dir: Path, *, force: bool = False) -> pd.DataFrame:
    """Long sample x entity miRNA CNV (locus+arm+family), with subproject strata."""
    cache = out_dir / "tables" / "sample_entity_cnv.tsv.gz"
    if cache.exists() and not force:
        print(f"[mirna_locus_cnv] reusing cached long CNV: {cache}")
        return pd.read_csv(cache, sep="\t", low_memory=False)

    from analysis.cohort_landscapes.cnv.dosage_landscape_cnv import (
        _add_cn_state_flags,
        assemble_mirna_cnv_long_table,
        build_entity_catalog,
        load_cohort_mirna_cnv_from_segments,
        load_mirna_family_map,
        resolve_mirna_cnv_source,
    )

    family_map = load_mirna_family_map(C.MIR_FAMILY_INFO)
    catalog = build_entity_catalog(
        loci_path=C.MIRNA_PRECURSOR_LOCI,
        arms_path=C.MIRNA_MATURE_LOCI,
        family_map=family_map,
        panel_arm_ids=None,  # full miRNA universe
    )
    cn_mode = resolve_mirna_cnv_source()
    if cn_mode == "ascat3_segment":
        print(f"[mirna_locus_cnv] extracting cohort CNV from ASCAT3 segments ({cn_mode}) ...")
        raw_long, _ = load_cohort_mirna_cnv_from_segments(
            catalog,
            cnv_dir=C.CNV_DIR,
            manifest_path=C.CNV_MANIFEST,
            primary_only=True,
        )
    else:
        from analysis.cohort_landscapes.cnv.dosage_landscape_cnv import (
            _per_sample_ascat3_paths,
            _reference_gene_table,
            build_entity_ensg_map,
            load_cohort_mirna_cnv,
            resolve_ascat3_cohort_source,
        )

        src = resolve_ascat3_cohort_source(C.CNV_GENE_TABLES_DIR)
        ref_path = src["path"] if src.get("mode") == "combined" else (
            src.get("paths", _per_sample_ascat3_paths(C.CNV_GENE_TABLES_DIR))[0]
        )
        ref = _reference_gene_table(Path(ref_path))
        entity_ensg = build_entity_ensg_map(catalog, ref)
        if entity_ensg.empty:
            raise ValueError("Empty entity->ENSG map; check miRNA coordinates vs CNV gene table")

        print(f"[mirna_locus_cnv] extracting cohort CNV for {len(entity_ensg)} mapped entities ...")
        raw_long, _ = load_cohort_mirna_cnv(
            entity_ensg, gene_tables_dir=C.CNV_GENE_TABLES_DIR, primary_only=True
        )
    long = assemble_mirna_cnv_long_table(raw_long, catalog, loci_path=C.MIRNA_PRECURSOR_LOCI)
    long = long.drop_duplicates(subset=["participant", "entity_id", "entity_level"], keep="first")

    # subproject strata (participant-keyed) + CN flags
    clin = D.load_clinical_strata()
    long = long.merge(clin, on="participant", how="left")
    long = _add_cn_state_flags(long)

    cache.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(cache, sep="\t", index=False, compression="gzip")
    print(f"[mirna_locus_cnv] cached long CNV ({len(long):,} rows) -> {cache}")
    return long


def _participant_arm_cn(long: pd.DataFrame) -> pd.DataFrame:
    """Long participant × arm with copy_number and cn_state."""
    arms = long.loc[long["entity_level"] == "arm"].copy()
    arms["copy_number"] = pd.to_numeric(arms["copy_number"], errors="coerce")
    arms = arms.dropna(subset=["copy_number"])
    arms = arms.drop_duplicates(subset=["participant", "entity_label"], keep="first")
    return arms[["entity_label", "participant", "copy_number", "cn_state"]]


def _iter_arm_cn_expr_merged(
    long: pd.DataFrame,
    *,
    min_n: int = 20,
) -> Iterator[Tuple[str, pd.DataFrame]]:
    """Yield (arm, df) with columns copy_number, cn_state, expr per participant."""
    cn = _participant_arm_cn(long)
    expr = D.load_mirna_arms()
    expr_arms = set(expr.index)
    for arm, sub in cn.groupby("entity_label"):
        if arm not in expr_arms:
            continue
        merged = (
            sub.set_index("participant")[["copy_number", "cn_state"]]
            .join(expr.loc[arm].rename("expr"), how="inner")
            .dropna(subset=["copy_number", "expr"])
        )
        if len(merged) < min_n:
            continue
        yield str(arm), merged


def _concordance_on_merged(
    merged: pd.DataFrame,
    clin: pd.DataFrame,
    conf_num: List[str],
) -> Dict[str, object]:
    cov = clin.reindex(merged.index)[conf_num] if conf_num else None
    stats = S.correlation_pair(merged["expr"], merged["copy_number"], cov)
    return stats


def cnv_expression_concordance(long: pd.DataFrame) -> pd.DataFrame:
    """Per mature-arm: Spearman/Pearson(copy number, arm expression) across participants.

    Confounder-adjusted partial Spearman and partial Pearson (CPE, HRD) where
    available; marginal coefficients are always reported alongside.
    """
    clin = D.load_clinical_strata().set_index("participant")
    conf_num = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]

    rows = []
    for arm, merged in _iter_arm_cn_expr_merged(long):
        st = _concordance_on_merged(merged, clin, conf_num)
        rows.append({
            "arm": arm,
            "n": st["n"],
            "spearman_rho": round(float(st["spearman_rho"]), 4),
            "spearman_p": st["spearman_p"],
            "pearson_r": round(float(st["pearson_r"]), 4),
            "pearson_p": st["pearson_p"],
            "partial_rho": (
                round(float(st["partial_rho"]), 4) if np.isfinite(st["partial_rho"]) else np.nan
            ),
            "partial_p": (
                float(st["partial_p"]) if np.isfinite(st["partial_p"]) else np.nan
            ),
            "partial_pearson_r": (
                round(float(st["partial_pearson_r"]), 4)
                if np.isfinite(st["partial_pearson_r"])
                else np.nan
            ),
            "partial_pearson_p": (
                float(st["partial_pearson_p"]) if np.isfinite(st["partial_pearson_p"]) else np.nan
            ),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["spearman_q"] = S.bh_fdr(df["spearman_p"].values)
        df["pearson_q"] = S.bh_fdr(df["pearson_p"].values)
        partial_p = df["partial_p"].dropna()
        if len(partial_p):
            df["partial_q"] = np.nan
            df.loc[df["partial_p"].notna(), "partial_q"] = S.bh_fdr(
                df.loc[df["partial_p"].notna(), "partial_p"].values
            )
        partial_pp = df["partial_pearson_p"].dropna()
        if len(partial_pp):
            df["partial_pearson_q"] = np.nan
            df.loc[df["partial_pearson_p"].notna(), "partial_pearson_q"] = S.bh_fdr(
                df.loc[df["partial_pearson_p"].notna(), "partial_pearson_p"].values
            )
        df = df.sort_values("spearman_rho", ascending=False)
    return df


def concordance_sample_scope_audit(long: pd.DataFrame, *, min_n: int = 20) -> pd.DataFrame:
    """Compare CN→expr ρ on full cohort vs CN>2-only vs gain|amp-only subsets."""
    clin = D.load_clinical_strata().set_index("participant")
    conf_num = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    rows: List[dict] = []
    for arm, merged in _iter_arm_cn_expr_merged(long, min_n=min_n):
        n_full = len(merged)
        n_diploid = int((merged["copy_number"] <= 2).sum())
        n_gt2 = int((merged["copy_number"] > 2).sum())
        st = merged["cn_state"].astype(str)
        n_gain_amp = int(st.isin(["gain", "amp"]).sum())
        st_full = _concordance_on_merged(merged, clin, conf_num)
        row: Dict[str, object] = {
            "arm": arm,
            "n_full": n_full,
            "n_cn_le_2": n_diploid,
            "n_cn_gt_2": n_gt2,
            "n_gain_or_amp": n_gain_amp,
            "pct_cn_gt_2": round(100.0 * n_gt2 / n_full, 2),
            "spearman_rho_full": round(float(st_full["spearman_rho"]), 4),
            "pearson_r_full": round(float(st_full["pearson_r"]), 4),
            "partial_rho_full": (
                round(float(st_full["partial_rho"]), 4)
                if np.isfinite(st_full["partial_rho"])
                else np.nan
            ),
            "partial_pearson_r_full": (
                round(float(st_full["partial_pearson_r"]), 4)
                if np.isfinite(st_full["partial_pearson_r"])
                else np.nan
            ),
        }
        for label, sub in (
            ("cn_gt_2", merged.loc[merged["copy_number"] > 2]),
            ("gain_or_amp", merged.loc[st.isin(["gain", "amp"])]),
        ):
            if len(sub) < min_n:
                row[f"n_{label}"] = int(len(sub))
                row[f"spearman_rho_{label}"] = np.nan
                row[f"pearson_r_{label}"] = np.nan
                row[f"partial_rho_{label}"] = np.nan
                row[f"partial_pearson_r_{label}"] = np.nan
                continue
            st_sub = _concordance_on_merged(sub, clin, conf_num)
            row[f"n_{label}"] = st_sub["n"]
            row[f"spearman_rho_{label}"] = round(float(st_sub["spearman_rho"]), 4)
            row[f"pearson_r_{label}"] = round(float(st_sub["pearson_r"]), 4)
            row[f"partial_rho_{label}"] = (
                round(float(st_sub["partial_rho"]), 4)
                if np.isfinite(st_sub["partial_rho"])
                else np.nan
            )
            row[f"partial_pearson_r_{label}"] = (
                round(float(st_sub["partial_pearson_r"]), 4)
                if np.isfinite(st_sub["partial_pearson_r"])
                else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def plot_concordance_rho_boxplot(
    conc: pd.DataFrame,
    out_png: Path,
    *,
    fdr_alpha: float = C.FDR_ALPHA,
) -> None:
    """Boxplot of per-arm marginal vs partial Spearman ρ and Pearson r (CN vs expression)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)
    marginal_s = pd.to_numeric(conc["spearman_rho"], errors="coerce").dropna()
    partial_s = pd.to_numeric(conc["partial_rho"], errors="coerce").dropna()
    marginal_p = pd.to_numeric(conc["pearson_r"], errors="coerce").dropna()
    partial_p = pd.to_numeric(conc["partial_pearson_r"], errors="coerce").dropna()
    sig = conc.loc[
        (conc["partial_q"] < fdr_alpha) & (pd.to_numeric(conc["partial_rho"], errors="coerce") > 0)
    ]
    n_sig = len(sig)

    fig, ax = plt.subplots(figsize=(8.5, 5))
    bp = ax.boxplot(
        [marginal_s.values, partial_s.values, marginal_p.values, partial_p.values],
        tick_labels=[
            "Marginal\nSpearman ρ",
            "Partial ρ\n(CPE+HRD)",
            "Marginal\nPearson r",
            "Partial r\n(CPE+HRD)",
        ],
        widths=0.45,
        patch_artist=True,
        showfliers=True,
        medianprops={"color": "black", "linewidth": 1.5},
    )
    colors = ("#9ecae1", "#fdae6b", "#c7e9c0", "#fdd0a2")
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)

    if n_sig:
        sig_partial = pd.to_numeric(sig["partial_rho"], errors="coerce").dropna()
        ax.scatter(
            [2] * len(sig_partial),
            sig_partial.values,
            s=14,
            c="#d62728",
            alpha=0.55,
            zorder=3,
            label=f"partial Spearman q<{fdr_alpha}, ρ>0 (n={n_sig})",
        )
        ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

    ax.axhline(0, color="0.45", linewidth=0.8, linestyle="--")
    ax.set_ylabel("Per-arm correlation (CN vs log2 RPM+1)")
    ax.set_title(
        f"TCGA CN→miRNA expression concordance ({len(conc)} arms tested)\n"
        "all participants per arm; not gain-filtered"
    )
    med_ms, med_ps = float(marginal_s.median()), float(partial_s.median())
    med_mp, med_pp = float(marginal_p.median()), float(partial_p.median())
    ax.text(
        0.02, 0.02,
        f"median Spearman ρ = {med_ms:.3f} / {med_ps:.3f}; "
        f"median Pearson r = {med_mp:.3f} / {med_pp:.3f}",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
    )
    plt.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def run(*, out_dir: Path = C.MIRNA_LOCUS_CNV_DIR, force: bool = False) -> pd.DataFrame:
    from analysis.cohort_landscapes.cnv.dosage_landscape_cnv import summarize_by_strata

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "tables").mkdir(parents=True, exist_ok=True)

    long = _build_long_cnv(out_dir, force=force)

    summary = summarize_by_strata(long, C.STRATUM_SPECS, include_cohort_all=True, min_samples=5)
    summary.to_csv(out_dir / "mirna_cnv_by_stratum.tsv", sep="\t", index=False)
    print(f"[mirna_locus_cnv] stratum summary rows: {len(summary):,}")

    conc = cnv_expression_concordance(long)
    conc.to_csv(out_dir / "mirna_cnv_expr_concordance.tsv", sep="\t", index=False)

    scope_audit = concordance_sample_scope_audit(long)
    scope_audit.to_csv(out_dir / "mirna_cnv_expr_concordance_scope_audit.tsv", sep="\t", index=False)
    fig_dir = out_dir / "figures"
    rho_fig = fig_dir / "mirna_cnv_expr_concordance_rho_boxplot.png"
    if not conc.empty:
        plot_concordance_rho_boxplot(conc, rho_fig)
        n_gt2_tested = int((scope_audit["n_cn_gt_2"] >= 20).sum()) if not scope_audit.empty else 0
        print(
            f"[mirna_locus_cnv] concordance scope audit: "
            f"median pct CN>2 per arm = "
            f"{scope_audit['pct_cn_gt_2'].median():.1f}% "
            f"({n_gt2_tested} arms with n≥20 CN>2 for sensitivity)"
        )
        print(f"[mirna_locus_cnv] wrote {rho_fig}")

    if not conc.empty:
        sig = conc.loc[
            (conc["partial_q"] < C.FDR_ALPHA) & (conc["partial_rho"] > 0)
        ]
        print(
            f"[mirna_locus_cnv] arms with partial CN->expr concordance "
            f"(q<{C.FDR_ALPHA}, partial_rho>0): {len(sig)}/{len(conc)}"
        )
        sig_m = conc.loc[(conc["spearman_q"] < C.FDR_ALPHA) & (conc["spearman_rho"] > 0)]
        print(
            f"[mirna_locus_cnv] (marginal Spearman q<{C.FDR_ALPHA}, rho>0): "
            f"{len(sig_m)}/{len(conc)}"
        )

    from mirna_hallmark.analyses.misc.concordance_target_join import build_concordance_target_join

    target_join = build_concordance_target_join(
        long, conc, top_n_arms=C.CONCORDANCE_TARGET_TOP_N_ARMS
    )
    target_join.to_csv(out_dir / "concordance_top_arm_target_join.tsv", sep="\t", index=False)
    if not target_join.empty:
        print(f"[mirna_locus_cnv] target-resolved join rows: {len(target_join):,} "
              f"({target_join['mirna_arm'].nunique()} arms)")

    from mirna_hallmark.analyses.cnv_locus.mirna_cnv_genome_maps import write_genome_maps

    map_paths = write_genome_maps(long, out_dir)

    from analysis.config import ANALYSIS as A_CFG
    cn_mode = A_CFG.mirna_cnv_dosage.mirna_cnv_source
    n_participants = int(long["participant"].nunique())
    n_primary_vials = int(long["sample_vial"].nunique())
    n_entities = long.groupby("entity_level")["entity_id"].nunique().to_dict()
    manifest = {
        "module": "mirna_hallmark.mirna_locus_cnv",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reuses": "analysis.cohort_landscapes.cnv.dosage_landscape_cnv (panel_arm_ids=None)",
        "cnv_source": f"ASCAT3 {cn_mode}, primary tumors only",
        "n_participants": n_participants,
        "n_primary_sample_vials": n_primary_vials,
        "mimat_cn_aggregate": A_CFG.mirna_cnv_dosage.mimat_cn_aggregate,
        "entity_levels": {"locus": "pre-miRNA hairpin", "arm": "mature arm (MIMAT)", "family": "TargetScan seed family"},
        "n_entities_by_level": {k: int(v) for k, v in n_entities.items()},
        "strata": [layer for _, layer in C.STRATUM_SPECS],
        "concordance": (
            "Spearman + Pearson(median CN, arm log2(RPM+1)); "
            "partial adjusts CPE+HRD; min n=20; all CN states"
        ),
        "concordance_scope_audit": "mirna_cnv_expr_concordance_scope_audit.tsv (full vs CN>2 vs gain|amp)",
        "concordance_rho_boxplot": "figures/mirna_cnv_expr_concordance_rho_boxplot.png",
        "n_concordance_arms": int(len(conc)),
        "n_concordance_target_join_rows": int(len(target_join)),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[mirna_locus_cnv] wrote outputs under {out_dir}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force-rebuild", action="store_true", help="Re-extract cohort CNV even if cached")
    ap.add_argument(
        "--concordance-figures-only",
        action="store_true",
        help="Rebuild concordance TSV, scope audit, and rho boxplot from cached long CNV",
    )
    ap.add_argument("--out-dir", type=Path, default=C.MIRNA_LOCUS_CNV_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.concordance_figures_only:
        cache = args.out_dir / "tables" / "sample_entity_cnv.tsv.gz"
        long = pd.read_csv(cache, sep="\t", low_memory=False)
        conc = cnv_expression_concordance(long)
        conc.to_csv(args.out_dir / "mirna_cnv_expr_concordance.tsv", sep="\t", index=False)
        scope_audit = concordance_sample_scope_audit(long)
        scope_audit.to_csv(
            args.out_dir / "mirna_cnv_expr_concordance_scope_audit.tsv", sep="\t", index=False
        )
        plot_concordance_rho_boxplot(
            conc, args.out_dir / "figures" / "mirna_cnv_expr_concordance_rho_boxplot.png"
        )
        print(f"[mirna_locus_cnv] concordance figures refreshed under {args.out_dir}")
        return
    run(out_dir=args.out_dir, force=args.force_rebuild)


if __name__ == "__main__":
    main()
