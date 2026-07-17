"""Deepen GSE270497 EV miRNA analysis: novel plasma indicators + TCGA linkage.

Builds on ``ev_mirna_screening`` outputs and TCGA tumor spine:

1. **Novel screening+ indicators** — significant EV elevation excluding well-known
   circulating oncomiRs / panel arms.
2. **TCGA joint table** — EV Δ vs tissue median log2(RPM+1), CN→expr partial ρ,
   Hallmark high-evidence target burden.
3. **Quadrant classification** — export vs dosage-concordance joint patterns.
4. **Scatter figure** — EV Δ vs TCGA partial ρ (novel vs known highlighted).

Outputs under ``mirna_hallmark/output/plasma_ev/deep/``.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import load_ev_log2_matrix, load_ev_samples
from pipeline.config import PANEL_MIRNA_TIER_ARM_IDS, PATHS
from pipeline.miRNA_exp.loaders import map_panel_arms_to_mimat_ids
from pipeline.sample_ids import tcga_sample_type_two_digit

DEEP_DIR = C.PLASMA_EV_DIR / "deep"

# Literature-common circulating breast / cancer miRNAs (exclude from "novel").
WELL_KNOWN_CIRCULATING_PATTERNS: Tuple[re.Pattern[str], ...] = tuple(
    re.compile(p, re.I)
    for p in (
        r"^hsa-miR-21",
        r"^hsa-miR-200[abc]",
        r"^hsa-let-7",
        r"^hsa-miR-155",
        r"^hsa-miR-10b",
        r"^hsa-miR-125[ab]",
        r"^hsa-miR-145",
        r"^hsa-miR-148a",
        r"^hsa-miR-34[abc]",
        r"^hsa-miR-126",
        r"^hsa-miR-210",
        r"^hsa-miR-1246",
        r"^hsa-miR-9",
        r"^hsa-miR-223",
        r"^hsa-miR-17",
        r"^hsa-miR-92[ab]",
        r"^hsa-miR-18a",
        r"^hsa-miR-19[ab]",
        r"^hsa-miR-29[abc]",
        r"^hsa-miR-93",
        r"^hsa-miR-196a",
        r"^hsa-miR-205",
        r"^hsa-miR-122",
        r"^hsa-miR-486",
        r"^hsa-miR-451a",
    )
)

PANEL_SET = set(PANEL_MIRNA_TIER_ARM_IDS)


def _is_well_known(arm: str) -> bool:
    return any(p.search(arm) for p in WELL_KNOWN_CIRCULATING_PATTERNS)


def _load_ev_contrasts(path: Optional[Path] = None) -> pd.DataFrame:
    p = Path(path or C.PLASMA_EV_DIR / "gse270497_screening_arm_contrasts.tsv")
    if not p.is_file():
        raise FileNotFoundError(f"Missing {p}; run: .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.ev_mirna_screening")
    return pd.read_csv(p, sep="\t")


def _load_tcga_concordance() -> pd.DataFrame:
    p = C.MIRNA_LOCUS_CNV_DIR / "mirna_cnv_expr_concordance.tsv"
    if not p.is_file():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def tcga_primary_median_log2rpm(arms: Sequence[str]) -> pd.DataFrame:
    """Map mature arms → MIMAT and return cohort median log2(RPM+1) in primary tumors."""
    tcga_path = PATHS.mirna_expression
    if not tcga_path.is_file():
        return pd.DataFrame(columns=["arm", "tcga_primary_median_log2rpm", "mimat"])

    header = pd.read_csv(tcga_path, sep="\t", nrows=0)
    cols = [c for c in header.columns if c != PATHS.mirna_gene_col]
    primary_cols = [c for c in cols if tcga_sample_type_two_digit(c) == "01"]
    if not primary_cols:
        return pd.DataFrame(columns=["arm", "tcga_primary_median_log2rpm", "mimat"])

    usecols = [PATHS.mirna_gene_col] + primary_cols
    tcga = pd.read_csv(tcga_path, sep="\t", usecols=usecols)
    mimat_ids = tcga[PATHS.mirna_gene_col].astype(str).tolist()
    arm_map = map_panel_arms_to_mimat_ids(arms, mimat_ids)

    loci = pd.read_csv(
        PATHS.mirna_mature_loci_csv,
        usecols=["mirbase_mature_id", "mature_accession"],
    )
    loci["mirbase_mature_id"] = loci["mirbase_mature_id"].astype(str).str.strip()
    loci["mature_accession"] = loci["mature_accession"].astype(str).str.strip()
    avail = set(mimat_ids)
    label_to_mimat = {
        lab: mim
        for lab, mim in zip(loci["mirbase_mature_id"], loci["mature_accession"])
        if mim in avail and lab
    }

    tcga_idx = tcga.set_index(PATHS.mirna_gene_col)
    rows = []
    for arm in arms:
        mimats = arm_map.get(arm) or []
        if not mimats:
            mim = label_to_mimat.get(arm)
            if mim:
                mimats = [mim]
        if not mimats:
            rows.append({"arm": arm, "tcga_primary_median_log2rpm": np.nan, "mimat": ""})
            continue
        sub = tcga_idx.loc[[m for m in mimats if m in tcga_idx.index]]
        med = float(sub[primary_cols].astype(float).median(axis=1).median()) if not sub.empty else np.nan
        rows.append({"arm": arm, "tcga_primary_median_log2rpm": med, "mimat": mimats[0]})
    return pd.DataFrame(rows)


def _hallmark_arm_summary(arms: Iterable[str]) -> pd.DataFrame:
    edges = D.load_hallmark_edges()
    he = D.high_evidence_edges(edges)
    he = he.loc[he["miRNA"].isin(arms)].copy()
    if he.empty:
        return pd.DataFrame(
            columns=[
                "arm",
                "hallmark_highev_target_n",
                "hallmark_highev_sets_n",
                "hallmark_highev_sets",
            ]
        )
    agg = (
        he.groupby("miRNA")
        .agg(
            hallmark_highev_target_n=("gene", "nunique"),
            hallmark_highev_sets_n=("hallmark_set", "nunique"),
            hallmark_highev_sets=("hallmark_set", lambda s: ";".join(sorted(set(s))[:8])),
        )
        .reset_index()
        .rename(columns={"miRNA": "arm"})
    )
    return agg


def build_tcga_linkage(
    contrasts: pd.DataFrame,
    *,
    delta_col: str = "delta_median_log2tpm_pos_minus_neg",
    q_col: str = "mannwhitney_q",
) -> pd.DataFrame:
    """Merge EV contrasts with TCGA tissue abundance, CN→expr concordance, Hallmark edges."""
    if delta_col not in contrasts.columns:
        raise KeyError(f"contrasts missing delta column {delta_col!r}")

    arms = contrasts["arm"].astype(str).tolist()
    tcga_med = tcga_primary_median_log2rpm(arms)
    conc = _load_tcga_concordance()
    hall = _hallmark_arm_summary(arms)

    out = contrasts.merge(tcga_med, on="arm", how="left")
    if not conc.empty:
        keep = [c for c in ("arm", "n", "spearman_rho", "partial_rho", "partial_p", "partial_q") if c in conc.columns]
        out = out.merge(conc[keep], on="arm", how="left")
    out = out.merge(hall, on="arm", how="left")

    out["well_known_circulating"] = out["arm"].map(_is_well_known)
    out["in_panel"] = out["arm"].isin(PANEL_SET)
    out["ev_case_higher_sig"] = (out[delta_col] > 0) & (out[q_col] < C.FDR_ALPHA)

    # Quadrants for arms significant in EV (case higher).
    def _quadrant(r: pd.Series) -> str:
        if not r["ev_case_higher_sig"]:
            return "not_ev_sig_case_higher"
        tcga_conc = (
            pd.notna(r.get("partial_rho"))
            and pd.notna(r.get("partial_q"))
            and r["partial_q"] < C.FDR_ALPHA
            and r["partial_rho"] > 0
        )
        tcga_abundant = pd.notna(r.get("tcga_primary_median_log2rpm")) and r["tcga_primary_median_log2rpm"] >= 8.0
        if tcga_conc and tcga_abundant:
            return "ev_and_tissue_dosage"
        if tcga_conc and not tcga_abundant:
            return "ev_plus_dosage_low_tissue"
        if not tcga_conc and tcga_abundant:
            return "ev_export_high_tissue"
        return "ev_circulation_led"

    out["tcga_joint_quadrant"] = out.apply(_quadrant, axis=1)
    sort_cols = [q_col, delta_col]
    return out.sort_values(sort_cols, ascending=[True, False])


def novel_positive_indicators(
    linkage: pd.DataFrame,
    *,
    delta_col: str = "delta_median_log2tpm_pos_minus_neg",
    min_delta: float = 0.8,
    max_well_known: bool = False,
) -> pd.DataFrame:
    """Less-known case-higher EV arms passing novelty + effect-size filters."""
    m = (
        linkage["ev_case_higher_sig"]
        & ~linkage["in_panel"]
        & (linkage[delta_col] >= min_delta)
    )
    if not max_well_known:
        m &= ~linkage["well_known_circulating"]
    out = linkage.loc[m].copy()
    out["novelty_rank"] = np.arange(1, len(out) + 1)
    return out


def cohort_level_correlations(
    linkage: pd.DataFrame,
    *,
    delta_col: str = "delta_median_log2tpm_pos_minus_neg",
) -> pd.DataFrame:
    """Spearman across arms: EV Δ vs TCGA tissue / concordance metrics."""
    from scipy.stats import spearmanr

    rows = []
    x = linkage[delta_col].astype(float)
    for label, ycol in (
        ("tcga_primary_median_log2rpm", "tcga_primary_median_log2rpm"),
        ("tcga_partial_rho", "partial_rho"),
        ("tcga_spearman_rho", "spearman_rho"),
        ("hallmark_highev_target_n", "hallmark_highev_target_n"),
    ):
        if ycol not in linkage.columns:
            continue
        y = pd.to_numeric(linkage[ycol], errors="coerce")
        mask = x.notna() & y.notna()
        if mask.sum() < 10:
            continue
        rho, p = spearmanr(x[mask], y[mask])
        rows.append({"comparison": f"ev_delta_vs_{label}", "n_arms": int(mask.sum()), "spearman_rho": rho, "p": p})

    for subset, name in (
        (linkage.loc[~linkage["well_known_circulating"]], "non_well_known"),
        (linkage.loc[linkage["ev_case_higher_sig"]], "ev_sig_case_higher"),
        (novel_positive_indicators(linkage, delta_col=delta_col), "novel_positive"),
    ):
        if subset.empty:
            continue
        xs = subset[delta_col].astype(float)
        ys = pd.to_numeric(subset["partial_rho"], errors="coerce")
        mask = xs.notna() & ys.notna()
        if mask.sum() < 8:
            continue
        rho, p = spearmanr(xs[mask], ys[mask])
        rows.append(
            {
                "comparison": f"ev_delta_vs_partial_rho_{name}",
                "n_arms": int(mask.sum()),
                "spearman_rho": rho,
                "p": p,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        m = out["p"].notna()
        out.loc[m, "q"] = S.bh_fdr(out.loc[m, "p"].values)
    return out


def plot_ev_vs_tcga_scatter(
    linkage: pd.DataFrame,
    out_path: Path,
    *,
    delta_col: str = "delta_median_log2tpm_pos_minus_neg",
    ylabel: str = "EV case − control Δ median log₂ abundance",
    title: str = "Plasma exosome vs tumor dosage concordance (sig. EV+ arms)",
) -> None:
    import matplotlib.pyplot as plt

    sub = linkage.loc[linkage["ev_case_higher_sig"]].copy()
    if sub.empty or "partial_rho" not in sub.columns:
        return

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for label, mask, color, marker in (
        ("novel", sub["novel_candidate"], "#C44E52", "o"),
        ("well-known", sub["well_known_circulating"], "#4C72B0", "s"),
        ("other", ~(sub["novel_candidate"] | sub["well_known_circulating"]), "#AAAAAA", "."),
    ):
        if "novel_candidate" not in sub.columns:
            break
        pts = sub.loc[mask]
        if pts.empty:
            continue
        ax.scatter(
            pts["partial_rho"],
            pts[delta_col],
            c=color,
            label=label,
            alpha=0.75,
            s=30 if marker != "." else 12,
            marker=marker,
        )

    # Label top novel points
    if "novel_candidate" in sub.columns:
        top = sub.loc[sub["novel_candidate"]].nlargest(8, delta_col)
        for r in top.itertuples():
            ax.annotate(
                r.arm.replace("hsa-miR-", ""),
                (r.partial_rho, getattr(r, delta_col)),
                fontsize=7,
                alpha=0.9,
            )

    ax.axhline(0, color="k", lw=0.5, alpha=0.4)
    ax.axvline(0, color="k", lw=0.5, alpha=0.4)
    ax.set_xlabel("TCGA primary CN→expr partial Spearman ρ (CPE+HRD adj.)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _write_cohort_deep_outputs(
    *,
    cohort: str,
    contrasts: pd.DataFrame,
    out_dir: Path,
    delta_col: str,
    q_col: str,
    min_novel_delta: float,
    ylabel: str,
    title: str,
) -> pd.DataFrame:
    linkage = build_tcga_linkage(contrasts, delta_col=delta_col, q_col=q_col)
    novel = novel_positive_indicators(linkage, delta_col=delta_col, min_delta=min_novel_delta)
    linkage["novel_candidate"] = linkage["arm"].isin(set(novel["arm"]))

    prefix = cohort.lower()
    linkage.to_csv(out_dir / f"{prefix}_tcga_linkage_all_arms.tsv", sep="\t", index=False)
    novel.to_csv(out_dir / f"{prefix}_novel_positive_indicators.tsv", sep="\t", index=False)

    quads = (
        linkage.loc[linkage["ev_case_higher_sig"]]
        .groupby("tcga_joint_quadrant", as_index=False)
        .agg(n_arms=("arm", "count"), median_ev_delta=(delta_col, "median"))
        .sort_values("n_arms", ascending=False)
    )
    quads.to_csv(out_dir / f"{prefix}_tcga_joint_quadrants.tsv", sep="\t", index=False)

    cors = cohort_level_correlations(linkage, delta_col=delta_col)
    cors.to_csv(out_dir / f"{prefix}_cohort_arm_correlations.tsv", sep="\t", index=False)

    plot_ev_vs_tcga_scatter(
        linkage,
        out_dir / "figures" / f"{prefix}_ev_delta_vs_tcga_partial_rho.png",
        delta_col=delta_col,
        ylabel=ylabel,
        title=title,
    )
    return linkage


def run(
    *,
    out_dir: Path = DEEP_DIR,
    contrasts_path: Optional[Path] = None,
    min_novel_delta: float = 0.8,
    cohort: str = "gse270497",
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    if cohort == "gse197020":
        contrasts_path = Path(
            contrasts_path or C.PLASMA_EV_DIR / "gse197020_case_arm_contrasts.tsv"
        )
        if not contrasts_path.is_file():
            raise FileNotFoundError(
                f"Missing {contrasts_path}; run: .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.ev_mirna_replication"
            )
        contrasts = pd.read_csv(contrasts_path, sep="\t")
        delta_col = "delta_median_log2_case_minus_control"
        ylabel = "EV cancer minus normal Δ median log₂(CPM+1)"
        title = "GSE197020 exosome vs tumor dosage concordance (sig. cancer-higher arms)"
        clinical_note = (
            "Per-sample stage/ER/HER2 not deposited in GEO; disease + author ML train/test only — "
            "see data/miRNA/ev/GSE197020/samples.tsv"
        )
    else:
        contrasts = _load_ev_contrasts(contrasts_path)
        delta_col = "delta_median_log2tpm_pos_minus_neg"
        ylabel = "EV screening+ minus − Δ median log₂(TPM+1)"
        title = "Plasma exosome vs tumor dosage concordance (sig. EV+ arms)"
        clinical_note = (
            "Per-sample clinical covariates not deposited; see data/miRNA/ev/GSE270497/clinical_aggregate.tsv"
        )

    linkage = _write_cohort_deep_outputs(
        cohort=cohort,
        contrasts=contrasts,
        out_dir=out_dir,
        delta_col=delta_col,
        q_col="mannwhitney_q",
        min_novel_delta=min_novel_delta,
        ylabel=ylabel,
        title=title,
    )
    novel = novel_positive_indicators(linkage, delta_col=delta_col, min_delta=min_novel_delta)

    print(f"[ev_mirna_deep] {cohort}: novel positive indicators (Δ≥{min_novel_delta}): {len(novel)}")
    if not novel.empty:
        show = novel.head(8)[
            ["arm", delta_col, "mannwhitney_q", "partial_rho", "tcga_primary_median_log2rpm", "tcga_joint_quadrant"]
        ]
        print(show.to_string(index=False))
    cors_path = out_dir / f"{cohort.lower()}_cohort_arm_correlations.tsv"
    if cors_path.is_file():
        cors = pd.read_csv(cors_path, sep="\t")
        if not cors.empty:
            print(f"[ev_mirna_deep] {cohort} arm-level correlations:")
            print(cors.to_string(index=False))

    manifest_path = out_dir / "method_manifest.json"
    manifest: Dict[str, object] = {}
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(
        {
            "module": "mirna_hallmark.analyses.dcis_ev.ev_mirna_deep",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            cohort: {
                "n_arms_linkage": int(len(linkage)),
                "n_novel_positive": int(len(novel)),
                "min_novel_delta": min_novel_delta,
                "clinical_note": clinical_note,
            },
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[ev_mirna_deep] wrote {cohort} outputs under {out_dir}")
    return linkage


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=DEEP_DIR)
    ap.add_argument("--contrasts", type=Path, default=None)
    ap.add_argument("--min-novel-delta", type=float, default=0.8)
    ap.add_argument(
        "--cohort",
        choices=("gse270497", "gse197020", "both"),
        default="gse270497",
        help="EV cohort for TCGA deep linkage (both runs screening then replication cohort)",
    )
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.cohort == "both":
        run(out_dir=args.out_dir, contrasts_path=args.contrasts, min_novel_delta=args.min_novel_delta, cohort="gse270497")
        run(out_dir=args.out_dir, min_novel_delta=args.min_novel_delta, cohort="gse197020")
    else:
        run(
            out_dir=args.out_dir,
            contrasts_path=args.contrasts,
            min_novel_delta=args.min_novel_delta,
            cohort=args.cohort,
        )


if __name__ == "__main__":
    main()
