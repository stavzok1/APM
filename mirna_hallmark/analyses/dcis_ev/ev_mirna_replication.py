"""Cross-cohort EV miRNA replication: GSE270497 screening vs GSE197020 cancer vs normal.

Runs case-control contrasts on both ingested plasma-exosome cohorts, computes
per-arm AUC, and joins with TCGA target-combined anticorr summaries for the
miR-29 family (and optional focus arms).

Outputs under ``mirna_hallmark/output/plasma_ev/``:
  ``gse197020_case_arm_contrasts.tsv``
  ``ev_cohort_replication_focus_arms.tsv`` — side-by-side EV + AUC + tissue spine
  ``mir29_ev_tissue_bridge.tsv`` — miR-29a/b/c EV × ECM anticorr summary
  ``figures/ev_cohort_replication_focus_arms.png``

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.ev_mirna_replication
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import (
    FEATURE_COL,
    MIR29_ARMS,
    case_control_arm_contrasts,
    compare_ev_delta_to_tcga_tissue,
    load_ev_log2_matrix,
    load_ev_samples,
    tcga_tissue_median_by_arm,
)
from pipeline.config import PANEL_MIRNA_TIER_ARM_IDS

# Arms to highlight in replication (panel + miR-29 + classic EV markers).
REPLICATION_FOCUS_ARMS: Tuple[str, ...] = tuple(
    dict.fromkeys(
        [
            *MIR29_ARMS,
            "hsa-miR-9-5p",
            "hsa-miR-125b-5p",
            "hsa-miR-148a-3p",
            "hsa-miR-342-5p",
            "hsa-miR-342-3p",
            "hsa-miR-21-5p",
            "hsa-miR-1246",
            "hsa-miR-155-5p",
            "hsa-miR-200c-3p",
        ]
    )
)

ECM_ANTICORR_GENES = frozenset(
    {
        "COL1A1",
        "COL3A1",
        "COL5A2",
        "LOX",
        "LOXL2",
        "MMP2",
        "SPARC",
        "FBN1",
        "SERPINH1",
        "TGFB2",
        "TGFB3",
        "ADAM12",
        "BMP1",
        "VEGFA",
        "PTEN",
        "MCL1",
    }
)


def load_gse197020_matrix(path: Optional[Path] = None) -> pd.DataFrame:
    p = Path(path or C.EV_MIRNA_GSE197020_LOG2_CPM)
    if not p.is_file():
        raise FileNotFoundError(
            f"Missing {p}. Run:\n"
            "  bash scripts/mirna/download_geo_ev_gse197020.sh\n"
            "  .venv/bin/python3 scripts/mirna/ingest_geo_ev_gse197020.py"
        )
    df = pd.read_csv(p, sep="\t")
    feat = df.pop(FEATURE_COL)
    df.index = feat.astype(str).str.strip()
    return df


def load_gse197020_samples(path: Optional[Path] = None) -> pd.DataFrame:
    p = Path(path or C.EV_MIRNA_GSE197020_SAMPLES)
    if not p.is_file():
        raise FileNotFoundError(f"Missing {p}; run ingest_geo_ev_gse197020.py first")
    return pd.read_csv(p, sep="\t")


def _arm_auc(expr: pd.DataFrame, samples: pd.DataFrame, label_col: str) -> Dict[str, float]:
    """Return arm -> AUC (1 = higher in case/boolean-True group)."""
    sample_ids = [c for c in expr.columns if c in set(samples["sample_id"])]
    labels = samples.set_index("sample_id").loc[sample_ids, label_col].astype(bool)
    y = labels.astype(int).values
    out: Dict[str, float] = {}
    for arm in expr.index:
        x = pd.to_numeric(expr.loc[arm, sample_ids], errors="coerce").values
        ok = np.isfinite(x)
        if ok.sum() < 10 or len(np.unique(y[ok])) < 2:
            out[arm] = np.nan
            continue
        out[arm] = float(roc_auc_score(y[ok], x[ok]))
    return out


def _anticorr_summary_for_arms(arms: Sequence[str]) -> pd.DataFrame:
    """Per-arm TCGA gated anticorr hit counts from target_combined outputs."""
    top_path = C.TARGET_COMBINED_ANTICORR_DIR / "target_combined_top_mirnas.tsv"
    by_pam_path = C.TARGET_COMBINED_ANTICORR_DIR / "target_combined_anticorr_by_pam50.tsv"
    if not top_path.is_file() or not by_pam_path.is_file():
        return pd.DataFrame()

    top = pd.read_csv(top_path, sep="\t")
    by_pam = pd.read_csv(by_pam_path, sep="\t")
    by_pam = by_pam[(by_pam["view"] == "gated") & (by_pam["partial_q"] < C.FDR_ALPHA)]
    by_pam = by_pam[by_pam["partial_rho"] < 0]

    rows: List[dict] = []
    for arm in arms:
        genes = set(top.loc[top["miRNA"] == arm, "gene"])
        sub = by_pam[by_pam["gene"].isin(genes)]
        ecm = sub[sub["gene"].isin(ECM_ANTICORR_GENES)]
        rows.append(
            {
                "arm": arm,
                "n_top_driver_genes": len(genes),
                "n_anticorr_gene_x_pam50": len(sub),
                "n_anticorr_basal": int((sub["pam50"] == "Basal").sum()),
                "n_anticorr_lumA": int((sub["pam50"] == "LumA").sum()),
                "n_anticorr_lumB": int((sub["pam50"] == "LumB").sum()),
                "n_anticorr_her2": int((sub["pam50"] == "Her2").sum()),
                "n_ecm_anticorr_hits": len(ecm),
                "median_abs_partial_rho": round(float(sub["partial_rho"].abs().median()), 4)
                if len(sub)
                else np.nan,
                "max_abs_partial_rho": round(float(sub["partial_rho"].abs().max()), 4) if len(sub) else np.nan,
                "top_ecm_targets": ";".join(
                    sorted(ecm.groupby("gene")["partial_rho"].median().nsmallest(5).index.astype(str))
                ),
            }
        )
    return pd.DataFrame(rows)


def build_replication_table(
    g270_contrasts: pd.DataFrame,
    g197_contrasts: pd.DataFrame,
    auc_270: Dict[str, float],
    auc_197: Dict[str, float],
    anticorr: pd.DataFrame,
    arms: Sequence[str],
) -> pd.DataFrame:
    """Side-by-side EV deltas, q-values, AUC, and tissue anticorr spine."""
    c270 = g270_contrasts.set_index("arm")
    c197 = g197_contrasts.set_index("arm")
    rows: List[dict] = []
    for arm in arms:
        r: dict = {"arm": arm}
        if arm in c270.index:
            r.update(
                {
                    "gse270497_delta": c270.loc[arm, "delta_median_log2tpm_pos_minus_neg"],
                    "gse270497_q": c270.loc[arm, "mannwhitney_q"] if "mannwhitney_q" in c270.columns else np.nan,
                    "gse270497_auc": auc_270.get(arm, np.nan),
                }
            )
        if arm in c197.index:
            r.update(
                {
                    "gse197020_delta": c197.loc[arm, "delta_median_log2_case_minus_control"],
                    "gse197020_q": c197.loc[arm, "mannwhitney_q"] if "mannwhitney_q" in c197.columns else np.nan,
                    "gse197020_auc": auc_197.get(arm, np.nan),
                }
            )
        d270 = r.get("gse270497_delta")
        d197 = r.get("gse197020_delta")
        if pd.notna(d270) and pd.notna(d197):
            r["ev_direction_concordant"] = bool(np.sign(d270) == np.sign(d197))
        else:
            r["ev_direction_concordant"] = np.nan
        rows.append(r)

    rep = pd.DataFrame(rows)
    if not anticorr.empty:
        rep = rep.merge(anticorr, on="arm", how="left")
    return rep.sort_values("gse270497_delta", ascending=False, na_position="last")


def plot_replication_bars(rep: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    sub = rep.dropna(subset=["gse270497_delta", "gse197020_delta"], how="all").copy()
    if sub.empty:
        return

    sub = sub.set_index("arm")
    arms = sub.index.tolist()
    x = np.arange(len(arms))
    w = 0.35
    fig, ax = plt.subplots(figsize=(max(8, len(arms) * 0.55), 4.5))
    m270 = sub["gse270497_delta"].fillna(0).values
    m197 = sub["gse197020_delta"].fillna(0).values
    ax.bar(x - w / 2, m270, w, label="GSE270497 screening+ vs −", color="#C44E52", alpha=0.85)
    ax.bar(x + w / 2, m197, w, label="GSE197020 cancer vs normal", color="#4C72B0", alpha=0.85)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(arms, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Δ median log2 abundance (case − control)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("EV miRNA replication — focus arms")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(*, out_dir: Path = C.PLASMA_EV_DIR) -> Tuple[pd.DataFrame, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    # GSE197020
    expr197 = load_gse197020_matrix()
    samp197 = load_gse197020_samples()
    contrasts197 = case_control_arm_contrasts(
        expr197,
        samp197,
        label_col="bc_case",
        case_label="cancer",
        control_label="normal",
    )
    contrasts197.to_csv(out_dir / "gse197020_case_arm_contrasts.tsv", sep="\t", index=False)

    tcga_med = tcga_tissue_median_by_arm()
    focus197 = contrasts197.loc[contrasts197["arm"].isin(PANEL_MIRNA_TIER_ARM_IDS)].copy()
    compare197 = compare_ev_delta_to_tcga_tissue(
        focus197,
        tcga_med,
        delta_col="delta_median_log2_case_minus_control",
    )
    compare197.to_csv(out_dir / "gse197020_tcga_tissue_rank_compare.tsv", sep="\t", index=False)
    if not compare197.empty and "spearman_rho" in compare197.attrs:
        print(
            f"[ev_mirna_replication] GSE197020 panel EV Δ vs TCGA tissue ρ="
            f"{compare197.attrs['spearman_rho']:.3f} (p={compare197.attrs['spearman_p']:.3g})"
        )

    from mirna_hallmark.analyses.dcis_ev.ev_mirna_deep import run as run_deep

    run_deep(out_dir=out_dir / "deep", cohort="gse197020")

    # GSE270497 (reload if contrasts file missing)
    path270 = out_dir / "gse270497_screening_arm_contrasts.tsv"
    if path270.is_file():
        contrasts270 = pd.read_csv(path270, sep="\t")
        expr270 = load_ev_log2_matrix()
        samp270 = load_ev_samples()
    else:
        from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import run as run270

        contrasts270 = run270(out_dir=out_dir)
        expr270 = load_ev_log2_matrix()
        samp270 = load_ev_samples()

    auc270 = _arm_auc(expr270, samp270, "screening_positive")
    auc197 = _arm_auc(expr197, samp197, "bc_case")

    arms = list(REPLICATION_FOCUS_ARMS) + [a for a in PANEL_MIRNA_TIER_ARM_IDS if a not in REPLICATION_FOCUS_ARMS]
    arms = list(dict.fromkeys(a for a in arms if a in set(contrasts270["arm"]) | set(contrasts197["arm"])))

    anticorr = _anticorr_summary_for_arms(arms)
    replication = build_replication_table(contrasts270, contrasts197, auc270, auc197, anticorr, arms)
    replication.to_csv(out_dir / "ev_cohort_replication_focus_arms.tsv", sep="\t", index=False)

    mir29 = replication[replication["arm"].isin(MIR29_ARMS)].copy()
    mir29.to_csv(out_dir / "mir29_ev_tissue_bridge.tsv", sep="\t", index=False)

    plot_replication_bars(
        replication[replication["arm"].isin(REPLICATION_FOCUS_ARMS)],
        out_dir / "figures" / "ev_cohort_replication_focus_arms.png",
    )

    n_conc = int(replication["ev_direction_concordant"].fillna(False).sum())
    n_both = int(replication[["gse270497_delta", "gse197020_delta"]].notna().all(axis=1).sum())
    print(
        f"[ev_mirna_replication] GSE197020: {len(contrasts197)} arms; "
        f"replication table n={len(replication)} ({n_both} with both cohorts, {n_conc} same direction)"
    )
    m29 = mir29[["arm", "gse270497_delta", "gse270497_auc", "n_ecm_anticorr_hits", "n_anticorr_gene_x_pam50"]]
    if not m29.empty:
        print("[ev_mirna_replication] miR-29 bridge:\n" + m29.to_string(index=False))

    manifest = {
        "module": "mirna_hallmark.analyses.dcis_ev.ev_mirna_replication",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohorts": {
            "GSE270497": {"design": "screening_positive vs screening_negative", "n": int(len(samp270))},
            "GSE197020": {"design": "breast_cancer vs normal", "n": int(len(samp197))},
        },
        "focus_arms": list(REPLICATION_FOCUS_ARMS),
        "mir29_arms": list(MIR29_ARMS),
        "n_direction_concordant": n_conc,
        "n_with_both_cohorts": n_both,
        "note": "log2(TPM+1) vs log2(CPM+1) — comparable rank/AUC, not absolute delta calibration.",
        "gse197020_metadata": {
            "per_sample": ["disease_status", "geo_source_label", "ml_split", "ml_label"],
            "not_in_geo": ["stage", "grade", "ER", "HER2", "age", "treatment"],
        },
        "tcga_linkage": [
            "gse197020_tcga_tissue_rank_compare.tsv",
            "deep/gse197020_tcga_linkage_all_arms.tsv",
            "deep/gse197020_cohort_arm_correlations.tsv",
        ],
    }
    (out_dir / "replication_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return replication, mir29


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.PLASMA_EV_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
