"""Independent GSE255660 serum miRNA screen: top AUC contrasts × TCGA pressure spines.

Contrasts (excl. dilution replicates):
  - cancer vs healthy
  - cancer vs benign
  - benign vs healthy

Ranks arms by discriminative AUC per contrast, joins:
  1. miRTarBase combined-pressure anticorr (``target_combined_anticorr``)
  2. TargetScan top predicted targets + cohort partial anticorr (novel vs curated)

Outputs under ``mirna_hallmark/output/plasma_ev/gse255660_screen/``:
  ``gse255660_all_contrast_auc.tsv``
  ``gse255660_top_auc_pressure_bridge.tsv``
  ``gse255660_targetscan_coherent_stories.tsv``
  ``figures/gse255660_top_auc_by_contrast.png``

Run:
  .venv/bin/python3 -m mirna_hallmark.ev_mirna_gse255660_screen
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from mirna_hallmark import config as C
from mirna_hallmark import stats as S
from mirna_hallmark.analyses.dcis_ev.ev_mirna_multi_cohort import _read_log2_wide
from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import case_control_arm_contrasts
from mirna_hallmark.analyses.dcis_ev.ev_mirna_ts_pressure import (
    TS_MIN_WEIGHT,
    build_ev_pressure_bridge,
    coherent_story_table,
)
from pipeline.config import PATHS

OUT_DIR = C.PLASMA_EV_DIR / "gse255660_screen"
TOP_N_PER_CONTRAST = 30
MIN_AUC = 0.58


def _prep_contrast(samples: pd.DataFrame, mode: str) -> pd.DataFrame:
    s = samples.loc[~samples["is_dilution_replicate"].astype(bool)].copy()
    if mode == "cancer_vs_healthy":
        s = s[s["bc_case"].astype(bool) | s["healthy_control"].astype(bool)]
        s["case_label"] = s["bc_case"].astype(bool)
    elif mode == "cancer_vs_benign":
        s = s[s["bc_case"].astype(bool) | s["bc_benign"].astype(bool)]
        s["case_label"] = s["bc_case"].astype(bool)
    elif mode == "benign_vs_healthy":
        s = s[s["bc_benign"].astype(bool) | s["healthy_control"].astype(bool)]
        s["case_label"] = s["bc_benign"].astype(bool)
    else:
        raise ValueError(mode)
    return s


def _discriminative_auc(auc: float) -> float:
    """AUC oriented to max separation (flip if <0.5)."""
    if pd.isna(auc):
        return np.nan
    return float(auc) if auc >= 0.5 else 1.0 - float(auc)


def run_all_contrast_aucs() -> pd.DataFrame:
    expr = _read_log2_wide(PATHS.mirna_ev_gse255660_log2_expr)
    samples = pd.read_csv(PATHS.mirna_ev_gse255660_samples, sep="\t")

    rows: List[dict] = []
    for mode, label in (
        ("cancer_vs_healthy", "cancer vs healthy"),
        ("cancer_vs_benign", "cancer vs breast_benign"),
        ("benign_vs_healthy", "benign vs healthy"),
    ):
        sub = _prep_contrast(samples, mode)
        contrasts = case_control_arm_contrasts(
            expr,
            sub,
            label_col="case_label",
            case_label="case",
            control_label="control",
            min_total_samples=10,
        )
        for _, r in contrasts.iterrows():
            arm = r["arm"]
            sample_ids = [c for c in expr.columns if c in set(sub["sample_id"])]
            labels = sub.set_index("sample_id").loc[sample_ids, "case_label"].astype(bool)
            x = pd.to_numeric(expr.loc[arm, sample_ids], errors="coerce")
            ok = x.notna()
            auc = np.nan
            if ok.sum() >= 10 and labels[ok].nunique() == 2:
                auc = float(roc_auc_score(labels[ok].astype(int), x[ok]))
            rows.append(
                {
                    "arm": arm,
                    "contrast": mode,
                    "contrast_label": label,
                    "n_case": int(r.get("n_case", labels.sum())),
                    "n_control": int(r.get("n_control", (~labels).sum())),
                    "delta_case_minus_control": r["delta_median_log2_case_minus_control"],
                    "mannwhitney_q": r.get("mannwhitney_q", np.nan),
                    "auc_case_higher": auc,
                    "discriminative_auc": _discriminative_auc(auc),
                    "direction": "case_higher" if r["delta_median_log2_case_minus_control"] > 0 else "control_higher",
                }
            )
    out = pd.DataFrame(rows)
    m = out["mannwhitney_q"].notna()
    out.loc[m, "mannwhitney_q"] = S.bh_fdr(out.loc[m, "mannwhitney_q"].values)
    out["sig_fdr05"] = out["mannwhitney_q"] < C.FDR_ALPHA
    return out.sort_values(["contrast", "discriminative_auc"], ascending=[True, False])


def _top_arms_per_contrast(auc_table: pd.DataFrame) -> pd.DataFrame:
    parts = []
    for contrast, sub in auc_table.groupby("contrast"):
        sub = sub.sort_values("discriminative_auc", ascending=False)
        sub = sub.loc[sub["discriminative_auc"] >= MIN_AUC].head(TOP_N_PER_CONTRAST)
        parts.append(sub)
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def plot_top_auc(top: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)
    for ax, (contrast, label) in zip(
        axes,
        [
            ("cancer_vs_healthy", "cancer vs healthy"),
            ("cancer_vs_benign", "cancer vs benign"),
            ("benign_vs_healthy", "benign vs healthy"),
        ],
    ):
        sub = top.loc[top["contrast"] == contrast].head(15)
        if sub.empty:
            ax.set_visible(False)
            continue
        y = np.arange(len(sub))
        colors = ["#C44E52" if d == "case_higher" else "#4C72B0" for d in sub["direction"]]
        ax.barh(y, sub["discriminative_auc"], color=colors, alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(sub["arm"], fontsize=7)
        ax.set_xlim(0.5, 1.0)
        ax.set_title(label, fontsize=9)
        ax.invert_yaxis()
    fig.suptitle("GSE255660 top serum miRNA discriminative AUC by contrast", fontsize=11)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(*, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = C.PLASMA_EV_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    auc_all = run_all_contrast_aucs()
    auc_all.to_csv(out_dir / "gse255660_all_contrast_auc.tsv", sep="\t", index=False)

    top = _top_arms_per_contrast(auc_all)
    top.to_csv(out_dir / "gse255660_top_auc_per_contrast.tsv", sep="\t", index=False)

    arms = sorted(top["arm"].unique())
    bridge, ts_detail = build_ev_pressure_bridge(auc_all, arms)
    if not ts_detail.empty:
        ts_detail.to_csv(out_dir / "gse255660_targetscan_edge_coupling.tsv", sep="\t", index=False)

    bridge, _ = build_ev_pressure_bridge(auc_all, arms)
    bridge.to_csv(out_dir / "gse255660_top_auc_pressure_bridge.tsv", sep="\t", index=False)

    stories = coherent_story_table(bridge, ts_detail, auc_all)
    stories.to_csv(out_dir / "gse255660_targetscan_coherent_stories.tsv", sep="\t", index=False)

    plot_top_auc(top, fig_dir / "gse255660_top_auc_by_contrast.png")

    manifest = {
        "module": "mirna_hallmark.ev_mirna_gse255660_screen",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "contrasts": ["cancer_vs_healthy", "cancer_vs_benign", "benign_vs_healthy"],
        "top_n_per_contrast": TOP_N_PER_CONTRAST,
        "min_auc": MIN_AUC,
        "ts_min_weight": TS_MIN_WEIGHT,
        "n_arms_tested": int(auc_all["arm"].nunique()),
        "n_top_union": len(arms),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[gse255660_screen] tested {manifest['n_arms_tested']} arms; top union n={len(arms)}")
    show = bridge[
        [
            c
            for c in [
                "arm",
                "best_discriminative_auc",
                "discriminative_auc_cancer_vs_healthy",
                "discriminative_auc_cancer_vs_benign",
                "discriminative_auc_benign_vs_healthy",
                "mirtar_n_anticorr_pam50",
                "n_ts_neg_sig_cohort",
                "n_ts_orphan_neg_sig",
                "top_ts_orphan_neg_sig",
                "pressure_coherence",
                "discovery_score",
            ]
            if c in bridge.columns
        ]
    ].head(20)
    print("\n[gse255660_screen] Top discovery bridge (AUC × pressure):\n")
    print(show.to_string(index=False))

    print("\n[gse255660_screen] TS orphan coherent stories:\n")
    print(stories.loc[stories["ts_orphan_story"]].head(12).to_string(index=False))

    return {"auc_all": auc_all, "bridge": bridge, "stories": stories}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
