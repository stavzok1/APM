"""GSE270497 EV screening cohort: discriminative AUC + TargetScan pressure bridge.

Design: screening-positive vs screening-negative plasma exosome miRNA (n≈180).

Outputs under ``mirna_hallmark/output/plasma_ev/gse270497_ts/``:
  ``gse270497_screening_auc.tsv``
  ``gse270497_top_auc_pressure_bridge.tsv``
  ``gse270497_targetscan_edge_coupling.tsv``
  ``gse270497_targetscan_coherent_stories.tsv``

Run:
  .venv/bin/python3 -m mirna_hallmark.ev_mirna_gse270497_ts
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from mirna_hallmark import config as C
from mirna_hallmark import stats as S
from mirna_hallmark.analyses.dcis_ev.ev_mirna_replication import REPLICATION_FOCUS_ARMS
from mirna_hallmark.analyses.dcis_ev.ev_mirna_screening import case_control_arm_contrasts, load_ev_log2_matrix, load_ev_samples
from mirna_hallmark.analyses.dcis_ev.ev_mirna_ts_pressure import (
    TS_MIN_WEIGHT,
    build_ev_pressure_bridge,
    coherent_story_table,
)

OUT_DIR = C.PLASMA_EV_DIR / "gse270497_ts"
TOP_N = 30
MIN_AUC = 0.58

# GSE255660 placental/oncofetal serum cluster + TS orphan follow-ups.
ONCOFETAL_ARMS: Sequence[str] = (
    "hsa-miR-642a-3p",
    "hsa-miR-208a-5p",
    "hsa-miR-187-5p",
    "hsa-miR-371a-5p",
    "hsa-miR-652-5p",
    "hsa-miR-135a-3p",
    "hsa-miR-92a-2-5p",
    "hsa-miR-370-3p",
    "hsa-miR-486-3p",
)
TS_ORPHAN_FOLLOWUP_ARMS: Sequence[str] = (
    "hsa-miR-493-3p",
    "hsa-miR-493-5p",
    "hsa-miR-301a-3p",
)


def _discriminative_auc(auc: float) -> float:
    if pd.isna(auc):
        return np.nan
    return float(auc) if auc >= 0.5 else 1.0 - float(auc)


def run_screening_auc() -> pd.DataFrame:
    expr = load_ev_log2_matrix()
    samples = load_ev_samples()
    contrasts = case_control_arm_contrasts(
        expr,
        samples,
        label_col="screening_positive",
        case_label="screening_positive",
        control_label="screening_negative",
    )
    sample_ids = [c for c in expr.columns if c in set(samples["sample_id"])]
    labels = samples.set_index("sample_id").loc[sample_ids, "screening_positive"].astype(bool)

    rows: List[dict] = []
    for _, r in contrasts.iterrows():
        arm = r["arm"]
        x = pd.to_numeric(expr.loc[arm, sample_ids], errors="coerce")
        ok = x.notna()
        auc = np.nan
        if ok.sum() >= 10 and labels[ok].nunique() == 2:
            auc = float(roc_auc_score(labels[ok].astype(int), x[ok]))
        rows.append(
            {
                "arm": arm,
                "contrast": "screening_pos_vs_neg",
                "n_case": int(r["n_screening_positive"]),
                "n_control": int(r["n_screening_negative"]),
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
    return out.sort_values("discriminative_auc", ascending=False)


def run(*, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)

    auc_all = run_screening_auc()
    auc_all.to_csv(out_dir / "gse270497_screening_auc.tsv", sep="\t", index=False)

    top = auc_all.loc[auc_all["discriminative_auc"] >= MIN_AUC].head(TOP_N)
    extra = [
        a
        for a in list(REPLICATION_FOCUS_ARMS) + list(ONCOFETAL_ARMS) + list(TS_ORPHAN_FOLLOWUP_ARMS)
        if a in set(auc_all["arm"])
    ]
    arms = sorted(set(top["arm"]) | set(extra))

    bridge, ts_detail = build_ev_pressure_bridge(auc_all, arms)
    bridge.to_csv(out_dir / "gse270497_top_auc_pressure_bridge.tsv", sep="\t", index=False)
    if not ts_detail.empty:
        ts_detail.to_csv(out_dir / "gse270497_targetscan_edge_coupling.tsv", sep="\t", index=False)

    stories = coherent_story_table(bridge, ts_detail, auc_all)
    stories.to_csv(out_dir / "gse270497_targetscan_coherent_stories.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.ev_mirna_gse270497_ts",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort": "GSE270497",
        "design": "screening_positive vs screening_negative",
        "top_n": TOP_N,
        "min_auc": MIN_AUC,
        "ts_min_weight": TS_MIN_WEIGHT,
        "n_arms_tested": int(auc_all["arm"].nunique()),
        "n_bridge_union": len(arms),
        "extra_arms": {"focus": list(REPLICATION_FOCUS_ARMS), "oncofetal": list(ONCOFETAL_ARMS)},
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[gse270497_ts] tested {manifest['n_arms_tested']} arms; bridge union n={len(arms)}")
    show_cols = [
        c
        for c in [
            "arm",
            "best_discriminative_auc",
            "discriminative_auc_screening_pos_vs_neg",
            "direction_screening_pos_vs_neg",
            "mirtar_n_anticorr_pam50",
            "n_ts_neg_sig_cohort",
            "n_ts_orphan_neg_sig",
            "top_ts_orphan_neg_sig",
            "pressure_coherence",
            "discovery_score",
        ]
        if c in bridge.columns
    ]
    print("\n[gse270497_ts] Top bridge:\n", bridge[show_cols].head(20).to_string(index=False))

    onc = bridge.loc[bridge["arm"].isin(ONCOFETAL_ARMS)]
    if not onc.empty:
        print("\n[gse270497_ts] Oncofetal cluster (screening EV × TCGA TS):\n")
        print(
            onc[
                [
                    c
                    for c in show_cols
                    if c in onc.columns
                ]
            ].to_string(index=False)
        )

    return {"auc_all": auc_all, "bridge": bridge, "stories": stories, "ts_detail": ts_detail}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
