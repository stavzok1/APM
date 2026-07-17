"""Deep miR-301a-3p target coupling: full TargetScan panel + PAM50 strata.

Separates evidence classes (do not bundle):
  - ``published_brca`` — luciferase / exosome literature (BTG1, PTEN, MEOX2, …)
  - ``mirtar_he`` — Hallmark high-evidence miRTar edge
  - ``mirtar_low`` — miRTar summary only (n_studies=1, score<=3)
  - ``ts_orphan`` — TargetScan only, not HE

Outputs under ``mirna_hallmark/output/tissue_reference/mir301a_depth/``:
  ``mir301a_full_ts_coupling_cohort.tsv`` — all TS targets (w>=0.25), cohort partials
  ``mir301a_coupling_by_pam50.tsv`` — same edges × (cohort + PAM50)
  ``mir301a_evidence_class_summary.tsv`` — per-class hit counts
  ``mir301a_published_vs_orphan_bridge.tsv`` — side-by-side cohort + subtype
  ``mir301a_target_combined_drivers.tsv`` — miRTar pressure spine rows

Run:
  .venv/bin/python3 -m mirna_hallmark.ev_mirna_301a_target_depth
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.dcis_ev.ev_mirna_ts_pressure import load_targetscan_all_for_arms
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.analyses.mir301.mir301_focus_genes import (
    MIR301A_ARM,
    PRIOR_ORPHAN_HITS,
    PUBLISHED_BRCA,
    evidence_class,
    he_pairs,
)
from mirna_hallmark.robustness_checks import (
    PAM50_SUBTYPES,
    _pam50_scope_iter,
    _partial_ladder,
    _proliferation_proxies,
)
from mirna_hallmark.eval.targetscan_orphan_coupling import _lit_spot_check

ARM = MIR301A_ARM
OUT_DIR = C.TISSUE_REFERENCE_DIR / "mir301a_depth"
# Luciferase / exosome-validated or strongly cited BRCA-relevant targets.
# (Defined in mir301_focus_genes; re-exported here for backward compatibility.)
from mirna_hallmark.analyses.mir301.mir301_focus_genes import PRIOR_ORPHAN_HITS, PUBLISHED_BRCA  # noqa: F401

TS_MIN_WEIGHT = 0.25
MIN_SCOPE_N = 20


def _he_pairs():
    return he_pairs()


def _evidence_class(gene: str, arm: str, he_pairs_set):
    return evidence_class(gene, arm, he_pairs_set)


def _coupling_rows(
    arm: str,
    genes: Sequence[str],
    *,
    scopes_only: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    clinical = D.load_clinical_strata()
    clin = clinical.set_index("participant")
    hs = HallmarkSets.load()
    proxies = _proliferation_proxies(rna, hs)
    cnv = D.load_cnv_target_genes(genes)
    he_pairs = _he_pairs()

    if arm not in mirna.index:
        return pd.DataFrame()

    x_all = mirna.loc[arm]
    rows: List[dict] = []
    for gene in genes:
        if gene not in rna.index:
            continue
        y_all = rna.loc[gene]
        cn_row = cnv.loc[gene] if gene in cnv.index else None
        for scope, samples in _pam50_scope_iter(clinical):
            if scopes_only is not None and scope not in scopes_only:
                continue
            cols = x_all.index.intersection(y_all.index).intersection(clin.index)
            if samples is not None:
                cols = cols.intersection(samples)
            cols = pd.Index(sorted(cols))
            if len(cols) < MIN_SCOPE_N:
                continue
            ladder = _partial_ladder(
                y_all,
                x_all,
                clin,
                cols,
                proxies,
                target_cn=cn_row,
            )
            rho_cn = ladder.get("rho_CPE_HRD_CN", ladder.get("rho_CPE_HRD"))
            p_cn = ladder.get("p_CPE_HRD_CN", ladder.get("p_CPE_HRD"))
            neg = bool(pd.notna(rho_cn) and rho_cn < 0 and pd.notna(p_cn) and p_cn < 0.05)
            rows.append(
                {
                    "arm": arm,
                    "gene": gene,
                    "scope": scope,
                    "n": len(cols),
                    "evidence_class": _evidence_class(gene, arm, he_pairs),
                    "mirtar_he": (arm, gene) in he_pairs,
                    "lit_spot_check": _lit_spot_check(gene, arm),
                    "neg_sig": neg,
                    **ladder,
                }
            )
    return pd.DataFrame(rows)


def _load_target_combined_for_arm(arm: str) -> pd.DataFrame:
    top_path = C.TARGET_COMBINED_ANTICORR_DIR / "target_combined_top_mirnas.tsv"
    by_pam_path = C.TARGET_COMBINED_ANTICORR_DIR / "target_combined_anticorr_by_pam50.tsv"
    if not top_path.is_file() or not by_pam_path.is_file():
        return pd.DataFrame()

    top = pd.read_csv(top_path, sep="\t")
    genes = set(top.loc[top["miRNA"] == arm, "gene"])
    by_pam = pd.read_csv(by_pam_path, sep="\t")
    by_pam = by_pam[(by_pam["view"] == "gated") & (by_pam["gene"].isin(genes))]
    by_pam = by_pam.assign(arm=arm)
    return by_pam.sort_values(["gene", "pam50"])


def _class_summary(cohort: pd.DataFrame, by_pam: pd.DataFrame) -> pd.DataFrame:
    rows: List[dict] = []
    for cls, sub in cohort.groupby("evidence_class"):
        neg = sub.loc[sub["neg_sig"]]
        rows.append(
            {
                "evidence_class": cls,
                "n_genes_tested": len(sub),
                "n_neg_sig_cohort": len(neg),
                "median_abs_rho_CN": round(float(neg["rho_CPE_HRD_CN"].abs().median()), 4)
                if "rho_CPE_HRD_CN" in neg.columns and len(neg)
                else np.nan,
                "top_neg_sig_genes": ";".join(
                    neg.reindex(neg["rho_CPE_HRD_CN"].abs().sort_values(ascending=False).index)["gene"]
                    .astype(str)
                    .head(8)
                    .tolist()
                )
                if len(neg)
                else "",
            }
        )
    # subtype counts for neg_sig
    if not by_pam.empty:
        for cls in cohort["evidence_class"].unique():
            genes = set(cohort.loc[cohort["evidence_class"] == cls, "gene"])
            sub = by_pam[(by_pam["gene"].isin(genes)) & (by_pam["partial_rho"] < 0) & (by_pam["partial_q"] < C.FDR_ALPHA)]
            for pam in list(PAM50_SUBTYPES) + ["cohort"]:
                if pam == "cohort":
                    continue
                n = int((sub["pam50"] == pam).sum()) if not sub.empty else 0
                rows.append({"evidence_class": cls, "pam50_neg_sig_gated": pam, "n": n})
    return pd.DataFrame(rows)


def _published_vs_orphan_bridge(cohort: pd.DataFrame, by_pam: pd.DataFrame) -> pd.DataFrame:
    focus = sorted(set(PUBLISHED_BRCA) | set(PRIOR_ORPHAN_HITS))
    c = cohort[cohort["gene"].isin(focus)].copy()
    if c.empty:
        return c

    wide_c = c.set_index("gene")
    rows: List[dict] = []
    for gene in focus:
        if gene not in wide_c.index:
            continue
        r = wide_c.loc[gene]
        row = {
            "gene": gene,
            "evidence_class": r["evidence_class"],
            "ts_weight": np.nan,
            "cohort_rho_CN": r.get("rho_CPE_HRD_CN"),
            "cohort_p_CN": r.get("p_CPE_HRD_CN"),
            "cohort_neg_sig": bool(r.get("neg_sig")),
            "lit_note": r.get("lit_spot_check", ""),
        }
        if not by_pam.empty:
            gs = by_pam[by_pam["gene"] == gene]
            for pam in PAM50_SUBTYPES:
                gp = gs[gs["pam50"] == pam]
                if len(gp):
                    row[f"{pam}_partial_rho"] = gp["partial_rho"].iloc[0]
                    row[f"{pam}_partial_q"] = gp["partial_q"].iloc[0]
                    row[f"{pam}_neg_sig"] = bool(
                        gp["partial_rho"].iloc[0] < 0 and gp["partial_q"].iloc[0] < C.FDR_ALPHA
                    )
        rows.append(row)

    out = pd.DataFrame(rows)
    # attach TS weight
    ts_path = C.PLASMA_EV_DIR / "gse255660_screen" / "hsa-miR-301a-3p_targetscan_targets.tsv"
    if ts_path.is_file():
        ts = pd.read_csv(ts_path, sep="\t")
        w = ts.set_index("gene")["ts_weight"]
        out["ts_weight"] = out["gene"].map(w)
    return out.sort_values(["evidence_class", "cohort_rho_CN"])


def run(*, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)

    rna = D.load_rna()
    ts_edges = load_targetscan_all_for_arms([ARM], rna.index, min_weight=TS_MIN_WEIGHT)
    ts_edges = ts_edges.drop_duplicates(["arm", "gene"])
    genes = sorted(set(ts_edges["gene"]) | set(PUBLISHED_BRCA) | set(PRIOR_ORPHAN_HITS))

    # merge ts weights
    ts_w = ts_edges.set_index("gene")["ts_weight"] if not ts_edges.empty else pd.Series(dtype=float)

    by_pam_all = _coupling_rows(ARM, genes)
    cohort = by_pam_all[by_pam_all["scope"] == "cohort"].copy()
    if not cohort.empty and len(ts_w):
        cohort["ts_weight"] = cohort["gene"].map(ts_w)
        cohort = cohort.sort_values("ts_weight", ascending=False, na_position="last")

    by_pam = by_pam_all[by_pam_all["scope"] != "cohort"].copy()
    if not by_pam.empty and len(ts_w):
        by_pam["ts_weight"] = by_pam["gene"].map(ts_w)

    tc = _load_target_combined_for_arm(ARM)

    cohort.to_csv(out_dir / "mir301a_full_ts_coupling_cohort.tsv", sep="\t", index=False)
    by_pam.to_csv(out_dir / "mir301a_coupling_by_pam50.tsv", sep="\t", index=False)
    tc.to_csv(out_dir / "mir301a_target_combined_drivers.tsv", sep="\t", index=False)

    bridge = _published_vs_orphan_bridge(cohort, tc)
    bridge.to_csv(out_dir / "mir301a_published_vs_orphan_bridge.tsv", sep="\t", index=False)

    summary_rows: List[dict] = []
    for cls in ("published_brca", "mirtar_he", "mirtar_low", "ts_orphan"):
        sub = cohort[cohort["evidence_class"] == cls] if not cohort.empty else pd.DataFrame()
        neg = sub.loc[sub["neg_sig"]] if not sub.empty else pd.DataFrame()
        row = {
            "evidence_class": cls,
            "n_genes_cohort": len(sub),
            "n_neg_sig_cohort": len(neg),
            "top_neg_sig": ";".join(neg["gene"].astype(str).head(8).tolist()) if len(neg) else "",
        }
        if not by_pam.empty:
            gset = set(sub["gene"]) if len(sub) else set()
            bp = by_pam[by_pam["gene"].isin(gset)]
            for pam in PAM50_SUBTYPES:
                ps = bp[(bp["scope"] == pam) & bp["neg_sig"]]
                row[f"n_neg_sig_{pam}"] = len(ps["gene"].unique()) if len(ps) else 0
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "mir301a_evidence_class_summary.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.ev_mirna_301a_target_depth",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "arm": ARM,
        "n_ts_genes": int(len(ts_edges)),
        "n_cohort_edges": int(len(cohort)),
        "published_brca": list(PUBLISHED_BRCA),
        "prior_orphan_hits": list(PRIOR_ORPHAN_HITS),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[mir301a_depth] cohort edges={len(cohort)}; by_pam rows={len(by_pam)}")
    print("\n[mir301a_depth] Evidence class summary:\n", summary.to_string(index=False))
    print("\n[mir301a_depth] Published vs orphan bridge:\n", bridge.to_string(index=False))

    # Highlight published targets that were missing from top-20 screen
    pub = cohort[cohort["evidence_class"] == "published_brca"][
        ["gene", "ts_weight", "rho_CPE_HRD_CN", "p_CPE_HRD_CN", "neg_sig"]
    ]
    if len(pub):
        print("\n[mir301a_depth] Published BRCA targets (full panel):\n", pub.to_string(index=False))

    return {
        "cohort": cohort,
        "by_pam": by_pam,
        "bridge": bridge,
        "summary": summary,
        "target_combined": tc,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
