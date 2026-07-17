"""EV (exosomal) miRNA across the isogenic MCF10A breast-progression series (GSE297448).

Why this cohort is special: MCF10A -> MCF10AT1 -> MCF10DCIS -> MCF10CA1a is an
**isogenic** progression model (Normal -> benign proliferative -> DCIS -> invasive
carcinoma), so the ordered trajectory `Normal->DCIS->invasive` that the cross-sectional
patient cohorts could NOT license (MH-48: no temporal order across patients) IS
legitimate here -- there is no patient heterogeneity or field-effect confound, only the
genetic step-series, and the readout is the **exosomal (EV) cargo**.

This lets us ask the EV-direction question MH-48 raised: the invasion loss-leaders
(miR-145 / miR-29c / miR-126 / miR-140) are *lost from cells* across progression -- are
they **EXPORTED into EVs** (EV abundance rises while the cellular pool falls = selective
shedding) or simply down everywhere? We cross-reference each miRNA's **EV progression
trend** against its **cellular tumour-vs-healthy direction** (TCGA GTEx-anchored
`gtex_to_tumor_pct`): cellular-DOWN + EV-UP = selective export.

n = 3 EV replicates x 4 stages = 12; readout is the monotone stage trend (Spearman vs
stage ordinal) and the DCIS->invasive step. Small isogenic cell-line cohort -> a
mechanistic-direction probe, not a population estimate.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.dcis_ev_progression
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.dcis_timing import _jonckheere, acquired_direction
from mirna_hallmark.stats import bh_fdr

# the source paper's headline exported tumour-suppressor miRNA (Rab27A-driven export);
# we cross-check our compositional selectivity against it.
PAPER_EXPORTED = ("miR-205",)

EV_MATRIX = (C.REPO_ROOT / "data" / "external" / "dcis_geo"
             / "GSE297448_Normalized_microRNA_readcounts_All_samples.txt.gz")
DE_LANDSCAPE = (C.TISSUE_REFERENCE_DIR / "cross_state_landscape"
                / "mirna_arm_de_landscape.tsv")
OUT_DIR = C.OUTPUT_ROOT / "plasma_ev" / "mcf10a_progression"
LOSS_LEADERS = ("miR-145", "miR-126", "miR-29c", "miR-140", "miR-497", "miR-143")
# isogenic step series -> ordinal stage
_STAGE = {"MCF10A": 0, "MCF10AT1": 1, "MCF10DCIS": 2, "MCF10CA1a": 3}
_STAGE_NAME = {0: "Normal", 1: "benign", 2: "DCIS", 3: "invasive"}


def load_ev() -> Tuple[pd.DataFrame, pd.Series]:
    """Return (ev: miRNA x sample log2 readcount, stage: sample -> ordinal 0..3)."""
    df = pd.read_csv(EV_MATRIX, sep="\t", index_col=0)
    # exact cell-line token (NB "MCF10AT1".startswith("MCF10A") -> must match the token, not a prefix)
    stage = pd.Series({c: _STAGE[c.split()[0]] for c in df.columns})
    ev = np.log2(df.astype(float) + 1.0)
    return ev, stage


def _clr(ev: pd.DataFrame) -> pd.DataFrame:
    """Centered-log-ratio: subtract each SAMPLE's mean log -> compositional EV *share*.

    Required here because the raw 'normalized' readcounts still carry a strong global
    EV-cargo increase along progression (~95% of all miRNA rise), so absolute 'EV-up'
    is the baseline and tells us nothing about *selective* export. CLR removes the
    per-sample global loading so the trend measures relative EV-share change.
    """
    return ev.sub(ev.mean(axis=0), axis=1)


def progression_trend(ev: pd.DataFrame, stage: pd.Series,
                      raw: pd.DataFrame = None) -> pd.DataFrame:
    """Per miRNA: Spearman of EV *share* (CLR) vs stage ordinal (12 pts) + DCIS->inv step.

    ``ev`` is the CLR (compositional) matrix; ``raw`` (log2) is used only for reporting
    absolute EV abundance so a high-share trend isn't read off a near-absent miRNA.
    """
    s = stage.reindex(ev.columns).to_numpy(float)
    cols = {k: ev.columns[stage == k] for k in (0, 1, 2, 3)}   # Normal/benign/DCIS/invasive
    raw = ev if raw is None else raw
    rows = []
    for arm in ev.index:
        y = ev.loc[arm].to_numpy(float)
        if np.unique(y).size < 2:
            continue
        rho, p = spearmanr(y, s)
        means = {k: float(ev.loc[arm, c].mean()) for k, c in cols.items()}
        steps = {"normal_to_benign": means[1] - means[0],
                 "benign_to_dcis": means[2] - means[1],
                 "dcis_to_invasive": means[3] - means[2]}
        jz, jp = _jonckheere(y, s)   # ordered-alternatives test across the 4 stages
        signs = [np.sign(v) for v in steps.values() if abs(v) > 1e-9]
        non_monotone = len(set(signs)) > 1   # at least one step reverses direction
        rows.append({"arm": arm, "ev_trend_rho": rho, "ev_trend_p": p,
                     "ev_jt_z": jz, "ev_jt_p": jp, "non_monotone": non_monotone,
                     "share_normal": round(means[0], 3), "share_benign": round(means[1], 3),
                     "share_dcis": round(means[2], 3), "share_invasive": round(means[3], 3),
                     "step_normal_to_benign": round(steps["normal_to_benign"], 3),
                     "step_benign_to_dcis": round(steps["benign_to_dcis"], 3),
                     "ev_share_dcis_to_invasive": round(steps["dcis_to_invasive"], 3),
                     "ev_abundance_log2": float(raw.loc[arm].mean())})
    res = pd.DataFrame(rows)
    res["ev_trend_q"] = bh_fdr(res["ev_trend_p"].to_numpy())
    return res


def export_direction(trend: pd.DataFrame) -> pd.DataFrame:
    """Join cellular tumour-vs-healthy direction (TCGA) -> classify export vs co-move."""
    de = pd.read_csv(DE_LANDSCAPE, sep="\t")
    acq = acquired_direction(de, None).rename(columns={"modern_arm": "arm"})[
        ["arm", "acq_dir", "gtex_to_tumor_pct", "trajectory"]]
    t = trend.merge(acq, on="arm", how="left")
    ev_up = t["ev_trend_rho"] > 0

    def _cls(r):
        if r.get("acq_dir", 0) == 0 or not np.isfinite(r["ev_trend_rho"]):
            return "unclassified"
        cell_down, cell_up = r["acq_dir"] < 0, r["acq_dir"] > 0
        up = r["ev_trend_rho"] > 0
        if cell_down and up:
            return "EXPORTED (cell-down, EV-up)"
        if cell_up and up:
            return "co-accumulated (cell-up, EV-up)"
        if cell_down and not up:
            return "depleted (cell-down, EV-down)"
        return "retained (cell-up, EV-down)"

    t["ev_class"] = t.apply(_cls, axis=1)
    t["loss_leader"] = t["arm"].str.contains("|".join(LOSS_LEADERS))
    return t.sort_values("ev_trend_rho", ascending=False).reset_index(drop=True)


def _summary(t: pd.DataFrame, raw_global_trend: float) -> dict:
    cls = t[t["ev_class"] != "unclassified"]
    losers = cls[cls["acq_dir"] < 0]
    gainers = cls[cls["acq_dir"] > 0]
    ll = t[t["loss_leader"] & t["ev_trend_rho"].notna()]
    conc = (np.sign(cls["ev_trend_rho"]) == np.sign(cls["acq_dir"])).mean()
    # the REAL export-selectivity test (on compositional CLR trend, global loading removed):
    # are cellular LOSERS preferentially gaining EV share vs cellular GAINERS?
    p_export = np.nan
    if len(losers) and len(gainers):
        _, p_export = mannwhitneyu(losers["ev_trend_rho"], gainers["ev_trend_rho"],
                                   alternative="greater")
    return {
        "n_ev_mirna": int(len(t)),
        "raw_global_median_trend_BEFORE_clr": round(float(raw_global_trend), 3),
        "note_on_normalization": ("raw EV trend is globally positive (global cargo increase) "
                                  "-> analysis uses CLR EV-SHARE; trends below are compositional"),
        "n_classified": int(len(cls)),
        "ev_share_class_counts": cls["ev_class"].value_counts().to_dict(),
        "ev_share_vs_cellular_direction_concordance": round(float(conc), 3),
        "EXPORT_selectivity_losers_gt_gainers_mwu_p": (
            float(f"{p_export:.3g}") if np.isfinite(p_export) else None),
        # FDR over the FULL cellular-loser set: how many losers are confirmed EV-share exporters
        "losers_FDR_confirmed_exported": int(
            ((losers["ev_trend_rho"] > 0) & (losers["ev_trend_q"] < C.FDR_ALPHA)).sum()),
        "cellular_LOSERS_gaining_EV_share": {
            "n": int(len(losers)),
            "frac_share_up": round(float((losers["ev_trend_rho"] > 0).mean()), 3)},
        "cellular_GAINERS_gaining_EV_share": {
            "n": int(len(gainers)),
            "frac_share_up": round(float((gainers["ev_trend_rho"] > 0).mean()), 3)},
        "n_non_monotone": int(t["non_monotone"].sum()) if "non_monotone" in t else None,
        "paper_exported_check_miR205": [
            {"arm": r.arm, "ev_share_trend_rho": round(r.ev_trend_rho, 3),
             "ev_trend_q": round(r.ev_trend_q, 4), "ev_class": r.ev_class,
             "steps": [round(r.step_normal_to_benign, 2), round(r.step_benign_to_dcis, 2),
                       round(r.ev_share_dcis_to_invasive, 2)]}
            for r in t[t["arm"].str.contains("|".join(PAPER_EXPORTED))].itertuples()],
        "loss_leaders": [
            {"arm": r.arm, "ev_share_trend_rho": round(r.ev_trend_rho, 3),
             "ev_trend_q": round(r.ev_trend_q, 4), "ev_jt_p": round(r.ev_jt_p, 4),
             "non_monotone": bool(r.non_monotone),
             "steps_N_b_D_i": [round(r.step_normal_to_benign, 2), round(r.step_benign_to_dcis, 2),
                               round(r.ev_share_dcis_to_invasive, 2)],
             "ev_class": r.ev_class}
            for r in ll.sort_values("ev_trend_rho", ascending=False).itertuples()],
        "top_EV_share_gain": [
            {"arm": r.arm, "rho": round(r.ev_trend_rho, 3), "q": round(r.ev_trend_q, 4)}
            for r in t.nlargest(8, "ev_trend_rho").itertuples()],
        "caveats": [
            "Rab27A-driven selective export is the source paper's mechanism (miR-205 headline); we test selectivity compositionally",
            "isogenic MCF10A cell-line EVs (n=3/stage) -> mechanistic direction, not population",
            "ordered trajectory legitimate ONLY because the series is isogenic (no patient/field confound)",
            "raw readcounts carry a global EV-cargo increase -> analysis is compositional (CLR share); absolute EV-up is uninformative",
            "EXPORT inferred from EV-SHARE gain vs the INDEPENDENT cellular (TCGA) direction; the cell-line's own cellular pool is not measured here",
        ],
    }


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[dcis_ev] loading MCF10A isogenic EV progression (GSE297448) ...")
    ev, stage = load_ev()
    print(f"[dcis_ev] {ev.shape[0]} EV miRNA x {ev.shape[1]} samples; stages "
          f"{ {(_STAGE_NAME[s]): int((stage==s).sum()) for s in sorted(stage.unique())} }")
    # raw global trend (documents the global-cargo confound that forces the CLR analysis)
    s = stage.reindex(ev.columns).to_numpy(float)
    raw_global = float(np.median([spearmanr(ev.loc[a].to_numpy(float), s)[0]
                                  for a in ev.index if np.unique(ev.loc[a]).size > 1]))
    trend = progression_trend(_clr(ev), stage, raw=ev)   # compositional EV-share trend
    t = export_direction(trend)
    t.to_csv(out_dir / "ev_progression_trend.tsv", sep="\t", index=False)

    summary = _summary(t, raw_global)
    manifest = {"module": "mirna_hallmark.analyses.dcis_ev.dcis_ev_progression",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "cohort": "GSE297448 (MCF10A isogenic EV: Normal->benign->DCIS->invasive)",
                "stages": _STAGE_NAME, **summary}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[dcis_ev] raw global trend BEFORE CLR={summary['raw_global_median_trend_BEFORE_clr']} "
          f"(global cargo increase -> using compositional EV-share)")
    print(f"[dcis_ev] EXPORT-selectivity (losers gain EV-share > gainers) MWU p="
          f"{summary['EXPORT_selectivity_losers_gt_gainers_mwu_p']}; "
          f"losers FDR-confirmed exported={summary['losers_FDR_confirmed_exported']}/"
          f"{summary['cellular_LOSERS_gaining_EV_share']['n']}; "
          f"non-monotone={summary['n_non_monotone']}")
    print(f"[dcis_ev] share-class counts {summary['ev_share_class_counts']}")
    if summary["paper_exported_check_miR205"]:
        print(f"[dcis_ev] miR-205 (paper headline): {summary['paper_exported_check_miR205']}")
    print("[dcis_ev] loss-leaders (compositional EV-share; steps N→b→D→i):")
    for r in summary["loss_leaders"]:
        print(f"   {r['arm']:18s} share-trend ρ={r['ev_share_trend_rho']:+.2f} q={r['ev_trend_q']:.3f} "
              f"steps={r['steps_N_b_D_i']} {'NONMONO ' if r['non_monotone'] else ''}{r['ev_class']}")
    print(f"[dcis_ev] wrote outputs under {out_dir}")
    return t


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
