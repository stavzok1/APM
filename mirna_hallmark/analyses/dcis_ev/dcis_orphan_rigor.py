"""Rigor retrofit of the TS/ENCORI ORPHAN discovery list (MH-36/38) — composition test.

The orphan edges (TargetScan∪ENCORI, miRTarBase-functional-absent) were validated at protein
(CPTAC) + TCGA + Buffa with purity/CIN adjustment, but **never composition-adjusted** (epi/immune/
stroma). MH-49/51/54 showed that exact gap turns "miR-29→collagen" into a stromal-composition
artifact — and the orphan flagship is literally **miR-29b-3p→COL6A3/FSTL1** (collagen/ECM). So:
re-test every orphan edge in the **within-patient DCIS cohort** (GSE59247 miRNA + GSE59246 mRNA,
40 paired) **raw vs composition+prolif+batch-adjusted** partial Spearman, and ask which orphan
nominations SURVIVE the adjustment vs which were composition-driven.

This applies the apm-rigor-protocol (B confounders, C cross-cohort, F two-sided) to the orphan list,
converting orphan *discovery* into orphan discovery *with compartment attribution*. NB the DCIS
cohort is small (n=40, bulk) → a screening retrofit; the definitive version re-runs the orphan
coupling on the TCGA spine with composition metagenes.

Run: ``.venv/bin/python3 -m mirna_hallmark.dcis_orphan_rigor``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr

from analysis.utils.common.loaders import partial_spearman
from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.dcis_mrna_coupling import _pair
from mirna_hallmark.analyses.dcis_ev.dcis_timing import bridge_v16_to_modern, composition_covariates
from mirna_hallmark.stats import bh_fdr

ORPHANS = C.OUTPUT_ROOT / "buffa_validation" / "orphan_triple_cohort_validation.tsv"
DE_LANDSCAPE = (C.TISSUE_REFERENCE_DIR / "cross_state_landscape" / "mirna_arm_de_landscape.tsv")
OUT_DIR = C.OUTPUT_ROOT / "dcis_orphan_rigor"
MIN_N = 15
ECM_TOK = ("COL", "FBN", "SPARC", "LOX", "FN1", "ELN", "LAMA", "LAMB", "FSTL", "POSTN", "VCAN", "MMP")


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    orph = pd.read_csv(ORPHANS, sep="\t")[["miRNA", "gene", "triple_validated", "triple_validated_sig"]].drop_duplicates()
    mirna, mrna, meta = _pair()
    de = pd.read_csv(DE_LANDSCAPE, sep="\t")
    bridge = bridge_v16_to_modern(mirna.index, de["arm"], dict(zip(de["arm"], de["tumor_median"])))
    m2g = dict(zip(bridge["modern_arm"], bridge["gse_arm"]))
    cov = composition_covariates(mrna, meta)
    print(f"[orphan_rigor] {len(orph)} orphan edges; {mirna.shape[1]} paired DCIS samples; "
          f"bridgeable arms {orph['miRNA'].isin(m2g).sum()}/{orph['miRNA'].nunique()}")

    rows = []
    for r in orph.itertuples(index=False):
        g = m2g.get(r.miRNA)
        if g is None or r.gene not in mrna.index:
            continue
        x = mirna.loc[g]; y = mrna.loc[r.gene]
        m = ~(x.isna() | y.isna())
        if m.sum() < MIN_N:
            continue
        rho, p = spearmanr(x[m], y[m])
        rho_adj, p_adj, _ = partial_spearman(y, x, cov)
        rows.append({"miRNA": r.miRNA, "gene": r.gene, "triple_validated": bool(r.triple_validated),
                     "is_ECM": any(t in r.gene.upper() for t in ECM_TOK),
                     "rho_raw": rho, "p_raw": p, "rho_adj": rho_adj, "p_adj": p_adj})
    res = pd.DataFrame(rows)
    res["q_raw"] = bh_fdr(res["p_raw"].fillna(1).to_numpy())
    ok = res["p_adj"].notna()
    res["q_adj"] = np.nan
    res.loc[ok, "q_adj"] = bh_fdr(res.loc[ok, "p_adj"].to_numpy())
    res.to_csv(out_dir / "orphan_composition_retrofit.tsv", sep="\t", index=False)

    def _frac_neg(d, col):
        s = d[col].dropna()
        return round(float((s < 0).mean()), 3) if len(s) else None
    neg_raw = res["rho_raw"] < 0
    surv = res[(res["rho_raw"] < 0) & (res["rho_adj"] < 0)]   # stayed negative after adjustment
    killed = res[(res["rho_raw"] < 0) & (res["rho_adj"] >= 0)]
    ecm = res[res["is_ECM"]]
    # is the median |rho| attenuated by adjustment among raw-negatives?
    raw_med = float(res.loc[neg_raw, "rho_raw"].median())
    adj_med = float(res.loc[neg_raw, "rho_adj"].median())
    summary = {
        "module": "mirna_hallmark.dcis_orphan_rigor",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_orphan_edges_tested": int(len(res)),
        "frac_negative_raw": _frac_neg(res, "rho_raw"),
        "frac_negative_adj": _frac_neg(res, "rho_adj"),
        "median_rho_among_raw_neg": {"raw": round(raw_med, 3), "adjusted": round(adj_med, 3),
                                     "attenuation": round(raw_med - adj_med, 3)},
        "raw_neg_kept_neg_after_adj": int(len(surv)),
        "raw_neg_killed_by_adj": int(len(killed)),
        "ECM_orphans": {
            "n": int(len(ecm)),
            "frac_neg_raw": _frac_neg(ecm, "rho_raw"),
            "frac_neg_adj": _frac_neg(ecm, "rho_adj"),
            "median_rho_raw": round(float(ecm["rho_raw"].median()), 3) if len(ecm) else None,
            "median_rho_adj": round(float(ecm["rho_adj"].median()), 3) if len(ecm) else None,
        },
        "top_survivors_adj": [
            {"edge": f"{x.miRNA}->{x.gene}", "rho_raw": round(x.rho_raw, 3),
             "rho_adj": round(x.rho_adj, 3), "q_adj": round(x.q_adj, 3) if np.isfinite(x.q_adj) else None,
             "ECM": bool(x.is_ECM)}
            for x in surv.dropna(subset=["rho_adj"]).nsmallest(10, "rho_adj").itertuples()],
        "caveats": ["DCIS within-patient n=40, bulk → screening retrofit; definitive = TCGA spine + composition metagenes",
                    "orphan list bridged v16→modern; coverage limited; orphans were curated-absent by design",
                    "tests whether composition adjustment changes orphan survival (the MH-49/54 lesson applied to orphans)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[orphan_rigor] tested {len(res)} edges; neg raw {summary['frac_negative_raw']:.0%} → adj {summary['frac_negative_adj']:.0%}; "
          f"median ρ(raw-neg) {raw_med:.3f}→{adj_med:.3f}")
    print(f"[orphan_rigor] raw-neg kept {len(surv)} / killed {len(killed)} by adjustment")
    print(f"[orphan_rigor] ECM orphans (n={len(ecm)}): neg raw {summary['ECM_orphans']['frac_neg_raw']} → adj {summary['ECM_orphans']['frac_neg_adj']}; "
          f"median ρ {summary['ECM_orphans']['median_rho_raw']}→{summary['ECM_orphans']['median_rho_adj']}")
    print(f"[orphan_rigor] top adj survivors: {[s['edge'] for s in summary['top_survivors_adj'][:6]]}")
    print(f"[orphan_rigor] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
