"""Within-sample miRNA<->mRNA coupling corroboration in the DCIS/IBC cohort.

The strongest layer GSE59247+GSE59246 offer: the two series are companion measurements
on the **same lesions** (paired by ``rep_group``), so for ~40 samples we have a sample's
miRNA *and* its mRNA. That lets us test the MH-48 finding **within-sample (within-patient)**
rather than by cross-platform direction-agreement:

  (1) COUPLING (within-sample): does each curated miRNA->target edge anti-correlate
      (miRNA up -> target mRNA down) across the paired DCIS/IBC samples? Negative ρ in this
      *independent* cohort corroborates the TCGA edge; we also report sign-concordance with
      the TCGA partial-ρ (``coupling_permutation/coupling_edge``) and the replication rate
      among TCGA-negative edges. Focused readout on the MH-48 invasion loss-leaders
      (miR-145 / miR-126 / miR-29c).
  (2) DE-REPRESSION (population): do the loss-leaders' target genes rise (de-repress)
      DCIS->IBC in the mRNA, as continued miRNA loss at invasion predicts?

This is association in a small independent cohort (~40 paired samples, Agilent log-signal,
rank-based), not proof of binding -- but it is genuine within-patient corroboration, which
the population timing analysis (MH-48, dcis_timing) could not provide.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.dcis_mrna_coupling
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, spearmanr

from analysis.utils.common.loaders import partial_spearman
from mirna_hallmark import config as C
from mirna_hallmark.analyses.dcis_ev.dcis_geo_loader import load_gse59246, load_gse59247
from mirna_hallmark.analyses.dcis_ev.dcis_timing import (
    _MARKER_SET, acquired_direction, bridge_v16_to_modern, composition_covariates,
)
from mirna_hallmark.stats import bh_fdr

OUT_DIR = C.OUTPUT_ROOT / "dcis_mrna_coupling"
EDGES = C.OUTPUT_ROOT / "edges" / "mirna_hallmark_edges.tsv.gz"
TCGA_EDGE = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
DE_LANDSCAPE = (C.TISSUE_REFERENCE_DIR / "cross_state_landscape"
                / "mirna_arm_de_landscape.tsv")
LOSS_LEADERS = ("miR-145", "miR-126", "miR-29c", "miR-140", "miR-497", "miR-143")
MIN_N = 15


def _pair() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Paired (miRNA arm x sample, mRNA gene x sample, meta) on common rep_group."""
    mirna, mmeta = load_gse59247()
    mrna, _ = load_gse59246()
    common = mirna.columns.intersection(mrna.columns)
    return mirna[common], mrna[common], mmeta.loc[common]


def _spearman(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < MIN_N:
        return np.nan, np.nan
    return spearmanr(x[m], y[m])


def within_sample_coupling(mirna, mrna, meta) -> pd.DataFrame:
    """Per high-evidence edge: within-sample Spearman(arm, target) across paired samples."""
    edges = pd.read_csv(EDGES, sep="\t")
    he = edges[edges["high_evidence"]][["miRNA", "gene"]].drop_duplicates()
    de = pd.read_csv(DE_LANDSCAPE, sep="\t")
    bridge = bridge_v16_to_modern(mirna.index, de["arm"], dict(zip(de["arm"], de["tumor_median"])))
    m2g = dict(zip(bridge["modern_arm"], bridge["gse_arm"]))      # modern -> GSE v16 arm
    acq = acquired_direction(de, None).set_index("modern_arm")["acq_dir"]

    cov = composition_covariates(mrna, meta)   # epi/immune/stroma/prolif + site batch
    he = he[he["miRNA"].isin(m2g) & he["gene"].isin(mrna.index)].copy()
    rows = []
    for arm, gene in he.itertuples(index=False):
        g = m2g[arm]
        rho, p = _spearman(mirna.loc[g].to_numpy(float), mrna.loc[gene].to_numpy(float))
        if not np.isfinite(rho):
            continue
        # composition+prolif+batch-ADJUSTED partial Spearman (the headline)
        rho_adj, p_adj, _ = partial_spearman(mrna.loc[gene], mirna.loc[g], cov)
        rows.append({"modern_arm": arm, "gse_arm": g, "gene": gene,
                     "rho_paired": rho, "p_paired": p,
                     "rho_paired_adj": rho_adj, "p_paired_adj": p_adj,
                     "acq_dir": int(acq.get(arm, 0)),
                     "target_is_marker": gene in _MARKER_SET})
    res = pd.DataFrame(rows)
    res["q_paired"] = bh_fdr(res["p_paired"].to_numpy())
    adj_ok = res["p_paired_adj"].notna()
    res["q_paired_adj"] = np.nan
    if adj_ok.any():
        res.loc[adj_ok, "q_paired_adj"] = bh_fdr(res.loc[adj_ok, "p_paired_adj"].to_numpy())
    # sign concordance with TCGA partial-ρ
    if TCGA_EDGE.exists():
        tc = pd.read_csv(TCGA_EDGE, sep="\t").rename(columns={"Unnamed: 0": "key"})
        tc["modern_arm"] = tc["key"].str.split(r"\|\|").str[0]
        tc["gene"] = tc["key"].str.split(r"\|\|").str[1]
        res = res.merge(tc[["modern_arm", "gene", "rho"]].rename(columns={"rho": "rho_tcga"}),
                        on=["modern_arm", "gene"], how="left")
        res["sign_concordant"] = np.sign(res["rho_paired"]) == np.sign(res["rho_tcga"])
        res["sign_concordant_adj"] = np.sign(res["rho_paired_adj"]) == np.sign(res["rho_tcga"])
    res["loss_leader"] = res["modern_arm"].str.contains("|".join(LOSS_LEADERS))
    return res.sort_values("rho_paired").reset_index(drop=True)


def target_derepression(mirna, mrna, meta) -> pd.DataFrame:
    """Do loss-leader target genes rise (de-repress) DCIS->IBC in the mRNA?"""
    edges = pd.read_csv(EDGES, sep="\t")
    he = edges[edges["high_evidence"]]
    tgt = he[he["miRNA"].str.contains("|".join(LOSS_LEADERS))]["gene"].unique()
    tgt = [g for g in tgt if g in mrna.index]
    dcis = meta.index[meta["state"] == "DCIS"]; ibc = meta.index[meta["state"] == "IBC"]
    rows = []
    for g in tgt:
        a = mrna.loc[g, ibc].to_numpy(float); b = mrna.loc[g, dcis].to_numpy(float)
        a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
        if len(a) < 4 or len(b) < 4:
            continue
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append({"gene": g, "delta_ibc_minus_dcis": float(np.median(a) - np.median(b)),
                     "rank_biserial": 2 * u / (len(a) * len(b)) - 1, "mwu_p": p})
    res = pd.DataFrame(rows)
    if not res.empty:
        res["mwu_q"] = bh_fdr(res["mwu_p"].to_numpy())
    return res.sort_values("delta_ibc_minus_dcis", ascending=False).reset_index(drop=True)


def _block(df: pd.DataFrame, rhocol: str, qcol: str) -> dict:
    s = df[rhocol].dropna()
    nfdr = int(((df[rhocol] < 0) & (df[qcol] < C.FDR_ALPHA)).sum())
    return {"n": int(len(s)), "frac_negative": round(float((s < 0).mean()), 3),
            "n_FDR_negative": nfdr, "median_rho": round(float(s.median()), 4)}


def _replication(cpl: pd.DataFrame, rhocol: str, concol: str) -> dict:
    both = cpl.dropna(subset=["rho_tcga", rhocol])
    tn = both[both["rho_tcga"] < 0]
    k, n = int((tn[rhocol] < 0).sum()), int(len(tn))
    return {"sign_concordance": round(float(both[concol].mean()), 3),
            "tcga_neg_replication": round(k / n, 3) if n else None,
            "binom_p_vs_50pct": float(f"{binomtest(k, n, 0.5, alternative='greater').pvalue:.3g}") if n else None,
            "n_tcga_negative_tested": n}


def _ll_enrich(cpl: pd.DataFrame, rhocol: str) -> dict:
    ll, bg = cpl[cpl["loss_leader"]], cpl[~cpl["loss_leader"]]
    a, b = ll[rhocol].dropna(), bg[rhocol].dropna()
    p = mannwhitneyu(a, b, alternative="less")[1] if len(a) and len(b) else np.nan
    return {"n_edges": int(len(a)), "frac_negative": round(float((a < 0).mean()), 3) if len(a) else None,
            "median_rho": round(float(a.median()), 4) if len(a) else None,
            "more_negative_than_background_mwu_p": float(f"{p:.3g}") if np.isfinite(p) else None}


def _summary(cpl: pd.DataFrame, der: pd.DataFrame, n_paired: int, meta) -> dict:
    has_tcga = "rho_tcga" in cpl
    ll = cpl[cpl["loss_leader"]]
    out = {
        "n_paired_samples": n_paired,
        "n_dcis": int((meta["state"] == "DCIS").sum()),
        "n_ibc": int((meta["state"] == "IBC").sum()),
        "n_edges_tested": int(len(cpl)),
        "covariates": "epi+immune+stroma metagenes + proliferation + collection-site batch",
        # ADJUSTED is the headline (stromal/site-confound removed); raw kept alongside
        "ADJUSTED": {
            "all_edges": _block(cpl, "rho_paired_adj", "q_paired_adj"),
            "replication_vs_tcga": _replication(cpl, "rho_paired_adj", "sign_concordant_adj") if has_tcga else None,
            "loss_leaders": _ll_enrich(cpl, "rho_paired_adj"),
        },
        "raw": {
            "all_edges": _block(cpl, "rho_paired", "q_paired"),
            "replication_vs_tcga": _replication(cpl, "rho_paired", "sign_concordant") if has_tcga else None,
            "loss_leaders": _ll_enrich(cpl, "rho_paired"),
        },
        "loss_leader_top_anticorrelated_adj": [
            {"arm": r.modern_arm, "gene": r.gene,
             "rho_adj": round(r.rho_paired_adj, 3) if np.isfinite(r.rho_paired_adj) else None,
             "rho_raw": round(r.rho_paired, 3), "q_adj": round(r.q_paired_adj, 4) if np.isfinite(r.q_paired_adj) else None,
             "target_is_marker": bool(r.target_is_marker)}
            for r in ll.dropna(subset=["rho_paired_adj"]).nsmallest(10, "rho_paired_adj").itertuples()],
        "n_loss_leader_targets_that_are_composition_markers": int(ll["target_is_marker"].sum()),
    }
    if not der.empty:
        up = der["delta_ibc_minus_dcis"] > 0
        out["derepression"] = {
            "n_loss_leader_targets": int(len(der)),
            "frac_up_at_invasion": round(float(up.mean()), 3),
            "n_FDR_up": int((up & (der["mwu_q"] < C.FDR_ALPHA)).sum()),
            "top_up": [{"gene": r.gene, "delta": round(r.delta_ibc_minus_dcis, 3),
                        "q": round(r.mwu_q, 4)} for r in der.head(6).itertuples()],
        }
    out["caveats"] = [
        "ADJUSTED (composition+prolif+site) is the headline; raw is shown for comparison",
        "association in a small independent cohort (~40 paired samples), rank-based; not binding",
        "within-sample = within-patient coupling (each sample contributes its own miRNA+mRNA)",
        "target_is_marker edges: composition adjustment may OVER-correct (target is itself a fraction marker)",
        "Agilent log-signal (miRNA GPL15019, mRNA GPL13607); v16->modern arm bridge; combined DCIS+IBC for power",
    ]
    return out


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[dcis_mrna_coupling] pairing GSE59247 miRNA <-> GSE59246 mRNA ...")
    mirna, mrna, meta = _pair()
    print(f"[dcis_mrna_coupling] {mirna.shape[1]} paired samples "
          f"(DCIS {int((meta['state']=='DCIS').sum())}, IBC {int((meta['state']=='IBC').sum())})")

    cpl = within_sample_coupling(mirna, mrna, meta)
    cpl.to_csv(out_dir / "within_sample_edge_coupling.tsv", sep="\t", index=False)
    der = target_derepression(mirna, mrna, meta)
    if not der.empty:
        der.to_csv(out_dir / "loss_leader_target_derepression.tsv", sep="\t", index=False)

    summary = _summary(cpl, der, mirna.shape[1], meta)
    manifest = {"module": "mirna_hallmark.analyses.dcis_ev.dcis_mrna_coupling",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "cohort": "GSE59247 miRNA + GSE59246 mRNA (paired, DCIS+IBC)", **summary}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    adj, raw = summary["ADJUSTED"], summary["raw"]
    print(f"[dcis_mrna_coupling] {summary['n_edges_tested']} edges (covariates: composition+prolif+site)")
    for tag, blk in (("ADJUSTED", adj), ("raw", raw)):
        rep = blk["replication_vs_tcga"]
        print(f"  [{tag}] all: {blk['all_edges']['frac_negative']:.0%} neg, "
              f"{blk['all_edges']['n_FDR_negative']} FDR-neg, median ρ={blk['all_edges']['median_rho']}"
              + (f"; TCGA-neg replication {rep['tcga_neg_replication']:.0%} (p={rep['binom_p_vs_50pct']})" if rep else "")
              + f"; loss-leaders {blk['loss_leaders']['frac_negative']:.0%} neg, "
              f"enrich p={blk['loss_leaders']['more_negative_than_background_mwu_p']}")
    print(f"  loss-leader targets that are composition markers: "
          f"{summary['n_loss_leader_targets_that_are_composition_markers']}")
    if "derepression" in summary:
        print(f"  loss-leader targets up at invasion: {summary['derepression']['frac_up_at_invasion']:.0%}")
    print(f"[dcis_mrna_coupling] wrote outputs under {out_dir}")
    return cpl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
