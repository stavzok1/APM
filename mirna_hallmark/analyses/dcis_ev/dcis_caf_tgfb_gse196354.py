"""Direct TGFβ-axis test in matched breast CAF vs NF mRNA (GSE196354).

MH-54 showed miR-29c is NOT down in CAFs → the desmoplastic ECM is miR-29-independent; the single-cell
Risom test showed FAP⁺/αSMA⁺ myCAFs carry the collagen (TGFβ-myCAF *proxies*). This closes the proxy
gap at the mRNA level: in **7 paired primary breast NF vs CAF** RNA-seq, are the **TGFβ pathway** and
the **ECM/collagen** program (the miR-29 targets) co-up-regulated in CAF — i.e. is the desmoplastic
collagen a TGFβ program rather than a miR-29-loss program?

Rigor moves (apm-rigor-protocol): A paired design (within-patient NF↔CAF); D direct measurement;
I positive control (ACTA2/FAP must rise) + null control (housekeeping must stay flat); G power-floor
(paired n=7 → min Wilcoxon p≈0.016); J known-vs-novel (TGFβ→CAF desmoplasia is textbook — this is
confirmation of the alternative driver, not a new mechanism).

DATA: ``data/external/caf_nf_mrna_gse196354/GSE196354_NFs_CAFs_RNAseq_FPKM.txt.gz`` (GEO).
Run: ``.venv/bin/python3 -m mirna_hallmark.dcis_caf_tgfb_gse196354``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr

FPKM = (C.REPO_ROOT / "data" / "external" / "caf_nf_mrna_gse196354"
        / "GSE196354_NFs_CAFs_RNAseq_FPKM.txt.gz")
OUT_DIR = C.OUTPUT_ROOT / "dcis_caf_tgfb_gse196354"
GENE_SETS = {
    "TGFb_pathway": ["TGFB1", "TGFB2", "TGFB3", "TGFBI", "SERPINE1", "CTGF", "CCN2", "CYR61",
                     "CCN1", "THBS1", "SKIL", "JUNB", "SMAD7", "LTBP1", "LTBP2", "TGFBR1", "TGFBR2"],
    "myCAF_activation": ["ACTA2", "FAP", "POSTN", "TAGLN", "PDPN", "MMP11"],
    "ECM_collagen_miR29targets": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL4A1",
                                  "SPARC", "FN1", "LOX", "LOXL2", "FBN1"],
    "positive_control": ["ACTA2", "FAP"],          # must be CAF-up if the cohort is real
    "null_control": ["ACTB", "GAPDH", "TBP", "B2M"],  # housekeeping → should be flat
}


def load_paired():
    """(genes x sample log2 FPKM, median-centred, NF cols, CAF cols) paired by patient P1..P7.

    NB the null-control check exposed a global CAF>NF FPKM shift, so we **median-centre each sample
    over expressed genes** (remove the global/library shift) before testing direction — otherwise
    everything (incl. housekeeping) drifts up. This is the compositional fix (apm-rigor-protocol B/F).
    """
    d = pd.read_csv(FPKM, sep="\t", index_col=0)
    d = np.log2(d.astype(float) + 1.0)
    expressed = d.index[d.mean(axis=1) > 1.0]            # genes expressed across the cohort
    d = d.sub(d.loc[expressed].mean(axis=0), axis=1)     # per-sample centring on expressed genes
    pats = sorted({c.split("_")[0] for c in d.columns})
    nf = [f"{p}_NFs" for p in pats if f"{p}_NFs" in d.columns]
    caf = [f"{p}_CAFs" for p in pats if f"{p}_CAFs" in d.columns]
    return d, nf, caf


def paired_genes(d, nf, caf, genes):
    rows = []
    for g in genes:
        if g not in d.index:
            continue
        a = d.loc[g, nf].to_numpy(float)
        b = d.loc[g, caf].to_numpy(float)
        if np.allclose(a, b):
            continue
        try:
            p = wilcoxon(b, a)[1]
        except ValueError:
            p = np.nan
        rows.append({"gene": g, "nf_mean": round(float(a.mean()), 2),
                     "caf_mean": round(float(b.mean()), 2),
                     "log2fc_caf_vs_nf": round(float(np.median(b - a)), 2),
                     "wilcoxon_p": p})
    return pd.DataFrame(rows)


def run(*, out_dir: Path = OUT_DIR):
    out_dir.mkdir(parents=True, exist_ok=True)
    d, nf, caf = load_paired()
    print(f"[caf_tgfb] {d.shape[0]} genes; {len(nf)} NF / {len(caf)} CAF (paired)")
    all_rows, summary = [], {}
    for setname, genes in GENE_SETS.items():
        r = paired_genes(d, nf, caf, genes)
        if r.empty:
            continue
        r["set"] = setname
        all_rows.append(r)
        up = r["log2fc_caf_vs_nf"] > 0
        summary[setname] = {
            "n_genes": int(len(r)), "frac_up_in_CAF": round(float(up.mean()), 2),
            "median_log2fc": round(float(r["log2fc_caf_vs_nf"].median()), 2),
            "n_sig_up_p<=0.05": int((up & (r["wilcoxon_p"] <= 0.05)).sum()),
            "genes": [{"g": x.gene, "log2fc": x.log2fc_caf_vs_nf, "p": round(float(x.wilcoxon_p), 3)}
                      for x in r.sort_values("log2fc_caf_vs_nf", ascending=False).itertuples()],
        }
    res = pd.concat(all_rows, ignore_index=True)
    res["q_bh"] = bh_fdr(res["wilcoxon_p"].fillna(1).to_numpy())
    res.to_csv(out_dir / "caf_vs_nf_tgfb_ecm.tsv", sep="\t", index=False)

    pos_ok = summary.get("positive_control", {}).get("frac_up_in_CAF", 0) >= 0.5
    null_flat = summary.get("null_control", {}).get("n_sig_up_p<=0.05", 9) == 0
    verdict = ("TGFβ pathway + ECM/collagen co-UP in CAF → desmoplastic ECM is a TGFβ program "
               "(miR-29-independent, consistent with MH-54)"
               if summary.get("TGFb_pathway", {}).get("frac_up_in_CAF", 0) >= 0.6
               and summary.get("ECM_collagen_miR29targets", {}).get("frac_up_in_CAF", 0) >= 0.6
               else "mixed — inspect per-gene")
    manifest = {"module": "mirna_hallmark.dcis_caf_tgfb_gse196354",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "cohort": "GSE196354 — 7 paired primary breast NF vs CAF RNA-seq (FPKM)",
                "by_set": summary,
                "controls_ok": {"positive(ACTA2/FAP)_up": bool(pos_ok), "null(housekeeping)_flat": bool(null_flat)},
                "verdict": verdict,
                "caveats": ["paired n=7 (Wilcoxon p-floor ~0.016); FPKM; primary CAF/NF cultures",
                            "TGFβ-pathway/ECM co-up is consistent-with but not proof of TGFβ causation",
                            "confirms the alternative driver to miR-29 (MH-54); textbook biology"]}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")

    for s in ("positive_control", "null_control", "TGFb_pathway", "myCAF_activation", "ECM_collagen_miR29targets"):
        if s in summary:
            v = summary[s]
            print(f"  [{s:26s}] up_in_CAF {v['frac_up_in_CAF']:.0%}  median log2FC {v['median_log2fc']:+.2f}  sig_up {v['n_sig_up_p<=0.05']}/{v['n_genes']}")
    print(f"  controls: positive_up={pos_ok}  null_flat={null_flat}")
    print(f"[caf_tgfb] verdict: {verdict}")
    print(f"[caf_tgfb] wrote outputs under {out_dir}")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
