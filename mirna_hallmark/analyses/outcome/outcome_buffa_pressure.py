"""Roadmap Phase 3 outcome — Buffa subset (GSE22216 miRNA + GSE22219 mRNA, 207) PRESSURE -> DRFS.

The Buffa METABRIC subset is the ONE independent cohort (besides TCGA) with **matched miRNA + mRNA +
clinical** — so it's the only place we can compute the framework's **pressure** construct and test
**pressure -> outcome independently of TCGA**. Computes Buffa-native gene/program pressure (same engine,
`compute_gene_pressure`), gated by Buffa's own AGO capacity, and Cox vs DRFS (overall + per-ER).
Validates MH-61: does pressure (esp. the HYPOXIA program — miR-210's program) predict relapse here too?

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.outcome.outcome_buffa_pressure``
"""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.eval.buffa_validation import load_buffa, BUFFA_MIRNA_SERIES
from mirna_hallmark.pressure_build import compute_gene_pressure
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.ago_gate import compute_capacity, gate_from_capacity
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z

OUT_DIR = C.OUTPUT_ROOT / "outcome_buffa_pressure"
TARGET_GENES = ["MKNK2", "KCNJ8"]
KEY = ["HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"]


def _clinical(samples) -> pd.DataFrame:
    lines = [l.rstrip("\n") for l in gzip.open(BUFFA_MIRNA_SERIES, "rt")]
    gsm = [x.strip('"') for x in [l for l in lines if l.lower().startswith("!sample_geo_accession")][0].split("\t")[1:]]

    def char(key):
        row = [l for l in lines if l.lower().startswith("!sample_characteristics") and key in l.lower()]
        if not row:
            return pd.Series(np.nan, index=gsm)
        vals = [x.strip('"').split(":", 1)[1].strip() if ":" in x else np.nan for x in row[0].split("\t")[1:]]
        return pd.Series(pd.to_numeric(vals, errors="coerce"), index=gsm)

    df = pd.DataFrame({"T": char("distant-relapse free survival"), "E": char("distant-relapse event"),
                       "er": char("er status")})
    return df.reindex(samples)


def _buffa_pressure():
    mi, rna = load_buffa()
    hs = HallmarkSets.load()
    gp = compute_gene_pressure(hs.universe, mirna=mi)                       # gene x sample (ungated)
    # Buffa 2011 Illumina array uses old aliases: EIF2C1-4 = AGO1-4, TNRC6 = TNRC6A
    rna = rna.rename(index={"EIF2C1": "AGO1", "EIF2C2": "AGO2", "EIF2C3": "AGO3", "EIF2C4": "AGO4", "TNRC6": "TNRC6A"})
    ago = [g for g in (list(C.AGO_GENES) + list(C.RISC_EXTRA_GENES)) if g in rna.index]
    sub = rna.loc[ago]                                                      # gene x sample (AGO genes as index)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=0).replace(0, np.nan), axis=0)
    gate = gate_from_capacity(compute_capacity(z)["ago_capacity_z"])
    gpg = gp.mul(gate.reindex(gp.columns), axis=1)                         # AGO-gated gene pressure
    prog = hallmark_pressure_matrix(gpg, hs)                               # program x sample
    return gpg, prog, mi


def _screen(mat, fam, clin):
    rows = []
    for feat in mat.index:
        x = mat.loc[feat].astype(float)
        for grp, ids in {"all": clin.index, "ER+": clin.index[clin["er"] == 1], "ER-": clin.index[clin["er"] == 0]}.items():
            c, p = _cox_one(clin.loc[ids, "T"], clin.loc[ids, "E"], x.reindex(ids), pd.DataFrame(index=ids))
            rows.append({"family": fam, "feature": feat, "stratum": grp, "hr": np.exp(c), "p": p,
                         "n_events": int(clin.loc[ids, "E"].sum())})
    r = pd.DataFrame(rows)
    if r.empty:
        return pd.DataFrame(columns=["family", "feature", "stratum", "hr", "p", "n_events", "q"])
    for st, idx in r.groupby("stratum").groups.items():
        ok = r.loc[idx, "p"].notna()
        r.loc[r.loc[idx].index[ok], "q"] = bh_fdr(r.loc[idx, "p"].dropna().to_numpy())
    return r


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    gpg, prog, mi = _buffa_pressure()
    clin = _clinical(prog.columns)
    print(f"[buffa-press] {prog.shape[1]} samples | DRFS events {int(clin['E'].sum())} "
          f"(ER+ {int(clin.loc[clin.er==1,'E'].sum())} / ER- {int(clin.loc[clin.er==0,'E'].sum())})")
    res = pd.concat([_screen(prog, "program_pressure", clin),
                     _screen(gpg.loc[gpg.index.intersection(TARGET_GENES)], "gene_pressure", clin),
                     _screen(mi.loc[mi.index.intersection(["hsa-miR-210-3p", "hsa-miR-210"])], "mir210_expr", clin)],
                    ignore_index=True)
    res.to_csv(out_dir / "buffa_pressure_drfs.tsv", sep="\t", index=False)

    prog_all = res[(res.family == "program_pressure") & (res.stratum == "all")]
    summary = {
        "module": "mirna_hallmark.analyses.outcome.outcome_buffa_pressure",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_samples": int(prog.shape[1]), "drfs_events": int(clin["E"].sum()),
        "program_pressure_DRFS_FDR_all": prog_all[prog_all.q < 0.1].sort_values("q")
            .apply(lambda r: f"{r['feature']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "HYPOXIA_pressure_DRFS": res[(res.feature == "HALLMARK_HYPOXIA")]
            .apply(lambda r: f"{r['stratum']} HR={r['hr']:.2f} p={r['p']:.3f}", axis=1).tolist(),
        "gene_pressure_DRFS": res[res.family == "gene_pressure"]
            .apply(lambda r: f"{r['feature']}/{r['stratum']} HR={r['hr']:.2f} p={r['p']:.3f}", axis=1).tolist(),
        "mir210_expr_DRFS": res[res.family == "mir210_expr"]
            .apply(lambda r: f"{r['feature']}/{r['stratum']} HR={r['hr']:.2f} p={r['p']:.3f}", axis=1).tolist(),
        "caveats": ["Buffa-native pressure (same engine + Buffa AGO gate); n=207, ~79 DRFS events",
                    "DRFS only (no OS in GSE22216); ER strata (no PAM50); independent of TCGA — the pressure-outcome cross-check"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[buffa-press] program-pressure DRFS FDR: {summary['program_pressure_DRFS_FDR_all']}")
    print(f"[buffa-press] HYPOXIA pressure DRFS: {summary['HYPOXIA_pressure_DRFS']}")
    print(f"[buffa-press] miR-210 expr DRFS: {summary['mir210_expr_DRFS']}")
    print(f"[buffa-press] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
