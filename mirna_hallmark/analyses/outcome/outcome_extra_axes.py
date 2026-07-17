"""Roadmap Phase 3 outcome — extra (non-pressure) axes vs TCGA survival, for completeness:
  #4 miRNA-CNV dosage (per-locus copy number + amplification/loss burden) — the genomic axis.
  #7 pressure RESIDUAL burden (mean |expression unexplained by miRNA-arm pressure|).
Cox vs TCGA DFI+OS, composition-adjusted, BH-FDR. (#8 acquired-miRNA reduces to a per-arm constant
shift of arm expression = already null MH-60; #5 NAT-referenced needs the NAT-pressure build — deferred.)

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_extra_axes``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, load_tcga_outcome

OUT_DIR = C.OUTPUT_ROOT / "outcome_extra_axes"
CNV = C.OUTPUT_ROOT / "mirna_locus_cnv" / "tables" / "sample_entity_cnv.tsv.gz"
RESID = C.OUTPUT_ROOT / "gene_expression_explainability" / "cohort" / "per_sample_residuals.tsv.gz"


def _screen(mat: pd.DataFrame, fam: str, cl: dict) -> pd.DataFrame:
    adj = pd.concat([cl["base"], cl["comp"]], axis=1)
    rows = []
    for f in mat.index:
        x = mat.loc[f].astype(float)
        if x.std(skipna=True) == 0 or x.notna().sum() < 50:
            continue
        dca, dpa = _cox_one(cl["dfi_t"], cl["dfi_e"], x, adj)
        oca, opa = _cox_one(cl["os_t"], cl["os_e"], x, adj)
        rows.append({"axis": fam, "feature": str(f), "dfi_hr": np.exp(dca), "dfi_p": dpa,
                     "os_hr": np.exp(oca), "os_p": opa})
    r = pd.DataFrame(rows)
    if r.empty:
        return r
    for col in ("dfi_p", "os_p"):
        ok = r[col].notna(); r.loc[ok, col.replace("_p", "_q")] = bh_fdr(r.loc[ok, col].to_numpy())
    return r


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cl = load_tcga_outcome()
    # #4 miRNA-CNV
    c = pd.read_csv(CNV, sep="\t", usecols=["participant", "entity_id", "copy_number", "has_gain", "has_loss"], low_memory=False)
    c["participant"] = c["participant"].astype(str)
    burden = c.groupby("participant").agg(cnv_mean=("copy_number", "mean"),
                                          frac_gain=("has_gain", "mean"), frac_loss=("has_loss", "mean")).T
    loci = c.pivot_table(index="entity_id", columns="participant", values="copy_number", aggfunc="mean")
    loci = loci[loci.notna().sum(axis=1) >= 200]                        # loci present in >=200 patients
    print(f"[extra] miRNA-CNV: {loci.shape[0]} loci + 3 burden features")
    # #7 residual burden
    d = pd.read_csv(RESID, sep="\t", index_col=0)
    rcol = "residual_core_arms" if "residual_core_arms" in d.columns else "residual_base_meth"
    resid = pd.DataFrame({"resid_unexplained_burden": d.groupby(level=0)[rcol].apply(lambda s: s.abs().mean())}).T

    res = pd.concat([_screen(burden, "cnv_burden", cl), _screen(loci, "cnv_locus", cl),
                     _screen(resid, "pressure_residual", cl)], ignore_index=True)
    res = res[~res.empty] if isinstance(res, pd.Series) else res
    res.to_csv(out_dir / "extra_axes_survival.tsv", sep="\t", index=False)

    fdr = res[(res.dfi_q < 0.1) | (res.os_q < 0.1)]
    summary = {
        "module": "mirna_hallmark.outcome_extra_axes",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "by_axis": {a: {"n": int(len(g)), "dfi_FDR": int((g.dfi_q < 0.1).sum()), "os_FDR": int((g.os_q < 0.1).sum())}
                    for a, g in res.groupby("axis")},
        "burden_results": res[res.axis == "cnv_burden"].apply(lambda r: f"{r['feature']} DFI HR={r['dfi_hr']:.2f} q={r['dfi_q']:.3f} | OS q={r['os_q']:.3f}", axis=1).tolist(),
        "residual_result": res[res.axis == "pressure_residual"].apply(lambda r: f"DFI HR={r['dfi_hr']:.2f} q={r['dfi_q']:.3f} | OS HR={r['os_hr']:.2f} q={r['os_q']:.3f}", axis=1).tolist(),
        "FDR_hits": fdr.sort_values(["dfi_q", "os_q"]).head(15).apply(lambda r: f"{r['axis']}:{r['feature']} DFI q={r['dfi_q']:.3f} OS q={r['os_q']:.3f}", axis=1).tolist(),
        "caveats": ["TCGA-only, composition-adjusted; #8 acquired-miRNA ~ arm-expression (null); #5 NAT-referenced deferred"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[extra] by axis: {summary['by_axis']}")
    print(f"[extra] CNV burden: {summary['burden_results']}")
    print(f"[extra] residual: {summary['residual_result']}")
    print(f"[extra] FDR hits: {summary['FDR_hits'][:8]}")
    print(f"[extra] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
