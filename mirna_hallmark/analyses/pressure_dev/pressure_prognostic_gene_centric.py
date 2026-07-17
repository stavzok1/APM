"""Design build 1b — the GENE-centric mirror of `pressure_prognostic_signature` (which controlled NODES).

Per target gene (cancer-relevant, coupled), per patient, readouts:
  received_pressure = framework gene-pressure it receives (Σ miRNA repression)
  realized          = -z(expression)                         (its repression state)
  decoupling        = resid(z(expr) ~ z(received_pressure))  (escape: expressed despite pressure)
  acquired_NAT      = -(log2TPM - NAT mean)                  (repressed vs matched NAT)
  acquired_GTEx     = -(log2TPM - GTEx breast)               (vs true-healthy)
Gene pool = BY-neg-coupled genes that are cancer-relevant (curated role OR high DepMap dependency).
Per-axis LASSO-Cox 5-fold nested-CV concordance vs full clinical; TCGA (DFI+OS) + Buffa (DRFS).

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.pressure_dev.pressure_prognostic_gene_centric``
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.gene_roles import load_gene_roles, load_gene_dependency
from mirna_hallmark.analyses.outcome.outcome_survival import _z, load_tcga_outcome
from mirna_hallmark.analyses.pressure_dev.pressure_prognostic_signature import (
    _cv_concordance, _log2tpm_tcga, _nat_mean, _gtex_mean, _edges_w,
)

OUT_DIR = C.OUTPUT_ROOT / "pressure_prognostic_gene_centric"
L1_GP = C.OUTPUT_ROOT / "pressure_layer_comparison" / "L1_coupling_spine" / "gene_pressure_per_sample.tsv.gz"


def _zrows(df):
    return df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1).replace(0, np.nan), axis=0)


def _gene_pool():
    coupled = set(_edges_w()["gene"])
    gr = set(load_gene_roles()["gene"])
    dep = load_gene_dependency().set_index("gene")["dep_role_weight"]
    hidep = set(dep[dep.abs() > dep.abs().quantile(0.7)].index)
    return sorted(coupled & (gr | hidep))                              # coupled AND cancer-relevant


def _readouts(l2tpm, gp, pool, nat_mean, gtex_mean):
    genes = [g for g in pool if g in l2tpm.index and g in gp.index]
    E = l2tpm.loc[genes]; P = gp.loc[genes]
    common = E.columns.intersection(P.columns)
    E, P = E[common], P[common]
    Ez, Pz = _zrows(E), _zrows(P)
    out = {"received_pressure": P, "realized": -Ez}
    # decoupling: per gene resid(Ez ~ Pz)
    A, B = Pz.to_numpy(), Ez.to_numpy()
    Ar = A - np.nanmean(A, 1, keepdims=True); Br = B - np.nanmean(B, 1, keepdims=True)
    slope = np.nansum(Ar * Br, 1) / np.where(np.nansum(Ar ** 2, 1) == 0, np.nan, np.nansum(Ar ** 2, 1))
    out["decoupling"] = pd.DataFrame(B - (A * slope[:, None] + (np.nanmean(B, 1) - slope * np.nanmean(A, 1))[:, None]), index=genes, columns=common)
    if nat_mean is not None:
        out["acquired_NAT"] = -(E.sub(nat_mean.reindex(genes).values, axis=0))
    if gtex_mean is not None:
        out["acquired_GTEx"] = -(E.sub(gtex_mean.reindex(genes).fillna(E.mean(axis=1)).values, axis=0))
    return out


def _buffa_full_clin(bmi):
    from mirna_hallmark.eval.buffa_validation import BUFFA_MIRNA_SERIES
    lines = [l.rstrip("\n") for l in gzip.open(BUFFA_MIRNA_SERIES, "rt")]
    gsm = [x.strip('"') for x in [l for l in lines if l.lower().startswith("!sample_geo_accession")][0].split("\t")[1:]]
    def ch(k):
        rr = [l for l in lines if l.lower().startswith("!sample_characteristics") and k in l.lower()]
        return pd.Series(pd.to_numeric([x.strip('"').split(":", 1)[1].strip() if ":" in x else np.nan for x in rr[0].split("\t")[1:]], errors="coerce"), index=gsm) if rr else pd.Series(np.nan, index=gsm)
    full = pd.DataFrame({"er": ch("er status"), "age": _z(ch("patient age")), "grade": _z(ch("tumour grade")),
                         "nodes": _z(ch("nodes involved")), "size": _z(ch("tumour size"))})
    return ch("distant-relapse free survival"), ch("distant-relapse event"), full


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pool = _gene_pool()
    print(f"[gene-sig] gene pool: {len(pool)} cancer-relevant coupled genes")
    results = {}
    # TCGA
    l2 = _log2tpm_tcga(); gp = pd.read_csv(L1_GP, sep="\t", index_col=0); cl = load_tcga_outcome()
    tro = _readouts(l2, gp, pool, _nat_mean(), _gtex_mean())
    cov = pd.concat([cl["base"], cl["comp"]], axis=1)
    for axis, fm in tro.items():
        for ep in ("dfi", "os"):
            r = _cv_concordance(fm, cov, cl[f"{ep}_t"], cl[f"{ep}_e"])
            if r:
                results[f"TCGA/{axis}/{ep.upper()}"] = r
    # Buffa (received_pressure/realized/decoupling)
    from mirna_hallmark.eval.buffa_validation import load_buffa
    from mirna_hallmark.analyses.outcome.outcome_buffa_pressure import _buffa_pressure
    bmi, brna = load_buffa(); gpg, _, _ = _buffa_pressure()
    bT, bE, bfull = _buffa_full_clin(bmi)
    bro = _readouts(brna, gpg, pool, None, None)
    for axis, fm in bro.items():
        r = _cv_concordance(fm, bfull, bT, bE)
        if r:
            results[f"Buffa/{axis}/DRFS"] = r

    summary = {
        "module": "mirna_hallmark.analyses.pressure_dev.pressure_prognostic_gene_centric",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_genes": len(pool),
        "results": results,
        "beats_clinical": {k: v for k, v in results.items() if v.get("delta", -1) > 0.01},
        "caveats": ["gene-centric mirror; cancer-relevant coupled gene pool; LASSO nested-CV vs full clinical",
                    "TCGA sparse for outcome (overfit-prone); METABRIC(Buffa)/METABRIC-full is the arbiter"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print("[gene-sig] CV concordance (pressure+clin vs clin) by axis/cohort:")
    for k, v in results.items():
        print(f"   {k:30s} {v['cv_concordance_pressure+clin']} vs {v['cv_concordance_clin_only']} (delta {v['delta']:+.3f})  feats={v['stable_features'][:4]}")
    print(f"[gene-sig] beats clinical: {list(summary['beats_clinical'])}")
    print(f"[gene-sig] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
