"""Roadmap Phase 3 outcome — METABRIC (partial/public, ~2,000) RFS + OS, the POWERED per-subtype cohort.

METABRIC partial (cBioPortal) has **1,002 RFS events + 1,144 OS deaths** across CLAUDIN subtypes with
hundreds each — ~10x TCGA's recurrence power, and genuinely per-subtype. It is **mRNA-only** (no miRNA),
so we validate the **mRNA-side** of the framework's prognostic signals:
  * the prognostic **target genes** from the gene-pressure hits (MKNK2/KCNJ8; MH-61) — expression vs outcome
  * the Hallmark **program scores** (mean-z), esp. **HYPOXIA** (miR-210's program), EMT, proliferation
vs **RFS (primary) + OS**, overall and per **CLAUDIN_SUBTYPE**, age-adjusted, BH-FDR over the 50 programs.

This is the powered, independent (different platform/cohort) outcome check; the miRNA/pressure side stays
TCGA/GSE22216 (METABRIC-full miRNA = EGA pending [[pending-controlled-data-requests]]).

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_metabric``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z

MET = C.REPO_ROOT / "data" / "METABRIC"
OUT_DIR = C.OUTPUT_ROOT / "outcome_metabric"
TARGET_GENES = ["MKNK2", "KCNJ8"]                       # gene-pressure prognostic hits (MH-61)
SUBTYPES = ["LumA", "LumB", "Basal", "Her2", "claudin-low"]
KEY = ["HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"]


def _load():
    c = pd.read_csv(MET / "data_clinical_patient.txt", sep="\t", comment="#", low_memory=False).set_index("PATIENT_ID")
    out = pd.DataFrame(index=c.index)
    out["rfs_t"] = pd.to_numeric(c["RFS_MONTHS"], errors="coerce")
    out["rfs_e"] = c["RFS_STATUS"].astype(str).str.split(":").str[0].replace({"nan": np.nan}).astype(float)
    out["os_t"] = pd.to_numeric(c["OS_MONTHS"], errors="coerce")
    out["os_e"] = c["OS_STATUS"].astype(str).str.split(":").str[0].replace({"nan": np.nan}).astype(float)
    out["age"] = _z(pd.to_numeric(c["AGE_AT_DIAGNOSIS"], errors="coerce"))
    out["subtype"] = c["CLAUDIN_SUBTYPE"]
    m = pd.read_csv(MET / "data_mrna_illumina_microarray.txt", sep="\t", low_memory=False)
    m = m.drop(columns=["Entrez_Gene_Id"]).rename(columns={"Hugo_Symbol": "gene"}).dropna(subset=["gene"])
    m = m.groupby("gene").mean(numeric_only=True)
    z = m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1).replace(0, np.nan), axis=0)   # gene z-scores
    return out, z


def _program_scores(z: pd.DataFrame) -> pd.DataFrame:
    hs = HallmarkSets.load()
    sets = hs.sets if hasattr(hs, "sets") else hs
    rows = {}
    for prog, genes in (sets.items() if isinstance(sets, dict) else [(p, hs.genes(p)) for p in hs.names()]):
        g = z.index.intersection(list(genes))
        if len(g) >= 5:
            rows[prog] = z.loc[g].mean(axis=0)
    return pd.DataFrame(rows).T                                                       # program x sample


def _screen(mat: pd.DataFrame, clin: pd.DataFrame, fam: str) -> pd.DataFrame:
    rows = []
    for feat in mat.index:
        x = mat.loc[feat]
        for ep, (t, e) in {"RFS": ("rfs_t", "rfs_e"), "OS": ("os_t", "os_e")}.items():
            # overall (age-adjusted)
            c, p = _cox_one(clin[t], clin[e], x, clin[["age"]])
            rec = {"family": fam, "feature": feat, "endpoint": ep, "stratum": "all", "hr": np.exp(c), "p": p,
                   "n_events": int(clin[e].sum())}
            rows.append(rec)
            for sub in SUBTYPES:
                ids = clin.index[clin["subtype"] == sub]
                cs, ps = _cox_one(clin.loc[ids, t], clin.loc[ids, e], x.reindex(ids), clin.loc[ids, ["age"]])
                rows.append({"family": fam, "feature": feat, "endpoint": ep, "stratum": sub,
                             "hr": np.exp(cs), "p": ps, "n_events": int(clin.loc[ids, e].sum())})
    r = pd.DataFrame(rows)
    for (ep, st), idx in r.groupby(["endpoint", "stratum"]).groups.items():
        ok = r.loc[idx, "p"].notna()
        r.loc[r.loc[idx].index[ok], "q"] = bh_fdr(r.loc[idx, "p"].dropna().to_numpy())
    return r


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    clin, z = _load()
    print(f"[metabric] {len(clin)} patients | RFS events {int(clin['rfs_e'].sum())} | OS events {int(clin['os_e'].sum())}")
    prog = _program_scores(z)
    print(f"[metabric] {prog.shape[0]} Hallmark program scores")
    res = pd.concat([_screen(prog, clin, "program"),
                     _screen(z.loc[z.index.intersection(TARGET_GENES)], clin, "target_gene")], ignore_index=True)
    res.to_csv(out_dir / "metabric_outcome.tsv", sep="\t", index=False)

    def hits(fam, ep, st, n=8):
        g = res[(res.family == fam) & (res.endpoint == ep) & (res.stratum == st) & (res.q < 0.1)]
        return g.sort_values("q").head(n).apply(lambda r: f"{r['feature']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist()
    summary = {
        "module": "mirna_hallmark.outcome_metabric",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_patients": int(len(clin)), "rfs_events": int(clin["rfs_e"].sum()), "os_events": int(clin["os_e"].sum()),
        "program_RFS_all_FDR": hits("program", "RFS", "all"),
        "program_OS_all_FDR": hits("program", "OS", "all"),
        "HYPOXIA_RFS": res[(res.feature == "HALLMARK_HYPOXIA") & (res.endpoint == "RFS")]
            .apply(lambda r: f"{r['stratum']} HR={r['hr']:.2f} q={r.get('q', np.nan):.3f} (ev={r['n_events']})", axis=1).tolist(),
        "target_genes_RFS": res[(res.family == "target_gene") & (res.endpoint == "RFS")]
            .apply(lambda r: f"{r['feature']}/{r['stratum']} HR={r['hr']:.2f} p={r['p']:.3f}", axis=1).tolist(),
        "caveats": ["METABRIC mRNA microarray (no miRNA/pressure) — validates the mRNA-side correlates only",
                    "program score = mean-z of member genes; age-adjusted; FDR within endpoint×stratum",
                    "CLAUDIN subtype (PAM50-like)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[metabric] program RFS FDR (all): {summary['program_RFS_all_FDR'][:6]}")
    print(f"[metabric] HYPOXIA RFS by subtype: {summary['HYPOXIA_RFS']}")
    print(f"[metabric] target genes RFS: {summary['target_genes_RFS']}")
    print(f"[metabric] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
