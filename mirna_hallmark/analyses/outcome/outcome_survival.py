"""Roadmap Phase 3 (outcome component): do the framework's headline miRNA arms carry PROGNOSTIC
signal — and does it survive composition adjustment?

Two cohorts, different confound structures (apm-rigor-protocol triangulation):
  * **TCGA-BRCA** overall survival (`os_event`/`os_time_days`) — Cox per headline arm, crude
    (age+stage) vs composition-adjusted (+epi/immune/stroma/proliferation metagenes).
  * **GSE22216** (METABRIC-miRNA, 210) distant-relapse-free survival — independent replication
    (Cox + age/grade/ER/nodes/size). Arms matched by base name (v16-style, no arm suffix) — caveat.

This fills the framework's shallow-outcome gap (TCGA miRNA→survival was untested; GSE22216 DRFS unused).

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.outcome.outcome_survival --max-arms 200``
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cross_state.cross_state_coupling import EPI_MARKERS, IMMUNE_MARKERS, STROMA_MARKERS, _metagene, _prolif_metagene
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr

OUT_DIR = C.OUTPUT_ROOT / "outcome_survival"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
GSE = C.REPO_ROOT / "data" / "external_cache" / "geo" / "GSE22216_mirna_arm_matrix.tsv.gz"
GSE_SERIES = C.REPO_ROOT / "data" / "external" / "dcis_geo" / "GSE22216_series_matrix.txt.gz"


def _z(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / (s.std(ddof=0) or 1.0)


def _cox_one(dur, evt, x, cov: pd.DataFrame) -> tuple[float, float]:
    """Cox HR (per SD) + p for x, adjusting for cov. Returns (coef, p) for x."""
    df = pd.concat([dur.rename("T"), evt.rename("E"), _z(x).rename("x"), cov], axis=1).dropna()
    if df["E"].sum() < 8 or df.shape[0] < 25 or df["x"].std() == 0:
        return np.nan, np.nan
    df = df.loc[:, df.std(numeric_only=True) > 0]                # drop constant dummies
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph = CoxPHFitter(penalizer=0.1).fit(df, "T", "E")
        return float(cph.params_["x"]), float(cph.summary.loc["x", "p"])
    except Exception:
        return np.nan, np.nan


def _headline_arms(mirna_index, max_arms: int) -> list[str]:
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e = e[(e["rho"] < 0) & (e["q_by"] < 0.05)]
    arms = e["key"].str.split(r"\|\|").str[0]
    counts = arms.value_counts()
    ordered = [a for a in counts.index if a in set(mirna_index)]
    return ordered[:max_arms]


def load_tcga_outcome() -> dict:
    """Shared TCGA outcome covariates from the full raw clinical: DFI (recurrence) + OS endpoints,
    age + finer TNM (base), and composition metagenes (comp). Reused by arm- and pressure-survival."""
    rna = D.load_rna(); rna = rna[~rna.index.duplicated(keep="first")]
    c = pd.read_csv(C.REPO_ROOT / "annotations" / "BRCA_clinical", sep="\t", low_memory=False)
    c["participant"] = c["sampleID"].astype(str).str[:12]
    c = c.drop_duplicates("participant").set_index("participant")
    dtl = pd.to_numeric(c["days_to_last_followup"], errors="coerce")
    dfi_e = c["new_tumor_event_after_initial_treatment"].astype(str).str.upper().eq("YES").astype(float)
    dfi_e[c["new_tumor_event_after_initial_treatment"].isna()] = np.nan
    dfi_t = pd.to_numeric(c["days_to_new_tumor_event_after_initial_treatment"], errors="coerce").fillna(dtl)
    os_e = c["vital_status"].astype(str).str.upper().eq("DECEASED").astype(float)
    os_t = pd.to_numeric(c["days_to_death"], errors="coerce").fillna(dtl)
    age = _z(pd.to_numeric(c.get("Age_at_Initial_Pathologic_Diagnosis_nature2012"), errors="coerce"))

    def tnm_dum(col, n):
        s = c[col].astype(str).str.extract(rf"({col[-1]}[0-{n}])", expand=False)
        return pd.get_dummies(s, prefix=col[-1]).astype(float)
    tnm = pd.concat([tnm_dum("pathologic_T", 4), tnm_dum("pathologic_N", 3), tnm_dum("pathologic_M", 1)], axis=1)
    comp = pd.DataFrame({"epi": _metagene(rna, EPI_MARKERS), "imm": _metagene(rna, IMMUNE_MARKERS),
                         "str": _metagene(rna, STROMA_MARKERS), "prolif": _prolif_metagene(rna, HallmarkSets.load())})
    return {"dfi_t": dfi_t, "dfi_e": dfi_e, "os_t": os_t, "os_e": os_e,
            "base": pd.concat([age.rename("age"), tnm], axis=1), "comp": comp}


def _tcga(arms: list[str]) -> pd.DataFrame:
    mirna = D.load_mirna_arms(); mirna = mirna[~mirna.index.duplicated(keep="first")]
    cl = load_tcga_outcome()
    dfi_t, dfi_e, os_t, os_e = cl["dfi_t"], cl["dfi_e"], cl["os_t"], cl["os_e"]
    base = cl["base"]; adj = pd.concat([base, cl["comp"]], axis=1)
    rows = []
    for a in arms:
        x = mirna.loc[a]
        cb, pb = _cox_one(dfi_t, dfi_e, x, base)          # DFI crude (age+TNM)
        ca, pa = _cox_one(dfi_t, dfi_e, x, adj)           # DFI composition-adjusted
        oc, op = _cox_one(os_t, os_e, x, adj)             # OS (secondary)
        rows.append({"arm": a, "tcga_dfi_hr_base": np.exp(cb), "tcga_dfi_p_base": pb,
                     "tcga_dfi_hr_adj": np.exp(ca), "tcga_dfi_p_adj": pa,
                     "tcga_os_hr_adj": np.exp(oc), "tcga_os_p_adj": op})
    r = pd.DataFrame(rows)
    for col in ("tcga_dfi_p_adj", "tcga_os_p_adj"):
        ok = r[col].notna(); r.loc[ok, col.replace("_p_", "_q_")] = bh_fdr(r.loc[ok, col].to_numpy())
    print(f"[outcome] TCGA DFI events: {int(dfi_e.sum())} | OS events: {int(os_e.sum())}")
    return r


def _gse22216(arms: list[str]) -> pd.DataFrame:
    m = pd.read_csv(GSE, sep="\t", index_col=0)
    lines = [l.rstrip("\n") for l in gzip.open(GSE_SERIES, "rt")]
    gsm = [x.strip('"') for x in [l for l in lines if l.startswith("!Sample_geo_accession")][0].split("\t")[1:]]

    def char(key):
        row = [l for l in lines if l.startswith("!Sample_characteristics") and key in l][0].split("\t")[1:]
        vals = [x.strip('"').split(":", 1)[1].strip() if ":" in x else np.nan for x in row]
        return pd.Series(pd.to_numeric(vals, errors="coerce"), index=gsm)

    meta = pd.DataFrame({"T": char("distant-relapse free survival"), "E": char("distant-relapse event"),
                         "age": char("patient age"), "grade": char("tumour grade"),
                         "er": char("er status"), "nodes": char("nodes involved"), "size": char("tumour size")})
    cov = pd.concat([_z(meta["age"]).rename("age"), _z(meta["grade"]).rename("grade"),
                     meta["er"].rename("er"), _z(meta["nodes"]).rename("nodes"), _z(meta["size"]).rename("size")], axis=1)
    # base-name match (GSE22216 lacks -3p/-5p arm suffix): hsa-miR-29c-3p -> hsa-miR-29c
    def base(a): return re.sub(r"-(3p|5p)$", "", a)
    gse_base = {base(i): i for i in m.index}
    rows = []
    for a in arms:
        gi = gse_base.get(base(a))
        if gi is None:
            rows.append({"arm": a, "gse_hr": np.nan, "gse_p": np.nan}); continue
        x = m.loc[gi].reindex(meta.index)
        c, p = _cox_one(meta["T"], meta["E"], x, cov)
        rows.append({"arm": a, "gse_hr": np.exp(c), "gse_p": p})
    r = pd.DataFrame(rows)
    ok = r["gse_p"].notna()
    r.loc[ok, "gse_q"] = bh_fdr(r.loc[ok, "gse_p"].to_numpy())
    return r


def run(*, out_dir: Path = OUT_DIR, max_arms: int = 200) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    mirna = D.load_mirna_arms()
    arms = _headline_arms(mirna.index, max_arms)
    print(f"[outcome] testing {len(arms)} headline arms for survival in TCGA (OS) + GSE22216 (DRFS)")
    res = _tcga(arms).merge(_gse22216(arms), on="arm", how="outer")
    res.to_csv(out_dir / "arm_survival.tsv", sep="\t", index=False)

    sig_dfi = res[(res["tcga_dfi_q_adj"] < C.FDR_ALPHA)]
    # primary cross-cohort: matched RECURRENCE endpoints (TCGA DFI <-> GSE22216 distant-relapse), same direction
    conc = res[(res["tcga_dfi_p_adj"] < 0.05) & (res["gse_p"] < 0.05) &
               (np.sign(np.log(res["tcga_dfi_hr_adj"])) == np.sign(np.log(res["gse_hr"])))]
    summary = {
        "module": "mirna_hallmark.analyses.outcome.outcome_survival",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_arms": len(arms),
        "tcga_DFI_recurrence": {
            "endpoint": "new_tumor_event_after_initial_treatment (TNM+age adj, then +composition)",
            "n_FDR_sig_composition_adjusted": int(sig_dfi.shape[0]),
            "top": sig_dfi.reindex(sig_dfi["tcga_dfi_p_adj"].sort_values().index)
                   .head(10).apply(lambda r: f"{r['arm']} HR={r['tcga_dfi_hr_adj']:.2f} q={r['tcga_dfi_q_adj']:.3f}", axis=1).tolist(),
            "killed_by_composition": int(res[(res["tcga_dfi_p_base"] < 0.05) & (res["tcga_dfi_p_adj"] > 0.1)].shape[0]),
        },
        "tcga_OS_secondary": {"n_FDR_sig": int((res["tcga_os_q_adj"] < C.FDR_ALPHA).sum())},
        "gse22216_DRFS": {"n_matched": int(res["gse_p"].notna().sum()),
                          "n_nominal_sig": int((res["gse_p"] < 0.05).sum()),
                          "FDR_hits": res[res["gse_q"] < 0.1].sort_values("gse_q")
                              .apply(lambda r: f"{r['arm']} HR={r['gse_hr']:.2f} q={r['gse_q']:.3f}", axis=1).tolist()},
        "cross_cohort_concordant_recurrence": conc.sort_values("tcga_dfi_p_adj").apply(
            lambda r: f"{r['arm']} (TCGA DFI HR {r['tcga_dfi_hr_adj']:.2f} / GSE DRFS HR {r['gse_hr']:.2f})", axis=1).tolist(),
        "caveats": ["Cox per arm, penalizer=0.1; arm z-scored (HR per SD)",
                    "GSE22216 base-name match (no arm suffix) — approximate; DRFS in the 210 open subset",
                    "TCGA DFI from new_tumor_event (richer raw clinical); composition = metagene proxy (deconvolution = Phase 0/1-refined)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[outcome] TCGA DFI: {sig_dfi.shape[0]} arms FDR-sig (composition-adjusted); "
          f"{summary['tcga_DFI_recurrence']['killed_by_composition']} killed by composition")
    print(f"[outcome] GSE22216 DRFS: {summary['gse22216_DRFS']['n_nominal_sig']}/{summary['gse22216_DRFS']['n_matched']} nominal")
    print(f"[outcome] cross-cohort concordant (recurrence): {summary['cross_cohort_concordant_recurrence'][:8]}")
    print(f"[outcome] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--max-arms", type=int, default=200)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, max_arms=args.max_arms)


if __name__ == "__main__":
    main()
