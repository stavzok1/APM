"""Roadmap Phase 3 outcome — FOCUSED per-subtype confirmation of cohort-wide survival hits.

Per-subtype power is limited (TCGA LumA 43 DFI / 63 OS events; LumB/Basal ~20-38; Her2 too small —
see the event-count check). So we do NOT sweep per subtype. Instead we take the **cohort-wide hits**
(from `outcome_survival` arms + `outcome_pressure_survival` resolutions) and confirm/localize ONLY those:
  * TCGA: LumA (primary), LumB / Basal (exploratory) — Cox, age-adjusted only (minimal, given few events).
  * GSE22216: ER+ vs ER- (the better-powered split in the outcome-grade cohort) — for expression hits.

Run (after the two cohort-wide modules): ``.venv/bin/python3 -m mirna_hallmark.outcome_subtype``
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
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome, GSE, GSE_SERIES
from mirna_hallmark.analyses.outcome.outcome_pressure_survival import _feature_matrices

OUT_DIR = C.OUTPUT_ROOT / "outcome_subtype"
SWEEP = C.OUTPUT_ROOT / "outcome_pressure_survival" / "pressure_resolution_survival.tsv"
ARMS = C.OUTPUT_ROOT / "outcome_survival" / "arm_survival.tsv"
TCGA_SUBTYPES = ["LumA", "LumB", "Basal"]            # Her2 too few events (8 DFI)
MAX_HITS = 15


def _select_hits() -> pd.DataFrame:
    """Cohort-wide hits to confirm: FDR<0.1 in either endpoint, else top nominal DFI."""
    sw = pd.read_csv(SWEEP, sep="\t")
    fdr = sw[(sw["dfi_q_adj"] < 0.1) | (sw["os_q_adj"] < 0.1)]
    if len(fdr) < 5:
        fdr = pd.concat([fdr, sw.sort_values("dfi_p_adj").head(MAX_HITS)]).drop_duplicates(["resolution", "feature"])
    return fdr.head(MAX_HITS)[["resolution", "feature"]]


def _tcga_subtype(hits: pd.DataFrame) -> pd.DataFrame:
    cl = load_tcga_outcome()
    mats = _feature_matrices()
    strata = D.load_clinical_strata()
    strata = strata.set_index("participant") if "participant" in strata.columns else strata
    pam = strata["PAM50_final"]
    age = cl["base"]["age"]
    rows = []
    for h in hits.itertuples(index=False):
        if h.resolution not in mats or h.feature not in mats[h.resolution].index:
            continue
        x = mats[h.resolution].loc[h.feature].astype(float)
        for sub in TCGA_SUBTYPES:
            ids = pam.index[pam == sub]
            cov = age.reindex(ids).to_frame("age")
            for ep, (t, e) in {"DFI": (cl["dfi_t"], cl["dfi_e"]), "OS": (cl["os_t"], cl["os_e"])}.items():
                c, p = _cox_one(t.reindex(ids), e.reindex(ids), x.reindex(ids), cov)
                rows.append({"resolution": h.resolution, "feature": h.feature, "subtype": sub,
                             "endpoint": ep, "hr": np.exp(c), "p": p,
                             "n_events": int(e.reindex(ids).sum())})
    return pd.DataFrame(rows)


def _gse_subtype(hits: pd.DataFrame) -> pd.DataFrame:
    m = pd.read_csv(GSE, sep="\t", index_col=0)
    lines = [l.rstrip("\n") for l in gzip.open(GSE_SERIES, "rt")]
    gsm = [x.strip('"') for x in [l for l in lines if l.startswith("!Sample_geo_accession")][0].split("\t")[1:]]

    def char(key):
        row = [l for l in lines if l.startswith("!Sample_characteristics") and key in l][0].split("\t")[1:]
        return pd.Series(pd.to_numeric([x.strip('"').split(":", 1)[1].strip() if ":" in x else np.nan for x in row], errors="coerce"), index=gsm)

    er, t, e = char("er status"), char("distant-relapse free survival"), char("distant-relapse event")
    expr_hits = hits[hits["resolution"].astype(str).str.startswith("expr")]["feature"].tolist()
    expr_hits = list(dict.fromkeys(expr_hits + ["hsa-miR-210-3p", "hsa-miR-148a-3p", "hsa-miR-150-5p"]))  # ensure the GSE FDR hits
    base = {re.sub(r"-(3p|5p)$", "", i): i for i in m.index}
    rows = []
    for a in expr_hits:
        gi = base.get(re.sub(r"-(3p|5p)$", "", str(a)))
        if gi is None:
            continue
        x = m.loc[gi].reindex(gsm)
        for grp, sel in {"ER+": er == 1, "ER-": er == 0}.items():
            ids = er.index[sel]
            c, p = _cox_one(t.reindex(ids), e.reindex(ids), x.reindex(ids), pd.DataFrame(index=ids))
            rows.append({"feature": a, "stratum": grp, "hr": np.exp(c), "p": p, "n_events": int(e.reindex(ids).sum())})
    return pd.DataFrame(rows)


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = _select_hits()
    print(f"[subtype] confirming {len(hits)} cohort-wide hits per subtype")
    tc = _tcga_subtype(hits); tc.to_csv(out_dir / "tcga_per_subtype.tsv", sep="\t", index=False)
    gs = _gse_subtype(hits); gs.to_csv(out_dir / "gse22216_per_ER.tsv", sep="\t", index=False)

    tc_sig = tc[(tc["p"] < 0.05)].sort_values("p")
    gs_sig = gs[(gs["p"] < 0.05)].sort_values("p")
    summary = {
        "module": "mirna_hallmark.outcome_subtype",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "design": "focused confirmation of cohort-wide hits; NOT a per-subtype sweep (power-limited)",
        "hits_tested": hits.apply(lambda r: f"{r['resolution']}:{r['feature']}", axis=1).tolist(),
        "tcga_subtype_nominal": tc_sig.apply(lambda r: f"{r['feature']}[{r['resolution']}] {r['subtype']}/{r['endpoint']} HR={r['hr']:.2f} p={r['p']:.3f} (ev={r['n_events']})", axis=1).tolist(),
        "gse22216_ER_nominal": gs_sig.apply(lambda r: f"{r['feature']} {r['stratum']} HR={r['hr']:.2f} p={r['p']:.3f} (ev={r['n_events']})", axis=1).tolist(),
        "caveats": ["EXPLORATORY: per-subtype events few (LumA 43 DFI/63 OS, LumB/Basal ~20-38, Her2 excluded)",
                    "age-only adjustment within subtype (cannot afford composition+TNM); nominal p (no per-subtype FDR)",
                    "real per-subtype power needs METABRIC-full (PAM50 + follow-up; EGA pending [[pending-controlled-data-requests]])"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[subtype] TCGA nominal: {summary['tcga_subtype_nominal'][:8]}")
    print(f"[subtype] GSE22216 ER nominal: {summary['gse22216_ER_nominal']}")
    print(f"[subtype] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
