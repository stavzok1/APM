"""Batch covariates for the CPTAC cohorts (prospective + TCGA-105).

"Batch" is assay-specific and has no single column in CPTAC. The canonical proteomics batch is
the **MS multiplex plex** (TMT/iTRAQ), and it is the *literal* batch where the design is known:

  - **prospective** (LinkedOmics HS_CPTAC_BRCA_2018, X01BR* cols): TMT plex from PDC study
    **PDC000120**, cached by ``scripts/cptac/fetch_cptac_brca_tmt_plex.py`` ->
    ``data/CPTAC/cptac_brca_prospective_tmt_plex.tsv`` (15 plexes).
  - **tcga105** (Mertins iTRAQ): plex = "Experiment Order" in
    ``CPTAC_TCGA_Breast_Cancer_iTRAQ_Sample_Mapping.xlsx`` (41 4-plex experiments).

The enrolling-**site** / TSS code parsed from the sample id (X01BR001 -> "01") is an
always-available accrual-level proxy shared by the miRNA-seq *and* protein assays — the fallback
when no plex design is present, and the only option for cohorts without a plex file.

``kind``: ``none`` | ``site`` | ``plex`` | ``auto`` (plex if its file exists, else site).
All dummy frames drop the most common level as reference and pool rare levels into ``other``.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import List, Sequence, Tuple

import pandas as pd

from pipeline.config import PATHS

MIN_SITE_N = 5
_CPTAC_DIR = PATHS.cptac_proteome_cct.parent
PROSPECTIVE_PLEX_TSV = _CPTAC_DIR / "cptac_brca_prospective_tmt_plex.tsv"
TCGA_ITRAQ_XLSX = _CPTAC_DIR / "CPTAC_TCGA_Breast_Cancer_iTRAQ_Sample_Mapping.xlsx"


def _onehot(labels: pd.Series, prefix: str, *, min_count: int = 1) -> pd.DataFrame:
    vc = labels.value_counts()
    keep = set(vc[vc >= min_count].index)
    lab = labels.where(labels.isin(keep), "other")
    if lab.nunique() < 2:
        return pd.DataFrame(index=labels.index)
    ref = lab.value_counts().index[0]  # drop most common as reference
    d = pd.get_dummies(lab, prefix=prefix).drop(columns=[f"{prefix}_{ref}"], errors="ignore")
    return d.astype(float)


# --------------------------------------------------------------------------- #
# Site (always available; accrual/processing proxy)
# --------------------------------------------------------------------------- #
def _site_code(s: str) -> str:
    m = re.match(r"X?(\d+)[A-Za-z]{2}\d", str(s))
    return m.group(1) if m else "other"


def site_dummies(samples: Sequence[str], min_count: int = MIN_SITE_N) -> pd.DataFrame:
    s = pd.Series({x: _site_code(x) for x in samples})
    return _onehot(s, "site", min_count=min_count)


# --------------------------------------------------------------------------- #
# Literal MS plex
# --------------------------------------------------------------------------- #
def prospective_plex_dummies(samples: Sequence[str]) -> pd.DataFrame:
    pm = pd.read_csv(PROSPECTIVE_PLEX_TSV, sep="\t").set_index("sample")["plex"].astype(str)
    s = pd.Series({x: pm.get(x, "other") for x in samples})
    return _onehot(s, "plex")


def tcga105_plex_dummies(samples: Sequence[str]) -> pd.DataFrame:
    x = pd.read_excel(TCGA_ITRAQ_XLSX, sheet_name="Proteome")
    part_to_plex = {}
    for _, r in x.iterrows():
        order = str(r["Experiment Order"])
        for ch in ("114-Biospecimen", "115-Biospecimen", "116-Biospecimen"):
            bc = str(r.get(ch, ""))
            if bc.startswith("TCGA"):
                part_to_plex[bc[:12]] = order
    s = pd.Series({x_: part_to_plex.get(str(x_)[:12], "other") for x_ in samples})
    return _onehot(s, "plex")


# --------------------------------------------------------------------------- #
# Dispatch
# --------------------------------------------------------------------------- #
def batch_dummies(cohort: str, samples: Sequence[str], kind: str = "auto") -> pd.DataFrame:
    samples = list(samples)
    if kind == "none":
        return pd.DataFrame(index=pd.Index(samples))
    if cohort in ("prospective", "pancan122"):
        if kind in ("plex", "auto") and PROSPECTIVE_PLEX_TSV.exists():
            return prospective_plex_dummies(samples)
        return site_dummies(samples)
    if cohort == "tcga105":
        if kind in ("plex", "auto") and TCGA_ITRAQ_XLSX.exists():
            return tcga105_plex_dummies(samples)
        return site_dummies(samples)  # TCGA participants carry no BR-site -> mostly 'other'
    return pd.DataFrame(index=pd.Index(samples))


def augment_cov(base_cov: pd.DataFrame, cohort: str, kind: str) -> Tuple[pd.DataFrame, List[str]]:
    """Return (base_cov + batch dummies joined on its index, list of added batch columns)."""
    if kind == "none":
        return base_cov, []
    b = batch_dummies(cohort, list(base_cov.index), kind=kind)
    if b.empty or b.shape[1] == 0:
        return base_cov, []
    return base_cov.join(b, how="left"), list(b.columns)
