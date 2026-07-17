"""Loader for the GSE59247 DCIS/IBC miRNA cohort as a LABELLED STATE SET.

Design note (deliberate): DCIS is loaded as a *state with labels* (tissue type, grade,
PAM50, source/batch, matched-patient key, outcome) -- NOT as a point on a linear
``healthy->NAT->DCIS->tumour`` trajectory. That ordering is biologically wrong (NAT is a
field-effect normal, not a DCIS precursor; DCIS->IDC is non-obligate and grade/subtype-
heterogeneous; the cohorts are cross-sectional and cross-platform). The correct downstream
use is per-edge TIMING (is a healthy->tumour change already present in DCIS = early, or only
in IDC = late), grade/subtype-stratified, with matched DCIS+IBC pairs and the progression
outcome -- never an assumed order. This loader just exposes the states + labels for that.

Integration with the TCGA spine is **rank-based** (the Agilent log-signal here is not on the
RPM scale): align by rank/quantile and read DIRECTION agreement, not absolute level.

Data: ``data/external/dcis_geo/GSE59247_series_matrix.txt.gz`` (miRNA, GPL15019).
"""

from __future__ import annotations

import gzip
import io
import re
from pathlib import Path
from typing import Tuple

import pandas as pd

from mirna_hallmark import config as C

_DCIS_GEO = C.REPO_ROOT / "data" / "external" / "dcis_geo"
GSE59247 = _DCIS_GEO / "GSE59247_series_matrix.txt.gz"          # miRNA (GPL15019)
GSE59246 = _DCIS_GEO / "GSE59246_series_matrix.txt.gz"          # mRNA  (GPL13607)
GPL13607_MAP = _DCIS_GEO / "GPL13607_feature_to_gene.tsv.gz"    # numeric feature -> gene
_PATIENT_RE = re.compile(r"^(?:DCIS|IBC|Normal)[-_](.+?)(?:_rep\d+)?$", re.I)


def _parse_series_matrix(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (data: probe x GSM, meta: GSM x field). Pure-stdlib GEO parser."""
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    gsm = titles = None
    char_lines = []
    t0 = t1 = None
    for i, ln in enumerate(lines):
        if ln.startswith("!Sample_geo_accession"):
            gsm = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics"):
            char_lines.append([x.strip().strip('"') for x in ln.rstrip("\n").split("\t")[1:]])
        elif ln.startswith("!series_matrix_table_begin"):
            t0 = i
        elif ln.startswith("!series_matrix_table_end"):
            t1 = i

    meta = {}
    for vals in char_lines:
        field = next((v.split(":", 1)[0].strip() for v in vals if ":" in v), None)
        if not field:
            continue
        meta[field] = [v.split(":", 1)[1].strip() if ":" in v else pd.NA for v in vals]
    meta_df = pd.DataFrame(meta, index=gsm)
    meta_df["title"] = titles
    # matched-patient / lesion key (shared between a DCIS and its synchronous IBC)
    meta_df["patient_key"] = [
        (_PATIENT_RE.match(t).group(1) if _PATIENT_RE.match(t) else t) for t in titles
    ]

    data = pd.read_csv(io.StringIO("".join(lines[t0 + 1 : t1])), sep="\t", index_col=0)
    data.index = [str(i).strip().strip('"') for i in data.index]
    data.columns = [str(c).strip().strip('"') for c in data.columns]
    return data, meta_df


def _label_meta(meta: pd.DataFrame) -> pd.DataFrame:
    """Attach authoritative ``state`` + ``rep_group`` keys (shared by both series).

    State comes from the curated ``tissue type`` field, NOT the title (titles are
    unreliable: GSM1431408 is titled DCIS but is IBC). ``rep_group`` = lesion base
    (title minus prefix minus ``_repN``) + state; it is the **cross-platform sample
    key** that pairs a miRNA GSM to its mRNA GSM in the companion series.
    """
    keep = {
        "tissue type": "tissue_type", "grade (combined)": "grade", "pam50": "pam50",
        "source": "source", "er (ihc)": "er", "pr (ihc)": "pr", "her2 (ihc)": "her2",
        "recurrence event": "recurrence", "dmfs event": "dmfs", "follow-up time": "followup",
        "title": "title", "patient_key": "patient_key",
    }
    meta = meta.rename(columns={k: v for k, v in keep.items() if k in meta.columns})
    meta = meta[[v for v in keep.values() if v in meta.columns]].copy()
    meta["state"] = meta["tissue_type"].map(
        lambda s: "DCIS" if "situ" in str(s).lower()
        else ("IBC" if "invasive" in str(s).lower() else "control")
    )
    base = (meta["title"].str.replace(r"(?i)_rep\d+", "", regex=True)
                         .str.replace(r"(?i)^(DCIS|IBC|Normal)[-_]", "", regex=True).str.strip())
    meta["rep_group"] = base + "|" + meta["state"]
    return meta


def _collapse_and_filter(mat: pd.DataFrame, meta: pd.DataFrame,
                         collapse_replicates: bool, drop_controls: bool):
    """Collapse technical replicates (by ``rep_group``, mean) and drop controls."""
    if collapse_replicates:
        cg = meta["rep_group"].reindex(mat.columns)
        mat = mat.T.groupby(cg).mean().T
        meta = meta.reset_index().rename(columns={"index": "gsm"}).groupby(
            "rep_group", as_index=True).agg("first")
        meta = meta.loc[mat.columns]
    if drop_controls:
        keep_s = meta.index[meta["state"] != "control"]
        mat, meta = mat[keep_s], meta.loc[keep_s]
    return mat, meta


def load_gse59247(
    path: Path = GSE59247,
    *,
    collapse_replicates: bool = True,
    drop_controls: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load GSE59247 miRNA -> (mirna: arm x sample log-signal, meta: sample x labels).

    - probes restricted to human (``hsa-``) mature arms (viral/control/Blank dropped);
    - **state from the curated ``tissue type``, NOT the title** (titles unreliable);
    - **technical replicates collapsed** by mean (``_repN`` pairs); Ambion controls dropped;
    - meta carries tissue_type/grade/pam50/source/ER/PR/HER2/recurrence/DMFS + ``rep_group``.
      NOTE: **no matched synchronous DCIS+IBC pairs** within a patient -- confirmed two ways:
      (i) title analysis = 0 shared lesion ids across states (the only "overlap", BRC/FW06,
      are collection-site prefixes, not patients); (ii) the GEO record states "48 arrays,
      41 unique tissues, independent groups" (4 normal/1 tissue, 28 DCIS/26, 16 inv/14 --
      exactly our 26 DCIS + 14 IBC after replicate collapse). The only within-patient layer
      is therefore the miRNA<->mRNA pairing with GSE59246 (``rep_group``).
    """
    data, meta = _parse_series_matrix(path)
    hsa = data.loc[data.index.str.startswith("hsa-")].apply(pd.to_numeric, errors="coerce")
    mirna = hsa.groupby(hsa.index).mean()  # collapse duplicate probe rows
    meta = _label_meta(meta)
    return _collapse_and_filter(mirna, meta, collapse_replicates, drop_controls)


def load_gse59246(
    path: Path = GSE59246,
    *,
    collapse_replicates: bool = True,
    drop_controls: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load GSE59246 mRNA -> (mrna: gene x sample log-signal, meta: sample x labels).

    Companion mRNA series to GSE59247 (same patients/lesions). GPL13607 is the
    *Feature-Number* Agilent 8x60K design, so the matrix row id is a numeric feature
    number; we map feature -> ``GeneName`` via the persisted ``GPL13607_feature_to_gene``
    table and collapse probes per gene by mean. ``rep_group`` matches GSE59247, so the
    two series can be paired for **within-sample miRNA<->mRNA coupling**.
    """
    data, meta = _parse_series_matrix(path)
    f2g = pd.read_csv(GPL13607_MAP, sep="\t", dtype=str).set_index("feature_id")["gene"]
    data.index = data.index.astype(str)
    genes = data.loc[data.index.isin(f2g.index)].apply(pd.to_numeric, errors="coerce")
    genes.index = pd.Index(f2g.loc[genes.index].values, name="gene")
    mrna = genes.groupby(genes.index).mean()  # collapse probes per gene
    meta = _label_meta(meta)
    return _collapse_and_filter(mrna, meta, collapse_replicates, drop_controls)


if __name__ == "__main__":
    mat, md = load_gse59247()
    print(f"miRNA matrix (reps collapsed, controls dropped): "
          f"{mat.shape[0]} hsa arms x {mat.shape[1]} samples")
    print("state x grade:\n", pd.crosstab(md["state"], md["grade"]))
    print("\nstate x pam50:\n", pd.crosstab(md["state"], md["pam50"]))
    print("\nNOTE: no reliable matched synchronous DCIS+IBC pairs -> use population-level "
          "timing (early/late) + grade/subtype stratification + outcomes, not within-patient.")
