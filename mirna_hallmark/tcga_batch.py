"""TCGA-BRCA sequencing batch covariates for the subproject's TCGA spine (ready for full refresh).

The subproject's TCGA coupling is miRNA-arm (miRNA-Seq) vs target-gene RNA (RNA-Seq). The canonical
TCGA batch is the analyte **plate** (MBatch BatchId), encoded in the aliquot barcode field 6
(``TCGA-3C-AAAU-01A-11R-A41B-07`` -> plate ``A41B``). Each assay has its own plate, so a partial
correlation that wants to remove batch from *both* sides adjusts **both** the miRNA-seq and RNA-seq
plates (``kind="plate_both"``, the default). ``tss`` (tissue source site, barcode field 2) is the
coarser patient/accrual-level proxy shared by both assays.

Source map fetched by ``scripts/tcga/fetch_tcga_brca_seq_batch.py`` from the GDC files API ->
``data/TCGA/tcga_brca_seq_batch.tsv`` (1097 participants; 41 miRNA / 43 RNA plates). Plates with
``< min_count`` samples are pooled to ``other`` (default 30 -> ~14 levels per assay, <20).

``augment_cov(base_cov, kind)`` mirrors ``cptac_batch.augment_cov`` so the TCGA spine modules
(target_combined_anticorr, hallmark_interaction, robustness_checks, ...) can opt in via the same
``--batch-kind`` plumbing at the next full refresh. **Not yet wired into the spine** — this is the
ready, validated helper.

``kind`` ∈ ``none`` | ``tss`` | ``plate`` (RNA only) | ``mirna_plate`` | ``plate_both``.
"""
from __future__ import annotations

import functools
from typing import List, Sequence, Tuple

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.cptac_batch import _onehot

MIN_PLATE_N = 30
# Cross-state states (NAT n~110, GTEx n~350) have finer-grained, smaller batches than the
# ~1k-tumour spine, so they need a lower pooling threshold to retain any batch structure.
CROSS_STATE_MIN_N = 10
BATCH_TSV = C.REPO_ROOT / "data/TCGA/tcga_brca_seq_batch.tsv"


@functools.lru_cache(maxsize=1)
def _load_map() -> pd.DataFrame:
    if not BATCH_TSV.exists():
        raise FileNotFoundError(
            f"Missing {BATCH_TSV}. Build it: .venv/bin/python3 scripts/tcga/fetch_tcga_brca_seq_batch.py"
        )
    return pd.read_csv(BATCH_TSV, sep="\t").drop_duplicates("participant").set_index("participant")


@functools.lru_cache(maxsize=512)
def _batch_dummies_cached(participants: Tuple[str, ...], kind: str, min_count: int) -> pd.DataFrame:
    """Deterministic batch dummies for a fixed (participant-set, kind, min_count).

    Cached because partial-correlation loops call this once per gene with the *same*
    sample set (and the map is a disk read) -- without the cache the one-hot build and
    TSV read are repeated tens of thousands of times. The returned frame is treated as
    read-only by callers (``augment_cov`` joins onto a copy)."""
    m = _load_map().reindex(list(participants))
    cols = {"tss": ["tss"], "plate": ["rna_plate"], "mirna_plate": ["mirna_plate"],
            "plate_both": ["mirna_plate", "rna_plate"]}.get(kind)
    if cols is None:
        raise ValueError(f"kind must be none/tss/plate/mirna_plate/plate_both, got {kind!r}")
    frames = []
    for c in cols:
        prefix = {"rna_plate": "rplate", "mirna_plate": "mplate", "tss": "tss"}[c]
        d = _onehot(m[c].fillna("other"), prefix, min_count=min_count)
        if not d.empty:
            frames.append(d)
    return pd.concat(frames, axis=1) if frames else pd.DataFrame(index=m.index)


def tcga_batch_dummies(participants: Sequence[str], kind: str = "plate_both",
                       min_count: int = MIN_PLATE_N) -> pd.DataFrame:
    if kind == "none":
        return pd.DataFrame(index=pd.Index(list(participants)))
    return _batch_dummies_cached(tuple(participants), kind, min_count)


def augment_cov(base_cov: pd.DataFrame, kind: str,
                min_count: int = MIN_PLATE_N) -> Tuple[pd.DataFrame, List[str]]:
    """Return (base_cov + TCGA batch dummies joined on its participant index, added cols)."""
    if kind == "none":
        return base_cov, []
    b = tcga_batch_dummies(list(base_cov.index), kind=kind, min_count=min_count)
    if b.empty or b.shape[1] == 0:
        return base_cov, []
    return base_cov.join(b, how="left"), list(b.columns)


# --------------------------------------------------------------------------- #
# Cross-state (GTEx / NAT / tumor) per-state batch
#   tumor/NAT : analyte plate parsed from the TCGA barcode (NAT = the sample-11 aliquot),
#               keyed by participant *or* full barcode (the cross-state code uses both).
#   GTEx      : SMNABTCH (nucleic-acid) + SMGEBTCH (sequencing) from the sample-attribute
#               table, mapped donor -> batch on the small-RNA (SMLRNA) breast samples.
# --------------------------------------------------------------------------- #
@functools.lru_cache(maxsize=1)
def _tcga_plate_maps() -> dict:
    """{'01': {participant: plate}, '11': {participant: plate}} from the miRNA matrix
    barcodes (field 6 = analyte plate; field 4[:2] = sample type)."""
    cols = pd.read_csv(C.MIRNA_EXPRESSION, sep="\t", index_col=0, nrows=0).columns
    maps = {"01": {}, "11": {}}
    for bc in cols:
        parts = str(bc).split("-")
        if len(parts) < 7:
            continue
        st, plate, part = parts[3][:2], parts[5], "-".join(parts[:3])
        if st in maps:
            maps[st].setdefault(part, plate)
    return maps


@functools.lru_cache(maxsize=1)
def _gtex_batch_map() -> Tuple[dict, dict]:
    """(donor->SMNABTCH, donor->SMGEBTCH) over GTEx breast small-RNA (SMLRNA) samples."""
    df = pd.read_csv(C.GTEX_SAMPLE_ATTR, sep="\t", dtype=str)
    df = df[(df.get("SMTSD") == "Breast - Mammary Tissue") & (df.get("SMAFRZE") == "SMLRNA")].copy()
    df["donor"] = df["SAMPID"].str.split("-").str[:2].str.join("-")
    nab = dict(zip(df["donor"], df["SMNABTCH"].fillna("other")))
    geb = dict(zip(df["donor"], df["SMGEBTCH"].fillna("other")))
    return nab, geb


def _barcode_or_participant_plate(i: str, pmap: dict) -> str:
    s = str(i)
    if len(s) > 12 and s.count("-") >= 6:          # full aliquot barcode
        p = s.split("-")
        return p[5] if len(p) >= 7 else "other"
    return pmap.get(s, "other")                     # 12-char participant


def cross_state_batch_dummies(state: str, ids: Sequence[str],
                              min_count: int = CROSS_STATE_MIN_N) -> pd.DataFrame:
    """Per-state batch dummies for the cross-state coupling covariates.

    ``state`` in {tumor, nat, gtex}; ``ids`` = that state's sample/participant/donor index.
    Returns an ids-indexed 0/1 dummy frame (rare levels pooled to 'other')."""
    ids = list(ids)
    idx = pd.Index(ids)
    if not ids or not getattr(C, "TCGA_BATCH_KIND", "none") or C.TCGA_BATCH_KIND == "none":
        return pd.DataFrame(index=idx)

    if state == "gtex":
        nab, geb = _gtex_batch_map()
        frames = []
        for mp, prefix in ((nab, "gnab"), (geb, "ggeb")):
            d = _onehot(pd.Series({i: mp.get(str(i), "other") for i in ids}, dtype=object),
                        prefix, min_count=min_count)
            if not d.empty:
                frames.append(d)
        return pd.concat(frames, axis=1) if frames else pd.DataFrame(index=idx)

    if state in ("tumor", "nat"):
        pmap = _tcga_plate_maps()["01" if state == "tumor" else "11"]
        prefix = "tplate" if state == "tumor" else "nplate"
        lab = pd.Series({i: _barcode_or_participant_plate(i, pmap) for i in ids}, dtype=object)
        d = _onehot(lab, prefix, min_count=min_count)
        return d if not d.empty else pd.DataFrame(index=idx)

    return pd.DataFrame(index=idx)
