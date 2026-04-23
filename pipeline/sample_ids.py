from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Optional, Set, Tuple


_TCGA_PREFIX = "TCGA-"


@dataclass(frozen=True)
class TcgaIds:
    """
    Normalized TCGA identifiers at multiple granularities.

    Levels:
      - participant: TCGA-XX-YYYY
      - sample:      TCGA-XX-YYYY-SS   (SS = 2-digit sample type, e.g. 01/10/11)
      - sample_vial: TCGA-XX-YYYY-SSA  (adds vial letter, e.g. 01A/10A/11A)
      - aliquot:     TCGA-XX-YYYY-SSA-PP... (full aliquot barcode; variable tail)
    """

    raw: str
    participant: Optional[str]
    sample: Optional[str]
    sample_vial: Optional[str]
    aliquot: Optional[str]


_participant_re = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})$")
_sample_re = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-([0-9]{2})$")
_sample_vial_re = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-([0-9]{2})([A-Z])$")
_aliquot_re = re.compile(r"^(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})-([0-9]{2})([A-Z])-(.+)$")


def normalize_tcga_id(raw: str) -> TcgaIds:
    """
    Accepts any of the observed TCGA id shapes and derives joinable keys.

    Examples:
      - TCGA-A7-A0CH                       -> participant
      - TCGA-A7-A0CH-01                    -> sample
      - TCGA-A7-A0CH-01A                   -> sample_vial
      - TCGA-DD-A1EI-11A-11R-A131-07       -> aliquot (+ all higher levels)
    """
    raw = (raw or "").strip()
    if not raw.startswith(_TCGA_PREFIX):
        return TcgaIds(raw=raw, participant=None, sample=None, sample_vial=None, aliquot=None)

    m = _participant_re.match(raw)
    if m:
        return TcgaIds(raw=raw, participant=m.group(1), sample=None, sample_vial=None, aliquot=None)

    m = _sample_re.match(raw)
    if m:
        participant = m.group(1)
        sample = f"{participant}-{m.group(2)}"
        return TcgaIds(raw=raw, participant=participant, sample=sample, sample_vial=None, aliquot=None)

    m = _sample_vial_re.match(raw)
    if m:
        participant = m.group(1)
        sample = f"{participant}-{m.group(2)}"
        sample_vial = f"{participant}-{m.group(2)}{m.group(3)}"
        return TcgaIds(raw=raw, participant=participant, sample=sample, sample_vial=sample_vial, aliquot=None)

    m = _aliquot_re.match(raw)
    if m:
        participant = m.group(1)
        sample = f"{participant}-{m.group(2)}"
        sample_vial = f"{participant}-{m.group(2)}{m.group(3)}"
        return TcgaIds(raw=raw, participant=participant, sample=sample, sample_vial=sample_vial, aliquot=raw)

    # Fallback: parse by splitting (keeps behavior stable even if tail changes)
    parts = raw.split("-")
    if len(parts) >= 3 and parts[0] == "TCGA":
        participant = "-".join(parts[:3])
        sample = None
        sample_vial = None
        aliquot = None
        if len(parts) >= 4 and len(parts[3]) >= 2:
            ss = parts[3][:2]
            if ss.isdigit():
                sample = f"{participant}-{ss}"
                if len(parts[3]) >= 3 and parts[3][2].isalpha():
                    sample_vial = f"{participant}-{ss}{parts[3][2].upper()}"
        if len(parts) >= 5:
            aliquot = raw
        return TcgaIds(raw=raw, participant=participant, sample=sample, sample_vial=sample_vial, aliquot=aliquot)

    return TcgaIds(raw=raw, participant=None, sample=None, sample_vial=None, aliquot=None)


def tcga_sample_type_two_digit(raw: str) -> Optional[str]:
    """
    Two-digit TCGA sample type embedded in the barcode (e.g. 01 primary tumor,
    10 blood-derived normal, 11 solid-tissue normal).

    Parsed from ``sample_vial`` (…-01A → 01) or ``sample`` (…-01 → 01) when available.
    """
    tid = normalize_tcga_id(str(raw).strip())
    s = tid.sample_vial or tid.sample or tid.participant
    if not s or not str(s).startswith(_TCGA_PREFIX):
        return None
    parts = str(s).split("-")
    if len(parts) >= 4 and len(parts[3]) >= 2:
        two = parts[3][:2]
        if two.isdigit():
            return two
    return None


def normalize_tcga_to_sample_vial_or_sample(raw: str) -> Optional[str]:
    """Prefer sample_vial, else sample, for stable per-specimen keys when classifying types."""
    tid = normalize_tcga_id(str(raw).strip())
    v = tid.sample_vial or tid.sample
    if v and str(v).startswith(_TCGA_PREFIX):
        return str(v)
    return None


def count_unique_tcga_sample_types(
    raw_ids: Iterable[str],
) -> Tuple[int, int, int, int, int]:
    """
    Count unique specimen keys (sample_vial or sample) by TCGA sample type code.

    Returns:
        (n_total_unique, n_primary_01, n_normal_blood_10, n_normal_solid_11, n_other_or_unknown)
    """
    uniq: Set[str] = set()
    for x in raw_ids:
        key = normalize_tcga_to_sample_vial_or_sample(str(x).strip())
        if key:
            uniq.add(key)

    n01 = n10 = n11 = nother = 0
    for key in uniq:
        code = tcga_sample_type_two_digit(key)
        if code == "01":
            n01 += 1
        elif code == "10":
            n10 += 1
        elif code == "11":
            n11 += 1
        else:
            nother += 1
    return (len(uniq), n01, n10, n11, nother)


def tcga_best_join_key(raw: str) -> Optional[str]:
    """
    Recommended single key for joining across your pipeline:

    Prefer sample_vial (01A/10A/11A) when available, else fall back to sample (01/10/11),
    else participant.
    """
    ids = normalize_tcga_id(raw)
    return ids.sample_vial or ids.sample or ids.participant


def add_tcga_id_columns_inplace(
    df,
    *,
    raw_id_col: Optional[str] = None,
    raw_id: Optional[str] = None,
    overwrite: bool = False,
) -> None:
    """
    Add normalized TCGA id columns to a DataFrame in-place.

    Produces columns:
      - participant
      - sample
      - sample_vial
      - aliquot

    Source can be:
      - a scalar raw_id (applied to all rows), OR
      - a column raw_id_col in df (row-wise)
    """
    import pandas as pd  # local import to avoid hard dependency at module import time

    if raw_id is None and raw_id_col is None:
        return

    if raw_id_col is not None and raw_id_col in df.columns:
        src = df[raw_id_col].astype(str)
    elif raw_id is not None:
        src = pd.Series([raw_id] * len(df), index=df.index, dtype="string")
    else:
        return

    def _maybe_set(col: str, values) -> None:
        if overwrite or col not in df.columns:
            df[col] = values

    ids = src.map(normalize_tcga_id)
    _maybe_set("participant", ids.map(lambda x: x.participant))
    _maybe_set("sample", ids.map(lambda x: x.sample))
    _maybe_set("sample_vial", ids.map(lambda x: x.sample_vial))
    _maybe_set("aliquot", ids.map(lambda x: x.aliquot))


def add_tcga_id_columns_from_index_inplace(df, *, overwrite: bool = False) -> None:
    """
    Convenience: treat df.index as the raw TCGA identifier source.
    """
    import pandas as pd

    src = pd.Series(df.index.astype(str), index=df.index, dtype="string")
    tmp_col = "__tmp_tcga_raw__"
    df[tmp_col] = src.values
    try:
        add_tcga_id_columns_inplace(df, raw_id_col=tmp_col, overwrite=overwrite)
    finally:
        df.drop(columns=[tmp_col], inplace=True, errors="ignore")

