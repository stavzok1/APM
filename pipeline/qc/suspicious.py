from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _is_emptyish_scalar(x: Any) -> bool:
    if x is None:
        return True
    if isinstance(x, float) and np.isnan(x):
        return True
    if isinstance(x, str):
        s = x.strip()
        return s == "" or s == "[]" or s == "{}" or s.lower() == "nan"
    return False


def _is_emptyish_cell(x: Any) -> bool:
    if _is_emptyish_scalar(x):
        return True
    if isinstance(x, (list, tuple, set, dict)):
        return len(x) == 0
    return False


def summarize_suspicious_columns(
    df: pd.DataFrame,
    *,
    max_columns: int = 500,
) -> pd.DataFrame:
    """
    Column-level QC summary: NaN rate, emptyish rate, constant-ish.

    Returns a long DataFrame with one row per column.
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=["column", "n", "nan_frac", "emptyish_frac", "n_unique"])

    n = len(df)
    cols = list(df.columns)[:max_columns]
    out: List[Dict[str, Any]] = []
    for c in cols:
        s = df[c]
        nan_frac = float(s.isna().mean()) if n else 0.0
        try:
            emptyish_frac = float(s.map(_is_emptyish_cell).mean())
        except Exception:
            emptyish_frac = float("nan")
        try:
            n_unique = int(s.nunique(dropna=True))
        except Exception:
            n_unique = -1
        out.append(
            {
                "column": str(c),
                "n": int(n),
                "nan_frac": nan_frac,
                "emptyish_frac": emptyish_frac,
                "n_unique": n_unique,
            }
        )
    return pd.DataFrame(out)


def summarize_suspicious_nested_cells(
    df: pd.DataFrame,
    *,
    columns: Optional[Sequence[str]] = None,
    max_bad_examples: int = 5,
) -> pd.DataFrame:
    """
    Cell-level QC for nested columns that are expected to contain dict/list payloads.

    Flags:
    - non-empty strings that look like JSON/list but weren't decoded (e.g. '[{...}]' as str)
    - unexpected scalars where list/dict expected
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=["column", "n", "n_str_payload", "n_unexpected_scalar", "examples"])

    cols = list(columns) if columns else [c for c in df.columns if df[c].dtype == object]
    out: List[Dict[str, Any]] = []

    for c in cols:
        s = df[c]
        n = len(s)
        n_str_payload = 0
        n_unexpected_scalar = 0
        ex: List[str] = []
        for v in s.head(min(n, 50_000)).tolist():  # cap scan
            if isinstance(v, str):
                vv = v.strip()
                if vv.startswith("[") or vv.startswith("{"):
                    if vv not in ("[]", "{}"):
                        n_str_payload += 1
                        if len(ex) < max_bad_examples:
                            ex.append(vv[:160])
            elif isinstance(v, (list, dict)) or _is_emptyish_cell(v):
                continue
            else:
                n_unexpected_scalar += 1
                if len(ex) < max_bad_examples:
                    ex.append(repr(v)[:160])

        out.append(
            {
                "column": str(c),
                "n": int(n),
                "n_str_payload": int(n_str_payload),
                "n_unexpected_scalar": int(n_unexpected_scalar),
                "examples": "; ".join(ex),
            }
        )
    return pd.DataFrame(out)

