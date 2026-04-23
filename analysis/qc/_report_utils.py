from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


def ensure_out_dir(run_id: str) -> Path:
    out = Path("analysis") / "qc" / "output" / run_id
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def md_link(path: Path, label: Optional[str] = None) -> str:
    lab = label or path.name
    return f"[{lab}]({path.name})"


def md_table_preview(df: pd.DataFrame, *, max_rows: int = 8, max_cols: int = 12) -> str:
    """
    Small markdown table preview (not for huge data).
    """
    if df.empty:
        return "_(empty)_\n"
    d2 = df.copy()
    if d2.shape[1] > max_cols:
        d2 = d2.iloc[:, :max_cols]
    if len(d2) > max_rows:
        d2 = d2.head(max_rows)
    return d2.to_markdown(index=False) + "\n"


def group_summary(
    df: pd.DataFrame,
    *,
    group_col: str,
    value_cols: List[str],
    min_n: int = 3,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for g, sub in df.groupby(group_col, dropna=False):
        rec: Dict[str, object] = {group_col: g, "n": int(len(sub))}
        if int(len(sub)) < int(min_n):
            rows.append(rec)
            continue
        for c in value_cols:
            if c not in sub.columns:
                continue
            v = pd.to_numeric(sub[c], errors="coerce")
            rec[f"{c}__mean"] = float(v.mean()) if v.notna().any() else np.nan
            rec[f"{c}__median"] = float(v.median()) if v.notna().any() else np.nan
        rows.append(rec)
    return pd.DataFrame(rows).sort_values(["n", group_col], ascending=[False, True])

