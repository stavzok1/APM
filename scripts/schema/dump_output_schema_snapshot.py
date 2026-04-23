from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _shape(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): _shape(v) for k, v in sorted(x.items(), key=lambda kv: str(kv[0]))}
    if isinstance(x, list):
        if not x:
            return []
        return [_shape(x[0])]
    return type(x).__name__


def _maybe_decode_nested_cell(v: Any) -> Any:
    """
    Decode nested cells that were serialized as JSON strings (preferred) or repr(list/dict)
    in CSV outputs (some legacy paths). If it doesn't look like encoded nested data, return as-is.
    """
    if not isinstance(v, str):
        return v
    s = v.strip()
    if not s:
        return v
    if s == "[]" or s == "{}":
        try:
            return json.loads(s)
        except Exception:
            return v
    if s[0] in "[{":
        # Prefer JSON; fall back to literal_eval (repr-style) if needed.
        try:
            return json.loads(s)
        except Exception:
            try:
                return ast.literal_eval(s)
            except Exception:
                return v
    return v


def _first_non_null(series) -> Optional[Any]:
    for v in series:
        if v is None:
            continue
        if isinstance(v, float) and v != v:  # NaN
            continue
        return v
    return None


def _snapshot_parquet(path: Path, *, nested_cols: Iterable[str]) -> Dict[str, Any]:
    import pandas as pd
    import pyarrow.parquet as pq

    available = set(pq.ParquetFile(path).schema.names)
    wanted = [c for c in nested_cols if c in available]
    missing = [c for c in nested_cols if c not in available]

    df = pd.read_parquet(path, columns=wanted)
    snap: Dict[str, Any] = {}
    for col in missing:
        snap[col] = {"__missing_column__": True}

    for col in wanted:
        v = _first_non_null(df[col].tolist())
        v = _maybe_decode_nested_cell(v)
        snap[col] = _shape(v) if v is not None else {"__all_null__": True}
    return snap


def build_snapshot() -> Dict[str, Any]:
    """
    Snapshot nested key structures from actual on-disk outputs.

    This is intentionally "schema-only": it reads a small set of columns and
    uses only the first non-null row to infer nested key structure.
    """
    from pipeline.config import PATHS

    out: Dict[str, Any] = {}

    # 1) Element focus table (primary nested evidence contract)
    elem_pq = Path(PATHS.regulatory_elements_table_with_evidence_parquet)
    if not elem_pq.exists():
        raise FileNotFoundError(
            f"Missing element focus parquet: {elem_pq}. "
            f"Run the pipeline Step 12 to generate it."
        )
    out["elem_focus_parquet"] = {
        "path": str(elem_pq),
        "nested_columns": _snapshot_parquet(
            elem_pq,
            nested_cols=[
                "gene_links",
                "TAD_domains",
                "TAD_boundary_overlaps",
                "atac_peaks",
                "atac_peak_links",
                "chip_hits",
                "ABC_enhancers",
                "hichip",
                "screen_exp",
                "screen_comp",
            ],
        ),
    }

    # 2) ATAC peaks annotated table (if present)
    atac_pq = Path(PATHS.working_dir) / "atac_peaks" / "atac_peaks_annotated.parquet"
    if atac_pq.exists():
        out["atac_peaks_parquet"] = {
            "path": str(atac_pq),
            "nested_columns": _snapshot_parquet(
                atac_pq,
                nested_cols=[
                    "gene_links",
                    "genes_by_tier",
                    "lncrna_links",
                    "lncrnas_by_tier",
                    "ccre_links",
                    "ccre_types",
                    "TAD_domains",
                    "TAD_boundary_overlaps",
                ],
            ),
        }
    else:
        out["atac_peaks_parquet"] = {"path": str(atac_pq), "__missing__": True}

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Dump schema snapshot from actual pipeline outputs.")
    ap.add_argument(
        "--out",
        type=str,
        default="tests/fixtures/output_schema_snapshot.json",
        help="Output JSON path (relative to repo root by default).",
    )
    ap.add_argument("--check", action="store_true", help="Fail if snapshot differs; do not write.")
    args = ap.parse_args()

    out_path = (REPO / args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    snap = build_snapshot()
    rendered = json.dumps(snap, indent=2, sort_keys=True) + "\n"

    if out_path.exists():
        old = out_path.read_text(encoding="utf-8")
        if old == rendered:
            print("[OK] output schema snapshot already up to date")
            return
        if args.check:
            raise SystemExit(f"[FAIL] output schema snapshot differs: {out_path}")

    out_path.write_text(rendered, encoding="utf-8")
    print(f"[OK] wrote: {out_path}")


if __name__ == "__main__":
    main()

