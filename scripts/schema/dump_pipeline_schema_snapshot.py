from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _shape(x: Any) -> Any:
    """
    Convert an arbitrary nested python object (dict/list/scalar) into a
    stable "shape" representation (keys + nested key structure).
    """
    if isinstance(x, dict):
        out: Dict[str, Any] = {}
        for k in sorted(x.keys(), key=lambda s: str(s)):
            out[str(k)] = _shape(x[k])
        return out
    if isinstance(x, list):
        if not x:
            return []
        # Use the first element as a representative schema.
        return [_shape(x[0])]
    # Scalars: just record the type name (None -> "NoneType")
    return type(x).__name__


def build_snapshot() -> Dict[str, Any]:
    from pipeline import schemas

    # Keep these arguments tiny but representative.
    screen_exp_biosamples = ["MCF7", "HMEC"]
    screen_exp_assays = ["Intact-HiC", "RNAPII-ChIAPET"]
    screen_comp_biosamples = ["MCF7"]
    screen_comp_assays = ["ABC", "EPI"]
    abc_celltypes = ["MCF7_ENCODE", "MCF10A-Ji2017"]
    hichip_celltypes = ["MCF7", "HMEC"]

    snap: Dict[str, Any] = {}
    snap["screen_block"] = _shape(
        schemas.empty_screen_block(screen_exp_biosamples, screen_exp_assays)
    )
    snap["abc_celltype_entry"] = _shape(schemas.empty_abc_celltype_entry())
    snap["abc_enhancer_entry"] = _shape(schemas.empty_abc_enhancer_entry(abc_celltypes))
    snap["hichip_loop_entry"] = _shape(schemas.empty_hichip_loop_entry())
    snap["hichip_celltype_entry"] = _shape(schemas.empty_hichip_celltype_entry())
    snap["hichip_block"] = _shape(schemas.empty_hichip_block(hichip_celltypes))
    snap["gene_link_entry"] = _shape(
        schemas.empty_gene_link_entry(
            screen_exp_biosamples=screen_exp_biosamples,
            screen_exp_assays=screen_exp_assays,
            screen_comp_biosamples=screen_comp_biosamples,
            screen_comp_assays=screen_comp_assays,
            abc_celltypes=abc_celltypes,
            hichip_celltypes=hichip_celltypes,
        )
    )
    snap["cell_line_signals"] = _shape(schemas.empty_cell_line_signals(["H3K27ac", "DNase"]))
    snap["mirna_entry"] = _shape(schemas.empty_mirna_entry())
    snap["atac_body_overlap_entry"] = _shape(schemas.empty_atac_body_overlap_entry())
    return snap


def main() -> None:
    ap = argparse.ArgumentParser(description="Dump a stable schema snapshot for key nested structures.")
    ap.add_argument(
        "--out",
        type=str,
        default="tests/fixtures/pipeline_schema_snapshot.json",
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
            print("[OK] schema snapshot already up to date")
            return
        if args.check:
            raise SystemExit(f"[FAIL] schema snapshot differs: {out_path}")

    out_path.write_text(rendered, encoding="utf-8")
    print(f"[OK] wrote: {out_path}")


if __name__ == "__main__":
    main()

