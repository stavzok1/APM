from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import pandas as pd


def _block_for_col(c: str) -> str:
    x = str(c)
    xl = x.lower()
    if x in ("TCGA Participant Barcode", "TCGA Study", "Immune Subtype", "TCGA Subtype"):
        return "id_core"
    if xl.startswith("pathologic_") or "pam50" in xl or xl.endswith("_status_nature2012"):
        return "clinical_subtype"
    if "response" in xl or x in ("Wound Healing", "Proliferation", "Macrophage Regulation", "CTA Score"):
        return "program_scores"
    if "neoantigen" in xl or "mutation rate" in xl or "nonsilent" in xl or "silent" in xl:
        return "mut_burden_neoantigen"
    if x in ("Number of Segments", "Fraction Altered", "Aneuploidy Score", "Homologous Recombination Defects"):
        return "cnv_instability"
    if xl.startswith("bcr ") or xl.startswith("tcr "):
        return "receptor_diversity"
    if "cells" in xl:
        return "cell_fractions"
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, default="annotations/Thornsson_immune_table.tsv")
    ap.add_argument("--out", type=str, default="")
    args = ap.parse_args()

    inp = Path(args.input)
    df0 = pd.read_csv(inp, sep="\t", nrows=0, low_memory=False)
    cols = list(df0.columns)
    blocks: Dict[str, List[str]] = {}
    for c in cols:
        b = _block_for_col(c)
        blocks.setdefault(b, []).append(c)

    out = (
        Path(args.out)
        if args.out
        else Path("pipeline")
        / "md"
        / "contract_tables"
        / "Thorsson_immune_table.profile.md"
    )
    out.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    lines.append("## Thorsson immune table: column inventory (grouped)\n")
    lines.append(f"- **file**: `{inp}`\n")
    lines.append(f"- **n_columns**: {len(cols)}\n")
    lines.append("\n### Blocks\n")
    for b in sorted(blocks.keys()):
        lines.append(f"#### {b} ({len(blocks[b])})\n")
        for c in blocks[b]:
            lines.append(f"- `{c}`\n")
        lines.append("\n")

    out.write_text("".join(lines), encoding="utf-8")
    print(f"[OK] wrote: {out}")


if __name__ == "__main__":
    main()

