from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _dedupe_preserve_order(xs: Iterable[str]) -> List[str]:
    out: List[str] = []
    seen = set()
    for x in xs:
        x = str(x)
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _tick_list(symbols: List[str]) -> str:
    return ", ".join(f"`{s}`" for s in symbols)


def _replace_generated_block(doc: str, *, key: str, new_body: str) -> str:
    begin = f"<!-- GENERATED:{key} BEGIN -->"
    end = f"<!-- GENERATED:{key} END -->"
    if begin not in doc or end not in doc:
        raise ValueError(f"Missing generated markers for {key!r}")
    pattern = re.compile(
        re.escape(begin) + r"[\s\S]*?" + re.escape(end),
        flags=re.MULTILINE,
    )
    replacement = begin + "\n" + new_body.rstrip() + "\n" + end
    out, n = pattern.subn(replacement, doc, count=1)
    if n != 1:
        raise ValueError(f"Expected to replace exactly one block for {key!r}, replaced {n}")
    return out


def _counts_table_rows() -> List[Tuple[str, str]]:
    from pipeline import config as cfg

    primary = list(cfg.PRIMARY_GENES)
    extended = list(cfg.EXTENDED_PRIMARY_GENES)
    tier1_lncrna = list(cfg.TIER1_LNCRNA_GENES)
    full_integration = list(cfg.FULL_INTEGRATION_GENES)
    tier2 = list(cfg.TIER2_MEDIUM_GENES)
    tier3 = list(cfg.TIER3_CNV_ONLY_GENES)
    strat = list(cfg.TIER_SNV_CNV_STRATIFIER_GENES)
    tier4 = list(cfg.TIER4_READOUT_GENES)
    cnv = list(cfg.CNV_GENES)
    alias = list(cfg.RNA_EXPRESSION_ALIAS_SEED_SYMBOLS)

    def _extended_note() -> str:
        n_lncrna = len(tier1_lncrna)
        n_protein = len(extended) - n_lncrna
        return f"{len(extended)} ({n_protein} protein-coding + **{n_lncrna}** lncRNAs)"

    rows: List[Tuple[str, str]] = [
        ("`PRIMARY_GENES`", str(len(primary))),
        ("`EXTENDED_PRIMARY_GENES`", _extended_note()),
        ("`TIER1_LNCRNA_GENES`", f"{len(tier1_lncrna)} (subset of the lncRNAs above)"),
        ("`FULL_INTEGRATION_GENES`", f"{len(full_integration)} (= deduped union of primary + extended)"),
        ("`TIER2_MEDIUM_GENES`", str(len(tier2))),
        ("`TIER3_CNV_ONLY_GENES`", str(len(tier3))),
        ("`TIER_SNV_CNV_STRATIFIER_GENES`", str(len(strat))),
        ("`TIER4_READOUT_GENES`", str(len(tier4))),
        ("`CNV_GENES`", str(len(cnv))),
        ("`RNA_EXPRESSION_ALIAS_SEED_SYMBOLS`", str(len(alias))),
    ]
    return rows


def build_blocks() -> Dict[str, str]:
    from pipeline import config as cfg

    blocks: Dict[str, str] = {}

    # Counts table
    rows = _counts_table_rows()
    lines = ["| List | Count |", "|------|------:|"]
    for k, v in rows:
        lines.append(f"| {k} | {v} |")
    blocks["panel_counts"] = "\n".join(lines) + "\n"

    blocks["extended_primary_config_order"] = _tick_list(list(cfg.EXTENDED_PRIMARY_GENES)) + "\n"
    blocks["tier2_medium_genes_config_order"] = _tick_list(list(cfg.TIER2_MEDIUM_GENES)) + "\n"
    blocks["tier3_cnv_only_genes_config_order"] = _tick_list(list(cfg.TIER3_CNV_ONLY_GENES)) + "\n"
    blocks["stratifier_genes_config_order"] = _tick_list(list(cfg.TIER_SNV_CNV_STRATIFIER_GENES)) + "\n"
    blocks["tier4_readout_genes_config_order"] = _tick_list(list(cfg.TIER4_READOUT_GENES)) + "\n"

    return blocks


def main() -> None:
    ap = argparse.ArgumentParser(description="Regenerate generated blocks in research_plan/01_gene_panel_extended.md")
    ap.add_argument(
        "--doc",
        type=str,
        default="research_plan/01_gene_panel_extended.md",
        help="Target markdown doc to patch in-place.",
    )
    ap.add_argument("--check", action="store_true", help="Fail if doc differs; do not write.")
    args = ap.parse_args()

    doc_path = (REPO / args.doc).resolve()
    if not doc_path.exists():
        raise FileNotFoundError(doc_path)

    original = doc_path.read_text(encoding="utf-8")
    out = original
    blocks = build_blocks()
    for key, body in blocks.items():
        out = _replace_generated_block(out, key=key, new_body=body)

    if out == original:
        print("[OK] doc already up to date")
        return

    if args.check:
        raise SystemExit("[FAIL] generated blocks are stale; rerun without --check to rewrite")

    doc_path.write_text(out, encoding="utf-8")
    print(f"[OK] wrote: {doc_path}")


if __name__ == "__main__":
    main()

