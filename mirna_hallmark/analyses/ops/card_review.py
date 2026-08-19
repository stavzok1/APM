"""CARD REVIEW — a prioritised worksheet for reviewing every card column by hand.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.card_review               # the plan: blocks, in review order
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.card_review --block ctx_  # one block, in full
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.card_review --card arm    # one card's blocks
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.card_review --emit        # -> card_review.tsv (the worksheet)

⭐ WHY THIS EXISTS. 692 columns is not reviewable as a flat list, and reviewing them in alphabetical or
card order wastes the reviewer on inert columns before they reach load-bearing ones. This orders the work
by **LEVERAGE** — how much would it cost to be wrong about this column — and groups it into (card, block)
units, because a block shares a producer and a caveat and is therefore one decision, not fifteen.

THE LEVERAGE SCORE (all measured, none curated):
  `consumers`  how many .py modules mention the column     -> being wrong propagates
  `fill`       fraction of rows populated                  -> an empty column cannot mislead anyone
  `is_key`     the card's own grain
  `cross_card` the name appears on >1 card at DIFFERENT rungs -> the MH-179/187/188/191 failure class
  `flagged`    constant / near-empty / duplicated          -> a prune or documentation candidate
  `has_caveat` its description carries a ⛔ or ⚠           -> already known to be dangerous

⚠ **Leverage is NOT importance.** A high score means *an error here would be expensive*, which is why it
should be reviewed first. A column can be scientifically central and score low because nothing reads it yet.

THE VERDICT VOCABULARY — what a reviewed column should end up tagged as:
  KEEP        correct, described, and the description matches what the producer does
  RENAME      the name misleads (the `p_fam` class: it is the design dimension, not a p-value)
  REDESCRIBE  the description is wrong or thin; the column is fine
  GATE        a ratio/threshold that needs a denominator gate (the `top_identity` class)
  PRUNE       constant, duplicated, or superseded
  INVESTIGATE the producer needs reading before a verdict is possible
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import pandas as pd

from mirna_hallmark.learned import card_glossary as GL

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "learned"


def _consumer_counts() -> dict[str, int]:
    """How many .py modules mention each column name. One pass over the tree, not one grep per column."""
    text: list[str] = []
    for p in ROOT.rglob("*.py"):
        if "__pycache__" in p.parts:
            continue
        try:
            text.append(p.read_text(errors="ignore"))
        except OSError:
            continue
    counts: dict[str, int] = {}
    reg = pd.read_csv(OUT / "card_registry.tsv", sep="\t", dtype=str).fillna("")
    for col in reg["column"].unique():
        pat = re.compile(rf"\b{re.escape(col)}\b")
        counts[col] = sum(1 for t in text if pat.search(t))
    return counts


def _block(col: str) -> str:
    m = re.match(r"^([a-z]{2,12}?)_", col)
    return m.group(1) + "_" if m else "(bare)"


def build() -> pd.DataFrame:
    reg = pd.read_csv(OUT / "card_registry.tsv", sep="\t", dtype=str).fillna("")
    cov = json.loads((OUT / "card_coverage.json").read_text())
    fill = {(c["card"], c["col"]): c["frac"] for c in cov["columns"]}
    cons = _consumer_counts()
    rung_by_name: dict[str, set] = {}
    for r in reg.itertuples():
        rung_by_name.setdefault(r.column, set()).add(r.rung)

    rows = []
    for r in reg.itertuples():
        desc = GL.describe(r.column, r.card) or ""
        f = fill.get((r.card, r.column), float("nan"))
        cross = len(rung_by_name.get(r.column, set())) > 1
        caveat = ("⛔" in desc) or ("⚠" in desc)
        flagged = (f == f) and (f < 0.05)
        lev = (min(cons.get(r.column, 0), 12) * 3          # propagation cost, capped
               + (f if f == f else 0) * 4                   # an empty column cannot mislead
               + (6 if r.rung == "key" else 0)
               + (5 if cross else 0)
               + (3 if caveat else 0)
               + (2 if flagged else 0))
        rows.append({"card": r.card, "block": _block(r.column), "column": r.column,
                     "rung": r.rung, "agg_of": r.agg_of, "fill": round(f, 3) if f == f else None,
                     "consumers": cons.get(r.column, 0), "cross_card_rung": cross,
                     "has_caveat": caveat, "leverage": round(lev, 2),
                     "description": desc, "verdict": "", "note": ""})
    return pd.DataFrame(rows)


def plan(d: pd.DataFrame) -> pd.DataFrame:
    g = d.groupby(["card", "block"]).agg(
        cols=("column", "size"), leverage=("leverage", "sum"),
        max_lev=("leverage", "max"), med_fill=("fill", "median"),
        consumers=("consumers", "sum"), caveats=("has_caveat", "sum"),
        cross=("cross_card_rung", "sum")).reset_index()
    return g.sort_values("leverage", ascending=False)


def main() -> None:
    args = sys.argv[1:]
    d = build()

    if "--emit" in args:
        d.sort_values("leverage", ascending=False).to_csv(OUT / "card_review.tsv", sep="\t", index=False)
        print(f"wrote {OUT/'card_review.tsv'}  ({len(d)} rows, sorted by leverage; fill in `verdict` + `note`)")
        return

    if "--block" in args:
        b = args[args.index("--block") + 1]
        sub = d[d.block == b].sort_values(["card", "leverage"], ascending=[True, False])
        if sub.empty:
            print(f"no block {b!r}"); return
        for card, s in sub.groupby("card"):
            print(f"\n{'='*100}\n[{card}]  block {b}   {len(s)} columns\n{'='*100}")
            for r in s.itertuples():
                print(f"\n  {r.column}   rung={r.rung or '-'}"
                      f"{'  agg_of=' + r.agg_of if r.agg_of else ''}"
                      f"   fill={r.fill}   consumers={r.consumers}   leverage={r.leverage}")
                print(f"     {r.description or '(NO DESCRIPTION)'}")
        return

    if "--card" in args:
        c = args[args.index("--card") + 1]
        d = d[d.card == c]

    p = plan(d)
    print(f"REVIEW PLAN — {len(d)} columns in {len(p)} (card, block) units, hardest-to-be-wrong-about first\n")
    print(f"{'#':>3} {'card':<12} {'block':<14} {'cols':>5} {'lev':>7} {'fill':>6} {'cons':>5} {'cav':>4} {'xrung':>6}")
    for i, r in enumerate(p.itertuples(), 1):
        mf = f"{r.med_fill:.2f}" if r.med_fill == r.med_fill else "  -"
        print(f"{i:>3} {r.card:<12} {r.block:<14} {r.cols:>5} {r.leverage:>7.0f} {mf:>6} "
              f"{r.consumers:>5} {r.caveats:>4} {r.cross:>6}")
    top = p.head(12)
    print(f"\n⭐ The first {len(top)} units cover {top.cols.sum()} columns "
          f"({top.cols.sum()/len(d):.0%} of the review) and {top.leverage.sum()/p.leverage.sum():.0%} of total leverage.")
    print("   Review them with  --block <name>  and record verdicts in the --emit worksheet.")


if __name__ == "__main__":
    main()
