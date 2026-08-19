"""⭐ FIND EVERY BOOLEAN THAT CANNOT SAY "UNKNOWN" — the session's most repeated defect class.

⛔ **THE FAILURE.** In Python `NaN < NaN`, `NaN >= 0.5` and `NaN > 0.5` are all **False**, never NaN. So a
flag built as a bare comparison of possibly-missing numbers is defined on EVERY row, and every unmeasurable
row silently reads as the negative class. Five separate instances were found by hand in one session:

  `echim_any` · `fst_is_dominant_{HLY,TUM}` · `arb_max_identity`'s inert gate ·
  `cptac_*_agg_beats_abund_prot` (MH-253, deflated a headline 25.1% -> 31.7% and manufactured a
  cohort difference) · `identity_reliable` (MH-255, overstated the unreliable rate 2.3x)

Each was found only because someone printed that block's coverage. This module does it for all of them.

**THE TEST.** A boolean column is SUSPECT when it is materially better-populated than the numeric columns
of its own block — because a flag *derived* from those numbers cannot honestly know more than they do. That
is a heuristic, not a proof: a flag with an independent source legitimately outruns its neighbours (an
`*_measured` coverage flag is the canonical example, and is excluded by name). ⇒ output is a REVIEW QUEUE
ranked by the size of the gap, not a verdict. Read the producer before changing anything.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.nan_bool_audit          # ranked report
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.nan_bool_audit --all    # include the excluded
"""
from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd

from mirna_hallmark.learned.card_rungs import OUT

CARDS = {"edge": "realization/edge_card.tsv", "gene": "realization/gene_card.tsv",
         "gene_family": "gene_family_card.tsv", "arm": "arm_card.tsv",
         "seed_family": "seed_family_card.tsv"}
# flags whose coverage is legitimately independent of their block's numerics
EXEMPT = re.compile(r"(^cov_|_measured$|_meas$|^has_|_has_|_present$|^is_(gold|curated)|_src$)")
GAP = 0.02          # a flag must beat its block by >2 points of fill to be queued


def _block(col: str) -> str:
    m = re.match(r"^([a-z]+_)", col)
    return m.group(1) if m else "(bare)"


def audit(show_all: bool = False) -> pd.DataFrame:
    rows = []
    for card, rel in CARDS.items():
        p = OUT / rel
        if not p.exists():
            continue
        d = pd.read_csv(p, sep="\t", low_memory=False)
        fill = d.notna().mean()
        blocks: dict[str, list[str]] = {}
        for c in d.columns:
            blocks.setdefault(_block(c), []).append(c)
        for c in d.columns:
            v = d[c].dropna()
            if v.empty or not set(map(str, v.unique())) <= {"True", "False", "1.0", "0.0", "1", "0"}:
                continue
            if EXEMPT.search(c) and not show_all:
                continue
            # the numeric siblings this flag would have to be derived from
            sibs = [s for s in blocks[_block(c)]
                    if s != c and pd.to_numeric(d[s], errors="coerce").notna().any()
                    and d[s].dropna().nunique() > 2]
            if not sibs:
                continue
            ref = float(np.median([fill[s] for s in sibs]))
            # ⭐ THE SHARP TEST, added after the fill-median heuristic produced mostly false positives:
            # is the flag defined on any row where EVERY numeric sibling is missing? A flag genuinely
            # derived from this block cannot know anything there. This needs no knowledge of the producer
            # and does not care which sibling is the real input — `identified` (from a 99.9%-filled `z`)
            # scores 0 here while sitting 33 points above the block median.
            allnan = pd.concat([pd.to_numeric(d[s], errors="coerce").isna() for s in sibs],
                               axis=1).all(axis=1)
            orphan = int((d[c].notna() & allnan).sum())
            if fill[c] - ref > GAP or orphan or show_all:
                rows.append({"card": card, "column": c, "flag_fill": round(float(fill[c]), 4),
                             "block_fill": round(ref, 4), "gap": round(float(fill[c]) - ref, 4),
                             "orphan": orphan,
                             "n_rows_beyond": int(round((fill[c] - ref) * len(d))),
                             "block": _block(c), "n_sibs": len(sibs)})
    q = pd.DataFrame(rows)
    return q if q.empty else q.sort_values(["orphan", "gap"], ascending=False)


def main() -> int:
    show_all = "--all" in sys.argv
    q = audit(show_all)
    if q.empty:
        print("✅ no boolean column outruns its block's numerics by more than "
              f"{GAP:.0%} — nothing queued.")
        return 0
    print(f"⚠ {len(q)} boolean column(s) better-populated than their block's numerics "
          f"(gap > {GAP:.0%}) — REVIEW QUEUE, not a verdict:\n")
    print(f"  {'card':<12}{'column':<38}{'flag':>7}{'block':>8}{'gap':>8}{'ORPHAN':>8}")
    for r in q.itertuples():
        mark = " ⛔" if r.orphan else ""
        print(f"  {r.card:<12}{r.column:<38}{r.flag_fill:>7.1%}{r.block_fill:>8.1%}"
              f"{r.gap:>8.1%}{r.orphan:>8d}{mark}")
    print("\n  ORPHAN = rows where the flag has a value but EVERY numeric sibling in its block is NaN.")
    print("  ⛔ ORPHAN > 0 is close to proof. A large `gap` with ORPHAN = 0 usually means the flag's real")
    print("     input is one well-populated sibling and the block MEDIAN is a poor reference — check, then")
    print("     move on.")
    print("\n⇒ for each: read the producer. If it is a bare comparison of possibly-missing numbers,\n"
          "  mask it — `(a < b).where(a.notna() & b.notna())`. If its source is genuinely independent,\n"
          "  add it to EXEMPT with a reason.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
