"""ADVERSARIAL AUDIT of the block structure. Six ways a 'block' can be wrong, each tested mechanically."""
import re
import sys
from collections import Counter, defaultdict

import pandas as pd

sys.path.insert(0, "/sci/labs/michall/stavzok/APM")
from mirna_hallmark.analyses.ops import gen_column_reference as G  # noqa: E402

g = G._load()
print(f"{len(g)} card-columns · {g.groupby(['card','block']).ngroups} (card, block) groups\n")

# ─────────────────────────────────────────────────────────────── A. truncated multi-token prefixes
print("A. TRUNCATED PREFIX — every column in the block shares a LONGER common prefix,")
print("   so the block name is one token short of the real convention.")
hits = []
for (card, blk), sub in g.groupby(["card", "block"]):
    if blk == "(bare)" or len(sub) < 2:
        continue
    rests = [c[len(blk):] for c in sub.column]
    first = {r.split("_")[0] for r in rests if "_" in r}
    if len(first) == 1 and all("_" in r for r in rests):
        tok = first.pop()
        hits.append((card, blk, f"{blk}{tok}_", len(sub)))
for card, blk, real, n in sorted(hits):
    print(f"   ⛔ {card:<13}{blk:<12}-> every one of its {n} columns is really `{real}`")
print(f"   {len(hits)} block(s)\n" if hits else "   ✅ none\n")

# ─────────────────────────────────────────────────────────────── B. block spans multiple rungs
print("B. RUNG-INCOHERENT — one block whose columns live on DIFFERENT units.")
bad = []
for (card, blk), sub in g.groupby(["card", "block"]):
    rungs = {r for r in sub.rung if r}
    if len(rungs) > 1:
        bad.append((card, blk, len(sub), sorted(rungs)))
for card, blk, n, rungs in sorted(bad, key=lambda x: -x[2])[:14]:
    print(f"   ⚠ {card:<13}{blk:<14}{n:>3} cols · rungs {rungs}")
print(f"   {len(bad)} block(s) span >1 rung\n")

# ─────────────────────────────────────────────────────────────── C. bimodal coverage within a block
print("C. COVERAGE-INCOHERENT — a block whose columns are populated on very different row sets,")
print("   i.e. probably two things sharing a prefix.")
bad = []
for (card, blk), sub in g.groupby(["card", "block"]):
    f = sub.fill.dropna()
    if len(f) < 3:
        continue
    span = f.max() - f.min()
    if span > 0.45:
        bad.append((card, blk, len(sub), f.min(), f.max(), span))
for card, blk, n, lo, hi, span in sorted(bad, key=lambda x: -x[5])[:12]:
    print(f"   ⚠ {card:<13}{blk:<14}{n:>3} cols · fill {lo:.0%} – {hi:.0%}  (span {span:.0%})")
print(f"   {len(bad)} block(s) with >45pt coverage span\n")

# ─────────────────────────────────────────────────────────────── D. singleton blocks
print("D. SINGLETON BLOCKS — a prefix used by exactly one column is not a block.")
singles = [(c, b) for (c, b), s in g.groupby(["card", "block"]) if len(s) == 1]
by_card = Counter(c for c, _ in singles)
print(f"   {len(singles)} singleton block(s): " + " · ".join(f"{k} {v}" for k, v in by_card.items()))
orphan_pref = [(c, b) for c, b in singles
               if b != "(bare)" and not any(b == b2 and c != c2 for (c2, b2), _ in g.groupby(["card", "block"]))]
print(f"   of those, prefixes that exist on NO other card either: {len(orphan_pref)}")
print("   " + ", ".join(f"{c}.{b}" for c, b in sorted(orphan_pref)[:16]) + "\n")

# ─────────────────────────────────────────────────────────────── E. same prefix, different meaning per card
print("E. CROSS-CARD PREFIX COLLISION — same block name, DISJOINT column sets ⇒ two different things.")
byblk = defaultdict(dict)
for (card, blk), sub in g.groupby(["card", "block"]):
    byblk[blk][card] = set(sub.column)
for blk, per in sorted(byblk.items()):
    if len(per) < 2:
        continue
    cards = sorted(per)
    for i in range(len(cards)):
        for j in range(i + 1, len(cards)):
            a, b = per[cards[i]], per[cards[j]]
            inter = len(a & b)
            if inter == 0 and min(len(a), len(b)) >= 2:
                print(f"   ⛔ {blk:<14}{cards[i]}({len(a)}) vs {cards[j]}({len(b)}) — 0 shared names")
print()

# ─────────────────────────────────────────────────────────────── F. concept split across blocks
print("F. SPLIT CONCEPT — the same trailing token appearing under several prefixes on one card.")
for card in ("edge", "arm", "gene"):
    sub = g[g.card == card]
    tail = defaultdict(set)
    for r in sub.itertuples():
        rest = r.column[len(r.block):] if r.block != "(bare)" else r.column
        if rest:
            tail[rest].add(r.block)
    multi = {t: bs for t, bs in tail.items() if len(bs) > 2}
    if multi:
        print(f"   {card}:")
        for t, bs in sorted(multi.items(), key=lambda kv: -len(kv[1]))[:6]:
            print(f"      `…{t}` appears under {len(bs)} prefixes: {sorted(bs)}")
