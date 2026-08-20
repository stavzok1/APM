"""Browse the card system: five rungs, ~700 columns, and what you must not do with them.

    python -m mirna_hallmark.learned.cards                  overview of all five cards
    python -m mirna_hallmark.learned.cards arm              one card, grouped into blocks
    python -m mirna_hallmark.learned.cards arm --block real_   one block, with per-column stats
    python -m mirna_hallmark.learned.cards --col beta       where a column lives, on every card
    python -m mirna_hallmark.learned.cards --find promisc   search names and descriptions
    python -m mirna_hallmark.learned.cards --warn           every caveat attached to a live column

Reads `card_registry.tsv` (rung + agg_of + domain, written by `card_rungs --check`) and the card
files themselves (coverage, dtype, examples). Read-only.

WHY THE `family` / `seed_family` DISTINCTION KEEPS BITING — this tool prints the KEY, always:

    CARDS["gene_family"]  key = [gene, family]     -> a GENE-BY-FAMILY card. One row per (gene, family):
                                                      "how family F behaves AT gene G".
    CARDS["seed_family"]  key = [seed_family]      -> the FAMILY rung. One row per family:
                                                      "what family F IS", independent of any gene.

The name `family` predates the fifth rung and is a misnomer; `gene_family` would be accurate. Renaming
the FILE touches 13 modules, so the honest interim fix is this: the key is shown next to the name every
time, because the key is what actually disambiguates them.
"""
from __future__ import annotations

import argparse
import re
import shutil

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.learned.card_rungs import CARDS

REG = C.REPO_ROOT / "mirna_hallmark/output/learned/card_registry.tsv"
W = min(shutil.get_terminal_size((120, 40)).columns, 140)
RULE = "─" * W


def _reg() -> pd.DataFrame:
    if not REG.exists():
        raise SystemExit("no card_registry.tsv — run: python -m mirna_hallmark.learned.card_rungs --check")
    r = pd.read_csv(REG, sep="\t")
    for c in ("agg_of", "domain"):
        if c not in r.columns:
            r[c] = ""
    return r.fillna({"agg_of": "", "domain": ""})


def _load(card: str) -> pd.DataFrame:
    return pd.read_csv(CARDS[card]["path"], sep="\t", low_memory=False)


def _block_of(col: str) -> str:
    """The block a column belongs to = its prefix, or '(bare)' for the un-prefixed lifted columns.

    ⛔ CASE-INSENSITIVE (fixed 2026-08-19, user-caught): a camelCase prefix such as `dGlobal_` / `dShare_`
    does not match `[a-z]+_`, so those columns were filed as un-prefixed. The `bc_|hx_|cov_` alternation
    stays — those are genuine two-token prefixes that would otherwise truncate to their first token.
    """
    m = re.match(r"^((?:bc_|hx_|cov_)?[A-Za-z]+_)", col)
    return m.group(1) if m else "(bare)"


def _wrap(s: str, indent: int, width: int = None) -> str:
    import textwrap
    return textwrap.fill(str(s), width or (W - indent), initial_indent=" " * indent,
                         subsequent_indent=" " * indent)


def _short(s: str, n: int) -> str:
    s = str(s)
    return s if len(s) <= n else s[: n - 1] + "…"


# ── views ─────────────────────────────────────────────────────────────────
def overview() -> None:
    r = _reg()
    print(f"\n{RULE}\nCARD SYSTEM — 5 rungs, {len(r):,} card-columns\n{RULE}")
    print(f"{'card':14s} {'key':24s} {'rows':>7s} {'cols':>6s}  blocks")
    for name, spec in CARDS.items():
        if not spec["path"].exists():
            print(f"{name:14s} {'[' + ', '.join(spec['key']) + ']':24s} {'—':>7s} {'—':>6s}  (not built)")
            continue
        d = _load(name)
        blocks = sorted({_block_of(c) for c in d.columns})
        print(f"{name:14s} {'[' + ', '.join(spec['key']) + ']':24s} {len(d):>7,d} {d.shape[1]:>6d}  "
              f"{_short(' '.join(blocks), W - 56)}")
    print(f"\n{'NOTE':14s} `family` is keyed [gene, family] — a GENE-BY-FAMILY card, not the family rung.")
    print(f"{'':14s} `seed_family` is the family rung, one row per family.")
    print(f"\nrung distribution across all cards:")
    for rung, n in r.rung.value_counts().items():
        print(f"   {rung:16s} {n:4d}")
    agg = r[r.agg_of != ""]
    if len(agg):
        print(f"\n{len(agg)} column(s) SUMMARISE a lower unit (agg_of): "
              + ", ".join(f"{k}×{v}" for k, v in agg.agg_of.value_counts().items()))
    print(f"\ntry:  cards <name>   |   cards --col <column>   |   cards --find <text>   |   cards --warn\n")


def show_card(card: str, block: str | None = None) -> None:
    r = _reg()
    r = r[r.card == card]
    d = _load(card)
    spec = CARDS[card]
    key = spec["key"]
    print(f"\n{RULE}\n{card.upper()} CARD   key=[{', '.join(key)}]   {len(d):,} rows × {d.shape[1]} cols")
    print(f"{spec['path']}\n{RULE}")

    cols = [c for c in d.columns if not block or c.startswith(block)]
    if block and not cols:
        print(f"no columns with prefix '{block}'"); return

    groups: dict[str, list[str]] = {}
    for c in cols:
        groups.setdefault(_block_of(c), []).append(c)

    dom = dict(zip(r.column, r.domain))
    rung = dict(zip(r.column, r.rung))
    agg = dict(zip(r.column, r.agg_of))

    for b in sorted(groups, key=lambda k: (-len(groups[k]), k)):
        cs = groups[b]
        cov = d[cs].notna().any(axis=1).mean()
        rungs = sorted({rung.get(c, "?") for c in cs})
        print(f"\n  {b:<16s} {len(cs):>3d} cols   coverage {100*cov:5.1f}%   rung {'/'.join(rungs)}")
        dd = {dom.get(c, "") for c in cs} - {""}
        for x in sorted(dd):
            print(_wrap(f"· {x}", 4))
        if block:                                    # detail mode: one line per column
            print()
            for c in cs:
                s = d[c]
                a = f" agg_of={agg[c]}" if agg.get(c) else ""
                stat = ""
                if pd.api.types.is_numeric_dtype(s) and s.notna().any():
                    stat = f"med {s.median():>10.4g}   [{s.min():.4g}, {s.max():.4g}]"
                elif s.notna().any():
                    top = s.value_counts().head(3)
                    stat = "top " + ", ".join(f"{k}={v}" for k, v in top.items())
                print(f"    {c:<32s} {100*s.notna().mean():5.1f}%  {_short(stat, W - 48)}{a}")
        else:
            print(_wrap(", ".join(cs), 4))
    if not block:
        print(f"\n  (use --block <prefix> for per-column detail)\n")


def show_column(col: str) -> None:
    r = _reg()
    hit = r[r.column == col]
    if hit.empty:
        near = sorted({c for c in r.column if col.lower() in c.lower()})[:12]
        print(f"\nno column named '{col}'." + (f" did you mean: {', '.join(near)}" if near else ""))
        return
    print(f"\n{RULE}\nCOLUMN  {col}\n{RULE}")
    if hit.rung.nunique() > 1:
        print("  ⚠ SAME NAME, DIFFERENT RUNG ACROSS CARDS — never join these without re-checking the unit.\n")
    for _, h in hit.iterrows():
        d = _load(h.card)
        s = d[col] if col in d.columns else None
        print(f"  on {h.card.upper():12s} key=[{', '.join(CARDS[h.card]['key'])}]  rung={h.rung}"
              + (f"  agg_of={h.agg_of}" if h.agg_of else ""))
        if s is not None:
            line = f"coverage {100*s.notna().mean():.1f}% of {len(s):,} rows   dtype {s.dtype}"
            if pd.api.types.is_numeric_dtype(s) and s.notna().any():
                line += f"   median {s.median():.4g}   range [{s.min():.4g}, {s.max():.4g}]"
            print(_wrap(line, 6))
            if not pd.api.types.is_numeric_dtype(s) and s.notna().any():
                print(_wrap("values: " + ", ".join(map(str, s.value_counts().head(5).index)), 6))
        if h.domain:
            print(_wrap(f"domain: {h.domain}", 6))
        print()


def find(text: str) -> None:
    r = _reg()
    t = text.lower()
    m = r[r.column.str.lower().str.contains(t, regex=False)
          | r.domain.str.lower().str.contains(t, regex=False)]
    if m.empty:
        print(f"\nnothing matches '{text}'"); return
    print(f"\n{RULE}\nMATCHES for '{text}' — {len(m)} card-column(s)\n{RULE}")
    print(f"{'column':<34s} {'card':<13s} {'rung':<14s} domain")
    for _, h in m.sort_values(["column", "card"]).iterrows():
        print(f"{_short(h.column, 33):<34s} {h.card:<13s} {h['rung']:<14s} {_short(h.domain, W - 64)}")
    print()


def warnings() -> None:
    """Every caveat on a live column, deduplicated — the accumulated 'do not' list.

    ⭐ READS BOTH SOURCES (2026-08-19). It used to read only `domain`, which is a ROW-APPLICABILITY
    statement — so the caveats that actually matter, the ones written into the per-column DESCRIPTIONS
    (`card_glossary`), never surfaced here. A warning nobody can list is a warning nobody reads:
    *"a coupling carried by a few spiking samples is fragile"* lived only in a docstring until this.
    """
    r = _reg()
    seen: dict[str, list] = {}
    for _, h in r[r.domain.str.contains("⛔|⚠", regex=True, na=False)].iterrows():
        seen.setdefault(h.domain, []).append((h.card, h.column, "domain"))
    try:
        from mirna_hallmark.learned import card_glossary as GL
        for _, h in r.iterrows():
            desc = GL.describe(h.column, h.card) or ""
            if "⛔" in desc or "⚠" in desc:
                seen.setdefault(desc, []).append((h.card, h.column, "description"))
    except Exception as ex:                                    # never let the glossary break the listing
        print(f"  ⚠ glossary caveats unavailable ({ex})")
    n_cols = len({(c, col) for v in seen.values() for c, col, _ in v})
    print(f"\n{RULE}\nCAVEATS ON LIVE COLUMNS — {len(seen)} distinct, over {n_cols} card-columns\n{RULE}")
    for dom, cols in sorted(seen.items(), key=lambda kv: -len(kv[1])):
        cards = sorted({c for c, _, _ in cols})
        src = sorted({s for _, _, s in cols})
        print(f"\n[{len(cols):>3d} cols on {', '.join(cards)}  · from {'+'.join(src)}]")
        print(_wrap(dom, 2))
    print()


def main() -> None:
    ap = argparse.ArgumentParser(description="browse the mirna_hallmark card system")
    ap.add_argument("card", nargs="?", choices=list(CARDS), help="show one card")
    ap.add_argument("--block", help="restrict to a column-prefix block, with per-column detail")
    ap.add_argument("--col", help="explain one column, on every card it appears on")
    ap.add_argument("--find", help="search column names and domain descriptions")
    ap.add_argument("--warn", action="store_true", help="list every caveat on a live column")
    a = ap.parse_args()
    if a.col:
        show_column(a.col)
    elif a.find:
        find(a.find)
    elif a.warn:
        warnings()
    elif a.card:
        show_card(a.card, a.block)
    else:
        overview()


if __name__ == "__main__":
    main()
