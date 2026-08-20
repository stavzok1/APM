"""⭐ GENERATE THE COLUMN REFERENCE — every card column, organised by CARD then BLOCK, print/PDF-ready.

> **User-directed 2026-08-19: *"we need a proper document that describes the columns. maybe even pdf. it
> needs to be organized into the card type and inner blocks and such."***

The glossary answers *"what is this column?"* one lookup at a time. It does not answer *"what is on the arm
card, and how is it grouped?"* — the flat 711-row TSV is a lookup table, not a document, and 14% of its
entries exceed 600 characters, which reads badly in a spreadsheet cell.

This builds the reading version: **one section per card, blocks within it, and every column with its rung,
coverage and meaning.** It is GENERATED, never hand-edited — the glossary and the registry stay the single
sources, so the document cannot drift from them.

⚠ **Print/PDF is a first-class target here, not an afterthought.** `@media print` gives each card a page
break, drops the nav and the interactive controls, and switches to a paper palette — so
*Print → Save as PDF* yields a usable offline reference rather than a screenshot of a dark web page.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_column_reference
"""
from __future__ import annotations

import datetime as dt
import html
import pathlib
import re
from collections import defaultdict

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "learned"
DEST = ROOT / "docs" / "derived" / "COLUMN_REFERENCE.html"

CARDS = [
    ("edge", "realization/edge_card.tsv", "one (miRNA arm, gene) pair",
     "The widest card and the one most claims are made from. ⚠ 84 of its columns come from the BASE card "
     "and 102 from the annotation layer — read the vintage banner before trusting the fitted ones."),
    ("gene", "realization/gene_card.tsv", "the gene's total incoming regulation",
     "One row per gene. Nothing here is checkable by invariance — there is only one row per unit — so "
     "`agg_of` carries the meaning instead."),
    ("gene_family", "gene_family_card.tsv", "one (seed family, gene) pair",
     "⚠ Keyed [gene, family] — it CANNOT express a property of the family itself. That is what the "
     "seed_family card is for."),
    ("arm", "arm_card.tsv", "the miRNA arm itself, gene-free",
     "The deepest card. An arm-rung column must be constant within arm ACROSS genes; if it varies by gene "
     "it is not an arm property, whatever its name says."),
    ("seed_family", "seed_family_card.tsv", "the seed family itself, gene-free",
     "⚠ 83.7% of families are SINGLETONS — for those, every 'family-level' statistic is trivially that one "
     "arm's value. Mask on `fam_single_member` before any family comparison."),
]
RUNG_NOTE = {
    "edge": "one (miRNA, gene) pair", "gene": "the gene", "family": "the seed family",
    "arm": "the arm, across genes", "arm-in-family": "the arm inside one family cell",
    "key": "an identifier, not a measurement", "seed_family": "the family itself",
}


def _block(col: str) -> str:
    m = re.match(r"^([a-z]+_)", col)
    return m.group(1) if m else "(bare)"


def _md(t: str) -> str:
    """The glossary writes light markdown; render the bits that carry meaning, escape the rest."""
    t = html.escape(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t, flags=re.S)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    return t


def _load() -> pd.DataFrame:
    g = pd.read_csv(OUT / "card_glossary.tsv", sep="\t", dtype=str).fillna("")
    fills = {}
    for card, rel, _, _ in CARDS:
        p = OUT / rel
        if p.exists():
            d = pd.read_csv(p, sep="\t", low_memory=False)
            fills[card] = (d.notna().mean(), len(d))
    g["fill"] = [round(float(fills.get(c, ({}, 0))[0].get(col, float("nan")) or 0), 4)
                 if c in fills and col in fills[c][0] else None
                 for c, col in zip(g.card, g.column)]
    g["nrows"] = [fills.get(c, ({}, 0))[1] for c in g.card]
    g["block"] = g.column.map(_block)
    return g


def _stamp_banner() -> str:
    try:
        from mirna_hallmark.learned import card_stamp as ST
        mixed, gap = ST.mixed_vintage()
        if mixed:
            return (f'<div class="call stop"><p class="h">Vintage mix &mdash; read this before quoting a '
                    f'fitted column</p><p>The delivered edge card is <strong>{gap:.0f} days newer</strong> '
                    f'than the base it layers over. Its <strong>84 base-owned columns</strong> '
                    f'(<code>coupling_</code>, <code>beta_</code>, <code>dose_</code>, <code>share_</code>, '
                    f'<code>rank_</code>) carry the older vintage; the other 102 are current. Re-run '
                    f'<code>canonical_card.build()</code> before trusting them.</p></div>')
    except Exception:
        pass
    return ""


CSS = """
:root{
  --paper:#F3F5F6; --panel:#FFFFFF; --ink:#15191C; --ink-2:#414D55; --ink-3:#6E7C86;
  --rule:#D3DBDF; --rule-2:#E4EAED; --accent:#0D6A6B; --accent-soft:#E2EFEF;
  --edge:#2C5E8A; --gene:#25683F; --fam:#7A4E8F; --arm:#8A5C10; --seed:#9E332C;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --cond:"IBM Plex Sans Condensed","Helvetica Neue",Arial,sans-serif;
  --serif:"IBM Plex Serif",Georgia,"Times New Roman",serif;
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){
  --paper:#0E1214; --panel:#151A1D; --ink:#E6ECEF; --ink-2:#AFBCC4; --ink-3:#7E8C95;
  --rule:#2A3438; --rule-2:#1F282C; --accent:#48B3AE; --accent-soft:#152B2C;
  --edge:#6FA6D6; --gene:#6FBE8C; --fam:#B98FD0; --arm:#D6A94E; --seed:#E08078;
}}
:root[data-theme="dark"]{
  --paper:#0E1214; --panel:#151A1D; --ink:#E6ECEF; --ink-2:#AFBCC4; --ink-3:#7E8C95;
  --rule:#2A3438; --rule-2:#1F282C; --accent:#48B3AE; --accent-soft:#152B2C;
  --edge:#6FA6D6; --gene:#6FBE8C; --fam:#B98FD0; --arm:#D6A94E; --seed:#E08078;
}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--serif);
     font-size:16px;line-height:1.6;-webkit-font-smoothing:antialiased}
.wrap{max-width:1080px;margin:0 auto;padding:0 26px 110px}
code{font-family:var(--mono);font-size:.86em;background:var(--accent-soft);
     padding:.08em .32em;border-radius:2px;color:var(--ink)}
h1{font-family:var(--cond);font-size:38px;line-height:1.1;margin:44px 0 6px;text-wrap:balance;
   letter-spacing:-.01em}
h2{font-family:var(--cond);font-size:25px;margin:0 0 4px;letter-spacing:-.005em}
h3{font-family:var(--mono);font-size:13px;text-transform:uppercase;letter-spacing:.1em;
   color:var(--ink-3);margin:30px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--rule)}
.sub{color:var(--ink-2);max-width:70ch}
.note{color:var(--ink-3);font-size:14px;max-width:74ch}
.card-head{border-left:4px solid var(--k);background:var(--panel);border:1px solid var(--rule);
           border-left-width:4px;border-radius:2px;padding:16px 18px;margin:0 0 6px}
.card-head .key{font-family:var(--mono);font-size:12px;color:var(--k);letter-spacing:.05em}
.meta{display:flex;flex-wrap:wrap;gap:18px;font-family:var(--mono);font-size:12px;
      color:var(--ink-3);margin-top:8px}
.meta b{color:var(--ink-2);font-weight:600}
.col{border-bottom:1px solid var(--rule-2);padding:11px 0;display:grid;
     grid-template-columns:minmax(210px,270px) 1fr;gap:20px;align-items:start;break-inside:avoid}
.col:last-child{border-bottom:none}
.cname{font-family:var(--mono);font-size:13.5px;color:var(--ink);word-break:break-word;font-weight:500}
.tags{margin-top:5px;display:flex;flex-wrap:wrap;gap:5px;align-items:center}
.tag{font-family:var(--mono);font-size:10px;letter-spacing:.05em;text-transform:uppercase;
     padding:2px 6px;border:1px solid var(--rule);border-radius:2px;color:var(--ink-3)}
.tag.rung{border-color:var(--k);color:var(--k)}
.bar{display:inline-block;height:5px;background:var(--rule-2);width:46px;border-radius:1px;
     vertical-align:middle;position:relative;overflow:hidden}
.bar i{position:absolute;inset:0 auto 0 0;background:var(--k);display:block}
.desc{color:var(--ink-2);font-size:14.5px;max-width:82ch}
.desc strong{color:var(--ink)}
.toc{background:var(--panel);border:1px solid var(--rule);border-radius:2px;padding:16px 18px;margin:26px 0}
.toc ul{margin:8px 0 0;padding-left:18px;columns:2;column-gap:34px}
.toc li{font-family:var(--mono);font-size:12.5px;margin:2px 0;break-inside:avoid}
a{color:var(--accent);text-decoration:none;border-bottom:1px solid transparent}
a:hover,a:focus-visible{border-bottom-color:var(--accent)}
a:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.call{border:1px solid var(--rule);border-left:3px solid var(--accent);background:var(--panel);
      padding:13px 16px;margin:18px 0;border-radius:2px}
.call.stop{border-left-color:var(--seed)}
.call .h{font-family:var(--cond);font-weight:700;margin:0 0 4px;font-size:15px}
.call p{margin:.25em 0;max-width:78ch}
@media print{
  :root{--paper:#fff;--panel:#fff;--ink:#000;--ink-2:#222;--ink-3:#555;--rule:#bbb;--rule-2:#ddd;
        --accent:#005f60;--accent-soft:#eee;
        --edge:#26506f;--gene:#1f5334;--fam:#5e3c6e;--arm:#6d4a0d;--seed:#7d2823}
  body{font-size:10.5pt;line-height:1.45}
  .wrap{max-width:none;padding:0}
  .toc{break-after:page}
  section.card{break-before:page}
  section.card:first-of-type{break-before:auto}
  h3{break-after:avoid}
  .col{grid-template-columns:minmax(160px,200px) 1fr;gap:14px;padding:7px 0}
  .cname{font-size:9pt}.desc{font-size:9.5pt;max-width:none}
  a{color:var(--ink);border:none}   /* token, not a literal — print redefines --ink to black above */
  @page{margin:15mm 14mm}
}
"""


def build() -> pathlib.Path:
    g = _load()
    n_total = len(g)
    parts = [f"<title>APM Card Column Reference</title>",
             '<link rel="preconnect" href="https://fonts.googleapis.com">',
             '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>',
             '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
             'family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans+Condensed:wght@600;700&'
             'family=IBM+Plex+Serif:ital,wght@0,400;0,600;1,400&display=swap">',
             f"<style>{CSS}</style>", '<div class="wrap">']
    parts.append(
        f"<h1>Card column reference</h1>"
        f'<p class="sub">Every column on every card in the miRNA&nbsp;Hallmark system &mdash; '
        f'<strong>{n_total}</strong> across five rungs &mdash; grouped by card, then by block. '
        f'Each entry carries the unit it lives on, how much of the card it covers, and what it means.</p>'
        f'<p class="note">Generated {dt.datetime.now():%Y-%m-%d} from <code>card_glossary.tsv</code> and '
        f'<code>card_registry.tsv</code> by <code>analyses/ops/gen_column_reference.py</code>. '
        f'<strong>Do not hand-edit</strong> &mdash; fix the glossary and regenerate, so this cannot drift '
        f'from the source. Print or Save-as-PDF is supported: each card starts a new page.</p>')
    parts.append(_stamp_banner())
    parts.append(
        '<div class="call"><p class="h">The one rule that governs every entry below</p>'
        '<p><strong>The rung is a property of (CARD, COLUMN) &mdash; never of the column name.</strong> '
        '<code>beta</code> is edge-rung on the edge card and family-rung on the gene_family card: same '
        'estimator, different unit. That is why this document is organised by card first, and why the same '
        'name can appear in two sections with different text.</p></div>')

    # TOC
    toc = ['<div class="toc"><strong>Contents</strong><ul>']
    for card, rel, _, _ in CARDS:
        sub = g[g.card == card]
        if not len(sub):
            continue
        toc.append(f'<li><a href="#{card}">{card}</a> &mdash; {len(sub)} cols, '
                   f'{sub.block.nunique()} blocks</li>')
    toc.append("</ul></div>")
    parts.append("".join(toc))

    for card, rel, keydesc, headnote in CARDS:
        sub = g[g.card == card]
        if not len(sub):
            continue
        var = {"edge": "--edge", "gene": "--gene", "gene_family": "--fam",
               "arm": "--arm", "seed_family": "--seed"}[card]
        nrows = int(sub.nrows.iloc[0]) if len(sub) else 0
        parts.append(f'<section class="card" id="{card}" style="--k:var({var})">')
        parts.append(
            f'<div class="card-head"><div class="key">{html.escape(rel)}</div>'
            f'<h2>{card}</h2><p class="sub">One row per {html.escape(keydesc)}. '
            f'{_md(headnote)}</p>'
            f'<div class="meta"><span><b>{len(sub)}</b> columns</span>'
            f'<span><b>{nrows:,}</b> rows</span>'
            f'<span><b>{sub.block.nunique()}</b> blocks</span>'
            f'<span>median coverage <b>{sub.fill.dropna().median():.0%}</b></span></div></div>')

        by_block: dict[str, list] = defaultdict(list)
        for r in sub.itertuples():
            by_block[r.block].append(r)
        for blk in sorted(by_block, key=lambda b: (-len(by_block[b]), b)):
            rows = sorted(by_block[blk], key=lambda r: r.column)
            parts.append(f'<h3>{html.escape(blk)} &nbsp;&middot;&nbsp; {len(rows)} column'
                         f'{"s" if len(rows) != 1 else ""}</h3>')
            for r in rows:
                tags = []
                if r.rung:
                    tags.append(f'<span class="tag rung" title="{html.escape(RUNG_NOTE.get(r.rung, ""))}">'
                                f'{html.escape(r.rung)}</span>')
                if r.agg_of:
                    tags.append(f'<span class="tag">agg of {html.escape(r.agg_of)}</span>')
                if r.fill is not None and r.fill == r.fill:
                    pct = f"{r.fill:.0%}"
                    tags.append(f'<span class="tag" title="share of rows populated">'
                                f'<span class="bar"><i style="width:{r.fill*100:.0f}%"></i></span> {pct}</span>')
                desc = _md(r.description) if r.description else \
                    '<em style="color:var(--ink-3)">no description &mdash; this is a gap, not a default</em>'
                dom = (f'<p class="note" style="margin-top:5px">Defined on: {_md(r.domain)}</p>'
                       if r.domain else "")
                parts.append(
                    f'<div class="col"><div><div class="cname">{html.escape(r.column)}</div>'
                    f'<div class="tags">{"".join(tags)}</div></div>'
                    f'<div class="desc">{desc}{dom}</div></div>')
        parts.append("</section>")

    parts.append("</div>")
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text("\n".join(parts))
    return DEST


def main() -> int:
    p = build()
    kb = p.stat().st_size / 1024
    print(f"✅ {p.relative_to(ROOT.parent)}  ({kb:.0f} KB)")
    print("   Print / Save-as-PDF: each card starts a new page, nav and controls drop out.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
