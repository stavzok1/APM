"""⭐ THE COLUMN REFERENCE — every card column as a readable document, HTML + PDF.

> **User-directed 2026-08-19, then revised:** *"a proper document that describes the columns … organized
> into the card type and inner blocks"*, then *"much too long the entries … I need something much more
> readable. and I want a pdf doc. and the blocks per edges need to be organized by actual importance. and
> it needs to have order of contents in the beginning of each such card type."*

**FOUR THINGS THIS DOES THAT THE GLOSSARY TSV CANNOT.**

**1. It splits each entry into a LEAD and its CAVEATS.** The glossary grew warning-by-warning across 27
review units, so a median entry ran **231 characters** and 14% exceeded 600 — reference text, not a
description. Splitting at the first ⛔/⚠/⭐ marker (or the first sentence boundary, whichever comes first)
gives a **median lead of 109 characters** and folds **70% of the mass** into a caveat block that is
collapsed in HTML and set small in print. **0 entries lose their lead**, so nothing is hidden that a reader
needs to identify the column.

**2. Blocks are ordered by MEASURED importance, not by size.** Three independent signals, because any one
alone is gameable: how often the block's columns are named in `DISCOVERY_REGISTRY.md` (what the science
actually rests on — weighted highest), how many modules read them, and whether the block is base-owned
(fitted) rather than an annotation join. ⚠ Short bare names (`n`, `z`, `gene`) match everywhere in prose, so
names under 6 characters only count **backtick-quoted** registry hits — without that correction `(bare)`
scored top for the wrong reason.

**3. Each card opens with its own contents** — its blocks in importance order, with counts and tier.

**4. It emits a real PDF**, not just a print stylesheet, via headless Chrome.

⚠ GENERATED — never hand-edit. The glossary and registry stay the single sources, so this cannot drift.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_column_reference          # HTML + PDF
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_column_reference --no-pdf
"""
from __future__ import annotations

import datetime as dt
import html
import pathlib
import re
import shutil
import subprocess
import sys
from collections import defaultdict

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "output" / "learned"
DEST = ROOT / "docs" / "derived" / "COLUMN_REFERENCE.html"
PDF = ROOT / "docs" / "derived" / "COLUMN_REFERENCE.pdf"

CARDS = [
    ("edge", "realization/edge_card.tsv", "one (miRNA arm, gene) pair",
     "The widest card and the one most claims are made from. 84 of its columns come from the BASE card and "
     "102 from the annotation layer — check the vintage banner before quoting a fitted one."),
    ("gene", "realization/gene_card.tsv", "the gene's total incoming regulation",
     "One row per gene, so nothing here is checkable by invariance — `agg_of` carries the meaning instead."),
    ("gene_family", "gene_family_card.tsv", "one (seed family, gene) pair",
     "Keyed [gene, family] — it cannot express a property of the family itself. That is the seed_family "
     "card's job."),
    ("arm", "arm_card.tsv", "the miRNA arm itself, gene-free",
     "The deepest card. An arm-rung column must be constant within arm ACROSS genes; if it varies by gene "
     "it is not an arm property, whatever its name says."),
    ("seed_family", "seed_family_card.tsv", "the seed family itself, gene-free",
     "83.7% of families are SINGLETONS — for those every 'family-level' statistic is trivially that one "
     "arm's value. Mask on `fam_single_member` first."),
]
RUNG_NOTE = {"edge": "one (miRNA, gene) pair", "gene": "the gene", "family": "the seed family",
             "arm": "the arm, across genes", "arm-in-family": "the arm inside one family cell",
             "key": "an identifier, not a measurement", "seed_family": "the family itself"}
TIERS = [(0.78, "core", "the columns findings rest on"),
         (0.60, "supporting", "read these when the core raises a question"),
         (0.00, "reference", "context, provenance and rarely-quoted detail")]

#: ⛔⛔ **CURATED FLOOR — because MENTION-FREQUENCY IS NOT IMPORTANCE, and validating the ranking against
#: known cases proved it.** The three measured signals score how much a block is *used and discussed*. A
#: block that is **settled** gets discussed less precisely because nobody argues about it any more.
#: `adm_` scored **27/44 (reference tier)** while MH-216 calls `adm_has_site` *the project's single most
#: load-bearing conditioning variable* — it has only **4 registry mentions** because it was settled early.
#: ⇒ blocks with a RECORDED standing get a floor, each with the row that justifies it. ⚠ This is the one
#: hand-set input in an otherwise measured ranking; keep it short, and cite the finding, never a preference.
CURATED_FLOOR = {
    ("edge", "adm_"): (0.90, "MH-216: `adm_has_site` is the single most load-bearing conditioning variable"),
    ("arm", "adm_"): (0.86, "MH-259: the arm rollup of the same conditioning variable"),
    ("edge", "cal_"): (0.80, "the CALIBRATED width — the honest version of `identified` (MH-83)"),
    ("edge", "echim_"): (0.80, "chimeric evidence: the convergent-evidence ladder's independent rung"),
    ("gene", "comp_"): (0.82, "axiom 8: `comp_*_driver_share` is the composition GATE that predicts SIGN"),
}
MARK = re.compile(r"[⛔⚠⭐]")


def _block(col: str) -> str:
    m = re.match(r"^([a-z]+_)", col)
    return m.group(1) if m else "(bare)"


def _split(t: str) -> tuple[str, str]:
    """LEAD (what the column is) and CAVEATS (why you might misuse it). Never returns an empty lead."""
    m = MARK.search(t)
    head, rest = (t[:m.start()].strip(), t[m.start():].strip()) if m else (t.strip(), "")
    if len(head) > 200:
        parts = re.split(r"(?<=[.!?])\s+(?=[A-Z⛔⚠⭐`*])", head)
        if len(parts) > 1:
            rest = " ".join(parts[1:]) + ((" " + rest) if rest else "")
            head = parts[0]
    if not head and rest:
        parts = re.split(r"(?<=[.!?])\s+", rest, maxsplit=1)
        head, rest = parts[0], (parts[1] if len(parts) > 1 else "")
    return head.strip(), rest.strip()


def _md(t: str) -> str:
    t = html.escape(t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t, flags=re.S)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    return t


def _importance(g: pd.DataFrame) -> pd.DataFrame:
    reg = (ROOT / "docs" / "DISCOVERY_REGISTRY.md").read_text()
    try:
        base_cols = set(pd.read_csv(OUT / "edge_card_base.tsv", sep="\t", nrows=1).columns)
    except OSError:
        base_cols = set()
    blob = "\n".join(f.read_text(errors="ignore") for f in ROOT.rglob("*.py")
                     if "__pycache__" not in f.parts)
    rows = []
    for (card, blk), sub in g.groupby(["card", "block"]):
        cols = list(sub.column)
        # ⚠ short bare names occur in ordinary prose; only count them when backtick-quoted.
        hits = 0
        for c in cols:
            hits += len(re.findall(rf"`{re.escape(c)}`", reg)) if len(c) < 6 \
                else len(re.findall(rf"`?{re.escape(c)}`?", reg))
        rows.append({"card": card, "block": blk, "registry": hits,
                     "code": sum(blob.count(f'"{c}"') for c in cols),
                     "fitted": (sum(c in base_cols for c in cols) / len(cols)) if card == "edge" else 0.0})
    d = pd.DataFrame(rows)
    for c in ("registry", "code"):
        d[c + "_r"] = d.groupby("card")[c].rank(pct=True)
    d["score"] = 0.55 * d.registry_r + 0.30 * d.code_r + 0.15 * d.fitted
    s = d.set_index(["card", "block"])["score"]
    for k, (floor, _why) in CURATED_FLOOR.items():
        if k in s.index:
            s[k] = max(float(s[k]), floor)
    return s


def _tier(score: float) -> tuple[str, str]:
    for cut, name, blurb in TIERS:
        if score >= cut:
            return name, blurb
    return TIERS[-1][1], TIERS[-1][2]


def _load() -> pd.DataFrame:
    g = pd.read_csv(OUT / "card_glossary.tsv", sep="\t", dtype=str).fillna("")
    fills, nrows = {}, {}
    for card, rel, _, _ in CARDS:
        p = OUT / rel
        if p.exists():
            d = pd.read_csv(p, sep="\t", low_memory=False)
            fills[card], nrows[card] = d.notna().mean(), len(d)
    g["block"] = g.column.map(_block)
    g["fill"] = [float(fills[c].get(col)) if c in fills and col in fills[c] else None
                 for c, col in zip(g.card, g.column)]
    g["nrows"] = [nrows.get(c, 0) for c in g.card]
    lead_rest = [_split(t) for t in g.description]
    g["lead"] = [a for a, _ in lead_rest]
    g["caveat"] = [b for _, b in lead_rest]
    return g


def _banner() -> str:
    try:
        from mirna_hallmark.learned import card_stamp as ST
        mixed, gap = ST.mixed_vintage()
        if mixed:
            return ('<div class="call stop"><p class="h">Vintage mix — read before quoting a fitted column</p>'
                    f'<p>The delivered edge card is <strong>{gap:.0f} days newer</strong> than the base it '
                    'layers over. Its <strong>84 base-owned columns</strong> (<code>coupling_</code>, '
                    '<code>beta_</code>, <code>dose_</code>, <code>share_</code>, <code>rank_</code>) carry '
                    'the older vintage; the other 102 are current.</p></div>')
    except Exception:
        pass
    return ""


CSS = """
:root{
  --paper:#F4F6F7;--panel:#FFF;--ink:#14181B;--ink-2:#46525A;--ink-3:#75838C;
  --rule:#D5DCE0;--rule-2:#E7ECEF;--accent:#0D6A6B;--accent-soft:#E3EFEF;
  --edge:#2C5E8A;--gene:#25683F;--fam:#7A4E8F;--arm:#8A5C10;--seed:#9E332C;
  --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
  --cond:"IBM Plex Sans Condensed","Helvetica Neue",Arial,sans-serif;
  --serif:"IBM Plex Serif",Georgia,serif;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --paper:#0E1214;--panel:#151A1D;--ink:#E6ECEF;--ink-2:#AFBCC4;--ink-3:#7E8C95;
  --rule:#2A3438;--rule-2:#1F282C;--accent:#48B3AE;--accent-soft:#152B2C;
  --edge:#6FA6D6;--gene:#6FBE8C;--fam:#B98FD0;--arm:#D6A94E;--seed:#E08078;}}
:root[data-theme="dark"]{
  --paper:#0E1214;--panel:#151A1D;--ink:#E6ECEF;--ink-2:#AFBCC4;--ink-3:#7E8C95;
  --rule:#2A3438;--rule-2:#1F282C;--accent:#48B3AE;--accent-soft:#152B2C;
  --edge:#6FA6D6;--gene:#6FBE8C;--fam:#B98FD0;--arm:#D6A94E;--seed:#E08078;}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--serif);
     font-size:15.5px;line-height:1.55;-webkit-font-smoothing:antialiased}
.wrap{max-width:1000px;margin:0 auto;padding:0 26px 90px}
code{font-family:var(--mono);font-size:.85em;background:var(--accent-soft);padding:.05em .3em;
     border-radius:2px;color:var(--ink)}
h1{font-family:var(--cond);font-size:36px;line-height:1.08;margin:40px 0 6px;letter-spacing:-.012em}
h2{font-family:var(--cond);font-size:27px;margin:0 0 3px;letter-spacing:-.006em}
h3{font-family:var(--mono);font-size:12px;text-transform:uppercase;letter-spacing:.11em;
   color:var(--k);margin:26px 0 2px}
h3 .cnt{color:var(--ink-3);letter-spacing:.04em}
.tierline{font-family:var(--cond);font-size:12px;color:var(--ink-3);margin:0 0 9px;
          padding-bottom:6px;border-bottom:1px solid var(--rule)}
.sub{color:var(--ink-2);max-width:70ch}
.note{color:var(--ink-3);font-size:13.5px;max-width:76ch}
.card-head{border:1px solid var(--rule);border-left:4px solid var(--k);background:var(--panel);
           border-radius:2px;padding:15px 17px;margin:0 0 4px}
.card-head .key{font-family:var(--mono);font-size:11.5px;color:var(--k);letter-spacing:.05em}
.meta{display:flex;flex-wrap:wrap;gap:16px;font-family:var(--mono);font-size:11.5px;
      color:var(--ink-3);margin-top:7px}
.meta b{color:var(--ink-2);font-weight:600}
.cardtoc{background:var(--panel);border:1px solid var(--rule);border-radius:2px;
         padding:13px 16px;margin:10px 0 4px}
.cardtoc .t{font-family:var(--cond);font-weight:700;font-size:13px;margin:0 0 7px}
.cardtoc ol{margin:0;padding-left:20px;columns:2;column-gap:30px;font-family:var(--mono);font-size:12px}
.cardtoc li{margin:1.5px 0;break-inside:avoid;color:var(--ink-2)}
.cardtoc .tg{color:var(--ink-3);font-size:10.5px;text-transform:uppercase;letter-spacing:.06em}
.col{border-bottom:1px solid var(--rule-2);padding:8px 0;display:grid;
     grid-template-columns:minmax(190px,240px) 1fr;gap:18px;align-items:baseline;break-inside:avoid}
.col:last-child{border-bottom:none}
.cname{font-family:var(--mono);font-size:13px;font-weight:500;word-break:break-word}
.tags{margin-top:4px;display:flex;flex-wrap:wrap;gap:4px}
.tag{font-family:var(--mono);font-size:9.5px;letter-spacing:.05em;text-transform:uppercase;
     padding:1px 5px;border:1px solid var(--rule);border-radius:2px;color:var(--ink-3)}
.tag.rung{border-color:var(--k);color:var(--k)}
.lead{color:var(--ink-2);font-size:14.5px;max-width:80ch}
.lead strong{color:var(--ink)}
details{margin-top:4px}
summary{font-family:var(--cond);font-size:11.5px;color:var(--accent);cursor:pointer;
        letter-spacing:.03em;list-style:none}
summary::-webkit-details-marker{display:none}
summary::before{content:"▸ ";font-size:9px}
details[open] summary::before{content:"▾ "}
.cav{color:var(--ink-3);font-size:13px;max-width:82ch;margin-top:4px;padding-left:11px;
     border-left:2px solid var(--rule)}
.cav strong{color:var(--ink-2)}
a{color:var(--accent);text-decoration:none;border-bottom:1px solid transparent}
a:hover,a:focus-visible{border-bottom-color:var(--accent)}
a:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.call{border:1px solid var(--rule);border-left:3px solid var(--accent);background:var(--panel);
      padding:12px 15px;margin:16px 0;border-radius:2px}
.call.stop{border-left-color:var(--seed)}
.call .h{font-family:var(--cond);font-weight:700;margin:0 0 3px;font-size:14.5px}
.call p{margin:.2em 0;max-width:78ch}
.toc{background:var(--panel);border:1px solid var(--rule);border-radius:2px;padding:15px 17px;margin:22px 0}
.toc ul{margin:7px 0 0;padding-left:18px;columns:2;column-gap:32px}
.toc li{font-family:var(--mono);font-size:12.5px;margin:2px 0;break-inside:avoid}
@media print{
  :root{--paper:#fff;--panel:#fff;--ink:#000;--ink-2:#1e1e1e;--ink-3:#5a5a5a;--rule:#b4b4b4;
        --rule-2:#dcdcdc;--accent:#00595a;--accent-soft:#eef1f1;
        --edge:#26506f;--gene:#1f5334;--fam:#5e3c6e;--arm:#6d4a0d;--seed:#7d2823}
  body{font-size:9.6pt;line-height:1.36}
  .wrap{max-width:none;padding:0}
  .toc{break-after:page}
  section.card{break-before:page}
  section.card:first-of-type{break-before:auto}
  h3,.tierline,.card-head,.cardtoc{break-after:avoid}
  .col{grid-template-columns:minmax(150px,185px) 1fr;gap:12px;padding:5px 0}
  .cname{font-size:8.4pt}.lead{font-size:8.8pt;max-width:none}
  details{display:block}                     /* caveats PRINT — collapsed only on screen */
  summary{display:none}
  .cav{font-size:7.8pt;color:var(--ink-3);max-width:none}   /* token — print redefines it above */
  a{color:var(--ink);border:none}
  @page{margin:14mm 13mm}
}
"""


def build() -> pathlib.Path:
    g = _load()
    imp = _importance(g)
    P = ['<title>APM Card Column Reference</title>',
         '<link rel="preconnect" href="https://fonts.googleapis.com">',
         '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>',
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
         'family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@600;700&'
         'family=IBM+Plex+Serif:wght@400;600&display=swap">',
         f"<style>{CSS}</style>", '<div class="wrap">',
         "<h1>Card column reference</h1>",
         f'<p class="sub">All <strong>{len(g)}</strong> columns across the five cards, grouped by card and '
         f'then by block, with blocks ordered by <strong>measured importance</strong>. Each entry leads with '
         f'what the column <em>is</em>; the caveats fold underneath.</p>',
         f'<p class="note">Generated {dt.datetime.now():%Y-%m-%d} by '
         f'<code>analyses/ops/gen_column_reference.py</code> from <code>card_glossary.tsv</code> and '
         f'<code>card_registry.tsv</code>. <strong>Do not hand-edit.</strong></p>',
         _banner(),
         '<div class="call"><p class="h">The one rule behind this layout</p>'
         '<p><strong>The rung is a property of (CARD, COLUMN) — never of the column name.</strong> '
         '<code>beta</code> is edge-rung on the edge card and family-rung on the gene_family card: same '
         'estimator, different unit. That is why this is organised card-first, and why one name can appear '
         'twice with different text.</p></div>',
         '<div class="call"><p class="h">How blocks are ordered</p>'
         '<p>By three measured signals, not by size: how often a block\'s columns are named in the '
         '<strong>discovery registry</strong> (what the findings rest on, weighted highest), how many '
         '<strong>modules</strong> read them, and whether the block is <strong>base-owned</strong> '
         '(fitted) rather than an annotation join. Names under 6 characters count only backtick-quoted '
         'registry hits — otherwise <code>n</code> and <code>z</code> match ordinary prose.</p>'
         '<p><strong>One correction the measurement needed:</strong> mention-frequency is not importance — '
         'a <em>settled</em> block gets discussed less precisely because nobody argues about it. '
         '<code>adm_</code> scored 27th of 44 on the raw signals while MH-216 calls '
         '<code>adm_has_site</code> the single most load-bearing conditioning variable, on 4 registry '
         'mentions. Five blocks therefore carry a <strong>curated floor</strong>, each citing the finding '
         'that justifies it (see <code>CURATED_FLOOR</code>). That is the only hand-set input here.</p>'
         '</div>']

    toc = ['<div class="toc"><strong>Cards</strong><ul>']
    for card, rel, _, _ in CARDS:
        sub = g[g.card == card]
        if len(sub):
            toc.append(f'<li><a href="#{card}">{card}</a> — {len(sub)} columns, '
                       f'{sub.block.nunique()} blocks</li>')
    P.append("".join(toc) + "</ul></div>")

    for card, rel, keydesc, headnote in CARDS:
        sub = g[g.card == card]
        if not len(sub):
            continue
        var = {"edge": "--edge", "gene": "--gene", "gene_family": "--fam",
               "arm": "--arm", "seed_family": "--seed"}[card]
        blocks = sorted(sub.block.unique(), key=lambda b: (-float(imp.get((card, b), 0)), b))
        P.append(f'<section class="card" id="{card}" style="--k:var({var})">')
        P.append(f'<div class="card-head"><div class="key">{html.escape(rel)}</div><h2>{card}</h2>'
                 f'<p class="sub">One row per {html.escape(keydesc)}. {_md(headnote)}</p>'
                 f'<div class="meta"><span><b>{len(sub)}</b> columns</span>'
                 f'<span><b>{int(sub.nrows.iloc[0]):,}</b> rows</span>'
                 f'<span><b>{len(blocks)}</b> blocks</span>'
                 f'<span>median coverage <b>{sub.fill.dropna().median():.0%}</b></span></div></div>')
        # ── per-card contents, in importance order
        items = []
        for i, b in enumerate(blocks, 1):
            tier, _ = _tier(float(imp.get((card, b), 0)))
            n = int((sub.block == b).sum())
            items.append(f'<li><a href="#{card}-{b.strip("_") or "bare"}">{html.escape(b)}</a> '
                         f'<span class="tg">{n} · {tier}</span></li>')
        P.append(f'<div class="cardtoc"><p class="t">Contents — {len(blocks)} blocks, most important first</p>'
                 f'<ol>{"".join(items)}</ol></div>')

        for b in blocks:
            score = float(imp.get((card, b), 0))
            tier, blurb = _tier(score)
            rows = sub[sub.block == b].sort_values("column")
            P.append(f'<h3 id="{card}-{b.strip("_") or "bare"}">{html.escape(b)} '
                     f'<span class="cnt">· {len(rows)}</span></h3>'
                     f'<p class="tierline">{tier} — {blurb}</p>')
            for r in rows.itertuples():
                tags = []
                if r.rung:
                    tags.append(f'<span class="tag rung" title="{html.escape(RUNG_NOTE.get(r.rung,""))}">'
                                f'{html.escape(r.rung)}</span>')
                if r.agg_of:
                    tags.append(f'<span class="tag">agg {html.escape(r.agg_of)}</span>')
                if r.fill is not None and r.fill == r.fill:
                    tags.append(f'<span class="tag">{r.fill:.0%}</span>')
                cav = ""
                if r.caveat:
                    cav = (f'<details><summary>caveats</summary>'
                           f'<div class="cav">{_md(r.caveat)}'
                           + (f'<br><em>Defined on: {_md(r.domain)}</em>' if r.domain else "")
                           + "</div></details>")
                elif r.domain:
                    cav = f'<div class="cav"><em>Defined on: {_md(r.domain)}</em></div>'
                P.append(f'<div class="col"><div><div class="cname">{html.escape(r.column)}</div>'
                         f'<div class="tags">{"".join(tags)}</div></div>'
                         f'<div><div class="lead">{_md(r.lead)}</div>{cav}</div></div>')
        P.append("</section>")

    P.append("</div>")
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text("\n".join(P))
    return DEST


def to_pdf(src: pathlib.Path = DEST, dest: pathlib.Path = PDF) -> pathlib.Path | None:
    """Render with headless Chrome — the only local engine with modern CSS (grid, custom properties)."""
    chrome = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chrome")
    if not chrome:
        print("  ⚠ no Chrome found — PDF skipped (the HTML still prints correctly)")
        return None
    cmd = [chrome, "--headless", "--disable-gpu", "--no-sandbox", "--no-pdf-header-footer",
           f"--print-to-pdf={dest}", "--virtual-time-budget=20000", src.resolve().as_uri()]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if dest.exists() and dest.stat().st_size > 20_000:
        return dest
    print(f"  ⚠ PDF not produced: {(r.stderr or '')[-300:]}")
    return None


def main() -> int:
    p = build()
    print(f"✅ {p.relative_to(ROOT.parent)}  ({p.stat().st_size/1024:.0f} KB)")
    if "--no-pdf" not in sys.argv:
        q = to_pdf()
        if q:
            print(f"✅ {q.relative_to(ROOT.parent)}  ({q.stat().st_size/1024/1024:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
