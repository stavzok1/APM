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

import numpy as np
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


def _anchor(card: str, block: str) -> str:
    """A safe fragment id. `(bare)` used to produce `edge-(bare)` — parentheses in an id are legal HTML5
    but break naive selectors and were silently unlinkable from the contents list."""
    return f"{card}-" + (re.sub(r"[^A-Za-z0-9]+", "", block) or "bare")


def _block(col: str) -> str:
    """The column's block = its first underscore-delimited token.

    ⛔ **CASE-INSENSITIVE, and it was not (user-caught).** `^([a-z]+_)` fails on a camelCase prefix — the
    `G` in `dGlobal_HLY_NAT` breaks the match — so **8 columns (`dGlobal_*`, `dShare_*`) fell into
    `(bare)`**. Worse, three of them share one description, which then became the MODAL text for `(bare)`
    and was printed as that block's SUMMARY: a specific sentence about cohort-wide dose rank presented as
    the description of every unprefixed column. ⇒ a naming convention the code half-recognises is worse
    than one it ignores entirely.
    """
    m = re.match(r"^([A-Za-z]+_)", col)
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


BOOLS = {"True", "False", "true", "false", "1.0", "0.0", "1", "0"}

#: ⭐⭐ TOKEN DECODER — so **every row carries a line of its own**, not a pointer to the block.
#:
#: ⛔ The problem it solves: de-duplicating the block text (MH-272) left **425 columns with nothing but
#: *"described by the block above"*** — which is worse for a lookup than the repetition was. But the
#: glossary genuinely has no per-column text for them, so there is nothing to render. What IS knowable is
#: what distinguishes a column from its block siblings: a systematic token. `share_NAT` differs from
#: `share_TUM` by exactly one, and that token is the column's summary.
#:
#: ✅ Measured: **346 of the 425 (81%)** carry at least one known token. The remaining **79 are a REAL
#: GLOSSARY GAP** and are listed by `--gaps`, not papered over — among them `ctx_ceiling`, one of the most
#: cited columns in the program. Those fall back to the block lead, labelled as block-level.
#:
#: ⚠ Every rendered line says where it came from — **authored / derived / block-level** — because a
#: mechanically composed line must never be mistaken for one somebody wrote and checked.
TOKENS = {
    "hly": "in the GTEx-healthy state", "nat": "in matched normal (NAT)", "tum": "in tumour",
    "prosp": "CPTAC prospective cohort (n≈101)", "t105": "CPTAC TCGA-overlap cohort (n=105)",
    "cptac": "CPTAC", "tcga": "TCGA", "buffa": "the Buffa cohort", "gtex": "GTEx",
    "prot": "against protein", "rna": "against mRNA", "disc": "against the discovery layer",
    "rho": "Spearman ρ", "r": "correlation", "p": "p-value", "q": "BH-adjusted q", "z": "z-score",
    "sd": "standard deviation", "med": "median", "median": "median", "mean": "mean",
    "max": "maximum", "min": "minimum", "frac": "fraction", "pct": "percent", "n": "count",
    "sum": "sum", "total": "total", "iqr": "interquartile range", "rank": "rank",
    "raw": "UNADJUSTED — no confounder block", "adj": "composition-adjusted",
    "deconv": "composition-adjusted (deconvolution block)", "core": "core confounder block",
    "agg": "β-weighted aggregate", "abund": "unweighted abundance-sum reference",
    "gated": "gated", "perm": "permutation null", "oof": "out-of-fold", "meas": "measured only (un-imputed)",
    "fam": "per seed family", "arm": "per arm", "gene": "per gene", "cell": "per family cell",
    "edge": "per edge", "genes": "genes", "arms": "arms", "members": "family members",
    "top": "the top-ranked", "best": "the best", "dominant": "the dominant one",
    "8mer": "8mer sites", "7mer": "7mer sites", "6mer": "6mer sites",
    "canonical": "canonical sites", "noncanon": "non-canonical sites",
    "mti": "miRNA–target interaction evidence", "functional": "functional assays",
    "weak": "weak evidence", "studies": "study count", "assay": "by assay class",
    "pmid": "distinct publications", "chimeric": "chimeric duplex evidence",
    "share": "share", "shift": "shift", "gap": "gap", "delta": "change", "d": "change",
    "dose": "dose", "budget": "pressure budget", "retention": "retention",
    "class": "class label", "src": "source", "beta": "β", "identity": "Shapley identity",
}
_THRESH = re.compile(r"^(c|b)(\d{2,3})$")


def _decode(column: str, block: str) -> str:
    """Compose a line from the tokens that distinguish this column from its block siblings."""
    rest = column[len(block):] if block != "(bare)" else column
    parts = []
    for t in re.split(r"[_.]", rest):
        if not t:
            continue
        tl = t.lower()
        m = _THRESH.match(tl)
        if m:
            parts.append(f"at the {'−0.' if m.group(1) == 'c' else '0.'}{m.group(2)} threshold")
        elif tl in TOKENS:
            parts.append(TOKENS[tl])
    return " · ".join(dict.fromkeys(parts))


def _vtype(s: pd.Series, name: str) -> tuple[str, str]:
    """The VALUE TYPE a reader needs before using a column, inferred from the data. (label, range hint).

    ⚠ **ORDER MATTERS — three rules were mis-ordered on the first pass and validation caught all three:**
    `constant` must precede `boolean` (`pip_dense` is 1.0 everywhere and read as boolean); the p-value and
    percent NAME-hints must precede the integer check (`arm_pct_above_floor` is 0–100 integers and read as
    a count); and `identifier` keys on the distinct COUNT, not a share of rows (`gene` has 1,420 distinct
    in 5,649 rows = 25%, under a 50% threshold). ✅ Validated 10/10 against columns whose type is known.
    """
    v = s.dropna()
    if v.empty:
        return "empty", "no values"
    if v.nunique() == 1:
        return "constant", f"= {v.iloc[0]}"
    if set(map(str, v.unique())) <= BOOLS and v.nunique() <= 2:
        t = int(v.astype(str).isin({"True", "true", "1", "1.0"}).sum())
        return "boolean", f"{t:,} true / {len(v)-t:,} false"
    n = pd.to_numeric(v, errors="coerce")
    if n.notna().mean() < 0.9:
        k = v.nunique()
        return ("identifier", f"{k:,} distinct") if k > 50 else \
               ("categorical", f"{k} level" + ("s" if k != 1 else ""))
    n = n.dropna()
    lo, hi = float(n.min()), float(n.max())
    if re.search(r"(^|_)(p|q)(_|$)|_p$|_q$|_pval", name) and 0 <= lo and hi <= 1.0001:
        return "p-value", f"{lo:.3g} – {hi:.3g}"
    if re.search(r"pct|percent", name) and 0 <= lo and hi <= 100.001:
        return "percent", f"{lo:.3g} – {hi:.3g}"
    if bool(np.allclose(n, n.round())) and lo >= 0:
        return "count", f"{lo:g} – {hi:g}"
    if -1.0001 <= lo and hi <= 1.0001:
        return ("fraction 0–1", f"{lo:.3g} – {hi:.3g}") if lo >= -1e-9 else \
               ("signed −1…1", f"{lo:+.3g} – {hi:+.3g}")
    return ("real ≥ 0" if lo >= 0 else "real, signed"), f"{lo:.3g} – {hi:.3g}"


def _load() -> pd.DataFrame:
    g = pd.read_csv(OUT / "card_glossary.tsv", sep="\t", dtype=str).fillna("")
    fills, nrows, frames = {}, {}, {}
    for card, rel, _, _ in CARDS:
        p = OUT / rel
        if p.exists():
            d = pd.read_csv(p, sep="\t", low_memory=False)
            fills[card], nrows[card], frames[card] = d.notna().mean(), len(d), d
    g["block"] = g.column.map(_block)
    g["fill"] = [float(fills[c].get(col)) if c in fills and col in fills[c] else None
                 for c, col in zip(g.card, g.column)]
    g["nrows"] = [nrows.get(c, 0) for c in g.card]
    tv = [_vtype(frames[c][col], col) if c in frames and col in frames[c] else ("—", "")
          for c, col in zip(g.card, g.column)]
    g["vtype"] = [a for a, _ in tv]
    g["vhint"] = [b for _, b in tv]
    lead_rest = [_split(t) for t in g.description]
    g["lead"] = [a for a, _ in lead_rest]
    g["caveat"] = [b for _, b in lead_rest]

    # ⛔⛔ DE-DUPLICATION. `describe()` falls back to a PREFIX BLOCK when a column has no exact entry, so
    # every column in a block inherits the same text. Measured: **437 of 711 rows (61%)** repeated a
    # sibling's description — **159,456 characters** printed over and over — and **all 74 repeat groups were
    # exactly one block**. ⇒ the shared text is a property of the BLOCK, so it is hoisted into a block
    # summary and shown ONCE; a column then carries only what is TRUE OF IT ALONE. A column with no own
    # text is not blank — it still shows its rung, type, range and coverage, and the block summary above
    # explains the family.
    # ⛔ FIXED (user-caught): this used to flag only the block's MODAL description, so a text shared by a
    # NON-modal pair slipped through and still printed twice — `fam_ctx_composition` and its sibling being
    # the reported case. ANY description occurring more than once within a (card, block) is shared.
    g["shared"] = False
    for (card, blk), sub in g.groupby(["card", "block"]):
        counts = sub.description.value_counts()
        rep = set(counts[counts > 1].index)
        if rep:
            g.loc[sub.index[sub.description.isin(rep)], "shared"] = True
    g["_modal"] = False
    for (card, blk), sub in g.groupby(["card", "block"]):
        counts = sub.description.value_counts()
        if len(counts) and counts.iloc[0] > 1:
            g.loc[sub.index[sub.description == counts.index[0]], "_modal"] = True
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
.blocksum{background:var(--panel);border:1px solid var(--rule);border-left:3px solid var(--k);
          border-radius:2px;padding:11px 14px;margin:0 0 10px}
.blead{margin:0 0 6px;color:var(--ink-2);font-size:14px;max-width:82ch}
.blead strong{color:var(--ink)}
.bstats{margin:2px 0 0;font-family:var(--mono);font-size:11px;color:var(--ink-3)}
.bstats b{color:var(--ink-2)}
.tag.ty{border-color:var(--accent);color:var(--accent)}
.lead.sh{color:var(--ink-3)}
.prov{font-family:var(--mono);font-size:9px;letter-spacing:.06em;text-transform:uppercase;
      margin-left:7px;padding:1px 4px;border-radius:2px;vertical-align:1px;white-space:nowrap}
.prov.auth{color:var(--gene);background:color-mix(in srgb,var(--gene) 12%,transparent)}
.prov.der{color:var(--ink-3);background:var(--rule-2)}
.prov.gap{color:var(--seed);background:color-mix(in srgb,var(--seed) 12%,transparent)}
.rng{font-family:var(--mono);font-style:normal;font-size:12px}
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
  .blocksum{padding:7px 10px;margin-bottom:7px;break-inside:avoid}
  .blead{font-size:8.6pt;max-width:none}.bstats{font-size:7.4pt}
  .prov{font-size:6.4pt;border:1px solid var(--rule);background:none}
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
            items.append(f'<li><a href="#{_anchor(card, b)}">{html.escape(b)}</a> '
                         f'<span class="tg">{n} · {tier}</span></li>')
        P.append(f'<div class="cardtoc"><p class="t">Contents — {len(blocks)} blocks, most important first</p>'
                 f'<ol>{"".join(items)}</ol></div>')

        for b in blocks:
            score = float(imp.get((card, b), 0))
            tier, blurb = _tier(score)
            rows = sub[sub.block == b].sort_values("column")
            P.append(f'<h3 id="{_anchor(card, b)}">{html.escape(b)} '
                     f'<span class="cnt">· {len(rows)}</span></h3>'
                     f'<p class="tierline">{tier} — {blurb}</p>')

            # ── BLOCK SUMMARY: the text every column in this block shared, hoisted and shown ONCE,
            #    plus what the block looks like as a whole.
            shared_rows = rows[rows._modal]
            types = rows.vtype.value_counts()
            typemix = " · ".join(f"{n}&times; {html.escape(t)}" for t, n in types.items())
            cov = rows.fill.dropna()
            bits = [f"<b>{len(rows)}</b> columns"]
            if len(cov):
                bits.append(f"coverage <b>{cov.min():.0%}–{cov.max():.0%}</b> (median {cov.median():.0%})")
            if len(shared_rows):
                bits.append(f"<b>{len(shared_rows)}</b> share the block description")
            summ = ""
            if len(shared_rows):
                lead, cav = _split(shared_rows.description.iloc[0])
                summ = (f'<p class="blead">{_md(lead)}</p>'
                        + (f'<details><summary>block caveats</summary>'
                           f'<div class="cav">{_md(cav)}</div></details>' if cav else ""))
            P.append(f'<div class="blocksum">{summ}'
                     f'<p class="bstats">{" &nbsp;·&nbsp; ".join(bits)}</p>'
                     f'<p class="bstats">{typemix}</p></div>')

            for r in rows.itertuples():
                tags = []
                if r.rung:
                    tags.append(f'<span class="tag rung" title="{html.escape(RUNG_NOTE.get(r.rung,""))}">'
                                f'{html.escape(r.rung)}</span>')
                tags.append(f'<span class="tag ty" title="{html.escape(r.vhint)}">'
                            f'{html.escape(r.vtype)}</span>')
                if r.agg_of:
                    tags.append(f'<span class="tag">agg {html.escape(r.agg_of)}</span>')
                if r.fill is not None and r.fill == r.fill:
                    tags.append(f'<span class="tag">{r.fill:.0%}</span>')
                # ⭐ EVERY row gets a line, and every line says where it came from.
                if not r.shared:
                    body = (f'<div class="lead">{_md(r.lead)}'
                            f'<span class="prov auth" title="written and checked by hand">authored</span>'
                            f'</div>')
                else:
                    dec = _decode(r.column, r.block)
                    if dec:
                        body = (f'<div class="lead">{_md(dec)}'
                                f'<span class="prov der" title="composed from the tokens that distinguish '
                                f'this column from its block siblings — mechanical, not authored">'
                                f'derived</span></div>')
                    else:
                        # ⚠ a REAL glossary gap — say so rather than implying a column-specific entry
                        body = (f'<div class="lead sh">{_md(_split(r.description)[0])}'
                                f'<span class="prov gap" title="no column-specific entry exists yet — '
                                f'this is the BLOCK description">block-level</span></div>')
                cav = ""
                if not r.shared and r.caveat:
                    cav = (f'<details><summary>caveats</summary><div class="cav">{_md(r.caveat)}'
                           + (f'<br><em>Defined on: {_md(r.domain)}</em>' if r.domain else "")
                           + "</div></details>")
                elif r.domain:
                    cav = f'<div class="cav"><em>Defined on: {_md(r.domain)}</em></div>'
                P.append(f'<div class="col"><div><div class="cname">{html.escape(r.column)}</div>'
                         f'<div class="tags">{"".join(tags)}</div></div>'
                         f'<div>{body}{cav}</div></div>')
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
    if "--gaps" in sys.argv:
        g = _load()
        rows = [r for r in g.itertuples() if r.shared and not _decode(r.column, r.block)]
        print(f"⚠ {len(rows)} column(s) have NO column-specific glossary entry AND no decodable token.")
        print("  These fall back to the BLOCK description in the reference. They are the real writing queue:\n")
        for r in sorted(rows, key=lambda x: (x.card, x.column)):
            print(f"  {r.card:<13}{r.column}")
        return 0
    p = build()
    print(f"✅ {p.relative_to(ROOT.parent)}  ({p.stat().st_size/1024:.0f} KB)")
    if "--no-pdf" not in sys.argv:
        q = to_pdf()
        if q:
            print(f"✅ {q.relative_to(ROOT.parent)}  ({q.stat().st_size/1024/1024:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
