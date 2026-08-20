"""Card viewer — a searchable, self-contained HTML browser for the five-rung card system.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_cards --build      # generate the viewer
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_cards --check      # registry<->card gate
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_cards --roundtrip  # encode+decode, lossless test
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_cards --verify-js  # run the pages' OWN js vs the TSV
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_cards --links      # cross-card referential audit
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_cards --stamp      # provenance footer only

Output: `docs/derived/cards/{index,arm,gene,edge,family,seed_family}.html` — six standalone pages, no
server, no network, no dependencies. Search an arm / gene / edge / family, read its card.

The audits (`--check`, `--links`, `--roundtrip`, `--verify-js`) are NOT scaffolding for the build; they
shipped first and stand alone. Two of them are checks this subproject otherwise lacks — card integrity
and cross-card referential integrity — and each of the four has caught a real defect (see below).

## What the reader must never be able to misread

  * ⭐ **THE RUNG TRAVELS WITH EVERY VALUE.** `card_registry.tsv` gives `(card, column) -> rung`. The
    rung is a property of the (CARD, COLUMN) PAIR, never of the column name — `beta` is edge-rung on one
    card and family-rung on another. A value that is not at the card's own grain gets a ⟨rung⟩ chip,
    because it repeats across rows BY CONSTRUCTION and reading it as a per-row fact is MH-179/187/188/191.
  * ⭐ **A BLANK IS NOT A ZERO, AND THERE ARE THREE OF THEM:** a measured `0`; `— not scanned`, where a
    `cov_*` flag says the scan never covered this row; and `— no value`, scanned but empty. ⛔ The
    flag->block mapping is MEASURED, never curated (`cov_map`) — a guessed mapping would make the page
    assert "not scanned" about data that was scanned.
  * ⭐ **THE PERCENTILE'S DENOMINATOR IS PRINTED** (`p68 of 320`, not `of 1,959`) and is computed over
    non-NaN rows only. A percentile without its denominator is the axiom-5 failure in a tooltip.
  * ⛔ **The four targetome universes get four separate axes.** One grouped chart is the VISUAL ACT of
    blending them; the registry forbids it in words, `charts()` forbids it in pixels.
  * ⛔ **Dangling cross-card links render as muted TEXT, never as a link to a blank page** — and a key
    present under a different label form says so, rather than claiming absence.

## Defects these checks caught (all mine, all before a reader saw them)

  * `--roundtrip` — group compression applied on the card's OWN grain corrupted 21 of 61 family columns,
    worst relative error 2.5e+05.
  * `--links` — 42 of 50 "missing" keys were present under a different label form; reporting them as
    absent would have been a false negative.
  * `--verify-js` — booleans reach the browser as the STRINGS "True"/"False", so the coverage chart's
    `=== true` matched nothing: every flag rendered OFF and the "not scanned" state was dead. The page
    LOOKED correct. Also: `— not scanned` never rendered at all, because the hide-empty-blocks rule was
    hiding exactly the blocks whose emptiness was the message (1,860 of 2,450 arms).
  * step 5 — `distribution()` was called after its loop on a stale variable, so `dist` shipped nearly
    empty for weeks; the round-trip could not see it because it never consumed the descriptors.

⚠ Note the pattern: every one of these passed the checks that existed at the time. A test that does not
CONSUME an output cannot defend it.

## Why the gate exists — it is not hypothetical

`card_registry.tsv` is BUILT FROM the cards (`card_rungs.build()` reads each card's actual columns), so
registry and card can only disagree if **a card changed after the registry was written**. That is exactly
what a clobber looks like, and it happened on 2026-08-03:

    gene_family_card.tsv (22:07) became byte-identical to readouts_edges.tsv — 27 columns where the registry
    (21:48) declared 61. It had silently lost its whole ctx_ (16) and comp_ (16) annotation blocks plus
    n_fam/w_max. Nothing noticed for ~14 hours.

Two independent signals catch it, and this module reports BOTH because either alone can miss:
  * **column conformance** — registry-declared columns absent from the card (the clobber's fingerprint),
    and card columns the registry has never seen (the registry is stale and must be regenerated).
  * ⭐ **mtime ordering** — any card NEWER than the registry means the rung labels may not describe the
    data on disk. This fires even when the column sets happen to match, which is the case a pure column
    diff misses: a card can be rewritten with the same header and different semantics.

⚠ A passing gate is NOT proof the cards are correct — only that they match what the registry recorded.
Rung correctness is `card_rungs --check`'s job; this is the layer beneath it.

Conventions follow `analyses/ops/gen_architecture.py` (module docstring with a Run: line, `run()`,
`print("[gen-cards] …")`), and `CARDS` is imported from `learned/card_rungs.py` rather than re-declared —
paths have one home.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.learned.card_rungs import CARDS

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
REGISTRY = OUT / "card_registry.tsv"


def _stamp(p: Path) -> dict:
    """mtime + size + content hash. The hash is what makes the stamp a real identity, not a guess:
    a file can be rewritten with an identical size and a touched mtime."""
    b = p.read_bytes()
    return {"path": p, "exists": True, "mtime": p.stat().st_mtime,
            "size": len(b), "sha": hashlib.sha256(b).hexdigest()[:8]}


def sources() -> dict:
    """Every input the viewer depends on, stamped. Ordered registry-first — it is the schema."""
    out = {"registry": _stamp(REGISTRY) if REGISTRY.exists()
           else {"path": REGISTRY, "exists": False}}
    for name, spec in CARDS.items():
        p = spec["path"]
        out[name] = _stamp(p) if p.exists() else {"path": p, "exists": False}
    return out


def check(verbose: bool = True) -> list[str]:
    """Return a list of problems; empty means the cards match the registry. Never raises on data."""
    problems: list[str] = []
    src = sources()

    if not src["registry"]["exists"]:
        return [f"registry missing: {REGISTRY} — run card_rungs --check"]
    reg = pd.read_csv(REGISTRY, sep="\t")
    reg_mtime = src["registry"]["mtime"]

    if verbose:
        print(f"[gen-cards] registry {REGISTRY.name}  {len(reg):,} rows  "
              f"sha {src['registry']['sha']}  {_ts(reg_mtime)}")
        print(f"\n{'card':<13s}{'registry':>9s}{'on disk':>9s}{'missing':>9s}{'extra':>7s}  source")

    for name in CARDS:
        s = src[name]
        if not s["exists"]:
            problems.append(f"{name}: card file missing ({s['path']})")
            if verbose:
                print(f"{name:<13s}{'—':>9s}{'MISSING':>9s}")
            continue

        declared = set(reg.loc[reg.card == name, "column"])
        have = list(pd.read_csv(s["path"], sep="\t", nrows=1, low_memory=False).columns)
        missing, extra = sorted(declared - set(have)), sorted(set(have) - declared)

        if verbose:
            print(f"{name:<13s}{len(declared):>9d}{len(have):>9d}{len(missing):>9d}{len(extra):>7d}"
                  f"  {s['path'].name} sha {s['sha']} {_ts(s['mtime'])}")

        # ⛔ the clobber fingerprint: the registry saw columns the card no longer has
        if missing:
            problems.append(f"{name}: {len(missing)} declared column(s) ABSENT from the card — "
                            f"the card may have been overwritten by a less complete build. "
                            f"e.g. {', '.join(missing[:6])}")
        # the registry is behind the card — labels exist for fewer columns than shipped
        if extra:
            problems.append(f"{name}: {len(extra)} column(s) on the card the registry has never seen — "
                            f"regenerate it (card_rungs --check). e.g. {', '.join(extra[:6])}")
        # ⭐ catches a same-header rewrite, which a column diff cannot see
        if s["mtime"] > reg_mtime:
            problems.append(f"{name}: card is NEWER than the registry ({_ts(s['mtime'])} > "
                            f"{_ts(reg_mtime)}) — the rung labels may not describe this data.")

    if verbose:
        print()
        if problems:
            print(f"[gen-cards] {len(problems)} PROBLEM(S):")
            for p in problems:
                print(f"   ⛔ {p}")
        else:
            print("[gen-cards] CLEAN — every card matches the registry, and the registry is newest.")
        print("[gen-cards] ⚠ a passing gate means the cards match what the registry RECORDED; "
              "rung correctness is card_rungs --check.")
    return problems


def _ts(t: float) -> str:
    import datetime as _dt
    return _dt.datetime.fromtimestamp(t).strftime("%m-%d %H:%M")


# --------------------------------------------------------------------------- #
# STEP 2 — payload encoding
#
# The viewer must embed ~19k rows × 699 columns in self-contained HTML. Per-row
# dicts are hopeless: every NaN pays for its key AND the 4 bytes of `null`, and
# these cards are 25-70% NaN BY DESIGN (a blank means UNSCANNED, not zero).
# Column-major lets each column choose its own encoding, and the detail view
# (one entity, all columns) is a strided read either way.
#
# ⛔⛔ THE ONE INVARIANT: the NaN MASK IS SACRED. A blank must decode to a blank.
# NaN here means "this scan never covered this row" — inventing a value would be
# the viewer fabricating data, the same error class the cards exist to prevent.
# `roundtrip()` asserts the mask EXACTLY (no tolerance); only the numeric values
# are compared with tolerance, because rounding is deliberate and declared.
# --------------------------------------------------------------------------- #
SEP = "\x1f"            # unit separator — cannot occur in a gene symbol or arm name
DICT_MAX_FRAC = 0.60    # dictionary-encode a string column below this distinct-fraction


def _fmt(x, sigfig: int) -> str:
    """Round to significant figures and drop the noise. `0.6190476190476191` -> `0.619`.

    ⚠ LOSSY, and deliberately so — this is a viewer, not the source of truth. The TSV keeps full
    precision, the footer names the sigfig, and `roundtrip()` compares numerics with a matching
    tolerance rather than pretending the loss is not there.
    """
    if x is None or x != x:
        return ""
    f = float(x)
    if f == 0:
        return "0"
    r = round(f, -int(__import__("math").floor(__import__("math").log10(abs(f)))) + (sigfig - 1))
    return repr(int(r)) if r == int(r) and abs(r) < 1e15 else repr(r)


def _group_key(card: str, rung: str) -> str | None:
    """Which card column a `rung` is constant within — the basis for group compression.

    ⛔⛔ A RUNG THAT NAMES THE CARD'S OWN GRAIN IS NOT COMPRESSIBLE, and assuming otherwise CORRUPTED
    the family card until the round-trip caught it. `card_rungs`' docstring states the asymmetry:

        on a (gene, arm)-keyed card:    arm  => constant within arm ACROSS genes   [compressible]
        on a (gene, family)-keyed card: family => FREE (it IS the grain)           [NOT compressible]

    On the family card, a rung='family' column like `beta` is a property of THIS (gene, family) cell,
    not of the family — grouping by `family` merges rows from different genes and broadcasts one
    gene's β over all of them. Measured before the guard: 21 of 61 family columns corrupted, worst
    relative error 2.5e+05.

    ⇒ compressible only if the rung names a key column that is NOT the last (grain-defining) key.
    And even then `encode_card` VERIFIES constancy per column before using it — the declaration is
    tested, never trusted (the subproject's own rule for rung labels).
    """
    keys = CARDS[card]["key"]
    if len(keys) < 2:
        return None                       # single-key card: every row is its own group, no win
    if rung not in keys or rung == keys[-1]:
        return None                       # not a key, or it IS the grain
    return rung


def encode_column(s: pd.Series, *, sigfig: int, group: pd.Series | None = None) -> dict:
    """Encode one column. Four shapes, chosen by dtype/cardinality, all NaN-exact."""
    if group is not None:
        return _encode_grouped(s, group, sigfig)
    if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
        return {"t": "n", "v": ",".join(_fmt(x, sigfig) for x in s)}
    v = s.astype("object").where(s.notna(), None)
    vals = [None if x is None else str(x) for x in v]
    distinct = {x for x in vals if x is not None}
    if len(distinct) <= max(2, DICT_MAX_FRAC * max(len(vals), 1)):
        levels = sorted(distinct)
        ix = {lv: i for i, lv in enumerate(levels)}
        return {"t": "d", "l": levels, "i": ",".join("" if x is None else str(ix[x]) for x in vals)}
    return {"t": "s", "v": SEP.join("" if x is None else x for x in vals)}


def _encode_grouped(s: pd.Series, group: pd.Series, sigfig: int) -> dict:
    """⭐ RUNG-AWARE GROUP COMPRESSION — store one value per group, not per row.

    ⛔ AND PRESERVE THE ROW-LEVEL NaN MASK. A row that is blank stays blank even when its group
    carries a value: `x` lists exactly those rows (delta-encoded). Broadcasting the group value over
    them would fill in data nobody measured.
    """
    order = list(dict.fromkeys(group.astype(str)))       # first-appearance order, stable
    gi = {g: i for i, g in enumerate(order)}
    gvals: list = [None] * len(order)
    for g, sub in s.groupby(group.astype(str), sort=False):
        nn = sub.dropna()
        gvals[gi[g]] = nn.iloc[0] if len(nn) else None
    is_num = pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s)
    # rows blank despite their group having a value — the mask that must survive
    gser = group.astype(str).map(lambda g: gvals[gi[g]])
    exc = [i for i, (a, b) in enumerate(zip(s.isna(), gser.notna())) if a and b]
    out = {"g": ",".join(order), "x": ",".join(map(str, exc))}
    if is_num:
        out |= {"t": "gn", "v": ",".join(_fmt(x, sigfig) for x in gvals)}
    else:
        out |= {"t": "gs", "v": SEP.join("" if x is None else str(x) for x in gvals)}
    return out


def decode_column(enc: dict, n: int, group: pd.Series | None = None) -> list:
    """Inverse of `encode_column` — the reference decoder the JS mirrors, and what `roundtrip` uses."""
    t = enc["t"]
    if t == "n":
        return [None if x == "" else float(x) for x in enc["v"].split(",")]
    if t == "d":
        lv = enc["l"]
        return [None if x == "" else lv[int(x)] for x in enc["i"].split(",")]
    if t == "s":
        return [None if x == "" else x for x in enc["v"].split(SEP)]
    order = enc["g"].split(",")
    raw = (enc["v"].split(",") if t == "gn" else enc["v"].split(SEP))
    gv = [None if x == "" else (float(x) if t == "gn" else x) for x in raw]
    gi = {g: i for i, g in enumerate(order)}
    out = [gv[gi[g]] for g in group.astype(str)]
    for i in (int(i) for i in enc["x"].split(",") if i != ""):
        out[i] = None                                    # restore the row-level blank
    return out


def distribution(s: pd.Series) -> dict | None:
    """21 quantiles + a 16-bin histogram, over NON-NaN rows only.

    ⭐ `n` is shipped and is the point: a percentile is meaningless without its denominator.
    `fam_dominant_share` is defined on the 320 multi-member families, not all 1,959 — quoting
    `p68 of 320` vs `of 1,959` is the difference between a fact and a misreading.
    ⚠ Suppressed where a percentile would be noise: <8 observations or <3 distinct values.
    """
    import numpy as np
    if not pd.api.types.is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s):
        return None
    v = pd.to_numeric(s, errors="coerce").dropna().to_numpy(float)
    if len(v) < 8 or len(np.unique(v)) < 3:
        return None
    q = np.quantile(v, np.linspace(0, 1, 21))
    lo, hi = float(v.min()), float(v.max())
    counts, _ = np.histogram(v, bins=16, range=(lo, hi if hi > lo else lo + 1e-9))
    return {"n": int(len(v)), "q": [round(float(x), 6) for x in q],
            "h": [int(c) for c in counts], "lo": round(lo, 6), "hi": round(hi, 6)}


# --------------------------------------------------------------------------- #
# STEP 3 — cross-card links, and the referential audit that falls out of them
#
# `card_rungs --check` verifies rungs WITHIN a card. NOTHING checks that a key
# on one card resolves to a row on another — so a viewer that renders every
# cross-link as clickable would send readers to blank pages, and a blank page
# reads as "the data is broken" rather than "this card does not cover that arm".
#
# ⭐ Dangling is therefore a FIRST-CLASS STATE, not an error: the counts below
# are a true statement about card coverage, and reporting them is a referential
# audit this subproject does not currently have.
# --------------------------------------------------------------------------- #
#   (source card, source column, target card)  — target is matched on its single key
LINKS = [
    ("arm", "seed_family", "seed_family"),
    ("edge", "gene", "gene"),
    ("edge", "arm", "arm"),
    ("edge", "seed_family", "seed_family"),
    ("gene_family", "gene", "gene"),
    ("gene_family", "family", "seed_family"),
    ("seed_family", "fam_dominant_arm", "arm"),
]


def _has_label_variant(name: str, universe: set) -> bool:
    """Is `name` present in `universe` under a different LABEL FORM (arm-suffix / case)?

    ⚠ Detection only — never a join. A bare stem matching a suffixed family is precisely the 5p/3p
    conflation the arm key refuses; this reports the ambiguity so the viewer can say
    "present as miR-101-3p" instead of the false "not on the card".
    """
    import re
    n = str(name)
    if n.lower() in {u.lower() for u in universe}:
        return True
    pat = re.compile(rf"^{re.escape(n)}(-(3p|5p)(\.\d+)?)?(/|$)", re.I)
    return any(pat.match(u) or f"/{n}-" in u for u in universe)


def links(verbose: bool = True) -> pd.DataFrame:
    """Resolve every cross-card link; return one row per link with its dangling count.

    ⚠ The denominator is NON-NULL source values, not rows — a NaN link is "no link declared", which
    is a different fact from "declared and unresolvable", and conflating them would inflate coverage.
    """
    idx: dict[str, set] = {}
    for c, spec in CARDS.items():
        if len(spec["key"]) == 1 and spec["path"].exists():
            k = spec["key"][0]
            idx[c] = set(pd.read_csv(spec["path"], sep="\t", usecols=[k],
                                     low_memory=False)[k].dropna().astype(str))
    rows = []
    for src, col, tgt in LINKS:
        if not CARDS[src]["path"].exists() or tgt not in idx:
            continue
        head = pd.read_csv(CARDS[src]["path"], sep="\t", nrows=1, low_memory=False).columns
        if col not in head:
            continue
        s = pd.read_csv(CARDS[src]["path"], sep="\t", usecols=[col],
                        low_memory=False)[col].dropna().astype(str)
        ok = s.isin(idx[tgt])
        miss = sorted(set(s[~ok]))
        # ⭐ SPLIT "absent" FROM "present under a different LABEL FORM" — they mean opposite things.
        # Measured: 18 of 21 families dangling off the edge/family cards exist on the seed_family card
        # SUFFIXED (`miR-101` there is `miR-101-3p` here). Cause: `arm_card._arm_key` drops bare-stem
        # arms because they conflate 5p/3p, so the arm-derived seed_family card never sees their
        # families. Calling that "missing" would be a FALSE NEGATIVE in the viewer — the family is
        # there, under the resolved label.
        # ⛔ Reported, NOT auto-joined: matching a bare stem onto a suffixed family is exactly the
        # 5p/3p conflation `_arm_key` refuses. The reader is told; nothing is silently merged.
        variant = {m for m in miss if _has_label_variant(m, idx[tgt])}
        rows.append({"from": src, "column": col, "to": tgt, "declared": len(s),
                     "resolved": int(ok.sum()), "dangling": int((~ok).sum()),
                     "distinct_dangling": len(miss), "distinct_absent": len(miss) - len(variant),
                     "distinct_label_variant": len(variant),
                     "examples": ", ".join(sorted(set(miss) - variant)[:3] or miss[:3])})
    R = pd.DataFrame(rows)
    if verbose and len(R):
        print(f"{'link':<34s}{'declared':>9s}{'resolved':>9s}{'dangling':>9s}   abs/dist  examples")
        for _, r in R.iterrows():
            lab = f"{r['from']}.{r['column']} -> {r['to']}"
            pct = 100 * r["resolved"] / max(r["declared"], 1)
            print(f"{lab:<34s}{r['declared']:>9,d}{r['resolved']:>9,d}{r['dangling']:>9,d}"
                  f"  {pct:5.1f}%  {int(r['distinct_absent']):>3d}/{int(r['distinct_dangling']):>3d}"
                  f"  {r['examples'][:30]}")
        tot = int(R.dangling.sum())
        absent, variant = int(R.distinct_absent.sum()), int(R.distinct_label_variant.sum())
        print(f"\n[gen-cards] {tot:,} dangling value(s) over {len(R)} link type(s) — "
              f"{absent} distinct key(s) genuinely ABSENT, {variant} present under a DIFFERENT LABEL FORM.")
        print("[gen-cards] ⚠ dangling is a COVERAGE FACT, not a bug — muted text in the viewer, never a "
              "link to a blank page. ⛔ label variants are REPORTED, never auto-joined (that would be "
              "the 5p/3p conflation the arm key exists to refuse).")
    return R


def encode_card(card: str, *, sigfig: int = 4) -> dict:
    """Encode one whole card: columns + distributions + the group bases used."""
    reg = pd.read_csv(REGISTRY, sep="\t").query("card == @card")
    rung = dict(zip(reg.column, reg.rung))
    d = pd.read_csv(CARDS[card]["path"], sep="\t", low_memory=False)
    cols, dist, grouped = {}, {}, {}
    refused: list[str] = []
    for c in d.columns:
        gk = _group_key(card, str(rung.get(c, "")))
        g = d[gk] if gk and gk in d.columns and c not in CARDS[card]["key"] else None
        # ⭐ VERIFY THE DECLARATION BEFORE COMPRESSING ON IT. A declared rung that the data
        # contradicts is a mislabel (card_rungs finds those); here it would silently corrupt, so the
        # encoder simply declines and falls back to the plain column. Belt and braces with
        # `_group_key`'s grain rule — that rule is the reasoning, this is the proof.
        if g is not None and d.groupby(g, dropna=False)[c].nunique(dropna=True).max() > 1:
            refused.append(c)
            g = None
        cols[c] = encode_column(d[c], sigfig=sigfig, group=g)
        if g is not None:
            grouped[c] = gk
        dd = distribution(d[c])       # per column — see the fix note below
        if dd:
            dist[c] = dd
    if refused:
        print(f"  ⚠ {card}: {len(refused)} column(s) declared a rung they do not hold — "
              f"compression refused: {', '.join(refused[:5])}")
    return {"card": card, "key": CARDS[card]["key"], "n": len(d),
            "cols": list(d.columns), "d": cols, "dist": dist, "grouped": grouped}
    # ⚠ FIXED 2026-08-04 (step 5 surfaced it): the `distribution()` call used to sit AFTER the loop,
    # inside `if refused:`, reading the stale loop variable `c` — so `dist` shipped with at most ONE
    # column and usually zero. `roundtrip()` passed throughout because it tests the VALUES, not the
    # descriptors. A test that does not consume an output cannot defend it; the percentile strip is
    # the first consumer, and it found the bug on contact.


def roundtrip(card: str, *, sigfig: int = 4, verbose: bool = True) -> list[str]:
    """⭐ THE TEST THAT MAKES THE ENCODING TRUSTWORTHY: decode and compare against the source.

    Two standards, deliberately different:
      * **NaN mask — EXACT.** Any difference is a failure. This is the invariant.
      * **numeric values — tolerance 10^-(sigfig-1) relative.** Rounding is declared, not a bug.
      * **strings/categoricals — EXACT.**
    """
    import numpy as np
    d = pd.read_csv(CARDS[card]["path"], sep="\t", low_memory=False)
    pay = encode_card(card, sigfig=sigfig)
    bad: list[str] = []
    for c in d.columns:
        gk = pay["grouped"].get(c)
        got = decode_column(pay["d"][c], len(d), group=d[gk] if gk else None)
        src = d[c]
        gnull, snull = np.array([x is None for x in got]), src.isna().to_numpy()
        if not (gnull == snull).all():
            bad.append(f"{card}.{c}: NaN MASK differs on {int((gnull != snull).sum())} row(s)")
            continue
        if pd.api.types.is_numeric_dtype(src) and not pd.api.types.is_bool_dtype(src):
            a = np.array([np.nan if x is None else float(x) for x in got], float)
            b = pd.to_numeric(src, errors="coerce").to_numpy(float)
            m = ~np.isnan(b)
            if m.any() and not np.allclose(a[m], b[m], rtol=10.0 ** -(sigfig - 1), atol=0, equal_nan=True):
                worst = float(np.nanmax(np.abs((a[m] - b[m]) / np.where(b[m] == 0, 1, b[m]))))
                bad.append(f"{card}.{c}: numeric drift beyond {sigfig} sig-figs (worst rel {worst:.2e})")
        else:
            a = ["" if x is None else str(x) for x in got]
            b = ["" if x != x else str(x) for x in src.astype("object")]
            if a != b:
                n = sum(1 for x, y in zip(a, b) if x != y)
                bad.append(f"{card}.{c}: {n} value(s) differ after decode")
    if verbose:
        print(f"  {card:<13s} {len(d.columns):>4d} cols  "
              + ("OK" if not bad else f"⛔ {len(bad)} FAILURE(S)"))
        for b in bad[:5]:
            print(f"      {b}")
    return bad


# --------------------------------------------------------------------------- #
# STEPS 4-6 — the viewer: blocks, coverage semantics, HTML
#
# A card is 33-297 columns. Dumping them alphabetically is a spreadsheet, not a
# card; the reading order below is the EVIDENCE order the subproject argues in
# (identity -> what it is; abundance -> how much; targetome -> what it could
# touch; model -> what was fitted; realization -> what landed; attribution ->
# who gets credit; context -> what could explain it away; ... -> fame LAST).
#
# ⚠ THIS LIST WILL ROT as columns are added. It is built to rot LOUDLY: an
# unmatched prefix lands in a visible "other" block and `build()` prints the
# count, so a new column APPEARS in the wrong place rather than vanishing.
# --------------------------------------------------------------------------- #
#         key            label                                        prefixes                          exact names
BLOCKS: list[tuple[str, str, tuple, tuple]] = [
    ("identity",   "Identity & locus",
     ("arm_", "seed_", "family_", "gene_", "gctx_", "field_", "tier_", "n_"),
     ("detection", "gene_cancer_role", "arms", "n", "identified", "spiker")),
    ("admiss",     "Admissibility", ("adm_",), ()),
    ("abundance",  "Abundance & dose",
     ("abund_", "dose_", "grank_", "rank_", "healthy_", "hly_", "ago_", "aid_", "pressure_"),
     ("total_pressure", "concentration")),
    ("targetome",  "Targetome — four universes, NEVER blended",
     ("seq_", "site_", "ts_", "sdsz_", "kd_"), ()),
    ("chimeric",   "Chimeric / ligation evidence", ("chim_", "echim_"), ()),
    ("model",      "Learned model — β",
     ("beta_", "pip_", "coupling_", "cal_", "oof_", "term_", "prior_", "model_", "arb_",
      "p_", "m_", "w_", "sd_", "z_", "max_beta"),
     ("beta", "z", "identity", "m", "w", "sd", "p", "d")),
    ("realization", "Realization ladder",
     ("real_", "greal_", "realized_", "realization_"), ("is_realization_owner",)),
    ("attribution", "Attribution & share",
     ("attr_", "share_", "top_", "dominant_", "owner_", "static_", "identity_"),
     ("regulatory_handoff",)),
    ("context",    "Confounders & context",
     ("ctx_", "comp_", "cload_", "comv_", "sub_", "esub_", "retention_", "composition_",
      "surrogate_", "edge_rho"),
     ("retention", "median_retention", "n_composition_explained", "n_cell_intrinsic")),
    ("progression", "Progression & state",
     ("wshift_", "shift_", "acquired_", "acquisition_", "tcga_", "frac_", "d_", "net_",
      "wiring_", "dGlobal", "dShare", "delta_", "mean_d", "mean_own", "own_specific",
      "rho_gene_pressure"), ()),
    ("cnv",        "Copy number", ("cnv_", "cnvc_", "foc_"), ()),
    ("isomir",     "isomiR composition", ("iso_", "isoc_"), ()),
    ("protein",    "CPTAC protein", ("cptac_",), ()),
    ("outcome",    "Outcome / survival", ("surv_",), ()),
    ("famrole",    "Family role", ("famrole_", "fam_"), ()),
    ("broadcast",  "⚠ Broadcast from a coarser unit", ("bc_",), ()),
    ("literature", "Literature ground truth", ("lit_",), ()),
    ("fame",       "⛔ Fame axes — curation depth, NOT targetome", ("fame_", "cur_"), ()),
    ("coverage",   "Coverage — a blank here means UNSCANNED", ("cov_",), ()),
    ("retired",    "⛔ Retired / superseded", ("hx_",), ()),
    ("other",      "Unclassified — the block taxonomy has not caught up", (), ()),
]
COLLAPSED = {"fame", "retired", "other"}

# ⭐ The card's OWN rung — declared, because it is a NAMING fact, not a derivable one: the rung
# vocabulary calls the edge card's grain `edge`, not `arm`, and the family card carries MORE
# gene-rung columns (35) than family-rung ones (24), so "the modal rung" would answer `gene` and be
# wrong. `build()` asserts each entry is genuinely the grain (never group-compressible).
# ⛔ FIXED 2026-08-19: these keys said `"family"`, but `card_rungs.CARDS` renamed that card to
# `"gene_family"` in 6dbd321 — and THIS FILE WAS UNTRACKED at the time, so the rename could not ripple
# to it. `--build` has raised `KeyError: 'family'` ever since; the shipped HTML predates the rename.
# The exact failure mode the downstream-ripple axiom exists to catch, hidden by the file being untracked.
CARD_RUNG = {"arm": "arm", "gene": "gene", "edge": "edge",
             "gene_family": "family", "seed_family": "family"}


def _block_of(card: str, col: str) -> str:
    """Assign a column to a reading block. `seed_family` is matched on the stem BENEATH its
    universal `fam_` prefix — otherwise all 32 of its columns collapse into one undifferentiated
    block and `fam_ctx_composition` (a confounder) reads as a family-role column."""
    if col in CARDS[card]["key"]:
        return "identity"                     # the key IS the identity, whatever it is named
    for probe in ((col[4:], col) if card == "seed_family" and col.startswith("fam_") else (col,)):
        for key, _lab, prefixes, exact in BLOCKS:
            if probe in exact or any(probe.startswith(p) for p in prefixes):
                return key
    return "other"


def cov_map(d: pd.DataFrame, blocks: dict[str, list[str]], verbose: bool = False) -> dict:
    """⭐ Which `cov_*` flag governs which block — MEASURED, never curated.

    The viewer must distinguish three empty states, and two of them look identical on disk:
    a measured `0`, a blank because the scan never covered this row (**UNSCANNED**), and a blank
    despite coverage. Only the `cov_*` flags separate the last two.

    ⛔ A HAND-WRITTEN flag→block table would be a GUESS, and a wrong guess makes the viewer assert
    "not scanned" about data that was scanned — fabricating a provenance claim, which is the exact
    failure the cards exist to prevent. So the mapping is derived and must hold near-EXACTLY:
    `flag is False` ⟺ `every column in the block is blank on that row`, both directions ≥ 0.995.
    Implication alone is too weak — a rarely-covered block satisfies it against almost any flag.

    Blocks that fail simply get NO flag, and their blanks render as the honest, weaker
    "— no value" rather than a confident "— not scanned".
    """
    flags = [c for c in d.columns if c.startswith("cov_")]
    out: dict[str, dict] = {}
    if not flags:
        return out
    for bk, cols in blocks.items():
        cols = [c for c in cols if not c.startswith("cov_")]
        if not cols or bk in ("coverage", "identity"):
            continue
        blank = d[cols].isna().all(axis=1)
        if blank.all() or not blank.any():
            continue                                  # no contrast: nothing to attribute
        for f in flags:
            off = ~d[f].fillna(False).astype(bool)
            prec = float(blank[off].mean()) if off.any() else 0.0   # flag off ⇒ blank
            rec = float(off[blank].mean())                          # blank ⇒ flag off
            if prec >= 0.995 and rec >= 0.995:
                out[bk] = {"flag": f, "prec": round(prec, 4), "rec": round(rec, 4)}
                break
    if verbose:
        print(f"  coverage flags resolved for {len(out)}/{len(blocks)} block(s): "
              + ", ".join(f"{k}→{v['flag']}" for k, v in sorted(out.items())) or "  (none)")
    return out


def payload(card: str, *, sigfig: int = 4, verbose: bool = True) -> dict:
    """Everything one card page needs: values, descriptors, rung/domain labels, blocks, links."""
    reg = pd.read_csv(REGISTRY, sep="\t").query("card == @card")
    reg["domain"] = reg["domain"].fillna("")
    reg["agg_of"] = reg["agg_of"].fillna("")
    meta = {r.column: {"rung": r.rung, "agg": r.agg_of, "dom": r.domain}
            for r in reg.itertuples()}
    d = pd.read_csv(CARDS[card]["path"], sep="\t", low_memory=False)
    pay = encode_card(card, sigfig=sigfig)

    keys = CARDS[card]["key"]
    grouped_cols: dict[str, list[str]] = {}
    for c in d.columns:
        grouped_cols.setdefault(_block_of(card, c), []).append(c)
    order = [b[0] for b in BLOCKS]
    blocks = [{"k": k, "label": dict((b[0], b[1]) for b in BLOCKS)[k],
               "cols": grouped_cols[k], "collapsed": k in COLLAPSED}
              for k in order if k in grouped_cols]

    # cross-card links, resolved HERE so the page never renders a link to a blank page
    lk: dict[str, dict] = {}
    for src, col, tgt in LINKS:
        if src != card or col not in d.columns or tgt not in CARDS:
            continue
        tk = CARDS[tgt]["key"][0]
        universe = set(pd.read_csv(CARDS[tgt]["path"], sep="\t", usecols=[tk],
                                   low_memory=False)[tk].dropna().astype(str))
        vals = set(d[col].dropna().astype(str))
        miss = sorted(vals - universe)
        lk[col] = {"to": tgt, "miss": miss,
                   "variant": sorted(m for m in miss if _has_label_variant(m, universe))}
    if verbose:
        unk = len(grouped_cols.get("other", []))
        print(f"  {card:<13s} {len(d):>6,d} rows  {len(d.columns):>4d} cols  "
              f"{len(blocks):>2d} blocks  {len(pay['dist']):>3d} descriptors"
              + (f"  ⚠ {unk} unclassified" if unk else ""))
    return pay | {"meta": meta, "blocks": blocks, "links": lk,
                  "own_rung": CARD_RUNG[card],
                  "cov": cov_map(d, grouped_cols, verbose=verbose),
                  "label": [str(x) for x in d[keys[0]].astype(str)] if len(keys) == 1 else None}


def stamp_line() -> str:
    """The provenance footer the viewer will carry, so a reader can diff it against `ls -la`."""
    parts = []
    for name, s in sources().items():
        parts.append(f"{s['path'].name} {_ts(s['mtime'])} ({s['size']/1e6:.2f} MB, sha {s['sha']})"
                     if s["exists"] else f"{s['path'].name} MISSING")
    return " · ".join(parts)


def _island(obj) -> str:
    """A JSON data island. `</` is escaped to `<\\/` — valid JSON, and it makes it impossible for a
    data value to close the <script> tag early (the classic self-contained-HTML injection)."""
    import json
    return json.dumps(obj, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")


def _shell(title: str, head: str, body: str, js: str, foot: str) -> str:
    from mirna_hallmark.analyses.ops._card_assets import CSS
    return (f"<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<title>{esc(title)}</title>\n<style>{CSS}</style></head><body>\n"
            f"<div class=\"wrap\">\n{head}\n{body}\n"
            f"<footer>{foot}</footer>\n</div>\n"
            f"<button class=\"theme\" title=\"toggle theme\">◐ theme</button>\n{js}\n</body></html>\n")


def esc(s) -> str:
    import html
    return html.escape(str(s), quote=True)


def _page(card: str, pay: dict, foot: str) -> str:
    """One card page: sticky search, result rail, detail pane."""
    from mirna_hallmark.analyses.ops._card_assets import JS
    keys = ", ".join(pay["key"])
    ncav = sum(1 for m in pay["meta"].values() if "⛔" in m["dom"] or "⚠" in m["dom"])
    head = (f'<div class="crumb"><a href="index.html">← all cards</a></div>'
            f'<header><h1>{esc(card)} card</h1>'
            f'<p class="sub">One row per <b>{esc(keys)}</b>. Every value carries the rung it was '
            f'MEASURED at — a chip ⟨like this⟩ marks a value that is not a property of this row but '
            f'of a coarser unit, so it repeats across rows by construction.</p></header>'
            f'<div class="stats">'
            f'<div class="stat"><b>{pay["n"]:,}</b><span>rows</span></div>'
            f'<div class="stat"><b>{len(pay["cols"])}</b><span>columns</span></div>'
            f'<div class="stat"><b>{len(pay["blocks"])}</b><span>blocks</span></div>'
            f'<div class="stat"><b>{len(pay["dist"])}</b><span>with a distribution</span></div>'
            f'<div class="stat"><b>{ncav}</b><span>columns carrying ⛔/⚠</span></div>'
            f'<div class="stat"><b>{esc(pay["own_rung"])}</b><span>card rung</span></div></div>')
    body = ('<div class="searchbar"><input id="q" placeholder="search — / or Ctrl-K to focus, '
            '↑↓ to move, Enter to open" autocomplete="off" spellcheck="false"><span id="qn"></span></div>'
            '<div class="grid"><aside class="rail" id="rail"></aside><main id="detail"></main></div>')
    js = (f'<script type="application/json" id="payload">{_island(pay)}</script>\n'
          f'<script>{JS}</script>')
    return _shell(f"{card} card — miRNA×Hallmark", head, body, js, foot)


def _index(pays: dict, foot: str, link_rows: pd.DataFrame) -> str:
    from mirna_hallmark.analyses.ops._card_assets import INDEX_JS
    order = ["arm", "gene", "edge", "gene_family", "seed_family"]
    order = [c for c in order if c in pays]
    total = sum(p["n"] for p in pays.values())
    ncols = sum(len(p["cols"]) for p in pays.values())
    tiles = "".join(
        f'<a class="tile" href="{c}.html"><b>{esc(c)}</b>'
        f'<span>{pays[c]["n"]:,} rows · {len(pays[c]["cols"])} columns<br>'
        f'one row per {esc(", ".join(pays[c]["key"]))}</span></a>' for c in order)
    dang = ""
    if len(link_rows):
        absent = int(link_rows.distinct_absent.sum())
        variant = int(link_rows.distinct_label_variant.sum())
        dang = (f'<div class="note">⭐ <b>Cross-card links are resolved at build time.</b> '
                f'{int(link_rows.resolved.sum()):,} of {int(link_rows.declared.sum()):,} declared '
                f'links land on a real row and are clickable. The rest render as muted text, never '
                f'as a link to a blank page — <b>{absent}</b> distinct key(s) are genuinely absent '
                f'from the target card and <b>{variant}</b> are present there under a different '
                f'label form (a bare stem vs a 5p/3p-suffixed name). ⛔ Label variants are reported, '
                f'never auto-joined: silently matching them is exactly the arm conflation the card '
                f'keys exist to refuse.</div>')
    head = (f'<header><h1>miRNA×Hallmark — card browser</h1>'
            f'<p class="sub">The five-rung card system, searchable. Each card answers questions at '
            f'ONE unit; the same column name can mean different things on different cards, so the '
            f'rung travels with every value.</p></header>'
            f'<div class="stats">'
            f'<div class="stat"><b>{len(order)}</b><span>cards</span></div>'
            f'<div class="stat"><b>{total:,}</b><span>rows</span></div>'
            f'<div class="stat"><b>{ncols}</b><span>columns</span></div></div>')
    body = ('<div class="searchbar"><input id="q" placeholder="search any arm, gene, edge, family — '
            'or a column name, to find which card carries it" autocomplete="off" spellcheck="false">'
            '<span id="qn"></span></div>'
            '<div class="grid"><aside class="rail" id="rail"></aside><main>'
            f'<div class="tiles">{tiles}</div>{dang}'
            '<div class="note">⚠ <b>A blank is not a zero.</b> Three empty states render differently '
            'and mean different things: a measured <span class="mono">0</span>; '
            '<span class="mono">— not scanned</span>, where a coverage flag says the scan never '
            'covered this row; and <span class="mono">— no value</span>, scanned but empty. '
            'Reading a blank as a zero invents a result.</div>'
            '</main></div>')
    idx = {"order": order, "total": total, "ncols": ncols,
           "cards": {c: {"n": pays[c]["n"], "cols": pays[c]["cols"],
                         "keys": _labels(c, pays[c])} for c in order}}
    js = (f'<script type="application/json" id="payload">{_island(idx)}</script>\n'
          f'<script>{INDEX_JS}</script>')
    return _shell("miRNA×Hallmark — card browser", head, body, js, foot)


def _labels(card: str, pay: dict) -> list[str]:
    """The searchable label per row — the key column(s) joined, matching what the page shows."""
    d = pd.read_csv(CARDS[card]["path"], sep="\t", usecols=CARDS[card]["key"], low_memory=False)
    cols = [d[k].astype(str).fillna("") for k in CARDS[card]["key"]]
    return [" · ".join(t) for t in zip(*cols)] if len(cols) > 1 else list(cols[0])


def build(out: Path | None = None, *, sigfig: int = 4) -> Path:
    """Generate the whole viewer in ONE step — no intermediate JSON to drift from.

    ⭐ That is the point, and it is a deliberate correction of the precedent: `gen_architecture.py`
    emits JSON and its HTML is pasted by hand, so its published `architecture.html` has been measured
    running HOURS behind its own data file. A generator whose output can silently disagree with its
    input is a staleness machine; this one writes the pages it stamps.
    """
    out = out or (C.REPO_ROOT / "mirna_hallmark/docs/derived/cards")
    out.mkdir(parents=True, exist_ok=True)

    problems = check(verbose=False)
    if problems:
        print("[gen-cards] ⛔ conformance gate FAILED — refusing to build a viewer that would "
              "label data with rungs that do not describe it:")
        for p in problems:
            print(f"   {p}")
        raise SystemExit(1)

    # ⭐ assert the declared card rung really IS the grain — a compressible "own rung" would mean
    # the header chip is lying about what a row is.
    for c, r in CARD_RUNG.items():
        assert _group_key(c, r) is None, f"{c}: declared own rung {r} is group-compressible"

    print(f"[gen-cards] building → {out}")
    pays, foot = {}, _footer()
    for card in CARDS:
        if not CARDS[card]["path"].exists():
            continue
        pays[card] = payload(card, sigfig=sigfig)
    lr = links(verbose=False)
    written = []
    for card, pay in pays.items():
        p = out / f"{card}.html"
        p.write_text(_page(card, pay, foot), encoding="utf-8")
        written.append(p)
    p = out / "index.html"
    p.write_text(_index(pays, foot, lr), encoding="utf-8")
    written.append(p)

    print()
    for p in written:
        print(f"  {p.name:<20s}{p.stat().st_size/1e6:>7.2f} MB")
    print(f"  {'TOTAL':<20s}{sum(x.stat().st_size for x in written)/1e6:>7.2f} MB")
    bad = verify(out)
    print("\n[gen-cards] " + ("self-containment VERIFIED — no external reference in any page"
                              if not bad else f"⛔ {len(bad)} external reference(s): {bad[:3]}"))
    print(f"[gen-cards] open {out/'index.html'}")
    return out


_HARNESS_HEAD = """
/* minimal DOM so the page's OWN javascript runs unmodified under node */
function el(){return{textContent:'',innerHTML:'',value:'',dataset:{},style:{},
  addEventListener(){},closest(){return null},scrollIntoView(){},querySelector(){return el()}};}
const __store={payload:{textContent:__JSON__}};
const document={getElementById:id=>__store[id]||(__store[id]=el()),addEventListener(){},
  querySelector(){return el()},querySelectorAll(){return[]},activeElement:null};
const window={addEventListener(){}};
const location={hash:''};
const matchMedia=()=>({matches:false});
"""

_HARNESS_TAIL = """
/* exercise the REAL functions: decode every column, render sample rows, probe behaviour */
const __out={cols:{},rendered:{},n:P.n,behave:{}};
for(const c of P.cols) __out.cols[c]=dec(c);
for(const i of __SAMPLES__){ render(i); __out.rendered[i]=document.getElementById('detail').innerHTML; }
/* ⭐ behaviour, not just data: the row that decodes perfectly can still render unusably. */
const B=__out.behave;
B.self = rank(LBL[0]).length? LBL[rank(LBL[0])[0][2]] : null;   /* a name must find ITSELF first */
const __q = LBL[0].replace(/-(3p|5p)$/,'');                      /* and a stem must prefer its own arm */
B.stem_q = __q; B.stem_top = rank(__q).slice(0,4).map(x=>LBL[x[2]]);
/* the three empty states must be DISTINGUISHABLE somewhere in this card */
let ns=0, nv=0, chips=0, sparks=0, dang=0, blocks=0;
for(const i of __SAMPLES__){ const h=__out.rendered[i];
  ns+=(h.match(/— not scanned/g)||[]).length; nv+=(h.match(/— no value/g)||[]).length;
  chips+=(h.match(/chip rung">⟨/g)||[]).length; sparks+=(h.match(/<svg class="spark"/g)||[]).length;
  dang+=(h.match(/class="dangle"/g)||[]).length; blocks+=(h.match(/<details class="blk"/g)||[]).length; }
B.not_scanned=ns; B.no_value=nv; B.rung_chips=chips; B.sparks=sparks; B.dangling=dang; B.blocks=blocks;
require('fs').writeFileSync(__OUTPATH__, JSON.stringify(__out));
"""


def verify_js(card: str, *, out: Path | None = None, n_sample: int = 6,
              verbose: bool = True) -> list[str]:
    """⭐⭐ RUN THE PAGE'S OWN JAVASCRIPT AND COMPARE ITS DECODE TO THE SOURCE TSV.

    `roundtrip()` proves Python encodes and Python decodes agree — it cannot see a bug in the
    BROWSER's decoder, which is a hand-written mirror in a second language and is exactly the kind
    of thing that drifts. Everything a reader sees comes through that mirror, so it is the one that
    has to be right. This runs the shipped `<script>` verbatim under node against the shipped
    payload, then diffs all 187-297 decoded columns against the TSV.

    ⚠ It already earned its keep: booleans reach JS as the STRINGS "True"/"False", so the coverage
    chart's `=== true` test matched nothing and every flag rendered OFF. The page LOOKED correct.

    Returns a list of mismatches; empty means the browser sees what the TSV holds.
    """
    import json
    import shutil
    import subprocess
    import tempfile
    import numpy as np

    if not shutil.which("node"):
        if verbose:
            print("[gen-cards] ⚠ node not found — skipping the JS-decode check "
                  "(the Python round-trip still covers the encoder)")
        return []
    out = out or (C.REPO_ROOT / "mirna_hallmark/docs/derived/cards")
    page = out / f"{card}.html"
    if not page.exists():
        return [f"{card}: {page} not built"]

    import re
    from mirna_hallmark.analyses.ops._card_assets import JS
    html = page.read_text(encoding="utf-8")
    m = re.search(r'<script type="application/json" id="payload">(.*?)</script>', html, re.S)
    if not m:
        return [f"{card}: no payload island in {page.name}"]
    payload_json = m.group(1).replace("<\\/", "</")     # undo the tag-safety escape

    d = pd.read_csv(CARDS[card]["path"], sep="\t", low_memory=False)
    rng = np.random.default_rng(0)
    samples = rng.choice(len(d), size=min(n_sample, len(d)), replace=False).tolist()
    # ⭐ DELIBERATELY SAMPLE AN UNSCANNED ROW. A uniform sample kept missing the "not scanned" state
    # on cards where it is common (1,860 of 2,450 arms have no realization scan), so the state
    # rendered ZERO times across every card and read as "working" — when in fact the block was being
    # hidden entirely. A state that the test never exercises is a state nobody is defending.
    covmap = json.loads(payload_json).get("cov", {})
    expect_ns = False
    if covmap:
        flags = sorted({v["flag"] for v in covmap.values()})
        off = d.index[~d[flags].fillna(False).astype(bool).all(axis=1)].tolist()
        if off:
            samples += off[:2]
            expect_ns = True
    samples = sorted(set(int(x) for x in samples))

    with tempfile.TemporaryDirectory() as td:
        res = Path(td) / "res.json"
        script = (_HARNESS_HEAD.replace("__JSON__", json.dumps(payload_json))
                  + JS
                  + _HARNESS_TAIL.replace("__SAMPLES__", json.dumps(samples))
                                 .replace("__OUTPATH__", json.dumps(str(res))))
        src = Path(td) / "h.js"
        src.write_text(script, encoding="utf-8")
        r = subprocess.run(["node", str(src)], capture_output=True, text=True, timeout=900)
        if r.returncode != 0:
            return [f"{card}: node FAILED — {(r.stderr or r.stdout).strip().splitlines()[-1][:300]}"]
        got = json.loads(res.read_text())

    bad: list[str] = []
    for c in d.columns:
        g, s = got["cols"][c], d[c]
        gnull = np.array([x is None for x in g])
        if not (gnull == s.isna().to_numpy()).all():
            bad.append(f"{card}.{c}: JS NaN mask differs on {int((gnull != s.isna().to_numpy()).sum())} row(s)")
            continue
        if pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
            a = np.array([np.nan if x is None else float(x) for x in g], float)
            b = pd.to_numeric(s, errors="coerce").to_numpy(float)
            msk = ~np.isnan(b)
            if msk.any() and not np.allclose(a[msk], b[msk], rtol=1e-3, atol=0):
                bad.append(f"{card}.{c}: JS numeric drift")
        else:
            a = ["" if x is None else str(x) for x in g]
            b = ["" if x != x else str(x) for x in s.astype("object")]
            if a != b:
                bad.append(f"{card}.{c}: {sum(1 for x, y in zip(a, b) if x != y)} JS value(s) differ")
    # The rendered pages must be non-trivial and must not leak a raw NaN/undefined INTO A VALUE.
    # ⚠ Scoped to value cells on purpose: the first version grepped the whole page and fired on
    # `NaN` inside a registry DOMAIN string ("⛔ 83.7% single-member ⇒ shares/HHI are 1 or NaN") —
    # legitimate documentation prose. A test that flags its own subject's vocabulary is noise, and
    # noise is how a real failure gets waved through.
    cell = re.compile(r'<(?:td class="v[^"]*"|span class="val")>(.*?)</(?:td|span)>', re.S)
    for i, h in got["rendered"].items():
        if len(h) < 400:
            bad.append(f"{card}: row {i} rendered only {len(h)} chars — renderer produced nothing")
        for v in cell.findall(h):
            for token in ("undefined", "NaN", "[object Object]"):
                if token in v:
                    bad.append(f"{card}: row {i} rendered a raw '{token}' in a VALUE cell "
                               f"({v[:60]!r}) — a decode gap reached the page")
        for m in re.finditer(r'(?:width|x1|x2|y|height):?="?(-?[\d.]*(?:NaN|Infinity)[^"\s;]*)', h):
            bad.append(f"{card}: row {i} emitted a non-finite SVG/CSS length ({m.group(1)!r})")
    b = got["behave"]
    # ⭐ a name must find ITSELF as the top hit — the cheapest possible search regression test
    first = " · ".join(str(d[k].iloc[0]) for k in CARDS[card]["key"])
    if b["self"] is not None and b["self"] != first:
        bad.append(f"{card}: searching a row's own name returned {b['self']!r} first")
    # ⭐ and a stem must prefer its OWN arm. Before the name-boundary tier, `miR-21` ranked
    # hsa-miR-217 / -2110 / -2113 above hsa-miR-21-3p, because ties break on length.
    if b["stem_q"] and b["stem_top"]:
        top = b["stem_top"][0]
        if not top.lower().startswith(b["stem_q"].lower()):
            bad.append(f"{card}: stem {b['stem_q']!r} ranked {top!r} first — "
                       f"a DIFFERENT entity outranks the one asked for")
    if b["blocks"] == 0:
        bad.append(f"{card}: no blocks rendered on any sampled row")
    # ⭐ the third empty state must actually FIRE where the data says it should
    if expect_ns and b["not_scanned"] == 0:
        bad.append(f"{card}: a sampled row has a cov_ flag OFF, but '— not scanned' never "
                   f"rendered — the unscanned state is being hidden, not shown")
    if verbose:
        print(f"  {card:<13s} {len(d.columns):>4d} cols · {len(samples)} rows  "
              + ("OK  " if not bad else f"⛔ {len(bad)} FAIL  ")
              + f"[{b['blocks']} blocks · {b['sparks']} sparklines · {b['rung_chips']} rung-chips · "
                f"{b['not_scanned']} not-scanned · {b['no_value']} no-value · {b['dangling']} dangling]")
        for x in bad[:6]:
            print(f"      {x}")
    return bad


def verify(out: Path) -> list[str]:
    """⛔ The repo rule: these pages must render with no network. Grep, don't trust."""
    import re
    bad = []
    for p in sorted(out.glob("*.html")):
        t = p.read_text(encoding="utf-8")
        for pat in (r"https?://", r"<script\s+src", r"<link\s+[^>]*href", r"@import"):
            for m in re.finditer(pat, t):
                bad.append(f"{p.name}: {t[m.start():m.start()+60]!r}")
    return bad


def _footer() -> str:
    """Provenance: every input's mtime/size/sha + the generator's own commit, so a reader can diff
    the page against `ls -la` and tell whether it is describing the cards on disk today."""
    import subprocess
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=C.REPO_ROOT,
                             capture_output=True, text=True, timeout=10).stdout.strip() or "?"
    except Exception:
        sha = "?"
    return (f"generated by analyses/ops/gen_cards.py @ {sha} · sources: {esc(stamp_line())}"
            f" · ⚠ this page is a VIEW; the cards under output/learned/ are the source of truth")


def main() -> None:
    ap = argparse.ArgumentParser(description="card viewer generator (step 1: conformance gate)")
    ap.add_argument("--check", action="store_true", help="audit the cards against the registry")
    ap.add_argument("--stamp", action="store_true", help="print the source stamp and exit")
    ap.add_argument("--links", action="store_true",
                    help="resolve every cross-card link and report the dangling ones")
    ap.add_argument("--roundtrip", action="store_true",
                    help="encode + decode every card and assert it survives (lossless test)")
    ap.add_argument("--build", action="store_true", help="generate the HTML viewer")
    ap.add_argument("--verify-js", action="store_true",
                    help="run the built pages' own JS under node and diff its decode vs the TSV")
    ap.add_argument("--out", type=Path, default=None, help="output dir (default docs/derived/cards)")
    ap.add_argument("--sigfig", type=int, default=4, help="numeric rounding (default 4)")
    a = ap.parse_args()
    if a.build:
        build(a.out, sigfig=a.sigfig)
        return
    if a.verify_js:
        print("[gen-cards] running each page's OWN javascript under node, "
              "diffing its decode against the source TSV\n")
        bad = []
        for card in CARDS:
            bad += verify_js(card, out=a.out)
        print("\n[gen-cards] " + ("JS DECODE CLEAN — the browser sees exactly what the TSV holds"
                                  if not bad else f"⛔ {len(bad)} FAILURE(S)"))
        raise SystemExit(1 if bad else 0)
    if a.stamp:
        print(stamp_line())
        return
    if a.links:
        links()
        return
    if a.roundtrip:
        import json
        print(f"[gen-cards] round-trip @ {a.sigfig} sig-figs — NaN mask must be EXACT\n")
        bad, sizes = [], {}
        for card in CARDS:
            if not CARDS[card]["path"].exists():
                continue
            bad += roundtrip(card, sigfig=a.sigfig)
            sizes[card] = len(json.dumps(encode_card(card, sigfig=a.sigfig),
                                         separators=(",", ":")).encode())
        raw = sum(p.stat().st_size for p in (c["path"] for c in CARDS.values()) if p.exists())
        tot = sum(sizes.values())
        print(f"\n  {'card':<13s}{'payload':>10s}")
        for c, n in sizes.items():
            print(f"  {c:<13s}{n/1e6:>9.2f}M")
        print(f"  {'TOTAL':<13s}{tot/1e6:>9.2f}M   (source TSV {raw/1e6:.2f}M)")
        print("\n" + ("[gen-cards] ROUND-TRIP CLEAN — encoding is lossless within the declared rounding"
                      if not bad else f"[gen-cards] ⛔ {len(bad)} FAILURE(S)"))
        raise SystemExit(1 if bad else 0)
    problems = check()
    raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
