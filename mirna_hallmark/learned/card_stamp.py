"""⭐ FRESHNESS STAMPS FOR THE CARDS — the contract the base card never had (MH-270).

> **User-directed 2026-08-19: *"add stamps to card"*.**

⛔ **THE GAP THIS CLOSES.** The delivered edge card is assembled from two layers: a **base**
(`edge_card_base.tsv`, 84 of 186 columns — the FITTED quantities, a multi-worker per-gene build over 1,549
genes) and an **annotation** layer (102 columns, seconds-cheap joins). That split is correct — the base is a
cache of the expensive thing. But it was a cache **with no freshness contract**:
`edge_card_base_provenance.tsv` records what each column *is* (column → block → estimator) and **nothing
about WHEN it was built or FROM WHICH INPUTS**, so nothing could ever invalidate it.

**That one gap produced four separate recorded defects:**
  · **MH-252** — 129 genes carry model output with ZERO edge rows (gene card 1,549 vs edge card 1,420)
  · **MH-266** — prunes that printed `✅ DROPPED` and did not take, because `_annotate` re-adds base blocks
  · **MH-269** — the renames needed post-annotation machinery for the same reason
  · a delivered card silently mixing a **2026-08-04** base with same-day annotations

⇒ a stamp records **when a card was built, from which inputs (mtime + size), and over what universe**, and
`verdict()` compares that against the inputs' CURRENT state. **`_annotate` warns loudly when it is layering
fresh annotations over a stale base** — which is the exact condition that has been silently true for 16 days.

⚠ **A stamp is a claim about INPUTS, not about correctness.** `FRESH` means nothing it was built from has
changed since; it does not mean the numbers are right. And an `UNSTAMPED` card is not necessarily stale —
it is *unknown*, which is the honest verdict and the one this module reports.

    from mirna_hallmark.learned import card_stamp as ST
    ST.write(ST.BASE_CARD, inputs=[...], universe=genes)   # at build time
    ST.report()                                            # what every card's vintage is
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import pathlib
from typing import Iterable, Optional, Sequence

OUT = pathlib.Path(__file__).resolve().parents[1] / "output" / "learned"
BASE_CARD = OUT / "edge_card_base.tsv"
#: cards a reader is expected to consume directly
DELIVERED = ("realization/edge_card.tsv", "realization/gene_card.tsv", "gene_family_card.tsv",
             "arm_card.tsv", "seed_family_card.tsv")


def _stamp_path(card: pathlib.Path) -> pathlib.Path:
    return card.with_suffix(card.suffix + ".stamp.json")


def _fingerprint(p: pathlib.Path) -> Optional[dict]:
    if not p.exists():
        return None
    st = p.stat()
    return {"path": str(p), "mtime": round(st.st_mtime, 3), "size": st.st_size,
            "mtime_h": dt.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")}


def write(card: pathlib.Path, *, inputs: Sequence[pathlib.Path] = (),
          universe: Optional[Iterable[str]] = None, note: str = "",
          backfilled: bool = False) -> pathlib.Path:
    """Record what this card was built from. Called at BUILD time, beside the card it describes."""
    u = sorted(set(universe)) if universe is not None else None
    doc = {
        "card": str(card),
        "built": dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "rows": None, "cols": None,
        "universe_n": len(u) if u is not None else None,
        "universe_sha1": hashlib.sha1("\n".join(u).encode()).hexdigest()[:12] if u else None,
        "inputs": [f for f in (_fingerprint(p) for p in inputs) if f],
        "note": note,
        # ⚠ a BACKFILLED stamp describes a card built BEFORE stamping existed. It records the card's own
        # mtime as the build time, which is the best available answer and is NOT the same as having watched
        # the build. Flagged so nobody reads it as first-hand provenance.
        "backfilled": backfilled,
    }
    if card.exists():
        st = card.stat()
        doc["card_mtime_h"] = dt.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")
        if backfilled:
            doc["built"] = doc["card_mtime_h"] + ":00"
        try:
            with card.open() as fh:
                header = fh.readline().rstrip("\n").split("\t")
                doc["cols"] = len(header)
                doc["rows"] = sum(1 for _ in fh)
        except OSError:
            pass
    p = _stamp_path(card)
    p.write_text(json.dumps(doc, indent=2) + "\n")
    return p


def read(card: pathlib.Path) -> Optional[dict]:
    p = _stamp_path(card)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def verdict(card: pathlib.Path) -> tuple[str, list[str]]:
    """`FRESH` / `STALE` / `UNVERIFIED` / `UNSTAMPED`, plus the inputs that moved.

    ⛔⛔ **A BACKFILLED STAMP CAN NEVER RETURN `FRESH`, and the first version of this wrongly did.**
    A backfill captures the inputs' fingerprints *at backfill time*, so by construction nothing has moved
    since and every card reports clean — including the 16-day-old base this module exists to catch. That is
    the same failure as a `verify` that returns hardcoded True (memory `validate-the-control-and-the-tail`).
    ⇒ **only a stamp written AT BUILD TIME can support a freshness claim.** A backfilled one reports
    `UNVERIFIED`: its inputs have not moved since the backfill, which is a much weaker statement and is said
    as such.
    """
    doc = read(card)
    if doc is None:
        return "UNSTAMPED", []
    moved = []
    for rec in doc.get("inputs", []):
        now = _fingerprint(pathlib.Path(rec["path"]))
        if now is None:
            moved.append(f"{pathlib.Path(rec['path']).name} (GONE)")
        elif now["mtime"] > rec["mtime"] + 1 or now["size"] != rec["size"]:
            moved.append(f"{pathlib.Path(rec['path']).name} {rec['mtime_h']} -> {now['mtime_h']}")
    if moved:
        return "STALE", moved
    return ("UNVERIFIED" if doc.get("backfilled") else "FRESH"), []


def mixed_vintage(card: pathlib.Path = OUT / "realization/edge_card.tsv",
                  base: pathlib.Path = BASE_CARD, *, days: float = 1.0) -> tuple[bool, float]:
    """⭐ Is the delivered card layering fresh annotations over an OLDER base? Needs no stamp at all.

    This is the condition that was silently true for 16 days, and it is visible from mtimes alone — which is
    why it is a separate check from `verdict()`. A stamp can be missing or backfilled; a file's age cannot.
    """
    if not (card.exists() and base.exists()):
        return False, 0.0
    gap = (card.stat().st_mtime - base.stat().st_mtime) / 86400.0
    return gap > days, gap


def check_base(*, loud: bool = True) -> bool:
    """⭐ Called by `_annotate`. Returns True when it is safe to layer fresh annotations over the base.

    ⚠ WARNS, does not raise. A hard failure here would block every annotate run on a repo whose base is
    legitimately older than a cosmetic input change, and the recorded cost of the silent version is a
    16-day drift nobody saw — a loud, specific warning fixes that without making the pipeline unusable.
    """
    v, moved = verdict(BASE_CARD)
    # ⭐ the mtime gap is checked FIRST and unconditionally: it needs no stamp, and it is the condition that
    # was silently true for 16 days while every other signal looked fine.
    mixed, gap = mixed_vintage()
    if mixed and loud:
        print(f"  ⛔⛔ VINTAGE MIX: the delivered edge card is {gap:.0f} DAYS newer than the base it layers "
              f"over. Its 84 base-owned columns (coupling_*, beta_*, dose_*, share_*, rank_*) are "
              f"{BASE_CARD.stat().st_mtime and dt.datetime.fromtimestamp(BASE_CARD.stat().st_mtime).strftime('%Y-%m-%d')} "
              f"vintage; the other 102 are today's.")
    if v == "FRESH" and not mixed:
        return True
    if not loud:
        return False
    if v == "UNSTAMPED":
        print(f"  ⚠ base card is UNSTAMPED ({BASE_CARD.name}) — its vintage is UNKNOWN, not verified fresh. "
              f"Run `card_stamp --backfill` to record what is knowable.")
        return False
    if v == "UNVERIFIED":
        print(f"  ⚠ base card's stamp is BACKFILLED — its inputs have not moved SINCE THE BACKFILL, which is "
              f"not a freshness claim. Only a build-time stamp can support one.")
        return False
    print(f"  ⛔⛔ BASE CARD IS STALE — layering fresh annotations over it. {len(moved)} input(s) moved:")
    for m in moved[:6]:
        print(f"       {m}")
    print("     ⇒ the delivered card MIXES vintages. Re-run `canonical_card.build()` before trusting any "
          "column it owns (coupling_*, beta_*, dose_*, share_*, rank_*).")
    return False


def report() -> int:
    """Print every card's vintage. The thing a reader should be able to see without inferring from mtimes."""
    rows = [(BASE_CARD, "BASE (fitted layer)")] + [(OUT / r, "delivered") for r in DELIVERED]
    print(f"  {'card':<34}{'verdict':<11}{'built':<20}{'rows':>7}{'cols':>6}  inputs moved")
    for p, what in rows:
        if not p.exists():
            continue
        doc = read(p) or {}
        v, moved = verdict(p)
        mark = {"FRESH": "✅", "STALE": "⛔", "UNSTAMPED": "⚠ ", "UNVERIFIED": "· "}[v]
        built = doc.get("built", "—")
        if doc.get("backfilled"):
            built += " (backfilled)"
        print(f"  {p.name:<34}{mark + v:<12}{built:<32}"
              f"{doc.get('rows') or 0:>7}{doc.get('cols') or 0:>6}  {len(moved) or ''}")
    mixed, gap = mixed_vintage()
    print()
    if mixed:
        print(f"  ⛔⛔ VINTAGE MIX: the delivered edge card is {gap:.0f} days newer than its base.")
        print(f"     84 of its 186 columns are base vintage (the FITTED ones: coupling_*, beta_*, dose_*,")
        print(f"     share_*, rank_*); the other 102 are current. Re-run `canonical_card.build()`.")
    else:
        print(f"  ✅ no vintage mix — the delivered edge card is not materially newer than its base.")
    print("\n  ⚠ `UNVERIFIED` is a BACKFILLED stamp: its inputs have not moved since the backfill, which is")
    print("     NOT a freshness claim. Only a stamp written at BUILD time can support one.")
    return 0


def main() -> int:
    import sys
    if "--backfill" in sys.argv:
        # ⚠ BACKFILL records what is KNOWABLE about cards built before stamping existed: their own mtime,
        # and the current fingerprints of the inputs they are known to read. It is explicitly flagged so it
        # is never mistaken for first-hand provenance.
        from mirna_hallmark.learned import card_stamp as _self  # noqa: F401
        base_inputs = [OUT / "readouts_arm_edges.tsv", OUT / "attribution_edges.tsv",
                       OUT / "progression_edges.tsv"]
        write(BASE_CARD, inputs=[p for p in base_inputs if p.exists()], backfilled=True,
              note="backfilled 2026-08-19 (MH-270); built before stamping existed")
        for rel in DELIVERED:
            p = OUT / rel
            if p.exists():
                write(p, inputs=[BASE_CARD] if "edge_card" in rel else [], backfilled=True,
                      note="backfilled 2026-08-19 (MH-270)")
        print("✅ stamps backfilled — verdicts below are now computable rather than unknown:\n")
    return report()


if __name__ == "__main__":
    raise SystemExit(main())
