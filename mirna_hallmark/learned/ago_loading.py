"""Per-arm AGO LOADING — a soft (0,1] functional-abundance multiplier that discounts RISC-passenger arms.
Non-circular (external measurement, never TCGA X/Y; never coupling-derived; never the miRBase −5p/−3p label).

    loading(arm) = engagement(arm) / (engagement(5p) + engagement(3p))         # within-hairpin dominance

PRIMARY — **AGO-duplex dominance (Manakov chimeric eCLIP, HEK293T)**: engagement = # chimera reads whose miRNA
side is this arm (`data/external_cache/manakov_chimeric/*`). This MEASURES RISC loading directly (the arm is in
an AGO-crosslinked duplex), so a "star"/passenger arm that is transcribed but not loaded gets a low ratio — the
`Ndups/Nunique` give the count, no matched total needed because the ratio is within-hairpin.
FALLBACK — **abundance dominance (GTEx)**: engagement = median GTEx abundance. Arm ratios are tissue-conserved, so
this transfers; abundance-dominance ≈ loading is an approximation (usually true). Used where Manakov lacks the arm.
Uncovered by both → loading = 1.0 (no discount — never penalise what we cannot measure).

Soft multiplier: a validated-targeting passenger (e.g. miR-21-3p) keeps nonzero weight, just discounted. Enters
FUNCTIONAL ABUNDANCE downstream; orthogonal to seq-specificity (abundance side, not identity). See
[[ago-loading-arm-axis]].

CLI: .venv/bin/python3 -m mirna_hallmark.learned.ago_loading            # build cache + sanity pairs
"""
from __future__ import annotations

import glob
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

_MANAKOV = "data/external_cache/manakov_chimeric/*.tsv"
_CACHE_PARQUET = Path("mirna_hallmark/output/learned/ago_loading.parquet")
_CACHE: dict = {}


def _base(arm: str) -> str:
    return re.sub(r"-(5p|3p)$", "", str(arm))


def _manakov_engagement() -> pd.Series:
    """Total AGO-chimera reads per arm (Nunique), shared chimeras split across listed arms. Cached to parquet."""
    if "mano" in _CACHE:
        return _CACHE["mano"]
    if _CACHE_PARQUET.exists():
        s = pd.read_parquet(_CACHE_PARQUET)["ago_reads"]
        _CACHE["mano"] = s
        return s
    reads: dict = defaultdict(float)
    for f in glob.glob(_MANAKOV):
        df = pd.read_csv(f, sep="\t", usecols=["noncodingRNA", "noncodingRNA_type", "Nunique"])
        df = df[df["noncodingRNA_type"] == "miRNA"]
        for nc, nu in zip(df["noncodingRNA"], df["Nunique"]):
            arms = str(nc).split("|")
            for a in arms:
                reads[a] += nu / len(arms)
    s = pd.Series(reads, name="ago_reads").sort_values(ascending=False)
    _CACHE_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    s.to_frame().to_parquet(_CACHE_PARQUET)
    _CACHE["mano"] = s
    return s


_PAIR_PARQUET = Path("mirna_hallmark/output/learned/manakov_pair_reads.parquet")


def _pair_reads() -> pd.Series:
    """Per (SPECIFIC mature arm, gene) Manakov chimera reads — the EMPIRICAL arm resolver. Unlike loading (which
    is per-hairpin dominance), this keeps the exact mature name so co-seed members are distinguished by their
    actual ligated reads (chimeric carries the full arm sequence — what non-chimeric CLIP cannot). Cached."""
    if "pair" in _CACHE:
        return _CACHE["pair"]
    if _PAIR_PARQUET.exists():
        s = pd.read_parquet(_PAIR_PARQUET)["reads"]
        _CACHE["pair"] = s; return s
    from collections import defaultdict
    reads: dict = defaultdict(float)
    for f in glob.glob(_MANAKOV):
        df = pd.read_csv(f, sep="\t", usecols=["noncodingRNA", "noncodingRNA_type", "gene_name", "Nunique"])
        df = df[df["noncodingRNA_type"] == "miRNA"]
        for nc, gn, nu in zip(df["noncodingRNA"], df["gene_name"], df["Nunique"]):
            mem = str(nc).split("|")                                   # ambiguous paralog reads split across members
            for a in mem:
                reads[(a, gn)] += nu / len(mem)
    s = pd.Series(reads, name="reads")
    s.index = pd.MultiIndex.from_tuples(s.index, names=["arm", "gene"])
    _PAIR_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    s.to_frame().to_parquet(_PAIR_PARQUET)
    _CACHE["pair"] = s; return s


def empirical_arm(gene: str, members) -> tuple:
    """EMPIRICAL L2: which family member is actually seen in Manakov chimeras with `gene` (max reads). Returns
    (arm, reads) or (None, 0) if no chimera — caller then falls back to expression×loading."""
    pr = _pair_reads()
    best, br = None, 0.0
    for m in members:
        k = m if str(m).startswith("hsa-") else "hsa-" + str(m)
        r = float(pr.get((k, gene), 0.0))
        if r > br:
            best, br = m, r
    return best, br


def _dominance(engagement: pd.Series, min_total: float = 10.0) -> pd.Series:
    """Within-hairpin dominance ratio from an engagement Series keyed by arm name (has -5p/-3p)."""
    byhp: dict = defaultdict(dict)
    for a, v in engagement.items():
        m = re.search(r"-(5p|3p)$", str(a))
        if m:
            byhp[_base(a)][a] = v
    out = {}
    for hp, arms in byhp.items():
        tot = sum(arms.values())
        if tot < min_total:
            continue
        for a, v in arms.items():
            out[a] = v / tot
    return pd.Series(out, name="loading")


def _gtex_engagement() -> pd.Series:
    """Median GTEx abundance per arm (fallback engagement)."""
    if "gtex" in _CACHE:
        return _CACHE["gtex"]
    try:
        from mirna_hallmark import gtex_mirna_matrix as GX
        M = GX.gtex_mirna_matrix()                                # donors × arms (log2 TPM+1)
        s = (2 ** M - 1).clip(lower=0).median(axis=0)             # back to linear TPM, per-arm median
    except Exception:
        s = pd.Series(dtype=float)
    _CACHE["gtex"] = s
    return s


def loading(arms, *, floor: float = 0.0) -> pd.Series:
    """Per-arm loading ∈[floor,1] — Manakov AGO-dominance, GTEx abundance-dominance fallback, else 1.0 (no
    discount). `arms` may be family labels (a/b) — resolved per member, family value = the loaded member's max."""
    ago = _dominance(_manakov_engagement())
    gtex = _dominance(_gtex_engagement())

    def norm(a):
        return a if str(a).startswith("hsa-") else "hsa-" + str(a)

    def one(name):
        best = np.nan
        for m in str(name).split("/"):
            k = norm(m)
            v = ago.get(k, gtex.get(k, np.nan))
            if v == v and (best != best or v > best):
                best = v
        return best

    out = {a: one(a) for a in arms}
    s = pd.Series(out, name="loading")
    return s.fillna(1.0).clip(lower=floor, upper=1.0)             # uncovered → 1.0 (no discount)


def _report():
    ago = _dominance(_manakov_engagement())
    print(f"[ago_loading] Manakov: {len(ago)} arms with a paired-hairpin loading ratio")
    print(f"  distribution {ago.describe()[['min','25%','50%','75%','max']].round(2).to_dict()}")
    print(f"  passengers (loading<0.2): {int((ago<0.2).sum())}/{len(ago)}")
    print("\n  known guide/passenger pairs (guide should be high):")
    for hp in ["hsa-miR-21", "hsa-miR-17", "hsa-miR-19a", "hsa-miR-155", "hsa-miR-15a", "hsa-miR-140", "hsa-miR-30a"]:
        sub = ago[ago.index.str.startswith(hp + "-")]
        if len(sub):
            print("   " + "  ".join(f"{a}={v:.2f}" for a, v in sub.items()))


if __name__ == "__main__":
    _report()
