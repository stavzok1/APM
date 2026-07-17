"""Seed-rarity specificity — the COLLINEARITY-IMMUNE resolver (LEARNED_MODEL_WHATS_NEXT §"Ground STRUCTURAL
specificity in targetome RARITY"). The coalition/soup problem is an EXPRESSION-collinearity artifact, so every
expression axis (coupling, abundance, budget) is trapped inside it. Seed rarity is SEQUENCE-only → immune to the
confound: it can adjudicate a co-expressed cluster on grounds the collinearity cannot touch. Validated: within
17~92, miR-18a N=294 ≪ miR-17-5p 1259 / miR-19a 1164 → sequence alone nominates miR-18 as PTEN's specialist,
converging with the coupling result from an orthogonal axis.

N_eff(arm) = affinity-WEIGHTED effective targetome size = Σ_genes best-site-weight (8mer=1.0, 7mer=0.5) — seed-type
aware, not a binary strong-site count. rarity = 1/log(1+N_eff), gene-normalised to O(1). Because members of a seed
FAMILY share the seed, N_eff is a property of the seed → looked up per family via any member.

SOURCE (staged):
  STAGE 1 (here): TargetScan default predictions (`Predicted_Targets_Context_Scores`) — genome-wide GENES but only
    ~321 CONSERVED matures. The arms that CAUSE soups are all conserved & covered (miR-17/19/25/103/18a), so the
    soup-resolution prior is fully served now. Non-conserved specialists (miR-616/548) return NaN → conservation
    bias for THOSE (handled by STAGE 2).
  STAGE 2 (todo): genome-wide seed-match scan over full MANE 3'UTRs (utr_seed_scan.py, currently HE-restricted) →
    uniform N for EVERY arm, removing the conservation bias. Swap `_targetscan_N` for the scan output; API stable.

Consumers: `structural_identity` (sequence specificity + "designed specialist" label) and `priors` (slab-scale
multiplier on the collinear direction). One source, two lenses.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

_TS = Path("data/external_cache/targetscan/Predicted_Targets_Context_Scores.default_predictions.txt")
_GENOME_SCAN = Path("mirna_hallmark/output/learned/seed_targetome_size.parquet")   # STAGE-2 output (preferred if present)
_STRONG = {2, 3}          # TargetScan Site Type: 3=8mer, 2=7mer-m8 (canonical strong); 1=7mer-1a
_CACHE: dict = {}


def build_genome_scan(out: Path = _GENOME_SCAN) -> pd.DataFrame:
    """STAGE-2: uniform genome-wide seed rarity — scan EVERY hsa- arm's 7-8mer seed over ALL MANE 3'UTRs and
    count target genes. Removes the conservation bias of TargetScan-default (which covers only ~321 conserved
    matures, NaN for specialists like miR-616/548). Reuses the utr_seed_scan seed-match engine. Backgroundable.
    Writes [arm, key, N_7mer_plus, N_8mer]; `targetome_size` then prefers this over TargetScan automatically."""
    from mirna_hallmark.method_dev.site_ladder import utr_seed_scan as U
    annot = pd.read_parquet(U._ANNOT, columns=["feature", "gene_name", "is_mane_select"])
    genes = set(annot[(annot["is_mane_select"]) & (annot["feature"] == "UTR")]["gene_name"].unique())
    utrs = U._gene_3utr(genes)
    arms = U._arm_seq()                                                    # hsa- matures
    print(f"[seed-rarity STAGE-2] {len(utrs)} MANE 3'UTRs × {len(arms)} hsa- arms", flush=True)
    rows = []
    for i, (arm, seq) in enumerate(arms.items()):
        if len(seq) < 8:
            continue
        n7 = n8 = 0
        for utr in utrs.values():
            _, c7, c8, _ = U._scan(utr, seq)
            n7 += c7 >= 1; n8 += c8 >= 1
        rows.append((arm, _norm(re.sub(r"^hsa-", "", arm)), n7, n8))
        if (i + 1) % 250 == 0:
            print(f"  … {i + 1}/{len(arms)} arms", flush=True)
    df = pd.DataFrame(rows, columns=["arm", "key", "N_7mer_plus", "N_8mer"])
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out)
    print(f"[seed-rarity STAGE-2] wrote {out}  ({len(df)} arms; N_strong median {int(df['N_8mer'].median())})", flush=True)
    return df


def _norm(name: str) -> str:
    """Normalise a mature name to a seed-comparable key: drop species prefix + the paralog letter right after the
    miR number (miR-20a-5p → miR-20-5p), so family members map together (they share the seed → share N)."""
    s = re.sub(r"^[a-z]{3}-", "", name)                       # hsa- / mml- / ptr- …
    s = re.sub(r"(miR-\d+)[a-z](?=-|$)", r"\1", s)            # miR-20a-5p → miR-20-5p ; miR-103a-3p → miR-103-3p
    return s


_W8, _W7 = 1.0, 0.5       # seed-type efficacy weights: 8mer=1.0, 7mer(m8/A1)=0.5 (Bartel/TargetScan efficacy tiers)


def _seed_N() -> pd.Series:
    """N_eff(seed-key) = affinity-WEIGHTED effective targetome size = Σ_genes best-site-weight (8mer=1.0, 7mer=0.5),
    NOT a binary strong-site count. Seed-type-aware (a 7mer-heavy targetome is *more specific* than an 8mer-heavy
    one at equal gene count) and STABLE across the Stage-1→Stage-2 swap (both stages use the same two-tier weight).
    Prefers the STAGE-2 genome-wide scan (uniform, no conservation bias) when built; else STAGE-1 TargetScan-default
    (conserved matures only)."""
    if "N" in _CACHE:
        return _CACHE["N"]
    if _GENOME_SCAN.exists():                                             # STAGE 2: uniform genome-wide
        g = pd.read_parquet(_GENOME_SCAN)                                 # N_8mer ⊆ N_7mer_plus (8mer gene also has m8)
        g["N_eff"] = _W8 * g["N_8mer"] + _W7 * (g["N_7mer_plus"] - g["N_8mer"])   # 8mer genes + ½·(7mer-only genes)
        _CACHE["N"] = g.groupby("key")["N_eff"].max()
        return _CACHE["N"]
    if not _TS.exists():                                                  # STAGE 1: TargetScan conserved
        _CACHE["N"] = pd.Series(dtype=float); return _CACHE["N"]
    df = pd.read_csv(_TS, sep="\t", usecols=["Gene Symbol", "miRNA", "Site Type"])
    h = df[df["miRNA"].str.startswith("hsa-") & df["Site Type"].isin(_STRONG)].copy()
    h["key"] = h["miRNA"].map(_norm)
    best = h.groupby(["key", "Gene Symbol"])["Site Type"].max()          # per-gene best site (3=8mer > 2=7mer-m8)
    _CACHE["N"] = best.map({3: _W8, 2: _W7}).groupby("key").sum()        # same two-tier weighting as STAGE 2
    return _CACHE["N"]


def targetome_size(arms) -> pd.Series:
    """N per arm/family name (uses any member's seed key). NaN where the seed is not in the conserved set (STAGE 1)."""
    N = _seed_N()
    out = {}
    for a in arms:
        for member in str(a).split("/"):                      # a family name lists co-seed members
            m = member if member.startswith("miR-") or member.startswith("let-") else "miR-" + member
            k = _norm(m)
            if k in N.index:
                out[a] = round(float(N[k]), 1); break                     # N_eff is fractional (affinity-weighted)
        else:
            out[a] = np.nan
    return pd.Series(out, name="targetome_N")


def rarity(arms) -> pd.Series:
    """Seed rarity = 1/log(1+N), min-max scaled to (0,1] over the covered arms. Rare seed (small N) → →1;
    ubiquitous seed → →0. NaN (uncovered, STAGE 1) left NaN so callers can choose a neutral fallback."""
    N = targetome_size(arms)
    r = 1.0 / np.log1p(N.astype(float))
    lo, hi = r.min(), r.max()
    return ((r - lo) / (hi - lo)) if hi > lo else r.clip(upper=1.0)


if __name__ == "__main__":
    import sys
    from mirna_hallmark.learned import structural_identity as SI
    genes = [a for a in sys.argv[1:] if not a.startswith("-")] or ["PTEN", "ZEB1"]
    for g in genes:
        s = SI.structural_identity(g)
        if s.empty:
            continue
        tab = pd.DataFrame({"targetome_N": targetome_size(s.index), "rarity": rarity(s.index).round(2),
                            "ev_mass_spec": s["specificity"]})
        tab = tab.sort_values("rarity", ascending=False)
        print(f"\n=== {g} — seed RARITY (sequence) vs evidence-mass specificity ===")
        print(tab.to_string())
