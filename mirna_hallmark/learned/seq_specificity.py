"""Per-edge SEQUENCE specificity = affinity CONCENTRATION (the validated replacement for seed-rarity).

Rarity (`seed_rarity.py`) failed because it is a per-ARM property (targetome SIZE) that tracks obscurity, not
dedication, and can't even depend on which gene you're asking about. The right statistic is per-(arm,gene):

    seq_spec(arm, gene) = affinity(arm→gene) / Σ_{g' ∈ universe(arm)} affinity(arm→g')

i.e. the SHARE of the arm's total predicted repression that is aimed at this gene — an arm "designed" to hit g
concentrates its affinity on g. Affinity = TargetScan **context++** magnitude (already folds in site-type +
3′-supplementary + local AU/position context + — in the conserved-predictions file — conservation; `Predicted
relative KD` is the alt biophysical weight). `universe`:
  * **HE**     — denominator = the arm's high-evidence (curated) target set. Cleaner (drops the weak genome tail),
                 aligned to the model's gene space. Shares the curated universe with evidence-mass specificity, so
                 the two agree at ρ~0.6 in the numerators (sequence vs evidence-depth) = convergent validation.
  * **genome** — denominator = all TargetScan targets of the arm. Purer sequence (needed for UNCURATED orphans in
                 discovery), but noisier (weak-tail dominated).

Family = seed → concentration looked up per seed key (any member); a family's value = its best-aligned member.

Validated (2026-07-09): recovers textbook specialists rarity missed — miR-19b/19a & miR-26a→PTEN, miR-15/16/195→
CCND1, miR-200/429→ZEB1, with the abundant generalist miR-17 LAST; converges with curated evidence-mass
specificity at Spearman +0.49…+0.69 across 6/6 genes (`structural_identity.specificity`). See
LEARNED_MODEL_DISCOVERY_SYNTHESIS §6a.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark.learned import seed_rarity as SR   # reuse _norm (seed-key normaliser)

_TS = Path("data/external_cache/targetscan/Predicted_Targets_Context_Scores.default_predictions.txt")
_CACHE: dict = {}


_FAMILY_INFO = Path("data/miRNA/miR_Family_Info.txt")


def seed_family_map() -> dict:
    """CANONICAL hsa-mature → TargetScan **seed-family** label (== FAM.family_of; from `miR_Family_Info.txt` col
    'miR family', human=Species 9606). THE seed grouping to use anywhere seeds matter: co-seed matures (miR-17-5p/
    20-5p/93-5p/106-5p) map to ONE family, fixing the old miR-number keying (`SR._norm`) which split them into
    separate keys and forced a MAX-over-members hack that inflated multi-member families."""
    if "sfmap" in _CACHE:
        return _CACHE["sfmap"]
    fi = pd.read_csv(_FAMILY_INFO, sep="\t", usecols=["miR family", "Species ID", "MiRBase ID"], dtype=str)
    fi = fi[fi["Species ID"] == "9606"]
    _CACHE["sfmap"] = dict(zip(fi["MiRBase ID"], fi["miR family"]))
    return _CACHE["sfmap"]


def _fam_of(token: str) -> Optional[str]:
    """Map any arm/family token → its canonical seed family. Accepts hsa-miR-X, miR-X, or a family label
    (returned as-is if already a family). Uses the first member that resolves."""
    sf = seed_family_map()
    s = str(token)
    if s in _CACHE.get("_famset", set()):
        return s
    for m in s.split("/"):
        for cand in (m, "hsa-" + m if not m.startswith("hsa-") else m,
                     "hsa-miR-" + m if not m.startswith(("hsa-", "miR-", "let-")) else m):
            if cand in sf:
                return sf[cand]
    return s if s in _CACHE.get("_famset", set()) else None


def _affinity() -> pd.DataFrame:
    """Per (SEED-FAMILY, gene) context++ affinity — co-seed matures POOLED (mean of each mature's summed sites),
    keyed by the canonical seed family. `key` column = family label. Precomputes `pct` (rank within family). Cached."""
    if "aff" in _CACHE:
        return _CACHE["aff"]
    df = pd.read_csv(_TS, sep="\t", usecols=["Gene Symbol", "miRNA", "weighted context++ score"])
    df = df[df["miRNA"].str.startswith("hsa-")].copy()
    df["aff"] = (-df["weighted context++ score"]).clip(lower=0)          # ≤0 context = ≥0 repression magnitude
    df["key"] = df["miRNA"].map(seed_family_map())                       # mature → SEED FAMILY (not miR-number)
    df = df.dropna(subset=["key"])
    per_mat = df.groupby(["key", "miRNA", "Gene Symbol"], observed=True)["aff"].sum().reset_index()
    g = per_mat.groupby(["key", "Gene Symbol"], observed=True)["aff"].mean().reset_index()  # POOL co-seed matures
    g["pct"] = g.groupby("key")["aff"].rank(pct=True)
    _CACHE["aff"] = g
    _CACHE["genome_total"] = g.groupby("key")["aff"].sum()
    _CACHE["_famset"] = set(g["key"].unique())
    _CACHE["_sizes"] = g.groupby("key").size()
    return g


def _he_genes() -> dict:
    """Per SEED-FAMILY set of HE (curated) target genes — the HE-universe denominator support. Cached."""
    if "he" in _CACHE:
        return _CACHE["he"]
    from mirna_hallmark.learned.evidence import ledger as LG
    he = LG.pooled_he_edges().copy()
    he["key"] = he["miRNA"].map(lambda m: _fam_of(str(m)))
    he = he.dropna(subset=["key"])
    _CACHE["he"] = he.groupby("key")["gene"].apply(set).to_dict()
    return _CACHE["he"]


def seq_spec(gene: str, candidates, *, universe: str = "HE") -> pd.Series:
    """Per-candidate affinity-CONCENTRATION share on `gene`, at SEED-FAMILY resolution. NaN if uncovered.
    universe='HE' → denominator over the family's curated targets; 'genome' → over all TargetScan targets."""
    g = _affinity()
    sub_gene = g[g["Gene Symbol"] == gene].set_index("key")["aff"]
    gtot = _CACHE["genome_total"]
    heg = _he_genes() if universe == "HE" else None
    out = {}
    for fam in candidates:
        k = _fam_of(fam)
        if k is None or k not in sub_gene.index:
            out[fam] = np.nan; continue
        num = float(sub_gene[k])
        if universe == "HE":
            den = float(g[(g["key"] == k) & (g["Gene Symbol"].isin(heg.get(k, set())))]["aff"].sum())
        else:
            den = float(gtot.get(k, np.nan))
        out[fam] = num / den if den and den > 0 else np.nan
    return pd.Series(out, name=f"seq_spec_{universe}").reindex(candidates)


def affinity_percentile(gene: str, candidates, *, min_targetome: int = 50) -> pd.Series:
    """Per-candidate: percentile of `gene`'s affinity within the SEED-FAMILY's genome targetome ∈[0,1] — "is g
    among the family's strongest targets". The robust IDENTITY signal (scale-free; captures a hub-family-with-a-
    famous-target, e.g. miR-200→ZEB1 = its #1). Now EXACT at the seed-family unit (co-seed matures pooled in
    `_affinity`) — no MAX-over-members inflation. NaN if uncovered or the family targetome < `min_targetome`."""
    g = _affinity()
    lut = g.set_index(["key", "Gene Symbol"])["pct"]
    sizes = _CACHE["_sizes"]
    out = {}
    for fam in candidates:
        k = _fam_of(fam)
        out[fam] = (float(lut[(k, gene)]) if k is not None and int(sizes.get(k, 0)) >= min_targetome
                    and (k, gene) in lut.index else np.nan)
    return pd.Series(out, name="aff_pctile").reindex(candidates)


def affinity_percentile_kd(gene: str, candidates) -> pd.Series:
    """PER-ARM analog of `affinity_percentile`, using the validated all-site GENOME-WIDE scanMiR K_D
    (`kd.genome_affinity_pct`; beats context++ at specialist recovery, kd_fair_bench 2026-07-10) at ARM
    resolution — the identity/selection specificity context++ (per-seed-family) structurally cannot give. NaN if
    the arm is uncovered in the genome scan. Drop-in `spec` input for `priors.inclusion_prior` and `attribution_identity`."""
    from mirna_hallmark.learned import kd as KD
    lut = KD.genome_affinity_pct()
    out = {a: (float(lut[(a, gene)]) if (a, gene) in lut.index else np.nan) for a in candidates}
    return pd.Series(out, name="aff_pctile_kd").reindex(candidates)


def zscore(spec: np.ndarray) -> np.ndarray:
    """Standardize a specificity vector over covered candidates; NaN (uncovered) → 0 (neutral)."""
    s = np.asarray(spec, dtype=float)
    cov = ~np.isnan(s)
    if cov.sum() < 2 or s[cov].std() < 1e-9:
        return np.zeros_like(s)
    z = np.zeros_like(s)
    z[cov] = (s[cov] - s[cov].mean()) / s[cov].std()
    return z


if __name__ == "__main__":
    import sys
    from mirna_hallmark.learned import data as LD, families as FAM
    for gene in (sys.argv[1:] or ["PTEN", "CCND1", "ZEB1"]):
        _, X, _, w = LD.assemble_gene(gene, w_prior_source="ledger")
        X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
        he = seq_spec(gene, X.columns, universe="HE")
        gen = seq_spec(gene, X.columns, universe="genome")
        d = pd.DataFrame({"seq_HE": he.round(3), "seq_genome": (gen * 1000).round(2)}).dropna(how="all")
        print(f"\n=== {gene} — top-8 by HE affinity-concentration ===")
        print(d.sort_values("seq_HE", ascending=False).head(8).to_string())
