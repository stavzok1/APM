"""scanMiR biochemical K_D → affinity-aware occupancy κ (Design §Decision C/§G; the "occupancy gauge").

The occupancy link occ(m,g,s)=a/(a+κ) needs a per-edge half-occupancy constant κ. A *per-arm* κ (the
self-anchoring default in occupancy.transform) ignores that a high-affinity site saturates at LOWER miRNA
abundance than a weak site. scanMiR (McGeary–Bartel KdModels, 2656 human) predicts, per (arm, gene), an
aggregate **repression** (integrating every site's log-K_D, site type, and background) — a principled
per-edge affinity. We map it to κ so stronger-affinity edges reach half-occupancy sooner:

    κ(m,g) = κ0(m) · 2^(β · repression(m,g)),   repression ≤ 0  ⇒  factor ≤ 1  ⇒  lower κ for stronger edges

κ0(m) = per-arm median gated abundance (the same self-anchor as the naive default), so a *typical*-affinity
edge keeps occ=0.5 at the arm's median; high-affinity edges saturate below it. Edges with no scanMiR model
or no site fall back to κ0(m). This is the biochemically-grounded fix for the coupling loss the naive κ cost.

Scan pipeline: `learned/scanmir_scan.R` (getKdModels → findSeedMatches → aggregateMatches) over the MANE
3'UTRs of the HE-edge genes (reuses method_dev.site_ladder.utr_seed_scan._gene_3utr) and the HE-edge arms.
Data: data/external_cache/scanmir/. CLI: `python -m mirna_hallmark.learned.kd` (build + report).
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

DIR = Path("data/external_cache/scanmir")
FASTA = DIR / "hallmark_3utr.fa"
ARMS = DIR / "model_arms.txt"
TSV = DIR / "kd_repression.tsv.gz"
CACHE = DIR / "kd_repression.parquet"
RSCRIPT = Path("mirna_hallmark/learned/scanmir_scan.R")

_CACHE: dict = {}


def prep_inputs(force: bool = False) -> None:
    """Write the MANE 3'UTR FASTA (HE-edge genes) + the HE-edge arm list the R scan consumes."""
    if FASTA.exists() and ARMS.exists() and not force:
        return
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.method_dev.site_ladder.utr_seed_scan import _gene_3utr
    he = D.high_evidence_edges()
    genes = set(he["gene"].dropna().astype(str))
    arms = sorted(set(he["miRNA"].dropna().astype(str)))
    utrs = _gene_3utr(genes)
    DIR.mkdir(parents=True, exist_ok=True)
    with open(FASTA, "w") as f:
        for g, s in utrs.items():
            f.write(f">{g}\n{s}\n")
    ARMS.write_text("\n".join(arms) + "\n")
    json.dump(arms, open(DIR / "model_arms.json", "w"))


def build(*, force: bool = False) -> pd.DataFrame:
    """Per-(arm, gene) scanMiR predicted repression. Loads parquet → tsv.gz → (else) runs the R scan."""
    if CACHE.exists() and not force:
        return pd.read_parquet(CACHE)
    if not TSV.exists() or force:
        prep_inputs(force=force)
        subprocess.run(["Rscript", str(RSCRIPT)], check=True)         # ~10 min; writes TSV
    df = pd.read_csv(TSV, sep="\t")
    df["repression"] = pd.to_numeric(df["repression"], errors="coerce")
    df = df.dropna(subset=["repression"])
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(CACHE)
    return df


def affinity() -> pd.DataFrame:
    if "aff" not in _CACHE:
        _CACHE["aff"] = build()
    return _CACHE["aff"]


def genome_affinity() -> pd.DataFrame:
    """Per-(arm, gene) scanMiR predicted repression from the GENOME-WIDE scan (`scanmir_genomewide.R`) — the
    UNBIASED (not HE-restricted) biochemical support for discovery/identity. Uses **all-site** repression
    (`repression_all`) — validated 2026-07-10 to recover canonical specialists better than context++ and than
    strong-site-only (kd_fair_bench.tsv). Concats every `genomewide_kd*.tsv.gz` shard, deduped per (arm, gene).
    Column `repression` for drop-in compatibility with `affinity()`."""
    if "gaff" not in _CACHE:
        import glob
        parts = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob("data/external_cache/scanmir/genomewide_kd*.tsv.gz"))]
        df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["arm", "gene", "repression_all"])
        df["repression"] = pd.to_numeric(df["repression_all"], errors="coerce")
        _CACHE["gaff"] = df.dropna(subset=["repression"]).drop_duplicates(["arm", "gene"])[["arm", "gene", "repression"]]
    return _CACHE["gaff"]


def genome_affinity_pct() -> pd.Series:
    """Per-(arm, gene) PERCENTILE of scanMiR K_D affinity within the arm's genome-wide targetome ∈[0,1] — the
    scale-free per-arm specificity ("is g among the arm's strongest targets"). One groupby-rank over the full
    genome-wide table, cached; the validated (kd_fair_bench) per-arm identity/candidacy signal context++ can't give."""
    if "gpct" not in _CACHE:
        ga = genome_affinity().copy()
        ga["aff"] = (-ga["repression"]).clip(lower=0.0)                   # ≥0 repression strength
        ga["pct"] = ga.groupby("arm")["aff"].rank(pct=True)
        _CACHE["gpct"] = ga.set_index(["arm", "gene"])["pct"]
    return _CACHE["gpct"]


def _family_affinity() -> pd.DataFrame:
    """Per (seed-family, gene) affinity = the STRONGEST (min repression) member — for family-mode occ."""
    if "fam_aff" not in _CACHE:
        from mirna_hallmark.learned import families as FAM
        aff = affinity().copy()
        fam = FAM.family_of(pd.Index(aff["arm"].unique()))
        aff["family"] = aff["arm"].map(fam)
        g = (aff.dropna(subset=["family"]).groupby(["family", "gene"])["repression"].min()
             .reset_index().rename(columns={"family": "arm"}))       # 'arm' col = family label (union key)
        _CACHE["fam_aff"] = g
    return _CACHE["fam_aff"]


def edge_kappa(X: pd.DataFrame, gene: str, *, beta: float = 1.0, floor: float = 0.1) -> pd.Series:
    """Affinity-aware per-column κ for `gene`'s occupancy transform. Columns of X are the candidate
    regulators — **arms** (arm mode) or **seed families** (family mode); both are looked up (family =
    strongest member). κ0 = per-column median gated abundance (self-anchor, matches occupancy.transform);
    scaled by 2^(β·repression) so stronger scanMiR affinity ⇒ lower κ. Unmatched columns keep κ0."""
    a = X.clip(lower=0.0)
    k0 = a[a > 0].median().fillna(1.0).replace(0.0, 1.0)
    aff = affinity()
    rg_arm = aff.loc[aff["gene"] == gene].set_index("arm")["repression"]
    fa = _family_affinity()
    rg_fam = fa.loc[fa["gene"] == gene].set_index("arm")["repression"]
    rg = pd.concat([rg_arm, rg_fam])
    rg = rg[~rg.index.duplicated()].reindex(X.columns).fillna(0.0).clip(upper=0.0)  # ≤0 (no positive boost)
    fac = np.power(2.0, beta * rg)
    return (k0 * fac).clip(lower=floor)


def _report() -> None:
    df = build()
    print(f"[kd] scanMiR repression rows: {len(df):,} | arms {df['arm'].nunique()} × genes {df['gene'].nunique()}")
    print(f"[kd] repression (predicted logFC, <0 = repressed): median {df['repression'].median():.3f} "
          f"min {df['repression'].min():.3f} | edges ≤ −1: {int((df['repression'] <= -1).sum()):,}")
    print("[kd] strongest predicted edges:")
    top = df.sort_values("repression").head(8)[["arm", "gene", "repression", "n_8mer", "n_7mer"]]
    print(top.to_string(index=False))
    for m, g in [("hsa-miR-124-3p", "PTEN"), ("hsa-miR-200c-3p", "ZEB1"), ("hsa-miR-21-5p", "PTEN")]:
        r = df[(df["arm"] == m) & (df["gene"] == g)]
        if len(r):
            print(f"    {m}→{g}: repression={r['repression'].iloc[0]:.3f}")


if __name__ == "__main__":
    _report()
