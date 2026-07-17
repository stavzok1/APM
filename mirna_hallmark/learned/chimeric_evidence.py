"""Formalized CHIMERIC (miRNA↔target ligation) evidence — the one source type that resolves the specific mature
arm (carries the full arm sequence, unlike seed-inferring AGO CLIP). Unifies multiple provenances so they can be
used SEPARATELY (per source/tissue) or POOLED:

  * manakov               — Manakov chimeric eCLIP, HEK293T (per-mature read counts, `Nunique`).
  * tarbase_qCLASH        — TarBase v9 qCLASH (Umbilical vein / Stomach / Intestine / Kidney … NO breast).
  * tarbase_CLASH         — TarBase v9 CLASH (Helwak HEK293, etc.).
  * tarbase_Chimeric fragments — TarBase v9 chimeric.

Each row: (source, tissue, mature, gene, weight). `weight` = Manakov reads, or 1 per TarBase record. Pooling is a
groupby-sum; per-source is a filter. NB (checked 2026-07-09): NO breast chimeric exists in any source, and TarBase
qCLASH shares ~52% of pairs with Manakov (not independent) — hence `sources=` and the per-source breakdown matter.

CLI: .venv/bin/python3 -m mirna_hallmark.learned.chimeric_evidence            # provenance summary
"""
from __future__ import annotations

import glob
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark.learned import seed_rarity as SR

_MANAKOV = "data/external_cache/manakov_chimeric/*.tsv"
_TARBASE = "data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz"
_TB_CHIM = {"qCLASH", "CLASH", "Chimeric fragments"}
_CACHE_PARQUET = Path("mirna_hallmark/output/learned/chimeric_evidence.parquet")
_CACHE: dict = {}


def _build() -> pd.DataFrame:
    """Unified (source, tissue, mature, gene, weight) chimeric evidence. Cached."""
    if "df" in _CACHE:
        return _CACHE["df"]
    if _CACHE_PARQUET.exists():
        _CACHE["df"] = pd.read_parquet(_CACHE_PARQUET); return _CACHE["df"]
    rows = []
    # Manakov (per-mature reads)
    man: dict = defaultdict(float)
    for f in glob.glob(_MANAKOV):
        d = pd.read_csv(f, sep="\t", usecols=["noncodingRNA", "noncodingRNA_type", "gene_name", "Nunique"])
        d = d[d["noncodingRNA_type"] == "miRNA"]
        for nc, gn, nu in zip(d["noncodingRNA"], d["gene_name"], d["Nunique"]):
            mem = str(nc).split("|")
            for a in mem:
                man[(a, gn)] += nu / len(mem)
    for (a, gn), w in man.items():
        rows.append(("manakov", "HEK293T", a, gn, w))
    # TarBase chimeric subsets (per-record weight 1, keep tissue)
    tb = pd.read_csv(_TARBASE, sep="\t", usecols=["mirna_name", "gene_name", "experimental_method", "tissue"])
    tb = tb[tb["experimental_method"].isin(_TB_CHIM)]
    agg = tb.groupby(["experimental_method", "tissue", "mirna_name", "gene_name"], observed=True).size().reset_index(name="w")
    for r in agg.itertuples(index=False):
        rows.append((f"tarbase_{r.experimental_method}", str(r.tissue), r.mirna_name, r.gene_name, float(r.w)))
    df = pd.DataFrame(rows, columns=["source", "tissue", "mature", "gene", "weight"])
    df["key"] = df["mature"].map(lambda m: SR._norm(str(m).replace("hsa-", "")))
    _CACHE_PARQUET.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(_CACHE_PARQUET)
    _CACHE["df"] = df
    return df


def evidence_matrix(gene: str, members, *, sources=None, breast_only: bool = False) -> pd.DataFrame:
    """Per-(member × source) chimeric weight for `gene` — the provenance breakdown (inspect overlap here)."""
    df = _build()
    sub = df[df["gene"] == gene]
    if sources:
        sub = sub[sub["source"].isin(sources)]
    if breast_only:
        sub = sub[sub["tissue"].str.contains("breast|mamm", case=False, na=False)]
    keys = {m: (m if str(m).startswith("hsa-") else "hsa-" + str(m)) for m in members}
    all_src = sorted(_build()["source"].unique())
    if len(sub) == 0:
        return pd.DataFrame(0.0, index=list(members), columns=all_src)
    piv = sub.groupby(["mature", "source"])["weight"].sum().unstack(fill_value=0.0)
    rows = {m: (piv.loc[keys[m]] if keys[m] in piv.index else pd.Series(0.0, index=piv.columns))
            for m in members}
    return pd.DataFrame(rows).T.fillna(0.0)


def evidence(gene: str, members, *, sources=None, breast_only: bool = False, combine: str = "mean_norm") -> pd.Series:
    """Pooled chimeric weight per family member for `gene`, OVERLAP-RESOLVED. Manakov (reads) and TarBase (record
    counts) are on different scales AND share ~52% of pairs, so a naive sum double-counts + lets Manakov dominate.
    `combine='mean_norm'` (default) normalises each source to within-member fractions then averages the sources that
    have signal → a pair in two sources is CORROBORATION (averaged), not additive inflation. 'sum'=naive (double-
    counts); 'corroboration'=# sources supporting the member. `sources` filters provenance (None = pool all)."""
    mat = evidence_matrix(gene, members, sources=sources, breast_only=breast_only)
    if mat.empty or mat.to_numpy().sum() == 0:
        return pd.Series(0.0, index=list(members), name="chim")
    if combine == "sum":
        pooled = mat.sum(axis=1)
    elif combine == "corroboration":
        pooled = (mat > 0).sum(axis=1).astype(float)
    else:                                                               # mean_norm: scale-fair, overlap = averaged
        frac = mat.div(mat.sum(axis=0).replace(0, np.nan), axis=1)      # each source → fraction within members
        active = (mat.sum(axis=0) > 0)
        pooled = frac.loc[:, active].mean(axis=1).fillna(0.0) if active.any() else mat.sum(axis=1) * 0.0
    return pooled.rename("chim")


def _report():
    df = _build()
    print("=== chimeric evidence provenance ===")
    print(df.groupby("source").agg(pairs=("gene", "size"), matures=("mature", "nunique"),
                                    genes=("gene", "nunique")).to_string())
    print("\ntissues:", sorted(df["tissue"].dropna().unique())[:20])
    print("breast rows:", int(df["tissue"].str.contains("breast|mamm", case=False, na=False).sum()))


if __name__ == "__main__":
    _report()
