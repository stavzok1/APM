"""The PMID-deduped, method-centric evidence ledger (Design §Decision D/E; BUILD_PLAN §1b).

Builds ONE long ledger keyed by (arm, gene, pmid, assay_class) from the two overlapping PMID-bearing
curated sources — **miRTarBase** (`References (PMID)`, `Experiments`) and **TarBase v9**
(`article_pubmed_id`, `experimental_method`) — then **deduplicates globally by (edge × PMID × assay_class)**
so a paper reported by both DBs counts once. The fused per-edge weight is method-centric:

    ledger_weight(m,g) = Σ_class  CLASS_WEIGHT[class] · log1p(#distinct PMIDs in that class)

This replaces the hand-weighted, source-additive, un-deduped `evidence_scoring.score_confidence_logclass`
(+ ENCORI boost). Non-Functional (validated non-target) miRTarBase rows are excluded. Sequence/K_D
(PMID-free) is a SEPARATE channel (not here). ENCORI/POSTAR3/Manakov (CLIP/chimeric, mostly PMID-poor
in-cache) are folded into the same assay classes in a later increment (union, never summed).

Scope: Hallmark edge-universe genes (what the model uses). Cache: output/learned/evidence_ledger.parquet.
CLI: `python -m mirna_hallmark.learned.evidence.ledger` builds + reports the dedup impact.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned.evidence.methods import CLASS_WEIGHT, classify

MIRTAR = Path("data/miRNA/mirtar.csv")
TARBASE = Path("data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz")
CACHE = Path("mirna_hallmark/output/learned/evidence_ledger.parquet")
WEIGHTS_CACHE = Path("mirna_hallmark/output/learned/edge_weights.parquet")

_LT_FUNC = {"reporter", "western", "proteomics"}   # LOW-THROUGHPUT functional (luciferase / blot / MS); excludes
                                                   # qpcr_rna, which lumps in high-throughput RNA-seq / microarray


def pooled_he_edges() -> pd.DataFrame:
    """(memoized — see `_MEM`) POOLED-HE regulator edges = miRTarBase-HE ∪ TarBase-v9 low-throughput functional (reporter/western/
    proteomics, non-weak). ≈5,940 = HE 5,300 + ~640 TarBase — the multi-source analog of
    `data_loaders.high_evidence_edges`. Columns [miRNA, gene, evidence_score]; score carried from miRTarBase-HE,
    NaN for TarBase-only edges (their adaptive prior comes from the ledger weight, w_prior_source='ledger').
    NOT the 43k 'any-functional' set (qpcr_rna swallows high-throughput RNA-seq/microarray), NOT prediction/CLIP."""
    if "pooled_he" in _MEM:
        return _MEM["pooled_he"].copy()
    he = D.high_evidence_edges()[["miRNA", "gene", "evidence_score"]].drop_duplicates(["miRNA", "gene"])
    lg = build_ledger()
    tb = lg.loc[(lg["source"] == "tarbase") & lg["assay_class"].isin(_LT_FUNC) & (~lg["weak"]),
                ["arm", "gene"]].drop_duplicates().rename(columns={"arm": "miRNA"})
    tb["evidence_score"] = np.nan
    pooled = _clean_edge_arms(pd.concat([he, tb], ignore_index=True))
    _MEM["pooled_he"] = pooled.drop_duplicates(["miRNA", "gene"], keep="first").reset_index(drop=True)
    return _MEM["pooled_he"].copy()


def _clean_edge_arms(df: pd.DataFrame) -> pd.DataFrame:
    """Source-level hygiene on curated arm names (parsing corruption in miRTarBase/TarBase rows): normalize unicode
    dashes (U+2010–2015 / U+2212 → ASCII '-'; caught `hsa-miR‑15a`), SPLIT comma-joined multi-arm rows (caught
    `hsa-miR-27a-3p,hsa-miR-30b-5p`), strip whitespace, keep only `hsa-` arms. (A mouse-style alias `hsa-miR-21a-5p`
    survives here — it is a curation, not a parsing, artifact — but is harmlessly dropped downstream: not in `X`.)"""
    df = df.copy()
    df["miRNA"] = (df["miRNA"].astype(str).str.strip()
                   .str.replace(r"[‐-―−]", "-", regex=True))
    df["miRNA"] = df["miRNA"].str.split(",")                        # explode comma-joined multi-arm entries
    df = df.explode("miRNA")
    df["miRNA"] = df["miRNA"].str.strip()
    return df[df["miRNA"].str.startswith("hsa-")]


def _hallmark_genes() -> set[str]:
    return set(D.load_hallmark_edges()["gene"].dropna().astype(str))


def _classmap(methods) -> dict:
    return {m: classify(m) for m in set(methods)}


def _ingest_mirtarbase(genes: set[str]) -> pd.DataFrame:
    cols = ["miRNA", "Target Gene", "Experiments", "Support Type", "References (PMID)"]
    keep = []
    for ch in pd.read_csv(MIRTAR, usecols=cols, dtype=str, chunksize=500_000, low_memory=False):
        ch = ch[ch["miRNA"].str.startswith("hsa", na=False) & ch["Target Gene"].isin(genes)]
        ch = ch[~ch["Support Type"].fillna("").str.startswith("Non-Functional")]  # drop validated non-targets
        if len(ch):
            keep.append(ch)
    if not keep:
        return pd.DataFrame(columns=["arm", "gene", "pmid", "assay_class", "source", "weak"])
    df = pd.concat(keep, ignore_index=True)
    df["Experiments"] = df["Experiments"].fillna("").str.split("//")
    df = df.explode("Experiments")
    cm = _classmap(df["Experiments"])
    df["assay_class"] = df["Experiments"].map(cm)
    df["pmid"] = pd.to_numeric(df["References (PMID)"], errors="coerce")
    df["weak"] = df["Support Type"].fillna("").str.contains("Weak")
    df = df.rename(columns={"miRNA": "arm", "Target Gene": "gene"})
    df["source"] = "mirtarbase"
    return df[["arm", "gene", "pmid", "assay_class", "source", "weak"]]


def _ingest_tarbase(genes: set[str]) -> pd.DataFrame:
    cols = ["mirna_name", "gene_name", "experimental_method", "regulation", "article_pubmed_id"]
    keep = []
    for ch in pd.read_csv(TARBASE, sep="\t", usecols=cols, dtype=str, chunksize=1_000_000):
        ch = ch[ch["gene_name"].isin(genes)]
        if len(ch):
            keep.append(ch)
    if not keep:
        return pd.DataFrame(columns=["arm", "gene", "pmid", "assay_class", "source", "weak"])
    df = pd.concat(keep, ignore_index=True)
    cm = _classmap(df["experimental_method"])
    df["assay_class"] = df["experimental_method"].map(cm)
    df["pmid"] = pd.to_numeric(df["article_pubmed_id"], errors="coerce")
    df["weak"] = df["regulation"].fillna("") == "Positive"  # (near-degenerate; Positives flagged elsewhere)
    df = df.rename(columns={"mirna_name": "arm", "gene_name": "gene"})
    df["source"] = "tarbase"
    return df[["arm", "gene", "pmid", "assay_class", "source", "weak"]]


_MEM: dict = {}          # in-memory memo. build_ledger/edge_weights/pooled_he_edges are PURE given the
                         # on-disk caches, but were re-reading parquet on EVERY call — and `assemble_gene`
                         # calls them per gene, so a 1.5k-gene sweep paid ~2 parquet reads per gene (72% of
                         # its runtime, cProfile). Callers get a COPY, so no shared-mutation ripple.


def build_ledger(*, force: bool = False) -> pd.DataFrame:
    if not force and "ledger" in _MEM:
        return _MEM["ledger"].copy()
    if CACHE.exists() and not force:
        _MEM["ledger"] = pd.read_parquet(CACHE)
        return _MEM["ledger"].copy()
    genes = _hallmark_genes()
    mt = _ingest_mirtarbase(genes)
    tb = _ingest_tarbase(genes)
    raw = pd.concat([mt, tb], ignore_index=True)
    raw = raw.dropna(subset=["pmid"])
    raw["pmid"] = raw["pmid"].astype("int64")
    # GLOBAL dedup: one (edge × PMID × assay_class) counts once, across both DBs and both tiers.
    ledger = raw.drop_duplicates(["arm", "gene", "pmid", "assay_class"]).reset_index(drop=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    ledger.to_parquet(CACHE)
    return ledger


def edge_weights(ledger: Optional[pd.DataFrame] = None, *, force: bool = False,
                 weights: Optional[dict] = None) -> pd.DataFrame:
    """Fused per-edge weight Σ_class w_class·log1p(#distinct PMIDs). `weights` overrides CLASS_WEIGHT (e.g.
    the task-matched CLASS_WEIGHT_MRNA from calibrate) — an override recomputes and does NOT touch the cache."""
    if weights is None and ledger is None and not force:
        if "weights" in _MEM:
            return _MEM["weights"].copy()
        if WEIGHTS_CACHE.exists():
            _MEM["weights"] = pd.read_parquet(WEIGHTS_CACHE)
            return _MEM["weights"].copy()
    if ledger is None:
        ledger = build_ledger(force=force)
    wmap = weights if weights is not None else CLASS_WEIGHT
    per = (ledger.groupby(["arm", "gene", "assay_class"])["pmid"].nunique()
           .reset_index(name="n_pmid"))
    per["contrib"] = per["assay_class"].map(wmap).fillna(0.0) * np.log1p(per["n_pmid"])
    w = per.groupby(["arm", "gene"])["contrib"].sum().reset_index(name="ledger_weight")
    if weights is None:                                    # only the canonical weights are cached
        WEIGHTS_CACHE.parent.mkdir(parents=True, exist_ok=True)
        w.to_parquet(WEIGHTS_CACHE)
    return w


def _report() -> None:
    genes = _hallmark_genes()
    print(f"[ledger] Hallmark edge-universe genes: {len(genes)}")
    mt = _ingest_mirtarbase(genes)
    tb = _ingest_tarbase(genes)
    raw = pd.concat([mt, tb], ignore_index=True).dropna(subset=["pmid"])
    raw["pmid"] = raw["pmid"].astype("int64")
    ledger = raw.drop_duplicates(["arm", "gene", "pmid", "assay_class"]).reset_index(drop=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    ledger.to_parquet(CACHE)
    print(f"[ledger] raw (m,g,pmid,class) rows: {len(raw):,}  ->  deduped: {len(ledger):,} "
          f"({100*(1-len(ledger)/len(raw)):.1f}% collapsed)")
    print(f"[ledger] distinct PMIDs: {ledger['pmid'].nunique():,} | distinct edges: "
          f"{ledger.groupby(['arm','gene']).ngroups:,}")
    print("[ledger] rows by assay class (deduped):")
    print(ledger["assay_class"].value_counts().to_string())
    # dedup impact concentrated in CLIP: raw vs deduped for ago_clip
    for cls in ["ago_clip", "reporter", "western", "qpcr_rna", "chimeric"]:
        r = (raw["assay_class"] == cls).sum(); d = (ledger["assay_class"] == cls).sum()
        if r:
            print(f"    {cls:10s} raw {r:>9,} -> dedup {d:>9,}  ({100*(1-d/r):.1f}% collapsed)")
    w = edge_weights(ledger, force=True)
    print(f"[ledger] edge_weights rows: {len(w):,}  (cached {WEIGHTS_CACHE})")
    print(w.sort_values("ledger_weight", ascending=False).head(5).to_string(index=False))


if __name__ == "__main__":
    _report()
