"""Breast-context literature prior + per-target arm re-ranking.

Problem: miRTarBase high-evidence (HE) status is **pan-context** — an edge is
"Functional MTI" if *any* paper in *any* cell line / tissue showed functional
repression (HeLa, HEK293, liver, …). The raw table carries no tissue column.
The only context handle is the **PMID** of each validating study.

This module:
  1. collects the Functional-MTI PMIDs backing each HE edge,
  2. fetches PubMed records (MeSH + title + abstract) for those PMIDs,
  3. scores each PMID for **breast context** (MeSH breast-neoplasm terms OR a
     breast cell-line / mammary mention in title/abstract),
  4. aggregates to a per-edge ``breast_context_score`` ∈ [0,1]
     (= fraction of an edge's validating PMIDs that are breast-context),
  5. uses it as a **within-target re-ranker** of each gene's regulators:

        pan_weight(m→g)    = log1p(evidence_score)
        breast_weight(m→g) = log1p(evidence_score) · (1 + β · breast_context_score)

     and reports, per gene, which arm is *credited* (top-1) under each scheme —
     the edges where breast context **reassigns** the top regulator are the
     headline output (non-circular: prior is literature-only, not TCGA coupling).

This is a literature PRIOR. It does not touch TCGA expression, so it cannot be
circular with the partial-ρ coupling analyses.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_breast_context
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_breast_context --beta 1.5
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_breast_context --no-fetch   # use cache only
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from pipeline.config import PATHS

OUT_DIR = C.TISSUE_REFERENCE_DIR / "edge_breast_context"
# data/ root derived from the (absolute) miRTarBase path: .../data/miRNA/mirtar.csv
_DATA_ROOT = Path(PATHS.mirtarbase_csv).resolve().parent.parent
CACHE_DIR = _DATA_ROOT / "external_cache" / "pubmed"
CACHE_FILE = CACHE_DIR / "mirtar_he_pmid_context.json"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
OPERATOR_EMAIL = "stav.zok@mail.huji.ac.il"   # NCBI etiquette; override with --email
TOOL_NAME = "APM_mirna_hallmark"

# ---- breast-context lexicon ---- #
BREAST_MESH = {
    "Breast Neoplasms",
    "Breast",
    "Mammary Neoplasms, Human",
    "Mammary Neoplasms, Experimental",
    "Triple Negative Breast Neoplasms",
    "Carcinoma, Ductal, Breast",
    "Inflammatory Breast Neoplasms",
    "Unilateral Breast Neoplasms",
    "Breast Neoplasms, Male",
    "Carcinoma, Lobular",
    "Hereditary Breast and Ovarian Cancer Syndrome",
}
# breast cell lines + generic breast terms (matched case-insensitively on title+abstract)
BREAST_CELLLINE_RX = re.compile(
    r"\b("
    r"MCF[- ]?7|MCF7|"
    r"MDA[- ]?MB[- ]?\d{2,3}|MDAMB\d{2,3}|"
    r"T[- ]?47D|T47D|"
    r"BT[- ]?(?:474|549|20|483|459)|"
    r"SK[- ]?BR[- ]?3|SKBR3|"
    r"ZR[- ]?75(?:[- ]?1)?|ZR75|"
    r"Hs[ ]?578T|"
    r"HCC[- ]?\d{3,4}|"
    r"CAMA[- ]?1|CAL[- ]?51|CAL[- ]?120|"
    r"HMEC|MCF[- ]?10A|MCF10A|"
    r"DU4475|HBL[- ]?100|EFM[- ]?19|"
    r"UACC[- ]?812|UACC[- ]?893"
    r")\b",
    re.IGNORECASE,
)
BREAST_TERM_RX = re.compile(
    r"\bbreast (?:cancer|carcinoma|tumou?r|neoplas|adenocarcinoma|epitheli)"
    r"|\bmammary (?:gland|carcinoma|tumou?r|epitheli)"
    r"|triple[- ]negative breast",
    re.IGNORECASE,
)


# --------------------------------------------------------------------------- #
# 1. Collect PMIDs per HE edge
# --------------------------------------------------------------------------- #

def _clean_pmid(x) -> Optional[str]:
    if pd.isna(x):
        return None
    s = str(x).strip()
    if s.endswith(".0"):
        s = s[:-2]
    return s if s.isdigit() else None


def collect_he_edge_pmids() -> Tuple[Dict[Tuple[str, str], Set[str]], pd.DataFrame]:
    """Return ({(miRNA,gene): set(pmid)}, he_edges_df).

    PMIDs are drawn from Functional MTI rows (the strong functional class that
    defines HE) of the raw miRTarBase table.
    """
    raw = pd.read_csv(PATHS.mirtarbase_csv)
    raw.columns = [c.strip() for c in raw.columns]
    fmti = raw[raw["Support Type"] == "Functional MTI"].copy()
    fmti = fmti.rename(columns={"Target Gene": "gene"})
    fmti["pmid"] = fmti["References (PMID)"].map(_clean_pmid)
    fmti = fmti.dropna(subset=["pmid"])

    edges = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY) if C.MIRTAR_HALLMARK_SUMMARY.exists() \
        else pd.read_csv(C.EDGES_DIR / "mirna_hallmark_edges.tsv.gz", sep="\t")
    if "high_evidence" not in edges.columns:
        edges = pd.read_csv(C.EDGES_DIR / "mirna_hallmark_edges.tsv.gz", sep="\t")
    he = edges[edges["high_evidence"]].copy() if "high_evidence" in edges.columns else edges
    he = (he.sort_values("evidence_score", ascending=False)
            .drop_duplicates(subset=["miRNA", "gene"]))

    he_pairs = set(zip(he["miRNA"], he["gene"]))
    edge_pmids: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    sub = fmti[fmti.set_index(["miRNA", "gene"]).index.isin(he_pairs)]
    for arm, gene, pmid in zip(sub["miRNA"], sub["gene"], sub["pmid"]):
        edge_pmids[(arm, gene)].add(pmid)

    print(f"[breast_ctx] HE edges: {len(he)}  "
          f"with >=1 Functional-MTI PMID: {len(edge_pmids)}")
    all_pmids = set().union(*edge_pmids.values()) if edge_pmids else set()
    print(f"[breast_ctx] unique PMIDs: {len(all_pmids)}")
    return edge_pmids, he


# --------------------------------------------------------------------------- #
# 2. Fetch PubMed records (cached)
# --------------------------------------------------------------------------- #

def _parse_pubmed_xml(xml_text: str) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    root = ET.fromstring(xml_text)
    for art in root.findall(".//PubmedArticle"):
        pmid_el = art.find(".//PMID")
        if pmid_el is None or not pmid_el.text:
            continue
        pmid = pmid_el.text.strip()
        title_el = art.find(".//ArticleTitle")
        title = "".join(title_el.itertext()).strip() if title_el is not None else ""
        abstract = " ".join(
            "".join(a.itertext()).strip()
            for a in art.findall(".//Abstract/AbstractText")
        )
        mesh = [
            d.text.strip()
            for d in art.findall(".//MeshHeading/DescriptorName")
            if d.text
        ]
        out[pmid] = {"title": title, "abstract": abstract, "mesh": mesh}
    return out


def _efetch_batch(pmids: Sequence[str], email: str, api_key: Optional[str]) -> Dict[str, dict]:
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "tool": TOOL_NAME,
        "email": email,
    }
    if api_key:
        params["api_key"] = api_key
    url = EUTILS + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": TOOL_NAME})
    with urllib.request.urlopen(req, timeout=60) as r:
        return _parse_pubmed_xml(r.read().decode("utf-8", errors="replace"))


def fetch_pubmed_context(
    pmids: Set[str],
    *,
    do_fetch: bool = True,
    email: str = OPERATOR_EMAIL,
    api_key: Optional[str] = None,
    batch: int = 200,
) -> Dict[str, dict]:
    """Fetch + cache PubMed records for the PMID set. Returns {pmid: record}."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache: Dict[str, dict] = {}
    if CACHE_FILE.exists():
        cache = json.loads(CACHE_FILE.read_text())
        print(f"[breast_ctx] cache: {len(cache)} PMIDs on disk")

    missing = sorted(pmids - set(cache))
    if missing and do_fetch:
        rate = 0.11 if api_key else 0.34   # 10/s with key, ~3/s without
        print(f"[breast_ctx] fetching {len(missing)} PMIDs "
              f"({(len(missing)+batch-1)//batch} batches) …")
        for i in range(0, len(missing), batch):
            chunk = missing[i:i + batch]
            try:
                rec = _efetch_batch(chunk, email, api_key)
                cache.update(rec)
                # PMIDs that returned nothing → mark as empty so we don't refetch
                for p in chunk:
                    cache.setdefault(p, {"title": "", "abstract": "", "mesh": []})
            except Exception as e:  # noqa: BLE001
                print(f"  [warn] batch {i//batch} failed: {type(e).__name__} {e}")
            CACHE_FILE.write_text(json.dumps(cache))
            time.sleep(rate)
            if (i // batch) % 5 == 0:
                print(f"  …{min(i+batch, len(missing))}/{len(missing)}")
        print(f"[breast_ctx] cache now {len(cache)} PMIDs")
    elif missing:
        print(f"[breast_ctx] --no-fetch: {len(missing)} PMIDs missing, scored as non-breast")

    return cache


# --------------------------------------------------------------------------- #
# 3. Score each PMID for breast context
# --------------------------------------------------------------------------- #

def score_pmid_breast(record: dict) -> Tuple[bool, List[str]]:
    """Return (is_breast, hit_reasons)."""
    hits: List[str] = []
    mesh = record.get("mesh", []) or []
    for m in mesh:
        if m in BREAST_MESH:
            hits.append(f"mesh:{m}")
    text = f"{record.get('title','')} {record.get('abstract','')}"
    cl = BREAST_CELLLINE_RX.findall(text)
    if cl:
        hits.append("cellline:" + ",".join(sorted(set(c if isinstance(c, str) else c[0] for c in cl))[:4]))
    if BREAST_TERM_RX.search(text):
        hits.append("term:breast/mammary")
    return (len(hits) > 0), hits


# --------------------------------------------------------------------------- #
# 4. Edge-level breast-context score
# --------------------------------------------------------------------------- #

def edge_breast_scores(
    edge_pmids: Dict[Tuple[str, str], Set[str]],
    pmid_cache: Dict[str, dict],
) -> pd.DataFrame:
    pmid_breast: Dict[str, Tuple[bool, List[str]]] = {}
    for p, rec in pmid_cache.items():
        pmid_breast[p] = score_pmid_breast(rec)

    rows = []
    for (arm, gene), pmids in edge_pmids.items():
        pmids = sorted(pmids)
        n_total = len(pmids)
        n_breast = sum(1 for p in pmids if pmid_breast.get(p, (False, []))[0])
        n_fetched = sum(1 for p in pmids if p in pmid_cache and
                        (pmid_cache[p].get("mesh") or pmid_cache[p].get("title")))
        reasons = sorted({
            r for p in pmids for r in pmid_breast.get(p, (False, []))[1]
        })
        rows.append({
            "miRNA": arm, "gene": gene,
            "n_pmids": n_total,
            "n_pmids_fetched": n_fetched,
            "n_breast_pmids": n_breast,
            "breast_context_score": (n_breast / n_total) if n_total else 0.0,
            "any_breast": n_breast > 0,
            "breast_reasons": "; ".join(reasons[:6]),
        })
    df = pd.DataFrame(rows)
    print(f"[breast_ctx] edges scored: {len(df)}  "
          f"with breast evidence: {int(df['any_breast'].sum())} "
          f"({100*df['any_breast'].mean():.1f}%)")
    return df


# --------------------------------------------------------------------------- #
# 5. Per-target arm re-ranking
# --------------------------------------------------------------------------- #

def per_target_reranking(
    he: pd.DataFrame,
    score_df: pd.DataFrame,
    *,
    beta: float = 1.0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Within each gene, rank regulators by pan vs breast-weighted weight.

    The breast boost is driven by the **absolute count of breast-context
    Functional-MTI PMIDs**, not the fraction — a single breast paper out of one
    must not outrank a heavily-replicated pan-context edge, while an edge with
    many breast papers (e.g. miR-21→PTEN, 7) should be credited:

        pan_weight    = log1p(evidence_score)
        breast_weight = log1p(evidence_score) + β · log1p(n_breast_pmids)

    Returns (edge_table, gene_flip_table).
    """
    df = he[["miRNA", "gene", "evidence_score"]].merge(
        score_df[["miRNA", "gene", "breast_context_score", "n_breast_pmids",
                  "n_pmids", "breast_reasons"]],
        on=["miRNA", "gene"], how="left",
    )
    df["breast_context_score"] = df["breast_context_score"].fillna(0.0)
    df["n_breast_pmids"] = df["n_breast_pmids"].fillna(0).astype(int)
    df["n_pmids"] = df["n_pmids"].fillna(0).astype(int)

    df["pan_weight"] = np.log1p(df["evidence_score"].clip(lower=0))
    df["breast_bonus"] = beta * np.log1p(df["n_breast_pmids"])
    df["breast_weight"] = df["pan_weight"] + df["breast_bonus"]

    # within-gene ranks (1 = top regulator)
    df["pan_rank"] = df.groupby("gene")["pan_weight"].rank(
        ascending=False, method="min")
    df["breast_rank"] = df.groupby("gene")["breast_weight"].rank(
        ascending=False, method="min")
    df["rank_shift"] = df["pan_rank"] - df["breast_rank"]   # >0 = promoted

    # within-gene normalized share (how much of the gene's credit each arm gets)
    df["pan_share"] = df["pan_weight"] / df.groupby("gene")["pan_weight"].transform("sum")
    df["breast_share"] = df["breast_weight"] / df.groupby("gene")["breast_weight"].transform("sum")
    df["share_delta"] = df["breast_share"] - df["pan_share"]

    df = df.sort_values(["gene", "breast_rank"]).reset_index(drop=True)

    # gene-level top-1 reassignment
    flips = []
    for gene, g in df.groupby("gene"):
        n_reg = len(g)
        if n_reg < 2:
            continue
        pan_top = g.sort_values("pan_weight", ascending=False).iloc[0]
        brt_top = g.sort_values("breast_weight", ascending=False).iloc[0]
        flipped = pan_top["miRNA"] != brt_top["miRNA"]
        flips.append({
            "gene": gene,
            "n_regulators": n_reg,
            "pan_top_arm": pan_top["miRNA"],
            "pan_top_evidence": pan_top["evidence_score"],
            "breast_top_arm": brt_top["miRNA"],
            "breast_top_evidence": brt_top["evidence_score"],
            "breast_top_score": brt_top["breast_context_score"],
            "breast_top_n_breast_pmids": int(brt_top["n_breast_pmids"]),
            "reassigned": flipped,
            "breast_top_reasons": brt_top.get("breast_reasons", ""),
        })
    flip_df = pd.DataFrame(flips).sort_values(
        ["reassigned", "breast_top_score"], ascending=[False, False]
    ).reset_index(drop=True)

    n_flip = int(flip_df["reassigned"].sum()) if not flip_df.empty else 0
    n_multi = len(flip_df)
    print(f"[breast_ctx] multi-regulator genes: {n_multi}  "
          f"top-arm reassigned by breast context: {n_flip} "
          f"({100*n_flip/max(n_multi,1):.1f}%)")
    return df, flip_df


# --------------------------------------------------------------------------- #
# 6. Figures
# --------------------------------------------------------------------------- #

def plot_summary(
    score_df: pd.DataFrame,
    flip_df: pd.DataFrame,
    out_path: Path,
    *,
    beta: float,
    dpi: int = 150,
    n_flip_show: int = 25,
) -> None:
    import matplotlib.pyplot as plt

    flips = (flip_df[flip_df["reassigned"]]
             .sort_values("breast_top_n_breast_pmids", ascending=False)
             .head(n_flip_show)).iloc[::-1]
    n_rows = max(len(flips), 1)
    fig_h = max(6.0, 0.34 * n_rows + 2.0)
    fig, axes = plt.subplots(1, 2, figsize=(16, fig_h),
                             gridspec_kw={"width_ratios": [1, 1.5]})

    # left: how much breast functional evidence backs each edge
    ax = axes[0]
    counts = score_df["n_breast_pmids"]
    bins = np.arange(0, max(counts.max() + 2, 6))
    ax.hist(counts[counts > 0], bins=bins, color="#C77CA6", alpha=0.85,
            edgecolor="#7B4F6B", align="left")
    n_any = int(score_df["any_breast"].sum())
    ax.set_xlabel("# breast-context Functional-MTI PMIDs per edge", fontsize=10)
    ax.set_ylabel("HE edges (with ≥1 breast PMID)", fontsize=10)
    ax.set_title(
        f"Breast-specific literature backing\n"
        f"{n_any}/{len(score_df)} HE edges ({100*score_df['any_breast'].mean():.0f}%) "
        "have ≥1 breast PMID",
        fontsize=10,
    )
    ax.grid(axis="y", alpha=0.2, linestyle="--")

    # right: top reassigned genes as a horizontal bar chart
    ax = axes[1]
    if flips.empty:
        ax.text(0.5, 0.5, "No top-arm reassignments at this β",
                ha="center", va="center", fontsize=11)
        ax.axis("off")
    else:
        y = np.arange(len(flips))
        ax.barh(y, flips["breast_top_n_breast_pmids"], color="#C77CA6",
                alpha=0.85, edgecolor="#7B4F6B")
        labels = [
            f"{r['gene']}: {r['pan_top_arm'].replace('hsa-','')} "
            f"(ev{r['pan_top_evidence']:.0f}) → "
            f"{r['breast_top_arm'].replace('hsa-','')} (ev{r['breast_top_evidence']:.0f})"
            for _, r in flips.iterrows()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.8)
        ax.set_xlabel("# breast PMIDs backing the breast-credited arm", fontsize=9)
        ax.set_title(
            f"Top-arm reassignments by breast context (β={beta})\n"
            f"pan-context arm → breast-credited arm   "
            f"({int(flip_df['reassigned'].sum())} genes reassigned)",
            fontsize=10,
        )
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        ax.margins(y=0.01)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[breast_ctx] figure → {out_path}")


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(
    *,
    out_dir: Path = OUT_DIR,
    beta: float = 1.0,
    do_fetch: bool = True,
    email: str = OPERATOR_EMAIL,
    api_key: Optional[str] = None,
    dpi: int = 150,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    edge_pmids, he = collect_he_edge_pmids()
    all_pmids = set().union(*edge_pmids.values()) if edge_pmids else set()

    cache = fetch_pubmed_context(all_pmids, do_fetch=do_fetch,
                                 email=email, api_key=api_key)

    score_df = edge_breast_scores(edge_pmids, cache)
    score_df.to_csv(out_dir / "edge_breast_context_scores.tsv", sep="\t", index=False)

    edge_tbl, flip_df = per_target_reranking(he, score_df, beta=beta)
    edge_tbl.to_csv(out_dir / "per_target_arm_reranking.tsv", sep="\t", index=False)
    flip_df.to_csv(out_dir / "per_target_arm_reassignments.tsv", sep="\t", index=False)

    plot_summary(score_df, flip_df,
                 out_dir / "figures" / f"breast_context_reranking_beta{beta}.png",
                 beta=beta, dpi=dpi)

    n_fetched = sum(1 for p in all_pmids if p in cache and
                    (cache[p].get("mesh") or cache[p].get("title")))
    manifest = {
        "module": "mirna_hallmark.analyses.edge_panels.edge_breast_context",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "beta": beta,
        "n_he_edges": len(he),
        "n_edges_with_pmids": len(edge_pmids),
        "n_unique_pmids": len(all_pmids),
        "n_pmids_fetched": n_fetched,
        "n_edges_any_breast": int(score_df["any_breast"].sum()),
        "frac_edges_any_breast": float(score_df["any_breast"].mean()),
        "n_multi_regulator_genes": int((flip_df["n_regulators"] >= 2).sum()) if not flip_df.empty else 0,
        "n_top_arm_reassigned": int(flip_df["reassigned"].sum()) if not flip_df.empty else 0,
        "weighting": "breast_weight = log1p(evidence) + beta * log1p(n_breast_pmids)",
        "breast_context_score": "fraction of edge's Functional-MTI PMIDs that are breast (diagnostic); re-ranker driven by absolute n_breast_pmids",
        "context_definition": "MeSH breast-neoplasm terms OR breast cell-line / mammary mention in title+abstract",
        "circularity": "literature prior only; independent of TCGA expression",
        "outputs": {
            "edge_scores": "edge_breast_context_scores.tsv",
            "edge_reranking": "per_target_arm_reranking.tsv",
            "gene_reassignments": "per_target_arm_reassignments.tsv",
        },
        "cache_file": str(CACHE_FILE),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[breast_ctx] done — β={beta}, "
          f"{manifest['n_top_arm_reassigned']} top-arm reassignments")
    return {"score_df": score_df, "edge_tbl": edge_tbl, "flip_df": flip_df,
            "manifest": manifest}


def run_beta_sweep(
    *,
    out_dir: Path = OUT_DIR,
    betas: Sequence[float] = (0.5, 1.0, 2.0),
    do_fetch: bool = False,
    email: str = OPERATOR_EMAIL,
    api_key: Optional[str] = None,
    dpi: int = 150,
) -> pd.DataFrame:
    """Re-rank at several β, report which top-arm reassignments are robust.

    A flip that already appears at the **most conservative** β (smallest, weakest
    breast boost) and keeps the same breast-credited arm across all β is robust;
    flips that only appear at large β are β-sensitive (typically single-PMID).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    betas = sorted(betas)

    edge_pmids, he = collect_he_edge_pmids()
    all_pmids = set().union(*edge_pmids.values()) if edge_pmids else set()
    cache = fetch_pubmed_context(all_pmids, do_fetch=do_fetch,
                                 email=email, api_key=api_key)
    score_df = edge_breast_scores(edge_pmids, cache)

    flips_by_beta: Dict[float, pd.DataFrame] = {}
    for b in betas:
        _, flip_df = per_target_reranking(he, score_df, beta=b)
        flips_by_beta[b] = flip_df.set_index("gene")

    # union of genes flipped at any β
    flipped_genes = sorted(set().union(*[
        set(fd.index[fd["reassigned"]]) for fd in flips_by_beta.values()
    ]))
    rows = []
    for g in flipped_genes:
        rec: dict = {"gene": g}
        arms_when_flipped = []
        for b in betas:
            fd = flips_by_beta[b]
            is_flip = bool(fd.loc[g, "reassigned"]) if g in fd.index else False
            rec[f"flip_b{b}"] = is_flip
            arm = fd.loc[g, "breast_top_arm"] if (g in fd.index and is_flip) else None
            rec[f"breast_top_arm_b{b}"] = arm
            if arm is not None:
                arms_when_flipped.append(arm)
        # n_breast_pmids of the breast-credited arm (from the largest β where flipped)
        last_arm = next((rec[f"breast_top_arm_b{b}"] for b in reversed(betas)
                         if rec[f"breast_top_arm_b{b}"]), None)
        nb = 0
        if last_arm is not None:
            m = score_df[(score_df["miRNA"] == last_arm) & (score_df["gene"] == g)]
            nb = int(m["n_breast_pmids"].iloc[0]) if not m.empty else 0
        rec["breast_top_n_breast_pmids"] = nb
        rec["flips_at_min_beta"] = rec[f"flip_b{betas[0]}"]
        rec["arm_stable"] = (len(set(arms_when_flipped)) == 1) if arms_when_flipped else False
        rec["robust"] = bool(rec["flips_at_min_beta"] and rec["arm_stable"])
        rows.append(rec)
    sweep = pd.DataFrame(rows)
    if not sweep.empty:
        sweep = sweep.sort_values(
            ["robust", "breast_top_n_breast_pmids"], ascending=[False, False]
        ).reset_index(drop=True)
    sweep.to_csv(out_dir / "per_target_arm_reassignment_beta_sweep.tsv",
                 sep="\t", index=False)

    n_by_beta = {b: int(flips_by_beta[b]["reassigned"].sum()) for b in betas}
    n_robust = int(sweep["robust"].sum()) if not sweep.empty else 0
    print(f"[breast_ctx] β-sweep flips per β: {n_by_beta}")
    print(f"[breast_ctx] robust flips (appear at β={betas[0]} + stable arm): {n_robust}")

    # figure: flips-per-β (left) + robust flips with breast PMID counts (right)
    import matplotlib.pyplot as plt
    robust = sweep[sweep["robust"]].sort_values(
        "breast_top_n_breast_pmids", ascending=True) if not sweep.empty else sweep
    n_rows = max(len(robust), 1)
    fig, axes = plt.subplots(
        1, 2, figsize=(15, max(6.0, 0.34 * n_rows + 2.0)),
        gridspec_kw={"width_ratios": [1, 1.6]})

    ax = axes[0]
    ax.bar([str(b) for b in betas], [n_by_beta[b] for b in betas],
           color="#C77CA6", alpha=0.85, edgecolor="#7B4F6B")
    for i, b in enumerate(betas):
        ax.text(i, n_by_beta[b] + 0.5, str(n_by_beta[b]), ha="center", fontsize=9)
    ax.set_xlabel("β (breast-evidence boost)", fontsize=10)
    ax.set_ylabel("# genes with top-arm reassigned", fontsize=10)
    ax.set_title(f"Reassignment count vs β\nrobust (β≥{betas[0]} + stable arm): {n_robust}",
                 fontsize=10)
    ax.grid(axis="y", alpha=0.2, linestyle="--")

    ax = axes[1]
    if robust.empty:
        ax.text(0.5, 0.5, "No robust reassignments", ha="center", va="center")
        ax.axis("off")
    else:
        y = np.arange(len(robust))
        last_arm_col = f"breast_top_arm_b{betas[-1]}"
        ax.barh(y, robust["breast_top_n_breast_pmids"], color="#7B4F6B", alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{r['gene']} → {str(r[last_arm_col]).replace('hsa-','')}"
             for _, r in robust.iterrows()], fontsize=6.8)
        ax.set_xlabel("# breast PMIDs backing the breast-credited arm", fontsize=9)
        ax.set_title("Robust top-arm reassignments\n(flip at all β, consistent arm)",
                     fontsize=10)
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        ax.margins(y=0.01)

    fig.tight_layout()
    fig_path = out_dir / "figures" / "breast_context_beta_sweep.png"
    fig.savefig(fig_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[breast_ctx] β-sweep figure → {fig_path}")

    manifest = {
        "module": "mirna_hallmark.analyses.edge_panels.edge_breast_context (beta_sweep)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "betas": list(betas),
        "n_flips_per_beta": n_by_beta,
        "n_robust_flips": n_robust,
        "robust_definition": f"reassigned at β={betas[0]} (most conservative) AND same breast-credited arm across all β",
        "output": "per_target_arm_reassignment_beta_sweep.tsv",
        "figure": "figures/breast_context_beta_sweep.png",
    }
    (out_dir / "beta_sweep_manifest.json").write_text(json.dumps(manifest, indent=2))
    return sweep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beta", type=float, default=1.0,
                    help="breast-context boost strength in the re-ranker")
    ap.add_argument("--beta-sweep", type=str, default=None,
                    help="comma-separated β values for a robustness sweep, e.g. 0.5,1.0,2.0")
    ap.add_argument("--no-fetch", dest="do_fetch", action="store_false",
                    help="use only cached PMIDs (no network)")
    ap.add_argument("--email", type=str, default=OPERATOR_EMAIL,
                    help="contact email for NCBI E-utilities etiquette")
    ap.add_argument("--api-key", type=str, default=None,
                    help="NCBI API key (raises rate limit to 10/s)")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.beta_sweep:
        betas = [float(x) for x in args.beta_sweep.split(",")]
        run_beta_sweep(out_dir=args.out_dir, betas=betas, do_fetch=args.do_fetch,
                       email=args.email, api_key=args.api_key, dpi=args.dpi)
    else:
        run(out_dir=args.out_dir, beta=args.beta, do_fetch=args.do_fetch,
            email=args.email, api_key=args.api_key, dpi=args.dpi)


if __name__ == "__main__":
    main()
