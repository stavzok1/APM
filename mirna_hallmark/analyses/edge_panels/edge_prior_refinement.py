"""Multi-axis edge-prior refinement — transparent, auditable, honest.

The landscape spine weights each edge by `log1p(evidence_score)` only (literature
study count). That is one axis of "how much should we believe m→g". This module
layers the **other defensible, non-circular axes** on top, keeping every axis as a
**separate visible term** so nothing is hidden, and attributes every per-target
re-ranking to the axis that caused it.

Axes (each a log-additive term; default weight 1.0, all tunable):

  base       = log1p(evidence_score)            literature study count (the spine anchor)
  + w_ts  ·  log1p(ts_weight)                   A. sequence site strength
                                                   (TargetScan Σ|weighted context++|;
                                                    "weighted" already folds conservation/PCT)
  + w_brst·  log1p(n_breast_pmids)              breast-context literature (from edge_breast_context)
  +          log(reliability_factor)            B. negative-evidence penalty (≤0)
                                                   reliability = n_func / (n_func + λ·n_nonfunc)
  +          log(arm_factor)                     C. guide/passenger penalty (≤0)
                                                   passenger arm → γ (<1), else neutral 1.0

  refined_weight = sum of the above.

**Axis C — guide vs passenger (external, miRGeneDB).** Fetched live from
mirgenedb.org (cached `data/external_cache/mirgenedb/`): the `mat=1` set is the
**guide** (functional/dominant) arms, `star=1` is the **passenger/star** arms.
miRGeneDB uses paralog-cluster names, so arms are mapped to miRBase names by
**exact mature-sequence match** against the repo `mature.fa` — sequence-defined,
not name-guessed, and **non-circular** (miRGeneDB curation, never TCGA arm ratios).
Only arms with *positive passenger evidence* are down-weighted; arms absent from
miRGeneDB's high-confidence set are left **neutral, not penalised**.

**Axis D — APA / 3'UTR shortening (DIAGNOSTIC, three sources, best→weakest).**
  (1) **Breast TUMOUR** (headline): Xia 2014 ΔPDUI, TCGA breast tumour vs matched
      normal (`data/miRNA/XIA_APA.xlsx`). `brca_dpdui` < 0 = 3'UTR **shortened in
      breast tumour** ⇒ distal miRNA sites excised in tumour. Directional, tumour-
      specific, significance-filtered (427 events, 106 HE targets). Genes absent =
      not recurrent-APA in breast tumour.
  (2) **Breast NORMAL**: APAatlas GTEx-breast PDUI (218 samples) — baseline 3'UTR
      length in breast tissue (`breast_pdui_median`, `breast_apa_shortened`).
  (3) **Pan-tissue potential**: PolyASite 2.0 terminal-exon poly-A cluster count
      (`n_polyA_TE`, `apa_shortenable`) — near-universal, weak discriminator.
All three are **gene-level**, so they apply equally to every regulator of a gene and
therefore **cancel in the within-gene arm re-ranking** — supplied as per-edge
diagnostic flags, NOT weight terms. The arm-discriminating, *site-resolved* version
(does THIS arm's TargetScan site, `UTR_start`/`UTR_end`, fall distal to a proximal
poly-A site and get excised) needs transcript→genomic coordinate conversion and is
the documented next step.

**Deliberately NOT included (honesty):**
  - any TCGA-expression-derived term — that is the outcome; a prior must not see it
    (this is why guide/passenger comes from miRGeneDB sequence curation, not the
    TCGA 5p/3p read ratio).

**What this is / is not:** a per-target NOMINATION / AUDIT overlay (which arm to
credit for a gene). It is NOT the landscape budget — the spine still runs on the
pan-context S1 weight. No external breast-miRNA cohort exists to validate utility
(METABRIC has no miRNA; CPTAC miRNA = TCGA-reprocessed), so no predictive claim is made.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_prior_refinement
  .venv/bin/python3 -m mirna_hallmark.analyses.edge_panels.edge_prior_refinement --w-ts 1.0 --w-breast 1.0 --lambda-neg 1.0
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.edge_panels.edge_breast_context import (
    collect_he_edge_pmids,
    edge_breast_scores,
    fetch_pubmed_context,
)
from mirna_hallmark.robustness_checks import _load_targetscan_weights
from pipeline.config import PATHS

OUT_DIR = C.TISSUE_REFERENCE_DIR / "edge_prior_refinement"
CORR_TABLE = (C.TISSUE_REFERENCE_DIR / "edge_partial_corr_panels"
              / "edge_partial_corr_panel_corr_table.tsv")

_DATA_ROOT = Path(PATHS.mirtarbase_csv).resolve().parent.parent
MIRGENEDB_DIR = _DATA_ROOT / "external_cache" / "mirgenedb"
MIRGENEDB_GUIDE_URL = "https://mirgenedb.org/fasta/hsa?mat=1"
MIRGENEDB_STAR_URL = "https://mirgenedb.org/fasta/hsa?star=1"
MATURE_FA = _DATA_ROOT / "miRNA" / "mature.fa"
GAMMA_PASSENGER = 0.5   # down-weight factor for a confirmed passenger arm

AXIS_COLORS = {
    "sequence (TargetScan)":     "#4E79A7",
    "breast literature":         "#C77CA6",
    "negative evidence":         "#E15759",
    "guide/passenger (miRGeneDB)": "#59A14F",
    "base only (tie)":           "#999999",
}


# --------------------------------------------------------------------------- #
# Axis C — guide vs passenger arm (miRGeneDB, sequence-matched)
# --------------------------------------------------------------------------- #

def _parse_fasta(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    name = None
    seq = ""
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = seq
            name = line[1:].strip()
            seq = ""
        else:
            seq += line.strip().upper().replace("U", "T")
    if name is not None:
        out[name] = seq
    return out


def _fetch_mirgenedb(do_fetch: bool = True) -> Tuple[Dict[str, str], Dict[str, str]]:
    """(guide_arms, star_arms) as {mirgenedb_name: DNA_seq}, cached on disk."""
    MIRGENEDB_DIR.mkdir(parents=True, exist_ok=True)
    gp = MIRGENEDB_DIR / "hsa_guide_mat.fas"
    sp = MIRGENEDB_DIR / "hsa_star.fas"

    def _get(url: str, path: Path) -> str:
        if path.exists():
            return path.read_text()
        if not do_fetch:
            return ""
        import urllib.request
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 APM"})
        txt = urllib.request.urlopen(req, timeout=30).read().decode("utf-8", "replace")
        path.write_text(txt)
        return txt

    return _parse_fasta(_get(MIRGENEDB_GUIDE_URL, gp)), _parse_fasta(_get(MIRGENEDB_STAR_URL, sp))


def _arm_end(name: str) -> Optional[str]:
    if name.endswith("-5p") or "_5p" in name:
        return "5p"
    if name.endswith("-3p") or "_3p" in name:
        return "3p"
    return None


def _seed(seq: str) -> Optional[str]:
    return seq[1:8] if len(seq) >= 8 else None


def mirgenedb_arm_status(do_fetch: bool = True) -> pd.DataFrame:
    """Per miRBase hsa arm: guide / passenger / ambiguous (else absent → neutral upstream).

    Maximises coverage with a three-tier, non-circular mapping of miRGeneDB
    (guide=`mat=1`, passenger=`star=1`) onto miRBase arms:
      1. exact mature-sequence match;
      2. **seed (nt 2-8) + arm-end (5p/3p)** match — robust to miRGeneDB's end
         re-trimming, which exact matching misses;
      3. **hairpin propagation** — within a hairpin (one guide + one star arm),
         an unlabelled arm whose partner is unambiguously classified takes the
         opposite label.
    `arm_class_source` records which tier assigned each arm. Arms still absent
    from miRGeneDB's high-confidence set stay unknown (neutral, never penalised).
    """
    guide, star = _fetch_mirgenedb(do_fetch=do_fetch)
    if not guide and not star:
        print("[prior] WARNING: miRGeneDB unavailable and no cache — axis C neutral for all")
        return pd.DataFrame(columns=["miRNA", "arm_class", "arm_class_source"])

    guide_seq = set(guide.values())
    star_seq = set(star.values())
    guide_se = {(_seed(s), _arm_end(n)) for n, s in guide.items() if _seed(s) and _arm_end(n)}
    star_se = {(_seed(s), _arm_end(n)) for n, s in star.items() if _seed(s) and _arm_end(n)}

    # miRBase hsa mature arms
    rows = []
    name = None
    with open(MATURE_FA) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif name and name.startswith("hsa-"):
                seq = line.upper().replace("U", "T")
                end = _arm_end(name)
                stem = re.sub(r"-(5p|3p)$", "", name)
                rows.append({"miRNA": name, "seq": seq, "end": end,
                             "seed": _seed(seq), "stem": stem})
    mb = pd.DataFrame(rows).drop_duplicates("miRNA")

    def _classify(r) -> Tuple[str, str]:
        if r.seq in guide_seq and r.seq in star_seq:
            return "ambiguous", "seq"
        if r.seq in guide_seq:
            return "guide", "seq"
        if r.seq in star_seq:
            return "passenger", "seq"
        k = (r.seed, r.end)
        ing, ins = k in guide_se, k in star_se
        if ing and ins:
            return "ambiguous", "seed"
        if ing:
            return "guide", "seed"
        if ins:
            return "passenger", "seed"
        return "unknown", "none"

    cs = mb.apply(lambda r: pd.Series(_classify(r), index=["arm_class", "arm_class_source"]), axis=1)
    mb = pd.concat([mb, cs], axis=1)

    # tier 3: hairpin propagation
    stem_defs = (mb[mb.arm_class.isin(["guide", "passenger"])]
                 .groupby("stem")["arm_class"].apply(set))
    n_prop = 0
    for i, r in mb.iterrows():
        if r.arm_class != "unknown":
            continue
        others = stem_defs.get(r.stem, set())
        if others == {"guide"}:
            mb.at[i, "arm_class"], mb.at[i, "arm_class_source"] = "passenger", "propagated"
            n_prop += 1
        elif others == {"passenger"}:
            mb.at[i, "arm_class"], mb.at[i, "arm_class_source"] = "guide", "propagated"
            n_prop += 1

    out = mb.loc[mb.arm_class.isin(["guide", "passenger", "ambiguous"]),
                 ["miRNA", "arm_class", "arm_class_source"]].reset_index(drop=True)

    # legacy **bare names** in the edge table (e.g. 'hsa-miR-375', no -3p/-5p) denote
    # the major/guide arm — map each guide arm's stem to a bare-name 'guide' row so
    # these resolve instead of falling through to 'unknown' (false exoneration).
    have = set(out["miRNA"])
    bare = []
    for arm in out.loc[out.arm_class == "guide", "miRNA"]:
        stem = re.sub(r"-(5p|3p)$", "", arm)
        if stem != arm and stem not in have:
            bare.append({"miRNA": stem, "arm_class": "guide",
                         "arm_class_source": "bare_name_guide"})
            have.add(stem)
    if bare:
        out = pd.concat([out, pd.DataFrame(bare)], ignore_index=True)

    vc = out.arm_class.value_counts()
    print(f"[prior] miRGeneDB axis C: guide={int(vc.get('guide',0))} "
          f"passenger={int(vc.get('passenger',0))} ambiguous={int(vc.get('ambiguous',0))} "
          f"(seq+seed+propagated {n_prop}); arms unclassified stay neutral")
    return out


# --------------------------------------------------------------------------- #
# Axis B — negative (contradicting) evidence
# --------------------------------------------------------------------------- #

def negative_evidence_counts(include_weak: bool = False) -> pd.DataFrame:
    """Per (miRNA, gene): n_func (Functional MTI PMIDs) and n_nonfunc (Non-Functional)."""
    raw = pd.read_csv(PATHS.mirtarbase_csv)
    raw.columns = [c.strip() for c in raw.columns]
    raw = raw.rename(columns={"Target Gene": "gene"})

    func = raw[raw["Support Type"] == "Functional MTI"]
    if include_weak:
        nonfunc = raw[raw["Support Type"].str.startswith("Non-Functional", na=False)]
    else:
        nonfunc = raw[raw["Support Type"] == "Non-Functional MTI"]

    nf = func.groupby(["miRNA", "gene"])["References (PMID)"].nunique().rename("n_func")
    nn = nonfunc.groupby(["miRNA", "gene"])["References (PMID)"].nunique().rename("n_nonfunc")
    out = pd.concat([nf, nn], axis=1).fillna(0).astype(int).reset_index()
    return out


# --------------------------------------------------------------------------- #
# Axis D (diagnostic) — target 3'UTR APA potential from PolyASite mass-statistics
# --------------------------------------------------------------------------- #

XIA_APA_XLSX = _DATA_ROOT / "miRNA" / "XIA_APA.xlsx"
APAATLAS_DIR = _DATA_ROOT / "external_cache" / "apaatlas"
APAATLAS_PDUI_ZIP = APAATLAS_DIR / "PDUI.txt.zip"
APAATLAS_BREAST_CACHE = APAATLAS_DIR / "breast_pdui_per_gene.tsv"
GTEX_SAMPLE_ATTR = _DATA_ROOT / "GTEx" / "GTEx_v10_SampleAttributesDS.txt"

POLYASITE_DIR = _DATA_ROOT / "external_cache" / "polyasite"
POLYASITE_BED = POLYASITE_DIR / "atlas.clusters.2.0.GRCh38.96.bed.gz"
POLYASITE_URL = ("https://polyasite.unibas.ch/download/atlas/2.0/"
                 "GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz")
GENCODE_SLIM = _DATA_ROOT / "gencode.v49.slim.parquet"


def target_apa_xia_brca(genes: Sequence[str]) -> pd.DataFrame:
    """Per target gene: **breast-TUMOUR** APA from Xia 2014 (TCGA tumour vs matched
    normal). `brca_dpdui` = ΔPDUI (tumour − normal); **negative = 3'UTR shortened in
    breast tumour** (distal poly-A usage drops → distal miRNA sites excised in tumour).

    This is the headline APA signal — breast + tumour + directional + significance-
    filtered (Xia lists only the recurrently-dynamic-APA events). Genes absent from
    the list are not recurrent-APA in breast tumour (→ neutral / not shortened).
    Still **gene-level** ⇒ cancels in within-gene arm re-ranking ⇒ diagnostic flag,
    not a weight (the arm-discriminating version needs site-vs-poly-A coordinates).
    """
    if not XIA_APA_XLSX.exists():
        print("[prior] XIA_APA.xlsx absent → breast-tumour APA diagnostic skipped")
        return pd.DataFrame(columns=["gene", "brca_dpdui", "brca_tumor_shortened"])
    df = pd.read_excel(XIA_APA_XLSX, sheet_name="TCGA_Tumor_Matched_Normal_APA")
    df["gene"] = df["Event_id"].astype(str).str.split("|").str[1]
    g = (df[["gene", "BRCA"]].dropna(subset=["BRCA"])
         .sort_values("BRCA").drop_duplicates("gene")
         .rename(columns={"BRCA": "brca_dpdui"}))
    g["brca_tumor_shortened"] = g["brca_dpdui"] < 0
    sub = g[g.gene.isin(set(genes))]
    print(f"[prior] Xia breast-TUMOUR APA: {len(sub)} HE target genes with significant "
          f"breast-tumour APA ({int(sub['brca_tumor_shortened'].sum())} shortened)")
    return g


def target_apa_breast(genes: Sequence[str]) -> pd.DataFrame:
    """Per target gene: median 3'UTR PDUI in **GTEx breast tissue** (APAatlas).

    PDUI (Percentage of Distal poly-A site Usage Index) ∈ [0,1]: 1 = long 3'UTR
    (distal poly-A used), low = shortened 3'UTR (proximal poly-A) — distal miRNA
    sites less available. Computed over the 218 GTEx breast samples present in the
    APAatlas v7 PDUI matrix; cached. **Breast-specific, but GTEx NORMAL breast**
    (tumour-acquired extra shortening not captured — no GET-able BRCA tumour PDUI).
    Gene-level ⇒ cancels in within-gene arm re-ranking ⇒ diagnostic flag, not a
    weight (same caveat as the PolyASite count, but breast-specific and continuous).
    """
    if APAATLAS_BREAST_CACHE.exists():
        g = pd.read_csv(APAATLAS_BREAST_CACHE, sep="\t")
    elif APAATLAS_PDUI_ZIP.exists() and GTEX_SAMPLE_ATTR.exists():
        import zipfile
        attr = pd.read_csv(GTEX_SAMPLE_ATTR, sep="\t", usecols=["SAMPID", "SMTSD"], dtype=str)
        breast = set(attr.loc[attr["SMTSD"] == "Breast - Mammary Tissue", "SAMPID"])
        z = zipfile.ZipFile(APAATLAS_PDUI_ZIP)
        with z.open("PDUI.txt") as fh:
            hdr = fh.readline().decode().rstrip("\n").split("\t")
        bcols = [c for c in hdr if c in breast]
        with z.open("PDUI.txt") as fh:
            df = pd.read_csv(fh, sep="\t", usecols=["gene_name"] + bcols)
        vals = df[bcols].apply(pd.to_numeric, errors="coerce")
        df["breast_pdui_median"] = vals.median(axis=1)
        df["breast_n_obs"] = vals.notna().sum(axis=1)
        g = (df.groupby("gene_name")
             .agg(breast_pdui_median=("breast_pdui_median", "median"),
                  breast_n_obs=("breast_n_obs", "max")).reset_index())
        g = g[g.breast_n_obs >= 20]
        APAATLAS_DIR.mkdir(parents=True, exist_ok=True)
        g.to_csv(APAATLAS_BREAST_CACHE, sep="\t", index=False)
    else:
        print("[prior] APAatlas breast PDUI unavailable → breast-APA diagnostic skipped")
        return pd.DataFrame(columns=["gene", "breast_pdui_median", "breast_apa_shortened"])

    cut = g["breast_pdui_median"].quantile(0.25)
    g["breast_apa_shortened"] = g["breast_pdui_median"] < cut
    g = g.rename(columns={"gene_name": "gene"})[["gene", "breast_pdui_median", "breast_apa_shortened"]]
    sub = g[g.gene.isin(set(genes))]
    print(f"[prior] APAatlas breast PDUI: {len(sub)} target genes covered; "
          f"{int(sub['breast_apa_shortened'].sum())} breast-shortened (PDUI<{cut:.2f})")
    return g


def target_apa_potential(genes: Sequence[str], do_fetch: bool = True) -> pd.DataFrame:
    """Per target gene: count of terminal-exon (3'UTR) poly-A clusters from the
    PolyASite 2.0 atlas — aggregate, cross-tissue APA mass-statistics.

    n_polyA_TE ≥ 2 ⇒ tandem 3'UTR ⇒ APA-shortenable target. **Gene-level**, so it
    does NOT discriminate a gene's competing arms (cancels in within-gene
    re-ranking) — supplied as a diagnostic flag, not a re-ranking term. The
    arm-discriminating (site-resolved) version needs TargetScan-site→genomic
    coordinate conversion and is left as future work. Not breast-specific.
    """
    POLYASITE_DIR.mkdir(parents=True, exist_ok=True)
    if not POLYASITE_BED.exists():
        if not do_fetch:
            print("[prior] PolyASite atlas absent and --no fetch → APA diagnostic skipped")
            return pd.DataFrame(columns=["gene", "n_polyA_TE", "apa_usage_sum", "apa_shortenable"])
        import urllib.request
        print("[prior] downloading PolyASite 2.0 atlas …")
        urllib.request.urlretrieve(POLYASITE_URL, POLYASITE_BED)

    # gene spans (protein-coding) for the requested targets
    gc = pd.read_parquet(GENCODE_SLIM, columns=["seqname", "feature", "start", "end",
                                                "strand", "gene_name", "gene_type"])
    gc = gc[(gc.feature == "gene") & (gc.gene_name.isin(set(genes)))].copy()
    gc["seqname"] = gc["seqname"].astype(str)
    gc["start"] = gc["start"].astype(int)
    gc["end"] = gc["end"].astype(int)

    # PolyASite terminal-exon clusters (col10 == 'TE'); col5 = usage TPM
    cols = ["chrom", "start", "end", "cid", "usage", "strand", "f7", "n_prot", "tpm", "ann", "sig"]
    pa = pd.read_csv(POLYASITE_BED, sep="\t", header=None, names=cols,
                     dtype={"chrom": str}, compression="gzip")
    pa = pa[pa["ann"] == "TE"].copy()
    pa["chrom"] = "chr" + pa["chrom"].str.replace("^chr", "", regex=True)
    pa["pos"] = (pa["start"] + pa["end"]) // 2

    # count TE clusters per gene (same chrom + strand, midpoint within span)
    import bisect
    by_cs: Dict[tuple, list] = {}
    for ch, st, pos, use in zip(pa["chrom"], pa["strand"], pa["pos"], pa["usage"]):
        by_cs.setdefault((ch, st), []).append((pos, use))
    for k in by_cs:
        by_cs[k].sort()
    rows = []
    for _, g in gc.iterrows():
        arr = by_cs.get((g.seqname, g.strand))
        n, usum = 0, 0.0
        if arr:
            positions = [p for p, _ in arr]
            lo = bisect.bisect_left(positions, g.start)
            hi = bisect.bisect_right(positions, g.end)
            n = hi - lo
            usum = float(sum(u for _, u in arr[lo:hi]))
        rows.append({"gene": g.gene_name, "n_polyA_TE": n,
                     "apa_usage_sum": round(usum, 3), "apa_shortenable": n >= 2})
    df = pd.DataFrame(rows).drop_duplicates("gene")
    print(f"[prior] APA (PolyASite): {int(df['apa_shortenable'].sum())}/{len(df)} "
          f"target genes have ≥2 terminal-exon poly-A clusters (shortenable 3'UTR)")
    return df


# --------------------------------------------------------------------------- #
# Assemble all axes
# --------------------------------------------------------------------------- #

def build_refined_prior(
    *,
    w_ts: float = 1.0,
    w_breast: float = 1.0,
    lambda_neg: float = 1.0,
    gamma_passenger: float = GAMMA_PASSENGER,
    do_fetch: bool = False,
) -> Tuple[pd.DataFrame, dict]:
    """Returns (edge_table_with_all_axes, coverage_stats)."""
    edge_pmids, he = collect_he_edge_pmids()
    he = he[["miRNA", "gene", "evidence_score"]].drop_duplicates().copy()
    genes = sorted(he["gene"].unique())

    # ---- base axis ---- #
    he["base"] = np.log1p(he["evidence_score"].clip(lower=0))

    # ---- axis A: TargetScan sequence site strength ---- #
    print("[prior] loading TargetScan weighted context++ (axis A) …")
    ts = _load_targetscan_weights(genes)          # (miRNA, gene, ts_weight)
    he = he.merge(ts, on=["miRNA", "gene"], how="left")
    he["ts_weight"] = he["ts_weight"].fillna(0.0)
    he["has_ts_site"] = he["ts_weight"] > 0
    he["ts_term"] = np.log1p(he["ts_weight"])

    # ---- breast literature axis ---- #
    print("[prior] scoring breast-context PMIDs (literature axis) …")
    all_pmids = set().union(*edge_pmids.values()) if edge_pmids else set()
    cache = fetch_pubmed_context(all_pmids, do_fetch=do_fetch)
    bscore = edge_breast_scores(edge_pmids, cache)[["miRNA", "gene", "n_breast_pmids"]]
    he = he.merge(bscore, on=["miRNA", "gene"], how="left")
    he["n_breast_pmids"] = he["n_breast_pmids"].fillna(0).astype(int)
    he["has_breast_pmid"] = he["n_breast_pmids"] > 0
    he["breast_term"] = np.log1p(he["n_breast_pmids"])

    # ---- axis B: negative evidence penalty ---- #
    print("[prior] counting Non-Functional MTI contradictions (axis B) …")
    neg = negative_evidence_counts(include_weak=False)
    he = he.merge(neg, on=["miRNA", "gene"], how="left")
    he["n_func"] = he["n_func"].fillna(0).astype(int)
    he["n_nonfunc"] = he["n_nonfunc"].fillna(0).astype(int)
    he["has_negative_evidence"] = he["n_nonfunc"] > 0
    denom = he["n_func"] + lambda_neg * he["n_nonfunc"]
    he["reliability_factor"] = np.where(denom > 0, he["n_func"] / denom, 1.0)
    # if an edge somehow has 0 func PMIDs (shouldn't for HE) keep neutral
    he.loc[he["n_func"] == 0, "reliability_factor"] = 1.0
    he["reliability_term"] = np.log(he["reliability_factor"].clip(lower=1e-6))

    # ---- axis C: guide/passenger arm (miRGeneDB, sequence-matched) ---- #
    print("[prior] classifying guide/passenger arms (axis C, miRGeneDB) …")
    arm_status = mirgenedb_arm_status(do_fetch=do_fetch)
    arm_map = dict(zip(arm_status["miRNA"], arm_status["arm_class"])) if not arm_status.empty else {}
    he["arm_class"] = he["miRNA"].map(arm_map).fillna("unknown")
    arm_src = dict(zip(arm_status["miRNA"], arm_status["arm_class_source"])) if not arm_status.empty else {}
    he["arm_class_source"] = he["miRNA"].map(arm_src).fillna("none")
    # only confirmed passenger is penalised; guide/ambiguous/unknown stay neutral
    he["arm_factor"] = np.where(he["arm_class"].eq("passenger"), gamma_passenger, 1.0)
    he["arm_term"] = np.log(he["arm_factor"])
    # validated-universe flag: arm recognised by miRGeneDB (high-confidence miRNA)
    he["mirgenedb_validated"] = ~he["arm_class"].eq("unknown")

    # ---- axis D (diagnostic): target 3'UTR APA potential (PolyASite mass-stats) ---- #
    print("[prior] target 3'UTR APA potential (axis D, PolyASite — diagnostic only) …")
    apa = target_apa_potential(genes, do_fetch=do_fetch)
    he = he.merge(apa, on="gene", how="left")
    he["n_polyA_TE"] = he["n_polyA_TE"].fillna(0).astype(int)
    he["apa_shortenable"] = he["apa_shortenable"].fillna(False)
    # breast-TUMOUR APA (Xia 2014 ΔPDUI) — the HEADLINE APA diagnostic (directional)
    apa_x = target_apa_xia_brca(genes)
    if not apa_x.empty:
        he = he.merge(apa_x, on="gene", how="left")
        he["brca_tumor_shortened"] = he["brca_tumor_shortened"].fillna(False).astype(bool)
    else:
        he["brca_dpdui"] = np.nan
        he["brca_tumor_shortened"] = False
    # breast-NORMAL APA (APAatlas GTEx breast PDUI) — secondary
    apa_b = target_apa_breast(genes)
    if not apa_b.empty:
        he = he.merge(apa_b, on="gene", how="left")
        he["breast_apa_shortened"] = he["breast_apa_shortened"].fillna(False).astype(bool)
    else:
        he["breast_pdui_median"] = np.nan
        he["breast_apa_shortened"] = False
    # NOTE: all APA signals are gene-level → deliberately NOT in refined_weight (cancels within-gene)

    # ---- combined refined weight ---- #
    he["refined_weight"] = (
        he["base"]
        + w_ts * he["ts_term"]
        + w_breast * he["breast_term"]
        + he["reliability_term"]
        + he["arm_term"]
    )

    # ---- per-edge multi-axis confidence (4 breadth signals) ---- #
    he["n_supporting_axes"] = (
        he["has_ts_site"].astype(int)
        + he["has_breast_pmid"].astype(int)
        + (~he["has_negative_evidence"]).astype(int)
        + he["arm_class"].eq("guide").astype(int)
    )

    coverage = {
        "n_he_edges": int(len(he)),
        "axis_A_ts_site": {
            "n": int(he["has_ts_site"].sum()),
            "frac": float(he["has_ts_site"].mean()),
        },
        "axis_breast_pmid": {
            "n": int(he["has_breast_pmid"].sum()),
            "frac": float(he["has_breast_pmid"].mean()),
        },
        "axis_B_negative_evidence": {
            "n": int(he["has_negative_evidence"].sum()),
            "frac": float(he["has_negative_evidence"].mean()),
        },
        "axis_C_arm": {
            "passenger": int((he["arm_class"] == "passenger").sum()),
            "guide": int((he["arm_class"] == "guide").sum()),
            "ambiguous": int((he["arm_class"] == "ambiguous").sum()),
            "unknown_neutral": int((he["arm_class"] == "unknown").sum()),
            "frac_classified": float(he["arm_class"].isin(["guide", "passenger", "ambiguous"]).mean()),
            "frac_passenger": float((he["arm_class"] == "passenger").mean()),
        },
        "axis_D_apa_diagnostic": {
            "xia_brca_tumor_edges_with_apa": int(he["brca_dpdui"].notna().sum()),
            "xia_brca_tumor_edges_shortened": int(he["brca_tumor_shortened"].sum()),
            "xia_brca_tumor_frac_shortened": float(he["brca_tumor_shortened"].mean()),
            "apaatlas_breast_normal_edges_shortened": int(he["breast_apa_shortened"].sum()),
            "polyasite_pan_edges_shortenable": int(he["apa_shortenable"].sum()),
            "note": "3 gene-level diagnostics — Xia breast-TUMOUR ΔPDUI (headline, directional), "
                    "APAatlas GTEx-breast-NORMAL PDUI, PolyASite pan-tissue; all cancel in "
                    "within-gene re-ranking → flags, not weight terms",
        },
        "edges_no_extra_axis": int((he["n_supporting_axes"] == 1).sum()),
    }
    print(f"[prior] axis coverage of {len(he)} HE edges: "
          f"TS site {coverage['axis_A_ts_site']['frac']*100:.0f}%  ·  "
          f"breast PMID {coverage['axis_breast_pmid']['frac']*100:.0f}%  ·  "
          f"contradicted {coverage['axis_B_negative_evidence']['frac']*100:.1f}%  ·  "
          f"arm classified {coverage['axis_C_arm']['frac_classified']*100:.0f}% "
          f"(passenger {coverage['axis_C_arm']['frac_passenger']*100:.1f}%)  ·  "
          f"breast-tumour-APA-shortened target {coverage['axis_D_apa_diagnostic']['xia_brca_tumor_frac_shortened']*100:.0f}%")
    return he, coverage


# --------------------------------------------------------------------------- #
# miRGeneDB validated-universe variant + exoneration of discarded arms
# --------------------------------------------------------------------------- #

def mirgenedb_exoneration(he: pd.DataFrame, corr_path: Path = CORR_TABLE) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Do the arms miRGeneDB *discards* still couple with their targets?

    miRGeneDB excludes most of miRBase as dubious (artefacts / non-miRNA fragments).
    'Validated' = arm recognised by miRGeneDB (guide/passenger/ambiguous);
    'discarded' = arm absent (arm_class == unknown). We compare TCGA-BRCA target
    coupling (partial ρ, CPE+HRD+prolif-adjusted, from edge_partial_corr_panels)
    between the two groups. This is a **characterisation of the miRGeneDB filter**,
    NOT a prior term — coupling is TCGA-derived, so it must not feed back into the
    weight (it would be circular with the coupling analyses). The discarded arms
    that nonetheless couple (neg-sig) are **exoneration candidates**: possibly
    functional despite miRGeneDB exclusion — flagged for follow-up, not auto-trusted.
    """
    if not corr_path.exists():
        print(f"[prior] exoneration skipped — no partial-corr table at {corr_path}")
        return pd.DataFrame(), pd.DataFrame()
    corr = pd.read_csv(corr_path, sep="\t")
    m = (he[["miRNA", "gene", "arm_class"]]
         .merge(corr[["miRNA", "gene", "rho_adj", "q_adj"]], on=["miRNA", "gene"], how="inner")
         .dropna(subset=["rho_adj"]))
    m["mirgenedb"] = np.where(m["arm_class"].eq("unknown"),
                              "discarded", "validated")
    rows = []
    for grp, g in m.groupby("mirgenedb"):
        rows.append({
            "group": grp, "n_edges": len(g),
            "median_rho_adj": round(float(g["rho_adj"].median()), 4),
            "frac_negative": round(float((g["rho_adj"] < 0).mean()), 3),
            "frac_neg_sig_q05": round(float(((g["rho_adj"] < 0) & (g["q_adj"] < 0.05)).mean()), 3),
        })
    summary = pd.DataFrame(rows)
    # exoneration candidates: discarded arm but neg-sig coupling in BRCA
    cand = m[(m["mirgenedb"] == "discarded") & (m["rho_adj"] < 0) & (m["q_adj"] < 0.05)].copy()
    cand = cand.sort_values("rho_adj")[["miRNA", "gene", "arm_class", "rho_adj", "q_adj"]]
    print(f"[prior] exoneration — neg-sig coupling: "
          f"validated {summary.loc[summary.group=='validated','frac_neg_sig_q05'].iloc[0] if (summary.group=='validated').any() else float('nan'):.3f} "
          f"vs discarded {summary.loc[summary.group=='discarded','frac_neg_sig_q05'].iloc[0] if (summary.group=='discarded').any() else float('nan'):.3f}"
          f"  ·  {len(cand)} exoneration candidates")
    return summary, cand


# --------------------------------------------------------------------------- #
# Per-target re-ranking + axis attribution
# --------------------------------------------------------------------------- #

def _dominant_axis(pan_top: pd.Series, ref_top: pd.Series, w_ts: float, w_breast: float) -> str:
    """Which axis contributed most to ref_top overtaking pan_top."""
    d_ts = w_ts * (ref_top["ts_term"] - pan_top["ts_term"])
    d_brst = w_breast * (ref_top["breast_term"] - pan_top["breast_term"])
    # negative evidence / passenger help ref_top if they penalize pan_top more
    d_neg = ref_top["reliability_term"] - pan_top["reliability_term"]
    d_arm = ref_top["arm_term"] - pan_top["arm_term"]
    contribs = {
        "sequence (TargetScan)": d_ts,
        "breast literature": d_brst,
        "negative evidence": d_neg,
        "guide/passenger (miRGeneDB)": d_arm,
    }
    best = max(contribs, key=contribs.get)
    return best if contribs[best] > 1e-9 else "base only (tie)"


def per_target_reranking(
    he: pd.DataFrame, *, w_ts: float, w_breast: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = he.copy()
    df["pan_rank"] = df.groupby("gene")["base"].rank(ascending=False, method="min")
    df["refined_rank"] = df.groupby("gene")["refined_weight"].rank(ascending=False, method="min")

    flips = []
    for gene, g in df.groupby("gene"):
        if len(g) < 2:
            continue
        pan_top = g.sort_values("base", ascending=False).iloc[0]
        ref_top = g.sort_values("refined_weight", ascending=False).iloc[0]
        flipped = pan_top["miRNA"] != ref_top["miRNA"]
        rec = {
            "gene": gene,
            "n_regulators": len(g),
            "pan_top_arm": pan_top["miRNA"],
            "pan_top_evidence": pan_top["evidence_score"],
            "refined_top_arm": ref_top["miRNA"],
            "refined_top_evidence": ref_top["evidence_score"],
            "refined_top_ts_weight": round(float(ref_top["ts_weight"]), 3),
            "refined_top_n_breast_pmids": int(ref_top["n_breast_pmids"]),
            "refined_top_n_nonfunc": int(ref_top["n_nonfunc"]),
            "refined_top_arm_class": ref_top["arm_class"],
            "pan_top_arm_class": pan_top["arm_class"],
            "refined_top_n_supporting_axes": int(ref_top["n_supporting_axes"]),
            "reassigned": bool(flipped),
            "dominant_axis": _dominant_axis(pan_top, ref_top, w_ts, w_breast) if flipped else "",
        }
        flips.append(rec)
    flip_df = pd.DataFrame(flips)
    if not flip_df.empty:
        flip_df = flip_df.sort_values(
            ["reassigned", "refined_top_n_supporting_axes"], ascending=[False, False]
        ).reset_index(drop=True)

    n_flip = int(flip_df["reassigned"].sum()) if not flip_df.empty else 0
    print(f"[prior] multi-regulator genes: {len(flip_df)}  "
          f"top-arm reassigned: {n_flip}")
    if n_flip:
        by_axis = flip_df.loc[flip_df["reassigned"], "dominant_axis"].value_counts()
        for ax, n in by_axis.items():
            print(f"    {ax}: {n}")
    return df, flip_df


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #

def plot_summary(
    he: pd.DataFrame, flip_df: pd.DataFrame, coverage: dict, out_path: Path,
    *, dpi: int = 150, n_show: int = 25,
) -> None:
    import matplotlib.pyplot as plt

    flips = flip_df[flip_df["reassigned"]] if not flip_df.empty else flip_df
    top = flips.head(n_show).iloc[::-1] if not flips.empty else flips
    n_rows = max(len(top), 1)
    fig, axes = plt.subplots(
        1, 2, figsize=(16, max(6.5, 0.34 * n_rows + 2.5)),
        gridspec_kw={"width_ratios": [1, 1.5]})

    # left: axis coverage + flip attribution
    ax = axes[0]
    cov_labels = ["TS site\n(seq)", "breast\nPMID", "contradicted\n(neg ev)", "passenger\narm"]
    cov_vals = [coverage["axis_A_ts_site"]["frac"] * 100,
                coverage["axis_breast_pmid"]["frac"] * 100,
                coverage["axis_B_negative_evidence"]["frac"] * 100,
                coverage["axis_C_arm"]["frac_passenger"] * 100]
    bars = ax.bar(cov_labels, cov_vals,
                  color=["#4E79A7", "#C77CA6", "#E15759", "#59A14F"], alpha=0.85)
    for b, v in zip(bars, cov_vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 1, f"{v:.0f}%",
                ha="center", fontsize=9)
    ax.set_ylabel("% of HE edges touched", fontsize=10)
    ax.set_title(
        f"Axis coverage over {coverage['n_he_edges']} HE edges\n"
        "(edges a given axis has any information for)",
        fontsize=10)
    ax.set_ylim(0, max(cov_vals) * 1.25 + 5)
    ax.grid(axis="y", alpha=0.2, linestyle="--")

    # right: top reassignments colored by dominant axis
    ax = axes[1]
    if top.empty:
        ax.text(0.5, 0.5, "No top-arm reassignments", ha="center", va="center")
        ax.axis("off")
    else:
        y = np.arange(len(top))
        colors = [AXIS_COLORS.get(a, "#999999") for a in top["dominant_axis"]]
        ax.barh(y, top["refined_top_n_supporting_axes"], color=colors, alpha=0.85)
        labels = [
            f"{r['gene']}: {r['pan_top_arm'].replace('hsa-','')}"
            f"(ev{r['pan_top_evidence']:.0f}) → {r['refined_top_arm'].replace('hsa-','')}"
            f" [ts{r['refined_top_ts_weight']:.1f}, brst{r['refined_top_n_breast_pmids']}]"
            for _, r in top.iterrows()
        ]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.5)
        ax.set_xlabel("# supporting axes for the refined-credited arm (max 4: TS, breast, no-neg, guide)", fontsize=8.5)
        n_flip = int(flips["reassigned"].sum())
        ax.set_title(
            f"Top-arm reassignments, colored by driving axis  ({n_flip} genes)\n"
            "pan-context arm → refined-credited arm",
            fontsize=10)
        ax.set_xlim(0, 4.3)
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        from matplotlib.patches import Patch
        present = [a for a in AXIS_COLORS if a in set(top["dominant_axis"])]
        ax.legend(handles=[Patch(facecolor=AXIS_COLORS[a], label=a) for a in present],
                  fontsize=7, frameon=False, loc="lower right")
        ax.margins(y=0.01)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[prior] figure → {out_path}")


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(
    *,
    out_dir: Path = OUT_DIR,
    w_ts: float = 1.0,
    w_breast: float = 1.0,
    lambda_neg: float = 1.0,
    gamma_passenger: float = GAMMA_PASSENGER,
    do_fetch: bool = False,
    dpi: int = 150,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    he, coverage = build_refined_prior(
        w_ts=w_ts, w_breast=w_breast, lambda_neg=lambda_neg,
        gamma_passenger=gamma_passenger, do_fetch=do_fetch)
    he.to_csv(out_dir / "edge_prior_refined.tsv", sep="\t", index=False)

    edge_tbl, flip_df = per_target_reranking(he, w_ts=w_ts, w_breast=w_breast)
    flip_df.to_csv(out_dir / "per_target_arm_reassignments_multiaxis.tsv",
                   sep="\t", index=False)

    # miRGeneDB-validated universe variant (high-confidence miRNA edges only)
    validated = he[he["mirgenedb_validated"]]
    validated[["miRNA", "gene", "evidence_score", "arm_class",
               "refined_weight"]].to_csv(
        out_dir / "edges_mirgenedb_validated_universe.tsv", sep="\t", index=False)
    n_val_edges = int(he["mirgenedb_validated"].sum())
    n_val_genes = int(validated["gene"].nunique())
    print(f"[prior] miRGeneDB-validated universe: {n_val_edges}/{len(he)} edges "
          f"({100*he['mirgenedb_validated'].mean():.0f}%), {n_val_genes} target genes")

    # exoneration of discarded arms (characterisation, not a prior term)
    exo_summary, exo_cand = mirgenedb_exoneration(he)
    if not exo_summary.empty:
        exo_summary.to_csv(out_dir / "mirgenedb_exoneration_summary.tsv", sep="\t", index=False)
        exo_cand.to_csv(out_dir / "mirgenedb_exoneration_candidates.tsv", sep="\t", index=False)

    plot_summary(he, flip_df, coverage,
                 out_dir / "figures" / "edge_prior_refinement.png", dpi=dpi)

    by_axis = (flip_df.loc[flip_df["reassigned"], "dominant_axis"].value_counts().to_dict()
               if not flip_df.empty else {})
    manifest = {
        "module": "mirna_hallmark.analyses.edge_panels.edge_prior_refinement",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "weights": {"w_ts": w_ts, "w_breast": w_breast, "lambda_neg": lambda_neg,
                    "gamma_passenger": gamma_passenger},
        "formula": "refined = log1p(evidence) + w_ts·log1p(ts_weight) "
                   "+ w_breast·log1p(n_breast_pmids) + log(n_func/(n_func+λ·n_nonfunc)) "
                   "+ log(arm_factor)",
        "axes_included": {
            "base": "literature study count (spine anchor)",
            "A_sequence": "TargetScan Σ|weighted context++| (folds site-type + conservation)",
            "breast_literature": "breast-context Functional-MTI PMIDs (edge_breast_context)",
            "B_negative_evidence": "Non-Functional MTI contradictions (strong only)",
            "C_guide_passenger": "miRGeneDB guide vs passenger (mat=1 / star=1), mapped to "
                                 "miRBase by exact-seq + seed/arm-end + hairpin propagation; "
                                 "only confirmed passenger arms down-weighted (γ), rest neutral",
            "D_apa_diagnostic": "PolyASite 2.0 terminal-exon poly-A cluster count per target "
                                "(n_polyA_TE≥2 = shortenable); GENE-LEVEL → cancels in arm "
                                "re-ranking, supplied as diagnostic flag not weight; cross-tissue",
        },
        "axes_deferred": {
            "D_site_resolved_breast_APA": "arm-discriminating version (is this site distal to a "
                                          "proximal poly-A cluster) needs TargetScan-site→genomic "
                                          "conversion; breast-specific ΔPDUI needs a non-GET source "
                                          "(TC3A JS-only SPA, APAatlas query-UI)",
        },
        "axes_deliberately_excluded": {
            "tcga_expression": "outcome variable; a prior must not see it (circularity)",
        },
        "coverage": coverage,
        "n_multi_regulator_genes": int(len(flip_df)),
        "n_top_arm_reassigned": int(flip_df["reassigned"].sum()) if not flip_df.empty else 0,
        "reassignments_by_dominant_axis": by_axis,
        "status": "per-target NOMINATION/AUDIT overlay; NOT the landscape budget",
        "validation": "none — no external breast-miRNA cohort (METABRIC no miRNA; CPTAC=TCGA-reprocessed)",
        "mirgenedb_validated_universe": {
            "n_edges": int(he["mirgenedb_validated"].sum()),
            "frac_edges": float(he["mirgenedb_validated"].mean()),
            "note": "high-confidence miRNA subset (arm recognised by miRGeneDB); "
                    "optional restricted universe — filters dubious miRBase arms",
        },
        "exoneration": {
            "summary": (exo_summary.to_dict("records") if not exo_summary.empty else []),
            "n_candidates": int(len(exo_cand)),
            "note": "discarded (non-miRGeneDB) arms with neg-sig BRCA coupling; "
                    "characterisation of the filter, NOT folded into the prior (circularity)",
        },
        "outputs": {
            "edge_table": "edge_prior_refined.tsv",
            "reassignments": "per_target_arm_reassignments_multiaxis.tsv",
            "validated_universe": "edges_mirgenedb_validated_universe.tsv",
            "exoneration_summary": "mirgenedb_exoneration_summary.tsv",
            "exoneration_candidates": "mirgenedb_exoneration_candidates.tsv",
            "figure": "figures/edge_prior_refinement.png",
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[prior] done — {manifest['n_top_arm_reassigned']} reassignments, "
          f"by axis: {by_axis}")
    return {"edge_tbl": edge_tbl, "flip_df": flip_df, "coverage": coverage,
            "manifest": manifest}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--w-ts", type=float, default=1.0,
                    help="weight on the TargetScan sequence axis")
    ap.add_argument("--w-breast", type=float, default=1.0,
                    help="weight on the breast-literature axis")
    ap.add_argument("--lambda-neg", type=float, default=1.0,
                    help="penalty strength per contradicting Non-Functional MTI")
    ap.add_argument("--gamma-passenger", type=float, default=GAMMA_PASSENGER,
                    help="down-weight factor for a confirmed miRGeneDB passenger arm")
    ap.add_argument("--fetch", dest="do_fetch", action="store_true",
                    help="(re)fetch missing PMIDs + miRGeneDB; default uses the cache")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, w_ts=args.w_ts, w_breast=args.w_breast,
        lambda_neg=args.lambda_neg, gamma_passenger=args.gamma_passenger,
        do_fetch=args.do_fetch, dpi=args.dpi)


if __name__ == "__main__":
    main()
