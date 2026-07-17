"""End-to-end CN -> miRNA regulatory-network attribution.

The subproject measures two *halves* of the dosage chain separately:

  Link 1  CN -> miRNA expression      (``mirna_locus_cnv`` concordance, ~196/756 arms)
  Link 2  miRNA pressure -> target    (``target_combined_anticorr``, ~373/1349 genes)

Neither answers the question this module exists for: **how much of the bulk
repressive network is actually driven by the genetic (copy-number) layer, in
which context, and at which genomic granularity** (single locus / chromosome-arm
aneuploidy / dense polycistronic cluster / combined multi-locus load).

Two analyses, both reusing the canonical spine pressure
(``softmax_z_logrpm`` + ``evidence_mass`` on M0 miRTarBase edges, AGO-gated) and the
same partial-Spearman machinery (CPE + HRD + target CN) as ``target_combined_anticorr``.
Every CN-relevant coupling is *also* reported after adding a **broad genome-wide
aneuploidy** covariate (Taylor et al. arm-level aneuploidy score; ``*_aneu`` columns)
and, cumulatively, a **proliferation metagene** (E2F+G2M, as ``robustness_checks``;
``*_full`` columns), so a locus-specific miRNA-CN effect on the network is separated
from the genome-wide CIN tide and the proliferation/metabolic state that CPE+HRD only
partly remove:

A. **CN-driven vs CN-residual pressure (the "how much").** Each regulator arm's
   log2(RPM+1) is split per-arm into a CN-explained component ``cn_hat`` (OLS fit on
   the arm's own copy number) and a CN-independent residual. Three pressure tracks are
   built on the *same edges* — ``total`` (real expression), ``cn`` (cn_hat only), and
   ``resid`` (residual only, re-anchored to the cohort-median abundance) — and each is
   tested against target RNA. The share of the total repressive coupling that the
   ``cn`` track reproduces = the fraction of the network attributable to copy number.

B. **Combined multi-locus CN-load granularity contest.** For each target gene a single
   per-sample scalar — the evidence-weighted summed CN dosage of *all* its regulators —
   is built at four resolutions (``locus`` arm CN; ``chrom`` chromosome-mean aneuploidy;
   ``cluster`` dense-cluster co-dosage; ``single_best`` the most-amplified single
   regulator) and tested for repressor-consistent (negative) coupling with target RNA.
   Answers whether combining loci beats any single locus and which scale carries it.

Outputs under ``output/cn_dosage_attribution/``:
  - ``cn_attribution_cohort.tsv``         per-gene total/cn/resid partial-rho + attribution
  - ``cn_attribution_by_context.tsv``     same, within PAM50 / immune / TNBC strata
  - ``cn_attribution_summary.tsv``        headline counts + median CN-attribution per context
  - ``cn_load_granularity_cohort.tsv``    per-gene combined-load partial-rho per granularity
  - ``cn_load_granularity_summary.tsv``   neg-FDR counts per granularity (cohort + context)
  - ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges, method_blurb, pressure_kwargs
from mirna_hallmark.pressure_engine import (
    compute_gene_pressure as engine_pressure,
    filter_edges_by_abundance,
)

DIPLOID = 2.0
MIN_FIT = 25          # min paired (CN, expr) samples to fit a per-arm CN slope
MIN_N = 25            # min samples for a correlation test
# Contexts to stratify by, beyond the cohort. PAM50 is the rigorous primary
# (stratum-fixed, so excluded from its own partial); the others are reported with
# CPE+HRD+target-CN partials but are NOT PAM50-adjusted (composition caveat).
CONTEXT_COLS = ("PAM50_final", "thornsson_immune_subtype", "tnbc_subtype_4")
# Broad genome-wide aneuploidy / CIN tide (Taylor et al. arm-level aneuploidy
# score). Added as an *extra* covariate so a locus-specific miRNA-CN effect on the
# network can be separated from the genome-wide aneuploidy background that
# CPE+HRD only partly capture. Every CN-relevant coupling is reported both at
# baseline (CPE+HRD+target-CN) and aneuploidy-adjusted (``*_aneu``).
ANEUPLOIDY_COVAR = "thornsson_aneuploidy_score"
# Proliferation metagene (cell-cycle Hallmarks) — the third covariate layer. The
# ``*_full`` columns adjust for CPE+HRD+target-CN+aneuploidy+proliferation, the
# bulletproof control: amplified/aneuploid tumors are more proliferative, and
# metabolic/cell-cycle target programs co-vary with that state, not necessarily with
# miRNA repression. Same metagene as ``robustness_checks`` Aim-1.
PROLIFERATION_SETS = ("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")

_CACHE = C.MIRNA_LOCUS_CNV_DIR / "tables" / "sample_entity_cnv.tsv.gz"
_GENOME_MAP = C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_locus_genome_map.tsv"
_CLUSTER_MAP = C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_locus_cluster_summary.tsv"


# --------------------------------------------------------------------------- #
# CN panels: arm x participant copy number at several granularities
# --------------------------------------------------------------------------- #
def _arm_to_mimat(cache: pd.DataFrame) -> pd.Series:
    """arm label (hsa-*) -> MIMAT accession from the arm-level cache rows."""
    arm = cache.loc[cache["entity_level"] == "arm", ["entity_label", "mature_accession"]]
    arm = arm.dropna(subset=["mature_accession"]).drop_duplicates("entity_label")
    return arm.set_index("entity_label")["mature_accession"]


def _mimat_locus_map() -> pd.DataFrame:
    """Long MIMAT -> (locus_id, chrom, midpoint) from the locus genome map."""
    gm = pd.read_csv(_GENOME_MAP, sep="\t")
    gm = gm.dropna(subset=["mature_accessions"])
    gm = gm.assign(mature_accession=gm["mature_accessions"].str.split(",")).explode(
        "mature_accession"
    )
    gm["mature_accession"] = gm["mature_accession"].str.strip()
    return gm[["mature_accession", "locus_id", "chrom", "midpoint"]].dropna(
        subset=["mature_accession"]
    )


def _assign_clusters(mimat_loci: pd.DataFrame) -> pd.Series:
    """locus_id -> cluster_id for loci falling in a dense (>=3 hairpins/100kb) cluster."""
    if not _CLUSTER_MAP.exists():
        return pd.Series(dtype=object)
    cl = pd.read_csv(_CLUSTER_MAP, sep="\t")
    loci = mimat_loci.drop_duplicates("locus_id")[["locus_id", "chrom", "midpoint"]]
    out: Dict[str, str] = {}
    for _, c in cl.iterrows():
        lo, hi = c["cluster_start_mb"] * 1e6, c["cluster_end_mb"] * 1e6
        cid = f"{c['chrom']}:{c['cluster_start_mb']:.2f}-{c['cluster_end_mb']:.2f}"
        hit = loci[(loci["chrom"] == c["chrom"]) & loci["midpoint"].between(lo, hi)]
        for lid in hit["locus_id"]:
            out[str(lid)] = cid
    return pd.Series(out, name="cluster_id")


def build_cn_panels(arm_labels: Sequence[str]) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """Return (panels, arm_meta).

    ``panels`` maps granularity -> arm x participant copy-number frame:
      - ``locus``   : the arm's own (paralog-weighted) composite CN.
      - ``chrom``   : per-participant mean locus CN on the arm's chromosome (broad aneuploidy).
      - ``cluster`` : per-participant mean locus CN of the arm's dense cluster
                      (falls back to ``locus`` for arms outside any cluster).
    ``arm_meta`` carries arm -> mimat / chrom / locus_id / cluster_id / in_cluster.
    """
    cache = pd.read_csv(_CACHE, sep="\t", low_memory=False)
    cache["copy_number"] = pd.to_numeric(cache["copy_number"], errors="coerce")

    # --- locus (arm composite) panel ---
    arm_rows = cache.loc[cache["entity_level"] == "arm"].drop_duplicates(
        ["entity_label", "participant"]
    )
    locus_panel = arm_rows.pivot(index="entity_label", columns="participant", values="copy_number")
    keep = [a for a in arm_labels if a in locus_panel.index]
    locus_panel = locus_panel.loc[keep]

    # --- arm -> chrom / locus_id / cluster meta ---
    a2m = _arm_to_mimat(cache)
    mimat_loci = _mimat_locus_map()
    # dominant locus per MIMAT = first (paralog tie broken by map order)
    mimat_first = mimat_loci.drop_duplicates("mature_accession").set_index("mature_accession")
    cluster_of_locus = _assign_clusters(mimat_loci)

    meta_rows = []
    for arm in keep:
        mimat = a2m.get(arm)
        chrom = locus_id = cluster_id = None
        if mimat is not None and mimat in mimat_first.index:
            chrom = mimat_first.loc[mimat, "chrom"]
            locus_id = mimat_first.loc[mimat, "locus_id"]
            cluster_id = cluster_of_locus.get(str(locus_id))
        meta_rows.append(
            {"arm": arm, "mimat": mimat, "chrom": chrom, "locus_id": locus_id,
             "cluster_id": cluster_id, "in_cluster": cluster_id is not None}
        )
    arm_meta = pd.DataFrame(meta_rows).set_index("arm")

    # --- per-participant locus-level CN with chrom + cluster keys ---
    loc_rows = cache.loc[cache["entity_level"] == "locus"].copy()
    loc_rows["chrom"] = loc_rows["overlap_segment"].astype(str).str.split(":").str[0]
    loc_rows = loc_rows.dropna(subset=["copy_number"])

    chrom_mean = (
        loc_rows.groupby(["participant", "chrom"])["copy_number"].mean().unstack("participant")
    )  # chrom x participant

    loc_rows["cluster_id"] = loc_rows["entity_id"].astype(str).map(cluster_of_locus)
    clustered = loc_rows.dropna(subset=["cluster_id"])
    cluster_mean = (
        clustered.groupby(["participant", "cluster_id"])["copy_number"].mean().unstack("participant")
        if not clustered.empty
        else pd.DataFrame()
    )  # cluster x participant

    participants = locus_panel.columns
    # --- chrom panel: broadcast chrom-mean onto each arm ---
    chrom_panel = pd.DataFrame(index=keep, columns=participants, dtype=float)
    for arm in keep:
        ch = arm_meta.loc[arm, "chrom"]
        if ch is not None and ch in chrom_mean.index:
            chrom_panel.loc[arm] = chrom_mean.loc[ch].reindex(participants)
    # arms with no chrom mapping fall back to their own composite CN
    chrom_panel = chrom_panel.where(chrom_panel.notna(), locus_panel)

    # --- cluster panel: clustered arms get cluster-mean, else their own composite ---
    cluster_panel = locus_panel.copy()
    for arm in keep:
        cid = arm_meta.loc[arm, "cluster_id"]
        if cid is not None and not cluster_mean.empty and cid in cluster_mean.index:
            cluster_panel.loc[arm] = cluster_mean.loc[cid].reindex(participants)

    panels = {"locus": locus_panel, "chrom": chrom_panel, "cluster": cluster_panel}
    return panels, arm_meta


# --------------------------------------------------------------------------- #
# Per-arm CN / residual decomposition of expression
# --------------------------------------------------------------------------- #
def decompose_expression(expr: pd.DataFrame, cn: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split each arm's log2(RPM+1) into CN-driven (``cn_hat``) and CN-residual.

    Per arm with >= ``MIN_FIT`` paired finite (CN, expr) samples and CN variance:
    OLS ``expr ~ CN`` gives ``cn_hat = a + b*CN``; residual = expr - cn_hat. The
    residual track is re-anchored to the arm's cohort-median abundance so the
    ``logrpm`` term of the spine pressure is comparable across tracks. Arms with no
    usable CN contribute *all* their variation to the residual track (cn_hat = const).
    """
    samples = expr.columns
    cn = cn.reindex(index=expr.index, columns=samples)
    cn_hat = pd.DataFrame(index=expr.index, columns=samples, dtype=float)
    resid = pd.DataFrame(index=expr.index, columns=samples, dtype=float)

    for arm in expr.index:
        x = pd.to_numeric(expr.loc[arm], errors="coerce")
        med = float(x.median()) if np.isfinite(x.median()) else 0.0
        mean_x = float(x.mean()) if np.isfinite(x.mean()) else 0.0
        c = pd.to_numeric(cn.loc[arm], errors="coerce") if arm in cn.index else pd.Series(np.nan, index=samples)
        both = x.notna() & c.notna()
        if int(both.sum()) >= MIN_FIT and float(c[both].std()) > 0:
            cc, xx = c[both], x[both]
            b = float(np.cov(cc, xx)[0, 1] / np.var(cc))
            a = float(xx.mean() - b * cc.mean())
            xhat = a + b * c            # NaN where CN missing
            xhat = xhat.where(c.notna(), mean_x)
        else:
            xhat = pd.Series(mean_x, index=samples)
        cn_hat.loc[arm] = xhat.where(x.notna(), np.nan)
        resid.loc[arm] = (x - xhat) + med     # NaN where x missing (x term)
    return cn_hat, resid


# --------------------------------------------------------------------------- #
# Shared correlation helper
# --------------------------------------------------------------------------- #
def _matrix_row(mat: pd.DataFrame, row_id: str) -> pd.Series:
    row = mat.loc[row_id]
    if isinstance(row, pd.DataFrame):
        row = row.apply(pd.to_numeric, errors="coerce").median(axis=0)
    return pd.to_numeric(row, errors="coerce")


def _partial(y: pd.Series, x: pd.Series, conf: Optional[pd.DataFrame],
             target_cn: Optional[pd.Series], aneu: Optional[pd.Series] = None,
             prolif: Optional[pd.Series] = None):
    """correlation_pair with CPE+HRD (+ target CN, + broad aneuploidy, + proliferation)."""
    cov = conf.copy() if (conf is not None and len(conf.columns)) else None
    extras: Dict[str, pd.Series] = {}
    if target_cn is not None:
        extras["target_cn"] = target_cn
    if aneu is not None:
        extras["aneuploidy"] = aneu
    if prolif is not None:
        extras["proliferation"] = prolif
    if extras:
        cov = cov.copy() if cov is not None else pd.DataFrame(index=y.index)
        for k, v in extras.items():
            cov[k] = pd.to_numeric(v, errors="coerce").reindex(y.index)
    if cov is not None:
        cov = D.augment_tcga_batch(cov)
    # cn_dosage only consumes the Spearman outputs (rho_/p_/sp_rho_); skip the
    # unused Pearson half. Batch dummies are cached in tcga_batch (same sample set
    # recurs across genes), so augment_tcga_batch is ~free after the first call.
    return S.correlation_pair(y, x, cov, spearman_only=True)


def _proliferation_metagene(rna: pd.DataFrame) -> pd.Series:
    """Per-participant mean z of HALLMARK_E2F_TARGETS + G2M_CHECKPOINT (cell-cycle).

    Same metagene as ``robustness_checks`` Aim-1: the proliferation proxy a reviewer
    demands, since amplified/aneuploid tumors are more proliferative and metabolic /
    cell-cycle target programs co-vary with that state rather than with miRNA repression.
    """
    hs = HallmarkSets.load()
    genes = sorted({g for s in PROLIFERATION_SETS for g in hs.sets.get(s, [])})
    present = [g for g in genes if g in rna.index]
    if not present:
        return pd.Series(dtype=float)
    return S.zscore_rows(rna.loc[present]).mean(axis=0)


# --------------------------------------------------------------------------- #
# Analysis A: CN-driven vs CN-residual pressure
# --------------------------------------------------------------------------- #
def _gene_pressures(edges: pd.DataFrame, matrices: Dict[str, pd.DataFrame], genes: Sequence[str],
                    gate: Optional[pd.Series]) -> Dict[str, pd.DataFrame]:
    """Compute (and optionally AGO-gate) spine pressure for each dosage track."""
    kw = pressure_kwargs()
    out: Dict[str, pd.DataFrame] = {}
    for track, mat in matrices.items():
        p = engine_pressure(edges, mat, genes=list(genes), **kw)
        if gate is not None and not p.empty:
            shared = p.columns.intersection(gate.index)
            p = p[shared].mul(gate.reindex(shared), axis=1)
        out[track] = p
    return out


def attribution_analysis(genes: Sequence[str], *, gated: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Per-gene total / cn / resid pressure coupling with target RNA + CN attribution."""
    expr = D.load_mirna_arms()
    panels, _ = build_cn_panels(expr.index)
    cn_locus = panels["locus"]
    expr = expr.loc[[a for a in expr.index if a in cn_locus.index or True]]  # keep all arms

    edges = load_mirtar_edges(list(genes), resolve_arms=True)
    floor = C.PRESSURE_ABUNDANCE_FLOOR
    edges = filter_edges_by_abundance(edges, expr, floor)

    cn_hat, resid = decompose_expression(expr, cn_locus)
    matrices = {"total": expr, "cn": cn_hat, "resid": resid}

    rna = D.load_rna()
    clinical = D.load_clinical_strata().set_index("participant")
    conf_cols = [c for c in C.CONFOUNDER_NUMERIC if c in clinical.columns]
    aneu_all = clinical[ANEUPLOIDY_COVAR] if ANEUPLOIDY_COVAR in clinical.columns else None
    prolif_all = _proliferation_metagene(rna)
    prolif_all = prolif_all if not prolif_all.empty else None
    gate = compute_ago_gate(rna)["ago_gate"] if gated else None

    pressures = _gene_pressures(edges, matrices, genes, gate)
    common_genes = sorted(set.intersection(*[set(p.index) for p in pressures.values()]) & set(rna.index))
    target_cn = _safe_target_cn(common_genes)

    cohort_rows: List[dict] = []
    ctx_rows: List[dict] = []
    ctx_cols = [c for c in CONTEXT_COLS if c in clinical.columns]
    # Precompute level -> participant index per context once (was re-grouped per gene).
    ctx_levels: Dict[str, Dict[str, pd.Index]] = {}
    for col in ctx_cols:
        s = clinical[col].dropna()
        ctx_levels[col] = {str(level): pd.Index(idx) for level, idx in s.groupby(s).groups.items()}

    for gene in common_genes:
        e = _matrix_row(rna, gene)
        rec: Dict[str, object] = {"gene": gene}
        n_used = np.nan
        for track, p in pressures.items():
            samples = p.columns.intersection(e.dropna().index)
            pv = p.loc[gene, samples]
            ev = e.reindex(samples)
            conf = clinical.reindex(samples)[conf_cols] if conf_cols else None
            tcn = _matrix_row(target_cn, gene) if (target_cn is not None and gene in target_cn.index) else None
            st = _partial(ev, pv, conf, tcn)
            n_used = st["n"]
            rec[f"rho_{track}"] = _r4(st["partial_rho"])
            rec[f"p_{track}"] = st["partial_p"]
            rec[f"sp_rho_{track}"] = _r4(st["spearman_rho"])
            if track in ("total", "cn") and aneu_all is not None:
                st_a = _partial(ev, pv, conf, tcn, aneu=aneu_all.reindex(samples))
                rec[f"rho_{track}_aneu"] = _r4(st_a["partial_rho"])
                rec[f"p_{track}_aneu"] = st_a["partial_p"]
            if track in ("total", "cn") and (aneu_all is not None or prolif_all is not None):
                st_f = _partial(ev, pv, conf, tcn,
                                aneu=aneu_all.reindex(samples) if aneu_all is not None else None,
                                prolif=prolif_all.reindex(samples) if prolif_all is not None else None)
                rec[f"rho_{track}_full"] = _r4(st_f["partial_rho"])
                rec[f"p_{track}_full"] = st_f["partial_p"]
        rec["n"] = int(n_used) if np.isfinite(n_used) else 0
        if rec["n"] >= MIN_N:
            cohort_rows.append(rec)

        # --- within-context (sample subset; pressure already per-sample) ---
        tcn_full = _matrix_row(target_cn, gene) if (target_cn is not None and gene in target_cn.index) else None
        for col in ctx_cols:
            fixed = col == "PAM50_final"
            for level, samples0 in ctx_levels[col].items():
                sub_rec = {"gene": gene, "context": col, "level": level}
                ok = False
                for track, p in pressures.items():
                    samples = p.columns.intersection(samples0).intersection(e.dropna().index)
                    if len(samples) < MIN_N:
                        break
                    conf_s = clinical.reindex(samples)[conf_cols] if conf_cols else None
                    tcn_s = tcn_full.reindex(samples) if tcn_full is not None else None
                    st = _partial(e.reindex(samples), p.loc[gene, samples], conf_s, tcn_s)
                    sub_rec[f"rho_{track}"] = _r4(st["partial_rho"])
                    sub_rec[f"p_{track}"] = st["partial_p"]
                    sub_rec["n"] = int(st["n"])
                    sub_rec["pam50_fixed"] = fixed
                    if track in ("total", "cn") and aneu_all is not None:
                        st_a = _partial(e.reindex(samples), p.loc[gene, samples], conf_s, tcn_s,
                                        aneu=aneu_all.reindex(samples))
                        sub_rec[f"rho_{track}_aneu"] = _r4(st_a["partial_rho"])
                        sub_rec[f"p_{track}_aneu"] = st_a["partial_p"]
                    if track in ("total", "cn") and (aneu_all is not None or prolif_all is not None):
                        st_f = _partial(e.reindex(samples), p.loc[gene, samples], conf_s, tcn_s,
                                        aneu=aneu_all.reindex(samples) if aneu_all is not None else None,
                                        prolif=prolif_all.reindex(samples) if prolif_all is not None else None)
                        sub_rec[f"rho_{track}_full"] = _r4(st_f["partial_rho"])
                        sub_rec[f"p_{track}_full"] = st_f["partial_p"]
                    ok = True
                if ok and sub_rec.get("n", 0) >= MIN_N:
                    ctx_rows.append(sub_rec)

    cohort = pd.DataFrame(cohort_rows)
    ctx = pd.DataFrame(ctx_rows)
    cohort = _finalize_attribution(cohort, group=None)
    if not ctx.empty:
        ctx = pd.concat(
            [_finalize_attribution(g, group=(c, l)) for (c, l), g in ctx.groupby(["context", "level"])],
            ignore_index=True,
        )
    return cohort, ctx


def _finalize_attribution(df: pd.DataFrame, group) -> pd.DataFrame:
    """Add FDR per track and the per-gene CN-attribution metric."""
    if df.empty:
        return df
    df = df.copy()
    for track in ("total", "cn", "resid", "total_aneu", "cn_aneu", "total_full", "cn_full"):
        pcol = f"p_{track}"
        if pcol in df.columns:
            df[f"q_{track}"] = S.bh_fdr(df[pcol].fillna(1.0).values)
    # attribution = share of the (negative) total coupling reproduced by the CN track.
    rt, rc = df.get("rho_total"), df.get("rho_cn")
    if rt is not None and rc is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = rc / rt
        # only meaningful where total is repressor-consistent (negative) and non-trivial
        frac = frac.where((rt < -0.05), np.nan).clip(lower=-1.0, upper=2.0)
        df["cn_attribution"] = frac.round(4)
        df["cn_driven"] = (rc < 0) & (rc <= 0.5 * rt)   # cn track keeps >=half the total negative rho
    # aneuploidy-adjusted attribution: does CN-driven coupling survive removing the CIN tide?
    rta, rca = df.get("rho_total_aneu"), df.get("rho_cn_aneu")
    if rta is not None and rca is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_a = (rca / rta).where((rta < -0.05), np.nan).clip(lower=-1.0, upper=2.0)
        df["cn_attribution_aneu"] = frac_a.round(4)
    # fully-adjusted attribution: survives aneuploidy AND proliferation together.
    rtf, rcf = df.get("rho_total_full"), df.get("rho_cn_full")
    if rtf is not None and rcf is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            frac_f = (rcf / rtf).where((rtf < -0.05), np.nan).clip(lower=-1.0, upper=2.0)
        df["cn_attribution_full"] = frac_f.round(4)
    return df


# --------------------------------------------------------------------------- #
# Analysis B: combined multi-locus CN-load granularity contest
# --------------------------------------------------------------------------- #
def cn_load_granularity(genes: Sequence[str]) -> pd.DataFrame:
    """Per-gene evidence-weighted combined regulator CN load vs target RNA, per granularity."""
    expr = D.load_mirna_arms()
    panels, arm_meta = build_cn_panels(expr.index)
    edges = load_mirtar_edges(list(genes), resolve_arms=True)
    edges = filter_edges_by_abundance(edges, expr, C.PRESSURE_ABUNDANCE_FLOOR)
    edges = edges.assign(w=np.log1p(pd.to_numeric(edges["evidence_score"], errors="coerce").fillna(0.0)))

    rna = D.load_rna()
    clinical = D.load_clinical_strata().set_index("participant")
    conf_cols = [c for c in C.CONFOUNDER_NUMERIC if c in clinical.columns]
    aneu_all = clinical[ANEUPLOIDY_COVAR] if ANEUPLOIDY_COVAR in clinical.columns else None
    prolif_all = _proliferation_metagene(rna)
    prolif_all = prolif_all if not prolif_all.empty else None

    arm_set = set(panels["locus"].index)
    participants = panels["locus"].columns
    target_cn = _safe_target_cn(genes)
    ctx_cols = [c for c in CONTEXT_COLS if c in clinical.columns]
    ctx_levels: Dict[str, Dict[str, pd.Index]] = {}
    for col in ctx_cols:
        s = clinical[col].dropna()
        ctx_levels[col] = {str(level): participants.intersection(pd.Index(idx))
                           for level, idx in s.groupby(s).groups.items()}

    rows: List[dict] = []
    for gene, grp in edges.groupby("gene"):
        if gene not in rna.index:
            continue
        present = grp.loc[grp["miRNA"].isin(arm_set)]
        if present.empty:
            continue
        arms = present["miRNA"].values
        w = present["w"].values
        e = _matrix_row(rna, gene)
        conf = clinical.reindex(participants)[conf_cols] if conf_cols else None
        tcn = _matrix_row(target_cn, gene) if (target_cn is not None and gene in target_cn.index) else None

        loads: Dict[str, pd.Series] = {}
        for gran in ("locus", "chrom", "cluster"):
            sub = panels[gran].reindex(arms)                       # arms x participant
            delta = sub.sub(DIPLOID)
            loads[gran] = delta.mul(w, axis=0).sum(axis=0, min_count=1)   # combined weighted load
        # single most-amplified regulator (by mean CN delta), locus CN
        mean_delta = panels["locus"].reindex(arms).mean(axis=1) - DIPLOID
        best_arm = mean_delta.idxmax() if mean_delta.notna().any() else None
        loads["single_best"] = (
            panels["locus"].loc[best_arm] if best_arm is not None else pd.Series(np.nan, index=participants)
        )

        rec = {"gene": gene, "n_regulators": int(len(arms)),
               "n_in_cluster": int(arm_meta.reindex(arms)["in_cluster"].fillna(False).sum()),
               "best_arm": best_arm}
        n_used = 0
        for gran, load in loads.items():
            st = _partial(e.reindex(participants), load.reindex(participants), conf, tcn)
            rec[f"rho_{gran}"] = _r4(st["partial_rho"])
            rec[f"p_{gran}"] = st["partial_p"]
            rec[f"sp_rho_{gran}"] = _r4(st["spearman_rho"])
            n_used = max(n_used, st["n"])
            if gran in ("locus", "chrom") and aneu_all is not None:
                st_a = _partial(e.reindex(participants), load.reindex(participants), conf, tcn,
                                aneu=aneu_all.reindex(participants))
                rec[f"rho_{gran}_aneu"] = _r4(st_a["partial_rho"])
                rec[f"p_{gran}_aneu"] = st_a["partial_p"]
            if gran in ("locus", "chrom") and (aneu_all is not None or prolif_all is not None):
                st_f = _partial(e.reindex(participants), load.reindex(participants), conf, tcn,
                                aneu=aneu_all.reindex(participants) if aneu_all is not None else None,
                                prolif=prolif_all.reindex(participants) if prolif_all is not None else None)
                rec[f"rho_{gran}_full"] = _r4(st_f["partial_rho"])
                rec[f"p_{gran}_full"] = st_f["partial_p"]
        rec["n"] = int(n_used)

        # per-context (combined multi-locus load only, to keep table compact)
        for col in ctx_cols:
            for level, samples in ctx_levels[col].items():
                if len(samples) < MIN_N:
                    continue
                cv = conf.reindex(samples) if conf is not None else None
                st = _partial(e.reindex(samples), loads["locus"].reindex(samples), cv,
                              tcn.reindex(samples) if tcn is not None else None)
                rec[f"rho_locus__{col}={level}"] = _r4(st["partial_rho"])
                rec[f"p_locus__{col}={level}"] = st["partial_p"]
        if rec["n"] >= MIN_N:
            rows.append(rec)

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for gran in ("locus", "chrom", "cluster", "single_best",
                 "locus_aneu", "chrom_aneu", "locus_full", "chrom_full"):
        pcol = f"p_{gran}"
        if pcol in df.columns:
            df[f"q_{gran}"] = S.bh_fdr(df[pcol].fillna(1.0).values)
    return df.sort_values("rho_locus")


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #
def _attribution_summary(cohort: pd.DataFrame, ctx: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def _summ(df: pd.DataFrame, label: str):
        if df.empty:
            return
        neg_tot = df[(df["rho_total"] < 0) & (df["q_total"] < C.FDR_ALPHA)]
        neg_cn = df[(df["rho_cn"] < 0) & (df["q_cn"] < C.FDR_ALPHA)]
        neg_res = df[(df["rho_resid"] < 0) & (df["q_resid"] < C.FDR_ALPHA)]
        has_aneu = "q_cn_aneu" in df.columns
        has_full = "q_cn_full" in df.columns
        neg_cn_aneu = (
            df[(df["rho_cn_aneu"] < 0) & (df["q_cn_aneu"] < C.FDR_ALPHA)] if has_aneu else df.iloc[0:0]
        )
        neg_cn_full = (
            df[(df["rho_cn_full"] < 0) & (df["q_cn_full"] < C.FDR_ALPHA)] if has_full else df.iloc[0:0]
        )
        rows.append({
            "stratum": label,
            "n_genes_tested": int(len(df)),
            "n_neg_total": int(len(neg_tot)),
            "n_neg_cn": int(len(neg_cn)),
            "n_neg_resid": int(len(neg_res)),
            "n_neg_cn_aneu_adj": int(len(neg_cn_aneu)) if has_aneu else np.nan,
            "n_neg_cn_full_adj": int(len(neg_cn_full)) if has_full else np.nan,
            "n_neg_total_that_are_cn_driven": int(neg_tot["cn_driven"].sum()) if "cn_driven" in neg_tot else 0,
            "frac_total_neg_cn_driven": round(float(neg_tot["cn_driven"].mean()), 4) if len(neg_tot) else np.nan,
            "median_cn_attribution_among_total_neg": round(float(neg_tot["cn_attribution"].median()), 4)
            if len(neg_tot) else np.nan,
            "median_cn_attribution_aneu_adj": round(float(neg_tot["cn_attribution_aneu"].median()), 4)
            if (has_aneu and len(neg_tot) and "cn_attribution_aneu" in neg_tot) else np.nan,
            "median_cn_attribution_full_adj": round(float(neg_tot["cn_attribution_full"].median()), 4)
            if (has_full and len(neg_tot) and "cn_attribution_full" in neg_tot) else np.nan,
        })

    _summ(cohort, "cohort")
    if not ctx.empty:
        for (col, lvl), g in ctx.groupby(["context", "level"]):
            _summ(g, f"{col}={lvl}")
    return pd.DataFrame(rows)


def _load_summary(load: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grans = ["single_best", "locus", "chrom", "cluster",
             "locus_aneu", "chrom_aneu", "locus_full", "chrom_full"]
    for gran in grans:
        if f"rho_{gran}" not in load.columns or f"q_{gran}" not in load.columns:
            continue
        neg = load[(load[f"rho_{gran}"] < 0) & (load[f"q_{gran}"] < C.FDR_ALPHA)]
        rows.append({
            "granularity": gran,
            "aneuploidy_adjusted": gran.endswith("_aneu"),
            "n_genes_tested": int(len(load)),
            "n_neg_fdr": int(len(neg)),
            "median_rho_neg": round(float(neg[f"rho_{gran}"].median()), 4) if len(neg) else np.nan,
            "median_rho_all": round(float(load[f"rho_{gran}"].median()), 4),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# small utils
# --------------------------------------------------------------------------- #
def _r4(v) -> float:
    return round(float(v), 4) if v is not None and np.isfinite(v) else np.nan


def _safe_target_cn(genes: Sequence[str]) -> Optional[pd.DataFrame]:
    try:
        return D.load_cnv_target_genes(list(genes))
    except Exception as exc:  # noqa: BLE001
        print(f"[cn_dosage_attribution] target CN unavailable; CN partials skipped: {exc}")
        return None


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(*, out_dir: Path = C.CN_DOSAGE_ATTRIBUTION_DIR, genes: Optional[Sequence[str]] = None,
        gated: bool = True) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    gene_list = list(genes or hs.universe)

    print(f"[cn_dosage_attribution] A: CN-driven vs residual pressure for {len(gene_list):,} genes ...")
    cohort, ctx = attribution_analysis(gene_list, gated=gated)
    cohort.to_csv(out_dir / "cn_attribution_cohort.tsv", sep="\t", index=False)
    ctx.to_csv(out_dir / "cn_attribution_by_context.tsv", sep="\t", index=False)
    attr_summary = _attribution_summary(cohort, ctx)
    attr_summary.to_csv(out_dir / "cn_attribution_summary.tsv", sep="\t", index=False)

    print(f"[cn_dosage_attribution] B: combined multi-locus CN-load granularity contest ...")
    load = cn_load_granularity(gene_list)
    load.to_csv(out_dir / "cn_load_granularity_cohort.tsv", sep="\t", index=False)
    load_summary = _load_summary(load) if not load.empty else pd.DataFrame()
    load_summary.to_csv(out_dir / "cn_load_granularity_summary.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.cn_dosage_attribution",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pressure": method_blurb(),
        "gated": gated,
        "n_genes": len(gene_list),
        "confounders_baseline": list(C.CONFOUNDER_NUMERIC) + ["target_cn"],
        "broad_aneuploidy_covar": ANEUPLOIDY_COVAR,
        "proliferation_metagene_sets": list(PROLIFERATION_SETS),
        "aneuploidy_adjusted_columns": "*_aneu (CPE+HRD+target_cn+aneuploidy)",
        "fully_adjusted_columns": "*_full (CPE+HRD+target_cn+aneuploidy+proliferation)",
        "contexts": list(CONTEXT_COLS),
        "granularities": ["single_best", "locus", "chrom", "cluster"],
        "min_fit": MIN_FIT,
        "min_n": MIN_N,
        "attribution_cohort_summary": attr_summary.to_dict("records"),
        "load_granularity_summary": load_summary.to_dict("records") if not load_summary.empty else [],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\n=== CN attribution (cohort + contexts) ===")
    print(attr_summary.to_string(index=False))
    print("\n=== combined CN-load by granularity (cohort) ===")
    if not load_summary.empty:
        print(load_summary.to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=C.CN_DOSAGE_ATTRIBUTION_DIR)
    ap.add_argument("--no-gate", action="store_true", help="Use ungated pressure only")
    ap.add_argument("--max-genes", type=int, default=0, help="Smoke test: cap gene universe")
    args = ap.parse_args()
    genes = None
    if args.max_genes > 0:
        genes = sorted(HallmarkSets.load().universe)[: args.max_genes]
    run(out_dir=args.out_dir, genes=genes, gated=not args.no_gate)


if __name__ == "__main__":
    main()
