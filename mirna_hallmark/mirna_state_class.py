"""miRNA pressure state-trajectory classification (GTEx-healthy -> NAT -> tumor).

Stage 1 of the gain/loss formalisation. Two resolutions are delivered together:

A. EDGE resolution (miRNA -> gene): movement of the pressure a miRNA applies on a
   target gene, as the JOIN of two axes across the three states:
     - pressure-movement (potential): cross-state contribution ``c_X(m,g)`` under the
       no-z ``softmax_logrpm`` attribution mode. Because the edge weight ``w(m,g)`` is
       static across states, the per-edge pressure *direction* is inherited from the
       miRNA's abundance trajectory; the magnitude localises how much of the gained
       budget lands on ``g`` (``edge_share``).
     - coupling-movement (realised): composition+proliferation-adjusted partial
       Spearman per state on **two parallel ladders**:
         (i) miRNA **abundance** vs target expression (raw GTEx bundles);
         (ii) per-edge **realized pressure** ``c(m,g,s)`` vs target expression
         (GTEx miRNA matrix collapse-repaired, same as the pressure axis).
   -> ``joint_edge_class`` in {acquired_realized, field_established_realized,
      acquired_unrealized, constitutive, lost, nat_decoupled, stable, non_monotonic,
      uncoupled}.

B. miRNA resolution (global Hallmark-universe budget): per-state pressure mass /
   share / within-state percentile rank, the rank deltas dHN/dNT/dHT, a primary
   state-trajectory class, a Tau-based subtype-specificity secondary class, and
   special-case flags (incl. odd-NAT). A per-Hallmark-set breakdown is the secondary.

Gain/loss formalised. Pressure on the global budget is
  ``mass_X(m) = sum_{g in HallmarkUniverse} w(m,g) . f_softmax_logrpm(m;X)``
and (GTEx is cross-platform, so contrasts use within-state percentile RANK):
  dHN(m) = rank_nat   - rank_gtex   (healthy -> NAT;  field effect)
  dNT(m) = rank_tumor - rank_nat    (NAT -> tumor;    tumor-specific increment)
  dHT(m) = rank_tumor - rank_gtex   (healthy -> tumor; the acquired axis = headline)
A "gainer" has dHT > 0 (gained repressive pressure over healthy). The healthy
anchor is the **GTEx state, never NAT**.

Run:
  .venv/bin/python3 -m mirna_hallmark.mirna_state_class --healthy-anchor
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
from mirna_hallmark.analyses.cross_state.cross_state_coupling import (
    Q_ALPHA,
    RHO_FLOOR,
    _METAGENE_GENES,
    _couple_state,
    _couple_state_pressure,
    _edge_archetype_row,
    _mirna_for_edge_pressure,
    _nat_decoupled_row,
    _state_covariates,
)
from mirna_hallmark.analyses.cross_state.cross_state_deep_dive import _state_bundles
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _PAM50
from mirna_hallmark.family_normal_reference import FAMILY_ARMS, _participant
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import (
    compute_gene_pressure_contributions,
    load_mirtar_edges,
)
from mirna_hallmark.pressure_engine import compute_edge_pressure_map
from mirna_hallmark.stats import bh_fdr

OUT_DIR = C.TISSUE_REFERENCE_DIR / "mirna_state_class"
ATTR_MODE = C.PRESSURE_ATTRIBUTION_EXPR_MODE       # no-z, cross-state comparable
ANCHOR_MODE = C.PRESSURE_HEALTHY_ANCHOR_MODE        # softmax_devhealthy_logrpm (acquired)
STRUCT_MODE = C.PRESSURE_STRUCT_EXPR_MODE           # softmax: edge_w·sm (no abundance double-count)

NET_THR = 0.15      # |dHT| rank delta for a gain/loss call
MOVE_THR = 0.10     # |dHN| or |dNT| rank delta to count a leg as "moved"
MAG_STRONG = 2.0    # |log2fc tumor-vs-healthy| for a strong magnitude gain/loss (>4-fold)
MAG_MOD = 1.0       # |log2fc| for a moderate magnitude gain/loss (>2-fold)
EDGE_PRESSURE_EPS = 0.0  # edge-level c_tumor - c_ref must exceed this for edge_pressure_gain_*
TAU_PAN = 0.30      # subtype specificity below = pan
TAU_SPECIFIC = 0.50  # subtype specificity at/above = subtype-specific
EXCLUDE_FRAC = 0.25  # one subtype < EXCLUDE_FRAC * median(others) = excluded-in-that-subtype


# --------------------------------------------------------------------------- #
# State matrices + per-state realised contributions
# --------------------------------------------------------------------------- #
def _impute_gtex_state(gtex: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Index]:
    """Repair the GTEx **abundance** matrix in place for arms zeroed by the multi-mapping
    collapse, using the healthy-anchor baseline `gtex_median` (same log2(TPM+1) scale as the
    state matrix; for `mited_anchor` arms it is the miTED rank-transferred value mapped onto
    the GTEx scale). An arm is repaired when it is absent / collapsed (state-matrix median
    <= ABUND_FLOOR) but the baseline carries a healthy reference (`source` in
    {gtex, mited_anchor, gtex_family}). The repaired value is broadcast across GTEx columns
    (level/rank is what the cross-state contribution axis needs; coupling keeps the raw
    matrix). Returns (repaired_matrix, injected_arm_index)."""
    from mirna_hallmark.healthy_anchor import OUT_DIR as HA_OUT, ABUND_FLOOR, load_baseline

    load_baseline()  # ensure the cache exists
    cache = HA_OUT / "gtex_qn_baseline.tsv"
    if not cache.is_file() or gtex is None or gtex.empty:
        return gtex, pd.Index([])
    tab = pd.read_csv(cache, sep="\t", index_col=0)
    gm = pd.to_numeric(tab.get("gtex_median"), errors="coerce")
    src = tab.get("source")
    if gm is None or src is None:
        return gtex, pd.Index([])
    have_ref = src.isin(["gtex", "mited_anchor", "gtex_family"]) & gm.notna() & (gm > ABUND_FLOOR)
    cur_med = gtex.median(axis=1)
    g = gtex.copy()
    injected: List[str] = []
    for arm in gm.index[have_ref]:
        collapsed = (arm not in g.index) or (float(cur_med.get(arm, 0.0)) <= ABUND_FLOOR)
        if collapsed:
            g.loc[arm] = float(gm[arm])   # broadcast across GTEx donor columns
            injected.append(arm)
    return g, pd.Index(injected)


def _contrib_by_state(
    states: Dict[str, pd.DataFrame], genes: Sequence[str], edges: pd.DataFrame,
    *, mode: str = ATTR_MODE,
) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    for st, mat in states.items():
        c = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=mat, expr_mode=mode,  # type: ignore[arg-type]
            resolve_arms=False,
        )
        out[st] = c
    return out


def _mass_table(contrib_by_state: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Per-miRNA global pressure mass / share / within-state percentile rank."""
    masses: Dict[str, pd.Series] = {}
    for st, c in contrib_by_state.items():
        if c is None or c.empty:
            continue
        masses[st] = c.groupby("miRNA")["mean_abs_contribution"].sum()
    if not masses:
        return pd.DataFrame()
    idx = sorted(set().union(*[set(m.index) for m in masses.values()]))
    df = pd.DataFrame(index=pd.Index(idx, name="miRNA"))
    for st, m in masses.items():
        tot = float(m.sum())
        df[f"mass_{st}"] = m.reindex(idx)
        df[f"share_{st}"] = (m / tot).reindex(idx) if tot > 0 else np.nan
        df[f"rank_{st}"] = m.rank(pct=True).reindex(idx)
    return df


# --------------------------------------------------------------------------- #
# miRNA-resolution primary state-trajectory class
# --------------------------------------------------------------------------- #
def _primary_class(s1: float, s2: float, net: float) -> str:
    """s1=dHN, s2=dNT, net=dHT (within-state percentile-rank deltas)."""
    if pd.isna(net):
        return "healthy_unknown"
    if abs(net) < NET_THR:
        return "stable"
    if (not pd.isna(s1) and not pd.isna(s2) and abs(s1) >= MOVE_THR
            and abs(s2) >= MOVE_THR and np.sign(s1) != np.sign(s2)):
        return "non_monotonic"
    direction = "gain" if net > 0 else "loss"
    if pd.isna(s1) or pd.isna(s2):
        return f"net_{direction}"
    s1_move = abs(s1) >= MOVE_THR and np.sign(s1) == np.sign(net)
    s2_move = abs(s2) >= MOVE_THR and np.sign(s2) == np.sign(net)
    if s1_move and s2_move:
        return f"progressive_{direction}"
    if s2_move and not s1_move:
        return f"tumor_acquired_{direction}"
    if s1_move and not s2_move:
        return f"field_effect_{direction}"
    return f"net_{direction}"


def _mirna_trajectory(mass: pd.DataFrame) -> pd.DataFrame:
    df = mass.copy()
    has_g = "rank_gtex" in df.columns
    rg = df["rank_gtex"] if has_g else pd.Series(np.nan, index=df.index)
    rn = df.get("rank_nat", pd.Series(np.nan, index=df.index))
    rt = df.get("rank_tumor", pd.Series(np.nan, index=df.index))
    df["dHN"] = rn - rg
    df["dNT"] = rt - rn
    df["dHT"] = rt - rg
    df["primary_class"] = [
        _primary_class(df["dHN"].iloc[i], df["dNT"].iloc[i], df["dHT"].iloc[i])
        for i in range(len(df))
    ]
    # same-platform NAT axis (always available) for healthy_unknown fallback
    mt = df.get("mass_tumor", pd.Series(np.nan, index=df.index))
    mn = df.get("mass_nat", pd.Series(np.nan, index=df.index))
    nat_axis = np.where(
        (mt - mn) > 0, "tumor_vs_nat_up",
        np.where((mt - mn) < 0, "tumor_vs_nat_down", "tumor_vs_nat_flat"),
    )
    df["nat_axis"] = nat_axis
    # flags
    s1, net = df["dHN"], df["dHT"]
    df["odd_nat_abundance"] = (
        (np.sign(s1) == np.sign(net)) & (net.abs() >= 0.10) & (s1.abs() >= 0.5 * net.abs())
    ).fillna(False)
    df["non_monotonic"] = df["primary_class"].eq("non_monotonic")
    df["pressure_move"] = np.where(
        net >= NET_THR, "gain", np.where(net <= -NET_THR, "loss",
                                         np.where(net.isna(), "unknown", "flat")),
    )
    df["is_family_arm"] = df.index.isin(set(FAMILY_ARMS))
    return df


# --------------------------------------------------------------------------- #
# Subtype specificity (Tau) over PAM50, on the acquired (dev) share when anchored
# --------------------------------------------------------------------------- #
def _subtype_contribs(
    tumor_mat: pd.DataFrame, genes: Sequence[str], edges: pd.DataFrame,
    baseline: Optional[pd.Series],
) -> Dict[str, pd.DataFrame]:
    clin = D.load_clinical_strata()
    if "PAM50_final" not in clin.columns:
        return {}
    part_pam = clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
    col_pam = {c: part_pam.get(_participant(c)) for c in tumor_mat.columns}
    mode = ANCHOR_MODE if baseline is not None else ATTR_MODE
    kw = {"healthy_baseline": baseline} if baseline is not None else {}
    out: Dict[str, pd.DataFrame] = {}
    for sub in _PAM50:
        cols = [c for c in tumor_mat.columns if col_pam.get(c) == sub]
        if len(cols) < 20:
            continue
        contrib = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=tumor_mat[cols], expr_mode=mode,  # type: ignore[arg-type]
            resolve_arms=False, **kw,
        )
        if not contrib.empty:
            out[sub] = contrib
    return out


def _tau(a: pd.DataFrame) -> pd.Series:
    """Yanai tissue-specificity index over subtype columns of share matrix `a`."""
    mx = a.max(axis=1)
    k = a.shape[1]
    tau = (1 - a.div(mx, axis=0)).sum(axis=1) / max(k - 1, 1)
    return tau.where(mx > 0)


def _subtype_secondary(tau: float, n_present: int, minval: float, med: float,
                       dom: str, mn: str) -> str:
    if pd.isna(tau):
        return "na"
    excluded = (n_present >= 3) and (med > 0) and (minval < EXCLUDE_FRAC * med)
    if tau >= TAU_SPECIFIC:
        return f"subtype_specific:{dom}"
    if excluded:
        return f"subtype_excluded:{mn}"
    if tau < TAU_PAN:
        return "pan"
    return "intermediate"


def _subtype_mirna_table(contribs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    shares: Dict[str, pd.Series] = {}
    for sub, c in contribs.items():
        m = c.groupby("miRNA")["mean_abs_contribution"].sum()
        tot = float(m.sum())
        shares[sub] = m / tot if tot > 0 else m
    if not shares:
        return pd.DataFrame()
    a = pd.DataFrame(shares).fillna(0.0)
    sub_cols = list(a.columns)
    tau = _tau(a)
    dom = a.idxmax(axis=1)
    mn = a.idxmin(axis=1)
    minval = a.min(axis=1)
    med = a.median(axis=1)
    n_present = (a > 0).sum(axis=1)
    out = a.rename(columns={s: f"subshare_{s}" for s in sub_cols}).copy()
    out["subtype_tau"] = tau
    out["dominant_subtype"] = dom
    out["min_subtype"] = mn
    out["secondary_class"] = [
        _subtype_secondary(tau.get(i), int(n_present.get(i, 0)), float(minval.get(i, np.nan)),
                           float(med.get(i, np.nan)), dom.get(i), mn.get(i))
        for i in a.index
    ]
    out.index.name = "miRNA"
    return out


def _subtype_edge_table(contribs: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for sub, c in contribs.items():
        s = c.set_index(["miRNA", "gene"])["mean_abs_contribution"].rename(f"subc_{sub}")
        frames.append(s)
    if not frames:
        return pd.DataFrame()
    edge = pd.concat(frames, axis=1)
    sub_cols = [c for c in edge.columns if c.startswith("subc_")]
    a = edge[sub_cols].fillna(0.0)
    a.columns = [c.replace("subc_", "") for c in sub_cols]
    edge["edge_subtype_tau"] = _tau(a).values
    edge["edge_dominant_subtype"] = a.idxmax(axis=1).values
    return edge.reset_index()


# --------------------------------------------------------------------------- #
# Edge-resolution pressure movement
# --------------------------------------------------------------------------- #
def _struct_share_table(struct_by_state: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Per-(miRNA, gene) STRUCTURAL gene-share — the arm's de-double-counted fraction of
    a gene's incoming regulation. Built from the ``softmax`` (edge_w·sm) contribution, so
    absolute abundance (``logrpm``) is NOT spent twice the way it is in ``gene_share_*``.
    Answers who is the gene's preferred/specific regulator, not who is most abundant."""
    if not struct_by_state:
        return pd.DataFrame()
    frames = []
    for st, c in struct_by_state.items():
        if c is None or c.empty:
            continue
        frames.append(c.set_index(["miRNA", "gene"])["mean_abs_contribution"].rename(f"cstruct_{st}"))
    if not frames:
        return pd.DataFrame()
    s = pd.concat(frames, axis=1).reset_index()
    for st in ("gtex", "nat", "tumor"):
        col = f"cstruct_{st}"
        if col in s.columns:
            cv = pd.to_numeric(s[col], errors="coerce")
            gsum = cv.groupby(s["gene"]).transform("sum").replace(0, np.nan)
            s[f"gene_struct_share_{st}"] = cv / gsum
    if "cstruct_tumor" in s.columns:
        s["gene_struct_share_rank_tumor"] = (pd.to_numeric(s["cstruct_tumor"], errors="coerce")
                                             .groupby(s["gene"]).rank(ascending=False, method="min"))
        if "gene_struct_share_gtex" in s.columns:
            s["gene_struct_share_delta_tumor_gtex"] = (
                s["gene_struct_share_tumor"] - s["gene_struct_share_gtex"]
            )
    keep = ["miRNA", "gene"] + [c for c in s.columns if c.startswith("gene_struct_share")]
    return s[keep]


def _struct_vs_abundance_audit(edge: pd.DataFrame) -> pd.DataFrame:
    """Per-gene comparison of the dominant regulator under the abundance-driven
    ``gene_share_tumor`` (old) vs the structural ``gene_struct_share_tumor`` (new).
    Quantifies how often de-double-counting abundance changes the "who owns this gene"
    call. Output: ``struct_vs_abundance_resolution.tsv``."""
    need = {"gene", "miRNA", "gene_share_tumor", "gene_struct_share_tumor"}
    if not need.issubset(edge.columns):
        return pd.DataFrame()
    e = edge[["gene", "miRNA", "gene_share_tumor", "gene_struct_share_tumor"]].dropna(
        subset=["gene_share_tumor", "gene_struct_share_tumor"]
    ).copy()
    rows: List[dict] = []
    for g, grp in e.groupby("gene"):
        if len(grp) < 2:
            continue
        ab = grp.loc[grp["gene_share_tumor"].idxmax()]
        st = grp.loc[grp["gene_struct_share_tumor"].idxmax()]
        ab_struct_rank = int((grp["gene_struct_share_tumor"] > ab["gene_struct_share_tumor"]).sum() + 1)
        st_ab_rank = int((grp["gene_share_tumor"] > st["gene_share_tumor"]).sum() + 1)
        rho = grp["gene_share_tumor"].corr(grp["gene_struct_share_tumor"], method="spearman")
        rows.append({
            "gene": g, "n_regulators": int(len(grp)),
            "abund_dominant_mirna": ab["miRNA"], "struct_dominant_mirna": st["miRNA"],
            "dominant_flipped": bool(ab["miRNA"] != st["miRNA"]),
            "abund_dominant_struct_rank": ab_struct_rank,
            "struct_dominant_abund_rank": st_ab_rank,
            "within_gene_spearman": float(rho) if pd.notna(rho) else np.nan,
        })
    out = pd.DataFrame(rows)
    return out.sort_values("dominant_flipped", ascending=False) if not out.empty else out


def _struct_coupling_assessment(edge: pd.DataFrame) -> pd.DataFrame:
    """Does the STRUCTURAL pick track REALIZED repression (the anti-correlation) any better
    than the abundance pick? Tested at every available level — both coupling ladders
    (arm-abundance ``adj_rho`` and realized-pressure ``adj_rho_pressure``). The headline
    finding (2026-06-27): **NO — share of either flavour is ~orthogonal to realized coupling
    (within-gene ρ≈0)**, so the structural split is a correctness win on the *attribution*
    axis only; realized resolution is owned by ``adj_rho``, not the share. The small drop in
    #1-pick FDR-hit rate is a power artifact (abundance favours higher-dynamic-range arms).
    Output: ``struct_vs_abundance_coupling_assessment.tsv`` (one row per coupling ladder)."""
    from scipy.stats import spearmanr, wilcoxon
    need = {"gene", "miRNA", "gene_share_tumor", "gene_struct_share_tumor"}
    ladders = [("abundance_ladder", "adj_rho_tumor", "adj_q_tumor"),
               ("pressure_ladder", "adj_rho_pressure_tumor", "adj_q_pressure_tumor")]
    if not need.issubset(edge.columns):
        return pd.DataFrame()
    e = edge.copy()
    for _, rho, q in ladders:
        for c in (rho, q):
            if c in e.columns:
                e[c] = pd.to_numeric(e[c], errors="coerce")
    e["gene_share_tumor"] = pd.to_numeric(e["gene_share_tumor"], errors="coerce")
    e["gene_struct_share_tumor"] = pd.to_numeric(e["gene_struct_share_tumor"], errors="coerce")
    rows: List[dict] = []
    for name, rho, q in ladders:
        if rho not in e.columns:
            continue
        align_ab, align_st, pick, flip_ab, flip_st = [], [], [], [], []
        for g, grp in e.groupby("gene"):
            d = grp.dropna(subset=["gene_share_tumor", "gene_struct_share_tumor"])
            if len(d) < 2:
                continue
            ab = d.loc[d["gene_share_tumor"].idxmax()]
            st = d.loc[d["gene_struct_share_tumor"].idxmax()]
            def _hit(r): return bool(pd.notna(r[rho]) and pd.notna(r.get(q)) and r[rho] < 0 and r[q] < 0.05)
            def _neg(r): return bool(pd.notna(r[rho]) and r[rho] < 0)
            pick.append((_hit(ab), _hit(st), _neg(ab), _neg(st), ab[rho], st[rho]))
            if ab["miRNA"] != st["miRNA"] and pd.notna(ab[rho]) and pd.notna(st[rho]):
                flip_ab.append(ab[rho]); flip_st.append(st[rho])
            dd = d.dropna(subset=[rho])
            if len(dd) >= 5 and (-dd[rho]).nunique() >= 3:
                ra = spearmanr(dd["gene_share_tumor"], -dd[rho]).statistic
                rs = spearmanr(dd["gene_struct_share_tumor"], -dd[rho]).statistic
                if pd.notna(ra) and pd.notna(rs):
                    align_ab.append(ra); align_st.append(rs)
        P = pd.DataFrame(pick, columns=["ab_hit", "st_hit", "ab_neg", "st_neg", "ab_rho", "st_rho"])
        fa, fs = np.array(flip_ab), np.array(flip_st)
        wp = float(wilcoxon(fs, fa).pvalue) if len(fa) > 10 else np.nan
        rows.append({
            "ladder": name, "rho_col": rho,
            "n_genes_alignment": len(align_ab),
            "abund_align_median_rho": float(np.median(align_ab)) if align_ab else np.nan,
            "struct_align_median_rho": float(np.median(align_st)) if align_st else np.nan,
            "struct_aligns_better_frac": float(np.mean(np.array(align_st) > np.array(align_ab)))
            if align_ab else np.nan,
            "n_genes_pick": int(len(P)),
            "abund_pick_fdr_hit_frac": float(P["ab_hit"].mean()),
            "struct_pick_fdr_hit_frac": float(P["st_hit"].mean()),
            "abund_pick_neg_frac": float(P["ab_neg"].mean()),
            "struct_pick_neg_frac": float(P["st_neg"].mean()),
            "abund_pick_median_rho": float(P["ab_rho"].median()),
            "struct_pick_median_rho": float(P["st_rho"].median()),
            "n_flip_both_scored": int(len(fa)),
            "flip_abund_median_rho": float(np.median(fa)) if len(fa) else np.nan,
            "flip_struct_median_rho": float(np.median(fs)) if len(fs) else np.nan,
            "flip_struct_more_neg_frac": float(np.mean(fs < fa)) if len(fa) else np.nan,
            "flip_wilcoxon_p": wp,
        })
    return pd.DataFrame(rows)


def _edge_pressure(
    contrib_by_state: Dict[str, pd.DataFrame], mass: pd.DataFrame,
    struct_by_state: Optional[Dict[str, pd.DataFrame]] = None,
) -> pd.DataFrame:
    frames = []
    for st, c in contrib_by_state.items():
        if c is None or c.empty:
            continue
        s = c.set_index(["miRNA", "gene"])["mean_abs_contribution"].rename(f"c_{st}")
        frames.append(s)
    edge = pd.concat(frames, axis=1).reset_index()
    ev = None
    for st, c in contrib_by_state.items():
        if c is not None and not c.empty and "evidence_score" in c.columns:
            ev = c.set_index(["miRNA", "gene"])["evidence_score"]
            break
    if ev is not None:
        edge = edge.merge(ev.rename("evidence_score").reset_index(), on=["miRNA", "gene"], how="left")
    # per-state edge share (the miRNA's fraction of its OWN global budget that lands on g, in
    # each state) + structural degree (how many arms regulate the gene; how many genes the arm
    # targets). edge_share is share-among-the-arm's-targets; n_regulators_gene contextualises how
    # contested the gene is (a gene with 80 regulators vs 1 reads very differently).
    for st in ("gtex", "nat", "tumor"):
        col = f"c_{st}"
        if col in edge.columns and f"mass_{st}" in mass.columns:
            edge[f"mass_{st}_arm"] = edge["miRNA"].map(mass[f"mass_{st}"])
            edge[f"edge_share_{st}"] = (pd.to_numeric(edge[col], errors="coerce")
                                        / edge[f"mass_{st}_arm"])
    if "edge_share_tumor" not in edge.columns and "mass_tumor" in mass.columns:
        edge["mass_tumor_arm"] = edge["miRNA"].map(mass["mass_tumor"])
        edge["edge_share_tumor"] = edge["c_tumor"] / edge["mass_tumor_arm"]
    edge["n_regulators_gene"] = edge["gene"].map(edge.groupby("gene")["miRNA"].nunique())
    edge["n_targets_mirna"] = edge["miRNA"].map(edge.groupby("miRNA")["gene"].nunique())
    # GENE-side regulatory share: this arm's fraction of the GENE's total incoming miRNA pressure
    # in each state (complement of edge_share, which is arm-side). gene_share rising tumor>gtex =
    # the arm tightens its grip on the gene's regulation; gene_share_rank_tumor = the arm's rank
    # among the gene's regulators (1 = dominant regulator).
    for st in ("gtex", "nat", "tumor"):
        col = f"c_{st}"
        if col in edge.columns:
            cvals = pd.to_numeric(edge[col], errors="coerce")
            gsum = cvals.groupby(edge["gene"]).transform("sum").replace(0, np.nan)
            edge[f"gene_share_{st}"] = cvals / gsum
    if "c_tumor" in edge.columns:
        edge["gene_share_rank_tumor"] = (pd.to_numeric(edge["c_tumor"], errors="coerce")
                                         .groupby(edge["gene"]).rank(ascending=False, method="min"))
        if "gene_share_gtex" in edge.columns:
            edge["gene_share_delta_tumor_gtex"] = edge["gene_share_tumor"] - edge["gene_share_gtex"]
    if "c_tumor" in edge.columns and "c_nat" in edge.columns:
        edge["edge_dNT_abs"] = edge["c_tumor"].fillna(0.0) - edge["c_nat"].fillna(0.0)
    # STRUCTURAL gene-share (de-abundance-double-counted): the "who is the preferred
    # regulator" axis, complementary to the abundance-driven gene_share_* above.
    struct = _struct_share_table(struct_by_state) if struct_by_state else pd.DataFrame()
    if not struct.empty:
        edge = edge.merge(struct, on=["miRNA", "gene"], how="left")
    return edge


# --------------------------------------------------------------------------- #
# Edge-resolution coupling movement (NAT separated from GTEx)
# --------------------------------------------------------------------------- #
def _edge_archetype(row: pd.Series) -> str:
    return _edge_archetype_row(row, rho_prefix="adj_rho", q_prefix="adj_q")


def _edge_archetype_pressure(row: pd.Series) -> str:
    return _edge_archetype_row(row, rho_prefix="adj_rho_pressure", q_prefix="adj_q_pressure")


def _nat_decoupled(row: pd.Series) -> bool:
    return _nat_decoupled_row(row, rho_prefix="adj_rho", q_prefix="adj_q")


def _nat_decoupled_pressure(row: pd.Series) -> bool:
    return _nat_decoupled_row(row, rho_prefix="adj_rho_pressure", q_prefix="adj_q_pressure")


def _coupling_table(
    edge_pairs: Sequence[Tuple[str, str]],
    hs: HallmarkSets,
    edges: pd.DataFrame,
    states_mirna: Dict[str, pd.DataFrame],
    *,
    cov_source: str = "estimate_epi",
) -> pd.DataFrame:
    if not edge_pairs:
        return pd.DataFrame()
    genes_needed = sorted(
        {g for _, g in edge_pairs}
        | set(_METAGENE_GENES)
        | set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
        | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))
    )
    bundles = _state_bundles(genes_needed)
    states = list(bundles)
    cov: Dict[str, pd.DataFrame] = {}
    for s, (_mir, rna) in bundles.items():
        try:
            cov[s] = _state_covariates(s, rna, hs, cov_source)
        except Exception as exc:  # pragma: no cover - covariate source fallback
            print(f"[state_class] WARN {s} {cov_source} covariates failed ({exc}); metagene fallback")
            cov[s] = _state_covariates(s, rna, hs, "metagene")

    pair_set = set(edge_pairs)
    edges_sub = edges.loc[
        edges.apply(lambda r: (str(r["miRNA"]), str(r["gene"])) in pair_set, axis=1)
    ].copy()
    press_cache: Dict[str, Dict[tuple[str, str], pd.Series]] = {}
    mir_for_pressure = _mirna_for_edge_pressure(bundles, states_mirna)
    for s in states:
        mir_mat = mir_for_pressure.get(s)
        if mir_mat is None or mir_mat.empty:
            continue
        press_cache[s] = compute_edge_pressure_map(
            edges_sub,
            mir_mat,
            genes=genes_needed,
            expr_mode=ATTR_MODE,  # type: ignore[arg-type]
            target_norm=C.PRESSURE_TARGET_NORM,  # type: ignore[arg-type]
        )

    rows: List[dict] = []
    n = len(edge_pairs)
    for i, (arm, g) in enumerate(edge_pairs):
        row = {"miRNA": arm, "gene": g}
        key = (str(arm), str(g))
        for s in states:
            mir, rna = bundles[s]
            _raw_rho, _raw_p, adj_rho, adj_p, ncoup = _couple_state(mir, rna, cov[s], arm, g)
            row[f"adj_rho_{s}"] = adj_rho
            row[f"adj_p_{s}"] = adj_p
            row[f"n_{s}"] = ncoup
            pser = press_cache.get(s, {}).get(key)
            _, _, adj_pr, adj_pp, npres = _couple_state_pressure(rna, cov[s], g, pser)
            row[f"adj_rho_pressure_{s}"] = adj_pr
            row[f"adj_p_pressure_{s}"] = adj_pp
            row[f"n_pressure_{s}"] = npres
        rows.append(row)
        if (i + 1) % 2500 == 0:
            print(f"[state_class] coupling {i + 1}/{n} edges")
    df = pd.DataFrame(rows)
    for s in states:
        if f"adj_p_{s}" in df.columns:
            df[f"adj_q_{s}"] = bh_fdr(df[f"adj_p_{s}"].fillna(1.0).values)
        if f"adj_p_pressure_{s}" in df.columns:
            df[f"adj_q_pressure_{s}"] = bh_fdr(df[f"adj_p_pressure_{s}"].fillna(1.0).values)
    df["coupling_archetype"] = df.apply(_edge_archetype, axis=1)
    df["nat_decoupled"] = df.apply(_nat_decoupled, axis=1)
    df["coupling_archetype_pressure"] = df.apply(_edge_archetype_pressure, axis=1)
    df["nat_decoupled_pressure"] = df.apply(_nat_decoupled_pressure, axis=1)
    if {"adj_rho_pressure_tumor", "adj_rho_pressure_gtex"}.issubset(df.columns):
        df["delta_tumor_gtex_pressure"] = (
            df["adj_rho_pressure_tumor"] - df["adj_rho_pressure_gtex"]
        )
    from mirna_hallmark.analyses.cross_state.cross_state_coupling import MIN_N
    from mirna_hallmark.gtex_mirna_matrix import gtex_arm_map_table

    amap = gtex_arm_map_table()
    if not amap.empty:
        df = df.merge(
            amap.reset_index().rename(columns={"tcga_arm": "miRNA"}),
            on="miRNA",
            how="left",
        )
    if "gtex_feature" not in df.columns:
        df["gtex_feature"] = ""
    if "gtex_arm_map_kind" not in df.columns:
        df["gtex_arm_map_kind"] = "missing"
    df = df.rename(columns={"gtex_feature": "gtex_arm_mapped"})
    df["gtex_coupling_scored"] = (
        pd.to_numeric(df["n_gtex"], errors="coerce").fillna(0).ge(MIN_N)
        & pd.to_numeric(df["adj_rho_gtex"], errors="coerce").notna()
    )
    return df


def _coupling_state_letter(rho: float, q: float, *, sign: int) -> str:
    """Strict significance letter: R/P/0 at |rho|>RHO_FLOOR and q<Q_ALPHA."""
    if pd.isna(rho) or pd.isna(q):
        return "?"
    if sign < 0 and rho < -RHO_FLOOR and q < Q_ALPHA:
        return "R"
    if sign > 0 and rho > RHO_FLOOR and q < Q_ALPHA:
        return "P"
    return "0"


def _coupling_state_dir(rho: float) -> str:
    """Directional letter (no FDR): r/R/p/0 for strong negative/weak negative/positive/neutral."""
    if pd.isna(rho):
        return "?"
    if rho < -RHO_FLOOR:
        return "r"
    if rho > RHO_FLOOR:
        return "p"
    return "0"


def _annotate_coupling_trajectory(row: pd.Series, *, rho_prefix: str = "adj_rho",
                                  q_prefix: str = "adj_q") -> dict:
    """Per-edge strict + directional GTEx/NAT/tumor coupling patterns."""
    sig_chars: List[str] = []
    dir_chars: List[str] = []
    for st in ("gtex", "nat", "tumor"):
        rho = row.get(f"{rho_prefix}_{st}", np.nan)
        q = row.get(f"{q_prefix}_{st}", 1.0)
        neg = _coupling_state_letter(rho, q, sign=-1)
        pos = _coupling_state_letter(rho, q, sign=1)
        if neg == "R":
            sig_chars.append("R")
        elif pos == "P":
            sig_chars.append("P")
        elif neg == "?" or pos == "?":
            sig_chars.append("?")
        else:
            sig_chars.append("0")
        dir_chars.append(_coupling_state_dir(rho))
    sig = "".join(sig_chars)
    directional = "".join(dir_chars)
    rn = row.get(f"{rho_prefix}_nat", np.nan)
    rt = row.get(f"{rho_prefix}_tumor", np.nan)
    nat_dir_neg = (not pd.isna(rn)) and rn < -RHO_FLOOR
    tum_dir_neg = (not pd.isna(rt)) and rt < -RHO_FLOOR
    nat_sig_neg = sig_chars[1] == "R"
    tum_sig_neg = sig_chars[2] == "R"
    return {
        "coupling_sig_pattern": sig,
        "coupling_dir_pattern": directional,
        "nat_repressive_directional": bool(nat_dir_neg),
        "tumor_repressive_directional": bool(tum_dir_neg),
        "nat_tumor_repressive_concordant": bool(nat_dir_neg and tum_dir_neg),
        "nat_underpowered_repressive": bool(
            nat_dir_neg and not nat_sig_neg and tum_sig_neg
        ),
    }


def _arm_pressure_acquired(row: pd.Series) -> bool:
    """Arm-level acquired pressure: rank OR magnitude gainer (same OR as acquired_gainer)."""
    if bool(row.get("arm_acquired_gainer", False)):
        return True
    return str(row.get("pressure_move", "")) == "gain"


def _joint_edge_class_from(
    row: pd.Series,
    *,
    archetype_col: str = "coupling_archetype",
    nat_dec_col: str = "nat_decoupled",
) -> str:
    arche = row.get(archetype_col)
    if arche is None or (isinstance(arche, float) and pd.isna(arche)) or arche in ("unscored", None):
        return "uncoupled"
    if bool(row.get(nat_dec_col, False)):
        return "nat_decoupled"
    if bool(row.get("non_monotonic", False)):
        return "non_monotonic"
    acq = _arm_pressure_acquired(row)
    ploss = str(row.get("pressure_move", "")) == "loss"
    t_neg = arche in ("constitutive_repressor", "tumor_acquired_repressor", "repression_replaces_coexpr")
    sig_pat = str(row.get("coupling_sig_pattern", ""))
    if acq and arche == "tumor_acquired_repressor":
        return "acquired_realized"
    if acq and arche == "repression_replaces_coexpr":
        return "acquired_realized"
    if (
        acq
        and arche == "constitutive_repressor"
        and sig_pat == "0RR"
        and bool(row.get("gtex_coupling_scored", False))
    ):
        return "field_established_realized"
    if acq and t_neg and arche != "constitutive_repressor":
        return "acquired_realized"
    if acq and not t_neg:
        return "acquired_unrealized"
    if arche == "constitutive_repressor":
        return "constitutive"
    if ploss or arche == "lost_repressor":
        return "lost"
    return "stable"


def _joint_edge_class(row: pd.Series) -> str:
    return _joint_edge_class_from(row)


def _joint_edge_class_pressure(row: pd.Series) -> str:
    return _joint_edge_class_from(
        row,
        archetype_col="coupling_archetype_pressure",
        nat_dec_col="nat_decoupled_pressure",
    )


# --------------------------------------------------------------------------- #
# Per-Hallmark-set secondary breakdown
# --------------------------------------------------------------------------- #
def _by_hallmark(contrib_by_state: Dict[str, pd.DataFrame], hs: HallmarkSets) -> pd.DataFrame:
    states = list(contrib_by_state)
    rows: List[dict] = []
    for sname, members in hs.sets.items():
        mset = set(members)
        masses: Dict[str, pd.Series] = {}
        for st, c in contrib_by_state.items():
            if c is None or c.empty:
                continue
            sub = c.loc[c["gene"].isin(mset)]
            if sub.empty:
                continue
            masses[st] = sub.groupby("miRNA")["mean_abs_contribution"].sum()
        if "tumor" not in masses:
            continue
        ranks = {st: masses[st].rank(pct=True) for st in masses}
        idx = sorted(set().union(*[set(m.index) for m in masses.values()]))
        for m in idx:
            r = {"hallmark": sname, "miRNA": m}
            for st in states:
                r[f"mass_{st}"] = float(masses[st].get(m, np.nan)) if st in masses else np.nan
                r[f"rank_{st}"] = float(ranks[st].get(m, np.nan)) if st in ranks else np.nan
            s1 = r.get("rank_nat", np.nan) - r.get("rank_gtex", np.nan)
            s2 = r.get("rank_tumor", np.nan) - r.get("rank_nat", np.nan)
            net = r.get("rank_tumor", np.nan) - r.get("rank_gtex", np.nan)
            r["dHN"], r["dNT"], r["dHT"] = s1, s2, net
            r["set_class"] = _primary_class(s1, s2, net)
            rows.append(r)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Baseline source / confidence
# --------------------------------------------------------------------------- #
def _baseline_source() -> Tuple[pd.Series, set]:
    from mirna_hallmark.healthy_anchor import OUT_DIR as HA_OUT, BLOOD_MIRNAS, load_baseline

    load_baseline()  # ensure the cache table is built
    cache = HA_OUT / "gtex_qn_baseline.tsv"
    if cache.is_file():
        t = pd.read_csv(cache, sep="\t", index_col=0)
        return t["source"], set(BLOOD_MIRNAS)
    return pd.Series(dtype=object), set(BLOOD_MIRNAS)


def _acquired_magnitude() -> pd.Series:
    """Per-arm tumor-vs-healthy log2 fold-change on the SAME (QN'd TCGA) scale:
    `tcga_tumor_median - healthy_baseline_tcga`. This is the magnitude axis the within-state
    rank delta (`dHT`) saturates on: an arm already high in healthy cannot gain rank even
    when its absolute abundance rises several-fold. Cross-platform-safe because the baseline
    is quantile-normalised onto the TCGA scale. NaN for floor0 (no healthy reference)."""
    from mirna_hallmark.healthy_anchor import OUT_DIR as HA_OUT, load_baseline

    load_baseline()
    cache = HA_OUT / "gtex_qn_baseline.tsv"
    if not cache.is_file():
        return pd.Series(dtype=float)
    t = pd.read_csv(cache, sep="\t", index_col=0)
    tum = pd.to_numeric(t.get("tcga_tumor_median"), errors="coerce")
    heal = pd.to_numeric(t.get("healthy_baseline_tcga"), errors="coerce")
    if tum is None or heal is None:
        return pd.Series(dtype=float)
    return (tum - heal).dropna()


def _mag_class(v: float) -> str:
    if pd.isna(v):
        return "unknown"
    if v > MAG_STRONG:
        return "strong_gain"
    if v > MAG_MOD:
        return "moderate_gain"
    if v < -MAG_STRONG:
        return "strong_loss"
    if v < -MAG_MOD:
        return "moderate_loss"
    return "flat"


def _acq_axis(rank_gain: bool, mag_gain: bool) -> str:
    """Which acquired axis fires: rank (dHT, share-among-regulators, saturating) vs
    magnitude (log2fc, absolute abundance surge)."""
    if rank_gain and mag_gain:
        return "rank+magnitude"
    if rank_gain:
        return "rank_only"
    if mag_gain:
        return "magnitude_only"
    return "none"


# --------------------------------------------------------------------------- #
# Gene-centric acquired pressure (rank genes by acquired total pressure +
# the arms that contribute the gain)
# --------------------------------------------------------------------------- #
def _gene_acquired_table(edge: pd.DataFrame, hs: HallmarkSets, *, top_n: int = 6,
                         src: Optional[pd.Series] = None,
                         blood: Optional[set] = None) -> pd.DataFrame:
    """Per gene: total miRNA pressure per state + acquired gain (tumor vs healthy/NAT) +
    the arms that drive the gain. `acquired_vs_gtex` uses the GTEx healthy anchor (subject
    to the multi-mapping collapse); `acquired_vs_nat` is the collapse-free NAT-anchored
    alternative.

    The floor0 gain (arms with c_gtex<=0) is split by the *healthy-anchor baseline source*
    of the driving arm:
      - **imputable artifact** — arm has positive healthy evidence elsewhere
        (`gtex` / `mited_anchor` / `gtex_family`): the GTEx zero is a multi-mapping-collapse
        / name-join artifact, so the gain is FALSE and is removed in
        `acquired_vs_gtex_corrected`;
      - **true acquired** — arm is `floor0` (no healthy evidence anywhere, not blood): a
        genuine tumor-acquired arm whose gain is kept;
      - **blood** — blood-contaminant arm, gain treated as untrustworthy.
    """
    need = {"gene", "miRNA", "c_tumor", "c_gtex", "c_nat"}
    if not need.issubset(edge.columns):
        return pd.DataFrame()
    e = edge.dropna(subset=["gene", "miRNA"]).copy()
    for col in ("c_tumor", "c_nat", "c_gtex"):
        e[col] = pd.to_numeric(e[col], errors="coerce").fillna(0.0)
    e["gain_gtex"] = e["c_tumor"] - e["c_gtex"]
    e["gain_nat"] = e["c_tumor"] - e["c_nat"]

    # classify each arm's floor0 status from the healthy-anchor baseline source
    src = src if src is not None else pd.Series(dtype=object)
    blood = blood or set()
    e["base_src"] = e["miRNA"].map(src).fillna("floor0")
    e["is_blood"] = e["miRNA"].isin(blood)
    floor0 = e["c_gtex"] <= 0
    has_healthy = e["base_src"].isin(["gtex", "mited_anchor", "gtex_family"])
    e["floor0_imputable"] = floor0 & has_healthy & ~e["is_blood"]
    e["floor0_blood"] = floor0 & e["is_blood"]
    e["floor0_true_acq"] = floor0 & ~has_healthy & ~e["is_blood"]

    # gene -> hallmark membership (a gene may sit in several sets)
    gene2sets: Dict[str, List[str]] = {}
    for hset, members in hs.sets.items():
        for g in members:
            gene2sets.setdefault(g, []).append(hset.replace("HALLMARK_", ""))

    rows: List[dict] = []
    for gene, grp in e.groupby("gene"):
        pt, pn, pg = float(grp["c_tumor"].sum()), float(grp["c_nat"].sum()), float(grp["c_gtex"].sum())
        gain = grp["gain_gtex"]
        pos_gain = gain.clip(lower=0.0).sum()
        floor0_gain = grp.loc[grp["c_gtex"] <= 0, "gain_gtex"].clip(lower=0.0).sum()
        imp_gain = grp.loc[grp["floor0_imputable"], "gain_gtex"].clip(lower=0.0).sum()
        blood_gain = grp.loc[grp["floor0_blood"], "gain_gtex"].clip(lower=0.0).sum()
        true_gain = grp.loc[grp["floor0_true_acq"], "gain_gtex"].clip(lower=0.0).sum()
        top_gain = grp.sort_values("gain_gtex", ascending=False).head(top_n)
        top_tum = grp.sort_values("c_tumor", ascending=False).head(top_n)
        sets = gene2sets.get(gene, [])
        rows.append({
            "gene": gene,
            "n_arms": int(grp["miRNA"].nunique()),
            "pressure_tumor": pt, "pressure_nat": pn, "pressure_gtex": pg,
            "acquired_vs_gtex": pt - pg,
            # remove the imputable-artifact (false) gain: those arms are healthy~tumor expressed
            "acquired_vs_gtex_corrected": (pt - pg) - float(imp_gain),
            "acquired_vs_nat": pt - pn,
            "frac_gain_from_floor0": float(floor0_gain / pos_gain) if pos_gain > 0 else np.nan,
            "frac_gain_imputable_artifact": float(imp_gain / pos_gain) if pos_gain > 0 else np.nan,
            "frac_gain_true_acquired": float(true_gain / pos_gain) if pos_gain > 0 else np.nan,
            "frac_gain_blood": float(blood_gain / pos_gain) if pos_gain > 0 else np.nan,
            "top_gaining_arms": ";".join(f"{r.miRNA}({r.gain_gtex:+.2f})" for _, r in top_gain.iterrows()),
            "top_tumor_arms": ";".join(f"{r.miRNA}({r.c_tumor:.2f})" for _, r in top_tum.iterrows()),
            "n_hallmark_sets": len(sets),
            "hallmark_sets": ";".join(sorted(sets)[:6]),
        })
    out = pd.DataFrame(rows)
    return out.sort_values("acquired_vs_gtex", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def run(
    *, out_dir: Path = OUT_DIR, gene_limit: Optional[int] = None, healthy_anchor: bool = True,
    coupling_top_n: Optional[int] = None, coupling: bool = True, cov_source: str = "estimate_epi",
    impute_gtex: bool = True,
) -> Dict[str, pd.DataFrame]:
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices

    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    genes = sorted(hs.universe)
    if gene_limit:
        genes = genes[:gene_limit]

    states = _state_matrices()
    # Repair the GTEx abundance matrix for collapse-zeroed arms (default ON). This makes the
    # cross-state contribution axis (c_gtex, mass_gtex, dHT, acquired_vs_gtex) collapse-free
    # AT SOURCE; coupling re-derives its own raw GTEx matrix, so it is unaffected.
    gtex_injected: pd.Index = pd.Index([])
    if impute_gtex and "gtex" in states:
        states["gtex"], gtex_injected = _impute_gtex_state(states["gtex"])
        print(f"[state_class] GTEx collapse-repair: injected healthy baseline for "
              f"{len(gtex_injected)} arm(s) (e.g. {', '.join(list(gtex_injected[:6]))})")
    print(f"[state_class] states={list(states)}; "
          f"sizes={{ {', '.join(f'{k}:{v.shape[1]}' for k, v in states.items())} }}")

    edges = load_mirtar_edges(genes, resolve_arms=True)
    contrib_by_state = _contrib_by_state(states, genes, edges)
    # parallel STRUCTURAL contribution (edge_w·sm, no logrpm) for the de-abundance-double-
    # counted "preferred regulator" resolution axis (gene_struct_share_*).
    struct_by_state = _contrib_by_state(states, genes, edges, mode=STRUCT_MODE)

    # ---- B. miRNA resolution ------------------------------------------------ #
    mass = _mass_table(contrib_by_state)
    mirna = _mirna_trajectory(mass)

    baseline = None
    base_src: Optional[pd.Series] = None
    base_blood: set = set()
    if healthy_anchor:
        from mirna_hallmark.healthy_anchor import load_baseline
        baseline = load_baseline()

    sub_contribs = _subtype_contribs(states["tumor"], genes, edges, baseline)
    sub_mirna = _subtype_mirna_table(sub_contribs)
    if not sub_mirna.empty:
        mirna = mirna.join(sub_mirna, how="left")

    # baseline-source confidence flags (the dev/acquired interpretation axis)
    if healthy_anchor:
        src, blood = _baseline_source()
        base_src, base_blood = src, blood
        mirna["baseline_source"] = mirna.index.map(src).fillna("floor0")
        mirna["blood"] = mirna.index.isin(blood)
        mirna["imputed_baseline"] = mirna["baseline_source"].isin(["gtex_family", "mited_anchor"])
        mirna["low_healthy_confidence"] = (
            ~mirna["baseline_source"].eq("gtex") | mirna["blood"]
        )
        # ---- MAGNITUDE acquired axis (complements the saturating rank axis dHT) ---- #
        mag = _acquired_magnitude()
        mirna["log2fc_tumor_vs_healthy"] = mirna.index.map(mag)
        mirna["acquired_magnitude_class"] = mirna["log2fc_tumor_vs_healthy"].apply(_mag_class)
        rank_gain = (mirna["dHT"] > NET_THR).fillna(False)
        mag_gain = (mirna["log2fc_tumor_vs_healthy"] > MAG_STRONG).fillna(False)
        # combined gainer: gained SHARE among regulators (rank) OR rose in ABUNDANCE (magnitude)
        mirna["acquired_gainer"] = rank_gain | mag_gain
        mirna["acquired_axis"] = [_acq_axis(r, m) for r, m in zip(rank_gain, mag_gain)]

    # per-miRNA odd_nat_coupling depends on edge coupling -> filled in after edges
    # ---- A. edge resolution: pressure movement ------------------------------ #
    edge = _edge_pressure(contrib_by_state, mass, struct_by_state=struct_by_state)
    # attach miRNA-level pressure trajectory to each edge
    mcols = ["dHN", "dNT", "dHT", "primary_class", "pressure_move", "non_monotonic"]
    if "acquired_gainer" in mirna.columns:
        mcols += ["acquired_gainer", "acquired_axis"]
    edge = edge.merge(
        mirna[mcols].rename(columns={"dHN": "mirna_dHN", "dNT": "mirna_dNT", "dHT": "mirna_dHT",
                                     "primary_class": "mirna_primary_class",
                                     "acquired_gainer": "arm_acquired_gainer",
                                     "acquired_axis": "arm_acquired_axis"}),
        left_on="miRNA", right_index=True, how="left",
    )

    # ---- A. edge resolution: coupling movement ------------------------------ #
    coupling_df = pd.DataFrame()
    if coupling:
        pairs = list(map(tuple, edge[["miRNA", "gene"]].dropna().values))
        if coupling_top_n:
            sort_key = "c_tumor" if "c_tumor" in edge.columns else "evidence_score"
            top = (edge.dropna(subset=["miRNA", "gene"])
                   .sort_values(sort_key, ascending=False)
                   .groupby("miRNA", group_keys=False).head(coupling_top_n))
            pairs = list(map(tuple, top[["miRNA", "gene"]].values))
        print(f"[state_class] coupling on {len(pairs)} edges "
              f"({'all' if not coupling_top_n else f'top-{coupling_top_n}/arm'})")
        coupling_df = _coupling_table(
            pairs, hs, edges, states, cov_source=cov_source,
        )
        if not coupling_df.empty:
            edge = edge.merge(coupling_df, on=["miRNA", "gene"], how="left")

    if "coupling_archetype" not in edge.columns:
        edge["coupling_archetype"] = "unscored"
        edge["nat_decoupled"] = False
    for col, default in (
        ("gtex_arm_mapped", ""),
        ("gtex_arm_map_kind", "missing"),
        ("gtex_coupling_scored", False),
    ):
        if col not in edge.columns:
            edge[col] = default
    traj_cols = edge.apply(
        lambda r: pd.Series(_annotate_coupling_trajectory(r)), axis=1,
    )
    edge = pd.concat([edge, traj_cols], axis=1)
    if {"c_tumor", "c_gtex", "c_nat"}.issubset(edge.columns):
        edge["edge_pressure_gain_gtex"] = (
            pd.to_numeric(edge["c_tumor"], errors="coerce").fillna(0.0)
            - pd.to_numeric(edge["c_gtex"], errors="coerce").fillna(0.0)
        ) > EDGE_PRESSURE_EPS
        edge["edge_pressure_gain_nat"] = (
            pd.to_numeric(edge["c_tumor"], errors="coerce").fillna(0.0)
            - pd.to_numeric(edge["c_nat"], errors="coerce").fillna(0.0)
        ) > EDGE_PRESSURE_EPS
    edge["joint_edge_class"] = edge.apply(_joint_edge_class, axis=1)
    edge["joint_edge_class_pressure"] = edge.apply(_joint_edge_class_pressure, axis=1)

    # propagate odd_nat_coupling up to the miRNA table
    nat_dec = edge.groupby("miRNA")["nat_decoupled"].any()
    mirna["odd_nat_coupling"] = mirna.index.map(nat_dec).fillna(False)

    # ---- A'. edge x PAM50 (pressure only, tumor) ---------------------------- #
    edge_pam50 = _subtype_edge_table(sub_contribs)

    # ---- secondary: per-Hallmark-set ---------------------------------------- #
    by_hallmark = _by_hallmark(contrib_by_state, hs)

    # ---- gene-centric acquired pressure ------------------------------------- #
    gene_acq = _gene_acquired_table(edge, hs, src=base_src, blood=base_blood)

    # ---- write -------------------------------------------------------------- #
    mirna_out = mirna.reset_index().sort_values("dHT", ascending=False, na_position="last")
    mirna_out.to_csv(out_dir / "mirna_state_class.tsv", sep="\t", index=False)
    edge_sorted = edge.sort_values("edge_share_tumor", ascending=False, na_position="last") \
        if "edge_share_tumor" in edge.columns else edge
    edge_sorted.to_csv(out_dir / "mirna_gene_edge_class.tsv", sep="\t", index=False)
    if not edge_pam50.empty:
        edge_pam50.to_csv(out_dir / "mirna_gene_edge_class_by_pam50.tsv", sep="\t", index=False)
    if not by_hallmark.empty:
        by_hallmark.to_csv(out_dir / "mirna_state_class_by_hallmark.tsv", sep="\t", index=False)
    if not gene_acq.empty:
        gene_acq.to_csv(out_dir / "gene_acquired_pressure.tsv", sep="\t", index=False)
    struct_audit = _struct_vs_abundance_audit(edge)
    if not struct_audit.empty:
        struct_audit.to_csv(out_dir / "struct_vs_abundance_resolution.tsv", sep="\t", index=False)
        flip = float(struct_audit["dominant_flipped"].mean())
        print(f"[state_class] struct-vs-abundance: dominant regulator flips on "
              f"{flip:.0%} of {len(struct_audit)} multi-regulator genes "
              f"(median within-gene Spearman {struct_audit['within_gene_spearman'].median():.2f})")
    struct_coup = _struct_coupling_assessment(edge)
    if not struct_coup.empty:
        struct_coup.to_csv(out_dir / "struct_vs_abundance_coupling_assessment.tsv", sep="\t", index=False)
        for _, r in struct_coup.iterrows():
            print(f"[state_class] coupling-assessment [{r['ladder']}]: "
                  f"#1-pick FDR-sig-repressor abundance {r['abund_pick_fdr_hit_frac']:.1%} vs "
                  f"structural {r['struct_pick_fdr_hit_frac']:.1%}; "
                  f"within-gene align ρ abundance {r['abund_align_median_rho']:+.3f} vs "
                  f"structural {r['struct_align_median_rho']:+.3f} "
                  f"(share ⊥ coupling — structural is an attribution fix, not a coupling fix)")

    manifest = {
        "module": "mirna_hallmark.mirna_state_class",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "states": {k: int(v.shape[1]) for k, v in states.items()},
        "attribution_expr_mode": ATTR_MODE,
        "struct_expr_mode": STRUCT_MODE,
        "anchor_expr_mode": ANCHOR_MODE if healthy_anchor else None,
        "healthy_anchor": healthy_anchor,
        "target_norm": C.PRESSURE_TARGET_NORM,
        "n_genes": len(genes),
        "n_arms": int(mirna.shape[0]),
        "n_edges": int(edge.shape[0]),
        "coupling": bool(coupling),
        "coupling_top_n": coupling_top_n,
        "coupling_cov_source": cov_source if coupling else None,
        "coupling_abundance": "partial Spearman(miRNA abundance, target expr); raw GTEx bundles",
        "coupling_pressure": (
            f"partial Spearman(edge c(m,g), target expr); expr_mode={ATTR_MODE}; "
            f"target_norm={C.PRESSURE_TARGET_NORM}; GTEx miRNA imputed for pressure"
        ),
        "impute_gtex": impute_gtex,
        "thresholds": {"net": NET_THR, "move": MOVE_THR, "tau_pan": TAU_PAN,
                       "tau_specific": TAU_SPECIFIC, "exclude_frac": EXCLUDE_FRAC,
                       "rho_floor": RHO_FLOOR, "q_alpha": Q_ALPHA},
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    _print_summary(mirna, edge)
    return {"mirna": mirna_out, "edge": edge_sorted, "edge_pam50": edge_pam50,
            "by_hallmark": by_hallmark, "gene_acquired": gene_acq}


def _print_summary(mirna: pd.DataFrame, edge: pd.DataFrame) -> None:
    pc = mirna["primary_class"].value_counts().to_dict()
    print(f"\n[state_class] miRNA primary_class counts: {pc}")
    if "secondary_class" in mirna.columns:
        sc = mirna["secondary_class"].value_counts().to_dict()
        print(f"[state_class] subtype secondary_class counts: {sc}")
    if "odd_nat_coupling" in mirna.columns:
        print(f"[state_class] odd_nat_coupling arms: {int(mirna['odd_nat_coupling'].sum())}; "
              f"odd_nat_abundance arms: {int(mirna['odd_nat_abundance'].sum())}")
    top = mirna.sort_values("dHT", ascending=False).head(10)
    print("[state_class] top acquired gainers (dHT):",
          "; ".join(f"{i}({r.dHT:+.2f},{r.primary_class})" for i, r in top.iterrows() if pd.notna(r.dHT)))
    bot = mirna.sort_values("dHT", ascending=True).head(10)
    print("[state_class] top losers (dHT):",
          "; ".join(f"{i}({r.dHT:+.2f},{r.primary_class})" for i, r in bot.iterrows() if pd.notna(r.dHT)))
    if "log2fc_tumor_vs_healthy" in mirna.columns:
        axc = mirna["acquired_axis"].value_counts().to_dict()
        print(f"[state_class] acquired_axis counts: {axc} "
              f"(magnitude_only = abundance gains the saturating dHT rank misses)")
        magtop = mirna[mirna["acquired_axis"] == "magnitude_only"].sort_values(
            "log2fc_tumor_vs_healthy", ascending=False).head(10)
        print("[state_class] top MAGNITUDE-only gainers (log2fc, dHT-invisible):",
              "; ".join(f"{i}({r.log2fc_tumor_vs_healthy:+.1f}log2,dHT{r.dHT:+.2f})"
                        for i, r in magtop.iterrows()))
    if "joint_edge_class" in edge.columns:
        jc = edge["joint_edge_class"].value_counts().to_dict()
        print(f"[state_class] joint_edge_class counts: {jc}")
    if "joint_edge_class_pressure" in edge.columns:
        jcp = edge["joint_edge_class_pressure"].value_counts().to_dict()
        print(f"[state_class] joint_edge_class_pressure counts: {jcp}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--healthy-anchor", action="store_true",
                    help="load the GTEx/miTED healthy baseline for the acquired (dev) "
                         "subtype shares + baseline-source confidence flags")
    ap.add_argument("--gene-limit", type=int, default=None)
    ap.add_argument("--coupling-top-n", type=int, default=None,
                    help="cap coupling annotation to top-N targets per miRNA (default: all edges)")
    ap.add_argument("--no-coupling", action="store_true",
                    help="skip the (heavy) per-state partial-rho coupling annotation")
    ap.add_argument("--cov-source", choices=["estimate", "estimate_epi", "metagene"],
                    default="estimate_epi",
                    help="composition covariates (default: ESTIMATE immune+stromal + epi + prolif)")
    ap.add_argument("--no-impute-gtex", action="store_true",
                    help="disable GTEx collapse-repair (revert to raw collapsed GTEx for the "
                         "c_gtex/dHT/acquired axis; for audit only)")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    run(out_dir=args.out_dir, gene_limit=args.gene_limit, healthy_anchor=args.healthy_anchor,
        coupling_top_n=args.coupling_top_n, coupling=not args.no_coupling, cov_source=args.cov_source,
        impute_gtex=not args.no_impute_gtex)


if __name__ == "__main__":
    main()
