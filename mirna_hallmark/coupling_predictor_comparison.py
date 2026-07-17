"""Does the pressure *construction* change realized coupling, vs raw abundance?

Question (raised in review): per-edge coupling is run on the arm's RAW abundance,
while gene/Hallmark coupling is run on the aggregate PRESSURE (which folds in the
softmax share, the z dynamics, and the static evidence weight). Among those pieces,
only the softmax share (co-regulator normalization) and the z*level product can
change a rank-coupling — the static weight cannot (constant per edge) and z alone is
affine in abundance (rank-identical). This module measures, head-to-head, how the
realized coupling behaves under each predictor construction, at BOTH the edge level
and the gene level.

For every predictor we compute the SAME confounder-adjusted partial Spearman used by
the Hallmark spine: CPE (purity) + HRD + target-gene copy number, BH-FDR within each
(level, predictor) family. Headline = neg-sig (partial_rho < 0 AND q < 0.05).

Predictors
----------
EDGE level, predictor for target gene g vs candidate regulator m:
  - abundance        : the arm's raw log2(RPM+1)                       (current default)
  - z                : within-cohort z of the arm                      (affine -> rank-identical to abundance)
  - softmax          : gene-local softmax share of the arm
  - softmax_logrpm   : share x level
  - softmax_z        : share x dynamics
  - softmax_z_logrpm : share x dynamics x level                        (the spine contribution)
  (edge weight is constant per edge -> invisible to a rank test; target_norm fixed at the spine default)

GENE level, predictor for target gene g vs its aggregate incoming dose:
  - abundance_sum    : sum over regulators of max(log2(RPM+1),0)       (naive: no share, no z, no weight)
  - softmax          : sum over regulators of share
  - softmax_logrpm   : sum of share x level                           (the no-z gene-corepression mode)
  - softmax_z        : sum of share x dynamics
  - softmax_z_logrpm : sum of share x dynamics x level                (the spine pressure)

Outputs under output/coupling_predictor_comparison/.
Ungated (the gate is a per-sample rescale, orthogonal to the predictor question).
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
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import (
    compute_gene_pressure,
    compute_gene_pressure_contributions,
    load_mirtar_edges,
)
from mirna_hallmark.pressure_engine import (
    _softmax_rows,
    competition_share_map as _competition_share_map,
    compute_edge_pressure_map,
    compute_gene_pressure as _engine_gene_pressure,
    expression_multiplier,
    filter_edges_by_abundance,
    gene_pressure_from_share as _gene_pressure_from_share,
    logit_matrix as _logit_matrix,
    mirna_target_degree,
    promisc_discount as _promisc_discount,
    sparsemax_rows as _sparsemax_rows,
)

OUT_DIR = C.OUTPUT_ROOT / "coupling_predictor_comparison"

EDGE_MODES = ["abundance", "z", "softmax", "softmax_logrpm", "softmax_z", "softmax_z_logrpm"]
GENE_ENGINE_MODES = ["softmax", "softmax_logrpm", "softmax_z", "softmax_z_logrpm"]
TARGET_NORM = "evidence_mass"   # spine default; rank-irrelevant at edge level
MIN_N = 40


def _row(mat: pd.DataFrame, key: str) -> pd.Series:
    """One row as numeric Series, median-collapsing duplicate labels."""
    r = mat.loc[key]
    if isinstance(r, pd.DataFrame):
        r = r.apply(pd.to_numeric, errors="coerce").median(axis=0)
    return pd.to_numeric(r, errors="coerce")


def _cov_for_gene(cov_base: pd.DataFrame, target_cn: pd.DataFrame, gene: str) -> pd.DataFrame:
    cov = cov_base.copy()
    if not target_cn.empty and gene in target_cn.index:
        cov = cov.join(_row(target_cn, gene).rename("target_cn"), how="left")
    return cov


def _flags(df: pd.DataFrame) -> pd.DataFrame:
    """BH-FDR within each predictor on the CN-adjusted partial p; neg-sig flag."""
    df = df.copy()
    df["partial_q_cn"] = np.nan
    for pred, idx in df.groupby("predictor").groups.items():
        df.loc[idx, "partial_q_cn"] = S.bh_fdr(df.loc[idx, "partial_p_cn"].fillna(1.0).values)
    df["neg_sig_cn"] = (df["partial_rho_cn"] < 0) & (df["partial_q_cn"] < C.FDR_ALPHA)
    return df


# --------------------------------------------------------------------------- #
# Edge level
# --------------------------------------------------------------------------- #
def edge_level(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    rna: pd.DataFrame,
    cov_base: pd.DataFrame,
    target_cn: pd.DataFrame,
    genes: Sequence[str],
) -> pd.DataFrame:
    # Canonical edge set = the spine contribution map's keys (presence-filtered).
    spine_map = compute_edge_pressure_map(
        edges, mirna, genes=genes, expr_mode="softmax_z_logrpm", target_norm=TARGET_NORM
    )
    edge_keys = list(spine_map.keys())
    print(f"[edge] {len(edge_keys):,} (arm,gene) edges in the canonical set")

    rna_genes = set(rna.index)
    rna_cache: Dict[str, pd.Series] = {}
    cov_cache: Dict[str, pd.DataFrame] = {}

    rows: List[dict] = []
    for mode in EDGE_MODES:
        if mode == "abundance":
            emap = None
        elif mode == "z":
            zmat = expression_multiplier(mirna, "z")
            emap = "z"  # sentinel; predictor = zmat row of the arm
        else:
            emap = compute_edge_pressure_map(
                edges, mirna, genes=genes, expr_mode=mode, target_norm=TARGET_NORM
            )
        n_done = 0
        for (arm, gene) in edge_keys:
            if gene not in rna_genes or arm not in mirna.index:
                continue
            if gene not in rna_cache:
                rna_cache[gene] = _row(rna, gene)
            y = rna_cache[gene]
            if gene not in cov_cache:
                cov_cache[gene] = _cov_for_gene(cov_base, target_cn, gene)
            cov = cov_cache[gene]

            if mode == "abundance":
                x = _row(mirna, arm)
            elif mode == "z":
                x = zmat.loc[arm] if arm in zmat.index else None
                if x is None:
                    continue
                x = pd.to_numeric(x, errors="coerce")
            else:
                x = emap.get((arm, gene))
                if x is None:
                    continue
            st = S.correlation_pair(y, x, cov)
            if st["n"] < MIN_N:
                continue
            rows.append(
                {
                    "predictor": mode,
                    "miRNA": arm,
                    "gene": gene,
                    "n": st["n"],
                    "spearman_rho": st["spearman_rho"],
                    "partial_rho_cn": st["partial_rho"],
                    "partial_p_cn": st["partial_p"],
                }
            )
            n_done += 1
        print(f"[edge] {mode:18s} scored {n_done:,} edges")
    return _flags(pd.DataFrame(rows))


# --------------------------------------------------------------------------- #
# Gene level
# --------------------------------------------------------------------------- #
def _abundance_sum_pressure(edges: pd.DataFrame, mirna: pd.DataFrame, genes: Sequence[str],
                            weights: Optional[Dict[tuple, float]] = None) -> pd.DataFrame:
    """Naive gene-level dose: sum over regulators of max(log2RPM,0). ``weights`` optionally applies
    a per-(arm,gene) aggregation weight (co-movement redundancy discount)."""
    logrpm = mirna.astype(float).clip(lower=0.0)
    gene_set = set(genes)
    out = {}
    for gene, grp in edges.groupby("gene"):
        if gene not in gene_set:
            continue
        arms = [a for a in grp["miRNA"].unique() if a in logrpm.index]
        if not arms:
            continue
        if weights is None:
            out[gene] = logrpm.loc[arms].sum(axis=0)
        else:
            wv = pd.Series({a: float(weights.get((a, gene), 1.0)) for a in arms})
            out[gene] = logrpm.loc[arms].mul(wv, axis=0).sum(axis=0)
    return pd.DataFrame(out).T


def gene_level(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    rna: pd.DataFrame,
    cov_base: pd.DataFrame,
    target_cn: pd.DataFrame,
    genes: Sequence[str],
) -> pd.DataFrame:
    predictors: Dict[str, pd.DataFrame] = {
        "abundance_sum": _abundance_sum_pressure(edges, mirna, genes),
    }
    for mode in GENE_ENGINE_MODES:
        predictors[mode] = compute_gene_pressure(
            genes, edges=edges, expr_mode=mode, target_norm=TARGET_NORM, resolve_arms=False
        )

    rna_genes = set(rna.index)
    rows: List[dict] = []
    for pred, mat in predictors.items():
        n_done = 0
        for gene in mat.index:
            if gene not in rna_genes:
                continue
            y = _row(rna, gene)
            x = pd.to_numeric(mat.loc[gene], errors="coerce")
            cov = _cov_for_gene(cov_base, target_cn, gene)
            st = S.correlation_pair(y, x, cov)
            if st["n"] < MIN_N:
                continue
            rows.append(
                {
                    "predictor": pred,
                    "gene": gene,
                    "n": st["n"],
                    "spearman_rho": st["spearman_rho"],
                    "partial_rho_cn": st["partial_rho"],
                    "partial_p_cn": st["partial_p"],
                }
            )
            n_done += 1
        print(f"[gene] {pred:18s} scored {n_done:,} genes")
    return _flags(pd.DataFrame(rows))


# --------------------------------------------------------------------------- #
# Summaries
# --------------------------------------------------------------------------- #
def _summary(df: pd.DataFrame, level: str, key_cols: List[str], baseline: str) -> pd.DataFrame:
    """Per-predictor: counts, neg-sig, and concordance vs the baseline predictor."""
    base = df.loc[df["predictor"] == baseline].set_index(key_cols)
    base_rho = base["partial_rho_cn"]
    base_negsig = set(base.index[base["neg_sig_cn"].fillna(False)])

    out = []
    for pred, sub in df.groupby("predictor"):
        sub_i = sub.set_index(key_cols)
        joined = pd.concat([sub_i["partial_rho_cn"].rename("rho"), base_rho.rename("base")], axis=1).dropna()
        sign_conc = float((np.sign(joined["rho"]) == np.sign(joined["base"])).mean()) if len(joined) else np.nan
        rho_rankcorr = (
            float(joined["rho"].corr(joined["base"], method="spearman")) if len(joined) > 2 else np.nan
        )
        negsig = set(sub_i.index[sub_i["neg_sig_cn"].fillna(False)])
        recov = len(negsig & base_negsig) / len(base_negsig) if base_negsig else np.nan
        neg_genes = sub_i.index[sub_i["partial_rho_cn"] < 0]
        med_neg = float(sub_i.loc[neg_genes, "partial_rho_cn"].median()) if len(neg_genes) else np.nan
        out.append(
            {
                "level": level,
                "predictor": pred,
                "n_tested": int(len(sub)),
                "n_neg": int((sub["partial_rho_cn"] < 0).sum()),
                "n_neg_sig_cn": int(sub["neg_sig_cn"].sum()),
                "median_partial_rho": round(float(sub["partial_rho_cn"].median()), 4),
                "median_partial_rho_among_neg": round(med_neg, 4) if np.isfinite(med_neg) else np.nan,
                "sign_concordance_vs_base": round(sign_conc, 4) if np.isfinite(sign_conc) else np.nan,
                "rho_rankcorr_vs_base": round(rho_rankcorr, 4) if np.isfinite(rho_rankcorr) else np.nan,
                "recovery_of_base_neg_sig": round(recov, 4) if np.isfinite(recov) else np.nan,
                "baseline": baseline,
            }
        )
    order = {p: i for i, p in enumerate(df["predictor"].unique())}
    return pd.DataFrame(out).sort_values("predictor", key=lambda s: s.map(order))


# --------------------------------------------------------------------------- #
# Edge predictor × conditioning sweep (E4 + the MH-44 predictor head-to-head)
# --------------------------------------------------------------------------- #
COND_PRESSURE_MODES = ["softmax", "softmax_logrpm", "softmax_z", "softmax_z_logrpm"]
COND_PREDICTORS = ["abundance"] + COND_PRESSURE_MODES


def edge_conditioning(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    rna: pd.DataFrame,
    cov_base: pd.DataFrame,
    target_cn: pd.DataFrame,
    genes: Sequence[str],
    *,
    extra_predictors: Optional[Dict[str, Dict]] = None,
    include_pressure_modes: bool = True,
) -> pd.DataFrame:
    """Predictor × conditioning head-to-head at the edge (E4 + the MH-44 follow-up).

    PREDICTOR    ∈ {abundance} ∪ {softmax, softmax_logrpm, softmax_z, softmax_z_logrpm}
    CONDITIONING ∈ {none, precision, attribution} with an ASYMMETRIC covariate policy that
    treats the co-regulators as a TARGET covariate, never a focal-arm confounder:

        Z0      = CPE + HRD + target_cn              (spine confounder set)
        Z_mirna = Z0 \\ {target_cn}                   (the miRNA's own non-co-regulator confounders)
        P_other = gene_total - c(m,g,s)              (canonical spine field; FIXED across predictors)

        none        resid X on Z0,      resid Y on Z0               (reproduces the published baseline)
        precision   resid X on Z_mirna, resid Y on Z0 + P_other     (ALL target covariates on Y; only
                                                                     miRNA confounders on X — co-regs are
                                                                     target noise, not an X confound)
        attribution resid X on Z0+P_other, resid Y on Z0 + P_other  (both-sides; honest unique effect)

    Question (user-driven): does treating co-regulators as target noise (precision) — instead
    of letting the standard partial strip co-regulator-correlated variance from the predictor
    (the share itself) — let any softmax pressure flavour beat raw abundance? Pressure loses
    under the both-sides partial (MH-44); this asks whether the conditioning policy flips it.
    """
    from scipy.stats import spearmanr
    from analysis.utils.common.loaders import residualize

    spine_map = compute_edge_pressure_map(
        edges, mirna, genes=genes, expr_mode="softmax_z_logrpm", target_norm=TARGET_NORM
    )
    gene_total: Dict[str, pd.Series] = {}
    gene_nreg: Dict[str, int] = {}
    for (arm, gene), s in spine_map.items():
        gene_total[gene] = s.copy() if gene not in gene_total else gene_total[gene].add(s, fill_value=0.0)
        gene_nreg[gene] = gene_nreg.get(gene, 0) + 1
    pred_maps: Dict[str, object] = {"abundance": None}
    if include_pressure_modes:
        for mode in COND_PRESSURE_MODES:
            pred_maps[mode] = spine_map if mode == "softmax_z_logrpm" else compute_edge_pressure_map(
                edges, mirna, genes=genes, expr_mode=mode, target_norm=TARGET_NORM
            )
    if extra_predictors:
        pred_maps.update(extra_predictors)
    print(f"[cond] {len(spine_map):,} edges; {len(pred_maps)} predictors; "
          f"P_others=spine field (fixed); asymmetric Z (co-regs target-only)")

    def _sp(a: pd.Series, b: pd.Series) -> tuple:
        j = a.index.intersection(b.index)
        if len(j) < MIN_N:
            return (np.nan, np.nan)
        r, p = spearmanr(a.loc[j], b.loc[j])
        return (float(r), float(p))

    rna_genes = set(rna.index)
    rna_cache: Dict[str, pd.Series] = {}
    cov_cache: Dict[str, pd.DataFrame] = {}
    rows: List[dict] = []
    for (arm, gene), c_edge in spine_map.items():
        if gene not in rna_genes or arm not in mirna.index:
            continue
        if gene not in rna_cache:
            rna_cache[gene] = _row(rna, gene)
        if gene not in cov_cache:
            cov_cache[gene] = _cov_for_gene(cov_base, target_cn, gene)
        y, cov0 = rna_cache[gene], cov_cache[gene]
        p_others = gene_total[gene] - c_edge
        base = pd.concat(
            [y.rename("y"), cov0, p_others.rename("p_others")], axis=1, join="inner"
        ).dropna()
        if len(base) < MIN_N:
            continue
        z0 = list(cov0.columns)
        z_mirna = [c for c in z0 if c != "target_cn"] or z0
        zp = z0 + ["p_others"]
        ry0 = residualize(base["y"], base[z0])           # predictor-independent
        ryp = residualize(base["y"], base[zp])
        idx = base.index
        xa = _row(mirna, arm).reindex(idx)
        corr_xp = float(spearmanr(xa, base["p_others"])[0]) if xa.notna().sum() >= MIN_N else np.nan
        for pred_name, pmap in pred_maps.items():
            xv = _row(mirna, arm) if pred_name == "abundance" else pmap.get((arm, gene))
            if xv is None:
                continue
            dd = pd.concat(
                [pd.to_numeric(xv, errors="coerce").rename("x"), base[z0], base["p_others"]],
                axis=1,
            ).dropna()
            if len(dd) < MIN_N:
                continue
            rx0 = residualize(dd["x"], dd[z0])
            rxm = residualize(dd["x"], dd[z_mirna])
            rxp = residualize(dd["x"], dd[zp])
            trip = {
                "none": _sp(ry0, rx0),
                "precision": _sp(ryp, rxm),
                "attribution": _sp(ryp, rxp),
            }
            for variant, (rho, p) in trip.items():
                rows.append({
                    "predictor": pred_name,
                    "conditioning": variant,
                    "miRNA": arm,
                    "gene": gene,
                    "n": int(len(dd)),
                    "partial_rho_cn": rho,
                    "partial_p_cn": p,
                    "n_regulators": gene_nreg[gene],
                    "corr_x_pothers": corr_xp,
                })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["partial_q_cn"] = np.nan
    for _, idxg in df.groupby(["predictor", "conditioning"]).groups.items():
        df.loc[idxg, "partial_q_cn"] = S.bh_fdr(df.loc[idxg, "partial_p_cn"].fillna(1.0).values)
    df["neg_sig_cn"] = (df["partial_rho_cn"] < 0) & (df["partial_q_cn"] < C.FDR_ALPHA)
    return df


def _conditioning_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Per (predictor, conditioning) cell: counts, neg-sig, and overlap vs the global
    baseline (abundance + none = the published edge headline)."""
    key = ["miRNA", "gene"]
    base = df[(df["predictor"] == "abundance") & (df["conditioning"] == "none")].set_index(key)
    base_rho = base["partial_rho_cn"]
    base_negsig = set(base.index[base["neg_sig_cn"].fillna(False)])
    out = []
    for (pred, var), sub in df.groupby(["predictor", "conditioning"]):
        sub_i = sub.set_index(key)
        joined = pd.concat(
            [sub_i["partial_rho_cn"].rename("rho"), base_rho.rename("base")], axis=1
        ).dropna()
        sign_conc = float((np.sign(joined["rho"]) == np.sign(joined["base"])).mean()) if len(joined) else np.nan
        negsig = set(sub_i.index[sub_i["neg_sig_cn"].fillna(False)])
        recov = len(negsig & base_negsig) / len(base_negsig) if base_negsig else np.nan
        neg = sub_i["partial_rho_cn"] < 0
        out.append({
            "predictor": pred,
            "conditioning": var,
            "n_tested": int(len(sub_i)),
            "n_neg": int(neg.sum()),
            "n_neg_sig_cn": int(sub_i["neg_sig_cn"].sum()),
            "median_partial_rho": round(float(sub_i["partial_rho_cn"].median()), 4),
            "median_rho_among_neg": round(float(sub_i.loc[neg, "partial_rho_cn"].median()), 4) if neg.any() else np.nan,
            "sign_conc_vs_abund_none": round(sign_conc, 4) if np.isfinite(sign_conc) else np.nan,
            "recovery_of_abund_none": round(recov, 4) if np.isfinite(recov) else np.nan,
            "new_vs_abund_none": int(len(negsig - base_negsig)),
        })
    res = pd.DataFrame(out)
    pred_order = {p: i for i, p in enumerate(COND_PREDICTORS)}
    cond_order = {"none": 0, "precision": 1, "attribution": 2}
    res["_c"] = res["conditioning"].map(cond_order)
    res["_p"] = res["predictor"].map(pred_order)
    return res.sort_values(["_c", "_p"]).drop(columns=["_c", "_p"]).reset_index(drop=True)


def _conditioning_mechanism(df: pd.DataFrame) -> pd.DataFrame:
    """For the ABUNDANCE predictor: Δρ(precision − none) by terciles of gene crowding and
    focal⊥field independence — tests whether any precision gain concentrates in crowded
    genes whose competitive field is ~independent of the focal arm (more-negative Δρ =
    sharper anti-correlation)."""
    key = ["miRNA", "gene"]
    sub = df[df["predictor"] == "abundance"]
    base = sub.loc[sub["conditioning"] == "none"].set_index(key)["partial_rho_cn"]
    prec = sub.loc[sub["conditioning"] == "precision"].set_index(key)
    j = prec.join(base.rename("base"), how="inner").dropna(subset=["partial_rho_cn", "base"]).copy()
    if j.empty:
        return pd.DataFrame()
    j["delta"] = j["partial_rho_cn"] - j["base"]
    j["abs_corr_xp"] = j["corr_x_pothers"].abs()
    rows = []
    for col, label in [("n_regulators", "crowding"), ("abs_corr_xp", "focal_field_indep")]:
        if j[col].dropna().nunique() < 3:
            continue
        try:
            t = pd.qcut(j[col].rank(method="first"), 3, labels=["low", "mid", "high"])
        except ValueError:
            continue
        for lvl, g in j.groupby(t):
            rows.append({
                "split": label,
                "tercile": str(lvl),
                "n": int(len(g)),
                "median_delta_rho": round(float(g["delta"].median()), 4),
                "frac_more_negative": round(float((g["delta"] < 0).mean()), 4),
            })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Softmax-LOGIT grid (Axis-Group-I share flavours): reference × measurement ×
# promiscuity. Each is a bare gene-local share used as the edge predictor; per the
# coupling-invariance lemma the logit reference/scale move the edge coupling ONLY
# through the softmax, so each logit kind is a genuinely different predictor.
# --------------------------------------------------------------------------- #
LOGIT_REFS = ["cohort", "healthy", "nat"]   # full-cohort per-arm population baselines
LOGIT_MEASURES = ["dev", "z"]               # x-ref, (x-ref)/sd_tumor  (+ reference-free "rawx")


def _ref_baselines(mirna: pd.DataFrame, edges: pd.DataFrame) -> Dict[str, pd.Series]:
    x = mirna.astype(float)
    deg = _promisc_discount(edges, x.index)  # log1p(#targets)
    out: Dict[str, pd.Series] = {
        "cohort_med": x.median(axis=1),
        "cohort_sd": x.std(axis=1).replace(0, np.nan),
        "promisc": deg,                 # default promiscuity = degree (back-compat with pn flag)
        "promisc_degree": deg,
    }
    try:  # evidence-mass promiscuity definition (per-arm log1p evidence mass over its targets)
        from mirna_hallmark.pressure_engine import mirna_mass_denominators
        out["promisc_evidence_mass"] = mirna_mass_denominators(edges, "evidence_mass").reindex(x.index).fillna(0.0)
    except Exception:  # noqa: BLE001
        out["promisc_evidence_mass"] = deg
    try:  # GENOME-WIDE HE (+breast-expr) promiscuity — the corrected denominator (rank-faithful r≈0.92
        # to hallmark, ~1.5x scale; see genomewide_promiscuity + ATTRIBUTION_IDENTITY_VS_MAGNITUDE §11)
        from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
        out["promisc_he_expr"] = GP.load_promiscuity(arms=x.index, col="promisc_he_expr")
    except Exception as e:  # noqa: BLE001
        print(f"[logit] genome-wide HE promiscuity unavailable -> hallmark degree fallback: {e}")
        out["promisc_he_expr"] = deg
    try:  # healthy = GTEx-breast per-arm baseline (imputed)
        from mirna_hallmark import healthy_anchor as HA
        out["healthy_med"] = HA.load_baseline().reindex(x.index)
    except Exception as e:  # noqa: BLE001
        print(f"[logit] healthy baseline unavailable -> cohort fallback: {e}")
        out["healthy_med"] = pd.Series(index=x.index, dtype=float)
    try:  # NAT-cohort per-arm median
        from mirna_hallmark.analyses.cross_state.cross_state_coupling import _states_mirna_for_pressure
        nat = _states_mirna_for_pressure(impute_gtex=False).get("nat")
        out["nat_med"] = nat.median(axis=1).reindex(x.index) if nat is not None else pd.Series(index=x.index, dtype=float)
    except Exception as e:  # noqa: BLE001
        print(f"[logit] NAT matrix unavailable -> cohort fallback: {e}")
        out["nat_med"] = pd.Series(index=x.index, dtype=float)
    return out


def build_logit_predictors(edges: pd.DataFrame, mirna: pd.DataFrame, genes: Sequence[str],
                           baselines: Dict[str, pd.Series]) -> Dict[str, Dict]:
    """Bare-share predictors over reference × measurement × promiscuity (full-cohort)."""
    preds: Dict[str, Dict] = {}
    for promisc_on in (False, True):
        tag = "pn" if promisc_on else "p0"
        preds[f"share|rawx|{tag}"] = _competition_share_map(
            edges, _logit_matrix(mirna, baselines, "rawx", None, promisc_on), genes, "softmax")
        for ref in LOGIT_REFS:
            for measure in LOGIT_MEASURES:
                preds[f"share|{ref}_{measure}|{tag}"] = _competition_share_map(
                    edges, _logit_matrix(mirna, baselines, ref, measure, promisc_on), genes, "softmax")
    return preds


def edge_logit_paired(edges: pd.DataFrame, rna: pd.DataFrame, cov_base: pd.DataFrame,
                      target_cn: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    """Paired-NAT reference logit on the matched tumor∩NAT subset: logit(m,s) =
    x_tumor(m,s) − x_NAT(m, patient(s)) (dev) or /sd_tumor (z), ± promiscuity. Reuses
    ``edge_conditioning`` with mirna = matched tumor (so abundance + P_others are matched-local)
    and the paired shares as extra predictors. Reduced n (~matched pairs) — flagged."""
    from mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels import _matched_pair_matrices
    tumor_mat, nat_mat, matched = _matched_pair_matrices()
    if tumor_mat is None or tumor_mat.empty:
        return pd.DataFrame()
    arms = tumor_mat.index.intersection(nat_mat.index)
    dev = tumor_mat.loc[arms].sub(nat_mat.loc[arms], fill_value=np.nan)
    sd = tumor_mat.loc[arms].std(axis=1).replace(0, np.nan)
    promisc = _promisc_discount(edges, arms)
    paired_preds: Dict[str, Dict] = {}
    for promisc_on in (False, True):
        tag = "pn" if promisc_on else "p0"
        d = dev.sub(promisc, axis=0) if promisc_on else dev
        z = dev.div(sd, axis=0)
        z = z.sub(promisc, axis=0) if promisc_on else z
        paired_preds[f"share|pairednat_dev|{tag}"] = _competition_share_map(edges, d, genes, "softmax")
        paired_preds[f"share|pairednat_z|{tag}"] = _competition_share_map(
            edges, z.replace([np.inf, -np.inf], np.nan), genes, "softmax")
    print(f"[logit] paired-NAT: {tumor_mat.shape[1]} matched pairs, {len(paired_preds)} paired-share predictors")
    return edge_conditioning(edges, tumor_mat, rna, cov_base, target_cn, genes,
                             extra_predictors=paired_preds, include_pressure_modes=False)


def _logit_grid_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Per (sample_set, predictor, conditioning): counts + overlap vs that sample_set's
    abundance+none baseline. Sorted so the strongest predictors surface per block."""
    key = ["miRNA", "gene"]
    cond_order = {"none": 0, "precision": 1, "attribution": 2}
    rows = []
    for sset, dd in df.groupby("sample_set"):
        base = dd[(dd["predictor"] == "abundance") & (dd["conditioning"] == "none")].set_index(key)
        base_negsig = set(base.index[base["neg_sig_cn"].fillna(False)])
        base_rho = base["partial_rho_cn"]
        for (pred, var), sub in dd.groupby(["predictor", "conditioning"]):
            sub_i = sub.set_index(key)
            negsig = set(sub_i.index[sub_i["neg_sig_cn"].fillna(False)])
            joined = pd.concat([sub_i["partial_rho_cn"].rename("r"), base_rho.rename("b")], axis=1).dropna()
            sign_conc = float((np.sign(joined["r"]) == np.sign(joined["b"])).mean()) if len(joined) else np.nan
            neg = sub_i["partial_rho_cn"] < 0
            rows.append({
                "sample_set": sset, "predictor": pred, "conditioning": var,
                "n_tested": int(len(sub_i)), "n_neg": int(neg.sum()),
                "n_neg_sig_cn": int(sub_i["neg_sig_cn"].sum()),
                "median_partial_rho": round(float(sub_i["partial_rho_cn"].median()), 4),
                "median_rho_among_neg": round(float(sub_i.loc[neg, "partial_rho_cn"].median()), 4) if neg.any() else np.nan,
                "recovery_of_abund_none": round(len(negsig & base_negsig) / len(base_negsig), 4) if base_negsig else np.nan,
                "new_vs_abund_none": int(len(negsig - base_negsig)),
                "sign_conc_vs_abund_none": round(sign_conc, 4) if np.isfinite(sign_conc) else np.nan,
            })
    res = pd.DataFrame(rows)
    res["_c"] = res["conditioning"].map(cond_order)
    return res.sort_values(["sample_set", "_c", "n_neg_sig_cn"],
                           ascending=[True, True, False]).drop(columns="_c").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Competition-FUNCTION + occupancy grid (non-softmax normalizers; AGO-gated dose).
# Share/logit construction primitives live in pressure_engine (shared, canonical).
# --------------------------------------------------------------------------- #
def build_function_predictors(edges: pd.DataFrame, mirna: pd.DataFrame, genes: Sequence[str],
                              baselines: Dict[str, pd.Series], gates: Dict[str, pd.Series]) -> Dict[str, Dict]:
    """Competition-function variants on a canonical input + AGO-gated abundance dose(s).
    Log-domain funcs use the cohort-deviance logit; linear/massaction use relative linear
    abundance. Each gate in ``gates`` (full / partial / purified) yields a gated-dose
    predictor = log2(RPM+1) × per-sample gate."""
    x = mirna.astype(float)
    logit = x.sub(baselines["cohort_med"], axis=0)
    rel = (np.power(2.0, x) - 1.0).clip(lower=0.0)
    rel = rel.div(rel.median(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    preds: Dict[str, Dict] = {
        "func|softmax": _competition_share_map(edges, logit, genes, "softmax", 1.0),
        "func|tempT0.5": _competition_share_map(edges, logit, genes, "temp", 0.5),
        "func|tempT2": _competition_share_map(edges, logit, genes, "temp", 2.0),
        "func|sparsemax": _competition_share_map(edges, logit, genes, "sparsemax"),
        "func|linearL1": _competition_share_map(edges, rel, genes, "linear"),
        "func|massaction": _competition_share_map(edges, rel, genes, "massaction"),
    }
    ag = edges.drop_duplicates(["miRNA", "gene"])
    for gname, gser in gates.items():
        g = gser.reindex(x.columns)
        gated: Dict[tuple, pd.Series] = {}
        for arm, gene in zip(ag["miRNA"], ag["gene"]):
            if arm in x.index:
                gated[(str(arm), str(gene))] = x.loc[arm] * g
        preds[gname] = gated
    return preds


# --------------------------------------------------------------------------- #
# GENE-level construction × prune × gating grid (the aggregation lemma flips:
# construction IS consequential at the aggregate). G2/G3/G8/G8b of the gene taxonomy.
# --------------------------------------------------------------------------- #
GENE_PRUNES = ["all", "evid_tertile", "evid_decile", "abund_floor", "topk"]


def _prune_edges(edges: pd.DataFrame, mirna: pd.DataFrame, method: str) -> pd.DataFrame:
    """Regulator-set prioritization for R(g) + the aggregate."""
    if method == "all":
        return edges
    ev = pd.to_numeric(edges["evidence_score"], errors="coerce")
    if method == "evid_tertile":
        return edges[ev >= ev.quantile(2 / 3)]
    if method == "evid_decile":
        return edges[ev >= ev.quantile(0.9)]
    if method == "abund_floor":
        return filter_edges_by_abundance(edges, mirna, 2.0)  # arm cohort-median log2(RPM+1) >= 2
    if method == "topk":
        med = mirna.median(axis=1)
        e = edges.copy()
        e["_ab"] = e["miRNA"].map(med).fillna(-np.inf)
        return e.sort_values("_ab", ascending=False).groupby("gene").head(5).drop(columns="_ab")
    raise ValueError(f"unknown prune: {method}")


def _rel_to_ref(mirna: pd.DataFrame, baselines: Dict[str, pd.Series], ref: str) -> pd.DataFrame:
    """Linear relative abundance for linear/mass-action functions, referenced to ``ref`` median."""
    lin = (np.power(2.0, mirna.astype(float)) - 1.0).clip(lower=0.0)
    if ref == "rawx":
        denom = lin.median(axis=1)
    else:
        center = baselines[f"{ref}_med"].reindex(lin.index).where(
            baselines[f"{ref}_med"].reindex(lin.index).notna(), baselines["cohort_med"])
        denom = (np.power(2.0, center) - 1.0).clip(lower=0.0)
    return lin.div(denom.replace(0, np.nan), axis=0).fillna(0.0)


def build_gene_constructions(edges: pd.DataFrame, mirna: pd.DataFrame, genes: Sequence[str],
                             baselines: Dict[str, pd.Series], z_mat: pd.DataFrame,
                             logrpm: pd.DataFrame, weights: Optional[Dict[tuple, float]] = None,
                             cartesian: bool = False) -> Dict[str, pd.DataFrame]:
    """Carry the edge function × logit-reference × promiscuity space to the GENE aggregate:
    per-arm competition share aggregated as `Σ w_m·share·z·logrpm`. ``weights`` = optional
    co-movement aggregation weight. ``cartesian`` = FULL function × reference × promiscuity-def
    cross (for the overnight grid); else the lighter marginal sweep. The engine's
    `softmax_z_logrpm` ≡ `shr|softmax|cohort|p0`."""
    out: Dict[str, pd.DataFrame] = {}
    if not cartesian:
        for ref in ("cohort", "healthy", "nat", "rawx"):
            for promisc in (False, True):
                tag = "pn" if promisc else "p0"
                logit = _logit_matrix(mirna, baselines, ref, "dev", promisc)
                out[f"shr|softmax|{ref}|{tag}"] = _gene_pressure_from_share(
                    _competition_share_map(edges, logit, genes, "softmax"), z_mat, logrpm, genes, weights)
        logit_c = _logit_matrix(mirna, baselines, "cohort", "dev", False)
        rel = _rel_to_ref(mirna, baselines, "cohort")
        for func, val, T in [("linear", rel, 1.0), ("sparsemax", logit_c, 1.0),
                             ("temp", logit_c, 0.5), ("massaction", rel, 1.0)]:
            out[f"shr|{func}|cohort|p0"] = _gene_pressure_from_share(
                _competition_share_map(edges, val, genes, func, T), z_mat, logrpm, genes, weights)
        return out
    # FULL CARTESIAN: function × reference × promiscuity-def × magnitude × aggregate × contrib.
    # Each share is built ONCE (func×ref×promisc), then aggregated every way (cheap vector ops).
    refs = ("cohort", "healthy", "nat", "rawx")
    pdefs = {"p0": None, "degree": "promisc_degree", "evidence_mass": "promisc_evidence_mass"}
    mags, aggs, contribs = ("logrpm", "z", "z_logrpm"), ("sum", "mean"), ("signed", "pos", "abs")
    for ref in refs:
        base_logit = _logit_matrix(mirna, baselines, ref, "dev", False)  # log-domain centering
        rel_ref = _rel_to_ref(mirna, baselines, ref)                      # linear-domain
        for pdef, key in pdefs.items():
            disc = None if key is None else baselines.get(key)
            log_v = base_logit if disc is None else base_logit.sub(disc.reindex(base_logit.index), axis=0)
            lin_v = rel_ref if disc is None else rel_ref.div(np.exp(disc.reindex(rel_ref.index)).replace(0, np.nan), axis=0).fillna(0.0)
            shares = {f: _competition_share_map(edges, log_v, genes, f, (0.5 if f == "temp" else 1.0))
                      for f in ("softmax", "temp", "sparsemax")}
            shares.update({f: _competition_share_map(edges, lin_v, genes, f, 1.0) for f in ("linear", "massaction")})
            for func, sh in shares.items():
                for mag in mags:
                    for agg in aggs:
                        for con in contribs:
                            out[f"shr|{func}|{ref}|{pdef}|{mag}|{agg}|{con}"] = _gene_pressure_from_share(
                                sh, z_mat, logrpm, genes, weights, magnitude=mag, aggregate=agg, contrib=con)
    return out


def iter_cartesian_constructions(edges, mirna, genes, baselines, z_mat, logrpm, weights=None):
    """Generator form of the FULL cartesian: yields (name, P_agg) ONE at a time so each matrix is
    scored + GC'd before the next is built (bounds memory; the all-at-once dict OOMs — see
    gene-cartesian-oom-gotcha). At most 5 share maps + 1 matrix live at once."""
    refs = ("cohort", "healthy", "nat", "rawx")
    pdefs = {"p0": None, "degree": "promisc_degree", "evidence_mass": "promisc_evidence_mass"}
    mags, aggs, contribs = ("logrpm", "z", "z_logrpm"), ("sum", "mean"), ("signed", "pos", "abs")
    for ref in refs:
        base_logit = _logit_matrix(mirna, baselines, ref, "dev", False)
        rel_ref = _rel_to_ref(mirna, baselines, ref)
        for pdef, key in pdefs.items():
            disc = None if key is None else baselines.get(key)
            log_v = base_logit if disc is None else base_logit.sub(disc.reindex(base_logit.index), axis=0)
            lin_v = rel_ref if disc is None else rel_ref.div(np.exp(disc.reindex(rel_ref.index)).replace(0, np.nan), axis=0).fillna(0.0)
            shares = {f: _competition_share_map(edges, log_v, genes, f, (0.5 if f == "temp" else 1.0))
                      for f in ("softmax", "temp", "sparsemax")}
            shares.update({f: _competition_share_map(edges, lin_v, genes, f, 1.0) for f in ("linear", "massaction")})
            for func, sh in shares.items():
                for mag in mags:
                    for agg in aggs:
                        for con in contribs:
                            yield (f"shr|{func}|{ref}|{pdef}|{mag}|{agg}|{con}",
                                   _gene_pressure_from_share(sh, z_mat, logrpm, genes, weights,
                                                             magnitude=mag, aggregate=agg, contrib=con))
            del shares


def _config_summary_row(cname: str, prune: str, gname: str, rhos: list, ps: list) -> dict:
    """Per-config summary from streamed per-gene (rho, p): BH-FDR over genes → neg-sig count."""
    rho_a = np.asarray(rhos)
    q = S.bh_fdr(np.nan_to_num(np.asarray(ps), nan=1.0))
    neg = rho_a < 0
    return {"predictor": cname, "prune": prune, "gating": gname,
            "n_tested": int(rho_a.size), "n_neg": int(neg.sum()),
            "n_neg_sig_cn": int(((rho_a < 0) & (q < C.FDR_ALPHA)).sum()),
            "median_partial_rho": round(float(np.median(rho_a)), 4),
            "median_rho_among_neg": round(float(np.median(rho_a[neg])), 4) if neg.any() else np.nan}


def _score_one_prune_cartesian(prune: str, edges, mirna, rna, cov_base, target_cn, genes,
                               baselines, z_mat, logrpm, gate_variants, comov: bool = False) -> List[dict]:
    """Memory-bounded cartesian scorer: STREAM constructions, score per gene via the cached fast
    partial, BH-FDR per config over genes, return per-CONFIG SUMMARY rows (NOT 91M per-gene rows)."""
    ped = _prune_edges(edges, mirna, prune)
    print(f"[gene-grid-cart] prune={prune}{' (comov)' if comov else ''}: {len(ped):,} edges")
    weights = _comov_weights(ped, mirna) if comov else None
    cache = _build_gene_cache(rna, cov_base, target_cn, genes, mirna.columns)
    rows: List[dict] = []

    def emit(cname: str, mat) -> None:
        if mat is None or getattr(mat, "empty", True):
            return
        for gname, gser in gate_variants.items():
            matg = mat if gser is None else mat.mul(gser.reindex(mat.columns), axis=1)
            rhos: list = []
            ps: list = []
            for gene in matg.index:
                e = cache.get(gene)
                if e is None:
                    continue
                rho, p, n = _fast_partial(matg.loc[gene], e)
                if not n or not np.isfinite(rho):
                    continue
                rhos.append(rho)
                ps.append(p)
            if rhos:
                rows.append(_config_summary_row(cname, prune, gname, rhos, ps))

    emit("abundance_sum", _abundance_sum_pressure(ped, mirna, genes, weights))
    if not comov:
        for mode in GENE_ANCHOR_MODES:
            emit(mode, _engine_gene_pressure(ped, mirna, genes=list(genes), expr_mode=mode,
                                             target_norm=TARGET_NORM, aggregate="sum", contrib_transform="signed"))
    for cname, mat in iter_cartesian_constructions(ped, mirna, genes, baselines, z_mat, logrpm, weights):
        emit(cname, mat)
    return rows


GENE_ANCHOR_MODES = ["z", "softmax_z", "softmax_z_logrpm"]


def _family_map() -> pd.Series:
    """Arm → TargetScan seed-family label (from the canonical annotated edges). Same-family arms
    share a seed → target the same sites → redundant on a gene."""
    path = C.OUTPUT_ROOT / "edges" / "mirna_hallmark_edges.tsv.gz"
    fam = pd.read_csv(path, sep="\t", usecols=["miRNA", "miRNA_family"]).drop_duplicates("miRNA")
    return fam.set_index("miRNA")["miRNA_family"]


def _collapse_to_families(mirna: pd.DataFrame, edges: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Seed-family collapse (count-once): family expression = abundance-weighted mean of member
    arms' log2(RPM+1); family edges dedup on (family,gene) with max evidence. Arms with no family
    label stay singletons. Treats each seed family as ONE effective regulator throughout."""
    a2f = _family_map()
    fam = pd.Series(mirna.index.map(lambda a: a2f.get(a, a)), index=mirna.index)  # singletons keep arm name
    w = (np.power(2.0, mirna.astype(float)) - 1.0).clip(lower=0.0).mean(axis=1) + 1e-9  # member abundance weight
    num = mirna.astype(float).mul(w, axis=0).groupby(fam).sum()
    den = w.groupby(fam).sum()
    fam_mirna = num.div(den, axis=0)
    e = edges.copy()
    e["miRNA"] = e["miRNA"].map(lambda a: a2f.get(a, a))
    fam_edges = e.groupby(["miRNA", "gene"], as_index=False)["evidence_score"].max()
    return fam_mirna, fam_edges


def _comov_weights(edges: pd.DataFrame, mirna: pd.DataFrame) -> Dict[tuple, float]:
    """Per-(arm,gene) co-movement redundancy weight `w = n_eff_F / |F|`, where F = the arm's
    same-seed-family members within R(g), `n_eff_F = 1 + (|F|−1)(1−r̄)`, and r̄ = mean pairwise
    Spearman among F across samples. Co-moving same-seed arms → down-weighted toward count-once;
    independent same-seed arms → w→1 (count-all). Singletons w=1. Used to weight the gene SUM
    (NOT a collapse) — the only form that expresses partial redundancy."""
    a2f = _family_map()
    have = set(mirna.index)
    ranks = mirna.loc[sorted(have)].rank(axis=1)  # Spearman via Pearson on ranks
    w: Dict[tuple, float] = {}
    for gene, grp in edges.groupby("gene"):
        arms = [a for a in grp["miRNA"].unique() if a in have]
        byfam: Dict[str, list] = {}
        for a in arms:
            byfam.setdefault(a2f.get(a, a), []).append(a)
        for members in byfam.values():
            k = len(members)
            if k == 1:
                w[(members[0], gene)] = 1.0
                continue
            cm = np.corrcoef(ranks.loc[members].values)
            iu = np.triu_indices(k, 1)
            rbar = float(np.nanmean(cm[iu])) if iu[0].size else 0.0
            rbar = min(max(rbar, 0.0), 1.0)  # anti-/uncorrelated same-seed arms → treat as independent
            wm = (1.0 + (k - 1) * (1.0 - rbar)) / k
            for a in members:
                w[(a, gene)] = wm
    return w


def _build_gene_cache(rna: pd.DataFrame, cov_base: pd.DataFrame, target_cn: pd.DataFrame,
                      genes: Sequence[str], samples: pd.Index) -> Dict[str, dict]:
    """Per-gene precompute for the fast partial-Spearman: complete-case samples (y+cov non-NaN,
    restricted to the construction sample space ``samples``), the cov design Z=[1|cov], its
    pseudo-inverse, and the CENTERED ranked y-residual. cov + y are construction-independent, so
    this is built ONCE per job and reused across all ~thousands of (construction × gating) configs
    — removing the dominant per-config OLS re-solve."""
    from scipy.stats import rankdata
    rna_genes = set(rna.index)
    cache: Dict[str, dict] = {}
    for gene in genes:
        if gene not in rna_genes:
            continue
        y = _row(rna, gene)
        cov = _cov_for_gene(cov_base, target_cn, gene)
        df = pd.concat([y.rename("y"), cov], axis=1).dropna()
        df = df.loc[df.index.intersection(samples)]  # align to the construction sample space
        if len(df) < MIN_N:
            continue
        Z = np.column_stack([np.ones(len(df)), df.drop(columns="y").values.astype(float)])
        Zpinv = np.linalg.pinv(Z)
        yv = df["y"].values.astype(float)
        ry = rankdata(yv - Z @ (Zpinv @ yv))
        ry = ry - ry.mean()
        cache[gene] = {"idx": df.index, "Z": Z, "Zpinv": Zpinv, "ry": ry,
                       "ry_norm": float(np.sqrt((ry ** 2).sum())), "n": len(df)}
    return cache


def _fast_partial(x: pd.Series, e: dict) -> tuple:
    """Partial Spearman(x, y | cov) reusing the gene's cached cov-projection + ranked y-residual
    (matches ``correlation_pair``'s residualize-then-Spearman). Dense-x fast path; falls back to
    the general scorer if x has NaNs on the gene's complete-case samples."""
    from scipy.stats import rankdata, t as _t
    xv = pd.to_numeric(x, errors="coerce").reindex(e["idx"]).values.astype(float)
    if np.isnan(xv).any():
        return (np.nan, np.nan, 0)  # caller falls back to correlation_pair (rare; pressures dense)
    xr = xv - e["Z"] @ (e["Zpinv"] @ xv)
    rx = rankdata(xr)
    rx = rx - rx.mean()
    denom = float(np.sqrt((rx ** 2).sum())) * e["ry_norm"]
    if denom == 0:
        return (np.nan, np.nan, e["n"])
    rho = float((rx * e["ry"]).sum() / denom)
    n = e["n"]
    if abs(rho) >= 1.0:
        p = 0.0
    else:
        tstat = rho * np.sqrt((n - 2) / (1.0 - rho * rho))
        p = float(2.0 * _t.sf(abs(tstat), n - 2))
    return (rho, p, n)


def _score_one_prune(prune: str, edges: pd.DataFrame, mirna: pd.DataFrame, rna: pd.DataFrame,
                     cov_base: pd.DataFrame, target_cn: pd.DataFrame, genes: Sequence[str],
                     baselines: Dict[str, pd.Series], z_mat: pd.DataFrame, logrpm: pd.DataFrame,
                     gate_variants: Dict[str, Optional[pd.Series]], comov: bool = False,
                     cartesian: bool = False) -> List[dict]:
    """One prune's full work (build constructions + score all construction × gating × gene).
    Self-contained so it can run in a joblib worker. ``comov`` → co-movement-weighted aggregation
    (P_agg = Σ w_m·c_m, w_m the same-seed redundancy discount on R(g)); skips the engine anchors
    (they aggregate internally and can't take per-arm weights). ``cartesian`` → full function ×
    reference × promiscuity-def construction cross."""
    ped = _prune_edges(edges, mirna, prune)
    print(f"[gene-grid] prune={prune}{' (comov)' if comov else ''}: {len(ped):,} edges")
    weights = _comov_weights(ped, mirna) if comov else None
    constructions: Dict[str, pd.DataFrame] = {"abundance_sum": _abundance_sum_pressure(ped, mirna, genes, weights)}
    if not comov:
        for mode in GENE_ANCHOR_MODES:
            constructions[mode] = _engine_gene_pressure(ped, mirna, genes=list(genes), expr_mode=mode,
                                                        target_norm=TARGET_NORM, aggregate="sum",
                                                        contrib_transform="signed")
    constructions.update(build_gene_constructions(ped, mirna, genes, baselines, z_mat, logrpm, weights, cartesian))
    gene_cache = _build_gene_cache(rna, cov_base, target_cn, genes, mirna.columns)  # built ONCE, reused across all configs
    rows: List[dict] = []
    for cname, mat in constructions.items():
        if mat is None or mat.empty:
            continue
        for gname, gser in gate_variants.items():
            matg = mat if gser is None else mat.mul(gser.reindex(mat.columns), axis=1)
            common = matg.index.intersection(list(gene_cache.keys()))
            for gene in common:
                x = matg.loc[gene]
                rho, p, n = _fast_partial(x, gene_cache[gene])
                if n == 0:  # NaN-x fallback to the general scorer (rare; pressures are dense)
                    e = gene_cache[gene]
                    st = S.correlation_pair(_row(rna, gene), pd.to_numeric(x, errors="coerce"),
                                            _cov_for_gene(cov_base, target_cn, gene), spearman_only=True)
                    rho, p, n = st["partial_rho"], st["partial_p"], st["n"]
                if not n or n < MIN_N or not np.isfinite(rho):
                    continue
                rows.append({"predictor": cname, "prune": prune, "gating": gname, "gene": gene,
                             "n": int(n), "partial_rho_cn": rho, "partial_p_cn": p})
    return rows


GENE_REG_UNITS = ["arm", "seed_family", "seed_comov"]


def gene_construction_grid(edges: pd.DataFrame, mirna: pd.DataFrame, rna: pd.DataFrame,
                           cov_base: pd.DataFrame, target_cn: pd.DataFrame, genes: Sequence[str],
                           *, gates: Dict[str, pd.Series], jobs: int = 1,
                           units: Sequence[str] = ("arm", "seed_family", "seed_comov"),
                           cartesian: bool = False) -> pd.DataFrame:
    """Score `partial-ρ(P_agg(g), y_g | CPE+HRD+target_cn)` over regulator-unit × construction ×
    prune × gating. ``regulator_unit`` spans the seed-redundancy spectrum: **arm** (count-all),
    **seed_family** (collapse same-seed arms → one regulator: share+sum over families, count-once),
    **seed_comov** (arm-level shares, gene SUM co-movement-weighted `Σ w_m·c_m`, w_m the same-seed
    redundancy discount — interpolates count-once↔count-all by within-family co-movement; engine
    anchors skipped for this unit). ``jobs`` parallelizes over prunes."""
    gate_variants: Dict[str, Optional[pd.Series]] = {"ungated": None, **gates}
    # Build all (unit × prune) job specs up front (per-unit inputs prepared once), then parallelize
    # over the FLATTENED list. Cartesian uses the STREAMING summary scorer (memory-bounded; returns
    # per-config summaries, not 91M per-gene rows); marginal uses the per-gene scorer.
    scorer = _score_one_prune_cartesian if cartesian else _score_one_prune
    specs: List[tuple] = []  # (unit, args)
    for unit in units:
        comov = unit == "seed_comov"
        if unit == "seed_family":
            u_mirna, u_edges = _collapse_to_families(mirna, edges)
            print(f"[gene-grid] unit=seed_family: {u_mirna.shape[0]:,} families (from {mirna.shape[0]:,} arms)")
        else:
            u_mirna, u_edges = mirna, edges
        baselines = _ref_baselines(u_mirna, u_edges)
        z_mat = expression_multiplier(u_mirna, "z")
        logrpm = u_mirna.astype(float)
        base_args = (u_edges, u_mirna, rna, cov_base, target_cn, genes, baselines, z_mat, logrpm, gate_variants, comov)
        for p in GENE_PRUNES:
            specs.append((unit, (p,) + base_args + (() if cartesian else (cartesian,))))
    if jobs and jobs != 1:
        try:
            from joblib import Parallel, delayed
            print(f"[gene-grid] parallel over {len(specs)} (unit×prune) jobs (jobs={jobs}; cartesian={cartesian})")
            results = Parallel(n_jobs=jobs)(delayed(scorer)(*spec[1]) for spec in specs)
        except Exception as e:  # noqa: BLE001
            print(f"[gene-grid] parallel failed ({e}); sequential")
            results = [scorer(*spec[1]) for spec in specs]
    else:
        results = [scorer(*spec[1]) for spec in specs]
    frames: List[pd.DataFrame] = []
    for (unit, _), res in zip(specs, results):
        if res:
            u_df = pd.DataFrame(res)
            u_df["regulator_unit"] = unit
            frames.append(u_df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    if cartesian:
        return df  # already per-config summary (n_neg_sig_cn computed per config via streamed FDR)
    df["partial_q_cn"] = np.nan
    for _, idx in df.groupby(["regulator_unit", "predictor", "prune", "gating"]).groups.items():
        df.loc[idx, "partial_q_cn"] = S.bh_fdr(df.loc[idx, "partial_p_cn"].fillna(1.0).values)
    df["neg_sig_cn"] = (df["partial_rho_cn"] < 0) & (df["partial_q_cn"] < C.FDR_ALPHA)
    return df


def _gene_grid_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Per (regulator_unit, predictor, prune, gating): counts + Δ vs the GLOBAL baseline
    (arm / abundance_sum / all / ungated), so seed-family constructions compare to the same naive
    abundance-sum baseline as the arm constructions."""
    key = "gene"
    has_unit = "regulator_unit" in df.columns
    bmask = (df["predictor"] == "abundance_sum") & (df["prune"] == "all") & (df["gating"] == "ungated")
    if has_unit:
        bmask = bmask & (df["regulator_unit"] == "arm")
    base = df[bmask].set_index(key)
    base_neg = set(base.index[base["neg_sig_cn"].fillna(False)])
    group_cols = (["regulator_unit"] if has_unit else []) + ["predictor", "prune", "gating"]
    out = []
    for keyvals, sub in df.groupby(group_cols):
        si = sub.set_index(key)
        neg = si["partial_rho_cn"] < 0
        nss = set(si.index[si["neg_sig_cn"].fillna(False)])
        rec = dict(zip(group_cols, keyvals if isinstance(keyvals, tuple) else (keyvals,)))
        rec.update({
            "n_tested": int(len(si)), "n_neg": int(neg.sum()), "n_neg_sig_cn": int(si["neg_sig_cn"].sum()),
            "median_partial_rho": round(float(si["partial_rho_cn"].median()), 4),
            "median_rho_among_neg": round(float(si.loc[neg, "partial_rho_cn"].median()), 4) if neg.any() else np.nan,
            "recovery_of_base": round(len(nss & base_neg) / len(base_neg), 4) if base_neg else np.nan,
            "new_vs_base": int(len(nss - base_neg)),
        })
        out.append(rec)
    return pd.DataFrame(out).sort_values("n_neg_sig_cn", ascending=False).reset_index(drop=True)


def _gene_role_coherence(edges: pd.DataFrame, mirna: pd.DataFrame, genes: Sequence[str],
                         grid_df: pd.DataFrame) -> pd.DataFrame:
    """Per-gene role (TSG/oncogene) + stack coherence joined to the spine
    (softmax_z_logrpm, all, ungated) net-repression call. Coherence = |Σ signed|/Σ|·| over the
    SIGNED ``softmax_z_logrpm`` contributions (z can flip sign, so a canceling stack ⇒ <1)."""
    from mirna_hallmark import gene_roles
    contrib = compute_gene_pressure_contributions(list(genes), edges=edges, expr_mode="softmax_z_logrpm",
                                                  target_norm=TARGET_NORM, mirna=mirna, resolve_arms=False)
    g = contrib.groupby("gene").agg(total_signed=("mean_signed_contribution", "sum"),
                                    total_abs=("mean_abs_contribution", "sum"),
                                    n_regulators=("miRNA", "size"))
    g["coherence"] = (g["total_signed"].abs() / g["total_abs"].replace(0, np.nan))
    smask = ((grid_df["predictor"] == "softmax_z_logrpm") & (grid_df["prune"] == "all") &
             (grid_df["gating"] == "ungated"))
    if "regulator_unit" in grid_df.columns:
        smask = smask & (grid_df["regulator_unit"] == "arm")
    spine = grid_df[smask].set_index("gene")
    g = g.join(spine[["partial_rho_cn", "neg_sig_cn"]], how="left")
    roles = gene_roles.load_gene_roles(list(genes))
    if "gene" in roles.columns:           # load_gene_roles returns a gene COLUMN, not a gene index
        roles = roles.set_index("gene")
    rcols = [c for c in ("role", "malignancy_sign") if c in roles.columns]
    g = g.join(roles[rcols], how="left")
    return g.reset_index()


def _build_ago_gates(rna: pd.DataFrame, cov_base: pd.DataFrame, hs) -> Dict[str, pd.Series]:
    """Per-sample AGO gate variants for gene-aggregate gating: full / partial-λ / purified
    (capacity residualized on proliferation + CPE + HRD, removing the proliferation confound)."""
    gates: Dict[str, pd.Series] = {}
    try:
        from mirna_hallmark.ago_gate import compute_ago_gate, gate_from_capacity
        from mirna_hallmark.analyses.cross_state.cross_state_coupling import _prolif_metagene
        from analysis.utils.common.loaders import residualize
        agt = compute_ago_gate(rna)
        full_gate = agt["ago_gate"]
        gates["full"] = full_gate
        gates["partial"] = 1.0 + 0.5 * (full_gate - 1.0)
        prolif = _prolif_metagene(rna, hs)
        covp = pd.concat([prolif.rename("prolif"), cov_base], axis=1)
        cap_pure = residualize(agt["ago_capacity_z"], covp.reindex(agt.index))
        if not cap_pure.empty:
            gates["purified"] = gate_from_capacity(cap_pure)
    except Exception as e:  # noqa: BLE001
        print(f"[gene-grid] AGO gates unavailable -> ungated only: {e}")
    return gates


def _run_gene_grid(out_dir: Path, edges, mirna, rna, cov_base, target_cn, pressure_genes,
                   gg_gates: Dict[str, pd.Series], jobs: int, cartesian: bool = False) -> pd.DataFrame:
    """Run + persist the gene construction grid (comparison, summary, role/coherence overlay)."""
    gene_grid_df = gene_construction_grid(edges, mirna, rna, cov_base, target_cn, pressure_genes,
                                          gates=gg_gates, jobs=jobs, cartesian=cartesian)
    if gene_grid_df.empty:
        return gene_grid_df
    if cartesian:
        # gene_grid_df is ALREADY the per-config summary (streamed; no 91M-row per-gene table,
        # no role/coherence overlay which needs per-gene). Sort + write directly.
        gg_summ = gene_grid_df.sort_values(["regulator_unit", "n_neg_sig_cn"], ascending=[True, False]).reset_index(drop=True)
        gg_summ.to_csv(out_dir / "gene_construction_grid_summary.tsv", sep="\t", index=False)
        print(f"\n===== CARTESIAN GENE GRID DONE — {len(gg_summ)} configs (predictor×prune×gating×unit) =====")
        with pd.option_context("display.max_columns", None, "display.width", 240):
            print(gg_summ.nlargest(20, "n_neg_sig_cn").to_string(index=False))
        return gene_grid_df
    gene_grid_df.to_csv(out_dir / "gene_construction_grid_comparison.tsv.gz", sep="\t", index=False)
    gg_summ = _gene_grid_summary(gene_grid_df)
    gg_summ.to_csv(out_dir / "gene_construction_grid_summary.tsv", sep="\t", index=False)
    try:
        rc = _gene_role_coherence(edges, mirna, pressure_genes, gene_grid_df)
        rc.to_csv(out_dir / "gene_construction_role_coherence.tsv", sep="\t", index=False)
    except Exception as e:  # noqa: BLE001
        print(f"[gene-grid] role/coherence overlay skipped: {e}")
    print("\n===== GENE CONSTRUCTION GRID — top configs by neg-sig (baseline=abundance_sum/all/ungated) =====")
    with pd.option_context("display.max_columns", None, "display.width", 240):
        print(gg_summ.head(20).to_string(index=False))
    return gene_grid_df


def _load_grid_inputs(genes: Sequence[str] | None):
    """Shared data load for the grids: edges, mirna, rna, cov_base (CPE+HRD+batch), target_cn."""
    hs = HallmarkSets.load()
    gene_list = list(genes or hs.universe)
    edges = load_mirtar_edges(gene_list, resolve_arms=True)
    mirna = D.load_mirna_arms()
    rna = D.load_rna()
    clinical = D.load_clinical_strata().set_index("participant")
    conf_cols = [c for c in C.CONFOUNDER_NUMERIC if c in clinical.columns]
    cov_base = D.augment_tcga_batch(clinical[conf_cols].apply(pd.to_numeric, errors="coerce"))
    pressure_genes = sorted(set(edges["gene"]) & set(rna.index))
    target_cn = D.load_cnv_target_genes(pressure_genes)
    return hs, edges, mirna, rna, cov_base, conf_cols, pressure_genes, target_cn


def run_gene_grid_only(out_dir: Path = OUT_DIR, genes: Sequence[str] | None = None, jobs: int = 1,
                       cartesian: bool = False) -> None:
    """Fast path: ONLY the gene construction grid — skips the base edge_level/gene_level recompute
    + the edge conditioning/logit/function grids (which are already in canonical). ``cartesian`` runs
    the FULL function × reference × promiscuity-def construction cross (the heavy overnight grid)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    hs, edges, mirna, rna, cov_base, _conf, pressure_genes, target_cn = _load_grid_inputs(genes)
    print(f"[gene-grid-only] {len(edges):,} edges, {len(pressure_genes):,} target genes, "
          f"jobs={jobs}, cartesian={cartesian}")
    gg_gates = _build_ago_gates(rna, cov_base, hs)
    _run_gene_grid(out_dir, edges, mirna, rna, cov_base, target_cn, pressure_genes, gg_gates, jobs, cartesian)
    print(f"\n[done] gene-grid-only wrote {out_dir}")


def run(out_dir: Path = OUT_DIR, genes: Sequence[str] | None = None, conditioning: bool = True,
        logit_grid: bool = False, function_grid: bool = False, gene_grid: bool = False,
        jobs: int = 1) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    gene_list = list(genes or hs.universe)

    print(f"[load] edges / mirna / rna / clinical / target CN ({len(gene_list):,} universe genes)")
    edges = load_mirtar_edges(gene_list, resolve_arms=True)
    mirna = D.load_mirna_arms()
    rna = D.load_rna()
    clinical = D.load_clinical_strata().set_index("participant")
    conf_cols = [c for c in C.CONFOUNDER_NUMERIC if c in clinical.columns]
    cov_base = D.augment_tcga_batch(clinical[conf_cols].apply(pd.to_numeric, errors="coerce"))
    pressure_genes = sorted(set(edges["gene"]) & set(rna.index))
    target_cn = D.load_cnv_target_genes(pressure_genes)
    print(f"[load] {len(edges):,} edges, {mirna.shape[0]:,} arms, {len(pressure_genes):,} target genes, "
          f"confounders={conf_cols}+target_cn")

    edge_df = edge_level(edges, mirna, rna, cov_base, target_cn, pressure_genes)
    gene_df = gene_level(edges, mirna, rna, cov_base, target_cn, pressure_genes)

    edge_df.to_csv(out_dir / "edge_coupling_by_predictor.tsv.gz", sep="\t", index=False)
    gene_df.to_csv(out_dir / "gene_coupling_by_predictor.tsv.gz", sep="\t", index=False)

    summ = pd.concat(
        [
            _summary(edge_df, "edge", ["miRNA", "gene"], "abundance"),
            _summary(gene_df, "gene", ["gene"], "abundance_sum"),
        ],
        ignore_index=True,
    )
    summ.to_csv(out_dir / "predictor_comparison_summary.tsv", sep="\t", index=False)

    cond_df = pd.DataFrame()
    if conditioning:
        cond_df = edge_conditioning(edges, mirna, rna, cov_base, target_cn, pressure_genes)
        if not cond_df.empty:
            cond_df.to_csv(out_dir / "edge_conditioning_comparison.tsv.gz", sep="\t", index=False)
            cond_summ = _conditioning_summary(cond_df)
            cond_summ.to_csv(out_dir / "edge_conditioning_summary.tsv", sep="\t", index=False)
            cond_mech = _conditioning_mechanism(cond_df)
            cond_mech.to_csv(out_dir / "edge_conditioning_mechanism.tsv", sep="\t", index=False)
            print("\n===== EDGE PREDICTOR × CONDITIONING SWEEP — neg-sig counts =====")
            pivot = (cond_summ.pivot(index="conditioning", columns="predictor", values="n_neg_sig_cn")
                     .reindex(index=["none", "precision", "attribution"], columns=COND_PREDICTORS))
            with pd.option_context("display.max_columns", None, "display.width", 220):
                print(pivot.to_string())
                print()
                print(cond_summ.to_string(index=False))
                if not cond_mech.empty:
                    print(cond_mech.to_string(index=False))

    logit_df = pd.DataFrame()
    if logit_grid:
        print("\n[logit] building reference baselines (cohort/healthy/NAT) + promiscuity discount")
        baselines = _ref_baselines(mirna, edges)
        logit_preds = build_logit_predictors(edges, mirna, pressure_genes, baselines)
        full = edge_conditioning(edges, mirna, rna, cov_base, target_cn, pressure_genes,
                                 extra_predictors=logit_preds, include_pressure_modes=False)
        full["sample_set"] = "tumor_full"
        parts = [full]
        try:
            paired = edge_logit_paired(edges, rna, cov_base, target_cn, pressure_genes)
            if not paired.empty:
                paired["sample_set"] = "matched_nat"
                parts.append(paired)
        except Exception as e:  # noqa: BLE001
            print(f"[logit] paired-NAT block skipped: {e}")
        logit_df = pd.concat(parts, ignore_index=True)
        logit_df.to_csv(out_dir / "edge_logit_grid_comparison.tsv.gz", sep="\t", index=False)
        logit_summ = _logit_grid_summary(logit_df)
        logit_summ.to_csv(out_dir / "edge_logit_grid_summary.tsv", sep="\t", index=False)
        print("\n===== EDGE LOGIT GRID — neg-sig by (sample_set, predictor, conditioning) =====")
        with pd.option_context("display.max_columns", None, "display.width", 240,
                               "display.max_rows", None):
            print(logit_summ.to_string(index=False))

    func_df = pd.DataFrame()
    if function_grid:
        print("\n[func] building competition-function variants + AGO-gated dose (full/partial/purified)")
        baselines = _ref_baselines(mirna, edges)
        gates: Dict[str, pd.Series] = {}
        try:
            from mirna_hallmark.ago_gate import compute_ago_gate, gate_from_capacity
            from mirna_hallmark.analyses.cross_state.cross_state_coupling import _prolif_metagene
            from analysis.utils.common.loaders import residualize
            agt = compute_ago_gate(rna)
            full_gate = agt["ago_gate"]
            gates["abundance_agogated"] = full_gate
            gates["abundance_agogated_partial"] = 1.0 + 0.5 * (full_gate - 1.0)
            # purified: AGO capacity net of proliferation + CPE + HRD, then re-gate
            prolif = _prolif_metagene(rna, hs)
            covp = pd.concat([prolif.rename("prolif"), cov_base], axis=1)
            cap = agt["ago_capacity_z"]
            cap_pure = residualize(cap, covp.reindex(cap.index))
            if not cap_pure.empty:
                gates["abundance_agogated_purified"] = gate_from_capacity(cap_pure)
        except Exception as e:  # noqa: BLE001
            print(f"[func] AGO gate variants unavailable: {e}")
        # full (=HE) universe + a STRICTER evidence-pruned universe (prioritized regulators
        # for the field + share denominator + tested edges)
        cut = float(pd.to_numeric(edges["evidence_score"], errors="coerce").quantile(2 / 3))
        universes = {"all_he": edges, "strict_evid": edges[edges["evidence_score"] >= cut]}
        parts = []
        for uname, ued in universes.items():
            fp = build_function_predictors(ued, mirna, pressure_genes, baselines, gates)
            d = edge_conditioning(ued, mirna, rna, cov_base, target_cn, pressure_genes,
                                  extra_predictors=fp, include_pressure_modes=False)
            d["sample_set"] = uname
            parts.append(d)
        func_df = pd.concat(parts, ignore_index=True)
        func_df.to_csv(out_dir / "edge_function_grid_comparison.tsv.gz", sep="\t", index=False)
        func_summ = _logit_grid_summary(func_df)
        func_summ.to_csv(out_dir / "edge_function_grid_summary.tsv", sep="\t", index=False)
        print(f"\n===== EDGE FUNCTION/OCCUPANCY GRID — neg-sig (strict_evid cut evidence>={cut:.2f}) =====")
        with pd.option_context("display.max_columns", None, "display.width", 240,
                               "display.max_rows", None):
            print(func_summ.to_string(index=False))

    gene_grid_df = pd.DataFrame()
    if gene_grid:
        print("\n[gene-grid] GENE construction × prune × AGO-gating grid")
        gg_gates = _build_ago_gates(rna, cov_base, hs)
        gene_grid_df = _run_gene_grid(out_dir, edges, mirna, rna, cov_base, target_cn,
                                      pressure_genes, gg_gates, jobs)

    manifest = {
        "module": "mirna_hallmark.coupling_predictor_comparison",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": "does the pressure construction (share/z/level) change realized coupling vs raw abundance",
        "confounders": conf_cols + ["target_cn"],
        "target_norm": TARGET_NORM,
        "fdr_alpha": C.FDR_ALPHA,
        "min_n": MIN_N,
        "gated": False,
        "edge_modes": EDGE_MODES,
        "gene_modes": ["abundance_sum"] + GENE_ENGINE_MODES,
        "n_edges_scored_spine": int((edge_df["predictor"] == "softmax_z_logrpm").sum()),
        "n_genes_scored_spine": int((gene_df["predictor"] == "softmax_z_logrpm").sum()),
        "conditioning_sweep": {
            "enabled": bool(conditioning),
            "predictors": COND_PREDICTORS,
            "conditioning_modes": ["none", "precision", "attribution"],
            "p_others_def": "gene_total(g,s) - c(m,g,s); c = spine softmax_z_logrpm/evidence_mass field, FIXED across predictors",
            "covariate_policy": "Z0=CPE+HRD+target_cn; Z_mirna=Z0 minus target_cn (miRNA's own confounders)",
            "none": "resid X,Y on Z0 (reproduces published baseline)",
            "precision": "resid Y on Z0+P_others (all target covariates); resid X on Z_mirna (miRNA confounders only, co-regs NOT a focal confound)",
            "attribution": "resid X,Y on Z0+P_others (both-sides honest unique effect)",
            "n_edges": int(((cond_df["predictor"] == "abundance") & (cond_df["conditioning"] == "none")).sum()) if not cond_df.empty else 0,
        },
        "logit_grid": {
            "enabled": bool(logit_grid),
            "question": "does a softmax SHARE under any logit flavour (reference x measurement x promiscuity) beat raw abundance for edge coupling, incl. under precision conditioning",
            "references": ["rawx"] + LOGIT_REFS + ["pairednat"],
            "measures": ["dev", "z"],
            "promiscuity": ["p0", "pn"],
            "predictor": "bare gene-local softmax share (no post-multiply)",
            "sample_sets": sorted(logit_df["sample_set"].unique().tolist()) if not logit_df.empty else [],
            "n_predictors": int(logit_df["predictor"].nunique()) if not logit_df.empty else 0,
        },
        "function_grid": {
            "enabled": bool(function_grid),
            "question": "does a non-softmax competition function (linear/temperature/sparsemax/mass-action) "
                        "or AGO-gated abundance beat raw abundance for edge coupling",
            "functions": ["softmax", "tempT0.5", "tempT2", "sparsemax", "linearL1", "massaction", "abundance_agogated"],
            "note": "total miRNA is RPM-constant (CV~0.05) so global competition enters via the AGO-capacity gate, "
                    "not a total-load term; mass-action keeps a free-capacity term (not sum=1) so reduces to "
                    "saturating abundance at single-regulator genes",
            "n_edges": int(((func_df["predictor"] == "abundance") & (func_df["conditioning"] == "none")).sum()) if not func_df.empty else 0,
        },
        "gene_construction_grid": {
            "enabled": bool(gene_grid),
            "question": "GENE level: which construction x prune x AGO-gating best detects net repression "
                        "partial-rho(P_agg, y_g | CPE+HRD+target_cn) vs abundance_sum? (aggregation lemma: "
                        "construction is consequential at the aggregate, unlike the edge)",
            "constructions": ["abundance_sum"] + GENE_ANCHOR_MODES + sorted(
                gene_grid_df["predictor"].unique().tolist()) if not gene_grid_df.empty else ["abundance_sum"] + GENE_ANCHOR_MODES,
            "construction_space": "engine anchors + (function {softmax/linear/sparsemax/temperature/mass-action}) × "
                                  "(logit-reference {cohort/healthy/NAT/rawx}) × promiscuity, aggregated as Σ share·z·logrpm",
            "prunes": GENE_PRUNES,
            "gatings": ["ungated", "full", "partial", "purified"],
            "aggregate": "sum", "contrib_transform": "signed", "target_norm": TARGET_NORM,
            "n_configs": int(gene_grid_df.groupby(["predictor", "prune", "gating"]).ngroups) if not gene_grid_df.empty else 0,
        },
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("\n================ SUMMARY ================")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(summ.to_string(index=False))
    print(f"\n[done] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--max-genes", type=int, default=0, help="cap universe genes (debug)")
    ap.add_argument("--no-conditioning", action="store_false", dest="conditioning",
                    help="skip the E4 edge conditioning sweep (none/precision/attribution)")
    ap.add_argument("--logit-grid", action="store_true",
                    help="deep-dive: grid the softmax-share LOGIT (reference x measurement x "
                         "promiscuity, incl. healthy/NAT/paired-NAT) as predictors vs abundance "
                         "under each conditioning (opt-in; ~slow)")
    ap.add_argument("--function-grid", action="store_true",
                    help="deep-dive: grid the competition FUNCTION (softmax/temperature/sparsemax/"
                         "linear-L1/mass-action occupancy) + AGO-gated abundance vs abundance "
                         "under each conditioning (opt-in)")
    ap.add_argument("--gene-grid", action="store_true",
                    help="deep-dive: GENE-level construction × prune × AGO-gating grid — "
                         "partial-ρ(aggregate pressure, gene expr) vs abundance_sum (opt-in)")
    ap.add_argument("--gene-grid-only", action="store_true",
                    help="FAST PATH: run ONLY the gene construction grid (skip base edge/gene "
                         "recompute + edge conditioning/logit/function grids)")
    ap.add_argument("--jobs", type=int, default=1,
                    help="parallelism for the gene grid (over unit×prune jobs; joblib). -1 = all cores")
    ap.add_argument("--full-cartesian", action="store_true",
                    help="gene grid: FULL function × reference × promiscuity-def construction cross "
                         "(heavy; intended as a standalone overnight run)")
    args = ap.parse_args()
    genes = None
    if args.max_genes > 0:
        genes = list(HallmarkSets.load().universe)[: args.max_genes]
    if args.gene_grid_only or args.full_cartesian:
        run_gene_grid_only(out_dir=args.out_dir, genes=genes, jobs=args.jobs, cartesian=args.full_cartesian)
        return
    run(out_dir=args.out_dir, genes=genes, conditioning=args.conditioning,
        logit_grid=args.logit_grid, function_grid=args.function_grid, gene_grid=args.gene_grid,
        jobs=args.jobs)


if __name__ == "__main__":
    main()
