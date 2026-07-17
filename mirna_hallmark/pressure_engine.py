"""Gene-level miRNA pressure with configurable expression and target-count weighting.

Core formula (per gene *g*, sample *s*):

    pressure(g, s) = AGGREGATE_m [ expr_mult(m, s) · edge_w(m, g) ]

``edge_w`` starts from ``log1p(evidence_score)`` (or hybrid ``edge_weight``) and may
be divided by a per-miRNA normalization denominator (degree / evidence mass / TS mass).
"""

from __future__ import annotations

from typing import Dict, List, Literal, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ExprMode = Literal[
    "z",
    "softmax",
    "softmax_logrpm",
    "softmax_absratio",
    "z_softmax",
    "softmax_z",
    "softmax_z_logrpm",
    "softmax_z_logrpm_spec",
    "softmax_z_logrpm_enc",
    "softmax_z_absratio",
    "softmax_z_blend",
    # healthy-anchored: replace the within-tumor z with deviation from the healthy
    # (GTEx) median, x(m,s) - median_healthy(m), in log2 units (no SD division, so no
    # variability gating). Requires a `healthy_baseline` per-arm Series; arms absent
    # in healthy get baseline 0 -> full expression counts as tumor-acquired.
    "devhealthy",
    "softmax_devhealthy",
    "softmax_devhealthy_logrpm",
]
TargetNormMode = Literal[
    "none",
    "unity",
    "degree",
    "degree_only",
    "evidence_mass",
    "evidence_sum_log",
    "evidence_outer_log",
    "evidence_log_only",
    "ts_mass",
    "combined_mass",
]
AggregateMode = Literal["sum", "mean"]
ContribTransform = Literal["signed", "pos", "abs"]

BLEND_ALPHA = 0.5  # softmax_z_blend: weight on gene-local share vs uniform


def cohort_median_expression(mirna_log2: pd.DataFrame) -> pd.Series:
    """Per-arm cohort-median log2(RPM+1) — shared abundance background."""
    return mirna_log2.median(axis=1)


def filter_edges_by_abundance(
    edges: pd.DataFrame,
    mirna_log2: pd.DataFrame,
    floor: Optional[float],
) -> pd.DataFrame:
    """Drop edges whose miRNA arm is below the cohort-median abundance floor."""
    if edges.empty or floor is None:
        return edges
    med = cohort_median_expression(mirna_log2)
    ok = set(med.index[pd.to_numeric(med, errors="coerce").fillna(-np.inf) >= float(floor)])
    return edges.loc[edges["miRNA"].isin(ok)].copy()


def mirna_target_degree(edges: pd.DataFrame) -> pd.Series:
    """Number of distinct genes each miRNA targets in the edge table."""
    if edges.empty:
        return pd.Series(dtype=float)
    return edges.groupby("miRNA")["gene"].nunique().astype(float)


def _log_edge_evidence(edges: pd.DataFrame) -> pd.Series:
    return np.log1p(pd.to_numeric(edges["evidence_score"], errors="coerce").fillna(0.0))


def _log_ts_weight(edges: pd.DataFrame) -> pd.Series:
    if "ts_weight" not in edges.columns:
        return pd.Series(0.0, index=edges.index)
    return np.log1p(pd.to_numeric(edges["ts_weight"], errors="coerce").fillna(0.0))


def mirna_mass_denominators(edges: pd.DataFrame, target_norm: TargetNormMode) -> pd.Series:
    """Per-miRNA normalization denominators for edge weights."""
    if edges.empty or target_norm in ("none", "unity", "evidence_log_only"):
        return pd.Series(dtype=float)
    log_ev = _log_edge_evidence(edges)
    log_ts = _log_ts_weight(edges)
    e = edges.copy()
    e["_log_ev"] = log_ev
    e["_log_ts"] = log_ts
    e["_combined"] = e["_log_ev"] + e["_log_ts"]

    if target_norm == "degree":
        return mirna_target_degree(edges).map(lambda d: np.log1p(max(float(d), 1.0)))

    if target_norm == "degree_only":
        return mirna_target_degree(edges).map(lambda d: np.log1p(max(float(d), 1.0)))

    if target_norm == "evidence_mass":
        mass = e.groupby("miRNA")["_log_ev"].sum()
        return np.log1p(mass.clip(lower=0.0))

    if target_norm == "evidence_sum_log":
        # Σ_g log1p(ev) — inner log only, no outer log1p on the sum
        return e.groupby("miRNA")["_log_ev"].sum().clip(lower=1e-9)

    if target_norm == "evidence_outer_log":
        # log1p(Σ_g ev) — outer log on raw study-count sum
        raw = pd.to_numeric(edges["evidence_score"], errors="coerce").fillna(0.0)
        e["_raw_ev"] = raw
        mass = e.groupby("miRNA")["_raw_ev"].sum()
        return np.log1p(mass.clip(lower=0.0))

    if target_norm == "evidence_log_only":
        return pd.Series(dtype=float)

    if target_norm == "ts_mass":
        mass = e.groupby("miRNA")["_log_ts"].sum()
        return np.log1p(mass.clip(lower=0.0))

    if target_norm == "combined_mass":
        mass = e.groupby("miRNA")["_combined"].sum()
        return np.log1p(mass.clip(lower=0.0))

    raise ValueError(f"Unknown target_norm: {target_norm}")


def expression_multiplier(
    mirna_log2: pd.DataFrame,
    mode: Literal["z", "softmax"],
) -> pd.DataFrame:
    """Return arm × sample matrix used as logits or multipliers."""
    x = mirna_log2.astype(float)
    if mode == "z":
        mu = x.mean(axis=1)
        sigma = x.std(axis=1).replace(0, np.nan)
        return x.sub(mu, axis=0).div(sigma, axis=0).fillna(0.0)
    if mode == "softmax":
        med = cohort_median_expression(x)
        return x.sub(med, axis=0).fillna(0.0)
    raise ValueError(f"Unknown base expr mode: {mode}")


def _softmax_rows(sub: pd.DataFrame) -> pd.DataFrame:
    """Gene-local softmax: sub is (n_arms × n_samples) logits."""
    shifted = sub.sub(sub.max(axis=0), axis=1)
    ex = np.exp(shifted)
    denom = ex.sum(axis=0).replace(0, np.nan)
    return ex.div(denom, axis=1).fillna(0.0)


def _abs_ratio_matrix(mirna_log2: pd.DataFrame) -> pd.DataFrame:
    """log2(RPM+1) / cohort_median, floored — universal abundance anchor."""
    med = cohort_median_expression(mirna_log2).replace(0, np.nan)
    ratio = mirna_log2.astype(float).div(med, axis=0)
    return ratio.clip(lower=0.25).fillna(0.25)


# --------------------------------------------------------------------------- #
# Construction primitives for the edge/gene question-grids (shared, canonical).
# These build per-(arm,gene) competition shares under a chosen normalizer and a
# chosen logit (reference × measurement × promiscuity); the *spine* path does not
# use them, but the coupling-predictor grids do, so they live here (not duplicated
# in a grid module). NOTE: data-dependent reference baselines (healthy/NAT medians)
# are computed by the caller and passed in via ``baselines``.
# --------------------------------------------------------------------------- #
def promisc_discount(edges: pd.DataFrame, index: pd.Index) -> pd.Series:
    """Per-arm promiscuity discount = log1p(#targets); subtracted from a logit so a
    broadly-targeting arm captures less share per gene (dilution over its target set)."""
    deg = mirna_target_degree(edges)
    return np.log1p(deg.reindex(index).fillna(0.0))


def logit_matrix(mirna: pd.DataFrame, baselines: Dict[str, pd.Series],
                 ref: str, measure: Optional[str], promisc_on: bool) -> pd.DataFrame:
    """arm × sample logit for a softmax/competition share. ``ref='rawx'`` is reference-free
    (x itself); otherwise center on the per-arm reference median (``baselines[f'{ref}_med']``)
    and for ``measure='z'`` scale by ``baselines['cohort_sd']``. Arms lacking a ref baseline
    fall back to the cohort median. Promiscuity discount subtracted when ``promisc_on``."""
    x = mirna.astype(float)
    if ref == "rawx":
        logit = x.copy()
    else:
        center = baselines[f"{ref}_med"].reindex(x.index)
        center = center.where(center.notna(), baselines["cohort_med"])
        logit = x.sub(center, axis=0)
        if measure == "z":
            logit = logit.div(baselines["cohort_sd"], axis=0)
    if promisc_on:
        logit = logit.sub(baselines["promisc"].reindex(x.index), axis=0)
    return logit.replace([np.inf, -np.inf], np.nan)


def sparsemax_rows(sub: pd.DataFrame) -> pd.DataFrame:
    """Column-wise sparsemax (Euclidean projection onto the simplex) over the arm axis —
    a sparse alternative to softmax (small shares set to exactly 0; only top-few 'bind')."""
    Z = sub.values.astype(float)
    n = Z.shape[0]
    Zs = -np.sort(-Z, axis=0)
    css = np.cumsum(Zs, axis=0)
    k = np.arange(1, n + 1).reshape(-1, 1)
    cond = (1.0 + k * Zs) > css
    kz = np.maximum((cond * k).max(axis=0), 1)
    tau = (np.take_along_axis(css, (kz - 1).reshape(1, -1), axis=0)[0] - 1.0) / kz
    return pd.DataFrame(np.maximum(Z - tau.reshape(1, -1), 0.0), index=sub.index, columns=sub.columns)


def competition_share_map(edges: pd.DataFrame, value_mat: pd.DataFrame, genes: Sequence[str],
                          func: str, T: float = 1.0) -> Dict[tuple, pd.Series]:
    """Per-(arm,gene) gene-local competitive weight under a chosen normalizer. Log-domain
    funcs (softmax/temp/sparsemax) take a signed logit; linear/massaction take ≥0 relative
    abundance. ``massaction`` = a_m/(1+Σ a) keeps a free-capacity term (NOT Σ=1), so it
    reduces to (saturating) abundance for single-regulator genes."""
    gene_set = set(genes)
    have = set(value_mat.index)
    samples = value_mat.columns
    out: Dict[tuple, pd.Series] = {}
    for gene, grp in edges.groupby("gene"):
        if gene not in gene_set:
            continue
        arms = [a for a in grp["miRNA"].unique() if a in have]
        if not arms:
            continue
        sub = value_mat.loc[arms, samples]
        sub = sub.loc[sub.notna().any(axis=1)]
        if sub.empty:
            continue
        if func in ("softmax", "temp"):
            sh = _softmax_rows(sub.fillna(-30.0) / T)
        elif func == "sparsemax":
            sh = sparsemax_rows(sub.fillna(-30.0))
        elif func == "linear":
            p = sub.fillna(0.0).clip(lower=0.0)
            sh = p.div(p.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
        elif func == "massaction":
            p = sub.fillna(0.0).clip(lower=0.0)
            sh = p.div(1.0 + p.sum(axis=0), axis=1)
        else:
            raise ValueError(f"unknown competition func: {func}")
        for arm in sh.index:
            out[(str(arm), str(gene))] = sh.loc[arm]
    return out


def gene_pressure_from_share(share_map: Dict[tuple, pd.Series], z_mat: pd.DataFrame,
                             logrpm: pd.DataFrame, genes: Sequence[str],
                             weights: Optional[Dict[tuple, float]] = None,
                             magnitude: str = "z_logrpm", aggregate: str = "sum",
                             contrib: str = "signed") -> pd.DataFrame:
    """Aggregate a per-(arm,gene) competition share into ``P_agg`` over the aggregation axes:
    ``magnitude`` ∈ {logrpm, z, z_logrpm} (per-arm contribution = share × that term; the magnitude
    keeps the sum non-degenerate, since a bare share sums to 1 over R(g)); ``contrib`` ∈
    {signed, pos, abs} (coherence of the contribution before aggregating); ``aggregate`` ∈
    {sum, mean}. ``weights`` = optional per-(arm,gene) co-movement aggregation weight."""
    gene_set = set(genes)
    acc: Dict[str, pd.Series] = {}
    cnt: Dict[str, int] = {}
    for (arm, gene), sh in share_map.items():
        if gene not in gene_set or arm not in z_mat.index:
            continue
        c = sh
        if magnitude in ("z", "z_logrpm"):
            c = c.mul(z_mat.loc[arm])
        if magnitude in ("logrpm", "z_logrpm"):
            c = c.mul(logrpm.loc[arm])
        if weights is not None:
            c = c * float(weights.get((arm, gene), 1.0))
        if contrib == "pos":
            c = c.clip(lower=0.0)
        elif contrib == "abs":
            c = c.abs()
        acc[gene] = c if gene not in acc else acc[gene].add(c, fill_value=0.0)
        cnt[gene] = cnt.get(gene, 0) + 1
    if not acc:
        return pd.DataFrame()
    df = pd.DataFrame(acc).T
    if aggregate == "mean":
        df = df.div(pd.Series(cnt).reindex(df.index).replace(0, np.nan), axis=0)
    return df


def _edge_weights(
    present: pd.DataFrame,
    mass_denom: pd.Series,
    target_norm: TargetNormMode,
) -> np.ndarray:
    log_ev = _log_edge_evidence(present).values
    log_ts = _log_ts_weight(present).values
    hybrid = (
        "edge_weight" in present.columns
        and "evidence_score" in present.columns
        and target_norm in ("none", "degree")
    )
    if hybrid:
        w = pd.to_numeric(present["edge_weight"], errors="coerce").fillna(0.0).values
    elif target_norm == "unity":
        w = np.ones(len(present), dtype=float)
    elif target_norm == "degree_only":
        w = np.ones(len(present), dtype=float)
    elif target_norm == "combined_mass":
        w = log_ev + log_ts
    else:
        w = log_ev.copy()
    if target_norm in (
        "degree",
        "degree_only",
        "evidence_mass",
        "evidence_sum_log",
        "evidence_outer_log",
        "ts_mass",
        "combined_mass",
    ):
        denom = present["miRNA"].map(mass_denom).fillna(1.0).values
        denom = np.maximum(denom, 1e-9)
        w = w / denom
    return w


def _softmax_gene_logits(ab_sub: pd.DataFrame, log_ev: np.ndarray) -> pd.DataFrame:
    """Gene-local softmax with per-arm log1p(evidence) added to abundance logits."""
    logits = ab_sub.add(log_ev.reshape(-1, 1), axis=0)
    return _softmax_rows(logits)


def _expr_contribution(
    arms: np.ndarray,
    samples: pd.Index,
    expr_mode: ExprMode,
    z_mat: pd.DataFrame,
    ab_mat: pd.DataFrame,
    logrpm_mat: pd.DataFrame,
    abs_ratio: pd.DataFrame,
    dev_mat: Optional[pd.DataFrame] = None,
    log_ev_arms: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    if expr_mode in ("devhealthy", "softmax_devhealthy", "softmax_devhealthy_logrpm"):
        if dev_mat is None:
            raise ValueError(f"expr_mode '{expr_mode}' requires healthy_baseline")
        dev = dev_mat.loc[arms, samples]
        if expr_mode == "devhealthy":
            return dev
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        if expr_mode == "softmax_devhealthy":
            return sm.mul(dev)
        abs_term = logrpm_mat.loc[arms, samples].clip(lower=0.0)
        return sm.mul(dev).mul(abs_term)
    if expr_mode == "z":
        return z_mat.loc[arms, samples]
    if expr_mode == "softmax":
        return _softmax_rows(ab_mat.loc[arms, samples])
    if expr_mode == "softmax_logrpm":
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        return sm.mul(logrpm_mat.loc[arms, samples].clip(lower=0.0))
    if expr_mode == "softmax_absratio":
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        return sm.mul(abs_ratio.loc[arms, samples])
    if expr_mode == "z_softmax":
        return _softmax_rows(z_mat.loc[arms, samples])
    if expr_mode == "softmax_z":
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        return sm.mul(z_mat.loc[arms, samples])
    if expr_mode == "softmax_z_logrpm":
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        abs_term = logrpm_mat.loc[arms, samples].clip(lower=0.0)
        return sm.mul(z_mat.loc[arms, samples]).mul(abs_term)
    if expr_mode == "softmax_z_logrpm_spec":
        if log_ev_arms is None:
            raise ValueError("softmax_z_logrpm_spec requires log_ev_arms")
        sm = _softmax_gene_logits(ab_mat.loc[arms, samples], log_ev_arms)
        abs_term = logrpm_mat.loc[arms, samples].clip(lower=0.0)
        return sm.mul(z_mat.loc[arms, samples]).mul(abs_term)
    if expr_mode == "softmax_z_logrpm_enc":
        if log_ev_arms is None:
            raise ValueError("softmax_z_logrpm_enc requires share_logit_bonus (log_ev_arms slot)")
        sm = _softmax_gene_logits(ab_mat.loc[arms, samples], log_ev_arms)
        abs_term = logrpm_mat.loc[arms, samples].clip(lower=0.0)
        return sm.mul(z_mat.loc[arms, samples]).mul(abs_term)
    if expr_mode == "softmax_z_absratio":
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        return sm.mul(z_mat.loc[arms, samples]).mul(abs_ratio.loc[arms, samples])
    if expr_mode == "softmax_z_blend":
        sm = _softmax_rows(ab_mat.loc[arms, samples])
        z_sub = z_mat.loc[arms, samples]
        blended = sm.mul(BLEND_ALPHA).add(1.0 - BLEND_ALPHA)
        return blended.mul(z_sub)
    raise ValueError(f"Unknown expr mode: {expr_mode}")


def _dev_health_matrix(mirna_log2: pd.DataFrame, healthy_baseline: Optional[pd.Series]) -> Optional[pd.DataFrame]:
    """x(m,s) - median_healthy(m); arms with no healthy baseline -> baseline 0
    (their full expression counts as tumor-acquired)."""
    if healthy_baseline is None:
        return None
    base = healthy_baseline.reindex(mirna_log2.index).fillna(0.0)
    return mirna_log2.astype(float).sub(base, axis=0)


def _apply_contrib_transform(contrib: pd.DataFrame, mode: ContribTransform) -> pd.DataFrame:
    if mode == "signed":
        return contrib
    if mode == "pos":
        return contrib.clip(lower=0.0)
    if mode == "abs":
        return contrib.abs()
    raise ValueError(f"Unknown contrib_transform: {mode}")


def compute_gene_pressure(
    edges: pd.DataFrame,
    mirna_log2: pd.DataFrame,
    *,
    genes: Sequence[str],
    expr_mode: ExprMode = "z",
    target_norm: TargetNormMode = "none",
    aggregate: AggregateMode = "sum",
    contrib_transform: ContribTransform = "signed",
    healthy_baseline: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Per-(gene, sample) pressure with configurable weighting."""
    if edges.empty:
        return pd.DataFrame()

    samples = mirna_log2.columns
    mass_denom = mirna_mass_denominators(edges, target_norm)
    z_mat = expression_multiplier(mirna_log2, "z")
    ab_mat = expression_multiplier(mirna_log2, "softmax")
    logrpm_mat = mirna_log2.astype(float)
    abs_ratio = _abs_ratio_matrix(mirna_log2)
    dev_mat = _dev_health_matrix(mirna_log2, healthy_baseline)

    out = pd.DataFrame(0.0, index=list(genes), columns=samples)
    counts = pd.Series(0, index=list(genes), dtype=int)

    for gene, grp in edges.groupby("gene"):
        present = grp.loc[grp["miRNA"].isin(z_mat.index)]
        if present.empty:
            continue
        arms = present["miRNA"].values
        weights = _edge_weights(present, mass_denom, target_norm)
        log_ev_arms = _log_edge_evidence(present).values
        if "share_logit_bonus" in present.columns:
            share_bonus = pd.to_numeric(present["share_logit_bonus"], errors="coerce").fillna(0.0).values
        else:
            share_bonus = np.zeros(len(arms), dtype=float)
        if expr_mode == "softmax_z_logrpm_spec":
            logit_arms: Optional[np.ndarray] = log_ev_arms
        elif expr_mode == "softmax_z_logrpm_enc":
            logit_arms = share_bonus
        else:
            logit_arms = None
        n_arms = len(arms)

        contrib = _apply_contrib_transform(
            _expr_contribution(
                arms, samples, expr_mode, z_mat, ab_mat, logrpm_mat, abs_ratio, dev_mat,
                log_ev_arms=logit_arms,
            ).mul(weights, axis=0),
            contrib_transform,
        )

        vals = contrib.sum(axis=0).values
        if aggregate == "mean" and n_arms > 0:
            vals = vals / float(n_arms)
        out.loc[gene] = vals
        counts[gene] = n_arms

    keep = counts.loc[counts > 0].index
    return out.loc[keep]


def compute_gene_pressure_contributions(
    edges: pd.DataFrame,
    mirna_log2: pd.DataFrame,
    *,
    genes: Sequence[str],
    expr_mode: ExprMode = "z",
    target_norm: TargetNormMode = "none",
    healthy_baseline: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Per-(gene, miRNA) realized contribution summaries under several share metrics.

    Four complementary share definitions are exported (no single one is canonical):

    - ``mean_abs_share`` — per-sample |contribution| / Σ|contribution|, averaged.
      Magnitude of involvement; ignores direction.
    - ``mean_pos_share`` — per-sample contribution⁺ / Σ contribution⁺, averaged.
      Share of the *repressive* (positive-pressure) signal only; arms that are
      net-activating in a sample contribute 0 to that sample's denominator.
    - ``global_abs_share`` — mean|contribution| / Σ_arms mean|contribution|.
      Gene-level magnitude decomposition; stable (single ratio, no per-sample
      denominator instability).
    - ``global_signed_share`` — mean(contribution) / Σ_arms mean(contribution).
      Gene-level *net* share; can be negative or exceed 1 when arms cancel, which
      is itself diagnostic of an incoherent pressure stack.
    """
    if edges.empty:
        return pd.DataFrame()

    samples = mirna_log2.columns
    mass_denom = mirna_mass_denominators(edges, target_norm)
    z_mat = expression_multiplier(mirna_log2, "z")
    ab_mat = expression_multiplier(mirna_log2, "softmax")
    logrpm_mat = mirna_log2.astype(float)
    abs_ratio = _abs_ratio_matrix(mirna_log2)
    dev_mat = _dev_health_matrix(mirna_log2, healthy_baseline)
    cohort_med = cohort_median_expression(logrpm_mat)

    rows = []
    for gene, grp in edges.groupby("gene"):
        if gene not in set(genes):
            continue
        present = grp.loc[grp["miRNA"].isin(z_mat.index)].copy()
        if present.empty:
            continue
        arms = present["miRNA"].values
        weights = _edge_weights(present, mass_denom, target_norm)
        log_ev_arms = _log_edge_evidence(present).values
        if "share_logit_bonus" in present.columns:
            share_bonus = pd.to_numeric(present["share_logit_bonus"], errors="coerce").fillna(0.0).values
        else:
            share_bonus = np.zeros(len(arms), dtype=float)
        if expr_mode == "softmax_z_logrpm_spec":
            logit_arms: Optional[np.ndarray] = log_ev_arms
        elif expr_mode == "softmax_z_logrpm_enc":
            logit_arms = share_bonus
        else:
            logit_arms = None
        contrib = _expr_contribution(
            arms, samples, expr_mode, z_mat, ab_mat, logrpm_mat, abs_ratio, dev_mat,
            log_ev_arms=logit_arms,
        ).mul(weights, axis=0)

        denom_abs = contrib.abs().sum(axis=0).replace(0, np.nan)
        abs_share = contrib.abs().div(denom_abs, axis=1)
        pos = contrib.clip(lower=0.0)
        denom_pos = pos.sum(axis=0).replace(0, np.nan)
        pos_share = pos.div(denom_pos, axis=1)

        mean_signed = contrib.mean(axis=1)
        mean_abs = contrib.abs().mean(axis=1)
        total_signed = float(mean_signed.sum())
        total_abs = float(mean_abs.sum())

        z_sub = z_mat.loc[arms, samples]
        log_sub = logrpm_mat.loc[arms, samples]
        med = cohort_med.reindex(arms)
        for idx, arm in enumerate(arms):
            row = present.iloc[idx]
            c = pd.to_numeric(contrib.iloc[idx], errors="coerce")
            z = pd.to_numeric(z_sub.iloc[idx], errors="coerce")
            x = pd.to_numeric(log_sub.iloc[idx], errors="coerce")
            ms = float(mean_signed.iloc[idx])
            ma = float(mean_abs.iloc[idx])
            rows.append(
                {
                    "gene": gene,
                    "miRNA": arm,
                    "expr_mode": expr_mode,
                    "target_norm": target_norm,
                    "edge_weight_norm": float(weights[idx]),
                    "evidence_score": float(pd.to_numeric(row.get("evidence_score"), errors="coerce")),
                    "mean_signed_contribution": ms,
                    "mean_abs_contribution": ma,
                    "mean_abs_share": float(abs_share.iloc[idx].mean(skipna=True)),
                    "mean_pos_share": float(pos_share.iloc[idx].mean(skipna=True)),
                    "global_abs_share": float(ma / total_abs) if total_abs > 0 else np.nan,
                    "global_signed_share": float(ms / total_signed) if total_signed != 0 else np.nan,
                    "median_log2rpm": float(med.iloc[idx]) if pd.notna(med.iloc[idx]) else np.nan,
                    "mean_log2rpm": float(x.mean()),
                    "sd_log2rpm": float(x.std()),
                    "iqr_log2rpm": float(x.quantile(0.75) - x.quantile(0.25)),
                    "mean_z": float(z.mean()),
                    "sd_z": float(z.std()),
                    "frac_abs_z_lt_0_25": float((z.abs() < 0.25).mean()),
                    "frac_abs_z_lt_0_50": float((z.abs() < 0.50).mean()),
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["gene", "expr_mode", "global_abs_share"], ascending=[True, True, False])


def compute_edge_pressure_map(
    edges: pd.DataFrame,
    mirna_log2: pd.DataFrame,
    *,
    genes: Sequence[str],
    expr_mode: ExprMode = "softmax_logrpm",
    target_norm: TargetNormMode = "evidence_mass",
) -> Dict[tuple[str, str], pd.Series]:
    """Per-(miRNA, gene) realized pressure ``c(m,g,s)`` as a sample-indexed Series.

    Same engine path as ``compute_gene_pressure_contributions`` (edge weight ×
    ``expr_mode`` multiplier) but retains per-sample values for coupling tests.
    """
    if edges.empty:
        return {}

    samples = mirna_log2.columns
    mass_denom = mirna_mass_denominators(edges, target_norm)
    z_mat = expression_multiplier(mirna_log2, "z")
    ab_mat = expression_multiplier(mirna_log2, "softmax")
    logrpm_mat = mirna_log2.astype(float)
    abs_ratio = _abs_ratio_matrix(mirna_log2)

    out: Dict[tuple[str, str], pd.Series] = {}
    gene_set = set(genes)
    for gene, grp in edges.groupby("gene"):
        if gene not in gene_set:
            continue
        present = grp.loc[grp["miRNA"].isin(z_mat.index)]
        if present.empty:
            continue
        arms = present["miRNA"].values
        weights = _edge_weights(present, mass_denom, target_norm)
        contrib = _expr_contribution(
            arms, samples, expr_mode, z_mat, ab_mat, logrpm_mat, abs_ratio, dev_mat=None,
        ).mul(weights, axis=0)
        for idx, arm in enumerate(arms):
            s = pd.to_numeric(contrib.iloc[idx], errors="coerce")
            s.index = samples
            out[(str(arm), str(gene))] = s
    return out
