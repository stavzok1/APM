"""Shared gene-pressure builder for the mirna_hallmark spine.

Default weighting (``config.PRESSURE_*``):
  expr_mode=softmax_z_logrpm, target_norm=evidence_mass, aggregate=sum

Hybrid modes (M7/M8/M11) use ``config.PRESSURE_HYBRID_*`` (default M7:
  expr_mode=softmax_z, target_norm=combined_mass).
"""

from __future__ import annotations

from typing import Optional, Sequence

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna, write_arm_resolve_audit
from mirna_hallmark.pressure_engine import (
    AggregateMode,
    ExprMode,
    TargetNormMode,
    compute_gene_pressure as _engine_pressure,
    compute_gene_pressure_contributions as _engine_pressure_contributions,
    filter_edges_by_abundance,
)


def hybrid_pressure_kwargs(
    *,
    expr_mode: Optional[ExprMode] = None,
    target_norm: Optional[TargetNormMode] = None,
    aggregate: Optional[AggregateMode] = None,
) -> dict:
    return {
        "expr_mode": expr_mode or C.PRESSURE_HYBRID_EXPR_MODE,  # type: ignore[return-value]
        "target_norm": target_norm or C.PRESSURE_HYBRID_TARGET_NORM,  # type: ignore[return-value]
        "aggregate": aggregate or C.PRESSURE_AGGREGATE,  # type: ignore[return-value]
    }


def pressure_kwargs(
    *,
    expr_mode: Optional[ExprMode] = None,
    target_norm: Optional[TargetNormMode] = None,
    aggregate: Optional[AggregateMode] = None,
) -> dict:
    return {
        "expr_mode": expr_mode or C.PRESSURE_EXPR_MODE,  # type: ignore[return-value]
        "target_norm": target_norm or C.PRESSURE_TARGET_NORM,  # type: ignore[return-value]
        "aggregate": aggregate or C.PRESSURE_AGGREGATE,  # type: ignore[return-value]
    }


def load_mirtar_edges(genes: Sequence[str], *, resolve_arms: bool = True) -> pd.DataFrame:
    """M0 Hallmark edges with unified S1 evidence (scorer + optional ENCORI α boost)."""
    from mirna_hallmark.encori_edges import enc_depth_lookup
    from mirna_hallmark.evidence_scoring import apply_encori_boost_to_edges, build_m0_edges

    mirna_df = D.load_mirna_arms() if resolve_arms else None
    scorer = getattr(C, "PRESSURE_EVIDENCE_SCORER", "tiered_permissive")
    alpha = float(getattr(C, "PRESSURE_ENCORI_ALPHA", 0.0))
    if scorer in ("tiered_permissive", "tiered_strong"):
        from analysis.expression.mirna_target_integration import load_mirtar_edges as _legacy

        e = _legacy(
            C.MIRTAR_HALLMARK_SUMMARY,
            genes=list(genes),
            min_evidence=C.PRESSURE_MIN_EVIDENCE,
            weight_mode="permissive" if scorer == "tiered_permissive" else "strong",
        )
        if resolve_arms and not e.empty and mirna_df is not None:
            resolved, _audit = resolve_edges_mirna(e, mirna_df)
            return resolved
        return e

    e = build_m0_edges(
        list(genes),
        D.load_mirna_arms(),
        scorer_id=scorer,
        min_evidence=C.PRESSURE_MIN_EVIDENCE,
    )
    if alpha > 0 and not e.empty:
        enc_lookup = enc_depth_lookup(genes=list(genes))
        e = apply_encori_boost_to_edges(e, enc_lookup, alpha=alpha)
    return e


def resolve_pressure_edges(
    edges: pd.DataFrame,
    mirna: Optional[pd.DataFrame] = None,
    *,
    write_audit: bool = False,
) -> pd.DataFrame:
    """Map edge miRNA names onto the expression matrix (isoform aliases)."""
    if edges.empty:
        return edges
    mirna_df = mirna if mirna is not None else D.load_mirna_arms()
    weight_col = "evidence_score" if "evidence_score" in edges.columns else "edge_weight"
    resolved, audit = resolve_edges_mirna(edges, mirna_df, weight_col=weight_col)
    if write_audit and not audit.empty:
        write_arm_resolve_audit(audit)
    return resolved


def compute_gene_pressure(
    genes: Sequence[str],
    edges: Optional[pd.DataFrame] = None,
    *,
    expr_mode: Optional[ExprMode] = None,
    target_norm: Optional[TargetNormMode] = None,
    aggregate: Optional[AggregateMode] = None,
    abundance_floor: Optional[float] = None,
    mirna: Optional[pd.DataFrame] = None,
    resolve_arms: bool = True,
    write_arm_audit: bool = False,
    healthy_baseline: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Per-(gene, sample) pressure using subproject defaults unless overridden."""
    mirna_df = mirna if mirna is not None else D.load_mirna_arms()
    e = edges if edges is not None else load_mirtar_edges(genes, resolve_arms=False)
    keep_cols = ["miRNA", "gene", "evidence_score", "share_logit_bonus"]
    e = e[[c for c in keep_cols if c in e.columns]].copy()
    if "evidence_score" not in e.columns and "edge_weight" in e.columns:
        e = e.rename(columns={"edge_weight": "evidence_score"})
    if resolve_arms:
        e = resolve_pressure_edges(e, mirna_df, write_audit=write_arm_audit)
    floor = C.PRESSURE_ABUNDANCE_FLOOR if abundance_floor is None else abundance_floor
    if floor is not None:
        e = filter_edges_by_abundance(e, mirna_df, floor)
    return _engine_pressure(
        e,
        mirna_df,
        genes=list(genes),
        healthy_baseline=healthy_baseline,
        **pressure_kwargs(expr_mode=expr_mode, target_norm=target_norm, aggregate=aggregate),
    )


def compute_gene_pressure_contributions(
    genes: Sequence[str],
    edges: Optional[pd.DataFrame] = None,
    *,
    expr_mode: Optional[ExprMode] = None,
    target_norm: Optional[TargetNormMode] = None,
    abundance_floor: Optional[float] = None,
    mirna: Optional[pd.DataFrame] = None,
    resolve_arms: bool = True,
    write_arm_audit: bool = False,
    healthy_baseline: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """Per-(gene, miRNA) realized pressure share using subproject defaults."""
    mirna_df = mirna if mirna is not None else D.load_mirna_arms()
    e = edges if edges is not None else load_mirtar_edges(genes, resolve_arms=False)
    keep_cols = ["miRNA", "gene", "evidence_score", "share_logit_bonus"]
    e = e[[c for c in keep_cols if c in e.columns]].copy()
    if "evidence_score" not in e.columns and "edge_weight" in e.columns:
        e = e.rename(columns={"edge_weight": "evidence_score"})
    if resolve_arms:
        e = resolve_pressure_edges(e, mirna_df, write_audit=write_arm_audit)
    floor = C.PRESSURE_ABUNDANCE_FLOOR if abundance_floor is None else abundance_floor
    if floor is not None:
        e = filter_edges_by_abundance(e, mirna_df, floor)
    return _engine_pressure_contributions(
        e,
        mirna_df,
        genes=list(genes),
        expr_mode=expr_mode or C.PRESSURE_EXPR_MODE,  # type: ignore[arg-type]
        target_norm=target_norm or C.PRESSURE_TARGET_NORM,  # type: ignore[arg-type]
        healthy_baseline=healthy_baseline,
    )


def method_blurb() -> str:
    kw = pressure_kwargs()
    scorer = getattr(C, "PRESSURE_EVIDENCE_SCORER", C.PRESSURE_WEIGHT_MODE)
    alpha = float(getattr(C, "PRESSURE_ENCORI_ALPHA", 0.0))
    enc = f"; encori_alpha={alpha}" if alpha > 0 else ""
    return (
        f"softmax_z abundance-share × cohort-z; target_norm={kw['target_norm']}; "
        f"aggregate={kw['aggregate']}; min_evidence={C.PRESSURE_MIN_EVIDENCE}; "
        f"evidence_scorer={scorer}{enc}"
    )
