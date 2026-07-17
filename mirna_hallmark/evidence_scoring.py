"""Alternative miRTarBase evidence scores for (miRNA, gene) edges.

The default ``evidence_score`` in ``build_edges`` counts studies per experiment class:
  3·reporter + 3·binding + 2·protein + 1·rna + 1·perturbation

That ignores Functional vs Weak MTI and reporter/protein vs binding-only quality, which
collapses many edges to the same score (e.g. 8.5). This module recomputes ``evidence_score``
from the rich cross-count columns already present in ``mirtar_interaction_summary.csv``.

ENCORI boost (optional): merge AGO CLIP / target-program support for protein-coding genes.
APM's ENCORI client lives in ``pipeline/lncRNA_interactions/encori.py`` (lncRNA-first);
mRNA queries use ``geneType=mRNA`` — see ``load_encori_mrna_boost``.

Scorers are registered in ``SCORERS``; consumed by ``evidence_scoring_sensitivity``.
"""

from __future__ import annotations

from typing import Callable, Dict, Sequence

import numpy as np
import pandas as pd
from pathlib import Path

from mirna_hallmark import config as C
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna

ScorerFn = Callable[[pd.DataFrame], pd.Series]

# Mirrors ``analysis.expression.mirna_target_integration._tiered_study_weights`` —
# this is what ``load_mirtar_edges`` applies when ``n_studies`` is present (spine default).
_TIERED_PAIRS = [
    ("n_reporter__functional_mti_studies", 3.0, 3.0),
    ("n_protein__functional_mti_studies", 3.0, 3.0),
    ("n_rna__functional_mti_studies", 2.5, 2.5),
    ("n_binding__functional_mti_studies", 2.0, 2.0),
    ("n_perturbation__functional_mti_studies", 2.0, 2.0),
    ("n_reporter__functional_mti_weak_studies", 0.0, 1.0),
    ("n_protein__functional_mti_weak_studies", 0.0, 1.0),
    ("n_rna__functional_mti_weak_studies", 0.0, 0.8),
    ("n_binding__functional_mti_weak_studies", 0.0, 0.6),
    ("n_perturbation__functional_mti_weak_studies", 0.0, 0.6),
]

_COL = {
    "rf": "n_reporter__functional_mti_studies",
    "pf": "n_protein__functional_mti_studies",
    "rw": "n_reporter__functional_mti_weak_studies",
    "pw": "n_protein__functional_mti_weak_studies",
    "bf": "n_binding__functional_mti_studies",
    "bw": "n_binding__functional_mti_weak_studies",
    "pert_f": "n_perturbation__functional_mti_studies",
    "pert_w": "n_perturbation__functional_mti_weak_studies",
    "rna_f": "n_rna__functional_mti_studies",
    "rna_w": "n_rna__functional_mti_weak_studies",
}

# Assay-directness weights for the confidence score: does the assay demonstrate
# *repression* (reporter/protein) or merely *proximity* (CLIP binding)?
_CONFIDENCE_WEIGHTS = {
    "reporter": 3.0,   # luciferase / seed-mutant: direct functional readout
    "protein": 2.5,    # western: functional, one step removed
    "rna": 1.5,        # qPCR / knockdown: functional but indirect
    "perturbation": 1.5,
    "binding": 0.5,    # CLIP: physical proximity, not repression
}
_CONFIDENCE_WEAK_DISCOUNT = 0.3  # miRTarBase "Functional MTI (Weak)" curator flag


def _num(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series(0.0, index=df.index)
    return pd.to_numeric(df[col], errors="coerce").fillna(0.0)


def _tiered_row(row: pd.Series, *, strong: bool) -> float:
    total = 0.0
    for col, w_strong, w_perm in _TIERED_PAIRS:
        v = pd.to_numeric(row.get(col), errors="coerce")
        if pd.notna(v) and v > 0:
            total += float(v) * (w_strong if strong else w_perm)
    return total


def score_tiered_permissive(df: pd.DataFrame) -> pd.Series:
    """Spine default via ``load_mirtar_edges`` (functional-MTI cross-count tiers)."""
    return df.apply(lambda r: _tiered_row(r, strong=False), axis=1)


def score_tiered_strong(df: pd.DataFrame) -> pd.Series:
    """Strong tier — weak functional MTI cross-counts contribute 0."""
    return df.apply(lambda r: _tiered_row(r, strong=True), axis=1)


def score_canonical(df: pd.DataFrame) -> pd.Series:
    """``build_edges`` CSV ``evidence_score`` (study counts by experiment class)."""
    return (
        3 * _num(df, "n_reporter_studies")
        + 3 * _num(df, "n_binding_studies")
        + 2 * _num(df, "n_protein_studies")
        + 1 * _num(df, "n_rna_studies")
        + 1 * _num(df, "n_perturbation_studies")
    )


def score_functional_validation(df: pd.DataFrame) -> pd.Series:
    """Weight Functional MTI cross-flags; reporter/protein > binding/rna."""
    return (
        6 * _num(df, _COL["rf"])
        + 6 * _num(df, _COL["pf"])
        + 2 * _num(df, _COL["rw"])
        + 2 * _num(df, _COL["pw"])
        + 2 * _num(df, _COL["pert_f"])
        + 1 * _num(df, _COL["rna_f"])
        + 0.5 * _num(df, _COL["bf"])
        + 0.25 * _num(df, _COL["bw"])
    )


def score_he_proxy(df: pd.DataFrame) -> pd.Series:
    """Proxy for high-evidence: functional reporter or protein cross-counts."""
    return (
        10 * _num(df, _COL["rf"])
        + 10 * _num(df, _COL["pf"])
        + 3 * _num(df, _COL["rw"])
        + 3 * _num(df, _COL["pw"])
        + _num(df, "n_studies")
    )


def score_study_diversity(df: pd.DataFrame) -> pd.Series:
    """Unique study count (``log1p`` applied later in ``edge_w``)."""
    return _num(df, "n_studies")


def score_confidence_logclass(df: pd.DataFrame) -> pd.Series:
    """Confidence score: assay-directness weights × diminishing returns × weak discount.

    Design intent (see ``docs/PRESSURE_FUTURE_OPTIONS.md`` D):
    - ``log1p`` is applied *within each evidence class* before summing, so a pair with
      10 reporter studies is not 10x a pair with 1 — this defuses publication bias, the
      main reason fame (not biology) currently picks within-gene winners.
    - Binding-only (CLIP) is weighted 0.5: physical proximity is not demonstrated repression.
    - "Functional MTI (Weak)" curator-flagged studies are discounted by 0.3.

    This is a *confidence* weight (is the edge real?), not a hierarchy driver — within-gene
    ranking is meant to come from specificity normalization + expression dynamics.
    """
    w = _CONFIDENCE_WEIGHTS
    delta = _CONFIDENCE_WEAK_DISCOUNT
    func = (
        w["reporter"] * np.log1p(_num(df, _COL["rf"]))
        + w["protein"] * np.log1p(_num(df, _COL["pf"]))
        + w["rna"] * np.log1p(_num(df, _COL["rna_f"]))
        + w["perturbation"] * np.log1p(_num(df, _COL["pert_f"]))
        + w["binding"] * np.log1p(_num(df, _COL["bf"]))
    )
    weak = (
        w["reporter"] * np.log1p(_num(df, _COL["rw"]))
        + w["protein"] * np.log1p(_num(df, _COL["pw"]))
        + w["rna"] * np.log1p(_num(df, _COL["rna_w"]))
        + w["perturbation"] * np.log1p(_num(df, _COL["pert_w"]))
        + w["binding"] * np.log1p(_num(df, _COL["bw"]))
    )
    return func + delta * weak


def score_binding_penalized(df: pd.DataFrame) -> pd.Series:
    """Canonical-like but down-weight binding-only mass."""
    binding_only = (
        _num(df, "n_binding_studies")
        - _num(df, _COL["rf"])
        - _num(df, _COL["pf"])
        - _num(df, _COL["rw"])
        - _num(df, _COL["pw"])
    ).clip(lower=0.0)
    return (
        3 * _num(df, "n_reporter_studies")
        + 2 * _num(df, "n_protein_studies")
        + 1 * _num(df, "n_rna_studies")
        + 1 * _num(df, "n_perturbation_studies")
        + 0.5 * binding_only
    )


def apply_scorer(summary: pd.DataFrame, scorer_id: str) -> pd.DataFrame:
    """Return summary copy with ``evidence_score`` replaced by the chosen scorer."""
    if scorer_id not in SCORERS:
        raise KeyError(f"Unknown scorer {scorer_id!r}; choose from {list(SCORERS)}")
    out = summary.copy()
    out["evidence_score"] = SCORERS[scorer_id](out).astype(float)
    return out


def compute_enc_depth(df: pd.DataFrame) -> pd.Series:
    """ENCORI CLIP depth score from collapsed ENCORI rows."""
    clip = _num(df, "clipExpNum")
    degra = _num(df, "degraExpNum")
    panc = _num(df, "pancancerNum")
    ts = _num(df, "TargetScan") > 0
    mir = _num(df, "miRanda") > 0
    pita = _num(df, "PITA") > 0
    consensus = (ts & mir & pita).astype(float) * 0.5
    return (
        2.0 * np.log1p(clip)
        + 1.0 * np.log1p(degra)
        + consensus
        + 0.25 * np.log1p(panc)
    )


def compute_enc_reliability(df: pd.DataFrame) -> pd.Series:
    """ENCORI reliability score — less count-heavy than ``compute_enc_depth``."""
    clip = _num(df, "clipExpNum")
    degra = _num(df, "degraExpNum")
    panc = _num(df, "pancancerNum")
    pred = (
        _num(df, "TargetScan") + _num(df, "miRanda") + _num(df, "PITA")
    ) / 3.0
    return (
        1.0 * np.log1p(clip)
        + 0.5 * np.log1p(degra)
        + 0.25 * pred
        + 0.1 * np.log1p(panc)
    )


EncoriScoreMode = str  # "depth" | "reliability" | "two_channel"

ENCORI_SCORERS = {
    "depth": compute_enc_depth,
    "reliability": compute_enc_reliability,
}


def compute_encori_edge_score(df: pd.DataFrame, mode: str = "depth") -> pd.Series:
    """Per-edge ENCORI score; ``two_channel`` returns uniform 1.0 (ranking separate)."""
    if mode == "two_channel":
        return pd.Series(1.0, index=df.index)
    fn = ENCORI_SCORERS.get(mode)
    if fn is None:
        raise KeyError(f"Unknown ENCORI score mode {mode!r}; choose from {list(ENCORI_SCORERS)} + two_channel")
    return fn(df)


def apply_encori_boost_to_edges(
    edges: pd.DataFrame,
    enc_lookup: pd.DataFrame,
    *,
    alpha: float,
) -> pd.DataFrame:
    """S1: evidence_score_eff = base + alpha * enc_depth on matching (miRNA, gene)."""
    if edges.empty or enc_lookup.empty or alpha == 0:
        return edges.copy()
    out = edges.copy()
    merged = out.merge(
        enc_lookup.rename(columns={"enc_depth": "_enc_depth"}),
        on=["miRNA", "gene"],
        how="left",
    )
    merged["_enc_depth"] = pd.to_numeric(merged["_enc_depth"], errors="coerce").fillna(0.0)
    merged["evidence_score"] = (
        pd.to_numeric(merged["evidence_score"], errors="coerce").fillna(0.0)
        + float(alpha) * merged["_enc_depth"]
    )
    return merged[["miRNA", "gene", "evidence_score"]]


def apply_mirtar_boost_to_edges(
    enc_edges: pd.DataFrame,
    mirtar_summary: pd.DataFrame,
    *,
    beta: float,
    scorer_id: str = "tiered_permissive",
) -> pd.DataFrame:
    """S2: evidence_score_eff = enc_depth + beta * log1p(miRTar tier)."""
    if enc_edges.empty or beta == 0:
        return enc_edges.copy()
    scored = apply_scorer(mirtar_summary, scorer_id)[["miRNA", "gene", "evidence_score"]]
    scored = scored.rename(columns={"evidence_score": "_mirtar_score"})
    out = enc_edges.merge(scored, on=["miRNA", "gene"], how="left")
    out["_mirtar_score"] = pd.to_numeric(out["_mirtar_score"], errors="coerce").fillna(0.0)
    base = pd.to_numeric(out["evidence_score"], errors="coerce").fillna(0.0)
    out["evidence_score"] = base + float(beta) * np.log1p(out["_mirtar_score"])
    return out[["miRNA", "gene", "evidence_score"]]


def build_m0_edges(
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
    *,
    scorer_id: str = "tiered_permissive",
    min_evidence: float | None = None,
) -> pd.DataFrame:
    """M0 spine edges: scorer + min_evidence + arm resolve."""
    summary = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, low_memory=False)
    summary = summary.loc[summary["gene"].isin(set(genes))]
    summary = summary.loc[summary["miRNA"].astype(str).str.startswith("hsa-", na=False)]
    scored = apply_scorer(summary, scorer_id)
    floor = C.PRESSURE_MIN_EVIDENCE if min_evidence is None else min_evidence
    scored = scored.loc[pd.to_numeric(scored["evidence_score"], errors="coerce").fillna(0) >= floor]
    edges = scored[["miRNA", "gene", "evidence_score"]].drop_duplicates()
    resolved, _ = resolve_edges_mirna(edges, mirna_df, weight_col="evidence_score")
    return resolved


def load_encori_mrna_table(parquet_path: Path | None = None) -> pd.DataFrame:
    """Load consolidated ENCORI mRNA targets parquet if present."""
    from mirna_hallmark import config as C

    path = Path(parquet_path or C.OUTPUT_ROOT / "encori" / "encori_mirna_mrna_targets.parquet")
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_parquet(path)


def load_encori_mrna_boost(
    summary: pd.DataFrame,
    *,
    parquet_path: Path | None = None,
    cache_dir=None,
    clip_exp_num: int = 1,
    alpha: float = 0.5,
    min_clip_exp: int = 1,
) -> pd.Series:
    """Per-edge ENCORI CLIP boost in [0, alpha] for (miRNA, gene) pairs with AGO support.

    Reads ``mirna_hallmark/output/encori/encori_mirna_mrna_targets.parquet`` when
    available (built by ``fetch_encori_mrna``). Falls back to live fetch only if
    the parquet is missing.
    """
    from pathlib import Path

    from pipeline.config import PATHS
    from pipeline.lncRNA_interactions.encori import EncoriQuery, fetch_encori_miRNA_targets_with_symbol_fallback

    boost = pd.Series(0.0, index=summary.index)
    enc = load_encori_mrna_table(parquet_path)
    if enc.empty:
        cache = Path(cache_dir or PATHS.working_dir / "external_cache")
        genes = sorted(summary["gene"].dropna().astype(str).unique())
        query = EncoriQuery(gene_type="mRNA", clip_exp_num=int(clip_exp_num))
        enc, _diag = fetch_encori_miRNA_targets_with_symbol_fallback(
            targets=genes,
            miRNA="all",
            query=query,
            cache_dir=cache,
            clip_exp_num=clip_exp_num,
        )
    if enc.empty or "geneName" not in enc.columns:
        return boost

    mcol = "miRNAname" if "miRNAname" in enc.columns else ("miRNA" if "miRNA" in enc.columns else None)
    if mcol is None:
        return boost

    sub = enc.copy()
    if min_clip_exp > 0 and "clipExpNum" in sub.columns:
        sub = sub.loc[pd.to_numeric(sub["clipExpNum"], errors="coerce").fillna(0) >= min_clip_exp]

    hits = (
        sub.assign(
            _m=sub[mcol].astype(str).str.strip(),
            _g=sub["geneName"].astype(str).str.strip(),
        )
        .loc[lambda d: d["_m"].ne("") & d["_g"].ne("")]
        .drop_duplicates(["_m", "_g"])
        .set_index(["_m", "_g"])
    )
    hit_set = set(hits.index)
    hit_set_lower = {(m.lower(), g) for m, g in hit_set}

    for idx, row in summary.iterrows():
        m, g = str(row["miRNA"]).strip(), str(row["gene"]).strip()
        if (m, g) in hit_set:
            boost.at[idx] = alpha
        elif (m.lower(), g) in hit_set_lower:
            boost.at[idx] = alpha * 0.5
    return boost


SCORERS: Dict[str, ScorerFn] = {
    "tiered_permissive": score_tiered_permissive,
    "tiered_strong": score_tiered_strong,
    "confidence_logclass": score_confidence_logclass,
    "csv_canonical": score_canonical,
    "functional_validation": score_functional_validation,
    "he_proxy": score_he_proxy,
    "study_diversity": score_study_diversity,
    "binding_penalized": score_binding_penalized,
}

SCORER_DESCRIPTIONS: Dict[str, str] = {
    "tiered_permissive": "Spine default — functional-MTI cross-count tiers (permissive)",
    "tiered_strong": "Functional-MTI tiers; weak cross-counts excluded",
    "confidence_logclass": "Assay-directness weights × log1p within class × weak discount (binding=0.5)",
    "csv_canonical": "build_edges CSV: 3·rep + 3·bind + 2·prot + 1·rna + 1·pert",
    "functional_validation": "Functional MTI cross-flags; rep/prot weighted >> binding",
    "he_proxy": "Heavy weight on functional reporter/protein cross-counts",
    "study_diversity": "n_studies raw count (log1p in edge_w)",
    "binding_penalized": "Canonical minus excess binding-only studies",
}
