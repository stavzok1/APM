"""Hybrid miRNA pressure: miRTarBase + TargetScan edge weights.

Builds per-(gene, sample) pressure under configurable modes (M0–M11), replays
Hallmark coupling by PAM50, and writes a mode-comparison summary vs the
miRTarBase-only baseline (M0).

Default run modes (agent defaults): M7 (tiered fusion), M8 (sum tracks),
M11 (pure TargetScan, all pairs, per-gene cap).

Expression weighting uses ``pressure_engine`` (default: per-gene softmax over
targeting arms with logits = log2(RPM+1) − cohort_median). TS-only / orphan
edges additionally require cohort-median abundance ≥ ``PRESSURE_TS_ABUNDANCE_FLOOR``.

Outputs (``output/hybrid_pressure/``):
- ``{mode}/hybrid_edges.tsv.gz`` — (miRNA, gene, edge_weight, edge_tier)
- ``{mode}/gene_pressure_per_sample.tsv.gz`` — AGO-gated gene × sample
- ``{mode}/hallmark_pressure_per_sample.tsv.gz`` — Hallmark × sample
- ``{mode}/hallmark_coupling_by_pam50_prolif.tsv``
- ``{mode}/hub_gene_pressure_per_sample.tsv.gz`` — 8 hub targets only
- ``mode_comparison_summary.tsv`` — neg-sig counts vs M0
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import hybrid_pressure_kwargs, pressure_kwargs
from mirna_hallmark.pressure_engine import (
    AggregateMode,
    ExprMode,
    TargetNormMode,
    compute_gene_pressure,
    filter_edges_by_abundance,
)
from mirna_hallmark.robustness_checks import (
    HUB_ROUTES,
    KEY_HALLMARKS,
    PAM50_SUBTYPES,
    _load_targetscan_weights,
    _proliferation_proxies,
    gene_pressure_route_partial_corr,
    hallmark_coupling_by_subtype,
)
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna
from mirna_hallmark.eval.targetscan_orphan_coupling import _he_pairs, _mirtar_lookup

Pair = Tuple[str, str]

DEFAULT_MODES = ("M7", "M8", "M11")
MODE_ALIASES = {
    "mirtar_canonical": "M0",
    "ts_orphan": "M1",
    "hybrid_he_gap": "M2",
    "hybrid_he_orphan_downwt": "M3",
    "hybrid_he_orphan_mid": "M4",
    "he_ts_boost": "M5",
    "he_ts_boost_strong": "M6",
    "tiered_full": "M7",
    "sum_tracks_light": "M8",
    "sum_tracks_mid": "M9",
    "mirtar_any_ts_fallback": "M10",
    "ts_primary": "M11",
}


@dataclass(frozen=True)
class HybridPressureSpec:
    mode: str
    min_ts: float = 0.25
    top_k_program: int = 5
    top_k_hub: int = 15
    beta_ts_boost: float = 0.5
    gamma_orphan: float = 0.25
    lambda_orphan_track: float = 0.25
    min_evidence: int = 2


def _resolve_mode(mode: str) -> str:
    m = mode.strip().upper()
    if m in MODE_ALIASES.values():
        return m
    if mode.lower() in MODE_ALIASES:
        return MODE_ALIASES[mode.lower()]
    raise ValueError(f"Unknown hybrid pressure mode: {mode}")


def _spec_for_mode(mode: str) -> HybridPressureSpec:
    mode = _resolve_mode(mode)
    presets: Dict[str, HybridPressureSpec] = {
        "M0": HybridPressureSpec(mode="M0"),
        "M1": HybridPressureSpec(mode="M1"),
        "M2": HybridPressureSpec(mode="M2"),
        "M3": HybridPressureSpec(mode="M3", gamma_orphan=0.25),
        "M4": HybridPressureSpec(mode="M4", gamma_orphan=0.5),
        "M5": HybridPressureSpec(mode="M5", beta_ts_boost=0.5),
        "M6": HybridPressureSpec(mode="M6", beta_ts_boost=1.0),
        "M7": HybridPressureSpec(mode="M7", beta_ts_boost=0.5, gamma_orphan=0.25),
        "M8": HybridPressureSpec(mode="M8", lambda_orphan_track=0.25),
        "M9": HybridPressureSpec(mode="M9", lambda_orphan_track=0.5),
        "M10": HybridPressureSpec(mode="M10"),
        "M11": HybridPressureSpec(mode="M11"),
    }
    return presets[mode]


def _top_k(gene: str, hub_genes: Set[str], spec: HybridPressureSpec) -> int:
    return spec.top_k_hub if gene in hub_genes else spec.top_k_program


def _log1p(x: float) -> float:
    return float(np.log1p(max(float(x), 0.0)))


def _weight_he_ts_boost(evidence: float, ts_weight: float, beta: float) -> float:
    w_e = _log1p(evidence)
    if ts_weight <= 0:
        return w_e
    return w_e * (1.0 + beta * _log1p(ts_weight))


def _assign_edge_weight(
    *,
    is_he: bool,
    evidence: float,
    ts_weight: float,
    n_studies: int,
    spec: HybridPressureSpec,
    tier: str,
) -> float:
    """Return edge weight (passed to compute_pressure as evidence_score)."""
    mode = spec.mode
    if mode in ("M5", "M6") and is_he:
        return _weight_he_ts_boost(evidence, ts_weight, spec.beta_ts_boost)
    if mode == "M7":
        if is_he:
            if ts_weight > 0:
                return _weight_he_ts_boost(evidence, ts_weight, spec.beta_ts_boost)
            return _log1p(evidence)
        return spec.gamma_orphan * _log1p(ts_weight)
    if is_he:
        return _log1p(evidence)
    if mode in ("M3", "M4", "M8", "M9"):
        return spec.gamma_orphan * _log1p(ts_weight)
    if mode == "M10" and n_studies > 0 and evidence >= spec.min_evidence:
        return _log1p(evidence)
    if mode == "M11":
        return _log1p(ts_weight)
    # M2 gap / orphan fill-in, M1 orphan-only, M10 TS fallback
    return _log1p(ts_weight)


def _cap_orphan_per_gene(
    df: pd.DataFrame,
    hub_genes: Set[str],
    spec: HybridPressureSpec,
    *,
    min_ts: float,
) -> pd.DataFrame:
    if df.empty:
        return df
    sub = df.loc[~df["is_he"] & (df["ts_weight"] >= min_ts)].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(["gene", "ts_weight"], ascending=[True, False])
    rows = []
    for gene, grp in sub.groupby("gene", sort=False):
        k = _top_k(str(gene), hub_genes, spec)
        rows.append(grp.head(k))
    return pd.concat(rows, ignore_index=True) if rows else sub.iloc[:0]


def _cap_ts_per_gene(
    df: pd.DataFrame,
    hub_genes: Set[str],
    spec: HybridPressureSpec,
    *,
    min_ts: float,
) -> pd.DataFrame:
    """Top-k TargetScan arms per gene (all TS pairs, regardless of miRTarBase)."""
    if df.empty:
        return df
    sub = df.loc[df["ts_weight"] >= min_ts].copy()
    if sub.empty:
        return sub
    sub = sub.sort_values(["gene", "ts_weight"], ascending=[True, False])
    rows = []
    for gene, grp in sub.groupby("gene", sort=False):
        k = _top_k(str(gene), hub_genes, spec)
        rows.append(grp.head(k))
    return pd.concat(rows, ignore_index=True) if rows else sub.iloc[:0]


def build_hybrid_edges(
    genes: Sequence[str],
    *,
    spec: HybridPressureSpec,
    hallmark_edges: pd.DataFrame,
    mirtar_summary: Optional[Path] = None,
) -> pd.DataFrame:
    """Build (miRNA, gene, edge_weight, edge_tier) for one hybrid mode."""
    mode = spec.mode
    gene_list = sorted(set(genes))
    hub_genes = set(HUB_ROUTES.keys())
    he_pairs = _he_pairs(hallmark_edges)
    mirtar = _mirtar_lookup(Path(mirtar_summary or C.MIRTAR_HALLMARK_SUMMARY))
    ts = _load_targetscan_weights(gene_list)

    from mirna_hallmark.pressure_build import load_mirtar_edges

    mt_edges = load_mirtar_edges(gene_list)

    if mode == "M0":
        out = mt_edges.copy()
        out["edge_tier"] = "mirtar_evidence"
        if not ts.empty:
            out = out.merge(ts, on=["miRNA", "gene"], how="left")
        else:
            out["ts_weight"] = 0.0
        out["ts_weight"] = pd.to_numeric(out["ts_weight"], errors="coerce").fillna(0.0)
        return _resolve_hybrid_edge_table(
            out[["miRNA", "gene", "evidence_score", "ts_weight", "edge_tier"]]
        )

    if mode == "M1":
        from mirna_hallmark.eval.targetscan_orphan_coupling import build_orphan_edge_table

        orphan = build_orphan_edge_table(
            gene_list,
            he_pairs,
            min_ts=spec.min_ts,
            top_n_per_gene=None,
            mirtar=mirtar,
        )
        if orphan.empty:
            return pd.DataFrame(columns=["miRNA", "gene", "edge_weight", "edge_tier"])
        orphan = orphan.sort_values(["gene", "ts_weight"], ascending=[True, False])
        capped = []
        for gene, grp in orphan.groupby("gene", sort=False):
            capped.append(grp.head(_top_k(str(gene), hub_genes, spec)))
        orphan = pd.concat(capped, ignore_index=True) if capped else orphan.iloc[:0]
        out = orphan[["miRNA", "gene", "ts_weight"]].copy()
        out["edge_weight"] = out["ts_weight"].map(_log1p)
        out["edge_tier"] = "ts_orphan"
        return _resolve_hybrid_edge_table(out[["miRNA", "gene", "edge_weight", "edge_tier"]])

    he_df = hallmark_edges.loc[hallmark_edges["high_evidence"]].drop_duplicates(["miRNA", "gene"])
    he_df = he_df.loc[he_df["gene"].isin(gene_list), ["miRNA", "gene", "evidence_score"]].copy()
    if not ts.empty:
        he_df = he_df.merge(ts, on=["miRNA", "gene"], how="left")
    else:
        he_df["ts_weight"] = 0.0
    if not mirtar.empty:
        he_df = he_df.merge(mirtar[["miRNA", "gene", "n_studies"]], on=["miRNA", "gene"], how="left")
    else:
        he_df["n_studies"] = 0
    he_df["ts_weight"] = pd.to_numeric(he_df["ts_weight"], errors="coerce").fillna(0.0)
    he_df["n_studies"] = pd.to_numeric(he_df["n_studies"], errors="coerce").fillna(0).astype(int)
    he_df["evidence_score"] = pd.to_numeric(he_df["evidence_score"], errors="coerce").fillna(0.0)
    he_df["is_he"] = True

    ts_only = ts.loc[ts["gene"].isin(gene_list)].copy() if not ts.empty else pd.DataFrame(
        columns=["miRNA", "gene", "ts_weight"]
    )
    if not mirtar.empty and not ts_only.empty:
        ts_only = ts_only.merge(
            mirtar[["miRNA", "gene", "n_studies", "evidence_score"]],
            on=["miRNA", "gene"],
            how="left",
        )
    ts_only["n_studies"] = pd.to_numeric(ts_only.get("n_studies", 0), errors="coerce").fillna(0).astype(int)
    ts_only["evidence_score"] = pd.to_numeric(ts_only.get("evidence_score", 0), errors="coerce").fillna(0.0)
    ts_only["is_he"] = [pair in he_pairs for pair in zip(ts_only["miRNA"], ts_only["gene"])]

    if mode in ("M5", "M6"):
        selected = he_df.copy()
    elif mode == "M2":
        gap = ts_only.loc[~ts_only["is_he"] & (ts_only["ts_weight"] >= spec.min_ts) & (ts_only["n_studies"] == 0)]
        selected = pd.concat([he_df, gap], ignore_index=True)
    elif mode in ("M3", "M4", "M7", "M8", "M9"):
        orphan = _cap_orphan_per_gene(ts_only, hub_genes, spec, min_ts=spec.min_ts)
        selected = pd.concat([he_df, orphan], ignore_index=True)
    elif mode == "M10":
        ts_ok = ts_only.loc[(ts_only["ts_weight"] >= spec.min_ts) & (~ts_only["is_he"])].copy()
        mt_ok = mt_edges.copy()
        mt_ok["ts_weight"] = 0.0
        mt_ok["n_studies"] = 0
        mt_ok["is_he"] = False
        if not ts.empty:
            mt_ok = mt_ok.merge(ts[["miRNA", "gene", "ts_weight"]], on=["miRNA", "gene"], how="left", suffixes=("", "_ts"))
            if "ts_weight_ts" in mt_ok.columns:
                mt_ok["ts_weight"] = mt_ok["ts_weight_ts"].fillna(0.0)
        if not mirtar.empty:
            mt_ok = mt_ok.merge(mirtar[["miRNA", "gene", "n_studies"]], on=["miRNA", "gene"], how="left")
            mt_ok["n_studies"] = pd.to_numeric(mt_ok["n_studies"], errors="coerce").fillna(0).astype(int)
        selected = pd.concat([mt_ok, ts_ok], ignore_index=True).drop_duplicates(["miRNA", "gene"])
    elif mode == "M11":
        selected = _cap_ts_per_gene(ts_only, hub_genes, spec, min_ts=spec.min_ts)
    else:
        selected = he_df.copy()

    selected = selected.drop_duplicates(["miRNA", "gene"])

    rows = []
    for _, row in selected.iterrows():
        is_he = bool(row.get("is_he", False))
        ev = float(row.get("evidence_score", 0.0))
        tw = float(row.get("ts_weight", 0.0))
        ns = int(row.get("n_studies", 0))
        if is_he and tw >= spec.min_ts and mode in ("M5", "M6", "M7"):
            tier = "he_ts_boost" if tw > 0 else "he_only"
        elif mode == "M11":
            tier = "ts_pure"
        elif is_he:
            tier = "he"
        elif ns == 0:
            tier = "gap"
        else:
            tier = "orphan"
        w = _assign_edge_weight(
            is_he=is_he, evidence=ev, ts_weight=tw, n_studies=ns, spec=spec, tier=tier,
        )
        if w <= 0:
            continue
        rows.append({
            "miRNA": row["miRNA"], "gene": row["gene"],
            "edge_weight": w,
            "evidence_score": _log1p(ev) if is_he else 0.0,
            "ts_weight": tw,
            "edge_tier": tier,
        })
    return _resolve_hybrid_edge_table(pd.DataFrame(rows))


def _resolve_hybrid_edge_table(edges: pd.DataFrame) -> pd.DataFrame:
    """Map miRNA names onto GDC expression rows (isoform / alias)."""
    if edges.empty:
        return edges
    mirna = D.load_mirna_arms()
    weight_col = (
        "edge_weight" if "edge_weight" in edges.columns else "evidence_score"
    )
    resolved, _ = resolve_edges_mirna(edges, mirna, weight_col=weight_col)
    keep = [c for c in ("miRNA", "gene", "edge_weight", "evidence_score", "ts_weight", "edge_tier") if c in resolved.columns]
    return resolved[keep]


def _pressure_from_edges(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    genes: Sequence[str],
    gate: pd.Series,
    *,
    expr_mode: ExprMode,
    target_norm: TargetNormMode,
    aggregate: AggregateMode,
    abundance_floor: Optional[float],
    ts_abundance_floor: Optional[float],
    ts_only_mode: bool = False,
) -> pd.DataFrame:
    if edges.empty:
        return pd.DataFrame()
    e = edges.copy()
    if "evidence_score" not in e.columns:
        e = e.rename(columns={"edge_weight": "evidence_score"})
    keep = ["miRNA", "gene", "evidence_score"]
    if "edge_tier" in e.columns:
        keep.append("edge_tier")
    if "edge_weight" in e.columns:
        keep.append("edge_weight")
    if "ts_weight" in e.columns:
        keep.append("ts_weight")
    e = e[keep].copy()
    if "edge_tier" not in e.columns:
        e["edge_tier"] = "unknown"
    if abundance_floor is not None or ts_abundance_floor is not None:
        floor_he = abundance_floor if abundance_floor is not None else ts_abundance_floor
        floor_ts = ts_abundance_floor if ts_abundance_floor is not None else abundance_floor
        if ts_only_mode:
            e = filter_edges_by_abundance(e, mirna, floor_ts)
        else:
            he_mask = e["edge_tier"].isin({"mirtar_evidence", "he", "he_only", "he_ts_boost"})
            e_he = filter_edges_by_abundance(e.loc[he_mask], mirna, floor_he)
            e_ts = filter_edges_by_abundance(e.loc[~he_mask], mirna, floor_ts)
            e = pd.concat([e_he, e_ts], ignore_index=True)
    pass_cols = ["miRNA", "gene", "evidence_score"]
    if "edge_weight" in e.columns:
        pass_cols.append("edge_weight")
    if "ts_weight" in e.columns:
        pass_cols.append("ts_weight")
    e = e[pass_cols]
    gp = compute_gene_pressure(
        e, mirna, genes=sorted(set(genes)),
        expr_mode=expr_mode, target_norm=target_norm, aggregate=aggregate,
    )
    shared = gp.columns.intersection(gate.index)
    return gp[shared].mul(gate.reindex(shared), axis=1)


def compute_hybrid_gene_pressure(
    genes: Sequence[str],
    *,
    spec: HybridPressureSpec,
    mirna: pd.DataFrame,
    gate: pd.Series,
    hallmark_edges: pd.DataFrame,
    expr_mode: ExprMode = "z",
    target_norm: TargetNormMode = "none",
    aggregate: AggregateMode = "sum",
    abundance_floor: Optional[float] = None,
    ts_abundance_floor: Optional[float] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (gene_pressure_gated, hybrid_edges)."""
    pe_kwargs = {
        "expr_mode": expr_mode,
        "target_norm": target_norm,
        "aggregate": aggregate,
        "abundance_floor": abundance_floor,
        "ts_abundance_floor": ts_abundance_floor,
    }
    if spec.mode in ("M8", "M9"):
        he_spec = HybridPressureSpec(mode="M0")
        orph_spec = HybridPressureSpec(mode="M1", gamma_orphan=spec.gamma_orphan)
        he_e = build_hybrid_edges(genes, spec=he_spec, hallmark_edges=hallmark_edges)
        orph_e = build_hybrid_edges(genes, spec=orph_spec, hallmark_edges=hallmark_edges)
        if not orph_e.empty:
            orph_e = orph_e.copy()
            orph_e["edge_weight"] = spec.gamma_orphan * orph_e["edge_weight"]
        p_he = _pressure_from_edges(he_e, mirna, genes, gate, **pe_kwargs)
        p_or = _pressure_from_edges(
            orph_e, mirna, genes, gate, ts_only_mode=True, **pe_kwargs,
        )
        if p_he.empty and p_or.empty:
            return pd.DataFrame(), pd.concat([he_e, orph_e], ignore_index=True)
        if p_he.empty:
            gp = spec.lambda_orphan_track * p_or
        elif p_or.empty:
            gp = p_he
        else:
            shared = p_he.index.intersection(p_or.index)
            cols = p_he.columns.intersection(p_or.columns)
            gp = p_he.loc[shared, cols] + spec.lambda_orphan_track * p_or.loc[shared, cols]
        combined = pd.concat([
            he_e.assign(track="he"),
            orph_e.assign(track="orphan"),
        ], ignore_index=True)
        return gp, combined

    edges = build_hybrid_edges(genes, spec=spec, hallmark_edges=hallmark_edges)
    gp = _pressure_from_edges(
        edges, mirna, genes, gate,
        ts_only_mode=spec.mode in ("M1", "M11"),
        **pe_kwargs,
    )
    return gp, edges


def _count_negsig(coupling: pd.DataFrame, *, rho_col: str = "rho_prolif_cn_wsd_adj", q_col: str = "q_prolif_cn_wsd_adj") -> pd.DataFrame:
    rows = []
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub]
        if g.empty:
            continue
        neg = g.loc[(g[rho_col] < 0) & (g[q_col] < 0.10)]
        key = g.loc[g["key_hallmark"]]
        neg_key = neg.loc[neg["key_hallmark"]]
        rows.append({
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_hallmarks": int(len(g)),
            "n_key_neg_sig": int(len(neg_key)),
            "n_key_hallmarks": int(len(key)),
            "median_rho": float(pd.to_numeric(g[rho_col], errors="coerce").median()),
        })
    return pd.DataFrame(rows)


def mode_comparison_summary(
    results: Dict[str, pd.DataFrame],
    *,
    baseline: str = "M0",
) -> pd.DataFrame:
    """Summarize neg-sig counts per mode × subtype vs baseline."""
    base = _count_negsig(results[baseline]) if baseline in results else None
    rows = []
    for mode, coupling in results.items():
        summ = _count_negsig(coupling)
        for _, r in summ.iterrows():
            row = {"mode": mode, **r.to_dict()}
            if base is not None:
                b = base.loc[base["subtype"] == r["subtype"]]
                if not b.empty:
                    row["delta_n_neg_sig_vs_baseline"] = int(r["n_neg_sig"] - b.iloc[0]["n_neg_sig"])
                    row["delta_n_key_neg_sig_vs_baseline"] = int(r["n_key_neg_sig"] - b.iloc[0]["n_key_neg_sig"])
            rows.append(row)
    return pd.DataFrame(rows)


def run(
    *,
    modes: Optional[Sequence[str]] = None,
    out_dir: Optional[Path] = None,
    expr_mode: Optional[ExprMode] = None,
    target_norm: Optional[TargetNormMode] = None,
    aggregate: Optional[AggregateMode] = None,
    abundance_floor: Optional[float] = None,
    ts_abundance_floor: Optional[float] = None,
    skip_baseline: bool = False,
) -> None:
    out_root = Path(out_dir or C.HYBRID_PRESSURE_DIR)
    out_root.mkdir(parents=True, exist_ok=True)

    aggregate = aggregate or C.PRESSURE_AGGREGATE  # type: ignore[assignment]
    abundance_floor = C.PRESSURE_ABUNDANCE_FLOOR if abundance_floor is None else abundance_floor
    ts_abundance_floor = C.PRESSURE_TS_ABUNDANCE_FLOOR if ts_abundance_floor is None else ts_abundance_floor

    def _weighting_for_mode(mode: str) -> Tuple[ExprMode, TargetNormMode]:
        if expr_mode is not None or target_norm is not None:
            em = expr_mode or C.PRESSURE_EXPR_MODE  # type: ignore[assignment]
            tn = target_norm or C.PRESSURE_TARGET_NORM  # type: ignore[assignment]
            return em, tn
        if mode == "M0":
            kw = pressure_kwargs()
        else:
            kw = hybrid_pressure_kwargs()
        return kw["expr_mode"], kw["target_norm"]  # type: ignore[return-value]

    mode_list = [_resolve_mode(m) for m in (modes or DEFAULT_MODES)]
    if not skip_baseline and "M0" not in mode_list:
        mode_list = ["M0", *mode_list]

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    hallmark_edges = D.load_hallmark_edges()
    proxies = _proliferation_proxies(rna, hs)
    cnv_genes = sorted(set(HUB_ROUTES) | {g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(cnv_genes)

    gate_df = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )
    gate = gate_df["ago_gate"]

    coupling_by_mode: Dict[str, pd.DataFrame] = {}
    edge_counts: Dict[str, int] = {}
    weighting_by_mode: Dict[str, dict] = {}

    for mode in mode_list:
        spec = _spec_for_mode(mode)
        mode_em, mode_tn = _weighting_for_mode(mode)
        weighting_by_mode[mode] = {"expr_mode": mode_em, "target_norm": mode_tn}
        print(f"[hybrid_pressure] mode {mode} ({spec.mode}) expr={mode_em} target_norm={mode_tn} ...")
        mode_dir = out_root / mode.lower()
        mode_dir.mkdir(parents=True, exist_ok=True)

        gp, edges = compute_hybrid_gene_pressure(
            hs.universe,
            spec=spec,
            mirna=mirna,
            gate=gate,
            hallmark_edges=hallmark_edges,
            expr_mode=mode_em,
            target_norm=mode_tn,
            aggregate=aggregate,
            abundance_floor=abundance_floor,
            ts_abundance_floor=ts_abundance_floor,
        )
        edge_counts[mode] = len(edges)
        edges.to_csv(mode_dir / "hybrid_edges.tsv.gz", sep="\t", index=False, compression="gzip")
        if gp.empty:
            print(f"[hybrid_pressure]   WARNING: empty pressure for {mode}")
            continue
        gp.to_csv(mode_dir / "gene_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

        hp = hallmark_pressure_matrix(gp, hs)
        hp.to_csv(mode_dir / "hallmark_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

        hub_gp = gp.loc[[g for g in sorted(HUB_ROUTES) if g in gp.index]]
        hub_gp.to_csv(mode_dir / "hub_gene_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["mode"] = mode
        coupling.to_csv(mode_dir / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t", index=False)
        coupling_by_mode[mode] = coupling

        hub_routes = gene_pressure_route_partial_corr(rna, clinical, proxies, hub_gp, cnv)
        hub_routes["mode"] = mode
        hub_routes.to_csv(mode_dir / "hub_gene_pressure_route_partial_corr.tsv", sep="\t", index=False)

        print(f"[hybrid_pressure]   {len(edges):,} edges, {gp.shape[0]} genes, "
              f"Basal neg-sig={int(((coupling['subtype'] == 'Basal') & (coupling['rho_prolif_cn_wsd_adj'] < 0) & (coupling['q_prolif_cn_wsd_adj'] < 0.10)).sum())}/50")

    for mode_dir in sorted(out_root.iterdir()):
        if not mode_dir.is_dir():
            continue
        mode_key = mode_dir.name.upper()
        cpath = mode_dir / "hallmark_coupling_by_pam50_prolif.tsv"
        if cpath.exists() and mode_key not in coupling_by_mode:
            coupling_by_mode[mode_key] = pd.read_csv(cpath, sep="\t")
    cmp_tbl = mode_comparison_summary(coupling_by_mode, baseline="M0")
    cmp_tbl.to_csv(out_root / "mode_comparison_summary.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.hybrid_pressure",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "modes_run": mode_list,
        "default_modes": list(DEFAULT_MODES),
        "spine_weighting": {
            "expr_mode": C.PRESSURE_EXPR_MODE,
            "target_norm": C.PRESSURE_TARGET_NORM,
            "aggregate": aggregate,
            "modes": ["M0"],
        },
        "hybrid_weighting": {
            "expr_mode": C.PRESSURE_HYBRID_EXPR_MODE,
            "target_norm": C.PRESSURE_HYBRID_TARGET_NORM,
            "aggregate": aggregate,
            "modes": [m for m in mode_list if m != "M0"],
        },
        "weighting_by_mode": weighting_by_mode,
        "cli_expr_mode_override": expr_mode,
        "cli_target_norm_override": target_norm,
        "abundance_floor": abundance_floor,
        "ts_abundance_floor": ts_abundance_floor,
        "spec_defaults": {
            "min_ts": 0.25,
            "top_k_program": 5,
            "top_k_hub": 15,
            "beta_ts_boost": 0.5,
            "gamma_orphan": 0.25,
            "lambda_orphan_track": 0.25,
            "min_evidence": C.PRESSURE_MIN_EVIDENCE,
        },
        "mode_notes": {
            "M7": "tiered fusion: HE log1p(evidence) boosted by TS agreement; orphan TS downweighted",
            "M8": "sum tracks: P_mirtar + lambda * P_ts_orphan (capped, non-HE only)",
            "M11": "pure TargetScan: all TS pairs >= min_ts, top-k per gene, TS abundance floor",
        },
        "edge_counts": edge_counts,
        "outputs": [
            "{mode}/hybrid_edges.tsv.gz",
            "{mode}/gene_pressure_per_sample.tsv.gz",
            "{mode}/hallmark_pressure_per_sample.tsv.gz",
            "{mode}/hub_gene_pressure_per_sample.tsv.gz",
            "{mode}/hallmark_coupling_by_pam50_prolif.tsv",
            "{mode}/hub_gene_pressure_route_partial_corr.tsv",
            "mode_comparison_summary.tsv",
        ],
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[hybrid_pressure] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--modes",
        nargs="+",
        default=list(DEFAULT_MODES),
        help=f"Hybrid modes to run (default: {' '.join(DEFAULT_MODES)}; M0 baseline added automatically)",
    )
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument(
        "--expr-mode",
        choices=(
            "z", "softmax", "z_softmax", "softmax_z",
            "softmax_z_logrpm", "softmax_z_absratio", "softmax_z_blend",
        ),
        default=None,
        help="Override all modes (default: M0=PRESSURE_*, M7+=PRESSURE_HYBRID_*)",
    )
    ap.add_argument(
        "--abundance-floor",
        type=float,
        default=None,
        help=f"Cohort-median log2(RPM+1) floor for miRTarBase arms (default: {C.PRESSURE_ABUNDANCE_FLOOR})",
    )
    ap.add_argument(
        "--ts-abundance-floor",
        type=float,
        default=None,
        help=f"Stricter floor for TargetScan / orphan arms (default: {C.PRESSURE_TS_ABUNDANCE_FLOOR})",
    )
    ap.add_argument(
        "--target-norm",
        choices=("none", "degree", "evidence_mass", "ts_mass", "combined_mass"),
        default=None,
        help="Override all modes (default: M0=evidence_mass, M7+=combined_mass)",
    )
    ap.add_argument(
        "--aggregate",
        choices=("sum", "mean"),
        default=None,
        help=f"Aggregate over targeting arms (default: {C.PRESSURE_AGGREGATE})",
    )
    ap.add_argument(
        "--skip-baseline",
        action="store_true",
        help="Do not rerun M0; reuse existing m0/hallmark_coupling_by_pam50_prolif.tsv for comparison",
    )
    args = ap.parse_args()
    run(
        modes=args.modes,
        out_dir=args.out_dir,
        expr_mode=args.expr_mode,
        target_norm=args.target_norm,
        aggregate=args.aggregate,
        abundance_floor=args.abundance_floor,
        ts_abundance_floor=args.ts_abundance_floor,
        skip_baseline=args.skip_baseline,
    )


if __name__ == "__main__":
    main()
