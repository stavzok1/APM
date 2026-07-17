"""Upstream miRTarBase evidence-score sensitivity on the L1 coupling spine.

Recomputes ``evidence_score`` from ``mirtar_interaction_summary.csv`` cross-count
columns (see ``evidence_scoring.py``), rebuilds M0 edges, and runs the same
Basal/key hallmark coupling grid as ``pressure_evidence_sensitivity``.

Pressure formula is fixed: ``softmax_z_logrpm`` + ``evidence_mass`` + ``sum``.

Run:
    .venv/bin/python3 -m mirna_hallmark.evidence_scoring_sensitivity
    .venv/bin/python3 -m mirna_hallmark.evidence_scoring_sensitivity --with-encori  # network
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.evidence_scoring import (
    SCORER_DESCRIPTIONS,
    SCORERS,
    apply_scorer,
    load_encori_mrna_boost,
)
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.mirna_arm_resolve import resolve_edges_mirna
from mirna_hallmark.pressure_engine import (
    TargetNormMode,
    _edge_weights,
    compute_gene_pressure,
    mirna_mass_denominators,
)
from mirna_hallmark.pressure_build import compute_gene_pressure_contributions
from mirna_hallmark.robustness_checks import (
    KEY_HALLMARKS,
    PAM50_SUBTYPES,
    _proliferation_proxies,
    hallmark_coupling_by_subtype,
)

OUT_DIR = C.OUTPUT_ROOT / "evidence_scoring_sensitivity"
DEFAULT_TARGET_NORM: TargetNormMode = "evidence_mass"
HUB_GENES = ("PTEN", "IRF1", "CDKN1A", "E2F1")


def _load_summary() -> pd.DataFrame:
    path = C.MIRTAR_HALLMARK_SUMMARY
    if not path.is_file():
        raise FileNotFoundError(f"Missing Hallmark miRTar summary: {path}")
    return pd.read_csv(path, low_memory=False)


def _build_edges(
    summary: pd.DataFrame,
    scorer_id: str,
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
    *,
    encori_boost: bool = False,
    fixed_pairs: pd.DataFrame | None = None,
) -> pd.DataFrame:
    scored = apply_scorer(summary, scorer_id)
    if encori_boost:
        boost = load_encori_mrna_boost(scored)
        scored["evidence_score"] = scored["evidence_score"] + boost.reindex(scored.index).fillna(0.0)
    sub = scored.loc[scored["gene"].isin(set(genes))].copy()
    sub = sub.loc[sub["miRNA"].astype(str).str.startswith("hsa-", na=False)]
    if fixed_pairs is not None and not fixed_pairs.empty:
        sub = fixed_pairs.merge(
            sub[["miRNA", "gene", "evidence_score"]],
            on=["miRNA", "gene"],
            how="left",
        )
        sub["evidence_score"] = pd.to_numeric(sub["evidence_score"], errors="coerce").fillna(0.0)
    else:
        sub = sub.loc[
            pd.to_numeric(sub["evidence_score"], errors="coerce").fillna(0) >= C.PRESSURE_MIN_EVIDENCE
        ]
    edges = sub[["miRNA", "gene", "evidence_score"]].drop_duplicates()
    resolved, _audit = resolve_edges_mirna(edges, mirna_df, weight_col="evidence_score")
    return resolved


def _coupling_summary(coupling: pd.DataFrame, scorer_id: str) -> pd.DataFrame:
    rows = []
    if coupling.empty or "subtype" not in coupling.columns:
        for sub in PAM50_SUBTYPES:
            rows.append({
                "scorer_id": scorer_id,
                "subtype": sub,
                "n_neg_sig": 0,
                "n_key_neg_sig": 0,
                "median_rho_prolif_cn_wsd_adj": float("nan"),
            })
        return pd.DataFrame(rows)
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub]
        if g.empty:
            continue
        neg = g.loc[(g["rho_prolif_cn_wsd_adj"] < 0) & (g["q_prolif_cn_wsd_adj"] < 0.10)]
        key = g.loc[g["key_hallmark"] == True]  # noqa: E712
        neg_key = neg.loc[neg["key_hallmark"] == True]  # noqa: E712
        rows.append({
            "scorer_id": scorer_id,
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_key_neg_sig": int(len(neg_key)),
            "median_rho_prolif_cn_wsd_adj": float(
                pd.to_numeric(g["rho_prolif_cn_wsd_adj"], errors="coerce").median()
            ),
        })
    return pd.DataFrame(rows)


def _edge_weight_table(edges: pd.DataFrame, target_norm: TargetNormMode) -> pd.DataFrame:
    if edges.empty:
        return edges
    mass = mirna_mass_denominators(edges, target_norm)
    w = _edge_weights(edges, mass, target_norm)
    out = edges[["miRNA", "gene", "evidence_score"]].copy()
    out["edge_w"] = w
    return out


def _hierarchy_metrics(edges: pd.DataFrame, target_norm: TargetNormMode) -> dict:
    ew = _edge_weight_table(edges, target_norm)
    if ew.empty:
        return {
            "n_edges": 0,
            "median_within_gene_max_min_ratio": float("nan"),
            "frac_genes_spread_lt_2x": float("nan"),
            "n_unique_evidence_scores": 0,
        }
    ratios = []
    for _gene, grp in ew.groupby("gene"):
        vals = pd.to_numeric(grp["edge_w"], errors="coerce").dropna()
        if len(vals) < 2:
            continue
        lo = float(vals.min())
        hi = float(vals.max())
        if lo > 0:
            ratios.append(hi / lo)
    ratio_s = pd.Series(ratios, dtype=float)
    return {
        "n_edges": int(len(ew)),
        "median_within_gene_max_min_ratio": float(ratio_s.median()) if len(ratio_s) else float("nan"),
        "frac_genes_spread_lt_2x": float((ratio_s < 2.0).mean()) if len(ratio_s) else float("nan"),
        "n_unique_evidence_scores": int(edges["evidence_score"].nunique()),
    }


def _hub_rank_rho(
    edges: pd.DataFrame,
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
    target_norm: TargetNormMode,
) -> pd.DataFrame:
    rows = []
    for gene in HUB_GENES:
        if gene not in genes:
            continue
        ge = edges.loc[edges["gene"] == gene]
        if ge.empty:
            continue
        ew = _edge_weight_table(ge, target_norm)
        contrib = compute_gene_pressure_contributions(
            [gene],
            edges=edges,
            expr_mode="softmax_z_logrpm",
            target_norm=target_norm,
            mirna=mirna_df,
            resolve_arms=False,
        )
        if contrib.empty:
            continue
        merged = ew.merge(
            contrib[["miRNA", "gene", "mean_abs_contribution"]],
            on=["miRNA", "gene"],
            how="inner",
        )
        if len(merged) < 3:
            rho = float("nan")
        else:
            rho, _p = spearmanr(merged["edge_w"], merged["mean_abs_contribution"])
        rows.append({
            "gene": gene,
            "n_arms": int(len(merged)),
            "spearman_edge_w_vs_abs_contrib": float(rho) if np.isfinite(rho) else float("nan"),
        })
    return pd.DataFrame(rows)


def run(
    *,
    out_dir: Path | None = None,
    scorers: Sequence[str] | None = None,
    with_encori: bool = False,
    fixed_edge_set: bool = True,
    target_norm: TargetNormMode = DEFAULT_TARGET_NORM,
) -> None:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)

    scorer_ids = list(scorers or SCORERS.keys())
    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted({g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    summary = _load_summary()

    baseline_edges = _build_edges(summary, "tiered_permissive", hs.universe, mirna)
    fixed_pairs = baseline_edges[["miRNA", "gene"]].drop_duplicates() if fixed_edge_set else None
    if fixed_edge_set:
        print(
            f"[evidence_scoring_sensitivity] fixed M0 edge set: "
            f"{len(fixed_pairs):,} (miRNA,gene) pairs from tiered_permissive"
        )

    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )["ago_gate"]

    spec_rows = []
    coupling_all = []
    summary_all = []
    hierarchy_all = []
    hub_all = []
    basal_key_all = []

    for sid in scorer_ids:
        desc = SCORER_DESCRIPTIONS.get(sid, sid)
        print(f"[evidence_scoring_sensitivity] {sid}: {desc}")
        spec_rows.append({
            "scorer_id": sid,
            "description": desc,
            "expr_mode": "softmax_z_logrpm",
            "target_norm": target_norm,
            "with_encori_boost": with_encori,
        })

        edges = _build_edges(
            summary, sid, hs.universe, mirna,
            encori_boost=with_encori,
            fixed_pairs=fixed_pairs,
        )
        hier = _hierarchy_metrics(edges, target_norm)
        hier["scorer_id"] = sid
        hierarchy_all.append(hier)
        hub_all.append(_hub_rank_rho(edges, hs.universe, mirna, target_norm).assign(scorer_id=sid))

        gp = compute_gene_pressure(
            edges,
            mirna,
            genes=list(hs.universe),
            expr_mode="softmax_z_logrpm",
            target_norm=target_norm,
            aggregate="sum",
        )
        shared = gp.columns.intersection(gate.index)
        hp = hallmark_pressure_matrix(gp[shared].mul(gate.reindex(shared), axis=1), hs)

        spec_dir = out_root / sid
        spec_dir.mkdir(parents=True, exist_ok=True)
        _edge_weight_table(edges, target_norm).to_csv(spec_dir / "edge_weights.tsv", sep="\t", index=False)

        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["scorer_id"] = sid
        coupling.to_csv(spec_dir / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t", index=False)
        coupling_all.append(coupling)
        summary_all.append(_coupling_summary(coupling, sid))

        basal = coupling.loc[
            (coupling["subtype"] == "Basal") & (coupling["hallmark_set"].isin(KEY_HALLMARKS))
        ].copy()
        basal_key_all.append(basal)

        b = summary_all[-1].loc[summary_all[-1]["subtype"] == "Basal"]
        if not b.empty:
            pten = hub_all[-1].loc[hub_all[-1]["gene"] == "PTEN", "spearman_edge_w_vs_abs_contrib"]
            pten_r = float(pten.iloc[0]) if len(pten) else float("nan")
            print(
                f"  Basal: neg_sig={int(b.iloc[0]['n_neg_sig'])}/50, "
                f"key={int(b.iloc[0]['n_key_neg_sig'])}/8, "
                f"median_rho={b.iloc[0]['median_rho_prolif_cn_wsd_adj']:.3f}, "
                f"edge_spread={hier['median_within_gene_max_min_ratio']:.2f}x, "
                f"PTEN_rank_rho={pten_r:.2f}"
            )

    pd.DataFrame(spec_rows).to_csv(out_root / "scorer_manifest.tsv", sep="\t", index=False)
    pd.concat(coupling_all, ignore_index=True).to_csv(
        out_root / "coupling_by_scorer_subtype.tsv", sep="\t", index=False,
    )
    summary = pd.concat(summary_all, ignore_index=True)
    summary.to_csv(out_root / "coupling_summary_by_scorer.tsv", sep="\t", index=False)
    pd.DataFrame(hierarchy_all).to_csv(out_root / "edge_hierarchy_by_scorer.tsv", sep="\t", index=False)
    pd.concat(hub_all, ignore_index=True).to_csv(out_root / "hub_rank_correlation.tsv", sep="\t", index=False)

    basal_wide = pd.concat(basal_key_all, ignore_index=True)
    pivot_rows = []
    for h in KEY_HALLMARKS:
        row: dict = {"hallmark_set": h}
        for sid in scorer_ids:
            sub = basal_wide.loc[
                (basal_wide["hallmark_set"] == h) & (basal_wide["scorer_id"] == sid)
            ]
            if sub.empty:
                row[f"{sid}_rho"] = float("nan")
                row[f"{sid}_neg_sig"] = False
            else:
                r = sub.iloc[0]
                row[f"{sid}_rho"] = float(r["rho_prolif_cn_wsd_adj"])
                row[f"{sid}_neg_sig"] = bool(
                    (r["rho_prolif_cn_wsd_adj"] < 0) and (r["q_prolif_cn_wsd_adj"] < 0.10)
                )
        pivot_rows.append(row)
    pd.DataFrame(pivot_rows).to_csv(out_root / "basal_key_pivot.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.evidence_scoring_sensitivity",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "scorers": scorer_ids,
        "with_encori_boost": with_encori,
        "fixed_edge_set": fixed_edge_set,
        "pressure": f"softmax_z_logrpm; target_norm={target_norm}; aggregate=sum",
        "benchmark": "Basal 42/50 neg-sig, 8/8 key (tiered_permissive baseline)",
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[evidence_scoring_sensitivity] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--scorers", nargs="+", default=None, choices=sorted(SCORERS.keys()))
    ap.add_argument(
        "--with-encori",
        action="store_true",
        help="Add ENCORI mRNA CLIP boost (requires network on first fetch)",
    )
    ap.add_argument(
        "--no-fixed-edge-set",
        action="store_true",
        help="Allow each scorer to change which (miRNA,gene) pairs pass min_evidence",
    )
    ap.add_argument(
        "--target-norm",
        default=DEFAULT_TARGET_NORM,
        choices=[
            "evidence_mass", "evidence_sum_log", "evidence_outer_log",
            "evidence_log_only", "degree", "degree_only", "none",
        ],
        help="Edge-weight normalization (evidence_sum_log = specificity option A)",
    )
    args = ap.parse_args()
    run(
        out_dir=args.out_dir,
        scorers=args.scorers,
        with_encori=args.with_encori,
        fixed_edge_set=not args.no_fixed_edge_set,
        target_norm=args.target_norm,
    )


if __name__ == "__main__":
    main()
