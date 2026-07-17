"""ENCORI-only S2 evidence-formula sensitivity (beta=0, no miRTar boost).

Compares three ENCORI edge scores on the same gates/cap:
  - depth        : current enc_depth (CLIP-heavy)
  - reliability  : softer count weighting + prediction mean
  - two_channel  : rank by enc_depth, uniform evidence_score=1 for pressure

Run:
    .venv/bin/python3 -m mirna_hallmark.encori_evidence_sensitivity
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
from mirna_hallmark.encori_edges import EncoriM0Spec, build_encori_m0_edges
from mirna_hallmark.evidence_scoring import build_m0_edges
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_engine import (
    TargetNormMode,
    _edge_weights,
    compute_gene_pressure,
    mirna_mass_denominators,
)
from mirna_hallmark.pressure_build import compute_gene_pressure_contributions
from mirna_hallmark.robustness_checks import (
    HUB_ROUTES,
    KEY_HALLMARKS,
    PAM50_SUBTYPES,
    _proliferation_proxies,
    hallmark_coupling_by_subtype,
)

OUT_DIR = C.OUTPUT_ROOT / "encori_evidence_sensitivity"
DEFAULT_TARGET_NORM: TargetNormMode = "evidence_mass"
HUB_GENES = tuple(HUB_ROUTES.keys())
SCORE_MODES = ("depth", "reliability", "two_channel")


def _coupling_summary(coupling: pd.DataFrame, run_id: str, mode: str) -> pd.DataFrame:
    rows = []
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub] if not coupling.empty else coupling
        if g.empty:
            rows.append({
                "run_id": run_id, "score_mode": mode, "subtype": sub,
                "n_neg_sig": 0, "n_key_neg_sig": 0,
            })
            continue
        neg = g.loc[(g["rho_prolif_cn_wsd_adj"] < 0) & (g["q_prolif_cn_wsd_adj"] < 0.10)]
        neg_key = neg.loc[neg["key_hallmark"] == True]  # noqa: E712
        rows.append({
            "run_id": run_id,
            "score_mode": mode,
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_key_neg_sig": int(len(neg_key)),
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
            [gene], edges=edges, expr_mode="softmax_z_logrpm",
            target_norm=target_norm, mirna=mirna_df, resolve_arms=False,
        )
        if contrib.empty:
            continue
        merged = ew.merge(
            contrib[["miRNA", "gene", "mean_abs_contribution"]],
            on=["miRNA", "gene"], how="inner",
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
    score_modes: Sequence[str] | None = None,
    clip_min: int = 2,
    top_k: int = 15,
    prolif_key: str = "e2f_g2m",
) -> None:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)
    modes = list(score_modes or SCORE_MODES)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted({g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)

    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )["ago_gate"]

    s0_edges = build_m0_edges(hs.universe, mirna)
    print(f"[encori_evidence] S0 reference: {len(s0_edges):,} M0 pairs")

    spec_rows: List[dict] = []
    coupling_all = []
    summary_all = []
    hub_all = []

    for mode in modes:
        rid = f"S2_encori_only__score={mode}__clip={clip_min}__K={top_k}"
        print(f"[encori_evidence] {rid}")
        spec = EncoriM0Spec(
            clip_min=clip_min, top_k_program=top_k, top_k_hub=50, score_mode=mode,
        )
        spec_rows.append({
            "run_id": rid, "score_mode": mode, "clip_min": clip_min,
            "top_k": top_k, "beta": 0, "prolif_key": prolif_key,
        })

        edges = build_encori_m0_edges(hs.universe, mirna, spec=spec)
        print(f"  edges={len(edges):,}")

        hub_df = _hub_rank_rho(edges, hs.universe, mirna, DEFAULT_TARGET_NORM).assign(
            run_id=rid, score_mode=mode,
        )
        hub_all.append(hub_df)

        gp = compute_gene_pressure(
            edges, mirna, genes=list(hs.universe),
            expr_mode="softmax_z_logrpm", target_norm=DEFAULT_TARGET_NORM, aggregate="sum",
        )
        shared = gp.columns.intersection(gate.index)
        hp = hallmark_pressure_matrix(gp[shared].mul(gate.reindex(shared), axis=1), hs)

        coupling = hallmark_coupling_by_subtype(
            rna, clinical, proxies, hs, hp, cnv, prolif_key=prolif_key,
        )
        coupling["run_id"] = rid
        coupling["score_mode"] = mode
        coupling_all.append(coupling)
        summary_all.append(_coupling_summary(coupling, rid, mode))

        b = summary_all[-1].loc[summary_all[-1]["subtype"] == "Basal"]
        irf1 = hub_df.loc[hub_df["gene"] == "IRF1", "spearman_edge_w_vs_abs_contrib"]
        irf1_r = float(irf1.iloc[0]) if len(irf1) else float("nan")
        if not b.empty:
            print(
                f"  Basal: neg_sig={int(b.iloc[0]['n_neg_sig'])}/50, "
                f"key={int(b.iloc[0]['n_key_neg_sig'])}/8, IRF1_rank_rho={irf1_r:.2f}"
            )

    pd.DataFrame(spec_rows).to_csv(out_root / "spec_manifest.tsv", sep="\t", index=False)
    pd.concat(summary_all, ignore_index=True).to_csv(
        out_root / "coupling_summary.tsv", sep="\t", index=False,
    )
    pd.concat(coupling_all, ignore_index=True).to_csv(
        out_root / "coupling_detail.tsv", sep="\t", index=False,
    )
    pd.concat(hub_all, ignore_index=True).to_csv(
        out_root / "hub_rank_correlation.tsv", sep="\t", index=False,
    )

    manifest = {
        "module": "mirna_hallmark.encori_evidence_sensitivity",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "score_modes": modes,
        "clip_min": clip_min,
        "top_k": top_k,
        "beta": 0,
        "prolif_covariate": prolif_key,
        "note": "Registered coupling benchmark uses rho_prolif_cn_wsd_adj with e2f_g2m by default.",
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[encori_evidence] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--score-modes", nargs="+", default=None, choices=SCORE_MODES)
    ap.add_argument("--clip-min", type=int, default=2)
    ap.add_argument("--top-k", type=int, default=15)
    ap.add_argument(
        "--prolif-key",
        default="e2f_g2m",
        choices=["e2f_g2m", "mki67", "ortho_noE2F_MYC", "pc1"],
    )
    args = ap.parse_args()
    run(
        out_dir=args.out_dir,
        score_modes=args.score_modes,
        clip_min=args.clip_min,
        top_k=args.top_k,
        prolif_key=args.prolif_key,
    )


if __name__ == "__main__":
    main()
