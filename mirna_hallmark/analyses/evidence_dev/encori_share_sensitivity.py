"""Lane 2: ENCORI depth in softmax share on fixed M0 edges (β_share grid).

Keeps S0 edge_w (tiered_permissive + evidence_mass). Adds ENCORI depth only to
gene-local softmax logits:

    logit(m,g,s) = (x_m,s − cohort_median_m) + β_share · enc_depth(m,g)
    c(m,g,s) = softmax_g(logit) × z(m,s) × logRPM(m,s) × edge_w(m,g)

ENCORI depth is 0 when no parquet row exists (~75% of M0 pairs). Do **not** combine
with S1 edge-weight boost (α on evidence_score) at full strength.

Run:
    .venv/bin/python3 -m mirna_hallmark.encori_share_sensitivity
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
from mirna_hallmark.encori_edges import enc_depth_lookup
from itertools import product

from mirna_hallmark.evidence_scoring import (
    apply_encori_boost_to_edges,
    apply_scorer,
    build_m0_edges,
)
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

OUT_DIR = C.OUTPUT_ROOT / "encori_share_sensitivity"
EXPR_MODE = "softmax_z_logrpm_enc"
TARGET_NORM: TargetNormMode = "evidence_mass"
HUB_GENES = tuple(HUB_ROUTES.keys())
DEFAULT_BETA_GRID = (0.0, 0.25, 0.5, 1.0)
DEFAULT_ALPHA_GRID = (0.0,)
S1_BASE_SCORER = "confidence_logclass"


def build_lane_edges(
    m0_edges: pd.DataFrame,
    enc_lookup: pd.DataFrame,
    *,
    alpha: float,
    beta_share: float,
    base_scorer: str = S1_BASE_SCORER,
) -> pd.DataFrame:
    """M0 pairs; optional S1 edge boost (α) + Lane 2 share bonus (β_share)."""
    if alpha > 0:
        summary = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, low_memory=False)
        scored = apply_scorer(summary, base_scorer)
        base = m0_edges[["miRNA", "gene"]].merge(
            scored[["miRNA", "gene", "evidence_score"]],
            on=["miRNA", "gene"],
            how="left",
        )
        base["evidence_score"] = pd.to_numeric(base["evidence_score"], errors="coerce").fillna(0.0)
        edges = apply_encori_boost_to_edges(base, enc_lookup, alpha=alpha)
    else:
        edges = m0_edges[["miRNA", "gene", "evidence_score"]].copy()
    return attach_enc_share_bonus(edges, enc_lookup, beta_share=beta_share)


def attach_enc_share_bonus(
    edges: pd.DataFrame,
    enc_lookup: pd.DataFrame,
    *,
    beta_share: float,
) -> pd.DataFrame:
    out = edges.merge(enc_lookup, on=["miRNA", "gene"], how="left")
    depth = pd.to_numeric(out["enc_depth"], errors="coerce").fillna(0.0)
    out["share_logit_bonus"] = float(beta_share) * depth
    return out


def _run_id(alpha: float, beta_share: float) -> str:
    if alpha > 0:
        return f"S1_share__alpha={alpha}__beta_share={beta_share}"
    return f"Lane2__beta_share={beta_share}"


def _coupling_summary(
    coupling: pd.DataFrame, run_id: str, alpha: float, beta_share: float,
) -> pd.DataFrame:
    rows = []
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub] if not coupling.empty else coupling
        if g.empty:
            rows.append({
                "run_id": run_id, "alpha": alpha, "beta_share": beta_share, "subtype": sub,
                "n_neg_sig": 0, "n_key_neg_sig": 0,
                "median_rho_prolif_cn_wsd_adj": float("nan"),
            })
            continue
        neg = g.loc[(g["rho_prolif_cn_wsd_adj"] < 0) & (g["q_prolif_cn_wsd_adj"] < 0.10)]
        neg_key = neg.loc[neg["key_hallmark"] == True]  # noqa: E712
        rows.append({
            "run_id": run_id,
            "alpha": alpha,
            "beta_share": beta_share,
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_key_neg_sig": int(len(neg_key)),
            "median_rho_prolif_cn_wsd_adj": float(
                pd.to_numeric(g["rho_prolif_cn_wsd_adj"], errors="coerce").median()
            ),
        })
    return pd.DataFrame(rows)


def _hub_rank_rho(
    edges: pd.DataFrame,
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for gene in HUB_GENES:
        if gene not in genes:
            continue
        ge = edges.loc[edges["gene"] == gene]
        if ge.empty:
            continue
        mass = mirna_mass_denominators(edges, TARGET_NORM)
        ew = ge[["miRNA", "gene", "evidence_score"]].copy()
        ew["edge_w"] = _edge_weights(ge, mass, TARGET_NORM)
        contrib = compute_gene_pressure_contributions(
            [gene],
            edges=edges,
            expr_mode=EXPR_MODE,
            target_norm=TARGET_NORM,
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
    beta_grid: Sequence[float] | None = None,
    alpha_grid: Sequence[float] | None = None,
) -> None:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)
    betas = list(beta_grid or DEFAULT_BETA_GRID)
    alphas = list(alpha_grid or DEFAULT_ALPHA_GRID)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted({g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    enc_lookup = enc_depth_lookup(genes=hs.universe)

    m0_edges = build_m0_edges(hs.universe, mirna, scorer_id="tiered_permissive")
    n_enc = m0_edges.merge(enc_lookup, on=["miRNA", "gene"], how="inner").shape[0]
    print(f"[encori_share] M0 edges: {len(m0_edges):,}; with ENCORI depth: {n_enc:,}")

    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )["ago_gate"]

    spec_rows: List[dict] = []
    coupling_all = []
    summary_all = []
    hub_all = []
    basal_key_all = []

    for alpha, beta in product(alphas, betas):
        rid = _run_id(alpha, beta)
        print(f"[encori_share] {rid}")
        edges = build_lane_edges(m0_edges, enc_lookup, alpha=alpha, beta_share=beta)

        spec_rows.append({
            "run_id": rid,
            "alpha": alpha,
            "beta_share": beta,
            "base_scorer": S1_BASE_SCORER if alpha > 0 else "tiered_permissive",
            "expr_mode": EXPR_MODE,
            "target_norm": TARGET_NORM,
            "edge_set": "M0_tiered_permissive",
            "n_edges": len(edges),
            "n_pairs_with_encori": int((edges["share_logit_bonus"] > 0).sum()),
        })

        hub_df = _hub_rank_rho(edges, hs.universe, mirna).assign(
            run_id=rid, alpha=alpha, beta_share=beta,
        )
        hub_all.append(hub_df)

        gp = compute_gene_pressure(
            edges,
            mirna,
            genes=list(hs.universe),
            expr_mode=EXPR_MODE,
            target_norm=TARGET_NORM,
            aggregate="sum",
        )
        shared = gp.columns.intersection(gate.index)
        hp = hallmark_pressure_matrix(gp[shared].mul(gate.reindex(shared), axis=1), hs)

        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["run_id"] = rid
        coupling["alpha"] = alpha
        coupling["beta_share"] = beta
        coupling_all.append(coupling)
        summary_all.append(_coupling_summary(coupling, rid, alpha, beta))

        basal = coupling.loc[
            (coupling["subtype"] == "Basal") & (coupling["hallmark_set"].isin(KEY_HALLMARKS))
        ].copy()
        basal["run_id"] = rid
        basal_key_all.append(basal)

        b = summary_all[-1].loc[summary_all[-1]["subtype"] == "Basal"]
        med_rho = hub_df["spearman_edge_w_vs_abs_contrib"].median()
        irf1 = hub_df.loc[hub_df["gene"] == "IRF1", "spearman_edge_w_vs_abs_contrib"]
        irf1_r = float(irf1.iloc[0]) if len(irf1) else float("nan")
        if not b.empty:
            print(
                f"  Basal: neg_sig={int(b.iloc[0]['n_neg_sig'])}/50, "
                f"key={int(b.iloc[0]['n_key_neg_sig'])}/8, "
                f"hub_median_rho={med_rho:.3f}, IRF1_rho={irf1_r:.2f}"
            )

    pd.DataFrame(spec_rows).to_csv(out_root / "spec_manifest.tsv", sep="\t", index=False)
    summary = pd.concat(summary_all, ignore_index=True)
    summary.to_csv(out_root / "coupling_summary.tsv", sep="\t", index=False)
    pd.concat(coupling_all, ignore_index=True).to_csv(
        out_root / "coupling_detail.tsv", sep="\t", index=False,
    )
    pd.concat(hub_all, ignore_index=True).to_csv(
        out_root / "hub_rank_correlation.tsv", sep="\t", index=False,
    )

    basal_wide = pd.concat(basal_key_all, ignore_index=True)
    pivot_rows = []
    for h in KEY_HALLMARKS:
        row: dict = {"hallmark_set": h}
        for alpha, beta in product(alphas, betas):
            rid = _run_id(alpha, beta)
            sub = basal_wide.loc[
                (basal_wide["hallmark_set"] == h) & (basal_wide["run_id"] == rid)
            ]
            col = f"a{alpha}_b{beta}"
            if sub.empty:
                row[f"{col}_rho"] = float("nan")
                row[f"{col}_neg_sig"] = False
            else:
                r = sub.iloc[0]
                row[f"{col}_rho"] = float(r["rho_prolif_cn_wsd_adj"])
                row[f"{col}_neg_sig"] = bool(
                    (r["rho_prolif_cn_wsd_adj"] < 0) and (r["q_prolif_cn_wsd_adj"] < 0.10)
                )
        pivot_rows.append(row)
    pd.DataFrame(pivot_rows).to_csv(out_root / "basal_key_pivot.tsv", sep="\t", index=False)

    # README
    s0_row = summary.loc[
        (summary.alpha == 0) & (summary.beta_share == 0) & (summary.subtype == "Basal")
    ]
    if len(s0_row):
        s0_neg = int(s0_row.iloc[0]["n_neg_sig"])
    else:
        s0_neg = 0
    lines = [
        "# ENCORI in edge_w (S1 α) and/or softmax share (β_share) on M0 edges",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Formula",
        "",
        "```",
        "edge_w  from confidence_logclass + α·enc_depth  (α=0 → tiered_permissive S0)",
        "logit   = (x−med) + β_share·enc_depth",
        "c       = softmax_g(logit) × z × logRPM × edge_w",
        "```",
        "",
        f"M0 pairs: {len(m0_edges):,}; ENCORI hits: {n_enc:,} ({100*n_enc/len(m0_edges):.1f}%)",
        "",
        "## Basal coupling (α × β_share grid)",
        "",
    ]
    hub_concat = pd.concat(hub_all)
    for alpha, beta in product(alphas, betas):
        sub = summary.loc[
            (summary.alpha == alpha) & (summary.beta_share == beta) & (summary.subtype == "Basal")
        ]
        if sub.empty:
            continue
        r = sub.iloc[0]
        h = hub_concat.loc[(hub_concat.alpha == alpha) & (hub_concat.beta_share == beta)]
        med = h["spearman_edge_w_vs_abs_contrib"].median()
        irf1 = h.loc[h.gene == "IRF1", "spearman_edge_w_vs_abs_contrib"]
        irf1_r = float(irf1.iloc[0]) if len(irf1) else float("nan")
        tag = ""
        if alpha > 0 and beta == 0:
            tag = " [S1 only]"
        elif alpha > 0 and beta > 0:
            tag = " [S1+share]"
        lines.append(
            f"- α={alpha}, β={beta}{tag}: neg_sig={int(r.n_neg_sig)}/50, "
            f"key={int(r.n_key_neg_sig)}/8, hub_ρ={med:.3f}, IRF1_ρ={irf1_r:.2f}"
        )
    lines += [
        "",
        f"α=0, β=0 reference Basal: {s0_neg}/50 (M0 tiered_permissive edge_w).",
    ]
    (out_root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest = {
        "module": "mirna_hallmark.encori_share_sensitivity",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "beta_share_grid": betas,
        "alpha_grid": alphas,
        "expr_mode": EXPR_MODE,
        "target_norm": TARGET_NORM,
        "m0_edges": len(m0_edges),
        "m0_with_encori": n_enc,
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[encori_share] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--beta-share", nargs="+", type=float, default=None, dest="beta_share")
    ap.add_argument("--alpha", nargs="+", type=float, default=None, dest="alpha")
    args = ap.parse_args()
    run(out_dir=args.out_dir, beta_grid=args.beta_share, alpha_grid=args.alpha)


if __name__ == "__main__":
    main()
