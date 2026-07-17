"""Dual-spine comparison: miRTar M0 (S0/S1) vs ENCORI M0′ (S2).

Run:
    .venv/bin/python3 -m mirna_hallmark.dual_spine_comparison
    .venv/bin/python3 -m mirna_hallmark.dual_spine_comparison --spines S0 S1 --alpha 0 0.5
    .venv/bin/python3 -m mirna_hallmark.dual_spine_comparison --spines S2 --clip-min 2 --beta 0.5
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from itertools import product
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.encori_edges import EncoriM0Spec, build_encori_m0_edges, enc_depth_lookup, write_encori_pair_table
from mirna_hallmark.analyses.evidence_dev.encori_mirtar_comparison import run as refresh_hub_comparison
from mirna_hallmark.evidence_scoring import (
    apply_encori_boost_to_edges,
    apply_mirtar_boost_to_edges,
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

OUT_DIR = C.OUTPUT_ROOT / "dual_spine_comparison"
DEFAULT_TARGET_NORM: TargetNormMode = "evidence_mass"
HUB_GENES = tuple(HUB_ROUTES.keys())
S1_BASE_SCORER = "confidence_logclass"
S1_ALPHA_GRID = (0.0, 0.25, 0.5, 1.0)
S2_CLIP_MIN_GRID = (1, 2, 3)
S2_BETA_GRID = (0.0, 0.5, 1.0)
S2_K_GRID = (15, 25, 50)


def _run_id(spine: str, **params) -> str:
    parts = [spine]
    for k, v in sorted(params.items()):
        if v is None:
            continue
        parts.append(f"{k}={v}")
    return "__".join(parts)


def _pair_key(df: pd.DataFrame) -> set[Tuple[str, str]]:
    if df.empty:
        return set()
    return set(zip(df["miRNA"].astype(str), df["gene"].astype(str)))


def _edge_weight_table(edges: pd.DataFrame, target_norm: TargetNormMode) -> pd.DataFrame:
    if edges.empty:
        return edges
    mass = mirna_mass_denominators(edges, target_norm)
    w = _edge_weights(edges, mass, target_norm)
    out = edges[["miRNA", "gene", "evidence_score"]].copy()
    out["edge_w"] = w
    return out


def _coupling_summary(coupling: pd.DataFrame, run_id: str, spine: str) -> pd.DataFrame:
    rows = []
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub] if not coupling.empty else coupling
        if g.empty:
            rows.append({
                "run_id": run_id, "spine": spine, "subtype": sub,
                "n_neg_sig": 0, "n_key_neg_sig": 0,
                "median_rho_prolif_cn_wsd_adj": float("nan"),
            })
            continue
        neg = g.loc[(g["rho_prolif_cn_wsd_adj"] < 0) & (g["q_prolif_cn_wsd_adj"] < 0.10)]
        neg_key = neg.loc[neg["key_hallmark"] == True]  # noqa: E712
        rows.append({
            "run_id": run_id,
            "spine": spine,
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_key_neg_sig": int(len(neg_key)),
            "median_rho_prolif_cn_wsd_adj": float(
                pd.to_numeric(g["rho_prolif_cn_wsd_adj"], errors="coerce").median()
            ),
        })
    return pd.DataFrame(rows)


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
        lo, hi = float(vals.min()), float(vals.max())
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


def _build_s0_edges(genes: Sequence[str], mirna_df: pd.DataFrame) -> pd.DataFrame:
    return build_m0_edges(genes, mirna_df, scorer_id="tiered_permissive")


def _build_s1_edges(
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
    m0_edges: pd.DataFrame,
    enc_lookup: pd.DataFrame,
    *,
    alpha: float,
    base_scorer: str = S1_BASE_SCORER,
) -> pd.DataFrame:
    summary = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, low_memory=False)
    scored = apply_scorer(summary, base_scorer)
    base = m0_edges[["miRNA", "gene"]].merge(
        scored[["miRNA", "gene", "evidence_score"]],
        on=["miRNA", "gene"],
        how="left",
    )
    base["evidence_score"] = pd.to_numeric(base["evidence_score"], errors="coerce").fillna(0.0)
    return apply_encori_boost_to_edges(base, enc_lookup, alpha=alpha)


def _build_s2_edges(
    genes: Sequence[str],
    mirna_df: pd.DataFrame,
    mirtar_summary: pd.DataFrame,
    *,
    clip_min: int,
    beta: float,
    top_k: int,
) -> pd.DataFrame:
    spec = EncoriM0Spec(clip_min=clip_min, top_k_program=top_k, top_k_hub=50)
    enc_edges = build_encori_m0_edges(genes, mirna_df, spec=spec)
    return apply_mirtar_boost_to_edges(enc_edges, mirtar_summary, beta=beta)


def _export_disagreement(
    s0_edges: pd.DataFrame,
    s2_edges: pd.DataFrame,
    gene: str,
    out_path: Path,
) -> None:
    s0_g = s0_edges.loc[s0_edges["gene"] == gene, ["miRNA", "gene", "evidence_score"]]
    s2_g = s2_edges.loc[s2_edges["gene"] == gene, ["miRNA", "gene", "evidence_score"]]
    s0_keys = _pair_key(s0_g)
    s2_keys = _pair_key(s2_g)
    rows = []
    for m, g in sorted(s0_keys - s2_keys):
        sc = s0_g.loc[(s0_g["miRNA"] == m) & (s0_g["gene"] == g), "evidence_score"]
        rows.append({
            "miRNA": m, "gene": g, "status": "S0_only",
            "evidence_score_s0": float(sc.iloc[0]) if len(sc) else float("nan"),
            "evidence_score_s2": float("nan"),
        })
    for m, g in sorted(s2_keys - s0_keys):
        sc = s2_g.loc[(s2_g["miRNA"] == m) & (s2_g["gene"] == g), "evidence_score"]
        rows.append({
            "miRNA": m, "gene": g, "status": "S2_only",
            "evidence_score_s0": float("nan"),
            "evidence_score_s2": float(sc.iloc[0]) if len(sc) else float("nan"),
        })
    for m, g in sorted(s0_keys & s2_keys):
        s0_sc = s0_g.loc[(s0_g["miRNA"] == m) & (s0_g["gene"] == g), "evidence_score"]
        s2_sc = s2_g.loc[(s2_g["miRNA"] == m) & (s2_g["gene"] == g), "evidence_score"]
        rows.append({
            "miRNA": m, "gene": g, "status": "shared",
            "evidence_score_s0": float(s0_sc.iloc[0]) if len(s0_sc) else float("nan"),
            "evidence_score_s2": float(s2_sc.iloc[0]) if len(s2_sc) else float("nan"),
        })
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)


def _evaluate_decision(
    coupling_summary: pd.DataFrame,
    hub_rho: pd.DataFrame,
    *,
    s0_run_id: str,
) -> dict:
    """Apply plan decision rules; return structured verdict."""
    s0_basal = coupling_summary.loc[
        (coupling_summary["run_id"] == s0_run_id) & (coupling_summary["subtype"] == "Basal")
    ]
    s0_key = int(s0_basal.iloc[0]["n_key_neg_sig"]) if len(s0_basal) else 0
    s0_neg = int(s0_basal.iloc[0]["n_neg_sig"]) if len(s0_basal) else 0

    s0_hub = hub_rho.loc[hub_rho["run_id"] == s0_run_id]
    s0_median_rho = float(s0_hub["spearman_edge_w_vs_abs_contrib"].median()) if len(s0_hub) else float("nan")
    s0_irf1 = s0_hub.loc[s0_hub["gene"] == "IRF1", "spearman_edge_w_vs_abs_contrib"]
    s0_irf1_rho = float(s0_irf1.iloc[0]) if len(s0_irf1) else float("nan")

    candidates = []
    for rid in coupling_summary["run_id"].unique():
        if rid == s0_run_id:
            continue
        basal = coupling_summary.loc[
            (coupling_summary["run_id"] == rid) & (coupling_summary["subtype"] == "Basal")
        ]
        if basal.empty:
            continue
        key_neg = int(basal.iloc[0]["n_key_neg_sig"])
        neg_sig = int(basal.iloc[0]["n_neg_sig"])
        hub_sub = hub_rho.loc[hub_rho["run_id"] == rid]
        med_rho = float(hub_sub["spearman_edge_w_vs_abs_contrib"].median()) if len(hub_sub) else float("nan")
        irf1_sub = hub_sub.loc[hub_sub["gene"] == "IRF1", "spearman_edge_w_vs_abs_contrib"]
        irf1_rho = float(irf1_sub.iloc[0]) if len(irf1_sub) else float("nan")

        passes_key = key_neg >= 8
        passes_breadth = neg_sig >= 40
        rho_gain = med_rho - s0_median_rho if np.isfinite(med_rho) and np.isfinite(s0_median_rho) else float("nan")
        passes_hub = (np.isfinite(rho_gain) and rho_gain >= 0.05) or (np.isfinite(irf1_rho) and irf1_rho >= 0.25)
        promote = passes_key and passes_breadth and passes_hub

        candidates.append({
            "run_id": rid,
            "spine": coupling_summary.loc[coupling_summary["run_id"] == rid, "spine"].iloc[0],
            "basal_neg_sig": neg_sig,
            "basal_key_neg_sig": key_neg,
            "hub_median_rho": med_rho,
            "hub_median_rho_gain_vs_s0": rho_gain,
            "irf1_rho": irf1_rho,
            "promote": promote,
        })

    winner = [c for c in candidates if c["promote"]]
    if not winner:
        verdict = "S0_retained"
        winner_id = s0_run_id
    else:
        winner.sort(key=lambda c: (c["basal_neg_sig"], c["hub_median_rho_gain_vs_s0"]), reverse=True)
        winner_id = winner[0]["run_id"]
        verdict = f"{winner[0]['spine']}_promoted"

    return {
        "verdict": verdict,
        "winner_run_id": winner_id,
        "s0_baseline": {"basal_neg_sig": s0_neg, "basal_key_neg_sig": s0_key, "hub_median_rho": s0_median_rho, "irf1_rho": s0_irf1_rho},
        "candidates": candidates,
        "note": "DISCOVERY_REGISTRY and pressure_build default unchanged until manual review.",
    }


def _write_decision_readme(out_root: Path, decision: dict) -> None:
    lines = [
        "# Dual-spine comparison — decision memo",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        f"**Verdict:** `{decision['verdict']}`",
        f"**Winner run_id:** `{decision['winner_run_id']}`",
        "",
        "## S0 baseline",
        "",
        f"- Basal neg-sig: {decision['s0_baseline']['basal_neg_sig']}/50",
        f"- Basal key neg-sig: {decision['s0_baseline']['basal_key_neg_sig']}/8",
        f"- Hub median rank-ρ: {decision['s0_baseline']['hub_median_rho']:.3f}",
        f"- IRF1 rank-ρ: {decision['s0_baseline']['irf1_rho']:.3f}",
        "",
        "## Decision rules",
        "",
        "Promote only if ALL: key 8/8, Basal breadth ≥ 40/50, hub median ρ gain ≥ 0.05 OR IRF1 ρ ≥ 0.25.",
        "",
        decision["note"],
        "",
        "## Candidate runs",
        "",
    ]
    for c in decision.get("candidates", []):
        flag = "PROMOTE" if c["promote"] else "hold"
        lines.append(
            f"- `{c['run_id']}` ({c['spine']}): Basal {c['basal_neg_sig']}/50, key {c['basal_key_neg_sig']}/8, "
            f"hub ρ={c['hub_median_rho']:.3f} (Δ={c['hub_median_rho_gain_vs_s0']:.3f}), "
            f"IRF1 ρ={c['irf1_rho']:.3f} → **{flag}**"
        )
    (out_root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(
    *,
    out_dir: Path | None = None,
    spines: Sequence[str] | None = None,
    alpha_grid: Sequence[float] | None = None,
    clip_min_grid: Sequence[int] | None = None,
    beta_grid: Sequence[float] | None = None,
    k_grid: Sequence[int] | None = None,
    refresh_hub_table: bool = True,
) -> None:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)
    disagree_dir = out_root / "disagreement"
    disagree_dir.mkdir(parents=True, exist_ok=True)

    if refresh_hub_table:
        refresh_hub_comparison()

    write_encori_pair_table()

    spine_set = set(spines or ("S0", "S1", "S2"))
    alphas = list(alpha_grid or S1_ALPHA_GRID)
    clip_mins = list(clip_min_grid or S2_CLIP_MIN_GRID)
    betas = list(beta_grid or S2_BETA_GRID)
    k_vals = list(k_grid or S2_K_GRID)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted({g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    mirtar_summary = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, low_memory=False)
    enc_lookup = enc_depth_lookup(genes=hs.universe)

    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )["ago_gate"]

    s0_edges = _build_s0_edges(hs.universe, mirna)
    s0_keys = _pair_key(s0_edges)
    print(f"[dual_spine] S0 M0 edges: {len(s0_edges):,} pairs")

    specs: List[dict] = []
    if "S0" in spine_set:
        specs.append({"spine": "S0", "run_id": "S0"})
    if "S1" in spine_set:
        for alpha in alphas:
            if alpha == 0 and "S0" in spine_set:
                continue
            specs.append({
                "spine": "S1",
                "run_id": _run_id("S1", alpha=alpha, base=S1_BASE_SCORER),
                "alpha": alpha,
                "base_scorer": S1_BASE_SCORER,
            })
    if "S2" in spine_set:
        for clip_min, beta, top_k in product(clip_mins, betas, k_vals):
            specs.append({
                "spine": "S2",
                "run_id": _run_id("S2", clip_min=clip_min, beta=beta, top_k=top_k),
                "clip_min": clip_min,
                "beta": beta,
                "top_k": top_k,
            })

    coupling_all = []
    summary_all = []
    hierarchy_all = []
    hub_all = []
    basal_key_all = []
    jaccard_rows = []
    s2_edge_cache: dict[str, pd.DataFrame] = {}

    for spec in specs:
        rid = spec["run_id"]
        spine = spec["spine"]
        print(f"[dual_spine] {rid}")

        if spine == "S0":
            edges = s0_edges
        elif spine == "S1":
            edges = _build_s1_edges(
                hs.universe, mirna, s0_edges, enc_lookup,
                alpha=float(spec["alpha"]),
                base_scorer=spec.get("base_scorer", S1_BASE_SCORER),
            )
        else:
            edges = _build_s2_edges(
                hs.universe, mirna, mirtar_summary,
                clip_min=int(spec["clip_min"]),
                beta=float(spec["beta"]),
                top_k=int(spec["top_k"]),
            )
            s2_edge_cache[rid] = edges
            inter = len(s0_keys & _pair_key(edges))
            union = len(s0_keys | _pair_key(edges))
            jaccard_rows.append({
                "run_id": rid,
                "s0_pairs": len(s0_keys),
                "s2_pairs": len(edges),
                "intersection": inter,
                "union": union,
                "jaccard": inter / union if union else float("nan"),
            })

        hier = _hierarchy_metrics(edges, DEFAULT_TARGET_NORM)
        hier.update({"run_id": rid, "spine": spine, **{k: v for k, v in spec.items() if k not in ("run_id", "spine")}})
        hierarchy_all.append(hier)

        hub_df = _hub_rank_rho(edges, hs.universe, mirna, DEFAULT_TARGET_NORM).assign(run_id=rid, spine=spine)
        hub_all.append(hub_df)

        gp = compute_gene_pressure(
            edges, mirna, genes=list(hs.universe),
            expr_mode="softmax_z_logrpm",
            target_norm=DEFAULT_TARGET_NORM,
            aggregate="sum",
        )
        shared = gp.columns.intersection(gate.index)
        hp = hallmark_pressure_matrix(gp[shared].mul(gate.reindex(shared), axis=1), hs)

        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["run_id"] = rid
        coupling["spine"] = spine
        coupling_all.append(coupling)
        summary_all.append(_coupling_summary(coupling, rid, spine))

        basal = coupling.loc[
            (coupling["subtype"] == "Basal") & (coupling["hallmark_set"].isin(KEY_HALLMARKS))
        ].copy()
        basal["run_id"] = rid
        basal_key_all.append(basal)

        b = summary_all[-1].loc[summary_all[-1]["subtype"] == "Basal"]
        if not b.empty:
            print(
                f"  Basal: neg_sig={int(b.iloc[0]['n_neg_sig'])}/50, "
                f"key={int(b.iloc[0]['n_key_neg_sig'])}/8"
            )

    spec_manifest = pd.DataFrame([
        {k: v for k, v in s.items()} for s in specs
    ])
    spec_manifest.to_csv(out_root / "spine_spec_manifest.tsv", sep="\t", index=False)

    coupling_summary = pd.concat(summary_all, ignore_index=True)
    coupling_summary.to_csv(out_root / "coupling_summary_by_spine.tsv", sep="\t", index=False)
    pd.concat(coupling_all, ignore_index=True).to_csv(
        out_root / "coupling_detail_by_spine.tsv", sep="\t", index=False,
    )

    basal_wide = pd.concat(basal_key_all, ignore_index=True)
    pivot_rows = []
    run_ids = [s["run_id"] for s in specs]
    for h in KEY_HALLMARKS:
        row: dict = {"hallmark_set": h}
        for rid in run_ids:
            sub = basal_wide.loc[
                (basal_wide["hallmark_set"] == h) & (basal_wide["run_id"] == rid)
            ]
            if sub.empty:
                row[f"{rid}_rho"] = float("nan")
                row[f"{rid}_neg_sig"] = False
            else:
                r = sub.iloc[0]
                row[f"{rid}_rho"] = float(r["rho_prolif_cn_wsd_adj"])
                row[f"{rid}_neg_sig"] = bool(
                    (r["rho_prolif_cn_wsd_adj"] < 0) and (r["q_prolif_cn_wsd_adj"] < 0.10)
                )
        pivot_rows.append(row)
    pd.DataFrame(pivot_rows).to_csv(out_root / "basal_key_pivot.tsv", sep="\t", index=False)

    hub_concat = pd.concat(hub_all, ignore_index=True)
    hub_concat.to_csv(out_root / "hub_rank_correlation.tsv", sep="\t", index=False)
    pd.DataFrame(hierarchy_all).to_csv(out_root / "edge_hierarchy_by_spine.tsv", sep="\t", index=False)

    if jaccard_rows:
        pd.DataFrame(jaccard_rows).to_csv(out_root / "edge_set_jaccard.tsv", sep="\t", index=False)

    # Disagreement exports: best S2 by Basal breadth (or first S2)
    if s2_edge_cache:
        best_s2_rid = max(
            s2_edge_cache.keys(),
            key=lambda r: int(
                coupling_summary.loc[
                    (coupling_summary["run_id"] == r) & (coupling_summary["subtype"] == "Basal"),
                    "n_neg_sig",
                ].iloc[0]
                if len(coupling_summary.loc[
                    (coupling_summary["run_id"] == r) & (coupling_summary["subtype"] == "Basal")
                ]) else 0
            ),
        )
        s2_best = s2_edge_cache[best_s2_rid]
        for gene in HUB_GENES:
            _export_disagreement(
                s0_edges, s2_best, gene,
                disagree_dir / f"{gene}_edge_diff.tsv",
            )

    decision = _evaluate_decision(coupling_summary, hub_concat, s0_run_id="S0")
    _write_decision_readme(out_root, decision)

    manifest = {
        "module": "mirna_hallmark.dual_spine_comparison",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "spines": list(spine_set),
        "pressure": "softmax_z_logrpm; target_norm=evidence_mass; aggregate=sum",
        "s1_alpha_grid": alphas,
        "s2_clip_min_grid": clip_mins,
        "s2_beta_grid": betas,
        "s2_k_grid": k_vals,
        "decision": decision,
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[dual_spine] verdict={decision['verdict']} -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--spines", nargs="+", default=None, choices=["S0", "S1", "S2"])
    ap.add_argument("--alpha", nargs="+", type=float, default=None, help="S1 ENCORI boost grid")
    ap.add_argument("--clip-min", nargs="+", type=int, default=None, dest="clip_min")
    ap.add_argument("--beta", nargs="+", type=float, default=None, help="S2 miRTar boost grid")
    ap.add_argument("--top-k", nargs="+", type=int, default=None, dest="top_k")
    ap.add_argument("--no-refresh-hub-table", action="store_true")
    args = ap.parse_args()
    run(
        out_dir=args.out_dir,
        spines=args.spines,
        alpha_grid=args.alpha,
        clip_min_grid=args.clip_min,
        beta_grid=args.beta,
        k_grid=args.top_k,
        refresh_hub_table=not args.no_refresh_hub_table,
    )


if __name__ == "__main__":
    main()
