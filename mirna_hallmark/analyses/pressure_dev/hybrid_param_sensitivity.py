"""One-at-a-time sensitivity for M7 hybrid hyperparameters (+ optional weighting variants).

Grids ``beta_ts_boost``, ``gamma_orphan``, and ``min_ts`` around M7 defaults while
holding other knobs fixed. Also tests M7 under extended ``expr_mode`` /
``target_norm`` combinations (evidence-mass, TS mass, absolute anchors).

Outputs under ``output/hybrid_param_sensitivity/``:
  - ``param_sensitivity_summary.tsv`` — Basal Hallmark counts + IRF1 partial ρ
  - ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.hybrid_pressure import (
    HybridPressureSpec,
    build_hybrid_edges,
    compute_hybrid_gene_pressure,
)
from mirna_hallmark.pressure_engine import AggregateMode, ExprMode, TargetNormMode
from mirna_hallmark.robustness_checks import (
    HUB_ROUTES,
    _proliferation_proxies,
    gene_pressure_route_partial_corr,
    hallmark_coupling_by_subtype,
)

OUT_DIR = C.OUTPUT_ROOT / "hybrid_param_sensitivity"


@dataclass(frozen=True)
class HybridRunSpec:
    run_id: str
    mode: str = "M7"
    beta_ts_boost: float = 0.5
    gamma_orphan: float = 0.25
    min_ts: float = 0.25
    expr_mode: ExprMode = "softmax_z"
    target_norm: TargetNormMode = "degree"
    description: str = ""


def _default_grid() -> List[HybridRunSpec]:
    base = HybridRunSpec("m7_default", description="M7 defaults")
    rows = [base]
    for b in (0.25, 1.0):
        rows.append(HybridRunSpec(
            f"m7_beta_{b}", beta_ts_boost=b,
            description=f"M7 beta_ts_boost={b}",
        ))
    for g in (0.1, 0.5):
        rows.append(HybridRunSpec(
            f"m7_gamma_{g}", gamma_orphan=g,
            description=f"M7 gamma_orphan={g}",
        ))
    for t in (0.15, 0.5):
        rows.append(HybridRunSpec(
            f"m7_min_ts_{t}", min_ts=t,
            description=f"M7 min_ts={t}",
        ))
    rows.extend([
        HybridRunSpec(
            "m7_combined_mass", target_norm="combined_mass",
            description="M7 + combined_mass target norm (miRTar+TS)",
        ),
        HybridRunSpec(
            "m7_evidence_mass", target_norm="evidence_mass",
            description="M7 + evidence_mass target norm",
        ),
        HybridRunSpec(
            "m7_logrpm_degree", expr_mode="softmax_z_logrpm",
            description="M7 + softmax_z_logrpm absolute anchor",
        ),
        HybridRunSpec(
            "m7_absratio_combined", expr_mode="softmax_z_absratio",
            target_norm="combined_mass",
            description="M7 + abs ratio anchor + combined_mass norm",
        ),
        HybridRunSpec(
            "m0_default", mode="M0",
            description="M0 miRTar baseline (post arm-resolve)",
        ),
    ])
    return rows


def _hybrid_spec(run: HybridRunSpec) -> HybridPressureSpec:
    return HybridPressureSpec(
        mode=run.mode,
        min_ts=run.min_ts,
        beta_ts_boost=run.beta_ts_boost,
        gamma_orphan=run.gamma_orphan,
    )


def _irf1_basal_row(hub_routes: pd.DataFrame) -> Tuple[float, float]:
    sub = hub_routes.loc[(hub_routes["target"] == "IRF1") & (hub_routes["scope"] == "Basal")]
    if sub.empty:
        return float("nan"), float("nan")
    r = sub.iloc[0]
    return float(r["rho_e2f_g2m_CN"]), float(r["p_e2f_g2m_CN"])


def run(
    *,
    specs: Sequence[HybridRunSpec] = None,
    out_dir: Optional[Path] = None,
) -> pd.DataFrame:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)
    specs = list(specs or _default_grid())

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    hallmark_edges = D.load_hallmark_edges()
    proxies = _proliferation_proxies(rna, hs)
    cnv_genes = sorted(set(HUB_ROUTES) | {g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(cnv_genes)
    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )["ago_gate"]

    rows: List[dict] = []
    for rs in specs:
        print(f"[hybrid_param_sensitivity] {rs.run_id}: {rs.description}")
        hspec = _hybrid_spec(rs)
        gp, edges = compute_hybrid_gene_pressure(
            hs.universe,
            spec=hspec,
            mirna=mirna,
            gate=gate,
            hallmark_edges=hallmark_edges,
            expr_mode=rs.expr_mode,
            target_norm=rs.target_norm,
            aggregate=C.PRESSURE_AGGREGATE,  # type: ignore[arg-type]
        )
        hp = hallmark_pressure_matrix(gp, hs)
        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        basal = coupling.loc[coupling["subtype"] == "Basal"]
        neg = basal.loc[
            (basal["rho_prolif_cn_wsd_adj"] < 0) & (basal["q_prolif_cn_wsd_adj"] < 0.10)
        ]
        key_neg = neg.loc[neg["key_hallmark"] == True]  # noqa: E712

        hub = gene_pressure_route_partial_corr(rna, clinical, proxies, gp, cnv)
        irf1_rho, irf1_p = _irf1_basal_row(hub)

        rows.append({
            "run_id": rs.run_id,
            "mode": rs.mode,
            "beta_ts_boost": rs.beta_ts_boost,
            "gamma_orphan": rs.gamma_orphan,
            "min_ts": rs.min_ts,
            "expr_mode": rs.expr_mode,
            "target_norm": rs.target_norm,
            "n_edges": len(edges),
            "basal_neg_sig": int(len(neg)),
            "basal_key_neg_sig": int(len(key_neg)),
            "basal_median_rho": float(basal["rho_prolif_cn_wsd_adj"].median()),
            "irf1_basal_rho_cn_prolif": irf1_rho,
            "irf1_basal_p_cn_prolif": irf1_p,
            "description": rs.description,
        })
        print(f"  edges={len(edges):,} Basal={len(neg)}/50 IRF1 ρ={irf1_rho:.3f} p={irf1_p:.4g}")

    summary = pd.DataFrame(rows)
    summary.to_csv(out_root / "param_sensitivity_summary.tsv", sep="\t", index=False)
    manifest = {
        "module": "mirna_hallmark.hybrid_param_sensitivity",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "n_runs": len(specs),
        "note": "One-at-a-time M7 hyperparameter grid; includes M0 baseline",
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[hybrid_param_sensitivity] done -> {out_root}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
