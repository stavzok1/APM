"""Compare pressure weighting schemes: z vs softmax variants × target-count norm.

Baseline (legacy): ``expr_mode=z``, ``target_norm=none``, ``aggregate=sum`` on
miRTarBase M0 edges — matches ``pressure_build`` / registered MH claims.

Each spec rebuilds gene → Hallmark pressure (AGO-gated) and replays
``hallmark_coupling_by_subtype`` (CPE + HRD + E2F/G2M proliferation + target CN +
**TCGA sequencing batch** (``plate_both``, via ``robustness_checks._partial_batch``) +
within-subtype z-scored pressure/expression).

Outputs (``output/pressure_sensitivity/``):
- ``spec_manifest.tsv`` — grid definition
- ``coupling_by_spec_subtype.tsv`` — full Hallmark × subtype × spec
- ``coupling_summary_by_spec.tsv`` — neg-sig counts + median rho per spec × subtype
- ``basal_key_programs_by_spec.tsv`` — 8 key Hallmarks @ Basal
- ``method_manifest.json``

Run:
    .venv/bin/python3 -m mirna_hallmark.pressure_sensitivity
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.hybrid_pressure import HybridPressureSpec, build_hybrid_edges
from mirna_hallmark.pressure_engine import AggregateMode, ExprMode, TargetNormMode, compute_gene_pressure
from mirna_hallmark.robustness_checks import (
    KEY_HALLMARKS,
    PAM50_SUBTYPES,
    HUB_ROUTES,
    _proliferation_proxies,
    gene_pressure_route_partial_corr,
    hallmark_coupling_by_subtype,
)

OUT_DIR = C.OUTPUT_ROOT / "pressure_sensitivity"


@dataclass(frozen=True)
class PressureSpec:
    spec_id: str
    expr_mode: ExprMode
    target_norm: TargetNormMode
    aggregate: AggregateMode
    description: str


DEFAULT_SPECS: List[PressureSpec] = [
    PressureSpec("z_sum_baseline", "z", "none", "sum",
                 "Legacy: cohort z × log1p(evidence), sum over arms"),
    PressureSpec("softmax_z_degree_sum", "softmax_z", "degree", "sum",
                 "Canonical: softmax share × z, degree norm"),
    PressureSpec("softmax_z_evidence_mass", "softmax_z", "evidence_mass", "sum",
                 "softmax_z; edge_w / log1p(Σ_g log1p(evidence)) per miRNA"),
    PressureSpec("softmax_z_ts_mass", "softmax_z", "ts_mass", "sum",
                 "softmax_z; edge_w / log1p(Σ_g log1p(TS)) per miRNA"),
    PressureSpec("softmax_z_combined_mass", "softmax_z", "combined_mass", "sum",
                 "softmax_z; (log1p(ev)+log1p(TS)) / log1p(Σ combined) per miRNA"),
    PressureSpec("softmax_z_logrpm_degree", "softmax_z_logrpm", "degree", "sum",
                 "softmax_z × log2(RPM+1) absolute anchor, degree norm"),
    PressureSpec("softmax_z_absratio_degree", "softmax_z_absratio", "degree", "sum",
                 "softmax_z × (RPM/median) ratio anchor (floor 0.25), degree norm"),
    PressureSpec("softmax_z_blend_degree", "softmax_z_blend", "degree", "sum",
                 "softmax_z with 0.5·share + 0.5 uniform blend on z multiplier"),
    PressureSpec("softmax_z_logrpm_evidence_mass", "softmax_z_logrpm", "evidence_mass", "sum",
                 "absolute logrpm anchor + evidence-mass norm"),
    PressureSpec("softmax_z_absratio_combined_mass", "softmax_z_absratio", "combined_mass", "sum",
                 "abs ratio anchor + miRTar+TS combined mass norm"),
]


def _m0_edges(genes: Sequence[str], hallmark_edges: pd.DataFrame) -> pd.DataFrame:
    spec = HybridPressureSpec(mode="M0")
    e = build_hybrid_edges(genes, spec=spec, hallmark_edges=hallmark_edges)
    e = e.rename(columns={"edge_weight": "evidence_score"})[["miRNA", "gene", "evidence_score"]]
    from mirna_hallmark.robustness_checks import _load_targetscan_weights

    ts = _load_targetscan_weights(genes)
    if ts.empty:
        e["ts_weight"] = 0.0
        return e
    return e.merge(ts, on=["miRNA", "gene"], how="left").assign(
        ts_weight=lambda d: pd.to_numeric(d["ts_weight"], errors="coerce").fillna(0.0)
    )


def _coupling_summary(coupling: pd.DataFrame, spec_id: str) -> pd.DataFrame:
    rows = []
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub]
        if g.empty:
            continue
        neg = g.loc[(g["rho_prolif_cn_wsd_adj"] < 0) & (g["q_prolif_cn_wsd_adj"] < 0.10)]
        key = g.loc[g["key_hallmark"] == True]  # noqa: E712
        neg_key = neg.loc[neg["key_hallmark"] == True]  # noqa: E712
        rows.append({
            "spec_id": spec_id,
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_hallmarks": int(len(g)),
            "n_key_neg_sig": int(len(neg_key)),
            "n_key_hallmarks": int(len(key)),
            "median_rho_prolif_cn_wsd_adj": float(pd.to_numeric(g["rho_prolif_cn_wsd_adj"], errors="coerce").median()),
        })
    return pd.DataFrame(rows)


def run(
    *,
    specs: Sequence[PressureSpec] = DEFAULT_SPECS,
    out_dir: Path | None = None,
) -> None:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    hallmark_edges = D.load_hallmark_edges()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted({g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    cnv_hub = D.load_cnv_target_genes(list(HUB_ROUTES.keys()))
    edges_m0 = _m0_edges(hs.universe, hallmark_edges)

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
    basal_key_all = []
    irf1_all = []

    for ps in specs:
        print(f"[pressure_sensitivity] {ps.spec_id}: {ps.description}")
        spec_rows.append({
            "spec_id": ps.spec_id,
            "expr_mode": ps.expr_mode,
            "target_norm": ps.target_norm,
            "aggregate": ps.aggregate,
            "description": ps.description,
        })

        gp = compute_gene_pressure(
            edges_m0,
            mirna,
            genes=list(hs.universe),
            expr_mode=ps.expr_mode,
            target_norm=ps.target_norm,
            aggregate=ps.aggregate,
        )
        shared = gp.columns.intersection(gate.index)
        gp_gated = gp[shared].mul(gate.reindex(shared), axis=1)
        hp = hallmark_pressure_matrix(gp_gated, hs)

        spec_dir = out_root / ps.spec_id
        spec_dir.mkdir(parents=True, exist_ok=True)
        gp.to_csv(spec_dir / "gene_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")
        hp.to_csv(spec_dir / "hallmark_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["spec_id"] = ps.spec_id
        coupling.to_csv(spec_dir / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t", index=False)
        coupling_all.append(coupling)

        summ = _coupling_summary(coupling, ps.spec_id)
        summary_all.append(summ)

        basal = coupling.loc[
            (coupling["subtype"] == "Basal") & (coupling["hallmark_set"].isin(KEY_HALLMARKS))
        ].copy()
        basal_key_all.append(basal)

        hub = gene_pressure_route_partial_corr(rna, clinical, proxies, gp_gated, cnv_hub)
        irf1 = hub.loc[(hub["target"] == "IRF1") & (hub["scope"] == "Basal")]
        if not irf1.empty:
            r = irf1.iloc[0]
            irf1_all.append({
                "spec_id": ps.spec_id,
                "irf1_basal_rho_cn_prolif": r["rho_e2f_g2m_CN"],
                "irf1_basal_p_cn_prolif": r["p_e2f_g2m_CN"],
                "irf1_basal_survives": r.get("survives_e2f_g2m_CN"),
            })

        b = summ.loc[summ["subtype"] == "Basal"]
        if not b.empty:
            print(f"  Basal: neg_sig={int(b.iloc[0]['n_neg_sig'])}/50, "
                  f"key_neg_sig={int(b.iloc[0]['n_key_neg_sig'])}/8, "
                  f"median_rho={b.iloc[0]['median_rho_prolif_cn_wsd_adj']:.3f}")

    pd.DataFrame(spec_rows).to_csv(out_root / "spec_manifest.tsv", sep="\t", index=False)
    full_coupling = pd.concat(coupling_all, ignore_index=True)
    full_coupling.to_csv(out_root / "coupling_by_spec_subtype.tsv", sep="\t", index=False)

    summary = pd.concat(summary_all, ignore_index=True)
    baseline = summary.loc[summary["spec_id"] == "softmax_z_degree_sum", ["subtype", "n_neg_sig", "n_key_neg_sig"]]
    if baseline.empty:
        baseline = summary.loc[summary["spec_id"] == "z_sum_baseline", ["subtype", "n_neg_sig", "n_key_neg_sig"]]
    baseline = baseline.rename(columns={"n_neg_sig": "baseline_n_neg_sig", "n_key_neg_sig": "baseline_n_key_neg_sig"})
    summary = summary.merge(baseline, on="subtype", how="left")
    summary["delta_n_neg_sig_vs_baseline"] = summary["n_neg_sig"] - summary["baseline_n_neg_sig"]
    summary["delta_n_key_neg_sig_vs_baseline"] = summary["n_key_neg_sig"] - summary["baseline_n_key_neg_sig"]
    summary.to_csv(out_root / "coupling_summary_by_spec.tsv", sep="\t", index=False)

    pd.concat(basal_key_all, ignore_index=True).to_csv(
        out_root / "basal_key_programs_by_spec.tsv", sep="\t", index=False,
    )
    if irf1_all:
        pd.DataFrame(irf1_all).to_csv(out_root / "irf1_basal_by_spec.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.pressure_sensitivity",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "edge_universe": "M0 miRTarBase evidence≥2",
        "proliferation_proxy": "e2f_g2m",
        "coupling_covariates": "CPE + HRD + proliferation + mean target CN + TCGA sequencing batch (plate_both, via _partial_batch); within-subtype z pressure & expression",
        "specs": [s.spec_id for s in specs],
        "outputs": [
            "spec_manifest.tsv",
            "coupling_by_spec_subtype.tsv",
            "coupling_summary_by_spec.tsv",
            "basal_key_programs_by_spec.tsv",
            "irf1_basal_by_spec.tsv",
            "{spec_id}/gene_pressure_per_sample.tsv.gz",
            "{spec_id}/hallmark_coupling_by_pam50_prolif.tsv",
        ],
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[pressure_sensitivity] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
