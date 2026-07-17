"""Three-layer pressure reporting + positive/absolute coupling sensitivity.

Layers (§6 methodology):
  L1 — coupling spine: ``softmax_z_logrpm`` + ``evidence_mass``, signed sum
  L2 — attribution proxy: ``softmax_logrpm`` + ``evidence_mass`` (no z; co-abundance read)
  L3 — acquisition: ``softmax_devhealthy_logrpm`` + ``evidence_mass`` (GTEx-anchored dev)
  L4a — spine + positive contributions only (``contrib_transform=pos``)
  L4b — spine + |contribution| sum (``contrib_transform=abs``)

Each layer rebuilds gene → Hallmark pressure (AGO-gated) and replays proliferation-
adjusted subtype coupling (same ladder as ``pressure_sensitivity``).

Also emits hub-gene **cancellation** diagnostics on the L1 signed spine.

Run:
    .venv/bin/python3 -m mirna_hallmark.pressure_layer_comparison
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Literal, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.healthy_anchor import gtex_qn_baseline
from mirna_hallmark.hybrid_pressure import HybridPressureSpec, build_hybrid_edges
from mirna_hallmark.pressure_engine import (
    AggregateMode,
    ContribTransform,
    ExprMode,
    TargetNormMode,
    compute_gene_pressure,
    compute_gene_pressure_contributions,
)
from mirna_hallmark.robustness_checks import (
    HUB_ROUTES,
    KEY_HALLMARKS,
    PAM50_SUBTYPES,
    _proliferation_proxies,
    hallmark_coupling_by_subtype,
)

OUT_DIR = C.OUTPUT_ROOT / "pressure_layer_comparison"

LayerId = Literal["L1_coupling_spine", "L2_attribution_noz", "L3_acquisition_devhealthy",
                  "L4a_posclip_spine", "L4b_abs_spine"]


@dataclass(frozen=True)
class LayerSpec:
    layer_id: LayerId
    expr_mode: ExprMode
    target_norm: TargetNormMode
    aggregate: AggregateMode
    contrib_transform: ContribTransform
    description: str
    needs_healthy_baseline: bool = False


LAYER_SPECS: List[LayerSpec] = [
    LayerSpec(
        "L1_coupling_spine",
        "softmax_z_logrpm",
        "evidence_mass",
        "sum",
        "signed",
        "Coupling spine: share × z × logRPM, signed sum",
    ),
    LayerSpec(
        "L2_attribution_noz",
        "softmax_logrpm",
        "evidence_mass",
        "sum",
        "signed",
        "Attribution proxy: share × logRPM, no z (non-negative contributions)",
    ),
    LayerSpec(
        "L3_acquisition_devhealthy",
        "softmax_devhealthy_logrpm",
        "evidence_mass",
        "sum",
        "signed",
        "Acquisition: share × (x − healthy) × logRPM",
        needs_healthy_baseline=True,
    ),
    LayerSpec(
        "L4a_posclip_spine",
        "softmax_z_logrpm",
        "evidence_mass",
        "sum",
        "pos",
        "Spine with Σ max(contribution, 0) — active push only",
    ),
    LayerSpec(
        "L4b_abs_spine",
        "softmax_z_logrpm",
        "evidence_mass",
        "sum",
        "abs",
        "Spine with Σ |contribution| — total activity regardless of sign",
    ),
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


def _coupling_summary(coupling: pd.DataFrame, layer_id: str) -> pd.DataFrame:
    rows = []
    for sub in PAM50_SUBTYPES:
        g = coupling.loc[coupling["subtype"] == sub]
        if g.empty:
            continue
        neg = g.loc[(g["rho_prolif_cn_wsd_adj"] < 0) & (g["q_prolif_cn_wsd_adj"] < 0.10)]
        key = g.loc[g["key_hallmark"] == True]  # noqa: E712
        neg_key = neg.loc[neg["key_hallmark"] == True]  # noqa: E712
        rows.append({
            "layer_id": layer_id,
            "subtype": sub,
            "n_neg_sig": int(len(neg)),
            "n_hallmarks": int(len(g)),
            "n_key_neg_sig": int(len(neg_key)),
            "n_key_hallmarks": int(len(key)),
            "median_rho_prolif_cn_wsd_adj": float(
                pd.to_numeric(g["rho_prolif_cn_wsd_adj"], errors="coerce").median()
            ),
        })
    return pd.DataFrame(rows)


def _hub_cancellation(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    hub_genes: Sequence[str],
) -> pd.DataFrame:
    """Signed-spine cancellation metrics for hub targets."""
    contrib_df = compute_gene_pressure_contributions(
        edges,
        mirna,
        genes=list(hub_genes),
        expr_mode="softmax_z_logrpm",
        target_norm="evidence_mass",
    )
    if contrib_df.empty:
        return pd.DataFrame()

    rows = []
    for gene in hub_genes:
        g = contrib_df.loc[contrib_df["gene"] == gene]
        if g.empty:
            continue
        signed = pd.to_numeric(g["global_signed_share"], errors="coerce")
        abs_sh = pd.to_numeric(g["global_abs_share"], errors="coerce")
        pos_sh = pd.to_numeric(g["mean_pos_share"], errors="coerce")
        # share-stack coherence: signed vs abs divergence
        stack = g[["miRNA", "global_signed_share", "global_abs_share", "mean_pos_share"]].copy()
        stack["share_divergence"] = (
            pd.to_numeric(stack["global_abs_share"], errors="coerce")
            - pd.to_numeric(stack["global_signed_share"], errors="coerce").abs()
        )
        top_cancel = stack.sort_values("share_divergence", ascending=False).head(3)
        rows.append({
            "gene": gene,
            "n_regulators": int(len(g)),
            "sum_global_abs_share": float(abs_sh.sum(skipna=True)),
            "sum_global_signed_share": float(signed.sum(skipna=True)),
            "signed_abs_divergence": float(abs_sh.sum(skipna=True) - abs(signed.sum(skipna=True))),
            "mean_pos_vs_abs_gap": float((abs_sh - pos_sh).mean(skipna=True)),
            "top_cancelling_arms": "; ".join(
                f"{r.miRNA}({r.share_divergence:.3f})" for r in top_cancel.itertuples()
            ),
        })
    return pd.DataFrame(rows)


def run(*, out_dir: Path | None = None, specs: Sequence[LayerSpec] = LAYER_SPECS) -> None:
    out_root = Path(out_dir or OUT_DIR)
    out_root.mkdir(parents=True, exist_ok=True)

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    hallmark_edges = D.load_hallmark_edges()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted(HUB_ROUTES.keys())
    cnv = D.load_cnv_target_genes(sorted({g for gs in hs.sets.values() for g in gs}))
    edges_m0 = _m0_edges(hs.universe, hallmark_edges)

    healthy_baseline: Optional[pd.Series] = None
    if any(s.needs_healthy_baseline for s in specs):
        healthy_baseline = gtex_qn_baseline()

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

    for ls in specs:
        print(f"[pressure_layer_comparison] {ls.layer_id}: {ls.description}")
        spec_rows.append({
            "layer_id": ls.layer_id,
            "expr_mode": ls.expr_mode,
            "target_norm": ls.target_norm,
            "aggregate": ls.aggregate,
            "contrib_transform": ls.contrib_transform,
            "description": ls.description,
        })

        kw: dict = {}
        if ls.needs_healthy_baseline:
            kw["healthy_baseline"] = healthy_baseline

        gp = compute_gene_pressure(
            edges_m0,
            mirna,
            genes=list(hs.universe),
            expr_mode=ls.expr_mode,
            target_norm=ls.target_norm,
            aggregate=ls.aggregate,
            contrib_transform=ls.contrib_transform,
            **kw,
        )
        shared = gp.columns.intersection(gate.index)
        gp_gated = gp[shared].mul(gate.reindex(shared), axis=1)
        hp = hallmark_pressure_matrix(gp_gated, hs)

        layer_dir = out_root / ls.layer_id
        layer_dir.mkdir(parents=True, exist_ok=True)
        gp.to_csv(layer_dir / "gene_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")
        hp.to_csv(layer_dir / "hallmark_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["layer_id"] = ls.layer_id
        coupling.to_csv(layer_dir / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t", index=False)
        coupling_all.append(coupling)
        summary_all.append(_coupling_summary(coupling, ls.layer_id))

        basal = coupling.loc[
            (coupling["subtype"] == "Basal") & (coupling["hallmark_set"].isin(KEY_HALLMARKS))
        ].copy()
        basal_key_all.append(basal)

        b = summary_all[-1].loc[summary_all[-1]["subtype"] == "Basal"]
        if not b.empty:
            print(
                f"  Basal: neg_sig={int(b.iloc[0]['n_neg_sig'])}/50, "
                f"key_neg_sig={int(b.iloc[0]['n_key_neg_sig'])}/8, "
                f"median_rho={b.iloc[0]['median_rho_prolif_cn_wsd_adj']:.3f}"
            )

    pd.DataFrame(spec_rows).to_csv(out_root / "layer_manifest.tsv", sep="\t", index=False)
    pd.concat(coupling_all, ignore_index=True).to_csv(
        out_root / "coupling_by_layer_subtype.tsv", sep="\t", index=False,
    )
    summary = pd.concat(summary_all, ignore_index=True)
    summary.to_csv(out_root / "coupling_summary_by_layer.tsv", sep="\t", index=False)

    basal_wide = pd.concat(basal_key_all, ignore_index=True)
    basal_wide.to_csv(out_root / "basal_key_programs_by_layer.tsv", sep="\t", index=False)

    # Headline pivot: key Hallmarks × layer @ Basal
    pivot_rows = []
    for h in KEY_HALLMARKS:
        row: dict = {"hallmark_set": h}
        for ls in specs:
            sub = basal_wide.loc[
                (basal_wide["hallmark_set"] == h) & (basal_wide["layer_id"] == ls.layer_id)
            ]
            if sub.empty:
                row[f"{ls.layer_id}_rho"] = np.nan
                row[f"{ls.layer_id}_neg_sig"] = False
            else:
                r = sub.iloc[0]
                row[f"{ls.layer_id}_rho"] = float(r["rho_prolif_cn_wsd_adj"])
                row[f"{ls.layer_id}_neg_sig"] = bool(
                    (r["rho_prolif_cn_wsd_adj"] < 0) and (r["q_prolif_cn_wsd_adj"] < 0.10)
                )
        pivot_rows.append(row)
    pd.DataFrame(pivot_rows).to_csv(out_root / "basal_key_pivot.tsv", sep="\t", index=False)

    cancel = _hub_cancellation(edges_m0, mirna, hub_genes)
    cancel.to_csv(out_root / "hub_cancellation_L1_spine.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.pressure_layer_comparison",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "edge_universe": "M0 miRTarBase evidence≥2",
        "layers": [s.layer_id for s in specs],
        "outputs": [
            "layer_manifest.tsv",
            "coupling_summary_by_layer.tsv",
            "basal_key_pivot.tsv",
            "hub_cancellation_L1_spine.tsv",
            "{layer_id}/hallmark_coupling_by_pam50_prolif.tsv",
        ],
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[pressure_layer_comparison] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
