"""Evidence mass + spec-in-softmax sensitivity on the L1 spine.

Compares edge-weight denominators and putting log1p(evidence) into softmax logits:

  baseline      edge_w = log1p(ev) / log1p(Σ log1p ev)     [current evidence_mass]
  inner_sum_log edge_w = log1p(ev) / Σ log1p ev             [inner log only, no outer log]
  outer_log     edge_w = log1p(ev) / log1p(Σ ev)           [outer log on raw sum]
  log_only      edge_w = log1p(ev)                          [no mass denominator]
  spec_unity    logit = (x−med) + log1p(ev); edge_w = 1
  spec_degree   logit = (x−med) + log1p(ev); edge_w = 1/log1p(degree)

All use expr_mode softmax_z_logrpm (or softmax_z_logrpm_spec for spec rows).

Run:
    .venv/bin/python3 -m mirna_hallmark.pressure_evidence_sensitivity
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Sequence

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
    _proliferation_proxies,
    hallmark_coupling_by_subtype,
)

OUT_DIR = C.OUTPUT_ROOT / "pressure_evidence_sensitivity"


@dataclass(frozen=True)
class EvidenceSpec:
    spec_id: str
    expr_mode: ExprMode
    target_norm: TargetNormMode
    description: str


EVIDENCE_SPECS: List[EvidenceSpec] = [
    EvidenceSpec(
        "baseline_evidence_mass",
        "softmax_z_logrpm",
        "evidence_mass",
        "log1p(ev) / log1p(Σ log1p ev) — current spine",
    ),
    EvidenceSpec(
        "inner_sum_log",
        "softmax_z_logrpm",
        "evidence_sum_log",
        "log1p(ev) / Σ log1p ev — inner log sum, no outer log1p",
    ),
    EvidenceSpec(
        "outer_log_mass",
        "softmax_z_logrpm",
        "evidence_outer_log",
        "log1p(ev) / log1p(Σ ev) — outer log on raw evidence sum",
    ),
    EvidenceSpec(
        "log_only_no_denom",
        "softmax_z_logrpm",
        "evidence_log_only",
        "log1p(ev) only — no promiscuity denominator",
    ),
    EvidenceSpec(
        "spec_in_softmax_unity",
        "softmax_z_logrpm_spec",
        "unity",
        "softmax logit = (x−med)+log1p(ev); edge_w=1",
    ),
    EvidenceSpec(
        "spec_in_softmax_degree",
        "softmax_z_logrpm_spec",
        "degree_only",
        "softmax logit = (x−med)+log1p(ev); edge_w=1/log1p(degree)",
    ),
]


def _m0_edges(genes: Sequence[str], hallmark_edges: pd.DataFrame) -> pd.DataFrame:
    spec = HybridPressureSpec(mode="M0")
    e = build_hybrid_edges(genes, spec=spec, hallmark_edges=hallmark_edges)
    return e.rename(columns={"edge_weight": "evidence_score"})[["miRNA", "gene", "evidence_score"]]


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
            "n_key_neg_sig": int(len(neg_key)),
            "median_rho_prolif_cn_wsd_adj": float(
                pd.to_numeric(g["rho_prolif_cn_wsd_adj"], errors="coerce").median()
            ),
        })
    return pd.DataFrame(rows)


def run(*, out_dir: Path | None = None, specs: Sequence[EvidenceSpec] = EVIDENCE_SPECS) -> None:
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

    for es in specs:
        print(f"[pressure_evidence_sensitivity] {es.spec_id}: {es.description}")
        spec_rows.append({
            "spec_id": es.spec_id,
            "expr_mode": es.expr_mode,
            "target_norm": es.target_norm,
            "aggregate": "sum",
            "description": es.description,
        })

        gp = compute_gene_pressure(
            edges_m0,
            mirna,
            genes=list(hs.universe),
            expr_mode=es.expr_mode,
            target_norm=es.target_norm,
            aggregate="sum",
        )
        shared = gp.columns.intersection(gate.index)
        hp = hallmark_pressure_matrix(gp[shared].mul(gate.reindex(shared), axis=1), hs)

        spec_dir = out_root / es.spec_id
        spec_dir.mkdir(parents=True, exist_ok=True)
        coupling = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
        coupling["spec_id"] = es.spec_id
        coupling.to_csv(spec_dir / "hallmark_coupling_by_pam50_prolif.tsv", sep="\t", index=False)
        coupling_all.append(coupling)
        summary_all.append(_coupling_summary(coupling, es.spec_id))

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

    pd.DataFrame(spec_rows).to_csv(out_root / "spec_manifest.tsv", sep="\t", index=False)
    pd.concat(coupling_all, ignore_index=True).to_csv(
        out_root / "coupling_by_spec_subtype.tsv", sep="\t", index=False,
    )
    summary = pd.concat(summary_all, ignore_index=True)
    summary.to_csv(out_root / "coupling_summary_by_spec.tsv", sep="\t", index=False)

    basal_wide = pd.concat(basal_key_all, ignore_index=True)
    pivot_rows = []
    for h in KEY_HALLMARKS:
        row: dict = {"hallmark_set": h}
        for es in specs:
            sub = basal_wide.loc[
                (basal_wide["hallmark_set"] == h) & (basal_wide["spec_id"] == es.spec_id)
            ]
            if sub.empty:
                row[f"{es.spec_id}_rho"] = float("nan")
                row[f"{es.spec_id}_neg_sig"] = False
            else:
                r = sub.iloc[0]
                row[f"{es.spec_id}_rho"] = float(r["rho_prolif_cn_wsd_adj"])
                row[f"{es.spec_id}_neg_sig"] = bool(
                    (r["rho_prolif_cn_wsd_adj"] < 0) and (r["q_prolif_cn_wsd_adj"] < 0.10)
                )
        pivot_rows.append(row)
    pd.DataFrame(pivot_rows).to_csv(out_root / "basal_key_pivot.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.pressure_evidence_sensitivity",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "specs": [s.spec_id for s in specs],
    }
    (out_root / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[pressure_evidence_sensitivity] done -> {out_root}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=None)
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
