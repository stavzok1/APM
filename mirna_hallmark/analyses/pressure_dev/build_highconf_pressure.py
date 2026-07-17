"""High-confidence (evidence + TargetScan) Hallmark pressure and subtype coupling.

Motivation (see chat 2026-06): the default whole-miRNA Hallmark pressure sums
~515 miRNAs per program, and the bulk of low-confidence edges contribute a noisy,
near-zero-coupling component that *dilutes* the real repressive signal — worst in
Her2, where the full-aggregate coupling collapses to ~0 even though the curated
miR-17~92/106b~25 cluster couples at rho=-0.37. Restricting the edge universe to
**high-evidence AND TargetScan-predicted** (miRNA, gene) edges is a non-circular,
sequence-based prune that:
  - cuts the miRNA universe ~515 -> ~197 (edges 2447 -> 1077 over the 8 key genes),
  - turns the non-cluster bulk weakly coherent (Her2 +0.02 -> -0.02), and
  - roughly doubles the FULL-aggregate Her2 coupling (-0.05 -> -0.12) while
    preserving the Basal-maximal ordering.

This module rebuilds the Hallmark pressure on that universe and re-runs the exact
subtype coupling used for Fig. 1 (``hallmark_coupling_by_subtype``), writing
``hallmark_coupling_by_pam50_highconf.tsv`` next to the default table. Fig. 2
(cluster) is unchanged.

    .venv/bin/python3 -m mirna_hallmark.build_highconf_pressure
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import robustness_checks as R
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.pressure_build import compute_gene_pressure

OUT_TABLE = R.AIM1_DIR / "hallmark_coupling_by_pam50_highconf.tsv"
OUT_PRESSURE = C.OUTPUT_ROOT / "hallmark_interaction" / "hallmark_pressure_per_sample_highconf.tsv.gz"


def highconf_edges(hs: HallmarkSets) -> pd.DataFrame:
    """High-evidence AND TargetScan-predicted (miRNA, gene) edges for pressure."""
    raw = D.load_hallmark_edges()
    he = raw.loc[raw["high_evidence"] == True].copy()  # noqa: E712
    ts = R._load_targetscan_weights(list(hs.universe))
    ts_pairs = set(map(tuple, ts[["miRNA", "gene"]].itertuples(index=False, name=None)))
    he["pair"] = list(zip(he["miRNA"], he["gene"]))
    he = he.loc[he["pair"].isin(ts_pairs)]
    edges = he[["miRNA", "gene", "evidence_score"]].drop_duplicates()
    print(f"[highconf] high-evidence+TargetScan edges: {len(edges)} "
          f"over {edges['miRNA'].nunique()} miRNAs / {edges['gene'].nunique()} genes")
    return edges


def build() -> Path:
    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    proxies = R._proliferation_proxies(rna, hs)
    hub_genes = sorted(set(R.HUB_ROUTES) | {g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(include_tnrc6=False, gate_min=C.AGO_GATE.gate_min,
                             gate_k=C.AGO_GATE.gate_k, gate_midpoint=C.AGO_GATE.gate_midpoint),
    )["ago_gate"]

    edges = highconf_edges(hs)
    gp = compute_gene_pressure(hs.universe, edges=edges)
    sh = gp.columns.intersection(gate.index)
    hp = hallmark_pressure_matrix(gp[sh].mul(gate.reindex(sh), axis=1), hs)
    OUT_PRESSURE.parent.mkdir(parents=True, exist_ok=True)
    hp.to_csv(OUT_PRESSURE, sep="\t", compression="gzip")

    coupling = R.hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp, cnv)
    coupling.to_csv(OUT_TABLE, sep="\t", index=False)
    print(f"[highconf] wrote {OUT_TABLE}")
    print(f"[highconf] wrote {OUT_PRESSURE}")
    key = coupling.loc[coupling["key_hallmark"] == True]  # noqa: E712
    piv = key.pivot(index="hallmark_set", columns="subtype", values="rho_prolif_cn_wsd_adj")
    print("[highconf] key-program coupling (rho_prolif_cn_wsd_adj):")
    print(piv[["LumA", "LumB", "Her2", "Basal"]].round(3).to_string())
    return OUT_TABLE


if __name__ == "__main__":
    build()
