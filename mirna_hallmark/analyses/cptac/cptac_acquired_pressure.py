"""Acquired-pressure (tumor − NAT) version of the gene-centric CPTAC validation.

The spine ``cptac_validation`` correlates a tumor's **realized** gene pressure with
its mRNA / protein. This module instead uses **acquired** pressure — the pressure
the tumor gained over its *own* matched normal tissue:

    acq(g, p) = pressure_tumor(g, p) − pressure_nat(g, p)        (paired, per participant)

Why paired (and not a population normal baseline): the gene-level / pooled /
Hallmark correlations are all invariant to subtracting a per-gene **constant**
(Spearman/Pearson are shift- and scale-invariant per gene, and the pooled test
within-gene-centers). A population normal median is exactly such a per-gene
constant → it would reproduce the spine bit-for-bit. Only a **per-participant**
normal (matched NAT miRNA) changes the across-sample correlation. So this analysis
is restricted to participants with both tumor and NAT miRNA.

Pressure is built **cross-state-comparable** (``softmax_logrpm``, no within-cohort
z-score, identical edge set in both states) so the tumor−NAT difference is on a
common scale — the same attribution mode as ``edge_acquired_pressure_panels``.
Two variants: all-edge spine and high-evidence-only, both **ungated** (the AGO gate
needs per-state RISC/RNA which is not defined for the NAT side here).

Cohort: **tcga105** only. The matched-NAT TCGA miRNA participants (n≈101) give a
well-powered acquired-pressure → **TCGA mRNA** test; the subset that also has CPTAC
protein (n≈14) is **underpowered for gene-level** (``correlation_pair`` floors at
n<15 → NaN) but the **pooled** within-gene-centered test (D63-style) still resolves
across gene×participant. Read the protein rows as descriptive only.

Layers (response ~ acquired pressure):
  mrna_anticorr    TCGA tumor mRNA           expect negative   (powered, n≤101)
  protein_anticorr CPTAC protein z           expect negative   (n≈14, pooled only)
  protein_resid    protein beyond its mRNA   expect negative   (n≈14, pooled only)
  rna_protein_gap  cptac rna_z − protein_z   expect positive   (n≈14, pooled only)

Run::

    .venv/bin/python3 -m mirna_hallmark.cptac_acquired_pressure
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.eval import cptac_validation as V
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.analyses.edge_panels.edge_acquired_pressure_panels import _matched_pair_matrices
from mirna_hallmark.pressure_build import method_blurb, resolve_pressure_edges
from mirna_hallmark.pressure_engine import (
    compute_gene_pressure as engine_pressure,
    filter_edges_by_abundance,
)

OUT_DIR = V.OUT_DIR / "acquired"

# Cross-state-comparable attribution (matches edge_acquired_pressure_panels).
ACQ_EXPR_MODE = "softmax_logrpm"
ACQ_TARGET_NORM = "evidence_mass"

# (label, high-evidence-only) — both ungated.
ACQ_VARIANTS: Tuple[Tuple[str, bool], ...] = (
    ("all_edge|ungated", False),
    ("highev|ungated", True),
)

# Response layers for the acquired analysis (extends the spine LAYERS with mRNA).
ACQ_LAYERS: Tuple[Tuple[str, str, str], ...] = (
    ("mrna_anticorr", "rna_z", "negative"),         # acquired pressure -> lower tumor mRNA
    ("protein_anticorr", "protein_z", "negative"),  # -> lower protein
    ("protein_resid", "protein_resid", "negative"), # -> protein below its mRNA-predicted level
    ("rna_protein_gap", "gap", "positive"),         # -> bigger rna-protein gap
)


def _acquired_pressure_matrices(
    hs: HallmarkSets,
) -> Tuple[Dict[str, pd.DataFrame], List[str]]:
    """gene×participant acquired (tumor − NAT) pressure, keyed by variant label.

    Edges are resolved + abundance-filtered **once** (on the tumor matrix) so the
    edge universe is identical in both states; pressure is then computed ungated and
    z-free (``softmax_logrpm``) in each state and differenced on matched columns.
    """
    universe = list(hs.universe)
    tumor_mat, nat_mat, matched = _matched_pair_matrices()
    print(f"[acq] matched tumor+NAT participants: {len(matched)}")

    edges_all = D.load_hallmark_edges()
    edges_he = D.high_evidence_edges(edges_all)

    out: Dict[str, pd.DataFrame] = {}
    for label, high_ev in ACQ_VARIANTS:
        edges = edges_he if high_ev else edges_all
        e = resolve_pressure_edges(edges, tumor_mat)  # arm space shared by tumor/NAT
        e = filter_edges_by_abundance(e, tumor_mat, C.PRESSURE_ABUNDANCE_FLOOR)
        print(f"[acq] {label}: {len(e):,} edges -> pressure (tumor & NAT) ...")
        gp_t = engine_pressure(e, tumor_mat, genes=universe,
                               expr_mode=ACQ_EXPR_MODE, target_norm=ACQ_TARGET_NORM)
        gp_n = engine_pressure(e, nat_mat, genes=universe,
                               expr_mode=ACQ_EXPR_MODE, target_norm=ACQ_TARGET_NORM)
        cols = gp_t.columns.intersection(gp_n.columns)
        genes = gp_t.index.intersection(gp_n.index)
        out[label] = gp_t.loc[genes, cols] - gp_n.loc[genes, cols]
    return out, matched


def run(*, out_dir: Path = OUT_DIR) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()

    pressures, matched = _acquired_pressure_matrices(hs)

    # Responses: TCGA tumor mRNA (matched participants) + CPTAC tcga105 protein layers.
    rna = D.load_rna()
    rna_matched = rna[[c for c in matched if c in rna.columns]]
    # TCGA RNA has duplicate gene-symbol rows; collapse to one row per gene.
    if not rna_matched.index.is_unique:
        rna_matched = rna_matched.groupby(level=0).mean()
    prot_layers = V.load_cptac_layers("tcga105")  # protein_z, rna_z(cptac), gap, protein_resid
    layers = {
        "rna_z": rna_matched,                       # NOTE: TCGA tumor mRNA (not CPTAC rna)
        "protein_z": prot_layers["protein_z"],
        "protein_resid": prot_layers["protein_resid"],
        "gap": prot_layers["gap"],
    }

    clinical = D.load_clinical_strata()
    cov = V._covariates(clinical, "tcga105")
    enr = V._target_enrichment(hs)

    # Reuse the validated spine machinery with the extended layer set.
    saved = V.LAYERS
    V.LAYERS = ACQ_LAYERS
    try:
        print("[acq] gene-level associations ...")
        gene_tbl = V.gene_level_associations(pressures, layers, cov)
        print("[acq] pooled within-gene-centered ...")
        pooled_tbl = V.pooled_associations(pressures, layers, clinical)
        print("[acq] Hallmark-level (cohort + strata) ...")
        hm_tbl = V.hallmark_associations(pressures, layers, clinical, hs, cov, enr)
    finally:
        V.LAYERS = saved

    gene_tbl.to_csv(out_dir / "acquired_gene_level_associations.tsv.gz",
                    sep="\t", index=False, compression="gzip")
    pooled_tbl.to_csv(out_dir / "acquired_pooled_associations.tsv", sep="\t", index=False)
    hm_tbl.to_csv(out_dir / "acquired_hallmark_associations.tsv", sep="\t", index=False)

    n_mrna = int(len(set(rna_matched.columns) & set(next(iter(pressures.values())).columns)))
    n_prot = int(len(set(prot_layers["protein_z"].columns)
                     & set(next(iter(pressures.values())).columns)))
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "module": "mirna_hallmark.cptac_acquired_pressure",
        "cohort": "tcga105 (matched tumor+NAT miRNA participants)",
        "pressure": (
            f"acquired = tumor − NAT, paired per participant; {ACQ_EXPR_MODE} / "
            f"{ACQ_TARGET_NORM}; ungated; identical edge set across states; "
            f"base method: {method_blurb()}"
        ),
        "variants": [lbl for lbl, _ in ACQ_VARIANTS],
        "n_participants_acquired_and_mrna": n_mrna,
        "n_participants_acquired_and_protein": n_prot,
        "layers": {name: f"{key} ~ acquired (expect {sign})" for name, key, sign in ACQ_LAYERS},
        "confounders": "CPE + thornsson_hrd_score (partial Spearman/Pearson; raw reported)",
        "power_note": (
            "mrna_anticorr is well-powered (n<=101). protein layers overlap only "
            f"n={n_prot} participants — gene-level is NaN by the n>=15 floor in "
            "correlation_pair; read the pooled protein rows (within-gene-centered, "
            "D63-style) as descriptive only."
        ),
        "outputs": [
            "acquired_gene_level_associations.tsv.gz",
            "acquired_pooled_associations.tsv",
            "acquired_hallmark_associations.tsv",
        ],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n",
                                                  encoding="utf-8")

    print(f"[acq] wrote {out_dir} (mRNA n={n_mrna}, protein n={n_prot})")
    cohort_pool = pooled_tbl[pooled_tbl["stratum_group"] == "cohort"]
    for _, r in cohort_pool.iterrows():
        print(f"  {r['variant']:<18} {r['layer']:<16} pooled ρ={r['pooled_spearman']:+.3f} "
              f"(p={r['pooled_p']:.1e}, n_gp={r['n_gene_participant']})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
