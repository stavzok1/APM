"""Per-gene and per-miRNA resolution for actual PAM50 tumor subtypes.

This module answers a different question from the visibility-archetype analysis:
for the actual PAM50 tumor cohorts (LumA, LumB, Her2, Basal; Normal excluded),
which target genes are differentially miRNA-pressured, and which miRNAs drive
that pressure?

Outputs (``output/pam50_gene_resolution/``):

- ``gene_pressure_by_pam50.tsv``
  One row per target gene with raw and per-gene-z pressure means for LumA/LumB/
  Her2/Basal, target-expression means, and biologically useful contrasts.

- ``gene_pressure_one_vs_rest.tsv``
  One row per (subtype, gene): subtype vs the complement of the other three
  tumor PAM50 groups (Normal excluded). Includes pressure delta and target RNA
  delta.

- ``gene_mirna_drivers_by_pam50.tsv``
  Top high-evidence driver miRNAs per gene x subtype, using realized load:
  structural targeting weight x subtype mean miRNA z.

The key effect size for cross-gene comparison is ``zpressure_delta_*`` because
raw pressure magnitude scales with miRNA in-degree. Raw pressure is retained for
hub-gene biology (PTEN, BCL2L11, etc.).
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.hallmark_interaction import compute_gene_pressure


TUMOR_PAM50 = ("LumA", "LumB", "Her2", "Basal")


def subtype_members(clinical: pd.DataFrame) -> Dict[str, List[str]]:
    pam = clinical.set_index("participant")["PAM50_final"]
    return {s: pam.index[pam.eq(s)].tolist() for s in TUMOR_PAM50}


def _mean_by_samples(df: pd.DataFrame, samples: List[str]) -> pd.Series:
    cols = [c for c in samples if c in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    return df[cols].mean(axis=1)


def gene_pressure_by_pam50(
    gp_gated: pd.DataFrame,
    rna: pd.DataFrame,
    hs: HallmarkSets,
    members: Dict[str, List[str]],
) -> pd.DataFrame:
    """One row per gene, with raw/z pressure and target RNA per PAM50 subtype."""
    gpz = S.zscore_rows(gp_gated)
    rows = pd.DataFrame(index=gp_gated.index)
    rows["hallmark_sets"] = rows.index.map(lambda g: ";".join(hs.gene_to_sets.get(g, [])))

    for subtype, samples in members.items():
        rows[f"pressure_{subtype}"] = _mean_by_samples(gp_gated, samples)
        rows[f"zpressure_{subtype}"] = _mean_by_samples(gpz, samples).round(4)
        rows[f"target_log2tpm_{subtype}"] = _mean_by_samples(rna.reindex(gp_gated.index), samples).round(4)

    # Common biology contrasts. Normal-like is excluded everywhere here.
    rows["zpressure_Basal_minus_luminal"] = (
        rows["zpressure_Basal"] - rows[["zpressure_LumA", "zpressure_LumB"]].mean(axis=1)
    ).round(4)
    rows["zpressure_Her2_minus_luminal"] = (
        rows["zpressure_Her2"] - rows[["zpressure_LumA", "zpressure_LumB"]].mean(axis=1)
    ).round(4)
    rows["zpressure_LumB_minus_LumA"] = (rows["zpressure_LumB"] - rows["zpressure_LumA"]).round(4)
    rows["zpressure_Basal_minus_Her2"] = (rows["zpressure_Basal"] - rows["zpressure_Her2"]).round(4)

    rows["target_Basal_minus_luminal"] = (
        rows["target_log2tpm_Basal"] - rows[["target_log2tpm_LumA", "target_log2tpm_LumB"]].mean(axis=1)
    ).round(4)
    rows["target_Her2_minus_luminal"] = (
        rows["target_log2tpm_Her2"] - rows[["target_log2tpm_LumA", "target_log2tpm_LumB"]].mean(axis=1)
    ).round(4)
    return rows.rename_axis("gene").reset_index()


def one_vs_rest(
    gp_gated: pd.DataFrame,
    rna: pd.DataFrame,
    hs: HallmarkSets,
    members: Dict[str, List[str]],
) -> pd.DataFrame:
    """Subtype vs complement of the other three tumor PAM50 groups (Normal excluded)."""
    gpz = S.zscore_rows(gp_gated)
    rna_sub = rna.reindex(gp_gated.index)
    rows = []
    for subtype, samples in members.items():
        comp = [p for s, ids in members.items() if s != subtype for p in ids]
        p_sub = _mean_by_samples(gp_gated, samples)
        p_comp = _mean_by_samples(gp_gated, comp)
        z_sub = _mean_by_samples(gpz, samples)
        z_comp = _mean_by_samples(gpz, comp)
        r_sub = _mean_by_samples(rna_sub, samples)
        r_comp = _mean_by_samples(rna_sub, comp)
        for gene in gp_gated.index:
            rows.append(
                {
                    "subtype": subtype,
                    "gene": gene,
                    "hallmark_sets": ";".join(hs.gene_to_sets.get(gene, [])),
                    "n_subtype": len([p for p in samples if p in gp_gated.columns]),
                    "n_complement": len([p for p in comp if p in gp_gated.columns]),
                    "pressure_subtype": round(float(p_sub.loc[gene]), 4),
                    "pressure_complement": round(float(p_comp.loc[gene]), 4),
                    "pressure_delta_subtype_minus_complement": round(float(p_sub.loc[gene] - p_comp.loc[gene]), 4),
                    "zpressure_subtype": round(float(z_sub.loc[gene]), 4),
                    "zpressure_complement": round(float(z_comp.loc[gene]), 4),
                    "zpressure_delta_subtype_minus_complement": round(float(z_sub.loc[gene] - z_comp.loc[gene]), 4),
                    "target_log2tpm_subtype": round(float(r_sub.loc[gene]), 4),
                    "target_log2tpm_complement": round(float(r_comp.loc[gene]), 4),
                    "target_delta_subtype_minus_complement": round(float(r_sub.loc[gene] - r_comp.loc[gene]), 4),
                }
            )
    return pd.DataFrame(rows)


def gene_mirna_drivers_by_pam50(
    edges: pd.DataFrame,
    mirna_expr: pd.DataFrame,
    members: Dict[str, List[str]],
    *,
    top_n: int = 5,
) -> pd.DataFrame:
    """Top high-evidence driver miRNAs per gene x subtype by realized load."""
    he = edges.loc[edges["high_evidence"]].copy()
    he["w"] = np.log1p(pd.to_numeric(he["evidence_score"], errors="coerce").fillna(0.0))
    load = he.groupby(["gene", "miRNA"], as_index=False).agg(
        structural_weight=("w", "sum"),
        n_hallmark_contexts=("hallmark_set", "nunique"),
    )
    z = S.zscore_rows(mirna_expr).fillna(0.0)
    load = load.loc[load["miRNA"].isin(z.index)].copy()

    parts = []
    for subtype, samples in members.items():
        cols = [c for c in samples if c in z.columns]
        mean_z = z[cols].mean(axis=1) if cols else pd.Series(0.0, index=z.index)
        tmp = load.copy()
        tmp["subtype"] = subtype
        tmp["mean_z_mirna"] = tmp["miRNA"].map(mean_z).round(4)
        tmp["realized_load"] = (tmp["structural_weight"] * tmp["mean_z_mirna"]).round(4)
        tmp["abs_realized_load"] = tmp["realized_load"].abs()
        tmp = (
            tmp.sort_values(["gene", "abs_realized_load"], ascending=[True, False])
            .groupby("gene", as_index=False)
            .head(top_n)
        )
        parts.append(tmp)
    out = pd.concat(parts, ignore_index=True)
    out["structural_weight"] = out["structural_weight"].round(4)
    return out.drop(columns=["abs_realized_load"])


def run(*, out_dir: Path = C.PAM50_GENE_RESOLUTION_DIR, include_tnrc6: Optional[bool] = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    members = subtype_members(clinical)
    print("[pam50_gene] group sizes: " + ", ".join(f"{k}={len(v)}" for k, v in members.items()))

    # Some gene symbols can appear more than once in the processed RNA matrix
    # (aliases/paralogs). Collapse before reindexing to avoid duplicate-axis
    # failures; mean is consistent with the parent loaders' sample-collapse style.
    rna = D.load_rna().groupby(level=0).mean()
    edges = D.load_hallmark_edges()
    mirna_expr = D.load_mirna_arms()

    gp = compute_gene_pressure(hs.universe)
    from dataclasses import replace

    params = C.AGO_GATE if include_tnrc6 is None else replace(C.AGO_GATE, include_tnrc6=include_tnrc6)
    gate = compute_ago_gate(rna, params=params)["ago_gate"]
    shared = gp.columns.intersection(gate.index)
    gp_gated = gp[shared].mul(gate.reindex(shared), axis=1)

    gene_table = gene_pressure_by_pam50(gp_gated, rna, hs, members)
    gene_table.to_csv(out_dir / "gene_pressure_by_pam50.tsv", sep="\t", index=False)

    ovr = one_vs_rest(gp_gated, rna, hs, members)
    ovr.to_csv(out_dir / "gene_pressure_one_vs_rest.tsv", sep="\t", index=False)

    drivers = gene_mirna_drivers_by_pam50(edges, mirna_expr, members)
    drivers.to_csv(out_dir / "gene_mirna_drivers_by_pam50.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.pam50_gene_resolution",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "subtypes": list(TUMOR_PAM50),
        "excluded": ["Normal"],
        "group_sizes": {k: len(v) for k, v in members.items()},
        "effect_size_note": "Use zpressure deltas for cross-gene comparison; raw pressure retained for hub genes.",
        "outputs": {
            "gene_pressure_by_pam50": "raw/z pressure + target RNA for LumA/LumB/Her2/Basal",
            "gene_pressure_one_vs_rest": "subtype vs complement (other tumor PAM50 groups)",
            "gene_mirna_drivers_by_pam50": "top high-evidence miRNAs per gene x subtype by realized load",
        },
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[pam50_gene] wrote outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include-tnrc6", action="store_true", default=None)
    ap.add_argument("--out-dir", type=Path, default=C.PAM50_GENE_RESOLUTION_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, include_tnrc6=args.include_tnrc6)


if __name__ == "__main__":
    main()
