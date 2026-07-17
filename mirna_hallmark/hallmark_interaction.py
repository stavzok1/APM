"""Core analysis: AGO-gated miRNA pressure on Hallmark programs.

Pipeline:
1. **Per-(gene, sample) miRNA pressure** -- ``pressure_build.compute_gene_pressure``
   (default spine: ``softmax_z_logrpm`` × ``evidence_mass`` on M0 miRTarBase edges).
2. **AGO/RISC gate** -- multiply each sample's pressure by its bounded RISC
   capacity gate (``ago_gate.py``); gated AND ungated are carried throughout.
3. **Hallmark pressure** -- mean gated pressure over a Hallmark's member genes,
   per sample (hallmark x sample matrices).
4. **Anti-correlation** -- per Hallmark, Spearman/Pearson(pressure, mean target-gene
   expression) across samples; expected **negative**. Partial Spearman and partial
   Pearson adjust CPE + HRD. Reported cohort-wide and per PAM50 subtype (interaction view).
5. **High-evidence target enrichment** -- per Hallmark hypergeometric test of
   whether its genes are over-represented among high-evidence miRNA targets.

Hybrid extended pressure (M7/M8/M11): ``hybrid_pressure`` module, same weighting.

Outputs (``output/hallmark_interaction/``):
- ``hallmark_pressure_per_sample.tsv.gz``     -- gated hallmark x sample pressure
- ``hallmark_anticorrelation.tsv``            -- cohort gated/ungated rho,p,q
- ``hallmark_anticorrelation_by_pam50.tsv``   -- per-subtype anti-correlation
- ``hallmark_highev_target_enrichment.tsv``   -- hypergeometric enrichment
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import compute_gene_pressure, method_blurb


def hallmark_pressure_matrix(
    gene_pressure: pd.DataFrame,
    hs: HallmarkSets,
) -> pd.DataFrame:
    """Hallmark x sample = mean per-gene pressure over present member genes."""
    rows = {}
    for hset, members in hs.sets.items():
        present = [g for g in members if g in gene_pressure.index]
        if not present:
            continue
        rows[hset] = gene_pressure.loc[present].mean(axis=0)
    return pd.DataFrame(rows).T


def hallmark_expression_matrix(
    rna: pd.DataFrame,
    hs: HallmarkSets,
) -> pd.DataFrame:
    """Hallmark x sample = mean target-gene log2(TPM+1) over present member genes."""
    rows = {}
    for hset, members in hs.sets.items():
        present = [g for g in members if g in rna.index]
        if not present:
            continue
        rows[hset] = rna.loc[present].mean(axis=0)
    return pd.DataFrame(rows).T


def _anticorr(
    pressure: pd.DataFrame,
    expression: pd.DataFrame,
    clinical: pd.DataFrame,
    *,
    label: str,
) -> pd.DataFrame:
    """Per-Hallmark Spearman/Pearson(pressure, expression) + partial variants."""
    clin = clinical.set_index("participant")
    conf_num = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    samples = pressure.columns.intersection(expression.columns)

    rows = []
    for hset in pressure.index.intersection(expression.index):
        p = pressure.loc[hset, samples]
        e = expression.loc[hset, samples]
        df = pd.concat([p.rename("p"), e.rename("e")], axis=1).dropna()
        if len(df) < 20:
            continue
        cov = D.augment_tcga_batch(clin.reindex(df.index)[conf_num]) if conf_num else None
        st = S.correlation_pair(df["e"], df["p"], cov)
        rows.append({
            "hallmark_set": hset,
            "view": label,
            "n": st["n"],
            "spearman_rho": round(float(st["spearman_rho"]), 4),
            "spearman_p": st["spearman_p"],
            "pearson_r": round(float(st["pearson_r"]), 4),
            "pearson_p": st["pearson_p"],
            "partial_rho": (
                round(float(st["partial_rho"]), 4) if np.isfinite(st["partial_rho"]) else np.nan
            ),
            "partial_p": (
                float(st["partial_p"]) if np.isfinite(st["partial_p"]) else np.nan
            ),
            "partial_pearson_r": (
                round(float(st["partial_pearson_r"]), 4)
                if np.isfinite(st["partial_pearson_r"])
                else np.nan
            ),
            "partial_pearson_p": (
                float(st["partial_pearson_p"]) if np.isfinite(st["partial_pearson_p"]) else np.nan
            ),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["spearman_q"] = S.bh_fdr(out["spearman_p"].values)
        out["pearson_q"] = S.bh_fdr(out["pearson_p"].values)
        if out["partial_p"].notna().any():
            out["partial_q"] = S.bh_fdr(out["partial_p"].fillna(1.0).values)
        if out["partial_pearson_p"].notna().any():
            out["partial_pearson_q"] = S.bh_fdr(out["partial_pearson_p"].fillna(1.0).values)
        out = out.sort_values("spearman_rho")
    return out


def target_enrichment(edges: pd.DataFrame, hs: HallmarkSets) -> pd.DataFrame:
    """Per-Hallmark hypergeometric enrichment of high-evidence miRNA target genes."""
    universe = set(hs.universe)
    he_targets = set(edges.loc[edges["high_evidence"], "gene"]) & universe
    N, K = len(universe), len(he_targets)

    rows, pvals = [], []
    for hset, members in hs.sets.items():
        mem = set(members) & universe
        n = len(mem)
        k = len(mem & he_targets)
        fold, p = S.hypergeom_enrichment(k, n, K, N)
        rows.append({
            "hallmark_set": hset, "n_genes": n, "n_highev_target_genes": k,
            "frac_highev_targets": round(k / n, 4) if n else np.nan,
            "fold_enrichment": (round(fold, 3) if pd.notna(fold) else np.nan),
            "hypergeom_p": p,
        })
        pvals.append(p)
    out = pd.DataFrame(rows)
    out["hypergeom_q"] = S.bh_fdr(out["hypergeom_p"].values)
    return out.sort_values("fold_enrichment", ascending=False)


def run(
    *,
    out_dir: Path = C.INTERACTION_DIR,
    include_tnrc6: Optional[bool] = None,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    edges = D.load_hallmark_edges()

    print("[interaction] computing per-gene miRNA pressure ...")
    gene_pressure = compute_gene_pressure(hs.universe)
    print(f"[interaction]   pressure for {gene_pressure.shape[0]:,} genes x {gene_pressure.shape[1]:,} samples")

    print("[interaction] applying AGO/RISC gate ...")
    rna = D.load_rna()
    from dataclasses import replace
    params = C.AGO_GATE if include_tnrc6 is None else replace(C.AGO_GATE, include_tnrc6=include_tnrc6)
    gate_df = compute_ago_gate(rna, params=params)
    gate = gate_df["ago_gate"]
    shared = gene_pressure.columns.intersection(gate.index)
    gp_ungated = gene_pressure[shared]
    gp_gated = gp_ungated.mul(gate.reindex(shared), axis=1)

    hp_gated = hallmark_pressure_matrix(gp_gated, hs)
    hp_ungated = hallmark_pressure_matrix(gp_ungated, hs)
    hp_gated.to_csv(out_dir / "hallmark_pressure_per_sample.tsv.gz", sep="\t", compression="gzip")

    he = hallmark_expression_matrix(rna, hs)

    print("[interaction] anti-correlation (pressure vs target expression) ...")
    ac_gated = _anticorr(hp_gated, he, clinical, label="gated")
    ac_ungated = _anticorr(hp_ungated, he, clinical, label="ungated")
    ac = pd.concat([ac_gated, ac_ungated], ignore_index=True)
    ac.to_csv(out_dir / "hallmark_anticorrelation.tsv", sep="\t", index=False)
    neg_sig = ac_gated.loc[(ac_gated["spearman_q"] < C.FDR_ALPHA) & (ac_gated["spearman_rho"] < 0)]
    print(f"[interaction]   Hallmarks with gated anti-correlation (q<{C.FDR_ALPHA}, rho<0): "
          f"{len(neg_sig)}/{len(ac_gated)}")

    print("[interaction] per-PAM50 anti-correlation ...")
    pam = clinical.set_index("participant")["PAM50_final"]
    parts = []
    for subtype in C.PAM50_ORDER:
        keep = pam.index[pam == subtype]
        cols = hp_gated.columns.intersection(keep)
        if len(cols) < 20:
            continue
        sub = _anticorr(hp_gated[cols], he[he.columns.intersection(cols)], clinical, label=f"gated|{subtype}")
        sub["pam50"] = subtype
        parts.append(sub)
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(
            out_dir / "hallmark_anticorrelation_by_pam50.tsv", sep="\t", index=False
        )

    print("[interaction] high-evidence target enrichment ...")
    enr = target_enrichment(edges, hs)
    enr.to_csv(out_dir / "hallmark_highev_target_enrichment.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.hallmark_interaction",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "pressure": method_blurb(),
        "ago_gate": {"include_tnrc6": params.include_tnrc6, "model": "ago_gate.compute_ago_gate"},
        "hallmark_pressure": "mean per-gene gated pressure over member genes",
        "anticorrelation": (
            "Spearman + Pearson(hallmark pressure, mean target log2TPM); expect rho/r<0; "
            "partial adjusts CPE+HRD; FDR within view"
        ),
        "enrichment": "hypergeometric: Hallmark genes over-represented among high-evidence targets",
        "n_hallmarks_pressure": int(hp_gated.shape[0]),
        "n_samples": int(hp_gated.shape[1]),
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[interaction] wrote outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include-tnrc6", action="store_true", default=None,
                    help="Include TNRC6A/B/C in the AGO gate (default: follow config.AGO_GATE)")
    ap.add_argument("--out-dir", type=Path, default=C.INTERACTION_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, include_tnrc6=args.include_tnrc6)


if __name__ == "__main__":
    main()
