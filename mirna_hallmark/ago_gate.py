"""AGO / RISC availability as a rate-limiting gate on miRNA repression.

miRNA-mediated repression requires loading the guide into an Argonaute (AGO1-4)
RISC complex (with TNRC6/GW182 effectors). Even abundant, high-evidence miRNAs
cannot repress beyond the available RISC capacity. We summarize per-sample RISC
capacity from AGO (+ optional TNRC6) RNA and map it through a bounded saturating
gate that scales miRNA pressure in the interaction module.

This is a documented sensitivity layer, NOT a causal model: the interaction
module always reports gated AND ungated results side by side.

Outputs (``output/ago_gate/``):
- ``per_sample_ago_capacity.tsv.gz`` -- participant, per-gene z, ago_capacity_z, ago_gate
- ``ago_expression_by_stratum.tsv``  -- mean AGO/RISC log2(TPM+1) per stratum + KW p/q
- ``ago_cnv_by_stratum.tsv``         -- AGO/RISC CNV prevalence per stratum
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.config import AgoGateParams


def gate_from_capacity(capacity_z: pd.Series, params: AgoGateParams = C.AGO_GATE) -> pd.Series:
    """Map capacity z-score -> multiplicative gate in [gate_min, 1] (logistic)."""
    sig = 1.0 / (1.0 + np.exp(-params.gate_k * (capacity_z - params.gate_midpoint)))
    return params.gate_min + (1.0 - params.gate_min) * sig


def compute_capacity(z: pd.DataFrame, params: AgoGateParams = C.AGO_GATE) -> pd.DataFrame:
    """Decompose per-gene z into the two both-required RISC components + capacity.

    ``z`` is a (gene x participant) z-score matrix containing AGO1-4 and, when
    ``include_tnrc6``, TNRC6A/B/C. Returns a (participant x {ago_load_z,
    effector_z, ago_capacity_z}) frame. Capacity is the AGO-load / effector
    *minimum* (co-limitation) by default, or their mean (compensatory blend) when
    ``effector_colimit`` is off.
    """
    weights = dict(params.ago_load_weights)
    ago_genes = [g for g in C.AGO_GENES if g in z.index]
    if not ago_genes:
        raise ValueError("No AGO genes found in the RNA matrix")
    w = pd.Series({g: float(weights.get(g, 1.0)) for g in ago_genes})
    ago_load = z.loc[ago_genes].mul(w, axis=0).sum(axis=0) / w.sum()  # weighted mean z

    out = pd.DataFrame({"ago_load_z": ago_load})
    eff_genes = [g for g in C.RISC_EXTRA_GENES if g in z.index]
    if params.include_tnrc6 and eff_genes:
        effector = z.loc[eff_genes].mean(axis=0)
        out["effector_z"] = effector
        if params.effector_colimit:
            # whichever resource (loaded AGO or TNRC6 effector) is scarcer rate-limits
            out["ago_capacity_z"] = pd.concat([ago_load, effector], axis=1).min(axis=1)
        else:
            out["ago_capacity_z"] = 0.5 * (ago_load + effector)
    else:
        out["ago_capacity_z"] = ago_load
    return out


def compute_ago_gate(
    rna: Optional[pd.DataFrame] = None,
    *,
    params: AgoGateParams = C.AGO_GATE,
) -> pd.DataFrame:
    """Per-participant RISC capacity + gate.

    Returns a DataFrame indexed by participant with per-gene z columns
    (``z_AGO2`` ...), the two component scores (``ago_load_z``, ``effector_z``),
    the co-limited ``ago_capacity_z``, and the multiplicative ``ago_gate``.
    """
    ago = D.load_ago_expression(rna, include_tnrc6=params.include_tnrc6)
    if ago.empty:
        raise ValueError("No AGO/RISC genes found in the RNA matrix")
    z = S.zscore_rows(ago)  # gene x participant
    components = compute_capacity(z, params)  # participant x component scores
    gate = gate_from_capacity(components["ago_capacity_z"], params)

    out = z.T.copy()
    out.columns = [f"z_{g}" for g in out.columns]
    out = out.join(components)
    out["ago_gate"] = gate
    out.index.name = "participant"
    return out


def _summarize_expression_by_stratum(
    ago: pd.DataFrame,
    clinical: pd.DataFrame,
) -> pd.DataFrame:
    """Mean AGO/RISC log2(TPM+1) per stratum-level + Kruskal-Wallis across levels."""
    long = (
        ago.T.rename_axis("participant").reset_index().melt(
            id_vars="participant", var_name="gene", value_name="log2tpm"
        )
    )
    long = long.merge(clinical, on="participant", how="left")

    rows = []
    pvals = []
    keys = []
    for stratum_col, layer in C.STRATUM_SPECS:
        if stratum_col not in clinical.columns:
            continue
        for gene, gsub in long.groupby("gene"):
            h, p, ngroups = S.kruskal_across_strata(gsub["log2tpm"], gsub[stratum_col])
            for level, lsub in gsub.dropna(subset=[stratum_col]).groupby(stratum_col):
                rows.append(
                    {
                        "stratification_layer": layer,
                        "stratum": str(level),
                        "gene": gene,
                        "n": int(lsub["log2tpm"].notna().sum()),
                        "mean_log2tpm": round(float(lsub["log2tpm"].mean()), 4),
                        "median_log2tpm": round(float(lsub["log2tpm"].median()), 4),
                        "kw_H": h,
                        "kw_p": p,
                    }
                )
            keys.append((layer, gene))
            pvals.append(p)
    df = pd.DataFrame(rows)
    # FDR across (layer, gene) KW tests, broadcast back to rows
    if keys:
        q = S.bh_fdr(pvals)
        qmap = {k: qi for k, qi in zip(keys, q)}
        df["kw_q"] = [qmap.get((r.stratification_layer, r.gene), np.nan) for r in df.itertuples()]
    return df


def _summarize_cnv_by_stratum(
    clinical: pd.DataFrame,
    *,
    include_tnrc6: bool,
) -> pd.DataFrame:
    """AGO/RISC copy-number prevalence per stratum level."""
    genes = list(C.AGO_GENES) + (list(C.RISC_EXTRA_GENES) if include_tnrc6 else [])
    # Dedicated cache so we never clobber the Hallmark-universe CNV cache
    cnv = D.load_cnv_target_genes(
        genes, cache_path=C.OUTPUT_ROOT / "matrices" / "cnv_ago_risc.tsv.gz"
    )
    long = (
        cnv.T.rename_axis("participant").reset_index().melt(
            id_vars="participant", var_name="gene", value_name="copy_number"
        )
    )
    long["cn_state"] = D.classify_cn(long["copy_number"])
    long = long.merge(clinical, on="participant", how="left")

    rows = []
    for stratum_col, layer in C.STRATUM_SPECS:
        if stratum_col not in clinical.columns:
            continue
        for (level, gene), sub in long.dropna(subset=[stratum_col, "copy_number"]).groupby(
            [stratum_col, "gene"]
        ):
            n = len(sub)
            st = sub["cn_state"]
            rows.append(
                {
                    "stratification_layer": layer,
                    "stratum": str(level),
                    "gene": gene,
                    "n": n,
                    "mean_copy_number": round(float(sub["copy_number"].mean()), 4),
                    "pct_gain": round(100.0 * st.isin(["gain", "amp"]).mean(), 2),
                    "pct_loss": round(100.0 * (st == "loss").mean(), 2),
                    "pct_deep_del": round(100.0 * (st == "deep_del").mean(), 2),
                }
            )
    return pd.DataFrame(rows)


def run(
    *,
    out_dir: Path = C.AGO_GATE_DIR,
    include_tnrc6: Optional[bool] = None,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    params = C.AGO_GATE if include_tnrc6 is None else replace(C.AGO_GATE, include_tnrc6=include_tnrc6)

    print("[ago_gate] loading RNA + computing RISC capacity ...")
    # ⚠⚠ STATE-BLIND BY CONSTRUCTION (flagged MH-223). `load_rna()` collapses barcodes to 12-char
    # participants by MEAN over ALL sample types, so for the 103 participants with both a tumour and a NAT
    # sample the capacity below is a **tumour/NAT average**, not a tumour measurement. The emitted
    # `per_sample_ago_capacity.tsv.gz` has 1,095 participant keys and **zero NAT keys** — so it LOOKS
    # available for all 103 paired patients and is not state-resolved for any of them.
    # ⛔ DO NOT condition a NAT or a paired-Δ model on this file. A state-resolved RISC capacity must be
    # rebuilt from `learned.states.state_matrices(sample_type)` — AGO1-4 / TNRC6A-C / DICER1 / DROSHA /
    # DGCR8 / XPO5 are all mRNA rows present in BOTH states, so it costs no new data.
    rna = D.load_rna()
    gate_df = compute_ago_gate(rna, params=params)
    gate_df.to_csv(out_dir / "per_sample_ago_capacity.tsv.gz", sep="\t", compression="gzip")
    print(f"[ago_gate] capacity for {len(gate_df):,} participants; "
          f"gate range [{gate_df['ago_gate'].min():.3f}, {gate_df['ago_gate'].max():.3f}]")

    clinical = D.load_clinical_strata()
    ago = D.load_ago_expression(rna, include_tnrc6=params.include_tnrc6)

    expr_strata = _summarize_expression_by_stratum(ago, clinical)
    expr_strata.to_csv(out_dir / "ago_expression_by_stratum.tsv", sep="\t", index=False)

    cnv_strata = _summarize_cnv_by_stratum(clinical, include_tnrc6=params.include_tnrc6)
    cnv_strata.to_csv(out_dir / "ago_cnv_by_stratum.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.ago_gate",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "ago_genes": list(C.AGO_GENES),
        "risc_extra_genes": list(C.RISC_EXTRA_GENES),
        "include_tnrc6": params.include_tnrc6,
        "gate_model": {
            "ago_load": "AGO-weighted mean per-gene z of AGO1-4 (loadable miRISC core)",
            "ago_load_weights": dict(params.ago_load_weights),
            "effector": "mean per-gene z of TNRC6A/B/C (GW182 effector)" if params.include_tnrc6 else None,
            "capacity": (
                "min(ago_load, effector) -- co-limitation (Liebig minimum)"
                if (params.include_tnrc6 and params.effector_colimit)
                else ("0.5*(ago_load + effector) -- compensatory blend" if params.include_tnrc6 else "ago_load only")
            ),
            "effector_colimit": params.effector_colimit,
            "gate": "gate_min + (1-gate_min) * sigmoid(gate_k * (capacity_z - gate_midpoint))",
            "gate_min": params.gate_min,
            "gate_k": params.gate_k,
            "gate_midpoint": params.gate_midpoint,
            "note": "documented sensitivity layer; interaction module reports gated AND ungated",
        },
        "n_participants": int(len(gate_df)),
        "strata": [layer for _, layer in C.STRATUM_SPECS],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[ago_gate] wrote outputs under {out_dir}")
    return gate_df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include-tnrc6", action="store_true",
                    help="(default on) include TNRC6A/B/C effector co-limitation")
    ap.add_argument("--no-tnrc6", action="store_true",
                    help="AGO-load only, no TNRC6 effector term (sensitivity variant)")
    ap.add_argument("--out-dir", type=Path, default=C.AGO_GATE_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    include_tnrc6: Optional[bool] = False if args.no_tnrc6 else None
    run(out_dir=args.out_dir, include_tnrc6=include_tnrc6)


if __name__ == "__main__":
    main()
