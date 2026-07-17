"""Canonical per-edge coupling significance for the learned model — reuses the precursor's family-aware
Freedman-Lane permutation partial-Spearman (`coupling_inference.permutation_pvalues`), which is generic
(any predictor/target/covariate frames; no dependence on the precursor's pressure construction).

The ONE adaptation: `permutation_pvalues` takes a single covariate frame for all rows, but our confounder
block C carries a PER-GENE term (`target_cn`). So we call it **per gene** — the gene's regulator arms as
rows, the gene's mRNA as the (repeated) target, the gene's full C as covariates, seed family per arm — which
is exactly correct (C is shared *within* a gene) and preserves the seed-family permutation dependence. This
gives, per edge: permutation p, smooth p_z, BH/BY q, and **q_family** (Simes-per-family) — the seed family is
the hypothesis unit (Design §F). Genome-wide, we combine across genes with a second-level BY.

CLI: `python -m mirna_hallmark.learned.coupling PTEN` · `--all` (genome-wide, background-friendly).
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.coupling_inference import permutation_pvalues, benjamini_yekutieli


def edge_coupling(gene: str, *, n_perm: int = 2000, deconv: bool = False, orphans: bool = False,
                  seed: int = 0) -> pd.DataFrame:
    """Per-edge coupling significance for one gene: family-aware permutation partial-Spearman of each
    regulator arm vs the gene's mRNA | C. Returns rho, p_perm, p_z, q_bh, q_by, q_family (seed-family unit)."""
    Y, X, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=deconv, orphans=orphans)
    arms = list(X.columns)
    samples = X.index
    pred = X.T                                                       # arm × sample
    tgt = pd.DataFrame(np.repeat(Y.to_numpy(float)[None, :], len(arms), axis=0),
                       index=arms, columns=samples)                  # same gene mRNA for every arm row
    fam = FAM.family_of(pd.Index(arms)).reindex(arms).astype(str).tolist()
    out = permutation_pvalues(pred, tgt, covariates=C, families=fam, n_perm=n_perm, tail="neg", seed=seed)
    out.index.name = "arm"
    out = out.reset_index()
    out.insert(0, "gene", gene)
    out["family"] = fam
    return out


def run(genes=None, *, n_perm: int = 2000, deconv: bool = False,
        out: str = "mirna_hallmark/output/learned/edge_coupling.tsv", progress: int = 50) -> pd.DataFrame:
    """Genome-wide per-edge coupling significance. Per gene: family-aware permutation (within-gene C + seed
    families). Then a **second-level Benjamini-Yekutieli** across all edges (q_by_global) — the edge unit;
    q_family is the within-gene seed-family combination. Background-friendly."""
    from pathlib import Path
    from mirna_hallmark.learned.evidence import ledger as LG
    genes = genes or sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))  # POOLED-HE gene universe (migration)
    parts = []
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[edge_coupling] {i}/{len(genes)}", flush=True)
        try:
            parts.append(edge_coupling(g, n_perm=n_perm, deconv=deconv))
        except Exception:
            continue
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if len(df):
        df["q_by_global"] = benjamini_yekutieli(df["p_z"].values)     # edge-level, genome-wide
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False)
        sig = df[df["q_by_global"] < 0.05]
        famsig = df[df["q_family"] < 0.05]
        print(f"\n=== GENOME-WIDE EDGE COUPLING ({len(df)} edges / {df['gene'].nunique()} genes) ===")
        print(f"edge-significant (q_BY_global<0.05 & rho<0): {int((sig['rho']<0).sum())} "
              f"| family-significant (q_family<0.05 & rho<0): {int((famsig['rho']<0).sum())} "
              f"| distinct sig seed-family→gene: {famsig[famsig['rho']<0].assign(h=famsig['family']+'|'+famsig['gene'])['h'].nunique()}")
        print(f"wrote {out}")
    return df


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if "--all" in sys.argv:
        run(n_perm=int(next((a.split("=")[1] for a in sys.argv if a.startswith("--nperm=")), 2000)))
    else:
        for g in (args or ["PTEN"]):
            r = edge_coupling(g)
            with pd.option_context("display.width", 180):
                print(f"\n=== {g} — per-edge coupling (family-aware permutation) ===")
                print(r.sort_values("rho")[["arm", "family", "rho", "p_perm", "p_z", "q_by", "q_family"]]
                      .head(15).round(4).to_string(index=False))
