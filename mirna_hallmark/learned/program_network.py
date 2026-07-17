"""Program-level view (a) — the miRNA(family) × gene NETWORK of canonical repression weights, rolled up from
the per-gene canonical fits (states.canonical_M, bagged z-NNLS on fixed family support). This is the honest
program-level object: it AGGREGATES per-gene fits (no cross-gene pooling bias — the hierarchical Gibbs pooling
does NOT improve coupling here, gate 2026-07-05: Δρ −0.002 full-n / −0.045 at n=150, because the program's
genes barely share regulators and effects are gene-specific). See LEARNED_MODEL_ESTIMATOR_MAP §8.

Three levels of summary:
  1. per-FAMILY (miRNA) engagement — how many program genes it represses, total/mean weight, FOCUS
     (max/sum = single-target vs broad regulator), top target.
  2. per-GENE incoming — how many families regulate it, total incoming weight, top regulator.
  3. NETWORK — shape, density, hub families (co-target ≥2 genes = where a program view adds over per-gene),
     concentration (share of total program weight held by the top hubs), most-regulated genes.

CLI: `python -m mirna_hallmark.learned.program_network EMT` (or a Hallmark name / gene list).
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from mirna_hallmark.learned import states as ST


def program_matrix(genes, *, sample_type: str = "01") -> pd.DataFrame:
    """family × gene canonical-M matrix (union of the genes' HE seed families; 0 where a family doesn't
    regulate a gene). Family = the identified estimand (arm_level=False), so nodes aren't within-family dupes."""
    cols = {}
    for g in genes:
        try:
            m = ST.canonical_M(g, sample_type, arm_level=False)
        except Exception:
            continue                                          # gene with no HE regulators / not in matrix
        m = m[m > 0]
        if len(m):
            cols[g] = m
    if not cols:
        return pd.DataFrame()
    return pd.DataFrame(cols).fillna(0.0)                      # family × gene


def family_engagement(W: pd.DataFrame) -> pd.DataFrame:
    """Per-family (miRNA) program engagement — the 'program miRNA status'."""
    pos = W > 0
    sumw = W.sum(axis=1)
    return pd.DataFrame({
        "n_targets": pos.sum(axis=1).astype(int),
        "sum_weight": sumw.round(3),
        "mean_weight": (sumw / pos.sum(axis=1).clip(lower=1)).round(3),
        "focus": (W.max(axis=1) / sumw.replace(0, np.nan)).round(2),   # 1.0 = single-target; low = broad regulator
        "top_target": W.idxmax(axis=1),
    }).sort_values(["n_targets", "sum_weight"], ascending=False)


def gene_incoming(W: pd.DataFrame) -> pd.DataFrame:
    """Per-gene incoming repression — how many families, total weight, dominant regulator."""
    pos = W > 0
    return pd.DataFrame({
        "n_regulators": pos.sum(axis=0).astype(int),
        "sum_incoming": W.sum(axis=0).round(3),
        "top_regulator": W.idxmax(axis=0),
        "top_share": (W.max(axis=0) / W.sum(axis=0).replace(0, np.nan)).round(2),
    }).sort_values("sum_incoming", ascending=False)


def program_network(genes, *, sample_type: str = "01", label: str = "", top: int = 12):
    """Roll up the per-gene canonical M into the program network + the three summaries. Returns
    (W, family_engagement, gene_incoming)."""
    W = program_matrix(genes, sample_type=sample_type)
    if W.empty:
        print(f"[{label}] no canonical weights"); return W, pd.DataFrame(), pd.DataFrame()
    fam = family_engagement(W)
    gen = gene_incoming(W)
    n_fam, n_gene = W.shape
    n_edges = int((W > 0).sum().sum())
    density = n_edges / (n_fam * n_gene)
    hubs = fam[fam["n_targets"] >= 2]                          # co-targeting families = the shared regulators
    tot = float(W.to_numpy().sum())
    hub_share = float(W.loc[hubs.index].to_numpy().sum() / tot) if tot > 0 and len(hubs) else 0.0
    top_hub_share = float(W.loc[fam.head(top).index].to_numpy().sum() / tot) if tot > 0 else np.nan

    print(f"\n===== PROGRAM NETWORK{' — ' + label if label else ''}  ({n_gene} genes × {n_fam} families) =====")
    print(f"edges (family→gene, M>0): {n_edges}  |  density {density:.2f}  |  "
          f"co-targeting hub families (≥2 genes): {len(hubs)}  holding {hub_share:.0%} of program weight  |  "
          f"top-{top} families hold {top_hub_share:.0%}")
    print(f"\n-- per-family engagement (program miRNA status; focus 1.0=single-target, low=broad) — top {top} --")
    print(fam.head(top).to_string())
    print(f"\n-- HUB families (co-target ≥2 program genes — where the program view adds over per-gene) --")
    print((hubs.head(top).to_string() if len(hubs) else "  (none — genes don't share regulators; per-gene suffices)"))
    print(f"\n-- per-gene incoming repression — top {top} --")
    print(gen.head(top).to_string())
    return W, fam, gen


_ALIAS = {"EMT": "EPITHELIAL_MESENCHYMAL_TRANSITION", "P53": "P53_PATHWAY", "G2M": "G2M_CHECKPOINT",
          "IFN": "INTERFERON_GAMMA_RESPONSE", "IFNG": "INTERFERON_GAMMA_RESPONSE", "HYPOXIA": "HYPOXIA"}


def module_pooling_test(anchor: str, genes, *, sample_type: str = "01"):
    """(b) Does cross-gene pooling help on a CO-REGULATED module (genes sharing the `anchor` family)? Extracts
    the anchor's targets within `genes` and runs the hierarchical pool-vs-solo gate. Empirically (2026-07-05):
    pooling is harmful on the disjoint p53/cell-cycle set (Δρ −0.045 at n=150) but neutral/marginally helpful
    on the miR-29 ECM module (+0.001 at n=150) — sign tracks how many regulators the genes SHARE, not a
    hyperparameter. The HE network is sparse (density ~0.02) so pooling rarely wins; per-gene stays primary."""
    from mirna_hallmark.learned import hierarchical as H
    W = program_matrix(genes, sample_type=sample_type)
    if anchor not in W.index:
        print(f"{anchor} not in program network"); return None
    mod = [g for g in W.columns if W.loc[anchor, g] > 0]
    shared = (W[mod] > 0).sum(axis=1)
    print(f"module = {anchor}'s targets: {len(mod)} genes | families co-regulating ≥3 of them: "
          f"{int((shared >= 3).sum())} (shared-regulator density → whether pooling can help)")
    return H.oof(mod, n_sub=150)


def _resolve(arg: str):
    """A Hallmark short-name (EMT/P53/G2M/IFN/HYPOXIA) → its HE gene set, else treat args as gene symbols."""
    from mirna_hallmark.learned.programs import _program_genes
    progs = _program_genes()
    want = _ALIAS.get(arg.upper(), arg.upper())
    for k in progs:
        if want == k.upper() or want in k.upper():
            return progs[k], k
    return None, arg


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if not args:
        args = ["EMT"]
    genes, label = _resolve(args[0])
    if genes is None:                                         # explicit gene list
        genes, label = args, "custom"
    program_network(genes, label=label)
