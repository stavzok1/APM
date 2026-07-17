"""Within-subtype gene-pressure NMF — per-sample signature activation analysis + heatmaps.

The `mirna_comovement._nmf_gene_pressure` within-subtype block fits an NMF *inside* each PAM50 group
but persists only the top-gene loadings, not the per-sample factor activation (the H = factor x sample
matrix). This module refits that exact decomposition and answers a different question: **how do the
within-subtype signatures distribute across the samples of that subtype?** i.e. do tumours of a subtype
split into clean sub-groups each dominated by one signature, or is every tumour a graded mix of several?

For each PAM50 subtype it:
  1. refits the gene x sample total-pressure NMF on that subtype's columns (k = NMF_K),
  2. column-normalises H so each sample's factor activations sum to 1 (a "signature composition"),
  3. per sample: dominant factor, dominant share, Shannon entropy, effective #factors (exp(entropy)),
     and a dominant/mixed call (dominant share >= DOMINANT_SHARE),
  4. writes a per-sample table + a per-subtype summary, and renders a signature x sample heatmap
     (rows = signatures, columns = samples ordered by dominant factor then descending share).

Default = abundance pressure (matches the shipped within-subtype tables); `--acquired` uses the
dev/healthy-anchored gain pressure.
"""
from __future__ import annotations
import argparse
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from mirna_hallmark.analyses.misc.mirna_comovement import NMF_K, OUT_DIR

DOMINANT_SHARE = 0.50          # a sample is "dominant" if its top signature holds >= this share
FIG_DIR = OUT_DIR / "figures"


def _fit_subtype(P_sub: pd.DataFrame, prefix: str, k: int):
    """NMF on a (gene x sample) non-negative matrix -> W (gene x factor), H (factor x sample)."""
    from sklearn.decomposition import NMF
    m = P_sub.loc[P_sub.sum(axis=1) > 0]
    if min(m.shape) < k:
        return None, None
    model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=1000)
    W = pd.DataFrame(model.fit_transform(m.values), index=m.index,
                     columns=[f"{prefix}{j + 1}" for j in range(k)])
    H = pd.DataFrame(model.components_, index=W.columns, columns=m.columns)
    return W, H


def _sample_composition(H: pd.DataFrame) -> pd.DataFrame:
    """Column-normalise H so each sample's signature activations sum to 1; derive mixing metrics."""
    comp = H.div(H.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)   # factor x sample, columns sum 1
    dom_factor = comp.idxmax(axis=0)
    dom_share = comp.max(axis=0)
    # Shannon entropy over the composition; effective number of active signatures = exp(H)
    p = comp.where(comp > 0)
    ent = -(p * np.log(p)).sum(axis=0)
    eff = np.exp(ent)
    out = pd.DataFrame({
        "sample": comp.columns,
        "dominant_factor": dom_factor.values,
        "dominant_share": dom_share.values,
        "entropy": ent.values,
        "eff_n_signatures": eff.values,
        "call": np.where(dom_share.values >= DOMINANT_SHARE, "dominant", "mixed"),
    })
    return out, comp


def _heatmap(comp: pd.DataFrame, samp: pd.DataFrame, subtype: str, labels: Dict[str, str],
             out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = samp.sort_values(["dominant_factor", "dominant_share"],
                             ascending=[True, False])["sample"].tolist()
    M = comp[order]
    # order rows (factors) by how many samples they dominate (descending) for a block-diagonal look
    dom_counts = samp["dominant_factor"].value_counts()
    row_order = sorted(comp.index, key=lambda f: -dom_counts.get(f, 0))
    M = M.loc[row_order]
    ylabels = [f"{f}  {labels.get(f, '')}" for f in row_order]

    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.045 + 3), 0.42 * len(row_order) + 1.6))
    im = ax.imshow(M.values, aspect="auto", cmap="magma", vmin=0, vmax=min(1.0, float(M.values.max())))
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xticks([])
    frac_dom = (samp["call"] == "dominant").mean()
    ax.set_xlabel(f"{len(order)} {subtype} tumours (ordered by dominant signature)  "
                  f"|  {frac_dom:.0%} dominant, median top-share {samp['dominant_share'].median():.2f}, "
                  f"median eff #sigs {samp['eff_n_signatures'].median():.1f}")
    ax.set_title(f"{subtype}: within-subtype gene-pressure signatures x samples (column-normalised)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.01)
    cbar.set_label("signature share of sample", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def run(*, acquired: bool = False, out_dir: Path = OUT_DIR, k: int = NMF_K) -> Dict[str, pd.DataFrame]:
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices, _PAM50, ATTR_MODE
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    from mirna_hallmark.analyses.misc.mirna_comovement import (HallmarkSets as HS, _gene_factor_label,
                                                  _load_state_class, STATE_DIR)
    from mirna_hallmark.config import PRESSURE_HEALTHY_ANCHOR_MODE
    from mirna_hallmark.healthy_anchor import load_baseline

    tag = "_acquired" if acquired else ""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    hs = HS.load()
    _, edge = _load_state_class(STATE_DIR)
    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None
    tumor = _state_matrices()["tumor"]

    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    mode = PRESSURE_HEALTHY_ANCHOR_MODE if acquired else ATTR_MODE
    kw = {"healthy_baseline": load_baseline()} if acquired else {}
    P = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode=mode,  # type: ignore[arg-type]
                              resolve_arms=False, **kw)
    P = P.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    P = P.loc[P.sum(axis=1) > 0]

    clin = D.load_clinical_strata()
    part_pam = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
                if "PAM50_final" in clin.columns else {})
    col_pam = pd.Series({c: part_pam.get(_participant(c)) for c in P.columns})

    per_sample, summary = [], []
    for sub in _PAM50:
        cols = [c for c in col_pam.index[col_pam == sub] if c in P.columns]
        if len(cols) < 30:
            continue
        W, H = _fit_subtype(P[cols], f"{sub}_G", k)
        if W is None:
            continue
        labels = {f: (_gene_factor_label(list(W[f].sort_values(ascending=False).head(30).index), hs, ec)
                      .get("dominant_hallmark") or "").replace("HALLMARK_", "")[:22] for f in W.columns}
        samp, comp = _sample_composition(H)
        samp.insert(0, "subtype", sub)
        per_sample.append(samp)
        out_png = FIG_DIR / f"within_subtype_signatures{tag}_{sub}.png"
        _heatmap(comp, samp, sub, labels, out_png)
        n = len(samp)
        summary.append({
            "subtype": sub, "n_samples": n, "k_signatures": k,
            "frac_dominant": float((samp["call"] == "dominant").mean()),
            "median_dominant_share": float(samp["dominant_share"].median()),
            "median_eff_n_signatures": float(samp["eff_n_signatures"].median()),
            "n_signatures_that_dominate_any": int(samp["dominant_factor"].nunique()),
            "biggest_signature": samp["dominant_factor"].value_counts().idxmax(),
            "biggest_signature_hallmark": labels.get(samp["dominant_factor"].value_counts().idxmax(), ""),
            "biggest_signature_n": int(samp["dominant_factor"].value_counts().max()),
            "figure": str(out_png.relative_to(OUT_DIR.parents[2])),
        })
        print(f"[within_subtype] {sub}: n={n}  dominant={summary[-1]['frac_dominant']:.0%}  "
              f"median top-share={summary[-1]['median_dominant_share']:.2f}  "
              f"eff#sigs={summary[-1]['median_eff_n_signatures']:.1f}  -> {out_png.name}")

    ps = pd.concat(per_sample, ignore_index=True) if per_sample else pd.DataFrame()
    sm = pd.DataFrame(summary)
    if not ps.empty:
        ps.to_csv(out_dir / f"within_subtype_sample_signatures{tag}.tsv", sep="\t", index=False)
        sm.to_csv(out_dir / f"within_subtype_signature_summary{tag}.tsv", sep="\t", index=False)
    return {"per_sample": ps, "summary": sm}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--acquired", action="store_true",
                    help="use dev/healthy-anchored gain pressure instead of abundance")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    run(acquired=args.acquired, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
