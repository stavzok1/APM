"""Per-sample NMF signature composition analysis + heatmaps for every NMF with a sample axis.

The comovement module fits several NMFs but persists only top loadings / PAM50 means — not the
factor×sample activation matrix (H). This refits each decomposition on its **appropriate cohort**
and asks: do tumours split into discrete signature-dominated sub-groups, or carry graded mixtures?

Analyses (each → per-sample TSV + heatmap PNG under `mirna_comovement/nmf_sample_signatures/`):

| analysis_id | matrix | cohort | factors |
|---|---|---|---|
| `gene_pressure_cohort` | gene total-pressure × sample | all tumors | P1–P12 |
| `gene_pressure_acquired_cohort` | dev-anchored gain pressure | all tumors | Pa1–Pa12 |
| `gene_pressure_signed_cohort` | ± gain/loss channels × sample | all tumors | S1–S12 |
| `gene_pressure_within_{PAM50}` | gene total-pressure | one subtype | {sub}_G1–G12 |
| `gene_pressure_acquired_within_{PAM50}` | dev gain pressure | one subtype | {sub}_Ga1–Ga12 |
| `mirna_within_tumor_cohort` | arm abundance × sample | all tumors | T1–T12 |
| `mirna_within_{PAM50}` | arm abundance × sample | one subtype | {sub}_S1–S8 |

**Not included** (no per-sample axis): aggregate `c_tumor` NMF (`nmf_programs`), share NMF
(`nmf_share`), per-subtype *static* arm NMF (`nmf_subtype_static` — cells are miRNA×gene, not samples).

Cohort-wide heatmaps order samples by PAM50 block then dominant signature; within-subtype heatmaps
order by dominant signature only.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark.analyses.misc.mirna_comovement import (
    NMF_K, OUT_DIR, STATE_DIR, TOP_LOAD, HallmarkSets, MIN_EXPR,
    _gene_factor_label, _load_state_class,
)

DOMINANT_SHARE = 0.50
K_SUBTYPE_ARM = 8
OUT_ROOT = OUT_DIR / "nmf_sample_signatures"
FIG_DIR = OUT_ROOT / "figures"
SAMPLE_DIR = OUT_ROOT / "per_sample"


@dataclass
class NmfFit:
    W: pd.DataFrame          # features × factors
    H: pd.DataFrame          # factors × samples
    labels: Dict[str, str]   # factor → short annotation
    cohort_label: str        # e.g. "all tumors" or "Basal only"


def _fit_nmf(mat: pd.DataFrame, prefix: str, k: int) -> Optional[NmfFit]:
    from sklearn.decomposition import NMF
    m = mat.loc[mat.sum(axis=1) > 0]
    if min(m.shape) < k:
        return None
    model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=1000)
    W = pd.DataFrame(model.fit_transform(m.values), index=m.index,
                     columns=[f"{prefix}{j + 1}" for j in range(k)])
    H = pd.DataFrame(model.components_, index=W.columns, columns=m.columns)
    return NmfFit(W=W, H=H, labels={}, cohort_label="")


def _composition(H: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    comp = H.div(H.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    p = comp.where(comp > 0)
    ent = -(p * np.log(p)).sum(axis=0)
    samp = pd.DataFrame({
        "sample": comp.columns,
        "dominant_factor": comp.idxmax(axis=0).values,
        "dominant_share": comp.max(axis=0).values,
        "entropy": ent.values,
        "eff_n_signatures": np.exp(ent.values),
        "call": np.where(comp.max(axis=0).values >= DOMINANT_SHARE, "dominant", "mixed"),
    })
    return samp, comp


def _heatmap(comp: pd.DataFrame, samp: pd.DataFrame, fit: NmfFit, analysis_id: str,
             out_png: Path, col_pam: Optional[pd.Series] = None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order: List[str] = []
    if col_pam is not None and samp["sample"].isin(col_pam.index).any():
        subs = ["LumA", "LumB", "Her2", "Basal"]
        for sub in subs:
            sc = samp[samp["sample"].map(col_pam).eq(sub)]
            order += sc.sort_values(["dominant_factor", "dominant_share"],
                                    ascending=[True, False])["sample"].tolist()
        rest = [s for s in samp["sample"] if s not in order]
        order += rest
    else:
        order = samp.sort_values(["dominant_factor", "dominant_share"],
                                 ascending=[True, False])["sample"].tolist()

    M = comp[order]
    dom_counts = samp["dominant_factor"].value_counts()
    row_order = sorted(comp.index, key=lambda f: -dom_counts.get(f, 0))
    M = M.loc[row_order]
    ylabels = [f"{f}  {fit.labels.get(f, '')[:24]}" for f in row_order]

    fig_w = max(10, len(order) * 0.04 + 4)
    fig_h = 0.42 * len(row_order) + 2.4
    has_pam = col_pam is not None and any(col_pam.get(s) for s in order)
    if has_pam:
        import matplotlib.gridspec as gridspec
        fig = plt.figure(figsize=(fig_w, fig_h))
        gs = gridspec.GridSpec(2, 1, height_ratios=[len(row_order), 0.35], hspace=0.02)
        ax = fig.add_subplot(gs[0])
        ax_bar = fig.add_subplot(gs[1], sharex=ax)
    else:
        fig, ax = plt.subplots(figsize=(fig_w, fig_h - 0.4))
        ax_bar = None

    im = ax.imshow(M.values, aspect="auto", cmap="magma", vmin=0, vmax=min(1.0, float(M.values.max())))

    if has_pam and ax_bar is not None:
        pam_colors = {"LumA": "#4C72B0", "LumB": "#DD8452", "Her2": "#55A868", "Basal": "#C44E52"}
        bar = np.array([[pam_colors.get(col_pam.get(s), 0.5) for s in order]], dtype=object)
        for j, s in enumerate(order):
            ax_bar.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, 1,
                                           color=pam_colors.get(col_pam.get(s), "#888888")))
        ax_bar.set_xlim(-0.5, len(order) - 0.5)
        ax_bar.set_ylim(-0.5, 0.5)
        ax_bar.set_yticks([])
        ax_bar.set_xticks([])
        pos = 0
        for sub in ["LumA", "LumB", "Her2", "Basal"]:
            n = sum(1 for s in order if col_pam.get(s) == sub)
            if n and pos:
                ax.axvline(pos - 0.5, color="white", lw=0.8, alpha=0.7)
                ax_bar.axvline(pos - 0.5, color="white", lw=0.8, alpha=0.7)
            pos += n

    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xticks([])
    frac_dom = (samp["call"] == "dominant").mean()
    ax.set_xlabel(
        f"{len(order)} samples ({fit.cohort_label})  |  {frac_dom:.0%} dominant, "
        f"median top-share {samp['dominant_share'].median():.2f}, "
        f"median eff #sigs {samp['eff_n_signatures'].median():.1f}",
        fontsize=9,
    )
    ax.set_title(f"{analysis_id}: signatures × samples (column-normalised)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01).set_label("share", fontsize=8)
    if col_pam is not None:
        fig.legend(handles=[Patch(facecolor=c, label=s) for s, c in
                          {"LumA": "#4C72B0", "LumB": "#DD8452", "Her2": "#55A868", "Basal": "#C44E52"}.items()],
                   loc="upper right", fontsize=7, frameon=False)
    if has_pam and ax_bar is not None:
        fig.subplots_adjust(hspace=0.02)
    else:
        fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def _pam50_maps(cols: Sequence[str]) -> pd.Series:
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.family_normal_reference import _participant
    clin = D.load_clinical_strata()
    part_pam = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
                if "PAM50_final" in clin.columns else {})
    return pd.Series({c: part_pam.get(_participant(c)) for c in cols})


def _gene_labels(W: pd.DataFrame, hs: HallmarkSets, ec) -> Dict[str, str]:
    return {f: (_gene_factor_label(list(W[f].sort_values(ascending=False).head(30).index), hs, ec)
                  .get("dominant_hallmark") or "").replace("HALLMARK_", "")[:24]
            for f in W.columns}


def _signed_labels(W: pd.DataFrame, hs: HallmarkSets) -> Dict[str, str]:
    labels = {}
    for f in W.columns:
        gmass = W[f][W.index.str.endswith("|gain")].sum()
        lmass = W[f][W.index.str.endswith("|loss")].sum()
        d = "gain" if gmass >= lmass else "loss"
        top = W[f].sort_values(ascending=False).head(5).index
        genes = [g.split("|")[0] for g in top if "|" in g][:3]
        labels[f] = f"{d} {','.join(genes)}"
    return labels


def _load_signed_fit(comov_dir: Path) -> Optional[NmfFit]:
    """Reuse mirna_comovement's persisted signed-NMF fit (full W/H) instead of re-fitting a
    second, independent copy. Two fits return permutation-arbitrary factor indices, so the
    LIFR/ABCE1/IGFBP5 backbone could land on S9 in the loadings/factor_subtype tables but S5
    in the per-sample/summary tables — a bare-index join across them is then silently wrong.
    Reusing one fit makes the factor indices AND the stable `factor_label` identical everywhere.
    Returns None when comovement has not been run yet (caller falls back to a standalone re-fit)."""
    wf = comov_dir / "nmf_gene_pressure_signed_W.tsv"
    hf = comov_dir / "nmf_gene_pressure_signed_H.tsv"
    if not (wf.exists() and hf.exists()):
        return None
    W = pd.read_csv(wf, sep="\t").set_index("gene_channel")
    H = pd.read_csv(hf, sep="\t").set_index("factor")
    H.columns = [str(c) for c in H.columns]
    labels: Dict[str, str] = {}
    sf = comov_dir / "nmf_gene_pressure_signed_factor_subtype.tsv"
    if sf.exists():
        fsub = pd.read_csv(sf, sep="\t")
        if "factor_label" in fsub.columns:
            labels = dict(zip(fsub["factor"].astype(str), fsub["factor_label"].astype(str)))
    return NmfFit(W=W, H=H, labels=labels, cohort_label="all tumors (signed gain+loss)")


def _arm_labels(W: pd.DataFrame, mir_class: Optional[pd.DataFrame] = None) -> Dict[str, str]:
    labels = {}
    cls = mir_class.set_index("miRNA") if mir_class is not None and not mir_class.empty else None
    for f in W.columns:
        top = W[f].sort_values(ascending=False).head(3).index
        short = ";".join(a.replace("hsa-", "") for a in top)
        sec = ""
        if cls is not None and top[0] in cls.index:
            sec = str(cls.loc[top[0]].get("secondary_class", ""))[:12]
        labels[f] = f"{short} [{sec}]" if sec else short
    return labels


def _run_one(analysis_id: str, fit: NmfFit, col_pam: Optional[pd.Series],
             out_root: Path = OUT_ROOT) -> pd.DataFrame:
    samp, comp = _composition(fit.H)
    samp.insert(0, "analysis_id", analysis_id)
    if col_pam is not None:
        samp["PAM50"] = samp["sample"].map(col_pam).fillna("")
    out_png = FIG_DIR / f"{analysis_id}.png"
    _heatmap(comp, samp, fit, analysis_id, out_png, col_pam=col_pam)
    samp.to_csv(SAMPLE_DIR / f"{analysis_id}.tsv", sep="\t", index=False)
    summary = {
        "analysis_id": analysis_id,
        "cohort": fit.cohort_label,
        "n_samples": len(samp),
        "k_signatures": len(fit.H),
        "frac_dominant": float((samp["call"] == "dominant").mean()),
        "median_dominant_share": float(samp["dominant_share"].median()),
        "median_eff_n_signatures": float(samp["eff_n_signatures"].median()),
        "n_signatures_active": int(samp["dominant_factor"].nunique()),
        "top_signature": samp["dominant_factor"].value_counts().idxmax(),
        "top_signature_label": fit.labels.get(samp["dominant_factor"].value_counts().idxmax(), ""),
        "top_signature_n": int(samp["dominant_factor"].value_counts().max()),
        "figure": str(out_png.relative_to(OUT_DIR.parents[2])),
    }
    print(f"[nmf_sample] {analysis_id}: n={summary['n_samples']} cohort={fit.cohort_label}  "
          f"dominant={summary['frac_dominant']:.0%}  top-share={summary['median_dominant_share']:.2f}  "
          f"eff#={summary['median_eff_n_signatures']:.1f}  -> {out_png.name}")
    return pd.DataFrame([summary])


def run(*, acquired_within: bool = True, out_root: Path = OUT_ROOT) -> pd.DataFrame:
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices, _PAM50, ATTR_MODE
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    from mirna_hallmark.config import PRESSURE_HEALTHY_ANCHOR_MODE
    from mirna_hallmark.healthy_anchor import load_baseline

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SAMPLE_DIR.mkdir(parents=True, exist_ok=True)

    hs = HallmarkSets.load()
    mir_class, edge = _load_state_class(STATE_DIR)
    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None
    tumor = _state_matrices()["tumor"]
    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    col_pam_all = _pam50_maps(tumor.columns)

    summaries: List[pd.DataFrame] = []

    # --- gene-pressure cohort (abundance) ---
    P = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode=ATTR_MODE,  # type: ignore[arg-type]
                              resolve_arms=False)
    P = P.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0).loc[lambda d: d.sum(axis=1) > 0]
    fit = _fit_nmf(P, "P", NMF_K)
    if fit:
        fit.labels = _gene_labels(fit.W, hs, ec)
        fit.cohort_label = "all tumors"
        summaries.append(_run_one("gene_pressure_cohort", fit, col_pam_all, out_root))

    # --- gene-pressure acquired cohort ---
    base = load_baseline()
    Pa = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode=PRESSURE_HEALTHY_ANCHOR_MODE,  # type: ignore[arg-type]
                               resolve_arms=False, healthy_baseline=base)
    Pa = Pa.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0).loc[lambda d: d.sum(axis=1) > 0]
    fit_a = _fit_nmf(Pa, "Pa", NMF_K)
    if fit_a:
        fit_a.labels = _gene_labels(fit_a.W, hs, ec)
        fit_a.cohort_label = "all tumors (acquired gain)"
        summaries.append(_run_one("gene_pressure_acquired_cohort", fit_a, col_pam_all, out_root))

    # --- signed gene-pressure cohort ---
    # Reuse the comovement signed-NMF fit so the per-sample dominant factor indices + labels
    # match nmf_gene_pressure_signed_{loadings,factor_subtype}.tsv exactly (single source of
    # truth). Only fall back to a standalone re-fit if comovement has not been run.
    fit_s = _load_signed_fit(OUT_DIR)
    if fit_s is None:
        Pdev = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode=PRESSURE_HEALTHY_ANCHOR_MODE,  # type: ignore[arg-type]
                                     resolve_arms=False, healthy_baseline=base)
        Pdev = Pdev.apply(pd.to_numeric, errors="coerce").fillna(0.0).loc[lambda d: d.abs().sum(axis=1) > 0]
        pos = Pdev.clip(lower=0.0); pos.index = [f"{g}|gain" for g in Pdev.index]
        neg = (-Pdev).clip(lower=0.0); neg.index = [f"{g}|loss" for g in Pdev.index]
        M = pd.concat([pos, neg]).loc[lambda d: d.sum(axis=1) > 0]
        fit_s = _fit_nmf(M, "S", NMF_K)
        if fit_s:
            fit_s.labels = _signed_labels(fit_s.W, hs)
            fit_s.cohort_label = "all tumors (signed gain+loss)"
    if fit_s:
        summaries.append(_run_one("gene_pressure_signed_cohort", fit_s, col_pam_all, out_root))

    # --- mirna within-tumor cohort ---
    expr = tumor.apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[expr.median(axis=1) >= MIN_EXPR].dropna(how="any").clip(lower=0.0)
    fit_t = _fit_nmf(expr, "T", NMF_K)
    if fit_t:
        fit_t.labels = _arm_labels(fit_t.W, mir_class)
        fit_t.cohort_label = "all tumors (arm abundance)"
        summaries.append(_run_one("mirna_within_tumor_cohort", fit_t, col_pam_all, out_root))

    # --- per-subtype: gene-pressure + arm within-subtype ---
    for sub in _PAM50:
        cols = [c for c in col_pam_all.index[col_pam_all == sub] if c in P.columns]
        if len(cols) < 30:
            continue
        col_sub = col_pam_all.reindex(cols)

        fit_g = _fit_nmf(P[cols], f"{sub}_G", NMF_K)
        if fit_g:
            fit_g.labels = _gene_labels(fit_g.W, hs, ec)
            fit_g.cohort_label = f"{sub} only"
            summaries.append(_run_one(f"gene_pressure_within_{sub}", fit_g, None, out_root))

        if acquired_within:
            cols_a = [c for c in col_pam_all.index[col_pam_all == sub] if c in Pa.columns]
            if len(cols_a) >= 30:
                fit_ga = _fit_nmf(Pa[cols_a], f"{sub}_Ga", NMF_K)
                if fit_ga:
                    fit_ga.labels = _gene_labels(fit_ga.W, hs, ec)
                    fit_ga.cohort_label = f"{sub} only (acquired)"
                    summaries.append(_run_one(f"gene_pressure_acquired_within_{sub}", fit_ga, None, out_root))

        cols_e = [c for c in col_pam_all.index[col_pam_all == sub] if c in expr.columns]
        if len(cols_e) >= 30:
            fit_ms = _fit_nmf(expr[cols_e], f"{sub}_S", K_SUBTYPE_ARM)
            if fit_ms:
                fit_ms.labels = _arm_labels(fit_ms.W, mir_class)
                fit_ms.cohort_label = f"{sub} only (arm abundance)"
                summaries.append(_run_one(f"mirna_within_{sub}", fit_ms, None, out_root))

    sm = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    if not sm.empty:
        sm.to_csv(out_root / "summary_all.tsv", sep="\t", index=False)
    return sm


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-acquired-within", action="store_true",
                    help="skip per-subtype acquired gene-pressure analyses")
    ap.add_argument("--out-root", type=Path, default=OUT_ROOT)
    args = ap.parse_args()
    run(acquired_within=not args.no_acquired_within, out_root=args.out_root)


if __name__ == "__main__":
    main()
