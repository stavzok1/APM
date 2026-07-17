"""Advanced NMF lenses — beyond the gene×sample *pressure* decompositions of §6.

Five analyses, each on a DIFFERENT signal or decomposition structure (not another pressure variant):

  1. edge_sample  — NMF on the **edge(arm,gene) × sample** contribution matrix. Keeps arm AND gene AND
                    patient identity jointly (gene×sample sums over arms; miRNA×gene sums over samples).
                    A factor = a set of specific (arm→gene) edges that co-fire in a patient subset.
  2. coupling     — NMF on the **realized** edge × stratum coupling (partial-ρ), NOT pressure. Every NMF
                    elsewhere is on pressure (potential); this is on realized repression → co-repression
                    programs that turn on together across PAM50 / immune / stage strata.
  3. consensus    — cophenetic-stability NMF (subsample + connectivity) to pick a robust k and flag which
                    factors are stable vs fit artifacts (answers the "are the modules real" caveat).
  4. supervised   — factor ↔ OUTCOME association: refit the gene×sample pressure NMF, then test each
                    factor's per-sample activity against immune subtype / HRD / aneuploidy / stage →
                    which pressure programs track the immune-visibility / genomic-instability axes.
  5. cnv_locus    — NMF on the miRNA-locus **CNV** (sample × locus copy number) matrix → co-amplified /
                    co-deleted oncomiR-locus programs by subtype (a DNA-dosage signal, not pressure).

Outputs under output/tissue_reference/nmf_advanced/. Run:
  .venv/bin/python3 -m mirna_hallmark.nmf_advanced --which edge_sample,supervised,cnv_locus,coupling,consensus
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.analyses.nmf.within_pathway_nmf import _choose_k, _fit_k, KMIN

OUT_DIR = C.TISSUE_REFERENCE_DIR / "nmf_advanced"
TOPN = 10


def _pam50(cols: Sequence[str]) -> pd.Series:
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.family_normal_reference import _participant
    clin = D.load_clinical_strata()
    m = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
         if "PAM50_final" in clin.columns else {})
    return pd.Series({c: m.get(_participant(c)) for c in cols})


def _kw_by_pam50(activity: pd.Series, col_pam: pd.Series) -> Tuple[str, float, Dict[str, float]]:
    from scipy.stats import kruskal
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _PAM50
    means, groups = {}, []
    for sub in _PAM50:
        v = activity.loc[[c for c in col_pam.index[col_pam == sub] if c in activity.index]].dropna()
        means[sub] = float(v.mean()) if len(v) else np.nan
        if len(v) >= 10:
            groups.append(v.values)
    tot = np.nansum(list(means.values()))
    frac = {s: (means[s] / tot if tot > 0 and np.isfinite(means[s]) else np.nan) for s in _PAM50}
    dom = max(frac, key=lambda s: (frac[s] if np.isfinite(frac[s]) else -1))
    try:
        p = float(kruskal(*groups)[1]) if len(groups) >= 2 else np.nan
    except ValueError:
        p = np.nan
    return dom, p, frac


# --------------------------------------------------------------------------- #
# 1. edge × sample
# --------------------------------------------------------------------------- #
def run_edge_sample(*, expr_mode: str = "softmax_z_logrpm", k: Optional[int] = None,
                    kmax: int = 14, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
    from mirna_hallmark.pressure_build import load_mirtar_edges
    from mirna_hallmark.pressure_engine import compute_edge_pressure_map
    from mirna_hallmark.analyses.misc.mirna_comovement import _load_state_class, STATE_DIR

    hs = HallmarkSets.load()
    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    tumor = _state_matrices()["tumor"]
    print(f"[edge_sample] building edge×sample contribution map ({expr_mode}) ...")
    emap = compute_edge_pressure_map(edges, tumor, genes=genes, expr_mode=expr_mode,  # type: ignore[arg-type]
                                     target_norm="evidence_mass")
    M = pd.DataFrame({f"{m}|||{g}": s for (m, g), s in emap.items()}).T
    M = M.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    M = M.loc[M.sum(axis=1) > 0]
    print(f"[edge_sample] matrix {M.shape[0]} edges × {M.shape[1]} samples; fitting NMF ...")
    if k:
        W, H, _ = _fit_k(M, k)
    else:
        W, H, k = _choose_k(M, kmax)
    col_pam = _pam50(M.columns)
    _, edgetab = _load_state_class(STATE_DIR)
    ec = edgetab.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edgetab.columns else None

    rows = []
    for f in W.columns:
        top = W[f].sort_values(ascending=False).head(TOPN).index
        edge_pairs = [e.split("|||") for e in top]
        top_edges = ";".join(f"{a.replace('hsa-', '')}→{g}" for a, g in edge_pairs)
        top_arms = ";".join(pd.Series([a.replace("hsa-", "") for a, _ in edge_pairs]).drop_duplicates().head(6))
        top_genes = ";".join(pd.Series([g for _, g in edge_pairs]).drop_duplicates().head(8))
        dom, p, frac = _kw_by_pam50(H.loc[f], col_pam)
        rows.append({"factor": f, "k": int(W.shape[1]), "dominant_subtype": dom,
                     **{f"frac_{s}": round(v, 3) if np.isfinite(v) else np.nan for s, v in frac.items()},
                     "kw_p": p, "top_edges": top_edges, "top_arms": top_arms, "top_genes": top_genes})
    fac = pd.DataFrame(rows)
    fac["kw_q"] = _bh(fac["kw_p"].fillna(1.0).values)
    # per-sample dominance
    comp = H.div(H.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    ent = -(comp.where(comp > 0) * np.log(comp.where(comp > 0))).sum(axis=0)
    ps = pd.DataFrame({"sample": comp.columns, "PAM50": [col_pam.get(c, "") for c in comp.columns],
                       "dominant_factor": comp.idxmax(axis=0).values,
                       "dominant_share": comp.max(axis=0).round(3).values,
                       "eff_n_factors": np.exp(ent).round(2).values})
    out_dir.mkdir(parents=True, exist_ok=True)
    fac.to_csv(out_dir / "edge_sample_factors.tsv", sep="\t", index=False)
    ps.to_csv(out_dir / "edge_sample_per_sample.tsv.gz", sep="\t", index=False, compression="gzip")
    print(f"[edge_sample] k={k}; {len(fac)} factors; {(ps['dominant_share']>=0.5).mean():.0%} single-edge-module-dominant")
    return {"edge_sample_factors": fac, "edge_sample_per_sample": ps}


# --------------------------------------------------------------------------- #
# 2. coupling (realized, edge × stratum)
# --------------------------------------------------------------------------- #
def run_coupling(*, rho_floor: float = 0.10, min_stratum_n: int = 40, k: Optional[int] = None,
                 kmax: int = 10, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.analyses.misc.mirna_comovement import _load_state_class, STATE_DIR
    from mirna_hallmark.stats import correlation_pair

    _, edge = _load_state_class(STATE_DIR)
    # repressive-coupled edges only (keeps the matrix meaningful + tractable)
    if "adj_rho_tumor" in edge.columns:
        coupled = edge[pd.to_numeric(edge["adj_rho_tumor"], errors="coerce") <= -rho_floor]
    else:
        coupled = edge
    pairs = list(coupled[["miRNA", "gene"]].drop_duplicates().itertuples(index=False, name=None))
    tumor = _state_matrices()["tumor"]
    rna = D.load_rna()  # gene × sample
    clin = D.load_clinical_strata().set_index("participant")
    # strata panel: PAM50 ∪ immune ∪ stage
    strata: Dict[str, List[str]] = {}
    part_of = {c: _participant(c) for c in tumor.columns}
    for col, prefix in [("PAM50_final", "PAM50"), ("thornsson_immune_subtype", "IMM"),
                        ("pathologic_stage_collapsed", "STAGE")]:
        if col not in clin.columns:
            continue
        for lvl, g in clin.groupby(col):
            samples = [c for c in tumor.columns if part_of[c] in set(g.index)]
            if len(samples) >= min_stratum_n:
                strata[f"{prefix}:{lvl}"] = samples
    cov_all = clin[["CPE", "thornsson_hrd_score"]] if {"CPE", "thornsson_hrd_score"}.issubset(clin.columns) else None
    print(f"[coupling] {len(pairs)} repressive-coupled edges × {len(strata)} strata; computing partial-ρ ...")

    rna = rna[~rna.index.duplicated(keep="first")]      # collapse duplicate gene symbols
    tumor = tumor[~tumor.index.duplicated(keep="first")]
    M = pd.DataFrame(0.0, index=[f"{m}|||{g}" for m, g in pairs], columns=list(strata))
    for (m, g) in pairs:
        if m not in tumor.index or g not in rna.index:
            continue
        for sname, samples in strata.items():
            xs = [c for c in samples if c in tumor.columns]
            x = pd.to_numeric(tumor.loc[m, xs], errors="coerce")
            ycols = [c for c in xs if _participant(c) in rna.columns]
            if len(ycols) < min_stratum_n:
                continue
            y = pd.to_numeric(rna.loc[g, [_participant(c) for c in ycols]], errors="coerce")
            y.index = ycols
            cov = None
            if cov_all is not None:
                cov = cov_all.reindex([_participant(c) for c in ycols]); cov.index = ycols
            res = correlation_pair(y, x.reindex(ycols), covariates=cov)
            rho = res.get("partial_rho") if cov is not None else res.get("spearman_rho")
            if rho is None or not np.isfinite(rho):
                rho = res.get("spearman_rho")
            if rho is not None and np.isfinite(rho):
                M.loc[f"{m}|||{g}", sname] = max(-float(rho), 0.0)  # repression strength
    M = M.loc[M.sum(axis=1) > 0]
    if min(M.shape) < KMIN:
        print("[coupling] too sparse for NMF"); return {}
    print(f"[coupling] matrix {M.shape}; NMF ...")
    if k:
        W, H, _ = _fit_k(M, k)
    else:
        W, H, k = _choose_k(M, min(kmax, M.shape[1]))
    rows = []
    for f in W.columns:
        te = W[f].sort_values(ascending=False).head(TOPN).index
        ts = H.loc[f].sort_values(ascending=False).head(6).index
        rows.append({"factor": f, "k": int(W.shape[1]),
                     "top_strata": ";".join(ts),
                     "top_edges": ";".join(e.replace("|||", "→").replace("hsa-", "") for e in te)})
    fac = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    fac.to_csv(out_dir / "coupling_factors.tsv", sep="\t", index=False)
    M.to_csv(out_dir / "coupling_matrix.tsv.gz", sep="\t", compression="gzip")
    print(f"[coupling] k={k}; {len(fac)} realized-repression programs across strata")
    return {"coupling_factors": fac}


# --------------------------------------------------------------------------- #
# 3. consensus (cophenetic stability)
# --------------------------------------------------------------------------- #
def run_consensus(*, matrix: str = "gene_sample", krange: Sequence[int] = range(2, 11),
                  n_runs: int = 15, subsample: float = 0.8, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    from sklearn.decomposition import NMF
    from scipy.cluster.hierarchy import linkage, cophenet
    from scipy.spatial.distance import squareform

    M = _matrix_for_consensus(matrix)
    print(f"[consensus] matrix={matrix} {M.shape}; krange={list(krange)} × {n_runs} runs ...")
    rng = np.random.default_rng(0)
    rows = []
    ncol = M.shape[1]
    for kk in krange:
        consensus = np.zeros((ncol, ncol))
        counts = np.zeros((ncol, ncol))
        for _ in range(n_runs):
            idx = rng.choice(ncol, int(subsample * ncol), replace=False)
            sub = M.iloc[:, idx]
            sub = sub.loc[sub.sum(axis=1) > 0]
            if min(sub.shape) < kk:
                continue
            H = NMF(n_components=kk, init="nndsvda", random_state=int(rng.integers(1e6)),
                    max_iter=300).fit(sub.values).components_
            assign = H.argmax(axis=0)
            for a in range(len(idx)):
                same = idx[assign == assign[a]]
                consensus[np.ix_([idx[a]], same)] += 1
                counts[np.ix_([idx[a]], same)] += 1
        with np.errstate(invalid="ignore"):
            Cmat = np.divide(consensus, counts, out=np.zeros_like(consensus), where=counts > 0)
        keep = counts.sum(axis=1) > 0
        Cs = Cmat[np.ix_(keep, keep)]
        d = 1.0 - Cs
        np.fill_diagonal(d, 0.0)
        d = (d + d.T) / 2.0
        try:
            Z = linkage(squareform(d, checks=False), method="average")
            coph = float(cophenet(Z, squareform(d, checks=False))[0])
        except Exception:
            coph = np.nan
        disp = float(np.sum(4 * (Cs - 0.5) ** 2) / Cs.size)  # 1=perfectly bimodal/stable, 0=diffuse
        rows.append({"matrix": matrix, "k": kk, "cophenetic_corr": round(coph, 4),
                     "dispersion": round(disp, 4)})
        print(f"[consensus] k={kk}: cophenetic={coph:.3f} dispersion={disp:.3f}")
    df = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / f"consensus_{matrix}.tsv", sep="\t", index=False)
    return df


def _matrix_for_consensus(matrix: str) -> pd.DataFrame:
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices, ATTR_MODE
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    hs = HallmarkSets.load()
    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    tumor = _state_matrices()["tumor"]
    P = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode="softmax_z_logrpm",  # type: ignore[arg-type]
                              resolve_arms=False)
    return P.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0).loc[lambda d: d.sum(axis=1) > 0]


# --------------------------------------------------------------------------- #
# 4. supervised (factor ↔ outcome)
# --------------------------------------------------------------------------- #
def run_supervised(*, k: int = 12, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    from scipy.stats import spearmanr, kruskal
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges

    hs = HallmarkSets.load()
    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    tumor = _state_matrices()["tumor"]
    P = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode="softmax_z_logrpm",  # type: ignore[arg-type]
                              resolve_arms=False).apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    P = P.loc[P.sum(axis=1) > 0]
    W, H, _ = _fit_k(P, k)
    clin = D.load_clinical_strata().set_index("participant")
    part = {c: _participant(c) for c in H.columns}
    col_pam = _pam50(H.columns)

    def series(col):
        return pd.Series({c: clin[col].get(part[c]) for c in H.columns}) if col in clin.columns else None

    hrd = pd.to_numeric(series("thornsson_hrd_score"), errors="coerce")
    aneu = pd.to_numeric(series("thornsson_aneuploidy_score"), errors="coerce")
    imm = series("thornsson_immune_subtype")
    stage = series("pathologic_stage_collapsed")
    rows = []
    for f in H.index:
        a = H.loc[f].astype(float)
        top = ";".join(W[f].sort_values(ascending=False).head(8).index)
        dom, _, _ = _kw_by_pam50(a, col_pam)
        dom_samples = set(col_pam.index[col_pam == dom])
        row = {"factor": f, "dominant_subtype": dom, "top_genes": top}
        for name, v in [("HRD", hrd), ("aneuploidy", aneu)]:
            ok = v.dropna().index.intersection(a.index)
            rho, p = spearmanr(a.loc[ok], v.loc[ok]) if len(ok) > 20 else (np.nan, np.nan)
            row[f"{name}_rho"] = round(float(rho), 3); row[f"{name}_p"] = float(p)
            # within the factor's own dominant subtype (removes the PAM50↔HRD confound)
            okd = [c for c in ok if c in dom_samples]
            rd, pdv = spearmanr(a.loc[okd], v.loc[okd]) if len(okd) > 20 else (np.nan, np.nan)
            row[f"{name}_rho_within_dom"] = round(float(rd), 3) if np.isfinite(rd) else np.nan
            row[f"{name}_p_within_dom"] = float(pdv)
        for name, cat in [("immune", imm), ("stage", stage)]:
            if cat is None:
                row[f"{name}_kw_p"] = np.nan; continue
            groups = [a.loc[[c for c in g.index if c in a.index]].values
                      for _, g in cat.dropna().groupby(cat) if len(g) >= 15]
            try:
                row[f"{name}_kw_p"] = float(kruskal(*groups)[1]) if len(groups) >= 2 else np.nan
            except ValueError:
                row[f"{name}_kw_p"] = np.nan
        rows.append(row)
    df = pd.DataFrame(rows)
    for c in ["HRD_p", "aneuploidy_p", "immune_kw_p", "stage_kw_p"]:
        df[c.replace("_p", "_q")] = _bh(df[c].fillna(1.0).values)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "supervised_factor_outcome.tsv", sep="\t", index=False)
    print(f"[supervised] {len(df)} factors scored vs HRD/aneuploidy/immune/stage")
    return df


# --------------------------------------------------------------------------- #
# 5. cnv_locus
# --------------------------------------------------------------------------- #
def run_cnv_locus(*, k: Optional[int] = None, kmax: int = 12, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    path = C.OUTPUT_ROOT / "mirna_locus_cnv" / "tables" / "sample_entity_cnv.tsv.gz"
    if not path.exists():
        print(f"[cnv_locus] {path} missing — run mirna_locus_cnv first"); return {}
    df = pd.read_csv(path, sep="\t", usecols=["participant", "entity_id", "entity_level", "copy_number"])
    df = df[df["entity_level"].isin(["locus", "mirna_locus", "arm_locus"]) | True]  # keep all entity rows
    M = df.pivot_table(index="entity_id", columns="participant", values="copy_number", aggfunc="median")
    M = M.apply(pd.to_numeric, errors="coerce")
    # center on diploid 2 → deviation; split amp/del channels for signed structure
    dev = (M - 2.0).fillna(0.0)
    amp = dev.clip(lower=0.0); amp.index = [f"{i}|amp" for i in M.index]
    dele = (-dev).clip(lower=0.0); dele.index = [f"{i}|del" for i in M.index]
    X = pd.concat([amp, dele]).loc[lambda d: d.sum(axis=1) > 0]
    col_pam = _pam50(X.columns)
    print(f"[cnv_locus] {X.shape[0]} locus-channels × {X.shape[1]} participants; NMF ...")
    if k:
        W, H, _ = _fit_k(X, k)
    else:
        W, H, k = _choose_k(X, kmax)
    rows = []
    for f in W.columns:
        top = W[f].sort_values(ascending=False).head(TOPN).index
        amp_l = [i.replace("|amp", "") for i in top if i.endswith("|amp")][:6]
        del_l = [i.replace("|del", "") for i in top if i.endswith("|del")][:6]
        dom, p, frac = _kw_by_pam50(H.loc[f], col_pam)
        rows.append({"factor": f, "k": int(W.shape[1]), "dominant_subtype": dom,
                     **{f"frac_{s}": round(v, 3) if np.isfinite(v) else np.nan for s, v in frac.items()},
                     "kw_p": p, "top_amplified": ";".join(amp_l), "top_deleted": ";".join(del_l)})
    fac = pd.DataFrame(rows)
    fac["kw_q"] = _bh(fac["kw_p"].fillna(1.0).values)
    out_dir.mkdir(parents=True, exist_ok=True)
    fac.to_csv(out_dir / "cnv_locus_factors.tsv", sep="\t", index=False)
    print(f"[cnv_locus] k={k}; {len(fac)} co-CNV locus programs")
    return {"cnv_locus_factors": fac}


def _bh(p: np.ndarray) -> np.ndarray:
    from mirna_hallmark.stats import bh_fdr
    return bh_fdr(p)


def main() -> None:
    ap = argparse.ArgumentParser(description="Advanced NMF lenses")
    ap.add_argument("--which", type=str, default="edge_sample,supervised,cnv_locus,coupling,consensus")
    args = ap.parse_args()
    which = list(args.which.split(","))
    fns = {"edge_sample": run_edge_sample, "supervised": run_supervised, "cnv_locus": run_cnv_locus,
           "coupling": run_coupling, "consensus": run_consensus}
    for name in which:
        if name not in fns:
            continue
        try:
            fns[name]()
        except Exception as exc:  # noqa: BLE001 — isolate so one lens can't kill the rest
            import traceback
            print(f"[nmf_advanced] {name} FAILED: {type(exc).__name__}: {exc}")
            traceback.print_exc()


if __name__ == "__main__":
    main()
