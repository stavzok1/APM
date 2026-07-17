"""Within-pathway NMF — decompose miRNA gene-pressure INSIDE a single Hallmark set.

The cross-set gene-pressure NMF (`mirna_comovement`) factorises the whole ~1,410-gene universe and
labels each factor by its dominant Hallmark. The within-pathway view instead restricts the
gene x sample pressure matrix to ONE set's genes and decomposes *that*, surfacing co-pressure
sub-modules INSIDE a program (e.g. within E2F: a master-TF module vs an effector module) that the
cross-set NMF averages away — it competes whole programs against each other, so sub-program
structure is lost.

FIVE pressure flavours (each a different question; within a single set the gain backbone does NOT
swamp the fit the way it does cohort-wide, so all are informative):
  - abundance (softmax_logrpm)             : attribution / steady-state pressure share
  - spine     (softmax_z_logrpm)           : the CANONICAL M0 pressure the whole subproject couples on
  - hybrid    (softmax_z + combined_mass)  : the M7 hybrid pressure (relative share x z)
  - acquired  (softmax_devhealthy, clip>=0): which genes co-GAIN pressure vs healthy
  - signed    (softmax_devhealthy, +/-)    : gain AND relief in one decomposition (direction kept)

Two SCOPES per (set, flavour):
  - cohort                : all tumours
  - within-subtype        : NMF re-fit on each PAM50 column-slice (LumA/LumB/Her2/Basal)

A SECOND, genuinely-different NMF family (not pressure multipliers — different *matrices*, no sample
axis; the within-pathway analogue of the report's `_nmf_programs` / `nmf_share_*`): edge-level
miRNA x gene lenses restricted to the set's genes, each factor pairing an ARM module with a GENE module:
  - aggregate   (c_tumor)                 : co-regulating arm module + co-pressured gene module
  - share       (gene_share_tumor)        : which arms OWN the program's gene budget
  - grip_change (share_tumor - share_gtex): which arms co-GAIN vs co-LOSE regulatory grip (signed)

Per-set k is TUNED by a reconstruction-error elbow (k in [2, kmax]); small sets get fewer factors.

Outputs (under output/tissue_reference/within_pathway_nmf/):
  factor_summary.tsv  — (hallmark, mode, scope, factor): label, top genes, driver arms, per-PAM50
                        activation + KW q, gain_frac (signed), k_chosen
  set_dominance.tsv   — (hallmark, mode, scope): n_genes, k, frac_dominant, median top-share,
                        median eff #signatures, leading factor (is the structure discrete or a blend?)
  per_sample.tsv.gz   — (hallmark, mode, scope, sample): dominant_factor, dominant_share, eff_n_sigs, call
  figures/*.png       — heatmaps for the focal sets (cohort scope, abundance+spine)

Run:  .venv/bin/python3 -m mirna_hallmark.analyses.nmf.within_pathway_nmf
        [--sets E2F_TARGETS,EMT] [--flavours abundance,spine,...] [--kmax 6] [--min-genes 40]
        [--no-within-subtype] [--no-heatmaps]
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.analyses.misc.mirna_comovement import _stable_factor_label

OUT_DIR = C.TISSUE_REFERENCE_DIR / "within_pathway_nmf"
FIG_DIR = OUT_DIR / "figures"
KMAX_DEFAULT = 6       # upper bound for the per-set elbow
KMIN = 2
MIN_GENES = 40         # skip sets with fewer targetable genes
MIN_SUBTYPE_N = 40     # skip a within-subtype fit with fewer samples
ELBOW_GAIN = 0.02      # stop adding factors when marginal recon-error drop < 2%
TOP_GENES = 8
TOP_ARMS = 8
DOMINANT_SHARE = 0.50
_PAM50_SCOPES = ("LumA", "LumB", "Her2", "Basal")

# flavour -> (expr_mode, target_norm, channel, needs_dev_baseline)
FLAVOURS: Dict[str, Tuple[str, str, str, bool]] = {
    "abundance": ("softmax_logrpm", "evidence_mass", "pos", False),
    "spine": ("softmax_z_logrpm", "evidence_mass", "pos", False),
    "hybrid": ("softmax_z", "combined_mass", "pos", False),
    "acquired": ("softmax_devhealthy_logrpm", "evidence_mass", "pos", True),
    "signed": ("softmax_devhealthy_logrpm", "evidence_mass", "signed", True),
}

FOCAL = {
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY",
    "HALLMARK_HYPOXIA", "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE",
}


def _pressure_for_flavours(hs: HallmarkSets, flavours: List[str]) -> Tuple[Dict[str, pd.DataFrame], pd.Series]:
    """Full-universe gene-pressure (gene x sample, UNCLIPPED) per *distinct* (expr_mode, target_norm),
    computed once; per-set decompositions are row/column slices. Returns ({flavour: P}, col_pam)."""
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    from mirna_hallmark.healthy_anchor import load_baseline

    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    tumor = _state_matrices()["tumor"]
    base = load_baseline()

    cache: Dict[Tuple[str, str, bool], pd.DataFrame] = {}
    out: Dict[str, pd.DataFrame] = {}
    for fl in flavours:
        expr_mode, tnorm, _chan, dev = FLAVOURS[fl]
        key = (expr_mode, tnorm, dev)
        if key not in cache:
            kw = {"healthy_baseline": base} if dev else {}
            P = compute_gene_pressure(genes, edges=edges, mirna=tumor, expr_mode=expr_mode,
                                      target_norm=tnorm, resolve_arms=False, **kw)
            cache[key] = P.apply(pd.to_numeric, errors="coerce").fillna(0.0)
            print(f"[within-pathway] computed pressure: expr_mode={expr_mode} target_norm={tnorm} "
                  f"dev={dev}  shape={cache[key].shape}")
        out[fl] = cache[key]

    clin = D.load_clinical_strata()
    part_pam = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
                if "PAM50_final" in clin.columns else {})
    any_P = next(iter(out.values()))
    col_pam = pd.Series({c: part_pam.get(_participant(c)) for c in any_P.columns})
    return out, col_pam


def _build_matrix(P_sub: pd.DataFrame, members: List[str], channel: str) -> pd.DataFrame:
    if channel == "signed":
        pos = P_sub.clip(lower=0.0); pos.index = [f"{g}|gain" for g in members]
        neg = (-P_sub).clip(lower=0.0); neg.index = [f"{g}|loss" for g in members]
        M = pd.concat([pos, neg])
    else:
        M = P_sub.clip(lower=0.0)
    return M.loc[M.sum(axis=1) > 0]


def _fit_k(M: pd.DataFrame, k: int, max_iter: int = 1000) -> Tuple[pd.DataFrame, pd.DataFrame, float]:
    from sklearn.decomposition import NMF
    model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=max_iter)
    W = pd.DataFrame(model.fit_transform(M.values), index=M.index,
                     columns=[f"f{j + 1}" for j in range(k)])
    H = pd.DataFrame(model.components_, index=W.columns, columns=M.columns)
    return W, H, float(model.reconstruction_err_)


def _choose_k(M: pd.DataFrame, kmax: int) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """Reconstruction-error elbow: increase k until the marginal relative-error drop < ELBOW_GAIN."""
    kmax_eff = max(KMIN, min(kmax, M.shape[0] // 8, M.shape[1] // 8))
    base_norm = float(np.linalg.norm(M.values)) or 1.0
    prev_rel = None
    best: Optional[Tuple[pd.DataFrame, pd.DataFrame, int]] = None
    for k in range(KMIN, kmax_eff + 1):
        W, H, err = _fit_k(M, k, max_iter=400)
        rel = err / base_norm
        if best is None:
            best = (W, H, k)
        elif prev_rel is not None and (prev_rel - rel) < ELBOW_GAIN:
            break
        else:
            best = (W, H, k)
        prev_rel = rel
    W, H, k = best  # type: ignore[misc]
    W, H, _ = _fit_k(M, k, max_iter=1000)  # refit chosen k to convergence
    return W, H, k


def _summarise(hset: str, mode: str, scope: str, W: pd.DataFrame, H: pd.DataFrame,
               col_pam: pd.Series, ec: Optional[pd.Series]) -> Tuple[List[dict], dict, pd.DataFrame]:
    from scipy.stats import kruskal
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _PAM50

    n_genes = int(W.shape[0] // (2 if mode == "signed" else 1))
    rows: List[dict] = []
    for f in W.columns:
        load = W[f].sort_values(ascending=False)
        if mode == "signed":
            gmass = float(W.loc[[i for i in W.index if i.endswith("|gain")], f].sum())
            lmass = float(W.loc[[i for i in W.index if i.endswith("|loss")], f].sum())
            gain_frac = round(gmass / (gmass + lmass + 1e-9), 3)
            top_genes = [i.replace("|gain", "").replace("|loss", "") for i in load.index][:TOP_GENES]
        else:
            gain_frac = ""
            top_genes = list(load.index[:TOP_GENES])
        a = H.loc[f].astype(float)
        means, groups = [], []
        row: dict = {"hallmark": hset.replace("HALLMARK_", ""), "mode": mode, "scope": scope,
                     "factor": f, "k_chosen": int(W.shape[1]), "n_genes": n_genes,
                     "focal": hset in FOCAL, "gain_frac": gain_frac}
        for sub in _PAM50:
            cols = [c for c in col_pam.index[col_pam == sub] if c in a.index]
            vals = a.loc[cols].dropna()
            means.append(float(vals.mean()) if len(vals) else np.nan)
            if len(vals) >= 10:
                groups.append(vals.values)
        arr = np.array(means, float)
        for i, sub in enumerate(_PAM50):
            row[f"frac_{sub}"] = round(arr[i] / np.nansum(arr), 3) if np.nansum(arr) > 0 and np.isfinite(arr[i]) else np.nan
        row["dominant_subtype"] = _PAM50[int(np.nanargmax(arr))] if np.nansum(arr) > 0 else "NA"
        try:
            row["kw_p"] = float(kruskal(*groups)[1]) if len(groups) >= 2 else np.nan
        except ValueError:
            row["kw_p"] = np.nan
        row["top_genes"] = ";".join(top_genes)
        if ec is not None:
            sub_e = ec[ec.index.get_level_values("gene").isin(top_genes)]
            row["top_driver_mirnas"] = ";".join(
                sub_e.groupby(level="miRNA").sum().sort_values(ascending=False).head(TOP_ARMS).index)
        else:
            row["top_driver_mirnas"] = ""
        row["factor_label"] = _stable_factor_label(
            {"dominant_subtype": row["dominant_subtype"], "dominant_hallmark": hset, "top_genes": row["top_genes"]})
        rows.append(row)

    comp = H.div(H.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    top_share = comp.max(axis=0)
    p = comp.where(comp > 0)
    ent = -(p * np.log(p)).sum(axis=0)
    dom_factor = comp.idxmax(axis=0)
    lead = dom_factor.value_counts().idxmax()
    dom = {
        "hallmark": hset.replace("HALLMARK_", ""), "mode": mode, "scope": scope,
        "focal": hset in FOCAL, "n_genes": n_genes, "k": int(W.shape[1]), "n_samples": int(H.shape[1]),
        "frac_dominant": round(float((top_share >= DOMINANT_SHARE).mean()), 3),
        "median_top_share": round(float(top_share.median()), 3),
        "median_eff_n_sigs": round(float(np.exp(ent).median()), 2),
        "leading_factor": lead,
        "leading_factor_label": next(r["factor_label"] for r in rows if r["factor"] == lead),
        "leading_factor_n": int((dom_factor == lead).sum()),
    }
    persamp = pd.DataFrame({
        "hallmark": hset.replace("HALLMARK_", ""), "mode": mode, "scope": scope,
        "sample": comp.columns, "PAM50": [col_pam.get(c, "") for c in comp.columns],
        "dominant_factor": dom_factor.values, "dominant_share": top_share.round(3).values,
        "eff_n_sigs": np.exp(ent).round(2).values,
        "call": np.where(top_share.values >= DOMINANT_SHARE, "dominant", "mixed"),
    })
    return rows, dom, persamp


def _heatmap(comp: pd.DataFrame, samp: pd.DataFrame, labels: Dict[str, str], title: str, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order: List[str] = []
    for sub in _PAM50_SCOPES:
        sc = samp[samp["PAM50"] == sub]
        order += sc.sort_values(["dominant_factor", "dominant_share"], ascending=[True, False])["sample"].tolist()
    order += [s for s in samp["sample"] if s not in order]
    dom_counts = samp["dominant_factor"].value_counts()
    row_order = sorted(comp.index, key=lambda f: -dom_counts.get(f, 0))
    M = comp.loc[row_order, order]
    fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.035 + 3), 0.5 * len(row_order) + 2))
    im = ax.imshow(M.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_yticks(range(len(row_order)))
    ax.set_yticklabels([f"{f}  {labels.get(f, '')[:32]}" for f in row_order], fontsize=8)
    ax.set_xticks([]); ax.set_xlabel(f"{len(order)} tumours (PAM50-ordered)")
    ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="signature share")
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def run(*, sets: Optional[List[str]] = None, flavours: Optional[List[str]] = None,
        kmax: int = KMAX_DEFAULT, min_genes: int = MIN_GENES, within_subtype: bool = True,
        heatmaps: bool = True, out_dir: Path = OUT_DIR) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    if heatmaps:
        FIG_DIR.mkdir(parents=True, exist_ok=True)
    flavours = flavours or list(FLAVOURS)
    hs = HallmarkSets.load()
    from mirna_hallmark.analyses.misc.mirna_comovement import _load_state_class, STATE_DIR
    _, edge = _load_state_class(STATE_DIR)
    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None

    print("[within-pathway] computing full-universe pressure matrices ...")
    Pmap, col_pam = _pressure_for_flavours(hs, flavours)

    targets = sets or sorted(hs.sets)
    factor_rows: List[dict] = []
    dom_rows: List[dict] = []
    persamp_parts: List[pd.DataFrame] = []
    for hset in targets:
        members0 = [g for g in hs.sets.get(hset, []) if g in next(iter(Pmap.values())).index]
        if len(members0) < min_genes:
            continue
        scopes = [("cohort", None)]
        if within_subtype:
            scopes += [(sub, sub) for sub in _PAM50_SCOPES]
        for fl in flavours:
            _, _, channel, _ = FLAVOURS[fl]
            P = Pmap[fl]
            members = [g for g in members0 if g in P.index]
            for scope_name, sub in scopes:
                cols = list(P.columns) if sub is None else [c for c in P.columns if col_pam.get(c) == sub]
                if sub is not None and len(cols) < MIN_SUBTYPE_N:
                    continue
                M = _build_matrix(P.loc[members, cols], members, channel)
                if min(M.shape) < KMIN:
                    continue
                W, H, k = _choose_k(M, kmax)
                frows, drow, ps = _summarise(hset, fl, scope_name, W, H, col_pam, ec)
                factor_rows.extend(frows); dom_rows.append(drow); persamp_parts.append(ps)
                if heatmaps and hset in FOCAL and scope_name == "cohort" and fl in ("abundance", "spine"):
                    comp = H.div(H.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
                    labels = {r["factor"]: r["factor_label"] for r in frows}
                    _heatmap(comp, ps, labels, f"{hset.replace('HALLMARK_', '')} — {fl} (cohort, k={k})",
                             FIG_DIR / f"{hset.replace('HALLMARK_', '')}_{fl}_cohort.png")
        print(f"[within-pathway] {hset.replace('HALLMARK_', ''):<34} n_genes={len(members0)}")

    fs = pd.DataFrame(factor_rows)
    if not fs.empty:
        fs["kw_q"] = _bh(fs["kw_p"].fillna(1.0).values)
    dm = pd.DataFrame(dom_rows)
    ps_all = pd.concat(persamp_parts, ignore_index=True) if persamp_parts else pd.DataFrame()
    fs.to_csv(out_dir / "factor_summary.tsv", sep="\t", index=False)
    dm.to_csv(out_dir / "set_dominance.tsv", sep="\t", index=False)
    ps_all.to_csv(out_dir / "per_sample.tsv.gz", sep="\t", index=False, compression="gzip")

    # second family: edge-level (miRNA x gene) lenses — different matrices, not pressure variants
    print("[within-pathway] edge-level lenses (aggregate / share / grip_change) ...")
    lens = _edge_lens_nmf(edge, hs, targets, min_genes, kmax)
    lens.to_csv(out_dir / "lens_factors.tsv", sep="\t", index=False)
    print(f"[within-pathway] wrote {len(fs)} pressure-factor rows / {len(dm)} dominance rows / "
          f"{len(ps_all)} per-sample rows / {len(lens)} edge-lens rows -> {out_dir}")
    return {"factor_summary": fs, "set_dominance": dm, "per_sample": ps_all, "lens_factors": lens}


def _bh(p: np.ndarray) -> np.ndarray:
    from mirna_hallmark.stats import bh_fdr
    return bh_fdr(p)


# --------------------------------------------------------------------------- #
# Edge-level (miRNA x gene) lenses — genuinely DIFFERENT matrices, not pressure
# multipliers. These have no sample axis (cohort-structural): each factor pairs an
# ARM module with a GENE module. Restricting the edge table to gene in the set gives the
# within-pathway version of the report's `_nmf_programs` / `nmf_share_{dominant,grip_change}`.
#   aggregate   : c_tumor              -> co-regulating arm module + co-pressured gene module
#   share       : gene_share_tumor     -> which arms OWN the program's gene budget
#   grip_change : gene_share_t - _gtex -> which arms co-GAIN vs co-LOSE grip (signed)
# --------------------------------------------------------------------------- #
def _edge_lens_nmf(edge: pd.DataFrame, hs: HallmarkSets, sets: List[str],
                   min_genes: int, kmax: int) -> pd.DataFrame:
    rows: List[dict] = []
    has_share = {"gene_share_tumor", "gene_share_gtex"}.issubset(edge.columns)
    for hset in sets:
        members = set(hs.sets.get(hset, []))
        e = edge[edge["gene"].isin(members)]
        if e["gene"].nunique() < min_genes or e["miRNA"].nunique() < KMIN:
            continue
        short = hset.replace("HALLMARK_", "")
        lenses: List[Tuple[str, pd.DataFrame, str]] = []
        if "c_tumor" in e.columns:
            lenses.append(("aggregate", _pivot_edge(e, "c_tumor"), "pos"))
        if "gene_share_tumor" in e.columns:
            lenses.append(("share", _pivot_edge(e, "gene_share_tumor"), "pos"))
        if has_share:
            d = e.assign(d_share=pd.to_numeric(e["gene_share_tumor"], errors="coerce")
                         - pd.to_numeric(e["gene_share_gtex"], errors="coerce"))
            gain = d.pivot_table(index="miRNA", columns="gene", values="d_share",
                                 aggfunc=lambda s: np.clip(s, 0, None).mean()).fillna(0.0)
            loss = d.pivot_table(index="miRNA", columns="gene", values="d_share",
                                 aggfunc=lambda s: np.clip(-s, 0, None).mean()).fillna(0.0)
            gain.columns = [f"{c}|gain" for c in gain.columns]
            loss.columns = [f"{c}|loss" for c in loss.columns]
            lenses.append(("grip_change", pd.concat([gain, loss], axis=1).fillna(0.0), "signed"))

        for lens, M, channel in lenses:
            M = M.loc[M.sum(axis=1) > 0, M.sum(axis=0) > 0]
            if min(M.shape) < KMIN:
                continue
            W, H, k = _choose_k(M, kmax)            # W = arm x factor, H = factor x gene(channel)
            for f in W.columns:
                top_arms = ";".join(W[f].sort_values(ascending=False).head(TOP_ARMS).index)
                gload = H.loc[f].sort_values(ascending=False)
                if channel == "signed":
                    gmass = float(H.loc[f, [c for c in H.columns if c.endswith("|gain")]].sum())
                    lmass = float(H.loc[f, [c for c in H.columns if c.endswith("|loss")]].sum())
                    direction = "gain" if gmass >= lmass else "loss"
                    gain_frac = round(gmass / (gmass + lmass + 1e-9), 3)
                    top_genes = ";".join(g.rsplit("|", 1)[0] for g in gload.head(TOP_GENES).index)
                else:
                    direction, gain_frac = "", ""
                    top_genes = ";".join(gload.head(TOP_GENES).index)
                rows.append({"hallmark": short, "lens": lens, "factor": f, "k_chosen": int(W.shape[1]),
                             "n_genes": int(e["gene"].nunique()), "n_arms": int(e["miRNA"].nunique()),
                             "focal": hset in FOCAL, "direction": direction, "gain_frac": gain_frac,
                             "top_mirnas": top_arms, "top_genes": top_genes})
    return pd.DataFrame(rows)


def _pivot_edge(e: pd.DataFrame, value: str) -> pd.DataFrame:
    m = e.pivot_table(index="miRNA", columns="gene", values=value, aggfunc="mean")
    return m.fillna(0.0).clip(lower=0.0)


def main() -> None:
    ap = argparse.ArgumentParser(description="Within-pathway (single Hallmark set) gene-pressure NMF")
    ap.add_argument("--sets", type=str, default=None, help="comma-separated Hallmark set names (default: all)")
    ap.add_argument("--flavours", type=str, default=None,
                    help=f"comma-separated subset of {list(FLAVOURS)} (default: all)")
    ap.add_argument("--kmax", type=int, default=KMAX_DEFAULT)
    ap.add_argument("--min-genes", type=int, default=MIN_GENES)
    ap.add_argument("--no-within-subtype", action="store_true")
    ap.add_argument("--no-heatmaps", action="store_true")
    args = ap.parse_args()
    sets = None
    if args.sets:
        sets = [s if s.startswith("HALLMARK_") else f"HALLMARK_{s}" for s in args.sets.split(",")]
    flav = args.flavours.split(",") if args.flavours else None
    run(sets=sets, flavours=flav, kmax=args.kmax, min_genes=args.min_genes,
        within_subtype=not args.no_within_subtype, heatmaps=not args.no_heatmaps)


if __name__ == "__main__":
    main()
