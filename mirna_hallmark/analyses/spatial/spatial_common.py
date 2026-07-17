"""Roadmap Phase 3 (spatial): shared toolkit for multi-resolution spatial localization of the
miRNA-Hallmark programs + the spatial MH-56(d) program↔pressure DECOUPLING retest.

This module carries the machinery reused by every resolution layer (`spatial_mibi_anchor`,
`spatial_visium_wu`, `spatial_visium_public`, `spatial_xenium`):

  * `PROGRAM_SETS`              the 5 framework programs → MSigDB-Hallmark sets
  * `program_gene_sets`         resolve them against `HallmarkSets`
  * `score_programs`            per-unit z-metagene program score (spot / cell)
  * `deconvolve_matrix`         NNLS deconvolution of any genes×units matrix vs the Wu signature
                                (balanced markers — the same recipe as `brca_deconvolution`)
  * `compartment_fractions`     fine cell-type fractions → 5 coarse compartments (`COARSE`)
  * `framework_target_genes`    the headline BY-neg coupling target genes (genes under miRNA pressure)
  * `load_pressure_delta`       per-gene released/gained pressure (`gene_pressure_by_state.delta_tumor_gtex`)
  * `gene_level_pressure_retest`  **the spatial MH-56d engine** — gene-by-gene over TARGETED genes,
                                joined to roles + pressure delta, classified brake-release /
                                concordant-repression / discordant. Never markers/averages (the MH-56 lesson).
  * `morans_i`                  hand-rolled spatial autocorrelation (kNN weights, no squidpy)

Binding constraint (roadmap): no measured spatial miRNA for breast tumours. The measured side is
mRNA/protein (cell-resolved); the miRNA side is projected (route i) or inferred (route ii) — every
caller tags which.

No scanpy/squidpy/scanpy in the venv — loaders parse raw mtx/h5 with scipy/h5py and stream.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.brca_deconvolution import build_signature
from mirna_hallmark.analyses.misc.compartment_annotate_network import COARSE
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr, zscore_rows

# The 5 programs the roadmap names, each a union of MSigDB-Hallmark sets.
PROGRAM_SETS: Dict[str, List[str]] = {
    "proliferation": ["HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"],
    "hypoxia": ["HALLMARK_HYPOXIA"],
    "glycolysis": ["HALLMARK_GLYCOLYSIS"],
    "EMT": ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
    "immune": ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE"],
}

# Anchor framework-target genes whose MH-56d pressure call we re-test in situ at every resolution.
ANCHOR_GENES = ["SLC2A1", "FAP", "ERBB2", "ESR1", "CDH1", "VIM", "HIF1A"]

EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
PRESSURE_BY_STATE = C.TISSUE_REFERENCE_DIR / "cross_state_landscape" / "gene_pressure_by_state.tsv"


# --------------------------------------------------------------------------- programs
def program_gene_sets(hs: Optional[HallmarkSets] = None) -> Dict[str, List[str]]:
    """Resolve PROGRAM_SETS → {program: [genes]} via HallmarkSets (union of member sets)."""
    hs = hs or HallmarkSets.load()
    out: Dict[str, List[str]] = {}
    for prog, sets in PROGRAM_SETS.items():
        genes: set = set()
        for s in sets:
            genes |= set(hs.sets.get(s, []))
        out[prog] = sorted(genes)
    return out


def score_programs(expr_log: pd.DataFrame, program_genes: Dict[str, List[str]]) -> pd.DataFrame:
    """Per-unit program score = mean per-gene z (across units) over a program's present genes.

    `expr_log` is genes×units (log expression for RNA, arcsinh for protein). Returns units×programs.
    Same metagene recipe as `robustness_checks._metagene` (z rows, average) — keeps every program on a
    comparable z scale across spatial units.
    """
    z = zscore_rows(expr_log)
    cols = {}
    for prog, genes in program_genes.items():
        present = [g for g in genes if g in z.index]
        cols[prog] = z.loc[present].mean(axis=0) if present else pd.Series(np.nan, index=z.columns)
    return pd.DataFrame(cols)


# --------------------------------------------------------------------------- deconvolution
def deconvolve_matrix(B_linear: pd.DataFrame, sig: Optional[pd.DataFrame] = None,
                      *, n_sig_genes: int = 1500) -> pd.DataFrame:
    """NNLS-deconvolve a genes×units LINEAR matrix against the Wu celltype_major signature.

    Reuses `brca_deconvolution.build_signature()` and its **balanced** per-cell-type marker recipe (top
    fold-change genes per type, expressed above the type median) so sharp immune/stromal markers don't
    starve the epithelial fit. Returns units×celltype fractions (rows sum to 1).
    """
    sig = build_signature() if sig is None else sig
    shared = sig.index.intersection(B_linear.index)
    sig_s = sig.loc[shared]
    per_type = max(60, n_sig_genes // sig_s.shape[1])
    markers: set = set()
    for t in sig_s.columns:
        others = sig_s.drop(columns=t).mean(axis=1)
        fc = (sig_s[t] + 1.0) / (others + 1.0)
        cand = fc[sig_s[t] > sig_s[t].quantile(0.5)].sort_values(ascending=False)
        markers |= set(cand.head(per_type).index)
    markers = sorted(markers & set(B_linear.index))
    S = sig_s.loc[markers].to_numpy()
    Bm = B_linear.loc[markers].to_numpy()
    fr = np.zeros((Bm.shape[1], S.shape[1]))
    for j in range(Bm.shape[1]):
        x, _ = nnls(S, Bm[:, j])
        fr[j] = x / x.sum() if x.sum() > 0 else x
    return pd.DataFrame(fr, index=B_linear.columns, columns=sig.columns)


def compartment_fractions(celltype_fr: pd.DataFrame) -> pd.DataFrame:
    """Collapse fine celltype_major fractions (units×celltype) → units×coarse-compartment (`COARSE`)."""
    grp: Dict[str, List[str]] = {}
    for ct, comp in COARSE.items():
        grp.setdefault(comp, []).append(ct)
    out = {}
    for comp, cts in grp.items():
        cols = [c for c in cts if c in celltype_fr.columns]
        out[comp] = celltype_fr[cols].sum(axis=1) if cols else pd.Series(0.0, index=celltype_fr.index)
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- framework objects
def framework_target_genes() -> set:
    """Headline BY-neg coupling target genes — the genes carrying measurable miRNA pressure."""
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e["gene"] = e["key"].str.split(r"\|\|").str[1]
    neg = e[(e["rho"] < 0) & (e["q_by"] < 0.05)]
    return set(neg["gene"].dropna().unique())


def load_pressure_delta() -> pd.Series:
    """gene → delta_tumor_gtex (released pressure if <0 = brake-release; gained if >0). The MH-56d object."""
    d = pd.read_csv(PRESSURE_BY_STATE, sep="\t")
    return d.set_index("gene")["delta_tumor_gtex"]


# --------------------------------------------------------------------------- the spatial MH-56d engine
def gene_level_pressure_retest(readout: pd.DataFrame, tumour_score: pd.Series, *,
                               target_genes: Optional[Iterable[str]] = None,
                               min_units: int = 30) -> pd.DataFrame:
    """Spatial re-test of the MH-56(d) program↔pressure decoupling — **gene-by-gene, targeted genes only**.

    For each framework-target gene with a measured spatial readout, correlate the gene's per-unit readout
    against a per-unit *tumour-content* score (epithelial fraction for spots; tumour-niche proximity /
    tumour-lineage for single cells). Compare the sign of that spatial association to the gene's bulk
    pressure change `delta_tumor_gtex`:

        delta < 0 (pressure RELEASED) & rho > 0  → brake_release          (gene rises in tumour as brake lifts; concordant w/ pressure model)
        delta > 0 (pressure GAINED)   & rho < 0  → concordant_repression  (gene falls where pressure gained)
        delta > 0                     & rho > 0  → discordant_rise        (gene rises in tumour DESPITE gained pressure = decoupling / overridden brake)
        delta < 0                     & rho < 0  → discordant_loss

    `readout` genes×units (log RNA / arcsinh protein); `tumour_score` per unit. Returns one row per gene.
    """
    units = readout.columns.intersection(tumour_score.dropna().index)
    if len(units) < min_units:
        return pd.DataFrame()
    ts = tumour_score.reindex(units).astype(float)
    delta = load_pressure_delta()
    # The decoupling universe is the genes carrying a measured pressure change (gene_pressure_by_state) —
    # this includes RELEASED genes (e.g. SLC2A1, whose coupling is lost in tumour so it is NOT a current
    # BY-neg edge), which the brake-release arm of the test depends on. Callers may narrow via target_genes.
    fw = set(target_genes) if target_genes is not None else set(delta.index)
    genes = [g for g in readout.index if g in fw and g in delta.index]
    from mirna_hallmark.gene_roles import load_gene_roles
    roles = load_gene_roles(genes).set_index("gene")

    rows = []
    for g in genes:
        v = readout.loc[g, units].astype(float)
        if v.notna().sum() < min_units or v.std() == 0:
            continue
        rho, p = spearmanr(v, ts, nan_policy="omit")
        if not np.isfinite(rho):
            continue
        d = float(delta[g])
        released = d < 0
        if released and rho > 0:
            cls = "brake_release"
        elif (not released) and rho < 0:
            cls = "concordant_repression"
        elif (not released) and rho > 0:
            cls = "discordant_rise"
        else:
            cls = "discordant_loss"
        rows.append({"gene": g, "delta_tumor_gtex": round(d, 3),
                     "pressure": "released" if released else "gained",
                     "spatial_tumour_rho": round(float(rho), 3), "p": float(p),
                     "classification": cls, "n_units": int(v.notna().sum()),
                     "role": roles.loc[g, "role"] if g in roles.index else "unknown",
                     "is_anchor": g in ANCHOR_GENES})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q"] = bh_fdr(out["p"].to_numpy())
        out = out.sort_values(["is_anchor", "p"], ascending=[False, True]).reset_index(drop=True)
    return out


# --------------------------------------------------------------------------- 10x loaders (no scanpy)
def _sniff_compression(path: Path) -> Optional[str]:
    """'gzip' if the file really is gzipped (magic bytes), else None — tolerates `.gz`-named plain text
    (the Zenodo export quirk). pandas' own `compression='infer'` trusts the extension and would crash."""
    with open(path, "rb") as fh:
        return "gzip" if fh.read(2) == b"\x1f\x8b" else None


def _read_features(path: Path) -> List[str]:
    """Gene symbols from a 10x features/genes tsv (plain or gzip; 1-col=symbol, 3-col=id/symbol/type)."""
    df = pd.read_csv(path, sep="\t", header=None, compression=_sniff_compression(path))
    return (df[1] if df.shape[1] > 1 else df[0]).astype(str).tolist()


def read_10x_mtx(mtx_dir: Path) -> pd.DataFrame:
    """Read a 10x mtx triplet (matrix.mtx[.gz], features, barcodes) → genes×cells raw counts.

    Tolerates the Zenodo/10x quirk where `.gz` files are actually plain text (pandas `compression='infer'`
    sniffs the magic bytes, not the extension). Built on scipy.io.mmread; never densifies beyond the result.
    """
    from scipy.io import mmread
    import gzip as _gz

    def _open(p: Path):
        with open(p, "rb") as fh:
            magic = fh.read(2)
        return _gz.open(p, "rt") if magic == b"\x1f\x8b" else open(p, "rt")

    mtx = next(mtx_dir.glob("matrix.mtx*"))
    genes = _read_features(next(mtx_dir.glob("features.tsv*")) if list(mtx_dir.glob("features.tsv*"))
                           else next(mtx_dir.glob("genes.tsv*")))
    bpath = next(mtx_dir.glob("barcodes.tsv*"))
    barcodes = pd.read_csv(bpath, sep="\t", header=None,
                           compression=_sniff_compression(bpath))[0].astype(str).tolist()
    with _open(mtx) as fh:
        m = mmread(fh).tocsr()                                  # genes × cells
    df = pd.DataFrame.sparse.from_spmatrix(m, index=genes, columns=barcodes)
    return df[~df.index.duplicated(keep="first")]


def read_10x_h5(h5_path: Path) -> pd.DataFrame:
    """Read a 10x filtered_feature_bc_matrix.h5 → genes×spots raw counts (CSC under matrix/)."""
    import h5py
    from scipy.sparse import csc_matrix

    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        data, indices, indptr = g["data"][:], g["indices"][:], g["indptr"][:]
        shape = tuple(g["shape"][:])
        names = [x.decode() for x in g["features/name"][:]]
        barcodes = [x.decode() for x in g["barcodes"][:]]
    m = csc_matrix((data, indices, indptr), shape=shape).tocsr()    # genes × spots
    df = pd.DataFrame.sparse.from_spmatrix(m, index=names, columns=barcodes)
    return df[~df.index.duplicated(keep="first")]


def cpm_log(counts: pd.DataFrame) -> pd.DataFrame:
    """genes×units raw counts → log2(CPM/100 + 1) (dense). Keeps the per-gene z-metagene scoring stable."""
    dense = counts.sparse.to_dense() if hasattr(counts, "sparse") else counts.astype(float)
    cpm = dense.div(dense.sum(axis=0).replace(0, np.nan), axis=1) * 1e4
    return np.log2(cpm + 1.0)


def read_visium_positions(path: Path) -> pd.DataFrame:
    """tissue_positions_list.csv → spots×(in_tissue,x,y) using pixel coords (handles header/no-header)."""
    head = pd.read_csv(path, nrows=1, header=None)
    has_header = str(head.iloc[0, 0]).lower() in ("barcode",)
    cols = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    df = pd.read_csv(path, header=0 if has_header else None)
    df.columns = cols[:df.shape[1]]
    df = df.set_index("barcode")
    return pd.DataFrame({"in_tissue": df["in_tissue"], "x": df["pxl_col"], "y": df["pxl_row"]})


# --------------------------------------------------------------------------- shared spot analysis
def analyse_spots(log: pd.DataFrame, lin: pd.DataFrame, coords: pd.DataFrame,
                  sig: pd.DataFrame, program_genes: Dict[str, List[str]]) -> dict:
    """Per-section spot analysis shared by every Visium layer: program scores → compartment localization,
    Moran's I niche structure, and the brake-release co-localization of the anchor target genes.

    `log`/`lin` are genes×spots (log and ~linear); `coords` spots×(x,y). Returns a JSON-able dict.
    """
    prog = score_programs(log, program_genes)                  # spots × 5 programs
    cfr = deconvolve_matrix(lin, sig=sig)                      # spots × celltype
    comp = compartment_fractions(cfr)                         # spots × compartment
    delta = load_pressure_delta()

    loc = {p: {c: round(float(spearmanr(prog[p], comp[c])[0]), 3) for c in comp.columns} for p in prog.columns}
    moran = {p: round(morans_i(prog[p], coords), 3) for p in prog.columns}
    epi = comp.get("Epithelial", pd.Series(0.0, index=log.columns))
    anchor = {}
    for g in ANCHOR_GENES:
        if g not in log.index:
            continue
        gv = log.loc[g]
        prog_rho = {p: float(spearmanr(gv, prog[p])[0]) for p in prog.columns}
        best = max(prog_rho, key=lambda k: abs(prog_rho[k]))
        anchor[g] = {"pressure": "released" if delta.get(g, 0) < 0 else "gained",
                     "best_program": best, "rho_best_program": round(prog_rho[best], 3),
                     "rho_epithelial_frac": round(float(spearmanr(gv, epi)[0]), 3),
                     "morans_i": round(morans_i(gv, coords), 3)}
    return {"n_spots": int(log.shape[1]), "localization": loc, "morans_i": moran, "anchor": anchor}


# --------------------------------------------------------------------------- malignancy (inferCNV-style)
_GENCODE_SLIM = C.REPO_ROOT / "data" / "gencode.v49.slim.parquet"


def load_gene_positions() -> pd.DataFrame:
    """gene_name → (chrom, start) from GENCODE (feature=='gene'), autosomes ordered by position."""
    g = pd.read_parquet(_GENCODE_SLIM, columns=["seqname", "feature", "start", "gene_name"])
    g = g[g["feature"] == "gene"].copy()
    g = g[g["seqname"].isin([f"chr{i}" for i in range(1, 23)])]
    g["chrn"] = g["seqname"].str.replace("chr", "", regex=False).astype(int)
    g = g.drop_duplicates("gene_name").set_index("gene_name")[["chrn", "start"]]
    return g.sort_values(["chrn", "start"])


def infercnv_malignancy(log_expr: pd.DataFrame, ref_mask: pd.Series, gene_pos: pd.DataFrame,
                        *, window: int = 100, clip: float = 3.0, exclude_chr: Optional[int] = None) -> pd.Series:
    """Per-unit malignancy score: genome-wide CNV energy from a chromosome-smoothed expression profile.

    inferCNV-style — each gene is centred on a non-malignant REFERENCE mean (`ref_mask`), clipped, and
    smoothed along chromosome position (rolling mean over `window` genes, never crossing chromosome
    boundaries); the score is the mean squared smoothed deviation across the genome. Malignant units carry
    large coordinated chromosomal shifts (CNV) → high score; diploid (immune/stromal) units → low.

    `exclude_chr` drops one chromosome (used to remove a tested gene's OWN-CNV circularity).
    """
    from scipy.ndimage import uniform_filter1d

    pos = gene_pos
    if exclude_chr is not None:
        pos = pos[pos["chrn"] != exclude_chr]
    genes = [g for g in pos.index if g in log_expr.index]
    if len(genes) < window * 2:
        return pd.Series(np.nan, index=log_expr.columns)
    E = log_expr.loc[genes].to_numpy(float)                       # genes(ordered) × units
    ref = E[:, ref_mask.reindex(log_expr.columns).fillna(False).to_numpy()]
    if ref.shape[1] < 10:
        ref = E                                                   # fall back to all units as baseline
    centered = np.clip(E - ref.mean(axis=1, keepdims=True), -clip, clip)
    chrn = pos.loc[genes, "chrn"].to_numpy()
    smoothed = np.empty_like(centered)
    for c in np.unique(chrn):
        idx = np.where(chrn == c)[0]
        smoothed[idx] = uniform_filter1d(centered[idx], size=min(window, len(idx)), axis=0, mode="nearest")
    return pd.Series((smoothed ** 2).mean(axis=0), index=log_expr.columns)


# --------------------------------------------------------------------------- spatial autocorrelation
def _knn_weights(coords: pd.DataFrame, k: int):
    """Row-standardised kNN neighbour index (first neighbour = self, dropped)."""
    tree = cKDTree(coords.to_numpy(float))
    _, nbr = tree.query(coords.to_numpy(float), k=k + 1)
    return nbr[:, 1:]


def bivariate_morans_i(a: pd.Series, b: pd.Series, coords: pd.DataFrame, *, k: int = 6,
                       n_perm: int = 199) -> tuple:
    """Bivariate Moran's I: is high-`a` spatially adjacent to high-`b`? (neighbourhood co-occurrence,
    not just global co-variation). Returns (I_ab, perm_p). Permutation null shuffles `b`'s locations."""
    idx = a.dropna().index.intersection(b.dropna().index).intersection(coords.dropna().index)
    if len(idx) < k + 10:
        return np.nan, np.nan
    za = (a.reindex(idx) - a.reindex(idx).mean()) / (a.reindex(idx).std() or 1)
    zb = (b.reindex(idx) - b.reindex(idx).mean()) / (b.reindex(idx).std() or 1)
    nbr = _knn_weights(coords.reindex(idx), k)
    za_v, zb_v = za.to_numpy(), zb.to_numpy()
    n = len(idx)

    def _I(zb_arr):
        return float((za_v * zb_arr[nbr].mean(axis=1)).sum() / n)

    obs = _I(zb_v)
    rng = np.random.default_rng(0)
    null = np.array([_I(rng.permutation(zb_v)) for _ in range(n_perm)])
    p = float((1 + (null >= obs).sum()) / (1 + n_perm))
    return round(obs, 3), p


def morans_i(values: pd.Series, coords: pd.DataFrame, *, k: int = 6) -> float:
    """Moran's I spatial autocorrelation with row-standardised kNN weights (hand-rolled; no squidpy).

    `coords` is units×(x,y); `values` per unit. Returns I in ~[-1,1] (1 = smooth gradient, 0 = random).
    """
    idx = values.dropna().index.intersection(coords.dropna().index)
    if len(idx) < k + 5:
        return np.nan
    v = values.reindex(idx).to_numpy(float)
    xy = coords.reindex(idx).to_numpy(float)
    z = v - v.mean()
    denom = (z ** 2).sum()
    if denom == 0:
        return np.nan
    tree = cKDTree(xy)
    _, nbr = tree.query(xy, k=k + 1)            # first neighbour is self
    nbr = nbr[:, 1:]
    num = 0.0
    n = len(idx)
    for i in range(n):
        num += (z[i] * z[nbr[i]]).sum() / k     # row-standardised weights (1/k each)
    W = n                                       # sum of row-standardised weights = n
    return float((n / W) * (num / denom))
