"""Roadmap Phase 3 (spatial) — deconvolution upgrade: cell2location (spot-native) vs the NNLS baseline.

The NNLS deconvolution (`spatial_common.deconvolve_matrix`, inherited from Phase 0) treats each Visium
spot as an independent bulk mixture and ignores the count noise model + per-spot detection efficiency —
which is why it over-assigns epithelium (MH-58, epithelial ρ=0.39). cell2location is the spot-native
Bayesian alternative: it learns reference signatures from the Wu scRNA atlas (negative-binomial
regression) then maps each spot's cell abundances with a proper count model.

This module re-deconvolves the primary Visium sections with cell2location and **re-validates the MH-64
program→compartment localization** against the NNLS result — a robustness check (the spatial findings
deliberately avoid the fragile epithelial axis, so agreement here confirms NNLS was adequate; divergence
would flag where the better tool matters).

CPU-only environment (no GPU) → conservative epochs + a balanced reference subsample. Tunable via CLI.

Run: ``.venv/bin/python3 -m mirna_hallmark.spatial_deconv_c2l --samples CID4290``     # benchmark one
     ``.venv/bin/python3 -m mirna_hallmark.spatial_deconv_c2l``                        # default subset
"""

from __future__ import annotations

import argparse
import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.spatial import spatial_common as SP
from mirna_hallmark.brca_deconvolution import SC as WU_SC
from mirna_hallmark.analyses.spatial.spatial_malignancy_retest import _load_visium
from mirna_hallmark.analyses.spatial.spatial_visium_wu import _samples

warnings.filterwarnings("ignore")
OUT_DIR = C.OUTPUT_ROOT / "spatial_deconv_c2l"
DEFAULT_SAMPLES = ["CID4290", "CID44971", "10x_public"]          # representative robustness subset


def load_wu_anndata(per_type: int = 1500, seed: int = 0):
    """Balanced subsample of the Wu scRNA atlas → AnnData(cells×genes raw counts, obs.celltype_major).

    Streams the mtx keeping only the chosen cells' columns (never densifies the full 100k×30k)."""
    import anndata as ad
    from scipy.sparse import csr_matrix

    genes = pd.read_csv(WU_SC / "count_matrix_genes.tsv", header=None)[0].astype(str).values
    barcodes = pd.read_csv(WU_SC / "count_matrix_barcodes.tsv", header=None)[0].astype(str).values
    meta = pd.read_csv(WU_SC / "metadata.csv", index_col=0)
    ct = meta["celltype_major"].reindex(barcodes)
    rng = np.random.default_rng(seed)
    keep_cols = []
    for t in sorted(ct.dropna().unique()):
        idx = np.where(ct.values == t)[0]
        keep_cols.append(rng.choice(idx, size=min(per_type, len(idx)), replace=False))
    keep = np.sort(np.concatenate(keep_cols))                    # 0-based cell indices
    keep_set = set((keep + 1).tolist())                         # mtx columns are 1-based
    colmap = {c: i for i, c in enumerate(keep + 1)}

    rows, cols, vals = [], [], []
    reader = pd.read_csv(WU_SC / "count_matrix_sparse.mtx", sep=r"\s+", skiprows=3, header=None,
                         names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": np.float32},
                         chunksize=20_000_000)
    for chunk in reader:
        m = chunk["c"].isin(keep_set).to_numpy()
        if not m.any():
            continue
        cc = chunk["c"].to_numpy()[m]
        rows.append(np.fromiter((colmap[c] for c in cc), dtype=np.int32))   # cell row
        cols.append(chunk["g"].to_numpy()[m] - 1)                            # gene col
        vals.append(chunk["v"].to_numpy()[m])
    X = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(keep), len(genes)))
    adata = ad.AnnData(X=X)
    adata.var_names = genes
    adata.obs_names = barcodes[keep]
    adata.obs["celltype_major"] = ct.values[keep]
    adata.var_names_make_unique()
    return adata


def reference_signatures(per_type: int, ref_epochs: int) -> pd.DataFrame:
    """cell2location RegressionModel → genes×celltype reference signature (detection-corrected)."""
    import cell2location
    from cell2location.models import RegressionModel

    adata = load_wu_anndata(per_type=per_type)
    # cell2location gene filter (drops very low-count genes)
    from cell2location.utils.filtering import filter_genes
    sel = filter_genes(adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    adata = adata[:, sel].copy()
    RegressionModel.setup_anndata(adata, labels_key="celltype_major")
    mod = RegressionModel(adata)
    mod.train(max_epochs=ref_epochs, accelerator="cpu")
    adata = mod.export_posterior(adata, use_quantiles=True, add_to_varm=["q05"],
                                 sample_kwargs={"batch_size": 2500})
    keys = [k for k in adata.varm.keys() if "q05" in k and "mu_fg" in k] or \
           [k for k in adata.varm.keys() if "means_per_cluster_mu_fg" in k]
    inf = adata.varm[keys[0]].copy()
    inf.columns = [c.split("mu_fg_")[-1] for c in inf.columns]
    print(f"[c2l] reference signature {inf.shape[0]} genes × {inf.shape[1]} types")
    return inf


def map_section(sample: str, inf_aver: pd.DataFrame, sp_epochs: int) -> pd.DataFrame:
    """cell2location spatial mapping → spots×celltype abundance fractions."""
    import anndata as ad
    import cell2location

    log, lin, pos = _load_visium(sample)
    counts = np.round(np.power(2.0, log) - 1.0).clip(lower=0)    # back to ~counts for the NB model
    shared = counts.index.intersection(inf_aver.index)
    adata = ad.AnnData(X=counts.loc[shared].T.values)            # spots × genes
    adata.var_names = shared
    adata.obs_names = counts.columns
    cell2location.models.Cell2location.setup_anndata(adata)
    mod = cell2location.models.Cell2location(adata, cell_state_df=inf_aver.loc[shared],
                                             N_cells_per_location=8, detection_alpha=20)
    mod.train(max_epochs=sp_epochs, batch_size=None, train_size=1, accelerator="cpu")
    adata = mod.export_posterior(adata, use_quantiles=True, add_to_obsm=["q05"],
                                 sample_kwargs={"batch_size": adata.n_obs})
    ab = adata.obsm[[k for k in adata.obsm if "q05" in k and "abundance" in k][0]]
    ab.columns = inf_aver.columns
    frac = ab.div(ab.sum(axis=1).replace(0, np.nan), axis=0)
    frac.index = counts.columns
    return frac


def run(*, out_dir: Path = OUT_DIR, samples: List[str] | None = None,
        per_type: int = 1500, ref_epochs: int = 250, sp_epochs: int = 3000) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    samples = samples or DEFAULT_SAMPLES
    ref_cache = out_dir / "reference_signatures.parquet"
    if ref_cache.exists():
        inf_aver = pd.read_parquet(ref_cache)
        print(f"[c2l] reusing cached reference signature {inf_aver.shape} ({ref_cache})")
    else:
        print(f"[c2l] reference (per_type={per_type}, ref_epochs={ref_epochs}) — slow on CPU, cached after ...")
        inf_aver = reference_signatures(per_type, ref_epochs)
        inf_aver.to_parquet(ref_cache)
    pg = SP.program_gene_sets()

    rows = []
    for s in samples:
        t0 = datetime.now(timezone.utc)
        cfr = map_section(s, inf_aver, sp_epochs)
        comp = SP.compartment_fractions(cfr)
        log, lin, pos = _load_visium(s)
        prog = SP.score_programs(log, pg)
        for p in prog.columns:
            for c in comp.columns:
                rows.append({"sample": s, "program": p, "compartment": c,
                             "c2l_rho": round(float(spearmanr(prog[p].reindex(comp.index), comp[c])[0]), 3)})
        dt = (datetime.now(timezone.utc) - t0).total_seconds()
        print(f"   {s}: mapped in {dt:.0f}s | mean fracs {comp.mean().round(2).to_dict()}")
    long = pd.DataFrame(rows)
    long.to_csv(out_dir / "c2l_program_compartment_long.tsv", sep="\t", index=False)

    # dominant compartment per program (c2l) vs the NNLS baseline (spatial_visium_wu manifest)
    c2l_dom = (long.groupby(["program", "compartment"])["c2l_rho"].median().reset_index()
               .sort_values("c2l_rho", ascending=False).groupby("program").first()["compartment"].to_dict())
    nnls_manifest = C.OUTPUT_ROOT / "spatial_visium_wu" / "method_manifest.json"
    nnls_dom = {}
    if nnls_manifest.exists():
        loc = json.loads(nnls_manifest.read_text())["program_compartment_localization_median"]
        nnls_dom = {p: max(comps, key=comps.get) for p, comps in loc.items()}
    agree = {p: (c2l_dom.get(p) == nnls_dom.get(p)) for p in c2l_dom}

    summary = {
        "module": "mirna_hallmark.spatial_deconv_c2l",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "cell2location (spot-native Bayesian) vs NNLS baseline",
        "settings": {"samples": samples, "per_type": per_type, "ref_epochs": ref_epochs, "sp_epochs": sp_epochs},
        "c2l_dominant_compartment": c2l_dom, "nnls_dominant_compartment": nnls_dom,
        "agreement_with_nnls": agree, "n_programs_agree": int(sum(agree.values())),
        "caveats": ["CPU-only → reduced epochs + reference subsample (approximate posterior)",
                    "robustness check on the localization claim; spatial findings don't rest on epithelial axis",
                    "still no spatial miRNA (binding constraint)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[c2l] dominant compartment (cell2location): {c2l_dom}")
    print(f"[c2l] vs NNLS: {nnls_dom}")
    print(f"[c2l] programs agreeing: {sum(agree.values())}/{len(agree)} — {agree}")
    print(f"[c2l] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--per-type", type=int, default=1500)
    ap.add_argument("--ref-epochs", type=int, default=250)
    ap.add_argument("--sp-epochs", type=int, default=3000)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, samples=args.samples, per_type=args.per_type,
        ref_epochs=args.ref_epochs, sp_epochs=args.sp_epochs)


if __name__ == "__main__":
    main()
