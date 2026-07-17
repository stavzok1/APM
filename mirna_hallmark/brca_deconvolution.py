"""Roadmap Phase 0: deconvolve TCGA-BRCA bulk into cell-type fractions (primary-tumour composition).

Builds a cell-type signature from the Wu et al. 2021 primary-breast scRNA atlas (GSE176078; 100,064
cells, `celltype_major`: Cancer/Normal Epithelial, CAFs, T-cells, Myeloid, B-cells, Plasmablasts,
Endothelial, PVL) and NNLS-deconvolves each TCGA bulk sample → per-sample fractions. This is the
*refined* composition reference (actual cell-type fractions) that upgrades the marker-metagene proxy
used in Phase 1 (`core_coupling_composition_retest`) and feeds Phases 1-refined / 2 / 3.

Memory-safe: the 2.4 GB sparse mtx is **streamed** (chunked) into a genes×celltype pseudobulk
accumulator (never densified). NNLS on the top cross-cell-type-variable genes.

Run: ``.venv/bin/python3 -m mirna_hallmark.brca_deconvolution``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cross_state.cross_state_coupling import EPI_MARKERS, IMMUNE_MARKERS, STROMA_MARKERS, _metagene

SC = C.REPO_ROOT / "data" / "external" / "brca_scrna_gse176078" / "Wu_etal_2021_BRCA_scRNASeq"
OUT_DIR = C.OUTPUT_ROOT / "brca_deconvolution"
N_SIG_GENES = 1500


def build_signature() -> pd.DataFrame:
    """Stream the sc mtx -> genes x celltype_major pseudobulk CPM signature (linear)."""
    genes = pd.read_csv(SC / "count_matrix_genes.tsv", header=None)[0].astype(str).values
    barcodes = pd.read_csv(SC / "count_matrix_barcodes.tsv", header=None)[0].astype(str).values
    meta = pd.read_csv(SC / "metadata.csv", index_col=0)
    ct = meta["celltype_major"].reindex(barcodes)
    types = sorted(ct.dropna().unique())
    tcode = {t: i for i, t in enumerate(types)}
    cell_type_idx = ct.map(tcode).to_numpy()                       # per-column (cell) type code, NaN if unmapped
    acc = np.zeros((len(genes), len(types)), dtype=np.float64)     # genes x types count sums

    mtx = SC / "count_matrix_sparse.mtx"
    reader = pd.read_csv(mtx, sep=r"\s+", skiprows=3, header=None,
                         names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": np.float32},
                         chunksize=20_000_000)
    for chunk in reader:
        cidx = cell_type_idx[chunk["c"].to_numpy() - 1]            # 1-indexed -> 0-indexed cell
        ok = ~np.isnan(cidx)
        np.add.at(acc, (chunk["g"].to_numpy()[ok] - 1, cidx[ok].astype(int)), chunk["v"].to_numpy()[ok])
    sig = pd.DataFrame(acc, index=genes, columns=types)
    sig = sig.div(sig.sum(axis=0), axis=1) * 1e6                   # CPM per cell type (linear)
    return sig


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[deconv] building cell-type signature from GSE176078 (streaming) ...")
    sig = build_signature()
    print(f"[deconv] signature {sig.shape[0]} genes x {sig.shape[1]} types: {list(sig.columns)}")

    rna = D.load_rna()                                             # NB: load_rna double-logs (file is already log2(TPM+1))
    rna = rna[~rna.index.duplicated(keep="first")]
    rna_lin = np.power(2.0, np.power(2.0, rna) - 1.0) - 1.0        # undo BOTH logs -> linear TPM (deconvolution needs linear)

    shared = sig.index.intersection(rna_lin.index)
    sig_s = sig.loc[shared]
    # BALANCED markers: top fold-change genes PER cell type (else sharp immune/stromal markers
    # dominate and epithelial is under-fit — the failure the epi-metagene control caught).
    per_type = max(60, N_SIG_GENES // sig_s.shape[1])
    markers = set()
    for t in sig_s.columns:
        others = sig_s.drop(columns=t).mean(axis=1)
        fc = (sig_s[t] + 1.0) / (others + 1.0)
        # require the marker to be reasonably expressed in its type
        cand = fc[sig_s[t] > sig_s[t].quantile(0.5)].sort_values(ascending=False)
        markers |= set(cand.head(per_type).index)
    markers = sorted(markers)
    S = sig_s.loc[markers].to_numpy()
    B = rna_lin.loc[markers].to_numpy()                            # markers x samples
    print(f"[deconv] NNLS on {len(markers)} marker genes x {B.shape[1]} TCGA samples ...")

    fracs = np.zeros((B.shape[1], S.shape[1]))
    for j in range(B.shape[1]):
        x, _ = nnls(S, B[:, j])
        fracs[j] = x / x.sum() if x.sum() > 0 else x
    frac = pd.DataFrame(fracs, index=rna_lin.columns, columns=sig.columns)
    frac.to_csv(out_dir / "tcga_celltype_fractions.tsv", sep="\t")

    # validate vs the marker-metagene proxies used in Phase 1
    epi_mg = _metagene(rna, EPI_MARKERS); imm_mg = _metagene(rna, IMMUNE_MARKERS); str_mg = _metagene(rna, STROMA_MARKERS)
    epi_f = frac.get("Cancer Epithelial", pd.Series(0, index=frac.index)) + frac.get("Normal Epithelial", 0)
    val = {
        "epithelial_frac_vs_epi_metagene": round(float(spearmanr(epi_f, epi_mg.reindex(frac.index))[0]), 3),
        "immune_frac_vs_immune_metagene": round(float(spearmanr(
            frac.get("T-cells", 0) + frac.get("Myeloid", 0) + frac.get("B-cells", 0) + frac.get("Plasmablasts", 0),
            imm_mg.reindex(frac.index))[0]), 3),
        "CAF_frac_vs_stroma_metagene": round(float(spearmanr(frac.get("CAFs", 0), str_mg.reindex(frac.index))[0]), 3),
    }
    summary = {
        "module": "mirna_hallmark.brca_deconvolution",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "reference": "GSE176078 Wu 2021 primary-breast scRNA (celltype_major)",
        "n_tcga_samples": int(frac.shape[0]), "cell_types": list(sig.columns), "n_marker_genes": int(len(markers)),
        "mean_fraction_by_type": frac.mean().round(3).to_dict(),
        "validation_vs_metagene_spearman": val,
        "caveats": ["NNLS bulk deconvolution; CPM(sc) vs TPM(bulk) scale; marker-gene subset",
                    "celltype_major granularity; reference is primary tumours (Wu); feeds Phases 1-refined/2/3"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[deconv] mean fractions: {summary['mean_fraction_by_type']}")
    print(f"[deconv] validation vs metagenes (Spearman): {val}")
    print(f"[deconv] wrote {out_dir}")
    return frac


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
