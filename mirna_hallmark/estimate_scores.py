"""Actual ESTIMATE (Yoshihara et al. 2013) immune/stromal scores per tissue state.

Runs the canonical R `estimate` package (official 141-gene stromal + 141-gene
immune signatures, package ssGSEA) on the FULL transcriptome of each state
(TCGA tumor / TCGA NAT / GTEx healthy breast) and caches per-sample
StromalScore / ImmuneScore / ESTIMATEScore. These replace the hand-curated
metagene proxies as composition covariates in `cross_state_coupling`.

Why ESTIMATE rather than precomputed CIBERSORT/Thorsson: those exist for TCGA
tumor only; ESTIMATE is computed identically from each state's own expression
(no platform-specific reference matrix), so the immune/stromal adjustment is
comparable across tumor / NAT / GTEx. platform="illumina" (RNA-seq) yields the
three scores; TumorPurity (affymetrix-calibrated) is intentionally not produced
because it has no NAT/GTEx analog.

ssGSEA ranks genes within each sample, so the monotone log2 transform we feed is
equivalent to linear TPM/RPM for scoring.

Run:
  .venv/bin/python3 -m mirna_hallmark.estimate_scores --build      # all 3 states
"""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.family_normal_reference import (
    GTEX_BREAST_TPM,
    _gtex_available,
    _gtex_donor,
    _participant,
    _split_types,
)

OUT_DIR = C.TISSUE_REFERENCE_DIR / "estimate"
R_RUNNER = Path(__file__).resolve().parent / "run_estimate.R"
PLATFORM = "illumina"


# --------------------------------------------------------------------------- #
# Full-transcriptome loaders (ESTIMATE needs the whole matrix, not focus genes)
# --------------------------------------------------------------------------- #
def _full_rna_tcga() -> pd.DataFrame:
    """gene_symbol x FULL-barcode log2(TPM+1), all genes."""
    df = pd.read_csv(C.RNA_EXPRESSION, sep="\t", index_col=0)
    if "Ensembl_ID" in df.columns:
        df = df.drop(columns=["Ensembl_ID"])
    df = df.apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    return np.log2(df.clip(lower=0) + 1)


def _full_rna_gtex() -> pd.DataFrame:
    """gene_symbol x GTEx-breast-donor log2(TPM+1), all genes."""
    rows = []
    for chunk in pd.read_csv(GTEX_BREAST_TPM, sep="\t", skiprows=2, chunksize=8000, low_memory=False):
        rows.append(chunk)
    df = pd.concat(rows)
    df = df.drop(columns=["Name"]).set_index("Description")
    df = df.apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    df.columns = [_gtex_donor(c) for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return np.log2(df.clip(lower=0) + 1)


def _state_full_rna() -> Dict[str, pd.DataFrame]:
    tcga = _full_rna_tcga()
    split = _split_types(tcga)

    def _to_part(d: pd.DataFrame) -> pd.DataFrame:
        d = d.rename(columns={c: _participant(c) for c in d.columns})
        return d.loc[:, ~d.columns.duplicated()]

    states = {"tumor": _to_part(split["tumor"]), "nat": _to_part(split["nat"])}
    if _gtex_available():
        states["gtex"] = _full_rna_gtex()
    return states


# --------------------------------------------------------------------------- #
# GCT IO + R invocation
# --------------------------------------------------------------------------- #
def _write_expr_tsv(mat: pd.DataFrame, path: Path) -> None:
    """Plain tab-delimited matrix for `filterCommonGenes` (read.table header=TRUE,
    row.names=1): first column = gene symbol, header = GeneSymbol + sample names."""
    mat = mat.dropna(how="all")
    with path.open("w", encoding="utf-8") as fh:
        fh.write("GeneSymbol\t" + "\t".join(map(str, mat.columns)) + "\n")
        for gene, vals in mat.iterrows():
            fh.write(str(gene) + "\t" + "\t".join(f"{v:.4f}" for v in vals.values) + "\n")


def _parse_scores(scores_gct: Path, samples: pd.Index) -> pd.DataFrame:
    """ESTIMATE output GCT -> sample x {StromalScore, ImmuneScore, ESTIMATEScore}.

    R mangles sample names (- -> .), so reassign columns by position (order is
    preserved by estimateScore).
    """
    raw = pd.read_csv(scores_gct, sep="\t", skiprows=2, index_col=0)
    raw = raw.drop(columns=[c for c in ("Description",) if c in raw.columns])
    if raw.shape[1] != len(samples):
        raise SystemExit(f"[estimate] score columns {raw.shape[1]} != samples {len(samples)}")
    raw.columns = list(samples)
    return raw.T


def compute_estimate(state: str, rna: pd.DataFrame, *, out_dir: Path = OUT_DIR,
                     platform: str = PLATFORM, force: bool = False) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / f"{state}_estimate_scores.tsv"
    if cache.is_file() and not force:
        return pd.read_csv(cache, sep="\t", index_col=0)
    with tempfile.TemporaryDirectory() as td:
        gct = Path(td) / f"{state}_expr.txt"
        scores = Path(td) / f"{state}_scores.gct"
        _write_expr_tsv(rna, gct)
        res = subprocess.run(
            ["Rscript", str(R_RUNNER), str(gct), str(scores), platform],
            capture_output=True, text=True,
        )
        if "ESTIMATE_DONE" not in res.stdout or not scores.is_file():
            raise SystemExit(f"[estimate] R failed for {state}:\n{res.stdout}\n{res.stderr}")
        df = _parse_scores(scores, rna.columns)
    df.index.name = "sample"
    df.to_csv(cache, sep="\t")
    print(f"[estimate] {state}: n={len(df)}  "
          f"ImmuneScore[{df['ImmuneScore'].min():.0f},{df['ImmuneScore'].max():.0f}]  "
          f"StromalScore[{df['StromalScore'].min():.0f},{df['StromalScore'].max():.0f}]")
    return df


def load_estimate(state: str, *, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    cache = out_dir / f"{state}_estimate_scores.tsv"
    if not cache.is_file():
        raise SystemExit(f"[estimate] missing cache {cache}; run `--build` first")
    return pd.read_csv(cache, sep="\t", index_col=0)


def build_all(*, out_dir: Path = OUT_DIR, force: bool = False) -> Dict[str, pd.DataFrame]:
    states = _state_full_rna()
    return {s: compute_estimate(s, rna, out_dir=out_dir, force=force) for s, rna in states.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build", action="store_true", help="compute + cache ESTIMATE for all states")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    if args.build:
        build_all(out_dir=args.out_dir, force=args.force)


if __name__ == "__main__":
    main()
