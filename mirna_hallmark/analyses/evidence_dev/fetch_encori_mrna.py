"""Fetch ENCORI miRNATarget rows for protein-coding (mRNA) Hallmark genes.

Uses ``pipeline.lncRNA_interactions.encori`` with ``geneType=mRNA`` (the only
value that returns rows for coding genes). Per-gene TSVs are cached under
``data/external_cache/encori/miRNATarget/``; this module consolidates them into
``mirna_hallmark/output/encori/encori_mirna_mrna_targets.parquet``.

Run:
    .venv/bin/python3 -m mirna_hallmark.fetch_encori_mrna
    .venv/bin/python3 -m mirna_hallmark.fetch_encori_mrna --limit 50   # smoke
    .venv/bin/python3 -m mirna_hallmark.fetch_encori_mrna --resume     # skip cached genes
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from pipeline.config import PATHS
from pipeline.lncRNA_interactions.encori import (
    EncoriQuery,
    fetch_encori_miRNA_targets_with_symbol_fallback,
)

OUT_DIR = C.OUTPUT_ROOT / "encori"
DEFAULT_PARQUET = OUT_DIR / "encori_mirna_mrna_targets.parquet"
DEFAULT_DIAG = OUT_DIR / "encori_fetch_diagnostics.tsv"


def _gene_list(*, source: str) -> list[str]:
    if source == "hallmark":
        return HallmarkSets.load().universe
    if source == "mirtar":
        df = pd.read_csv(C.MIRTAR_HALLMARK_SUMMARY, usecols=["gene", "miRNA"], low_memory=False)
        df = df.loc[df["miRNA"].astype(str).str.startswith("hsa-", na=False)]
        return sorted(df["gene"].dropna().astype(str).unique())
    raise ValueError(f"Unknown gene source: {source!r}")


def _read_valid_cache_tsv(path: Path) -> pd.DataFrame:
    """Parse a cached ENCORI TSV; return empty if error stub or invalid."""
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False, comment="#", on_bad_lines="skip")
    except Exception:
        return pd.DataFrame()
    if df.empty or "miRNAname" not in df.columns:
        return pd.DataFrame()
    return df


def _cache_path_for_gene(cache_dir: Path, gene: str, query: EncoriQuery) -> Path:
    from pipeline.lncRNA_interactions.encori import _cache_key, _encori_url

    params = {
        "assembly": query.assembly,
        "geneType": query.gene_type,
        "miRNA": "all",
        "clipExpNum": int(query.clip_exp_num),
        "degraExpNum": int(query.degra_exp_num),
        "pancancerNum": int(query.pancancer_num),
        "programNum": int(query.program_num),
        "program": query.program,
        "target": gene,
        "cellType": query.cell_type,
    }
    url = _encori_url("miRNATarget", params)
    return cache_dir / "encori" / "miRNATarget" / f"{gene}.{_cache_key(url)}.tsv"


def _gene_has_valid_cache(cache_dir: Path, gene: str, query: EncoriQuery) -> bool:
    path = _cache_path_for_gene(cache_dir, gene, query)
    return path.is_file() and not _read_valid_cache_tsv(path).empty


def consolidate_cache(
    genes: Sequence[str],
    *,
    cache_dir: Path,
    query: EncoriQuery,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rebuild combined table from on-disk per-gene TSV cache (no network)."""
    frames: list[pd.DataFrame] = []
    diag_rows: list[dict] = []
    for gene in genes:
        path = _cache_path_for_gene(cache_dir, gene, query)
        df = _read_valid_cache_tsv(path) if path.is_file() else pd.DataFrame()
        n = len(df)
        diag_rows.append({
            "query_gene": gene,
            "cache_path": str(path) if path.is_file() else "",
            "n_rows_total": n,
            "from_cache": True,
        })
        if n:
            df = df.copy()
            df["__query_gene"] = gene
            frames.append(df)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return combined, pd.DataFrame(diag_rows)


def fetch_encori_mrna_targets(
    genes: Sequence[str],
    *,
    cache_dir: Path,
    query: EncoriQuery,
    resume: bool = False,
    progress_every: int = 25,
    max_retries: int = 3,
    retry_sleep_s: float = 5.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch miRNA-target rows for *genes*; return (combined, diagnostics)."""
    gene_list = [str(g).strip() for g in genes if str(g).strip()]

    frames: list[pd.DataFrame] = []
    diag_rows: list[dict] = []
    failed: list[str] = []

    for i, gene in enumerate(gene_list, start=1):
        if resume and _gene_has_valid_cache(cache_dir, gene, query):
            df, diag = consolidate_cache([gene], cache_dir=cache_dir, query=query)
        else:
            df = pd.DataFrame()
            diag = pd.DataFrame()
            last_err: Exception | None = None
            for attempt in range(1, max_retries + 1):
                try:
                    df, diag = fetch_encori_miRNA_targets_with_symbol_fallback(
                        targets=[gene],
                        miRNA="all",
                        query=query,
                        cache_dir=cache_dir,
                        clip_exp_num=query.clip_exp_num,
                    )
                    last_err = None
                    break
                except (urllib.error.URLError, TimeoutError, OSError) as exc:
                    last_err = exc
                    if attempt < max_retries:
                        wait = retry_sleep_s * attempt
                        print(
                            f"[fetch_encori_mrna] {gene}: attempt {attempt}/{max_retries} "
                            f"failed ({exc}); retry in {wait:.0f}s"
                        )
                        time.sleep(wait)
            if last_err is not None:
                print(f"[fetch_encori_mrna] {gene}: FAILED after {max_retries} attempts — skipping")
                failed.append(gene)
                diag_rows.append({
                    "query_gene": gene,
                    "n_rows_total": 0,
                    "error": str(last_err),
                })
                continue

        if not df.empty:
            df = df.copy()
            df["__query_gene"] = gene
            frames.append(df)
        if not diag.empty:
            row = diag.iloc[0].to_dict()
            row["query_gene"] = gene
            diag_rows.append(row)
        elif not df.empty:
            diag_rows.append({"query_gene": gene, "n_rows_total": len(df)})

        if progress_every and (i % progress_every == 0 or i == len(gene_list)):
            n_rows = sum(len(f) for f in frames)
            n_hit = sum(1 for r in diag_rows if int(r.get("n_rows_total", 0) or 0) > 0)
            print(
                f"[fetch_encori_mrna] {i}/{len(gene_list)} genes "
                f"({n_hit} with hits, {n_rows:,} rows so far"
                + (f", {len(failed)} failed" if failed else "")
                + ")"
            )

    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    diagnostics = pd.DataFrame(diag_rows)
    if failed:
        fail_path = OUT_DIR / "fetch_failed_genes.txt"
        fail_path.write_text("\n".join(failed) + "\n", encoding="utf-8")
        print(f"[fetch_encori_mrna] wrote {len(failed)} failed genes -> {fail_path}")
    return combined, diagnostics


def run(
    *,
    genes: Sequence[str] | None = None,
    gene_source: str = "mirtar",
    cache_dir: Path | None = None,
    out_parquet: Path = DEFAULT_PARQUET,
    out_diag: Path = DEFAULT_DIAG,
    clip_exp_num: int = 1,
    program_num: int = 1,
    limit: int | None = None,
    resume: bool = False,
    consolidate_only: bool = False,
) -> None:
    cache = Path(cache_dir or PATHS.working_dir / "external_cache")
    query = EncoriQuery(gene_type="mRNA", clip_exp_num=int(clip_exp_num), program_num=int(program_num))

    gene_list = list(genes) if genes is not None else _gene_list(source=gene_source)
    if limit is not None:
        gene_list = gene_list[: int(limit)]

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if consolidate_only:
        print(f"[fetch_encori_mrna] consolidating cache for {len(gene_list):,} genes (no network)")
        combined, diagnostics = consolidate_cache(gene_list, cache_dir=cache, query=query)
    else:
        print(
            f"[fetch_encori_mrna] fetching {len(gene_list):,} genes "
            f"(source={gene_source}, geneType=mRNA, clipExpNum={clip_exp_num})"
        )
        combined, diagnostics = fetch_encori_mrna_targets(
            gene_list,
            cache_dir=cache,
            query=query,
            resume=resume,
            progress_every=25,
        )

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(out_parquet, index=False)
    diagnostics.to_csv(out_diag, sep="\t", index=False)

    n_hit = int((diagnostics.get("n_rows_total", pd.Series(dtype=float)) > 0).sum()) if not diagnostics.empty else 0
    manifest = {
        "module": "mirna_hallmark.fetch_encori_mrna",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "gene_source": gene_source,
        "n_genes_queried": len(gene_list),
        "n_genes_with_hits": n_hit,
        "n_rows": int(len(combined)),
        "cache_dir": str(cache),
        "out_parquet": str(out_parquet),
        "query": {
            "gene_type": "mRNA",
            "clip_exp_num": clip_exp_num,
            "program_num": program_num,
            "program": query.program,
        },
    }
    (OUT_DIR / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(
        f"[fetch_encori_mrna] done: {len(combined):,} rows, "
        f"{n_hit}/{len(gene_list)} genes with hits -> {out_parquet}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gene-source", choices=["mirtar", "hallmark"], default="mirtar")
    ap.add_argument("--cache-dir", type=Path, default=None)
    ap.add_argument("--out-parquet", type=Path, default=DEFAULT_PARQUET)
    ap.add_argument("--out-diag", type=Path, default=DEFAULT_DIAG)
    ap.add_argument("--clip-exp-num", type=int, default=1)
    ap.add_argument("--program-num", type=int, default=1)
    ap.add_argument("--limit", type=int, default=None, help="Fetch first N genes only (smoke test)")
    ap.add_argument("--resume", action="store_true", help="Skip genes with valid on-disk cache")
    ap.add_argument(
        "--consolidate-only",
        action="store_true",
        help="Rebuild parquet from cache only (no network)",
    )
    args = ap.parse_args()
    run(
        gene_source=args.gene_source,
        cache_dir=args.cache_dir,
        out_parquet=args.out_parquet,
        out_diag=args.out_diag,
        clip_exp_num=args.clip_exp_num,
        program_num=args.program_num,
        limit=args.limit,
        resume=args.resume,
        consolidate_only=args.consolidate_only,
    )


if __name__ == "__main__":
    main()
