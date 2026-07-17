"""Hub-route gene methylation from pipeline per-sample probe betas.

Pipeline ``*_gene_meth.csv`` files cover ``PIPELINE_GENE_PANEL`` (133 genes) and
use promoter-centric aggregation. Most hub-route targets are outside that panel;
HM450 coverage is often gene-body-only (e.g. IRF1) or promoter-sparse (CDKN1A).

This module re-aggregates hub genes from ``data/Methylation/per_sample/*/*_probes.csv``
using GENCODE intervals and the full HM450 probe reference, preferring promoter
beta when ``promoter_n_probes >= min_probes`` and falling back to gene-body beta.
"""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from pipeline.Methylation.aggregation import aggregate_to_genes
from pipeline.Methylation.probe_loader import (
    annotate_probes_with_gene_bodies,
    annotate_probes_with_promoters,
)
from pipeline.config import PATHS, THRESHOLDS
from pipeline.genes.gene_loader import load_gencode_multifeature_subset

DEFAULT_MIN_PROBES = 5
HUB_METH_LONG = C.OUTPUT_ROOT / "matrices" / "hub_gene_methylation.tsv.gz"
HUB_METH_META = C.OUTPUT_ROOT / "matrices" / "hub_gene_methylation_meta.json"
HUB_METH_PROGRESS = C.OUTPUT_ROOT / "matrices" / "hub_gene_methylation_progress.json"
DEFAULT_PROGRESS_INTERVAL_SEC = 30


def _load_probe_reference_tsv() -> pd.DataFrame:
    raw = pd.read_csv(
        PATHS.methylation_probe_reference,
        sep="\t",
        usecols=["CpG_chrm", "CpG_beg", "CpG_end", "probeID"],
    )
    ref = raw.rename(columns={"CpG_chrm": "chrom", "CpG_beg": "start", "CpG_end": "end"})
    ref["chrom"] = ref["chrom"].astype(str).apply(
        lambda c: c if str(c).startswith("chr") else f"chr{c}"
    )
    ref["start"] = pd.to_numeric(ref["start"], errors="coerce")
    ref["end"] = pd.to_numeric(ref["end"], errors="coerce")
    ref["in_CGI"] = False
    ref["CGI_context"] = ""
    return ref.dropna(subset=["start", "end"])


def build_hub_probe_reference(genes: Sequence[str]) -> pd.DataFrame:
    """Annotate HM450 probes for hub genes (promoter + gene body lists)."""
    gene_list = sorted({g for g in genes if g})
    gencode = load_gencode_multifeature_subset(gene_list)
    gencode = gencode.loc[gencode["feature"].eq("gene")].drop_duplicates("gene_name")
    ref = _load_probe_reference_tsv()
    ref = annotate_probes_with_promoters(
        ref,
        gencode,
        gene_panel=gene_list,
        upstream_bp=THRESHOLDS.meth_promoter_upstream_bp,
        downstream_bp=THRESHOLDS.meth_promoter_downstream_bp,
    )
    ref = annotate_probes_with_gene_bodies(
        ref,
        gencode,
        gene_panel=gene_list,
        exclude_promoter=False,
        upstream_bp=THRESHOLDS.meth_promoter_upstream_bp,
        downstream_bp=THRESHOLDS.meth_promoter_downstream_bp,
    )
    return ref


def _hub_probe_ids(probe_ref: pd.DataFrame) -> set[str]:
    prom = probe_ref["promoter_genes"].apply(lambda xs: len(xs) > 0 if isinstance(xs, list) else False)
    body = probe_ref["gene_body_genes"].apply(lambda xs: len(xs) > 0 if isinstance(xs, list) else False)
    return set(probe_ref.loc[prom | body, "probeID"].astype(str))


def _write_progress(
    path: Path,
    *,
    phase: str,
    done: int,
    total: int,
    started_at: float,
    extra: Optional[dict] = None,
) -> None:
    elapsed = time.perf_counter() - started_at
    rate = done / elapsed if elapsed > 0 and done > 0 else 0.0
    remaining = total - done
    eta_sec = remaining / rate if rate > 0 else None
    payload = {
        "phase": phase,
        "done": done,
        "total": total,
        "pct": round(100.0 * done / total, 1) if total else 0.0,
        "elapsed_sec": round(elapsed, 1),
        "eta_sec": round(eta_sec, 1) if eta_sec is not None else None,
        "updated_utc": datetime.now(timezone.utc).isoformat(),
    }
    if extra:
        payload.update(extra)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _maybe_log_progress(
    *,
    phase: str,
    done: int,
    total: int,
    started_at: float,
    last_log_at: float,
    interval_sec: float,
    progress_path: Path,
    extra: Optional[dict] = None,
) -> float:
    now = time.perf_counter()
    if done >= total or (now - last_log_at) >= interval_sec:
        _write_progress(
            progress_path,
            phase=phase,
            done=done,
            total=total,
            started_at=started_at,
            extra=extra,
        )
        elapsed = now - started_at
        rate = done / elapsed if elapsed > 0 and done > 0 else 0.0
        eta = (total - done) / rate if rate > 0 else float("nan")
        eta_s = f"{eta:.0f}s" if rate > 0 else "?"
        print(
            f"[hub_meth] {phase}: {done}/{total} "
            f"({100.0 * done / total:.1f}%) elapsed={elapsed:.0f}s eta={eta_s}",
            flush=True,
        )
        return now
    return last_log_at


def effective_meth_row(
    row: pd.Series,
    *,
    min_probes: int = DEFAULT_MIN_PROBES,
    gene_min_probes: Optional[Dict[str, int]] = None,
) -> tuple[float, str, int]:
    """Return (beta, region, n_probes) using promoter-first, body-fallback."""
    gene = str(row["gene_name"])
    thr = (gene_min_probes or {}).get(gene, min_probes)
    prom_n = int(row.get("promoter_n_probes") or 0)
    body_n = int(row.get("gene_body_n_probes") or 0)
    if prom_n >= thr and pd.notna(row.get("promoter_beta_mean")):
        return float(row["promoter_beta_mean"]), "promoter", prom_n
    if body_n >= thr and pd.notna(row.get("gene_body_beta_mean")):
        return float(row["gene_body_beta_mean"]), "gene_body", body_n
    return np.nan, "NA", 0


def build_hub_gene_methylation(
    genes: Sequence[str],
    *,
    per_sample_dir: Path | None = None,
    out_long: Path = HUB_METH_LONG,
    min_probes: int = DEFAULT_MIN_PROBES,
    gene_min_probes: Optional[Dict[str, int]] = None,
    force: bool = False,
    progress_interval_sec: float = DEFAULT_PROGRESS_INTERVAL_SEC,
    progress_path: Path = HUB_METH_PROGRESS,
) -> pd.DataFrame:
    """Aggregate hub genes across pipeline methylation samples; cache long TSV."""
    gene_list = sorted({g for g in genes if g})
    if not gene_list:
        return pd.DataFrame()

    per_sample_dir = per_sample_dir or (PATHS.methylation_output_dir / "per_sample")
    out_long.parent.mkdir(parents=True, exist_ok=True)

    if out_long.exists() and not force:
        cached = pd.read_csv(out_long, sep="\t")
        if set(gene_list).issubset(set(cached["gene_name"].unique())):
            return cached.loc[cached["gene_name"].isin(gene_list)].copy()

    started_at = time.perf_counter()
    last_log_at = started_at
    _write_progress(
        progress_path,
        phase="annotate_probes",
        done=0,
        total=1,
        started_at=started_at,
        extra={"genes": gene_list},
    )
    print("[hub_meth] annotating HM450 probes for hub genes...", flush=True)
    probe_ref = build_hub_probe_reference(gene_list)
    hub_ids = _hub_probe_ids(probe_ref)
    paths = sorted(per_sample_dir.glob("*/*_probes.csv"))
    rows: List[dict] = []
    n_paths = len(paths)
    last_log_at = _maybe_log_progress(
        phase="aggregate_samples",
        done=0,
        total=n_paths,
        started_at=started_at,
        last_log_at=last_log_at,
        interval_sec=progress_interval_sec,
        progress_path=progress_path,
        extra={"n_hub_probes": len(hub_ids), "genes": gene_list},
    )

    for i, path in enumerate(paths, start=1):
        try:
            probes = pd.read_csv(path, usecols=["probeID", "beta", "participant"])
        except Exception:
            continue
        probes = probes.loc[probes["probeID"].astype(str).isin(hub_ids)].copy()
        if probes.empty:
            continue
        participant = probes["participant"].dropna().iloc[0]
        agg = aggregate_to_genes(
            probes,
            probe_ref,
            gene_panel=gene_list,
            include_gene_body=True,
        )
        for _, r in agg.iterrows():
            beta, region, n_probes = effective_meth_row(
                r, min_probes=min_probes, gene_min_probes=gene_min_probes
            )
            rows.append(
                {
                    "participant": participant,
                    "gene_name": r["gene_name"],
                    "meth_beta_mean": beta,
                    "meth_region": region,
                    "n_probes": n_probes,
                    "promoter_n_probes": int(r.get("promoter_n_probes") or 0),
                    "gene_body_n_probes": int(r.get("gene_body_n_probes") or 0),
                }
            )
        last_log_at = _maybe_log_progress(
            phase="aggregate_samples",
            done=i,
            total=n_paths,
            started_at=started_at,
            last_log_at=last_log_at,
            interval_sec=progress_interval_sec,
            progress_path=progress_path,
            extra={"n_hub_probes": len(hub_ids), "n_rows_so_far": len(rows)},
        )

    long = pd.DataFrame(rows)
    if not long.empty:
        _write_progress(
            progress_path,
            phase="write_cache",
            done=n_paths,
            total=n_paths,
            started_at=started_at,
            extra={"n_rows": len(long)},
        )
        print("[hub_meth] writing cache...", flush=True)
        long.to_csv(out_long, sep="\t", index=False, compression="gzip")
        meta = {
            "built_utc": datetime.now(timezone.utc).isoformat(),
            "n_samples": int(long["participant"].nunique()),
            "genes": gene_list,
            "min_probes": min_probes,
            "n_rows": int(len(long)),
            "elapsed_sec": round(time.perf_counter() - started_at, 1),
        }
        HUB_METH_META.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        _write_progress(
            progress_path,
            phase="done",
            done=n_paths,
            total=n_paths,
            started_at=started_at,
            extra=meta,
        )
        print(
            f"[hub_meth] wrote {out_long} ({meta['n_samples']} samples, "
            f"{meta['elapsed_sec']}s)",
            flush=True,
        )
    else:
        _write_progress(
            progress_path,
            phase="failed_empty",
            done=n_paths,
            total=n_paths,
            started_at=started_at,
        )
    return long


def load_hub_gene_methylation_matrix(
    genes: Sequence[str],
    *,
    per_sample_dir: Path | None = None,
    min_probes: int = DEFAULT_MIN_PROBES,
    gene_min_probes: Optional[Dict[str, int]] = None,
    force_rebuild: bool = False,
) -> pd.DataFrame:
    """Return gene x participant matrix of effective cis methylation beta."""
    long = build_hub_gene_methylation(
        genes,
        per_sample_dir=per_sample_dir,
        min_probes=min_probes,
        gene_min_probes=gene_min_probes,
        force=force_rebuild,
    )
    if long.empty:
        return pd.DataFrame()
    wide = long.pivot_table(
        index="gene_name",
        columns="participant",
        values="meth_beta_mean",
        aggfunc="median",
    )
    return wide.reindex(sorted({g for g in genes if g}))


def main(argv: Optional[Sequence[str]] = None) -> None:
    import argparse

    from mirna_hallmark.robustness_checks import HUB_ROUTES

    p = argparse.ArgumentParser(description="Build hub-route gene methylation cache.")
    p.add_argument("--genes", nargs="*", default=sorted(HUB_ROUTES))
    p.add_argument("--force", action="store_true")
    p.add_argument("--min-probes", type=int, default=DEFAULT_MIN_PROBES)
    p.add_argument(
        "--progress-interval-sec",
        type=float,
        default=DEFAULT_PROGRESS_INTERVAL_SEC,
        help="Log progress to stdout and hub_gene_methylation_progress.json every N seconds.",
    )
    args = p.parse_args(argv)
    build_hub_gene_methylation(
        args.genes,
        min_probes=args.min_probes,
        force=args.force,
        progress_interval_sec=args.progress_interval_sec,
    )


if __name__ == "__main__":
    main()
