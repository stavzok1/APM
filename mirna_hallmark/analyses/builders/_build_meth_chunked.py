"""Resumable, chunked Hallmark-target methylation builder.

The one-shot build of per-sample methylation over hundreds of genes exceeds the ~20-min process
wall in this environment (reading 796 full HM450 probe CSVs dominates). This script processes a
slice of samples per invocation and writes a part file, so a death never loses progress; run it
until all parts exist, then it concatenates into the final long table.

Usage:
    python -m mirna_hallmark._build_meth_chunked --batch 200      # process next 200 unprocessed
    python -m mirna_hallmark._build_meth_chunked --finalize       # concat parts -> final tsv.gz
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.config import PATHS
from mirna_hallmark.analyses.misc.hub_gene_methylation import (
    DEFAULT_MIN_PROBES,
    aggregate_to_genes,
    build_hub_probe_reference,
    effective_meth_row,
    _hub_probe_ids,
)

OUT = C.OUTPUT_ROOT / "matrices" / "hallmark_gene_methylation.tsv.gz"
PARTS = C.OUTPUT_ROOT / "matrices" / "_hallmark_meth_parts"
CURSOR = PARTS / "_cursor.txt"


def _genes() -> list[str]:
    e = pd.read_csv(C.TISSUE_REFERENCE_DIR / "mirna_state_class" / "mirna_gene_edge_class.tsv", sep="\t")
    dec = e[e["joint_edge_class"].isin(["lost", "nat_decoupled"])]
    return sorted(set(dec["gene"]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--finalize", action="store_true")
    args = ap.parse_args()

    PARTS.mkdir(parents=True, exist_ok=True)
    genes = _genes()

    if args.finalize:
        frames = [pd.read_csv(p, sep="\t") for p in sorted(PARTS.glob("part_*.tsv"))]
        if not frames:
            print("no parts to finalize")
            return
        long = pd.concat(frames, ignore_index=True)
        long = long.drop_duplicates(subset=["participant", "gene_name"], keep="last")
        long.to_csv(OUT, sep="\t", index=False, compression="gzip")
        print(f"FINALIZED {OUT} rows={len(long)} participants={long['participant'].nunique()} "
              f"genes={long['gene_name'].nunique()}")
        return

    print(f"[chunk] annotating probes for {len(genes)} decoupled-target genes...", flush=True)
    probe_ref = build_hub_probe_reference(genes)
    hub_ids = _hub_probe_ids(probe_ref)
    paths = sorted((PATHS.methylation_output_dir / "per_sample").glob("*/*_probes.csv"))

    # deterministic index cursor (the per-sample DIR name carries a vial suffix that does not match
    # the participant id, so we advance by position, not by participant set).
    offset = int(CURSOR.read_text().strip()) if CURSOR.exists() else 0
    todo = paths[offset: offset + args.batch]
    print(f"[chunk] offset={offset}/{len(paths)}; this batch={len(todo)}", flush=True)
    if not todo:
        print("ALL_DONE", flush=True)
        return

    rows = []
    t0 = time.perf_counter()
    for i, path in enumerate(todo, 1):
        try:
            probes = pd.read_csv(path, usecols=["probeID", "beta", "participant"])
        except Exception:
            continue
        probes = probes.loc[probes["probeID"].astype(str).isin(hub_ids)].copy()
        if probes.empty:
            continue
        participant = probes["participant"].dropna().iloc[0]
        agg = aggregate_to_genes(probes, probe_ref, gene_panel=genes, include_gene_body=True)
        for _, r in agg.iterrows():
            beta, region, n_probes = effective_meth_row(r, min_probes=DEFAULT_MIN_PROBES)
            rows.append({"participant": participant, "gene_name": r["gene_name"],
                         "meth_beta_mean": beta, "meth_region": region, "n_probes": n_probes,
                         "promoter_n_probes": int(r.get("promoter_n_probes") or 0),
                         "gene_body_n_probes": int(r.get("gene_body_n_probes") or 0)})
        if i % 50 == 0:
            rate = i / (time.perf_counter() - t0)
            print(f"[chunk] {i}/{len(todo)} rate={rate:.2f}/s eta={ (len(todo)-i)/rate:.0f}s",
                  flush=True)
    part = PARTS / f"part_{offset:04d}.tsv"
    pd.DataFrame(rows).to_csv(part, sep="\t", index=False)
    CURSOR.write_text(str(offset + len(todo)))
    print(f"BATCH_DONE wrote {part} rows={len(rows)} next_offset={offset + len(todo)}", flush=True)


if __name__ == "__main__":
    main()
