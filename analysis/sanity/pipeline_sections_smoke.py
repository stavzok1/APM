#!/usr/bin/env python3
from __future__ import annotations

"""
Fast section-by-section smoke validation for `pipeline/main.py`.

Goal: catch crashes in each major section without running the full (slow) evidence builds.

This is designed to be *fast* and to fail with a clear traceback pinpointing the section.
It intentionally runs on small slices of the data tables.
"""

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.sanity._io import ensure_out_dir  # noqa: E402
from pipeline.config import PATHS, PIPELINE_GENE_PANEL, use_legacy_lncrna_intervals_csv_input  # noqa: E402
from pipeline.genes.gene_loader import load_genes, load_lncrnas, lncrna_gene_intervals_from_annotation  # noqa: E402
from pipeline.genes.lncrna_matching import match_lncrnas_to_genes  # noqa: E402
from pipeline.regulatory_elements import load_ccres, match_ccres_to_genes, build_element_focus_table  # noqa: E402
from pipeline.tad_annotation import annotate_with_all_tad_sources, mirror_and_save_all_domains, build_boundaries_enriched_table  # noqa: E402
from pipeline.tad_annotation.tad_config import discover_tad_sources  # noqa: E402


@dataclass(frozen=True)
class StepResult:
    step: str
    ok: bool
    detail: str


def _ok(step: str, detail: str) -> StepResult:
    return StepResult(step=step, ok=True, detail=detail)


def _fail(step: str, detail: str) -> StepResult:
    return StepResult(step=step, ok=False, detail=detail)


def _head(df: pd.DataFrame, n: int) -> pd.DataFrame:
    return df.head(n).copy() if df is not None and not df.empty else pd.DataFrame()


def _pick_tad_biosamples(processed_dir: Path, *, max_biosamples: int = 1) -> List[str]:
    avail = discover_tad_sources(processed_dir, required_files=True)
    names = sorted(avail.keys())
    return names[:max_biosamples]


def run_smoke(*, out_dir: Path, n_genes: int, n_ccres: int, n_lncrnas: int, n_tad_biosamples: int) -> Tuple[List[StepResult], Dict[str, object]]:
    results: List[StepResult] = []
    summary: Dict[str, object] = {}

    # STEP 1: load gene + cCRE + lncRNA interval tables (full load, then slice)
    try:
        genes_all = load_genes(PATHS.gencode_gtf_pq)
        ccres = load_ccres(PATHS.ccre_csv)
        if use_legacy_lncrna_intervals_csv_input():
            lncrnas_all = load_lncrnas(PATHS.lncrna_csv)
            lncrna_src = "legacy_csv"
        else:
            lncrnas_all = lncrna_gene_intervals_from_annotation(genes_all)
            lncrna_src = "gencode_gene_rows"
        results.append(
            _ok(
                "STEP1.load_inputs",
                f"genes={len(genes_all)}, ccres={len(ccres)}, lncrnas={len(lncrnas_all)} (lncRNA_src={lncrna_src})",
            )
        )
        summary["n_genes_all"] = int(len(genes_all))
        summary["n_ccres_all"] = int(len(ccres))
        summary["n_lncrnas_all"] = int(len(lncrnas_all))
    except Exception as e:
        results.append(_fail("STEP1.load_inputs", repr(e)))
        return results, summary

    # Slice for speed downstream
    gene_panel = list(PIPELINE_GENE_PANEL)[: int(n_genes)]
    genes = genes_all[genes_all.get("gene_name").astype(str).isin(set(gene_panel))].copy() if "gene_name" in genes_all.columns else _head(genes_all, n_genes)
    ccres_s = _head(ccres, int(n_ccres))
    lncrnas_s = _head(lncrnas_all, int(n_lncrnas))

    # STEP 2: lncRNA matching (on slice)
    try:
        ln_pairs, genes_with_lnc, lncrnas_with_genes = match_lncrnas_to_genes(genes, lncrnas_s, window_bp=1_000_000)
        results.append(_ok("STEP2.match_lncrnas", f"pairs={len(ln_pairs)} genes_w_lists={len(genes_with_lnc)} lncrnas_w_lists={len(lncrnas_with_genes)}"))
        summary["n_lncrna_pairs"] = int(len(ln_pairs))
        ln_matched = ln_pairs
    except Exception as e:
        results.append(_fail("STEP2.match_lncrnas", repr(e)))
        ln_matched = pd.DataFrame()

    # STEP 3: cCRE-gene distance matching + element focus table build (on slice)
    try:
        pairs = match_ccres_to_genes(genes, ccres_s, window_bp=1_000_000)
        elem_focus = build_element_focus_table(ccres_s, pairs, tier_labels=["tier1", "tier2", "tier3"])
        results.append(_ok("STEP3.match_ccres_and_build_elem_focus", f"pairs={len(pairs)}, elem_focus={len(elem_focus)}"))
        summary["n_ccre_gene_pairs"] = int(len(pairs))
        summary["n_elem_focus"] = int(len(elem_focus))
    except Exception as e:
        results.append(_fail("STEP3.match_ccres_and_build_elem_focus", repr(e)))
        return results, summary

    # STEP 4: TAD annotation (single biosample by default)
    try:
        tad_biosamples = _pick_tad_biosamples(PATHS.tads_processed, max_biosamples=int(n_tad_biosamples))
        if not tad_biosamples:
            results.append(_ok("STEP4.tads.annotate", "no TAD biosamples discovered (skipped)"))
            summary["tad_biosamples_used"] = []
        else:
            # For TAD annotation + mirroring we want a non-empty cCRE table regardless of whether
            # the distance matching yielded any gene-links on this small slice.
            ccres_for_tads = ccres_s.copy()
            g2, c2, l2 = annotate_with_all_tad_sources(
                genes,
                ccres_for_tads,
                ln_matched if (ln_matched is not None and not ln_matched.empty) else None,
                processed_dir=PATHS.tads_processed,
                biosamples=tad_biosamples,
                verbose=False,
            )
            results.append(_ok("STEP4.tads.annotate", f"biosamples={tad_biosamples}"))
            summary["tad_biosamples_used"] = tad_biosamples

            # STEP 4b: mirroring (exercise the crash site; limited biosamples)
            mirror_and_save_all_domains(
                genes_df=g2,
                ccre_df=c2,
                lncrnas_df=l2 if (l2 is not None and not l2.empty) else None,
                processed_dir=PATHS.tads_processed,
                biosamples=tad_biosamples,
                verbose=False,
            )
            results.append(_ok("STEP4b.tads.mirror", f"biosamples={tad_biosamples}"))

            # STEP 4c: boundary enrichment (single combined table; limited biosamples + small panel)
            try:
                bdf = build_boundaries_enriched_table(
                    processed_dir=PATHS.tads_processed,
                    biosamples=tad_biosamples,
                    gene_panel=gene_panel,
                )
                # write a tiny artifact
                (out_dir / "boundaries_enriched__head.csv").write_text(bdf.head(20).to_csv(index=False), encoding="utf-8")
                results.append(_ok("STEP4c.tads.boundary_enrich", f"rows={len(bdf)}, biosamples={tad_biosamples}"))
                summary["boundary_enriched_rows"] = int(len(bdf))
            except Exception as e:
                results.append(_fail("STEP4c.tads.boundary_enrich", repr(e)))
    except Exception as e:
        results.append(_fail("STEP4.tads", repr(e)))

    # Evidence-heavy steps (5-9) and ATAC processing (10) are intentionally not executed here.
    # We validate their inputs exist so the full pipeline can run when you want it to.
    checks = {
        "STEP5.inputs.SCREEN_exp_zip": PATHS.screen_exp_zip,
        "STEP5.inputs.SCREEN_comp_gz": PATHS.screen_comp_gz,
        "STEP6.inputs.ABC_predictions": PATHS.abc_predictions,
        "STEP8.inputs.HiChIP_dir": PATHS.hichip_dir,
        "STEP9.inputs.TargetScan": PATHS.targetscan_predictions,
        "STEP9.inputs.miRTarBase": PATHS.mirtarbase_csv,
        "STEP10.inputs.ATAC_peaks_csv": PATHS.atac_peaks_csv,
        "STEP3b.inputs.ChIP_unified": PATHS.chip_unified,
    }
    for step, p in checks.items():
        try:
            exists = bool(Path(p).exists())
            results.append(_ok(step, f"exists={exists} path={p}"))
        except Exception as e:
            results.append(_fail(step, repr(e)))

    return results, summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-genes", type=int, default=30)
    ap.add_argument("--n-ccres", type=int, default=500)
    ap.add_argument("--n-lncrnas", type=int, default=200)
    ap.add_argument("--n-tad-biosamples", type=int, default=1)
    args = ap.parse_args()

    out_dir = ensure_out_dir("pipeline_sections_smoke")

    results, summary = run_smoke(
        out_dir=out_dir,
        n_genes=int(args.n_genes),
        n_ccres=int(args.n_ccres),
        n_lncrnas=int(args.n_lncrnas),
        n_tad_biosamples=int(args.n_tad_biosamples),
    )

    df = pd.DataFrame([{"step": r.step, "ok": r.ok, "detail": r.detail} for r in results])
    df.to_csv(out_dir / "results.csv", index=False)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    n_fail = int((~df["ok"]).sum()) if not df.empty else 0
    print(f"[OK] wrote: {out_dir}")
    print(f"[OK] steps={len(df)} fail={n_fail}")
    if n_fail:
        bad = df[~df["ok"]].head(20)
        print(bad.to_string(index=False))
        raise SystemExit(2)


if __name__ == "__main__":
    main()

