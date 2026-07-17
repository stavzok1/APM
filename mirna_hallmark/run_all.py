"""Orchestrate the full miRNA x Hallmark subproject end to end.

Run order (each step is independently runnable as ``python -m mirna_hallmark.<step>``):
  1. build_edges            -- regenerate Hallmark-scoped miRTarBase edges
  2. ago_gate               -- per-sample AGO/RISC capacity + gate + stratum tables
  3. mirna_locus_cnv        -- miRNA-universe CNV landscape (SLOW; cached) + concordance
  4. stratum_characterization -- per-Hallmark target CNV + miRNA expression by stratum
  5. hallmark_interaction   -- AGO-gated pressure, anti-correlation, enrichment
  6. subtype_contrasts      -- luminal vs non-luminal contrasts, differential gene
                               pressure, per-Hallmark miRNA impact, immune axis + IRF1
  7. visibility_archetype_contrasts -- actual cold_Basal / hot_Luminal (APM archetypes),
                               all 50 Hallmarks + tumor-suppressor gene route panel
  8. pam50_gene_resolution  -- actual LumA/LumB/Her2/Basal per-gene pressure and
                               top driver miRNAs (Normal excluded)
  9. robustness_checks      -- grant Aims 1-2: proliferation-adjusted coupling +
                               curation-bias (evidence/binary/TargetScan) hub robustness
                               (needs step 5's cached hallmark pressure matrix)
 10. coupling_predictor_comparison -- methods/sensitivity: does the pressure construction
                               (share/z/level) change realized coupling vs raw abundance,
                               at edge and gene resolution (MH-44); ~6 min, --skip- to omit

``mirna_locus_cnv`` does a one-time ~12 min cohort CNV extraction (then cached);
use ``--skip-mirna-cnv`` to skip it when the cache already exists and you only
need the Hallmark-interaction layer.

Writes ``output/run_manifest.json`` with per-step status + timing.
"""

from __future__ import annotations

import argparse
import json
import time
import traceback
from datetime import datetime, timezone
from typing import Optional
from pathlib import Path

from mirna_hallmark import config as C


def _step(name: str, fn, manifest: dict, *, fatal: bool = True, **kwargs) -> bool:
    print(f"\n{'=' * 70}\n[run_all] STEP: {name}\n{'=' * 70}")
    t0 = time.time()
    entry = {"started_utc": datetime.now(timezone.utc).isoformat()}
    try:
        fn(**kwargs)
        entry["status"] = "ok"
        ok = True
    except Exception as exc:  # noqa: BLE001
        entry["status"] = "error"
        entry["error"] = f"{type(exc).__name__}: {exc}"
        traceback.print_exc()
        ok = False
    entry["elapsed_s"] = round(time.time() - t0, 1)
    manifest["steps"][name] = entry
    if not ok and fatal:
        raise RuntimeError(f"[run_all] fatal failure in step '{name}'")
    return ok


def run_all(
    *,
    force_edges: bool = False,
    skip_mirna_cnv: bool = False,
    include_tnrc6: Optional[bool] = None,
    skip_predictor_comparison: bool = False,
) -> dict:
    C.ensure_output_dirs()
    # ``include_tnrc6=None`` -> follow the canonical config default (``C.AGO_GATE``,
    # the MH-46 TNRC6 co-limitation gate). Pass an explicit bool only to override.
    tnrc6 = C.AGO_GATE.include_tnrc6 if include_tnrc6 is None else include_tnrc6
    manifest = {
        "subproject": "mirna_hallmark",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "params": {"force_edges": force_edges, "skip_mirna_cnv": skip_mirna_cnv,
                   "include_tnrc6": tnrc6,
                   "skip_predictor_comparison": skip_predictor_comparison},
        "steps": {},
    }

    from mirna_hallmark import (
        build_edges,
        ago_gate,
        stratum_characterization,
        hallmark_interaction,
        subtype_contrasts,
        pam50_gene_resolution,
    )

    _step("build_edges", lambda: build_edges.run(force=force_edges), manifest)

    _step("ago_gate", lambda: ago_gate.run(include_tnrc6=tnrc6), manifest, fatal=False)

    if not skip_mirna_cnv:
        from mirna_hallmark import mirna_locus_cnv
        _step("mirna_locus_cnv", lambda: mirna_locus_cnv.run(force=False), manifest, fatal=False)
    else:
        manifest["steps"]["mirna_locus_cnv"] = {"status": "skipped"}

    _step("stratum_characterization", lambda: stratum_characterization.run(), manifest, fatal=False)
    _step("hallmark_interaction", lambda: hallmark_interaction.run(include_tnrc6=tnrc6),
          manifest, fatal=False)
    _step("subtype_contrasts", lambda: subtype_contrasts.run(include_tnrc6=tnrc6),
          manifest, fatal=False)
    _step("pam50_gene_resolution", lambda: pam50_gene_resolution.run(include_tnrc6=tnrc6),
          manifest, fatal=False)
    from mirna_hallmark import visibility_archetype_contrasts
    _step("visibility_archetype_contrasts",
          lambda: visibility_archetype_contrasts.run(include_tnrc6=tnrc6),
          manifest, fatal=False)

    from mirna_hallmark import robustness_checks
    _step("robustness_checks", lambda: robustness_checks.run(), manifest, fatal=False)

    from mirna_hallmark import hybrid_pressure
    _step("hybrid_pressure", lambda: hybrid_pressure.run(), manifest, fatal=False)

    if not skip_predictor_comparison:
        from mirna_hallmark import coupling_predictor_comparison
        _step("coupling_predictor_comparison", lambda: coupling_predictor_comparison.run(),
              manifest, fatal=False)
    else:
        manifest["steps"]["coupling_predictor_comparison"] = {"status": "skipped"}

    manifest["ended_utc"] = datetime.now(timezone.utc).isoformat()
    (C.OUTPUT_ROOT / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"\n[run_all] wrote run manifest -> {C.OUTPUT_ROOT / 'run_manifest.json'}")
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force-edges", action="store_true", help="Rebuild miRTarBase Hallmark edges")
    ap.add_argument("--skip-mirna-cnv", action="store_true", help="Skip the slow miRNA-universe CNV step")
    ap.add_argument("--include-tnrc6", action="store_true", default=None,
                    help="Force TNRC6A/B/C into the AGO gate (default: follow config.AGO_GATE)")
    ap.add_argument("--no-tnrc6", action="store_true",
                    help="Exclude TNRC6A/B/C from the AGO gate (override the config default)")
    ap.add_argument("--skip-predictor-comparison", action="store_true",
                    help="Skip the coupling predictor-construction sensitivity (MH-44; ~6 min)")
    args = ap.parse_args()
    include_tnrc6 = False if args.no_tnrc6 else args.include_tnrc6  # None -> config default
    run_all(
        force_edges=args.force_edges,
        skip_mirna_cnv=args.skip_mirna_cnv,
        include_tnrc6=include_tnrc6,
        skip_predictor_comparison=args.skip_predictor_comparison,
    )


if __name__ == "__main__":
    main()
