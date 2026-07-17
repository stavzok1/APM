#!/usr/bin/env python3
"""Parallel orchestrator for the Normal-like-discarded + batch-corrected full rerun.

Same 28 modules as ``rerun_normal_excluded_batch.py``, but the heavy downstream tail
runs CONCURRENTLY. The two-step PREFIX (healthy anchor + the core spine) must finish
first because almost everything reads its edges/pressure/gate outputs. After that, the
remaining work is organized into independent CHAINS: steps within a chain run serially
(they have real data dependencies), but the chains run in parallel up to ``MAX_WORKERS``.

Wall-clock therefore drops from sum(all steps) to ~ prefix + longest chain.

Like the serial version: a single long-lived python process (survives harness orphaning),
a PID lockfile (no concurrent writers), all output to a log file, failures are logged and
do not abort the rest.
"""
import datetime
import os
import pathlib
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

REPO = "/sci/labs/michall/stavzok/APM"
PY = os.path.join(REPO, ".venv/bin/python3")
LOGDIR = pathlib.Path(REPO, "mirna_hallmark/output/logs")
LOGDIR.mkdir(parents=True, exist_ok=True)
LOCK = LOGDIR / "rerun_all.lock"
LOG = LOGDIR / f"rerun_parallel_{datetime.datetime.utcnow():%Y%m%d_%H%M%S}.log"
MAX_WORKERS = int(os.environ.get("RERUN_WORKERS", str(min(12, (os.cpu_count() or 8)))))


def M(*a):  # module step
    return [PY, "-m", "mirna_hallmark." + a[0], *a[1:]]


# Must finish before any chain starts (everything reads the spine).
PREFIX = [
    M("healthy_anchor", "--build"),
    M("run_all", "--skip-mirna-cnv"),
]

# Independent chains (parallel across, serial within). Grouping is conservative:
# steps that share a data dependency stay in one chain.
CHAINS = [
    [M("target_combined_anticorr")],
    [M("cn_dosage_attribution")],
    [M("geneset_architecture", "--all-hallmarks"), M("within_pathway_nmf")],
    [M("mirna_comovement"), M("nmf_sample_signatures")],
    [M("mir301_family_depth"), M("mir301_family_network"), M("family_normal_reference")],
    [M("buffa_validation")],
    # cross-state: ONLY state_class -> coupling -> decoupling is a real chain (decoupling
    # reads both their outputs). landscape / deep-dive / panels / umap are independent leaves
    # (they read the spine + state-bundle loaders, not these modules' outputs) -> own chains.
    [
        M("mirna_state_class", "--healthy-anchor"),
        M("cross_state_coupling"),
        M("decoupling_validation"),
    ],
    [M("cross_state_landscape", "--healthy-anchor")],
    [M("cross_state_deep_dive", "--gene", "PTEN", "--mirna", "hsa-miR-200c-3p", "--healthy-anchor")],
    [M("cross_state_expression_panels", "--preset", "landscape_extended", "--out-name", "landscape_extended")],
    [M("cross_state_expression_panels", "--preset", "coupling_focus", "--out-name", "coupling_focus")],
    [M("cross_state_expression_panels", "--preset", "family", "--out-name", "family")],
    [M("nat_tumor_umap")],
    # CPTAC chain: the two proteome loads first, then everything that residualizes on them
    [
        M("cptac_validation", "--cohort", "tcga105"),
        M("cptac_validation", "--cohort", "prospective"),
        M("cptac_resid_survivors"),
        M("cptac_orphan_discovery"),
        M("cptac_target_specificity"),
        M("cptac_orphan_confound_pilot"),
        M("cptac_orphan_evidence_table"),
    ],
]


def main():
    # lock guard (same as serial version)
    if LOCK.exists():
        try:
            old = int(LOCK.read_text().split()[0])
            os.kill(old, 0)
            sys.exit(f"another rerun is live (pid {old}); exiting")
        except (ValueError, ProcessLookupError, PermissionError):
            pass
    LOCK.write_text(f"{os.getpid()} {LOG}\n")

    with open(LOG, "w", buffering=1) as f:
        def log(m):
            f.write(m + "\n"); f.flush()

        def run_step(step, tag):
            name = " ".join(step[2:])
            t0 = time.time()
            log(f"[{tag}] START {datetime.datetime.utcnow():%H:%M:%S} {name}")
            try:
                rc = subprocess.run(step, cwd=REPO, stdout=f, stderr=subprocess.STDOUT).returncode
            except Exception as e:  # noqa: BLE001
                log(f"[{tag}] EXCEPTION {type(e).__name__}: {e} (continuing)")
                return
            log(f"[{tag}] DONE rc={rc} elapsed={time.time()-t0:.0f}s {name}"
                + ("" if rc == 0 else "  !!!! NONZERO"))

        def run_chain(idx, chain):
            tag = f"chain{idx}"
            for step in chain:
                run_step(step, tag)
            return idx

        skip_prefix = "--skip-prefix" in sys.argv
        log(f"START {datetime.datetime.utcnow().isoformat()}Z pid={os.getpid()} workers={MAX_WORKERS}")
        log(f"PREFIX: {len(PREFIX)} serial steps{' (SKIPPED - reuse existing spine)' if skip_prefix else ''}; "
            f"then {len(CHAINS)} parallel chains")

        # ---- serial prefix (reuses an already-built spine when --skip-prefix) ----
        if not skip_prefix:
            for step in PREFIX:
                run_step(step, "prefix")

        # ---- parallel chains ----
        t0 = time.time()
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(run_chain, i, c): i for i, c in enumerate(CHAINS)}
            for fut in as_completed(futs):
                log(f"  chain {fut.result()} complete")
        log(f"\nALL CHAINS DONE in {time.time()-t0:.0f}s")
        log(f"ALL DONE {datetime.datetime.utcnow().isoformat()}Z")

    try:
        LOCK.unlink()
    except OSError:
        pass


if __name__ == "__main__":
    main()
