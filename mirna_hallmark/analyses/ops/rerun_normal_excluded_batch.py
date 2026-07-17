#!/usr/bin/env python3
"""Single-process orchestrator for the Normal-like-discarded + batch-corrected full rerun.

Why python (not a bash loop): the harness kills/retries/duplicates background *bash*
wrappers, orphaning the heavy child. A single long-lived python process survives the
orphaning and runs every step sequentially. A PID lockfile makes a duplicate launch a
no-op, so a harness retry cannot create concurrent writers. All output goes to a log
file (no terminal pipe -> no SIGPIPE stall).
"""
import os, sys, subprocess, time, datetime, pathlib

REPO = "/sci/labs/michall/stavzok/APM"
PY = os.path.join(REPO, ".venv/bin/python3")
LOGDIR = pathlib.Path(REPO, "mirna_hallmark/output/logs")
LOGDIR.mkdir(parents=True, exist_ok=True)
LOCK = LOGDIR / "rerun_all.lock"
LOG = LOGDIR / f"rerun_all_{datetime.datetime.utcnow():%Y%m%d_%H%M%S}.log"

# ---- lockfile guard: refuse to run a second concurrent instance ----
if LOCK.exists():
    try:
        old = int(LOCK.read_text().split()[0])
        os.kill(old, 0)            # raises if not alive
        sys.exit(f"another rerun is live (pid {old}); exiting")
    except (ValueError, ProcessLookupError, PermissionError):
        pass                       # stale lock -> take over
LOCK.write_text(f"{os.getpid()} {LOG}\n")

STEPS = [
    [PY, "-m", "mirna_hallmark.healthy_anchor", "--build"],
    [PY, "-m", "mirna_hallmark.run_all", "--skip-mirna-cnv"],
    [PY, "-m", "mirna_hallmark.target_combined_anticorr"],
    [PY, "-m", "mirna_hallmark.cn_dosage_attribution"],
    [PY, "-m", "mirna_hallmark.mirna_state_class", "--healthy-anchor"],
    [PY, "-m", "mirna_hallmark.analyses.cross_state.cross_state_coupling"],
    [PY, "-m", "mirna_hallmark.analyses.misc.mirna_comovement"],
    [PY, "-m", "mirna_hallmark.nmf_sample_signatures"],
    [PY, "-m", "mirna_hallmark.decoupling_validation"],
    [PY, "-m", "mirna_hallmark.analyses.cross_state.cross_state_landscape", "--healthy-anchor"],
    [PY, "-m", "mirna_hallmark.analyses.cross_state.cross_state_deep_dive", "--gene", "PTEN", "--mirna", "hsa-miR-200c-3p", "--healthy-anchor"],
    [PY, "-m", "mirna_hallmark.analyses.mir301.mir301_family_depth"],
    [PY, "-m", "mirna_hallmark.family_normal_reference"],
    [PY, "-m", "mirna_hallmark.mir301_family_network"],
    [PY, "-m", "mirna_hallmark.analyses.cross_state.cross_state_expression_panels", "--preset", "landscape_extended", "--out-name", "landscape_extended"],
    [PY, "-m", "mirna_hallmark.analyses.cross_state.cross_state_expression_panels", "--preset", "coupling_focus", "--out-name", "coupling_focus"],
    [PY, "-m", "mirna_hallmark.analyses.cross_state.cross_state_expression_panels", "--preset", "family", "--out-name", "family"],
    [PY, "-m", "mirna_hallmark.nat_tumor_umap"],
    [PY, "-m", "mirna_hallmark.geneset_architecture", "--all-hallmarks"],
    [PY, "-m", "mirna_hallmark.analyses.nmf.within_pathway_nmf"],
    [PY, "-m", "mirna_hallmark.eval.cptac_validation", "--cohort", "tcga105"],
    [PY, "-m", "mirna_hallmark.eval.cptac_validation", "--cohort", "prospective"],
    [PY, "-m", "mirna_hallmark.cptac_resid_survivors"],
    [PY, "-m", "mirna_hallmark.analyses.cptac.cptac_orphan_discovery"],
    [PY, "-m", "mirna_hallmark.cptac_target_specificity"],
    [PY, "-m", "mirna_hallmark.analyses.cptac.cptac_orphan_confound_pilot"],
    [PY, "-m", "mirna_hallmark.cptac_orphan_evidence_table"],
    [PY, "-m", "mirna_hallmark.eval.buffa_validation"],
]


def main():
    with open(LOG, "w", buffering=1) as f:
        def log(m):
            f.write(m + "\n"); f.flush()
        log(f"START {datetime.datetime.utcnow().isoformat()}Z  pid={os.getpid()}  LOG={LOG}")
        log(f"n_steps={len(STEPS)}")
        for i, step in enumerate(STEPS, 1):
            name = " ".join(step[2:]) if len(step) > 2 else " ".join(step)
            log(f"\n######## {datetime.datetime.utcnow():%H:%M:%S}  STEP {i}/{len(STEPS)}: {name} ########")
            t0 = time.time()
            try:
                rc = subprocess.run(step, cwd=REPO, stdout=f, stderr=subprocess.STDOUT).returncode
                log(f"---- step {i} rc={rc} elapsed={time.time()-t0:.0f}s"
                    + ("" if rc == 0 else "  !!!! NONZERO (continuing)"))
            except Exception as e:  # noqa: BLE001
                log(f"---- step {i} EXCEPTION {type(e).__name__}: {e} (continuing)")
        log(f"\nALL DONE {datetime.datetime.utcnow().isoformat()}Z")
    try:
        LOCK.unlink()
    except OSError:
        pass


if __name__ == "__main__":
    main()
