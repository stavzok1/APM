"""Sharded CLI driver for the paired-realization engine (`learned/realization.py`).

Independent-subprocess shards (NO multiprocessing.Pool — it deadlocked three ways on this NFS box;
CLAUDE.md axiom 3a). Each shard runs ladder (gene+edge) + null + NAT-adjudication on a gene slice and
writes its own part files; `--finalize` concats + summarizes.

  per shard:  python -m mirna_hallmark.analyses.progression.realize_run --shard i --nshards N --m-ref complement
  finalize:   python -m mirna_hallmark.analyses.progression.realize_run --finalize
Driver script: run the shards in parallel via a bash loop, then --finalize (see run ledger).
"""
from __future__ import annotations

import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")

import argparse
import glob
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

from mirna_hallmark import config as C
from mirna_hallmark.learned import realization as R

OUT = Path(C.OUTPUT_ROOT) / "learned" / "realization"


def _universe():
    from mirna_hallmark.learned import data as LD
    return sorted(set(LD.D.high_evidence_edges()["gene"].unique()))


def run_shard(shard: int, nshards: int, m_ref: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = _universe()[shard::nshards]                     # round-robin slice
    tag = f"{shard:02d}"
    R.realize_gene_edge(genes, m_ref=m_ref).to_csv(OUT / f"ladder_{tag}.tsv", sep="\t", index=False)
    R.realize_null(genes, m_ref=m_ref).to_csv(OUT / f"null_{tag}.tsv", sep="\t", index=False)
    R.nat_reference_adjudication(genes, m_ref=m_ref).to_csv(OUT / f"nat_{tag}.tsv", sep="\t", index=False)
    R.dose_shift_edge(genes, m_ref=m_ref).to_csv(OUT / f"dose_edge_{tag}.tsv", sep="\t", index=False)  # DOSE axis
    print(f"[shard {tag}] {len(genes)} genes -> ladder/null/nat/dose parts")


def run_shard_phase2(shard: int, nshards: int, m_ref: str) -> None:
    """Phase-2 shard: Res-3 family + Res-5 between-family Shapley on the paired Δ (the expensive M-fit + bootstrap
    part). Written to fam_/btwfam_ part files; owner-convergence + retention joins run once at --finalize-phase2."""
    OUT.mkdir(parents=True, exist_ok=True)
    genes = _universe()[shard::nshards]
    tag = f"{shard:02d}"
    R.realize_family(genes, m_ref=m_ref).to_csv(OUT / f"fam_{tag}.tsv", sep="\t", index=False)
    R.realize_between_family(genes, m_ref=m_ref).to_csv(OUT / f"btwfam_{tag}.tsv", sep="\t", index=False)
    print(f"[phase2 shard {tag}] {len(genes)} genes -> fam/btwfam parts")


def finalize_phase2() -> None:
    fam = _concat("fam_*.tsv"); bf = _concat("btwfam_*.tsv")
    fam.to_csv(OUT / "realization_family.tsv", sep="\t", index=False)
    bf.to_csv(OUT / "realization_between_family.tsv", sep="\t", index=False)
    print(f"=== Phase 2 — family + between-family ({fam.gene.nunique() if len(fam) else 0} genes) ===")
    if len(fam):
        fa = fam[fam.resolution == "family_agg"]; fm = fam[fam.resolution == "family"]
        print(f"family_agg (nonlinear-pooled aggregate): n={len(fa)} | median ρ_adj {fa.rho_adj.median():.3f}")
        print(f"family (per-family):                     n={len(fm)} | median ρ_adj {fm.rho_adj.median():.3f} "
              f"| ρ_adj<0 in {100*(fm.rho_adj<0).mean():.0f}% | median own_specific_frac {fm.own_specific_frac.median():.2f}")
    if len(bf):
        R.owner_convergence()
    R.retention_realization()
    # ⭐ integrated progression cards (edge + gene) — now with the Phase-2 realization/owner columns
    # annotate ONCE, after both are built (each call would otherwise re-run the full join)
    ec = R.edge_card(annotate=False); gc = R.gene_card()   # gene_card's annotate covers edge too
    print(f"\n[progression cards] edge {ec.shape} · gene {gc.shape}")
    R.class_realization()   # cross-state class -> within-patient realization (gene-clustered, dose-ctrl, decoy-arbitrated)
    R.patient_realization_efficiency()   # FU-1: within-sample per-patient efficiency (decoy-controlled)
    R.acquired_unrealized_buffers()      # FU-2: acquired-but-unrealized gene buffers (power-controlled)
    R.hallmark_realization()             # rigorous per-program realization (class+power adjusted, decoy-arbitrated)
    print(f"-> {OUT}/ : realization_family · realization_between_family · owner_convergence · retention_realization "
          f"· edge_card · gene_card")


def drive_phase2(nshards: int = 8, m_ref: str = "complement") -> None:
    import subprocess
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    for pat in ("fam_*.tsv", "btwfam_*.tsv"):
        for f in OUT.glob(pat):
            f.unlink()
    env = {**_os.environ, "PYTHONPATH": ".", "OMP_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    print(f"[drive-phase2] launching {nshards} independent shards (m_ref={m_ref})...")
    procs = [subprocess.Popen([sys.executable, "-m", "mirna_hallmark.analyses.progression.realize_run",
                               "--shard", str(i), "--nshards", str(nshards), "--m-ref", m_ref, "--phase2"], env=env)
             for i in range(nshards)]
    rc = [p.wait() for p in procs]
    print(f"[drive-phase2] shards done (exit codes {rc}); finalizing...")
    finalize_phase2()


def _concat(pat: str) -> pd.DataFrame:
    fs = sorted(glob.glob(str(OUT / pat)))
    dfs = [pd.read_csv(f, sep="\t") for f in fs if _os.path.getsize(f) > 5]
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def finalize() -> None:
    ladder = _concat("ladder_*.tsv"); nd = _concat("null_*.tsv"); nat = _concat("nat_*.tsv")
    ladder.to_csv(OUT / "realization_ladder.tsv", sep="\t", index=False)
    nd.to_csv(OUT / "realization_null_edges.tsv", sep="\t", index=False)
    nat.to_csv(OUT / "nat_reference.tsv", sep="\t", index=False)
    if len(nd):
        R.null_summary(nd).to_csv(OUT / "realization_null.tsv", sep="\t", index=False)
    # DOSE axis: edge-level (concat shards) + arm-level global (one genome-agnostic call)
    _concat("dose_edge_*.tsv").to_csv(OUT / "dose_shift_edge.tsv", sep="\t", index=False)
    da = R.dose_shift_arm(); da.to_csv(OUT / "dose_shift_arm.tsv", sep="\t", index=False)
    if len(da):
        print(f"\nDOSE-shift ARM: n={len(da)} arms | median own_specific_frac {da.own_specific_frac.median():.2f} "
              f"(dose-shift variance that is patient-specific NAT baseline)")
    # summary to stdout
    print(f"=== paired realization — genome-wide ({ladder.gene.nunique() if len(ladder) else 0} genes) ===")
    if len(ladder):
        gg = ladder[ladder.resolution == "gene"]
        print(f"GENE-aggregate: n={len(gg)} | median ρ_raw {gg.rho_raw.median():.3f} ρ_adj {gg.rho_adj.median():.3f} "
              f"| ρ_adj<0 in {100*(gg.rho_adj<0).mean():.0f}% | composition_explained {100*gg.composition_explained.mean():.0f}%")
    if len(nd):
        print("\nNULL (REAL vs DECOY, set-level, orientation × C):")
        print(R.null_summary(nd).to_string(index=False))
    if len(nat):
        print("\nNAT-reference adjudication (mean gene-agg):")
        print(nat.groupby("reference")[["resid_var_dY", "rho_raw", "rho_adj"]].mean().round(3).to_string())
    # DESCRIPTIVE PATTERN LAYER (cheap — reuses the outputs just written; no new heavy compute)
    R.edge_pattern_table(); R.mirna_nat_retention(); R.patient_nat_identity()
    print(f"\n-> {OUT}/ : realization_ladder · realization_null · nat_reference · dose_shift_{{arm,edge}} · "
          f"master_edge_patterns · mirna_nat_retention · patient_nat_identity")


def drive(nshards: int = 8, m_ref: str = "complement") -> None:
    """Orchestrate the full sharded run in-process (replaces the scratchpad shard driver): launch `nshards`
    INDEPENDENT single-process subprocesses (axiom 3a — no in-process Pool), wait, then finalize."""
    import subprocess
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    for pat in ("ladder_*.tsv", "null_*.tsv", "nat_*.tsv", "dose_edge_*.tsv"):
        for f in OUT.glob(pat):
            f.unlink()
    env = {**_os.environ, "PYTHONPATH": ".", "OMP_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    print(f"[drive] launching {nshards} independent shards (m_ref={m_ref})...")
    procs = [subprocess.Popen([sys.executable, "-m", "mirna_hallmark.analyses.progression.realize_run",
                               "--shard", str(i), "--nshards", str(nshards), "--m-ref", m_ref], env=env)
             for i in range(nshards)]
    rc = [p.wait() for p in procs]
    print(f"[drive] shards done (exit codes {rc}); finalizing...")
    finalize()


def decoy_control() -> None:
    """The expensive (~3 min) own-vs-cohort decoy control + the shift_class stratification (MH-158 refinement)."""
    from mirna_hallmark.learned import data as LD
    genes = sorted(set(LD.D.high_evidence_edges()["gene"].unique()))
    R.nat_decoy_control(genes).to_csv(OUT / "nat_decoy_control.tsv", sep="\t", index=False)
    R.decoy_stratify()
    print(f"-> {OUT}/nat_decoy_control.tsv · nat_decoy_stratified{{,_summary}}.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shard", type=int, default=None)
    ap.add_argument("--nshards", type=int, default=8)
    ap.add_argument("--m-ref", default="complement", choices=["complement", "full", "matched"])
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--drive", action="store_true", help="orchestrate the full sharded run + finalize")
    ap.add_argument("--decoy-control", action="store_true", help="the ~3min own-vs-cohort decoy control + stratify")
    ap.add_argument("--phase2", action="store_true", help="shard: run Res-3 family + Res-5 between-family instead")
    ap.add_argument("--drive-phase2", action="store_true", help="orchestrate the Phase-2 sharded run + finalize")
    ap.add_argument("--finalize-phase2", action="store_true", help="concat Phase-2 parts + owner-convergence/retention")
    a = ap.parse_args()
    if a.drive_phase2:
        drive_phase2(a.nshards, a.m_ref)
    elif a.finalize_phase2:
        finalize_phase2()
    elif a.drive:
        drive(a.nshards, a.m_ref)
    elif a.finalize:
        finalize()
    elif a.decoy_control:
        decoy_control()
    elif a.shard is not None and a.phase2:
        run_shard_phase2(a.shard, a.nshards, a.m_ref)
    elif a.shard is not None:
        run_shard(a.shard, a.nshards, a.m_ref)
    else:
        ap.error("give --drive[-phase2], or --shard i --nshards N [--phase2], or --finalize[-phase2], or --decoy-control")


if __name__ == "__main__":
    main()
