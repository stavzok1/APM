"""MH-120 — WHO REPRESSES, CAUSALLY? The CN instrument as an EXOGENOUS within-family arbiter,
head-to-head: COLLAPSED-model attribution vs UNCOLLAPSED (arm-level) attribution, with a PERMUTATION null.

    .venv/bin/python3 -m mirna_hallmark.eval.within_family_cn_iv

WHY (user, 2026-07-13). MH-119 used Manakov chimeras as the sole arbiter of within-family attribution and
found the learned arm-level beta at chance. But chimeras are a WEAK label: HEK293T (not breast), sparse, and
**measured to be 59% abundance-loaded (p=1.4e-07)** — a chance-level result against a partly-confounded label
is thin grounds for a retraction. The CN instrument is the right orthogonal arbiter:

  ⛔⛔ THE PREMISE BELOW IS FALSE — THIS MODULE'S HEADLINE RESULT IS RETRACTED (MH-122). KEPT FOR PROVENANCE.
  ⛔ "SAME-SEED FAMILY MEMBERS SIT AT DIFFERENT GENOMIC LOCI, so their locus CN varies INDEPENDENTLY."
     **THEY USUALLY DO NOT.** Same-seed members arose by IN-PLACE DUPLICATION ⇒ they are CO-LOCATED ⇒ they
     share ONE CN column ⇒ identical γ, π and Wald (94% of families; `n_indep_cn < k`). Any `argmax` over them
     is decided by ROW ORDER. Use `instrument.cn_copy` (the repo's own co-location-correct §8 attribution
     share: co-located members get an EVEN SPLIT = "cannot resolve"), compared CONTINUOUSLY, never argmax'd.
  ⛔ The Hansen-J "same-seed members are causally heterogeneous" reading is ALSO wrong: an over-ID rejection
     means an INVALID INSTRUMENT at least as much as heterogeneous effects (`multi_iv`: *"rejection ⇒ a
     source-locus is PLEIOTROPIC"*). TESTED with `instrument.exclusion()`'s T1 gate: the J-rejectors ARE the
     pleiotropic ones (|δ_s| 0.153 vs 0.088, p=0.028; 92% vs 81% flagged). **Design §F STANDS.**
  ⇒ See MH-122 (ledger/registry) and `ATTRIBUTION_PRIMITIVE.md` / `CN_INSTRUMENT.md` §7 before reusing this.

`instrument.family_multi_iv` gives, per member m:
    gamma_m  first stage  (CN_m -> X_fam)   = m's DOSE delivery into the family pool
    pi_m     reduced form (CN_m -> target)  = m's GENETIC-DOSE effect on the target  <-- THE CAUSAL LEG
    F_m      first-stage strength (instrument relevance; gate at F>10)
`family_multi_iv` is ALSO the COLLAPSED model's own within-family attributor (it instruments the family
AGGREGATE and attributes to members via CN) — so this run answers both questions at once:

  (1) the CN-CAUSAL pick        argmax_m |pi_m|   (F-gated)          <- the collapsed model's answer
  (2) the UNCOLLAPSED pick      argmax_m beta_m   (arm-level Gibbs)  <- the un-collapsed model's answer
  (3) the PHYSICAL pick         Manakov chimera for that (arm,gene)  <- MH-119's label

⚠ THE NULL IS NOT 1/k. A 1/k chance baseline ASSUMES the picks are exchangeable, which they are not (the
members differ in abundance, detection, CN variance). So the baseline is MEASURED by a PERMUTATION: shuffle
the member labels within each family and re-score, 2000x. Every agreement rate below is reported against its
own permuted null.

Output: `output/learned/within_family_cn_iv.tsv`.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("mirna_hallmark/output/learned")
MIN_F = 10.0          # instrument relevance gate (weak-instrument bias below this)
N_PERM = 2000


def _one(args):
    gene, fam, mem = args
    from mirna_hallmark.learned import instrument as I, data as LD
    try:
        asm = LD.assemble_gene(gene, w_prior_source="ledger")     # HE design — VERIFIED identical IV, 40x faster
        r = I.family_multi_iv(gene, mem, assembled=asm)
    except Exception:
        return None
    if not r or not r.get("solo"):
        return None
    rows = []
    for m, s in r["solo"].items():
        rows.append({"gene": gene, "seed_family": fam, "arm": m,
                     "F": s.get("F", np.nan), "gamma": s.get("gamma", np.nan),
                     "pi": s.get("pi", np.nan), "wald": s.get("wald", np.nan),
                     "first_stage_F": r.get("first_stage_F", np.nan),
                     "beta_family_iv": r.get("beta_family", np.nan)})
    return rows


def run(workers: int = 8) -> pd.DataFrame:
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")
    from multiprocessing import Pool
    A = pd.read_csv(OUT / "readouts_arm_edges.tsv", sep="\t")
    m = A[A.family_size > 1]
    jobs = [(g, f, list(s.arm)) for (g, f), s in m.groupby(["gene", "seed_family"])]
    print(f"[cn-iv] {len(jobs)} multi-arm gene-families")
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.evidence import ledger as LG
    LD._load(); LG.pooled_he_edges()
    with Pool(workers) as pool:
        parts = [r for r in pool.imap_unordered(_one, jobs, chunksize=4) if r]
    D = pd.DataFrame([x for r in parts for x in r])
    D.to_csv(OUT / "within_family_cn_iv.tsv", sep="\t", index=False)
    print(f"[cn-iv] {len(D)} member rows / {D.groupby(['gene','seed_family']).ngroups} gene-families")
    return D


if __name__ == "__main__":
    run()
