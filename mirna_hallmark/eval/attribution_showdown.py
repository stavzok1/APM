"""MH-122 — COLLAPSED vs UNCOLLAPSED ATTRIBUTION, using the SUBPROJECT'S OWN DOCTRINE AND ESTIMATORS.

    .venv/bin/python3 -m mirna_hallmark.eval.attribution_showdown

WHY THIS EXISTS. MH-118..121 compared the two models using estimators I INVENTED, while the subproject already
had a documented attribution doctrine I had not read. This module redoes the comparison with the CANONICAL
tools. The doctrine (`ATTRIBUTION_IDENTITY_VS_MAGNITUDE.md` §0, `LEARNED_MODEL_METHODS.md` §4/§5/§16,
`DESIGN_RESPONSE` §F):

  * **IDENTITY** = **Shapley/LMG credit for R²(X·M, Y)** (`attribution.shapley_identity`) — "who fairly OWNS
    the coupling, **splitting shared credit under collinearity**". This is THE readout for collinear units.
    ⚠ NOT `readouts.share` (= β_f/Σβ). For the ADDITIVE value function Shapley IS trivially β_f — it splits
    NOTHING under collinearity. Using it on near-collinear same-seed arms was my core estimator error.
  * **MAGNITUDE** = the realized budget `M(f)·X_fam` (`states.budget_shift`) — abundance-IN.
  * **CANONICAL COEFFICIENT** = **BAGGED NNLS** (`states._bagged_nnls` / `canonical_M`), *because*
    "single-lasso coefficients are unstable under collinearity (corr 0.03)"; bagged NNLS is corr ~0.99.
  * **DESIGN §F — THE CLAIM UNDER TEST:** *"Same-seed arms are near-collinear ⇒ **the identified estimand is
    family→gene** (arm = **nomination**)"*; the per-arm split is *"a nomination with uncertainty, **not a
    coefficient**"*, apportioned by abundance (`canonical_M(arm_level=True)`) or by chimeric provenance
    (`chimeric_evidence`, the L2 nominator). METHODS §16 already tested a sequence/chimeric arm nomination and
    found it **INERT** ("NO regression value ... identity_payoff_genome 0% flip").

SO THE SCIENTIFIC QUESTION IS SHARP: **is the within-family arm split really unidentified (Design §F), or can
an ARM-level fit resolve it beyond the abundance nomination?** MH-121 says the latter (arm β agrees with the
exogenous CN attribution at ρ=+0.285, p=0.005, GIVEN dose) — which contradicts §F and §16. Settle it with the
doctrine's own estimators.

WHAT THIS RUNS, per gene:
  COLLAPSED    M_fam = bagged NNLS on X_fam            → LMG-Shapley identity over FAMILIES
               arm split = the doctrine's NOMINATION (chimeric L2, else expression×loading)
  UNCOLLAPSED  M_arm = bagged NNLS on the ARM design   → LMG-Shapley identity over ARMS
  Then:  (A) BETWEEN-FAMILY ranking agreement (uncollapsed arms re-aggregated to families)
         (B) WITHIN-FAMILY arm ranking: nomination vs uncollapsed identity, scored against the exogenous
             CN attribution (`instrument.cn_copy`, dose-conditioned), AGO loading and scanMiR K_D
         (C) POSITIVE CONTROLS from the docs (see CONTROLS below)

Output: `output/learned/attribution_showdown_{arms,genes}.tsv`.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("mirna_hallmark/output/learned")
N_PERM = 300
N_BOOT = 40

# --- POSITIVE CONTROLS the docs themselves name (ATTRIBUTION_IDENTITY_VS_MAGNITUDE §0/§11/§12, METHODS §16,
#     DESIGN_RESPONSE §1h). Each is an expectation the attribution MUST reproduce, or explain.
CONTROLS = {
    "ESR1":  "quiet owner miR-18 (identity .54 / budget .20) vs LOUD PASSENGER miR-22 (budget .26 / identity .01)",
    "PTEN":  "identity≈magnitude (rank-corr +0.91); curated drivers miR-21/103a/181/182/96; §16 miR-19/26",
    "NCOA3": "§12 REALITY CHECK: surfaced specialist miR-137 does NOT couple (NS); the ABUNDANT miR-17-5p does (-0.25)",
    "CCND1": "§16 canonical specialist miR-15/16",
    "ZEB1":  "§16 canonical miR-200/141; the COLLINEAR case (miR-200bc/429 vs miR-141/200a)",
    "ABCA1": "headline attribution case → miR-33b",
    "CASP8": "headline attribution case → miR-376a",
    "SPRY2": "headline attribution case → miR-24",
}


def _resid(V, Cm):
    from sklearn.linear_model import LinearRegression
    return V - LinearRegression().fit(Cm, V).predict(Cm)


def _bagged(Xz: np.ndarray, yr: np.ndarray, *, n_boot=N_BOOT, seed=0) -> np.ndarray:
    """The doctrine's canonical coefficient (METHODS §5): bagged z-scored NNLS. Stable under collinearity
    where a single lasso is not (corr 0.03 → 0.99)."""
    from scipy.optimize import nnls
    rng = np.random.default_rng(seed)
    n, p = Xz.shape
    acc = np.zeros(p)
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        c, _ = nnls(Xz[idx], yr[idx])
        acc += c
    return acc / n_boot


def _prep(gene):
    from mirna_hallmark.learned import data as LD, families as FAM
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    if X is None or X.shape[1] < 2:
        return None
    fam = FAM.family_of(pd.Index(X.columns)).reindex(X.columns)
    Xf, _, _ = FAM.collapse_by_family(X, w, fam)
    Cm = C.to_numpy(float)
    yr = -_resid(Y.to_numpy(float), Cm)                       # sign: repression → POSITIVE weight (§5)
    def z(M):
        A = M.to_numpy(float)
        A = _resid(A, Cm)
        sd = A.std(0)
        A = (A - A.mean(0)) / (sd + 1e-9)
        A[:, sd < 0.1] = 0.0                                  # §5 VARIANCE FLOOR
        return A
    return dict(gene=gene, yr=yr, Xa=z(X), Xf=z(Xf), arms=list(X.columns),
                fams=list(Xf.columns), fam_of=fam)


def _one(gene):
    from mirna_hallmark.learned import attribution as AT, chimeric_evidence as CE, data as LD
    try:
        d = _prep(gene)
    except Exception:
        return None
    if d is None:
        return None
    out = []
    for level, Xz, cols in [("family", d["Xf"], d["fams"]), ("arm", d["Xa"], d["arms"])]:
        M = _bagged(Xz, d["yr"])                              # CANONICAL coefficient (bagged NNLS)
        phi = AT.shapley_identity(Xz, M, d["yr"], n_perm=N_PERM).clip(min=0)   # LMG identity (collinearity-fair)
        ident = phi / phi.sum() if phi.sum() > 0 else phi
        mag = M * np.abs(Xz).mean(0)                          # realized budget (abundance-IN)
        mag = mag / mag.sum() if mag.sum() > 0 else mag
        for c, m, i, g in zip(cols, M, ident, mag):
            out.append({"gene": gene, "level": level, "unit": c, "M": m, "identity": i, "magnitude": g})
    R = pd.DataFrame(out)
    # --- the doctrine's ARM NOMINATION inside each family (L2 chimeric, else expression×loading)
    noms = []
    for f, mem in pd.Series(d["fam_of"]).groupby(d["fam_of"]):
        members = list(mem.index)
        if len(members) < 2:
            continue
        try:
            ev = CE.evidence(gene, members)                   # pooled chimeric provenance
        except Exception:
            ev = pd.Series(dtype=float)
        src = "chimeric"
        if ev is None or len(ev) == 0 or float(np.nansum(ev.to_numpy(float))) <= 0:
            X = LD._load()["X"]
            ex = X.reindex(members).median(axis=1)
            ev = ex / ex.sum() if ex.sum() > 0 else pd.Series(1.0 / len(members), index=members)
            src = "expr_load"
        ev = ev.reindex(members).fillna(0.0)
        s = ev.sum()
        ev = ev / s if s > 0 else pd.Series(1.0 / len(members), index=members)
        for m in members:
            noms.append({"gene": gene, "seed_family": f, "arm": m,
                         "nomination": float(ev.get(m, 0.0)), "nom_src": src})
    return R, pd.DataFrame(noms)


def run(genes=None, workers: int = 8):
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")
    from multiprocessing import Pool
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.evidence import ledger as LG
    LD._load(); LG.pooled_he_edges()
    if genes is None:
        A = pd.read_csv(OUT / "readouts_arm_edges.tsv", sep="\t")
        genes = sorted(set(A[A.family_size > 1].gene) | set(CONTROLS))    # multi-arm genes + the controls
    print(f"[showdown] {len(genes)} genes (incl. {len(CONTROLS)} documented positive controls)")
    with Pool(workers) as p:
        res = [r for r in p.imap_unordered(_one, genes, chunksize=4) if r]
    U = pd.concat([r[0] for r in res], ignore_index=True)
    N = pd.concat([r[1] for r in res if len(r[1])], ignore_index=True)
    U.to_csv(OUT / "attribution_showdown_units.tsv", sep="\t", index=False)
    N.to_csv(OUT / "attribution_showdown_nominations.tsv", sep="\t", index=False)
    print(f"[showdown] {len(U)} unit rows / {U.gene.nunique()} genes | {len(N)} multi-arm nomination rows")
    return U, N


if __name__ == "__main__":
    run()
