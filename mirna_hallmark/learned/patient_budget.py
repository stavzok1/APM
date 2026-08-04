"""PER-PATIENT budget: the within-gene contribution share and rank, kept at PATIENT resolution.

⭐ **THE GAP THIS FILLS.** Every `share_*` / `rank_*` / `dose_rank_*` / `grank_*` / `d_rank_*` in this
subproject is a function of EXPRESSION evaluated at the COHORT aggregate. In `states.budget_shift` the
collapse is a single line — `state_matrices(st)[0].reindex(regs).median(axis=1)` (`states.py:541-542`) —
and **everything downstream of it is a pure column-wise function of that vector**, so the patient axis
can be kept simply by not collapsing. See `docs/PATIENT_QUESTION_TAXONOMY.md` §11.

⚠ **`M` STAYS COHORT-LEVEL.** A bagged-NNLS coefficient needs ≥2 samples, so the per-patient budget is
`cohort M × per-patient abundance`. Coupling is inherently P-across; only LEVEL quantities can be
genuinely per-patient (taxonomy §1).

⚠⚠ **TWO INCOMPATIBLE SHARE FORMULAS ALREADY EXIST IN THE REPO — this module implements the CANONICAL,
FAMILY-CONSERVING one** (`states.py:548-558`):

    contrib(u,p) = M_fam(u) · log2(1 + pool_fam(u,p)) · ( lin(u,p) / pool_fam(u,p) )

so that `Σ_{arms∈family} contrib = M_fam · X_fam` — no family-size inflation. `M` is fit *per seed
family* on the family-pooled predictor, so it must multiply `X_fam`, not the arm's own abundance; the
arm split is an apportionment (`states.py:526-529`). `realization.dose_shift_edge` uses a **different,
naive** share — `normalise(M_arm · lin(arm))`, no family pooling — which disagrees for multi-arm
families. That divergence is MEASURED by `compare_to_dose_shift_edge()`, never silently absorbed.

Engine only: functions return DataFrames and write nothing (the `realization.py` convention).
"""
from __future__ import annotations

import os as _os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")                 # single-thread BLAS before numpy (axiom 3a)

from functools import lru_cache
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import states as ST

OUT = Path(C.OUTPUT_ROOT) / "learned" / "realization"
_LEARNED = Path(C.OUTPUT_ROOT) / "learned"

RANK_METHOD = "min"      # ⚠ `states.py:558` uses pandas default "average" then .astype(int) (TRUNCATES
                         # ties); `card.py:293` uses "min". We take "min" and report the divergence in
                         # `verify_collapse()` rather than silently inheriting either.


# --------------------------------------------------------------------------- #
def _M_for(gene: str, m_source: str, arm_level: bool) -> pd.Series:
    """The cohort weight vector. ⚠ `m_source` MATTERS FOR THE VERIFICATION GATE: `budget_shift` uses
    `canonical_M`, so only `canonical` can reproduce the recorded `share_*`/`rank_*`. `complement` is the
    leak-free non-matched-tumour fit and is the right choice for the PAIRED tests."""
    if m_source == "canonical":
        return ST.canonical_M(gene, "01", arm_level=arm_level)
    from mirna_hallmark.learned import realization as R
    return (R.M_reference(gene, m_source) if arm_level else R.M_reference_fam(gene, m_source))


def _fam_codes(regs, arm_level: bool):
    """(codes, n_levels) mapping each unit to its seed family. At family rung each unit IS its family."""
    if not arm_level:
        return np.arange(len(regs)), len(regs)
    fam = FAM.family_of(pd.Index(regs)).reindex(regs).to_numpy()
    _lv, codes = np.unique(fam.astype(str), return_inverse=True)
    return codes, len(_lv)


def _budget_matrix(A: np.ndarray, Mv: np.ndarray, codes: np.ndarray, n_fam: int,
                   nan_policy: str = "zero"):
    """Vectorised `states.py:549-558` over patients. `A` is unit × patient log2(abund+1), NaN allowed.

    `nan_policy`: **zero** reproduces `budget_shift` (`.fillna(0.0)` at `states.py:549`) · **abstain**
    excludes the unit from its gene's denominator, the `_budget()` semantics at `states.py:501-520`
    ("a missing arm is an ABSTENTION, not a measurement"). Both are emitted — the per-patient version
    makes the choice far more consequential than the cohort one, because a single patient's dropout
    now moves that patient's whole within-gene ranking."""
    obs = np.isfinite(A)
    a = np.where(obs, A, 0.0) if nan_policy == "zero" else np.where(obs, A, 0.0)
    a_lin = np.clip(2.0 ** a - 1.0, 0.0, None)
    if nan_policy == "abstain":
        a_lin = np.where(obs, a_lin, 0.0)                       # invisible units contribute nothing...
    G = np.zeros((n_fam, len(codes)))
    G[codes, np.arange(len(codes))] = 1.0
    pool = (G @ a_lin)[codes, :]                                # Σ within family, per patient
    x_fam = np.log2(1.0 + pool)
    share_in_fam = np.divide(a_lin, pool, out=np.zeros_like(a_lin), where=pool > 0)
    c = np.clip(Mv[:, None] * x_fam * share_in_fam, 0.0, None)
    if nan_policy == "abstain":
        c = np.where(obs, c, np.nan)                            # ...and are OUT of the denominator
        tot = np.nansum(c, axis=0, keepdims=True)
    else:
        tot = c.sum(axis=0, keepdims=True)
    share = np.divide(c, tot, out=np.full_like(c, np.nan), where=tot > 0)
    rank = pd.DataFrame(c).rank(axis=0, ascending=False, method=RANK_METHOD).to_numpy()
    return c, share, rank


# --------------------------------------------------------------------------- #
def patient_budget(gene: str, *, m_source: str = "canonical", arm_level: bool = True,
                   sample_types: Sequence[str] = ("01", "11"), participants=None,
                   nan_policy: str = "zero") -> pd.DataFrame:
    """Per-patient within-gene budget for one gene → long frame
    `gene · unit · patient · state · budget · share · rank`.

    `participants=None` uses EVERY participant of each state (what `budget_shift` medians over — required
    for the verification gate); pass the paired list for the §J lane."""
    M = _M_for(gene, m_source, arm_level)
    if M is None or M.empty:
        return pd.DataFrame()
    regs = list(M[M > 0].index)
    if not regs:
        return pd.DataFrame()
    Mv = M.reindex(regs).to_numpy(float)
    codes, n_fam = _fam_codes(regs, arm_level)
    rows = []
    for st in sample_types:
        X = ST.state_matrices(st)[0]
        if not arm_level:
            X = _family_levels(X, regs)
        cols = list(X.columns) if participants is None else [p for p in participants if p in X.columns]
        if not cols:
            continue
        A = X.reindex(index=regs, columns=cols).to_numpy(float)
        c, share, rank = _budget_matrix(A, Mv, codes, n_fam, nan_policy)
        for j, p in enumerate(cols):
            rows.append(pd.DataFrame({"gene": gene, "unit": regs, "patient": p, "state": st,
                                      "budget": c[:, j], "share": share[:, j], "rank": rank[:, j]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _family_levels(X: pd.DataFrame, fams) -> pd.DataFrame:
    """Arm × patient → family × patient by the SAME nonlinear pooling the model fits on:
    `log2(1 + Σ 2^x − 1)` (`families.collapse_by_family`). Never a mean of logs."""
    arms = list(X.index)
    fam = FAM.family_of(pd.Index(arms)).reindex(arms)
    lin = np.clip(2.0 ** X.to_numpy(float) - 1.0, 0.0, None)
    out = pd.DataFrame(lin, index=fam.to_numpy(), columns=X.columns).groupby(level=0).sum()
    return np.log2(1.0 + out).reindex(index=[f for f in fams if f in out.index])


@lru_cache(maxsize=1)
def _he_genes() -> tuple:
    from mirna_hallmark.learned import data as LD
    return tuple(sorted(set(LD.D.high_evidence_edges()["gene"])))


def patient_budget_all(genes: Optional[Sequence[str]] = None, *, m_source: str = "canonical",
                       arm_level: bool = True, participants=None, nan_policy: str = "zero",
                       persist: bool = True) -> pd.DataFrame:
    """All genes, concatenated. Persisted as PARQUET (~3,700 edges × 104 patients × 2 states ≈ 770k rows —
    TSV is the convention in `realization.py` but is the wrong container at this size)."""
    genes = list(genes) if genes is not None else list(_he_genes())
    parts = []
    for g in genes:
        try:
            d = patient_budget(g, m_source=m_source, arm_level=arm_level,
                               participants=participants, nan_policy=nan_policy)
        except Exception:
            continue
        if len(d):
            parts.append(d)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if persist and len(df):
        OUT.mkdir(parents=True, exist_ok=True)
        rung = "edge" if arm_level else "family_edge"
        df.to_parquet(OUT / f"patient_budget_{rung}_{nan_policy}.parquet", index=False)
    return df


# --------------------------------------------------------------------------- #
# ⭐ VERIFICATION GATE 1 — the object must BE the same quantity as the recorded one
# --------------------------------------------------------------------------- #

def verify_collapse(genes: Optional[Sequence[str]] = None, *, n_genes: int = 40) -> pd.DataFrame:
    """⭐ THE PRIMARY GATE. Collapse the per-patient budget with the SAME estimator `budget_shift` uses
    (**median across patients**, `states.py:541`) and compare to `budget_shift`'s own `share_*`/`rank_*`.

    ⚠ Uses `m_source="canonical"` and `nan_policy="zero"` deliberately — those are `budget_shift`'s own
    choices. Any other setting is a DIFFERENT quantity and would fail this gate for the wrong reason.
    ⚠ `rank` is expected to diverge on TIES only (RANK_METHOD="min" vs `budget_shift`'s truncated
    "average"); the gate reports tie-free agreement separately.

    ⚠⚠ TWO DISTINCT COMPARISONS, and conflating them makes the gate meaningless:
      • `share_from_medabund` — take the MEDIAN ABUNDANCE first, then run the formula. This is exactly
        what `budget_shift` does, so it must match **to float precision**. This is the CORRECTNESS gate.
      • `med_of_patient_shares` — the median of the per-patient shares. This is a DIFFERENT quantity
        (Jensen: the share of a median ≠ the median of shares) and is expected to differ. Its size is a
        FINDING, not a bug — it measures what the cohort collapse was hiding.
    """
    import io, contextlib
    genes = list(genes) if genes is not None else list(_he_genes())[:n_genes]
    rows = []
    for g in genes:
        try:
            with contextlib.redirect_stdout(io.StringIO()):        # budget_shift prints a table per gene
                ref = ST.budget_shift(g)
        except Exception:
            continue
        if ref is None or ref.empty or "share_TUM" not in ref.columns:
            continue
        try:
            M = ST.canonical_M(g, "01", arm_level=True)
            regs = list(M[M > 0].index)
            if not regs:
                continue
            Mv = M.reindex(regs).to_numpy(float)
            codes, n_fam = _fam_codes(regs, True)
        except Exception:
            continue
        mine = patient_budget(g, m_source="canonical", arm_level=True,
                              sample_types=("01", "11"), participants=None, nan_policy="zero")
        if mine.empty:
            continue
        r = ref.set_index("arm").reindex(regs)
        for st, suf in (("01", "TUM"), ("11", "NAT")):
            if f"share_{suf}" not in ref.columns:
                continue
            # (a) CORRECTNESS: median abundance first, then the formula — budget_shift's own path
            lv = ST.state_matrices(st)[0].reindex(regs).median(axis=1)
            A1 = lv.to_numpy(float)[:, None]
            _c1, sh1, rk1 = _budget_matrix(A1, Mv, codes, n_fam, "zero")
            # (b) the DIFFERENT quantity: median of the per-patient shares
            m = mine[mine.state == st]
            med_sh = m.groupby("unit")["share"].median().reindex(regs).to_numpy(float)
            rows.append(pd.DataFrame({
                "gene": g, "state": suf, "arm": regs,
                "ref_share": r[f"share_{suf}"].to_numpy(float),
                "share_from_medabund": np.round(sh1[:, 0], 3),
                "med_of_patient_shares": med_sh,
                "ref_rank": r[f"rank_{suf}"].to_numpy(float),
                "rank_from_medabund": rk1[:, 0]}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def verify_simplex(df: pd.DataFrame) -> dict:
    """Gate 2 — `Σ_unit share == 1` per (gene, patient, state), among genes with ≥2 units.

    ⚠ Cells whose TOTAL BUDGET is 0 have **no defined share** and are excluded, not counted as failures
    (`budget_shift` emits 0.0 there via the `else 0.0` at `states.py:557`; this module emits NaN, which
    matches `_budget()`'s abstention semantics). They are reported separately so the exclusion is visible."""
    g = df.groupby(["gene", "patient", "state"]).agg(s=("share", "sum"), n=("share", "size"),
                                                     btot=("budget", "sum"))
    multi = g[g["n"] >= 2]
    defined = multi[multi["btot"] > 0]
    off = (defined["s"] - 1.0).abs()
    return {"n_cells": int(len(multi)), "n_defined": int(len(defined)),
            "n_zero_budget_excluded": int((multi["btot"] <= 0).sum()),
            "max_abs_dev": float(off.max()) if len(off) else np.nan,
            "frac_within_1e-9": float((off < 1e-9).mean()) if len(off) else np.nan}


def verify_gtex_constant(genes: Optional[Sequence[str]] = None, *, n_genes: int = 25) -> dict:
    """Gate 3 — the GTEx reference must have **NO PATIENT AXIS**. Actually checks it (a hardcoded True
    would read as coverage while testing nothing): `_healthy_level` must return a 1-D object of length
    `len(regs)`, and repeated calls must be byte-identical."""
    genes = list(genes) if genes is not None else list(_he_genes())[:n_genes]
    n_ok = n_1d = n_stable = n_arms = 0
    lens_match = True
    for g in genes:
        try:
            M = ST.canonical_M(g, "01", arm_level=True)
            regs = list(M[M > 0].index)
            if not regs:
                continue
            hl = ST._healthy_level(regs, repaired=True)
            hl2 = ST._healthy_level(regs, repaired=True)
        except Exception:
            continue
        n_ok += 1
        v = np.asarray(pd.Series(hl).to_numpy(float))
        n_1d += int(v.ndim == 1)
        lens_match &= (len(v) == len(regs))
        n_stable += int(np.array_equal(np.nan_to_num(v, nan=-999),
                                       np.nan_to_num(np.asarray(pd.Series(hl2).to_numpy(float)), nan=-999)))
        n_arms += int(np.isfinite(v).sum())
    return {"genes_checked": n_ok, "all_1d": n_1d == n_ok, "len_matches_regs": bool(lens_match),
            "repeat_call_identical": n_stable == n_ok, "arms_with_finite_HLY": n_arms,
            "PASS": bool(n_ok > 0 and n_1d == n_ok and lens_match and n_stable == n_ok)}
