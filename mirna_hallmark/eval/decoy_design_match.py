"""IS THE DECOY DESIGN MATCHED TO THE REAL ONE — in FAMILY STRUCTURE, not just in arm properties?

    .venv/bin/python3 -m mirna_hallmark.eval.decoy_design_match

WHY THIS EXISTS (2026-08-01, user-asked). `decoy_bench` matches ARMS — signed loadings + dose + variance,
Hungarian 1:1 — and only THEN collapses each side to seed families via `pool_family`. So two properties of
the DESIGN the estimator actually sees are never in the cost function:

  (A) WIDTH — the number of family columns. `_one` recorded only `Xr.shape[1]` (the REAL side), so the
      fake's width was not even in the artifact and this was unmeasurable until now (fixed: `n_fam_fake`).
      This matters because MH-147 measured that the GAP SCALES WITH DESIGN WIDTH
      (spearman(ceiling, n_fam)=+0.551; n_fam>=3 gap −0.0264 vs −0.0114 in the complement), so a systematic
      width asymmetry would contaminate the headline gap with a width effect.
  (B) WITHIN-DESIGN COLLINEARITY — the |r_S|<0.30 decorrelation filter constrains fake-vs-REAL correlation
      ONLY, never fake-vs-fake. Real regulators of one gene are often co-transcribed cluster mates
      (miR-17~92 is the canonical case), so the real design can be systematically more collinear, and
      collinearity moves the Gibbs's identifiability and hence rho directly.

⚠ THIS IS A CONTROL, NOT A FIX. If an asymmetry moved the gap, the answer would be to match on it — not to
  reweight the headline. The measured answer (below) is that neither does.

RESULT (2026-08-01, full universe, 1,395 genes with decoys / 1,349 scored)
-------------------------------------------------------------------------
BOTH ASYMMETRIES ARE REAL:
    width        real 3.096  fake 3.510   Δ +0.414   Wilcoxon p=4.6e-42
                 curation BUNDLES same-seed mates (1.087 arms/family real vs 1.004 fake); 81.8% of genes
                 are EQUAL-width — the mismatch is a 17.8% wide-design TAIL (+1.82 at real n_fam>=5).
    collinearity real 0.1663 fake 0.1547  Δ −0.0116  Wilcoxon p=6.7e-03
NEITHER MOVES THE GAP:
    b(Δn_fam)  = −0.00090 (p=0.722)  ⇒ contamination −0.0004 [−0.0025, +0.0017]
    b(Δcollin) = −0.02367 (p=0.618)  ⇒ contamination +0.0003 [−0.0008, +0.0014]
    b(Δdose)   = +0.01061 (p=0.020)  ⇒ the ONLY design axis with a detectable effect (already corrected
                 for post hoc — that is MH-139's b·Δ correction, gap@Δdose=0 = −0.0129).
⇒ the unmatched family structure contaminates the −0.0129 by <=~0.0025 rho at the CI edge (~19%),
  point estimate ~3%.

⚠ THE PRE-REGISTERED PREDICTION WAS WRONG — record it, do not quietly drop it (axiom 1a). Predicted: the
  asymmetries DEFLATE the gap (a wider, better-conditioned fake fits better), making −0.0129 conservative.
  Measured: UNDETECTED, not directional. **Do not carry a "conservative" framing from this control.** The
  defensible statement is a BOUND on contamination, nothing more.
"""
from __future__ import annotations

import os
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "mirna_hallmark/output/learned/decoy_design_match.tsv"
BENCH = ROOT / "mirna_hallmark/output/learned/decoy_bench.tsv"


def _ols(y: np.ndarray, X: np.ndarray) -> tuple:
    """b, se, p for an intercept-first design. Plain OLS — the regressors are gene-level and independent."""
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    n, k = X.shape
    sig2 = ((y - X @ b) ** 2).sum() / (n - k)
    se = np.sqrt(np.diag(sig2 * np.linalg.inv(X.T @ X)))
    return b, se, 2 * stats.t.sf(np.abs(b / se), n - k)


def geometry() -> pd.DataFrame:
    """Per-gene design geometry for BOTH sides. No Gibbs, no OOF — pure structure."""
    from mirna_hallmark.eval import decoy_bench as DB

    ctx = DB._ctx()
    P = DB.build_decoys()
    print(f"[design_match] decoys for {P.gene.nunique():,} genes / {len(P):,} pairs")
    idx = pd.Index(ctx["X"].columns)
    rows = []
    for g, sub in P.groupby("gene"):
        ra, fa = list(sub.real_arm), list(sub.fake_arm)
        Xr, Xf = DB.pool_family(ra, idx), DB.pool_family(fa, idx)
        rows.append({"gene": g, "n_arm": len(ra),
                     "n_fam_real": Xr.shape[1], "n_fam_fake": Xf.shape[1],
                     "collapse_real": len(ra) / Xr.shape[1], "collapse_fake": len(fa) / Xf.shape[1],
                     "collin_real": DB._collin(Xr), "collin_fake": DB._collin(Xf)})
    R = pd.DataFrame(rows)
    R["d_fam"] = R.n_fam_fake - R.n_fam_real
    R["d_collin"] = R.collin_fake - R.collin_real
    R.to_csv(OUT, sep="\t", index=False)
    return R


def report(R: pd.DataFrame) -> None:
    print(f"\n=== DECOY DESIGN-STRUCTURE MATCH — {len(R):,} genes ===\n")
    print("--- (A) DESIGN WIDTH ---")
    p = stats.wilcoxon(R.n_fam_fake, R.n_fam_real).pvalue
    print(f"  n_fam  real {R.n_fam_real.mean():.3f}  fake {R.n_fam_fake.mean():.3f}  "
          f"Δ {R.d_fam.mean():+.3f}  p={p:.2e}")
    print(f"  fake WIDER {(R.d_fam > 0).mean():.1%} · EQUAL {(R.d_fam == 0).mean():.1%} · "
          f"NARROWER {(R.d_fam < 0).mean():.1%}")
    print(f"  arms/family  real {R.collapse_real.mean():.3f}  fake {R.collapse_fake.mean():.3f}   "
          "(>1 = curation bundles same-seed mates)")
    print("\n--- (B) WITHIN-DESIGN COLLINEARITY ---")
    c = R.dropna(subset=["collin_real", "collin_fake"])
    pc = stats.wilcoxon(c.collin_fake, c.collin_real).pvalue
    print(f"  mean |rho| among family cols  real {c.collin_real.mean():.4f}  fake {c.collin_fake.mean():.4f}"
          f"  Δ {c.d_collin.mean():+.4f}  p={pc:.2e}  (n={len(c):,})")

    if not BENCH.exists():
        print("\n  ⚠ no decoy_bench.tsv — cannot test whether the asymmetries MOVE the gap")
        return
    D = pd.read_csv(BENCH, sep="\t").merge(R, on="gene", suffixes=("", "_g"))
    print(f"\n--- DO THEY MOVE THE GAP? (joined n={len(D):,}) ---")
    for lab, cols in [("d_fam", ["d_fam", "d_dose"]), ("d_collin", ["d_collin", "d_dose"])]:
        s = D.dropna(subset=["gap_core"] + cols)
        b, se, pv = _ols(s.gap_core.to_numpy(float), np.c_[np.ones(len(s)), s[cols].to_numpy(float)])
        i, mD = cols.index(lab) + 1, s[lab].mean()
        lo, hi = sorted(((b[i] - 1.96 * se[i]) * mD, (b[i] + 1.96 * se[i]) * mD))
        print(f"  b({lab:9s}) = {b[i]:+.5f}  p={pv[i]:.3g}   mean Δ {mD:+.4f}  "
              f"⇒ contamination {b[i] * mD:+.5f} [{lo:+.5f}, {hi:+.5f}]")
    s = D.dropna(subset=["gap_core", "d_dose"])
    b, se, pv = _ols(s.gap_core.to_numpy(float), np.c_[np.ones(len(s)), s.d_dose.to_numpy(float)])
    print(f"  b(d_dose   ) = {b[1]:+.5f}  p={pv[1]:.3g}   ⇒ gap@Δdose=0 = {b[0]:+.4f} "
          f"[{b[0] - 1.96 * se[0]:+.4f}, {b[0] + 1.96 * se[0]:+.4f}]   (the MH-139 correction)")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    report(geometry())
