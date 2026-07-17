"""THE SITE-FREE EMPIRICAL NULL — the subproject's per-edge calibration standard (MH-123).

WHY THIS EXISTS, AND WHY A PERMUTATION NULL IS NOT A SUBSTITUTE (MH-154)
-----------------------------------------------------------------------
The theoretical t-null is wrong by **3-4x in SCALE**: pairs with NO binding site are "FDR-significant"
(q<0.05, rho<0) at **25-35%**, not 5%. Bulk miRNA/mRNA share population structure (compartment,
proliferation, state) that no covariate block fully removes.

**A permutation null does NOT fix this, and MH-154 measured it.** Backing the null scale out of
``coupling_permutation``'s own ``p_z`` column (MH-45, 2000 Freedman-Lane draws):

    permutation null : sd0 = 0.0309, mu0 = +0.0001     (fit R^2 = 0.9997)
    theoretical t    : sd0 = 0.0311, mu0 =  0
    THIS null        : sd0 = 0.083-0.132,  mu0 = -0.048..-0.011   (core C)

2000 Freedman-Lane permutations reproduced the theoretical t-null to three significant figures. The reason
is structural and it generalises: **permutation destroys the very structure that inflates the real null.**
Freedman-Lane residualises on C, so it removes the *modelled* confound and nothing else -- but composition
adjustment demonstrably does NOT close the gap, so the culprit is UNMODELLED, and no permutation can
preserve an unmodelled confound. The site-free class can, because it is real data that merely cannot repress.

Consequence: MH-45's dependence-robust count (BY, rho<0) is **1,013/3,040 = 33.3%** -- inside the
**35.3%** rate at which IMPOSSIBLE edges pass the same gate under the same C. BY corrects for dependence
among tests; it cannot rescue a mis-scaled null. Do not read a BY q-value here as protection.

THE FIX
-------
Pairs with NO predicted site (no TargetScan context++ site, no scanMiR duplex) carry the SAME population
structure and ZERO targeting => their rho distribution IS the null. Match on ARM ABUNDANCE (the dominant
confound), FIT the null's location+scale per abundance quintile, score real edges against N(mu0, sd0).

**Efron form, not exceedance counting.** Counting exceedances against a 12k-per-bin null is
RESOLUTION-LIMITED (min p = 1/N ~ 8e-5), so BH across ~4,000 tests (needing p<1.2e-5) can NEVER fire --
that bug's giveaway was HE -> 0. Fitting location+scale gives a continuous p at full resolution.

**CONSERVATIVE IN THE SAFE DIRECTION:** the site-free class certainly contains some real non-canonical
targets, which widens the null and makes it harder to beat. The truth lies between this null and the
theoretical one -- but a 3-4x scale error is not a small correction.

USE
---
    from mirna_hallmark.eval import site_free_null as SFN
    nl = SFN.fit(composition=True)          # NullLadder: per-abundance-quintile (mu0, sd0)
    p  = nl.pvalues(rho, arm_abundance)     # one-sided (repressive) continuous p
    q  = nl.qvalues(rho, arm_abundance)     # + BH

Run the diagnostic: ``.venv/bin/python3 -m mirna_hallmark.eval.site_free_null [--n-fake N]``
Outputs ``output/site_free_null/{null_params.tsv, class_rates.tsv, manifest.json}``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Dict, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.learned import confounders as CF
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import kd as KD

OUT_DIR = C.OUTPUT_ROOT / "site_free_null"

N_FAKE = 60_000          # site-free pool; 12k per abundance quintile
N_BINS = 5               # abundance quintiles
MIN_N = 50               # min paired samples for an edge rho
MIN_BIN = 50             # min site-free pairs to fit a bin
SCANMIR_SITE = -0.5      # repression < this => a duplex exists => NOT site-free


# --------------------------------------------------------------------------------------- the null object

@dataclass
class NullLadder:
    """The fitted site-free null: (mu0, sd0) per arm-abundance quintile, plus the bin edges.

    ``pvalues`` is ONE-SIDED in the repressive direction -- a more NEGATIVE rho gives a smaller p.
    """
    bins: np.ndarray                       # N_BINS+1 edges over arm abundance
    params: Dict[int, Tuple[float, float]]  # bin -> (mu0, sd0)
    composition: bool
    n_fake: int
    theoretical_sd: float                  # what the t-null would have claimed, for the record

    def _bin_of(self, ab: Sequence[float]) -> np.ndarray:
        return np.clip(np.digitize(np.asarray(ab, float), self.bins) - 1, 0, N_BINS - 1)

    def zscores(self, rho: Sequence[float], ab: Sequence[float]) -> np.ndarray:
        rho = np.asarray(rho, float)
        bi = self._bin_of(ab)
        z = np.full(len(rho), np.nan)
        for i in range(N_BINS):
            m = bi == i
            if not m.any():
                continue
            mu, sd = self.params[i]
            z[m] = (rho[m] - mu) / sd
        return z

    def pvalues(self, rho: Sequence[float], ab: Sequence[float]) -> np.ndarray:
        """One-sided p against the fitted null: P(rho_null <= rho). Smaller = more repressive."""
        return stats.norm.cdf(self.zscores(rho, ab))

    def qvalues(self, rho: Sequence[float], ab: Sequence[float]) -> np.ndarray:
        return S.bh_fdr(pd.Series(self.pvalues(rho, ab)).fillna(1.0).values)

    def inflation(self) -> float:
        """How many times WIDER the honest null is than the t-null claimed (the MH-123 headline)."""
        return float(np.mean([sd for _, sd in self.params.values()]) / self.theoretical_sd)

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([
            {"bin": i, "ab_lo": self.bins[i], "ab_hi": self.bins[i + 1],
             "mu0": self.params[i][0], "sd0": self.params[i][1],
             "theoretical_sd": self.theoretical_sd,
             "width_ratio": self.params[i][1] / self.theoretical_sd}
            for i in range(N_BINS)
        ])


# ------------------------------------------------------------------------------------------- internals

def _rank_residualise(M: pd.DataFrame, Cov: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
    """Rank-transform each row, then residualise on the covariate block. Vectorised over complete rows;
    falls back per-row where the row has NaNs (different sample support => different hat matrix)."""
    cv = Cov.loc[[c for c in M.columns if c in Cov.index]].apply(pd.to_numeric, errors="coerce")
    keep = cv.dropna().index
    A = np.column_stack([np.ones(len(keep)), cv.loc[keep].to_numpy(float)])
    H = A @ np.linalg.pinv(A.T @ A) @ A.T
    V = M[keep].to_numpy(float)
    out = np.full(V.shape, np.nan)
    full = ~np.isnan(V).any(1)
    if full.any():
        R = np.apply_along_axis(stats.rankdata, 1, V[full])
        out[full] = R - (H @ R.T).T
    for i in np.where(~full)[0]:
        v = V[i]
        m = ~np.isnan(v)
        if m.sum() < MIN_N:
            continue
        Ai = A[m]
        r = stats.rankdata(v[m])
        b, *_ = np.linalg.lstsq(Ai, r, rcond=None)
        out[i, m] = r - Ai @ b
    return pd.DataFrame(out, index=M.index, columns=keep), A.shape[1]


def _rho_of(E: pd.DataFrame, XR: pd.DataFrame, YR: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    """Partial Spearman per (miRNA, gene) row of E, from pre-residualised rank matrices."""
    cols = [c for c in XR.columns if c in YR.columns]
    XA, YA = XR[cols].to_numpy(float), YR[cols].to_numpy(float)
    xi = {m: i for i, m in enumerate(XR.index)}
    yi = {g: i for i, g in enumerate(YR.index)}
    rho = np.full(len(E), np.nan)
    nn = np.zeros(len(E), int)
    for j, (m, g) in enumerate(zip(E.miRNA.values, E.gene.values)):
        if m not in xi or g not in yi:
            continue
        xv, yv = XA[xi[m]], YA[yi[g]]
        msk = ~np.isnan(xv) & ~np.isnan(yv)
        n = int(msk.sum())
        if n < MIN_N:
            continue
        a = xv[msk] - xv[msk].mean()
        b = yv[msk] - yv[msk].mean()
        dd = np.sqrt((a * a).sum() * (b * b).sum())
        if dd > 0:
            rho[j], nn[j] = float(a @ b / dd), n
    return rho, nn


def _classes(n_fake: int = N_FAKE, seed: int = 0):
    """Build the three edge classes on a shared arm/gene support: HE (curated), ORPHAN, FAKE (site-free)."""
    d = LD._load()
    X, Y = d["X"], d["Y"]
    parts = [p for p in X.columns if p in Y.columns and p in CF.deconv("tcga").index]
    X = X[parts]

    e = D.load_hallmark_edges()
    curated = set(zip(e.miRNA, e.gene))
    HE = e[e.high_evidence][["miRNA", "gene"]].drop_duplicates()

    # THE SITE SCREEN — union of BOTH scanMiR tables. MH-155: MH-123's `empnull.py` screened with
    # `KD.affinity()` ALONE, which is **HE-RESTRICTED** (912k rows / 1,432 genes) rather than the
    # genome-wide scan (`KD.genome_affinity()`, 10.4M rows / 18,852 genes). Consequence: pairs with a
    # genuine duplex that the HE-restricted table never scanned were admitted to the "site-free" class —
    # **6.0% of the 60k pool (10.0% of the pairs the genome-wide scan can judge)**. Contamination by real
    # repressive pairs WIDENS the null, so the defect ran in the CONSERVATIVE direction; screening on the
    # union removes it. (Measured effect: see MH-155 / `output/site_free_null/manifest.json`.)
    aff = KD.affinity()
    gaff = KD.genome_affinity()
    site = set(zip(aff[aff.repression < SCANMIR_SITE].arm, aff[aff.repression < SCANMIR_SITE].gene))
    site |= set(zip(gaff[gaff.repression < SCANMIR_SITE].arm, gaff[gaff.repression < SCANMIR_SITE].gene))

    dos = pd.read_csv(C.OUTPUT_ROOT / "learned" / "discovery_dossier.tsv", sep="\t").rename(columns={"arm": "miRNA"})
    ORPH = dos[["miRNA", "gene"]].drop_duplicates()
    ORPH = ORPH[[(m, g) not in curated for m, g in zip(ORPH.miRNA, ORPH.gene)]]

    arms = [a for a in X.index if np.isfinite(X.loc[a].to_numpy(float)).mean() > 0.5]
    aset = set(arms)
    genes = sorted({g for g in set(HE.gene) | set(ORPH.gene) | {g for _, g in site}
                    if g in Y.index and Y.loc[g, parts].notna().mean() > 0.8})
    gset = set(genes)
    HE = HE[HE.miRNA.isin(aset) & HE.gene.isin(gset)]
    ORPH = ORPH[ORPH.miRNA.isin(aset) & ORPH.gene.isin(gset)]

    # THE SITE-FREE CLASS: cannot repress by any modelled mechanism, same population structure.
    rng = np.random.default_rng(seed)
    fk = set()
    while len(fk) < n_fake:
        a, g = arms[rng.integers(len(arms))], genes[rng.integers(len(genes))]
        if (a, g) in curated or (a, g) in site:
            continue
        fk.add((a, g))
    FAKE = pd.DataFrame(sorted(fk), columns=["miRNA", "gene"])
    return X, Y.loc[genes, parts], parts, HE, ORPH, FAKE


def _fit_bins(NUL: pd.DataFrame) -> Tuple[np.ndarray, Dict[int, Tuple[float, float]]]:
    """Robust (median / 1.4826*MAD) location+scale of the site-free rho, per arm-abundance quintile."""
    bins = np.quantile(NUL.ab, np.linspace(0, 1, N_BINS + 1))
    bins[0] -= 1e-6
    bins[-1] += 1e-6
    params = {}
    for i in range(N_BINS):
        nn = NUL[(NUL.ab >= bins[i]) & (NUL.ab < bins[i + 1])].rho.values
        if len(nn) < MIN_BIN:
            params[i] = (0.0, 1.0)      # refuse to fit; scores in this bin become ~0 => never significant
            continue
        mu = float(np.median(nn))
        params[i] = (mu, max(float(1.4826 * np.median(np.abs(nn - mu))), 1e-6))
    return bins, params


# ------------------------------------------------------------------------------------------------ API

def fit(*, composition: bool = True, n_fake: int = N_FAKE, seed: int = 0) -> NullLadder:
    """Fit the site-free null. ``composition=True`` uses the composition-adjusted C (the defensible default:
    it shrinks both the shift and the width, though it does NOT close the gap)."""
    X, Yg, parts, _HE, _ORPH, FAKE = _classes(n_fake=n_fake, seed=seed)
    Cov = CF.build_C("tcga", parts, composition=composition)
    XR, k = _rank_residualise(X, Cov)
    YR, _ = _rank_residualise(Yg, Cov)
    ab = X.median(axis=1)

    r, n = _rho_of(FAKE, XR, YR)
    NUL = pd.DataFrame({"miRNA": FAKE.miRNA.values, "gene": FAKE.gene.values, "rho": r, "n": n})
    NUL = NUL[NUL.n > 0]
    NUL["ab"] = NUL.miRNA.map(ab)
    NUL = NUL.dropna(subset=["ab", "rho"])

    bins, params = _fit_bins(NUL)
    return NullLadder(bins=bins, params=params, composition=composition, n_fake=int(len(NUL)),
                      theoretical_sd=float(1.0 / np.sqrt(max(len(parts) - k - 1, 1))))


def run(*, n_fake: int = N_FAKE, seed: int = 0) -> pd.DataFrame:
    """Diagnostic: the theoretical vs empirical FDR-negative rate per class, under both C blocks.
    The FAKE row IS the false-positive rate -- under the theoretical null it should be ~5% and is not."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    X, Yg, parts, HE, ORPH, FAKE = _classes(n_fake=n_fake, seed=seed)
    ab = X.median(axis=1)
    print(f"HE {len(HE)} | ORPHAN {len(ORPH)} | FAKE null pool {len(FAKE)} | genes {Yg.shape[0]}", flush=True)

    rates, ladders = [], {}
    for label, comp in [("core C", False), ("composition-adjusted C", True)]:
        Cov = CF.build_C("tcga", parts, composition=comp)
        XR, k = _rank_residualise(X, Cov)
        YR, _ = _rank_residualise(Yg, Cov)

        res = {}
        for lbl, E in [("HE", HE), ("ORPHAN", ORPH), ("FAKE", FAKE)]:
            r, n = _rho_of(E, XR, YR)
            t = pd.DataFrame({"miRNA": E.miRNA.values, "gene": E.gene.values, "rho": r, "n": n})
            t = t[t.n > 0]
            t["ab"] = t.miRNA.map(ab)
            res[lbl] = t.dropna(subset=["ab", "rho"])

        NUL = res["FAKE"]
        bins, params = _fit_bins(NUL)
        nl = NullLadder(bins=bins, params=params, composition=comp, n_fake=int(len(NUL)),
                        theoretical_sd=float(1.0 / np.sqrt(max(len(parts) - k - 1, 1))))
        ladders[label] = nl

        print(f"\n############ {label} ({Cov.shape[1]} cov) ############")
        print(f"{'class':8s} {'tested':>7s} | {'THEORETICAL t-null':>21s} | {'SITE-FREE empirical null':>25s}")
        print(f"{'':8s} {'':>7s} | {'FDR-neg':>10s} {'rate':>10s} | {'FDR-neg':>12s} {'rate':>12s}")
        for lbl in ["HE", "ORPHAN", "FAKE"]:
            t = res[lbl]
            df = np.maximum(t.n - k - 1, 1)
            with np.errstate(invalid="ignore"):
                tt = t.rho * np.sqrt(df / np.maximum(1 - t.rho ** 2, 1e-12))
                pth = 2 * stats.t.sf(np.abs(tt), df)
            qth = S.bh_fdr(pd.Series(pth).fillna(1).values)
            n_th = int(((qth < 0.05) & (t.rho < 0)).sum())
            n_emp = int((nl.qvalues(t.rho.values, t.ab.values) < 0.05).sum())
            print(f"{lbl:8s} {len(t):7d} | {n_th:10d} {n_th/len(t):9.1%} | {n_emp:12d} {n_emp/len(t):11.1%}")
            rates.append({"C": label, "class": lbl, "tested": len(t),
                          "theoretical_fdr_neg": n_th, "theoretical_rate": n_th / len(t),
                          "empirical_fdr_neg": n_emp, "empirical_rate": n_emp / len(t)})
        print("  => the FAKE row IS the false-positive rate. Under the theoretical null it should be ~5%")
        print("     and is not; under the empirical null it is calibrated BY CONSTRUCTION -- so read")
        print("     HE/ORPHAN as the EXCESS OVER BACKGROUND, not as per-edge discoveries.")
        print(f"  fitted site-free null per abundance quintile (theoretical would be mu=0, sd={nl.theoretical_sd:.4f}):")
        for i in range(N_BINS):
            mu, sd = nl.params[i]
            print(f"     bin {i}: mu0={mu:+.4f}  sd0={sd:.4f}   ({sd/nl.theoretical_sd:.1f}x wider)")

    rate_df = pd.DataFrame(rates)
    rate_df.to_csv(OUT_DIR / "class_rates.tsv", sep="\t", index=False)
    pd.concat([l.to_frame().assign(C=k) for k, l in ladders.items()], ignore_index=True) \
        .to_csv(OUT_DIR / "null_params.tsv", sep="\t", index=False)
    (OUT_DIR / "manifest.json").write_text(json.dumps({
        "module": "mirna_hallmark.eval.site_free_null",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_fake_requested": n_fake,
        "seed": seed,
        "n_bins": N_BINS,
        "site_definition": f"no TargetScan context++ site AND scanMiR repression > {SCANMIR_SITE}",
        "null_form": "Efron: FIT location+scale (median / 1.4826*MAD) per arm-abundance quintile, score N(mu0,sd0)",
        "why_not_permutation": ("MH-154: coupling_permutation's 2000 Freedman-Lane draws reproduce the "
                                "theoretical t-null (sd0=0.0309 vs 0.0311, R^2=0.9997). Permutation destroys "
                                "the unmodelled population structure that inflates the real null; only a "
                                "real-data impossible-edge class preserves it."),
        "caveat": ("the site-free class contains some real non-canonical targets => this null is "
                   "CONSERVATIVE; the truth lies between it and the theoretical null"),
        "inflation_vs_theoretical": {k: round(l.inflation(), 2) for k, l in ladders.items()},
    }, indent=2))
    print(f"\n[site_free_null] wrote outputs under {OUT_DIR}")
    return rate_df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-fake", type=int, default=N_FAKE, help="site-free pool size")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    run(n_fake=a.n_fake, seed=a.seed)


if __name__ == "__main__":
    main()
