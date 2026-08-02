"""GENE AXES — the reusable multi-axis characterisation for ANY per-gene outcome.

    from mirna_hallmark.learned import gene_axes as GA
    A   = GA.build_axes(genes, cohort_rna=..., cohort_mirna_families=...)   # ~100 axes/gene
    S   = GA.scan(outcome, A)                                              # FDR modifier scan
    GA.contrast(outcome, A, k=80)                                          # top-vs-bottom profile
    GA.sign_analysis(reference, candidate)                                 # is the gain SIGN CORRECTION?

⭐ WHY THIS EXISTS (user-directed, 2026-08-02, generalised out of MH-201). The Buffa transfer test began as
one number — "does β transport?" — and the number alone (β −0.0186 vs abundance −0.0022) was nearly
uninterpretable. What made it a *finding* was asking the same question along many gene axes at once. Four
things came out that a single statistic could not have shown, and **three of the four came from axes I had
not thought to build**:

  1. **A degeneracy that governed the headline.** 43% of genes have ONE measurable regulator family, where
     `M = β·Z` and the reference is `Z`, so the two are IDENTICAL by construction (verified max|Δ| = 0.0).
     Pooling them made an abundance effect look like a β effect.
  2. **The strongest modifier was on the REGULATOR side, not the gene's** — `reg_dose_hhi` (q=2.1e-05).
  3. **Measurement profile beat expression level** — a gene's/regulator's DYNAMIC RANGE predicts, its MEAN
     does not (`buffa_sd` q=0.0016 vs `buffa_mean` q=0.45).
  4. **The gain was SIGN CORRECTION, asymmetrically applied** — where the reference had the wrong sign the
     candidate moved it (p=1.5e-16); where the reference was already right it moved it by EXACTLY 0.000.

⚠⚠ **TWO TRAPS THIS MODULE EXISTS TO STOP YOU REPEATING** (both cost real errors in MH-201):
  * **THE DEGENERACY TRAP.** Always ask where your candidate is *mathematically inert* relative to its
    reference and SPLIT on it. `mask_degenerate()` does this. A pooled statistic over a universe that is
    43% inert produced an incoherent "median Δ = +0.0000 with Wilcoxon p = 1.3e-27".
  * **THE MOVING-SUPPORT TRAP (axiom 5).** A concentration index like HHI is bounded below by 1/k, so
    `corr(k, HHI)` is strongly negative NO MATTER WHAT the biology does — I reported −0.667 as evidence
    that big designs are diffuse; **floor-corrected it is −0.075.** `hhi()` returns the normalised form by
    default for exactly this reason.

THE AXIS FAMILIES (all optional; build what the cohort supports)
---------------------------------------------------------------
  `card_*`    whatever the gene/family cards carry (design size, identifiability, composition, pressure)
  `self_*`    the TARGET's own expression profile — mean, sd, IQR ⇒ **sd/IQR are the informative ones**
  `reg_*`     ⭐ the REGULATOR ENSEMBLE — level, dynamic range, and how both DISTRIBUTE across regulators
              (`reg_dose_hhi`, `reg_var_hhi`, `reg_frac_flat`). Usually the strongest family. The
              mechanism: a weighted estimator can only beat an unweighted sum when the ensemble is SPREAD —
              if one member dominates the abundance, the sum already IS that member.
  `w_*`       the weight vector's own distribution (concentration, top share) — ⚠ in MH-201 this was INERT;
              what matters is how the DATA distributes, not how the weights do.
  `ident_*`   attribution/Shapley — ⚠ entered as **COVERAGE** (`ident_n_def`, q=4e-05) and NOT as shape
              (`ident_hhi` q=0.23, `ident_top1` q=0.44, agreement-with-β q=0.82).

⚠⚠ **THE LOADING-vs-DRIVER DISTINCTION — the subtlest trap here, and it looks like a refutation when it is
not.** A gene's own **expression loading** on a compartment (`spearman(gene, CAF_fraction)`) and that
compartment's **share of the gene's COUPLING** (`comp_*_driver_share`) are DIFFERENT QUANTITIES, and in
MH-201 **only the second predicted anything**. Driver-share separates the harm tails at **p=0.00999** and
gates real damage; **all 11 loading axes were null under FDR** (best q=0.18 — and CAFs themselves p=0.59,
with the sign running the *opposite* way). *"Is this gene expressed by stroma"* is not *"is this gene's
miRNA-target relationship an artifact of stroma"*. **Build the driver-share axis, not just the loading
axis — and never report a loading null as if it refuted a composition effect.**

⚠ **A STRUCTURALLY MISSING COVARIATE IS NOT AUTOMATICALLY THE CULPRIT.** Buffa cannot supply `target_cn`,
so CN-driven genes were the obvious suspects for the residual harm. **Null on every measure** — the
per-gene C-ablation cost vs the outcome p=0.72; `cn_var` p=0.87, `cn_amp_frac` p=0.30, `cn_absdev` p=0.54;
the harmed genes not CN-extreme (MWU p=0.998). Test the suspicion; do not assume the gap you know about is
the gap that hurts.

⛔ **SCAN THE CATEGORICAL CLASSES TOO — `scan()` IS NUMERIC-ONLY AND WILL SILENTLY SKIP THEM.** In MH-201
that dropped **6 class columns**, and one of them turned out to carry the mechanism: `ctx_apriori_class`
(Kruskal q=**0.0084**; every other categorical null at q≈0.97). Use `scan_categorical()`.

⭐⭐ **THE DEEPEST LESSON OF THE ARC, and it only appeared once the classes were tested: SOME AXES PREDICT
THE *MAGNITUDE OF ACTION* IN BOTH DIRECTIONS, OTHERS PREDICT ITS *SIGN*. DO NOT EXPECT THE FIRST KIND TO
GATE HARM.** `ctx_apriori_class = A_COMPETENT` (`n_fam≥3 ∧ w_max>median`) marks where β **acts at all**:
median Δ −0.0190 and genuine-harm **9.4%**, against **exactly +0.0000** and 3.0% for every other class —
and all 8 genuinely-harmed genes are A_COMPETENT. So the competence/ensemble axes select genes where β does
*more of everything*, help and harm alike; only the confounder/retention axes tell you which. That is why
gating on `reg_dose_hhi` moved the margin 3.7× and the harm rate not at all.
"""


from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr, wilcoxon

MIN_N_SCAN = 100


def hhi(v, *, normalise: bool = True) -> float:
    """Concentration of a non-negative vector.

    ⚠ **RAW HHI IS BOUNDED BELOW BY 1/k, SO IT FALLS AS k GROWS WHATEVER THE DATA DOES.** Correlating raw
    HHI against the item count therefore manufactures a strong negative correlation out of nothing — it
    cost me a wrong claim in MH-201 (−0.667 raw vs **−0.075** floor-corrected). `normalise=True` returns
    `(HHI − 1/k)/(1 − 1/k)` ∈ [0,1], which is comparable across k. Use the raw form only when k is fixed.
    """
    v = np.abs(np.asarray(v, float))
    v = v[np.isfinite(v)]
    k = len(v)
    if k < 2 or v.sum() <= 0:
        return np.nan
    h = float(((v / v.sum()) ** 2).sum())
    return float((h - 1 / k) / (1 - 1 / k)) if normalise else h


def regulator_axes(dose: np.ndarray) -> dict:
    """⭐ THE REGULATOR-ENSEMBLE AXES — usually the strongest family, and the one most often forgotten.

    `dose` is (n_regulators × n_samples) for ONE gene, on the cohort you are scoring in. Captures both the
    LEVEL and the DYNAMIC RANGE of each regulator, and — the part that matters — **how each DISTRIBUTES
    across the ensemble**.
    """
    if dose.ndim != 2 or dose.shape[0] < 2:
        return {}
    med = np.nanmedian(dose, axis=1)
    sd = np.nanstd(dose, axis=1)
    flat_cut = np.nanpercentile(sd, 10) + 1e-9
    return {"reg_n": float(dose.shape[0]),
            "reg_dose_med": float(np.nanmedian(med)), "reg_dose_max": float(np.nanmax(med)),
            "reg_dose_hhi": hhi(med),          # ⭐ one regulator dominating ABUNDANCE ⇒ a weighted
                                               #    estimator cannot beat the unweighted sum
            "reg_var_med": float(np.nanmedian(sd)), "reg_var_max": float(np.nanmax(sd)),
            "reg_var_min": float(np.nanmin(sd)), "reg_var_hhi": hhi(sd),
            "reg_frac_flat": float(np.mean(sd < flat_cut))}   # unmeasurable regulators = nothing to weight


def self_axes(y: np.ndarray) -> dict:
    """The TARGET's own measurement profile. ⚠ In MH-201 **sd and IQR predicted (q≈0.001–0.01) while the
    MEAN was nothing (q=0.45)** — a flat gene cannot correlate however highly expressed it is. Always carry
    a dispersion term; a level term alone will mislead you into thinking expression doesn't matter."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    if len(y) < 10:
        return {}
    return {"self_mean": float(y.mean()), "self_sd": float(y.std()),
            "self_iqr": float(np.percentile(y, 75) - np.percentile(y, 25))}


def weight_axes(w, identity=None) -> dict:
    """The weight vector's own shape, and the attribution/Shapley vector's.
    ⚠ Both were largely INERT in MH-201 — kept because a null on a plausible axis is worth recording, and
    because `ident_n_def` (COVERAGE) did predict while every SHAPE term did not."""
    out = {}
    w = np.abs(np.asarray(w, float))
    if np.isfinite(w).sum() >= 2:
        out |= {"w_hhi": hhi(w), "w_top1": float(np.nanmax(w) / np.nansum(w)) if np.nansum(w) > 0 else np.nan}
    if identity is not None:
        iv = np.asarray(identity, float)
        out |= {"ident_n_def": float(np.isfinite(iv).sum()),        # ⭐ the one that predicts
                "ident_frac_def": float(np.isfinite(iv).mean()),
                "ident_hhi": hhi(iv)}
        if np.isfinite(iv).sum() >= 3 and np.isfinite(w).sum() >= 3:
            out["ident_w_agree"] = float(spearmanr(np.abs(iv), w, nan_policy="omit").correlation)
    return out


def mask_degenerate(n_units: pd.Series, min_units: int = 2) -> pd.Series:
    """⚠⚠ **ALWAYS CALL THIS BEFORE COMPARING A WEIGHTED ESTIMATOR TO AN UNWEIGHTED ONE.**
    With a single unit, `Σ w·Z = w·Z` and the unweighted sum is `Z`; **Spearman is scale-invariant, so the
    two are IDENTICAL BY CONSTRUCTION** — the weights are mathematically inert, not merely weak. In MH-201
    that was **43% of the universe**, and pooling it produced a nonsense statistic and made an abundance
    effect read as a weight effect. Returns the KEEP mask."""
    return n_units >= min_units


def scan(outcome: pd.Series, axes: pd.DataFrame, *, min_n: int = MIN_N_SCAN,
         circular: set | None = None) -> pd.DataFrame:
    """Spearman every axis against one per-gene outcome, BH-FDR across axes.
    `circular` names axes that are the outcome's own source (flagged, not dropped — report them, don't lean
    on them)."""
    from mirna_hallmark.stats import bh_fdr
    circular = circular or set()
    rows = []
    for c in axes.columns:
        if not pd.api.types.is_numeric_dtype(axes[c]):
            continue
        m = pd.concat([outcome.rename("_y"), axes[c]], axis=1).dropna()
        if len(m) < min_n or m[c].nunique() < 5:
            continue
        r = spearmanr(m[c].to_numpy(float), m["_y"].to_numpy(float))
        rows.append({"axis": c, "n": len(m), "rho": float(r.correlation), "p": float(r.pvalue),
                     "circular": c in circular})
    S = pd.DataFrame(rows)
    if len(S):
        S["q"] = bh_fdr(S["p"].to_numpy())
        S = S.sort_values("p").reset_index(drop=True)
    return S


def contrast(outcome: pd.Series, axes: pd.DataFrame, k: int = 80) -> pd.DataFrame:
    """Profile the k most-helped vs k least-helped genes on every axis. The complement to `scan`: the scan
    says WHICH axes matter, this says WHAT the extremes actually look like — and it is where recognisable
    biology shows up (MH-201's top genes were MYC / CCND1 / HIF1A / FOXC1, its worst was FN1, a CAF/ECM
    gene from the composition-confounded class the scan had independently flagged)."""
    o = outcome.dropna().sort_values()
    top, bot = o.index[:k], o.index[-k:]
    rows = []
    for c in axes.columns:
        if not pd.api.types.is_numeric_dtype(axes[c]):
            continue
        a, b = axes.loc[top, c].dropna(), axes.loc[bot, c].dropna()
        if len(a) < 10 or len(b) < 10:
            continue
        rows.append({"axis": c, "low_outcome": float(a.median()), "high_outcome": float(b.median()),
                     "p": float(mannwhitneyu(a, b).pvalue)})
    return pd.DataFrame(rows).sort_values("p").reset_index(drop=True)


def sign_analysis(reference: pd.Series, candidate: pd.Series, *, expect_negative: bool = True) -> dict:
    """⭐ **IS THE GAIN A SIGN CORRECTION?** — the question that turned MH-201's modest margin into a claim.

    A weighted estimator can beat an unweighted one two ways: by sharpening an already-correct answer, or by
    **fixing answers the reference gets backwards**. They mean very different things. In MH-201 the
    unweighted reference had the correct sign in **51.1%** of genes — chance — while β had it in **58.5%**,
    rescuing 85 genes and breaking 36 (**net +49**); and the correction was **asymmetric**: where the
    reference was wrong β moved it (p=1.5e-16), where the reference was already right β moved it by
    **exactly 0.0000**. That is a corrective, not a rescaling, and it is invisible in any mean.
    """
    m = pd.concat([reference.rename("ref"), candidate.rename("cand")], axis=1).dropna()
    s = (lambda v: v < 0) if expect_negative else (lambda v: v > 0)
    ref_ok, cand_ok = s(m.ref), s(m.cand)
    tab = [[int((~ref_ok & cand_ok).sum()), int((~ref_ok & ~cand_ok).sum())],
           [int((ref_ok & cand_ok).sum()), int((ref_ok & ~cand_ok).sum())]]
    wrong, right = m[~ref_ok], m[ref_ok]
    out = {"n": len(m), "ref_sign_correct": float(ref_ok.mean()), "cand_sign_correct": float(cand_ok.mean()),
           "rescued": tab[0][0], "broken": tab[1][1], "net_rescue": tab[0][0] - tab[1][1],
           "fisher_p": float(fisher_exact(tab)[1]),
           "shift_where_ref_wrong": float((wrong.cand - wrong.ref).median()) if len(wrong) > 10 else np.nan,
           "shift_where_ref_right": float((right.cand - right.ref).median()) if len(right) > 10 else np.nan}
    if len(wrong) > 10:
        out["p_shift_where_ref_wrong"] = float(wilcoxon(wrong.cand - wrong.ref).pvalue)
    return out


def scan_categorical(outcome: pd.Series, cats: pd.DataFrame, *, min_level_n: int = 25) -> pd.DataFrame:
    """Kruskal-Wallis across the levels of each categorical class column, BH-FDR across columns.

    ⛔ **EXISTS BECAUSE `scan()` IS NUMERIC-ONLY AND SILENTLY SKIPS CLASS COLUMNS.** In MH-201 that dropped
    six of them, including `ctx_apriori_class` — the one that carried the mechanism (q=0.0084 while every
    other categorical was null at q≈0.97). A scan that quietly ignores a column type is worse than no scan,
    because its silence reads as coverage."""
    from mirna_hallmark.stats import bh_fdr
    from scipy.stats import kruskal
    rows = []
    for c in cats.columns:
        m = pd.concat([outcome.rename("_y"), cats[c]], axis=1).dropna()
        grps = [g["_y"].to_numpy(float) for _, g in m.groupby(c) if len(g) >= min_level_n]
        if len(grps) < 2:
            continue
        rows.append({"axis": c, "n": len(m), "levels": len(grps), "p": float(kruskal(*grps).pvalue)})
    S = pd.DataFrame(rows)
    if len(S):
        S["q"] = bh_fdr(S["p"].to_numpy())
        S = S.sort_values("p").reset_index(drop=True)
    return S
