"""THE DECOY GAP, STRATIFIED BY THE GENE'S REGULATOR-DESIGN COMPOSITION (user-asked, 2026-08-01).

    .venv/bin/python3 -m mirna_hallmark.eval.decoy_by_gene_strata

The gap is a PER-GENE quantity, so it can be read on the design each gene actually has, rather than only as
one pooled headline. The question this answers: **which genes is the real-vs-fake advantage coming from?**

THE DESIGN DESCRIPTORS (computed from the gene's REAL curated regulators — decoy-independent, a priori):
    n_arm            how many curated regulator arms
    n_fam            how many SEED FAMILIES they collapse to  (the design width the estimator sees)
    collapse         n_arm / n_fam  — >1 means curation bundles same-seed mates onto this gene
    n_fam_multi      families contributing >=2 arms
    n_abund          real arms above the canonical floor `log2(11)` (RPM>=10) — `_FLOOR`, the same cut
                     `eval/admissibility_bench` and the cn_guide_passenger panel use
    n_fam_abund      families containing at least one abundant arm
    frac_abund       n_abund / n_arm
    dose_max/med     max / median arm dose (log2 RPM+1) among the real regulators

WHY ABUNDANCE IS A DESIGN AXIS AND NOT A NUISANCE HERE (MH-166, user-corrected): abundant arms being the
dominant realized regulators IS the biology, not a tilt to be removed. This module therefore reads the gap
ACROSS abundance strata rather than adjusting abundance away.
⚠ Do not confuse this with `d_dose` = fake−real dose MISMATCH, which is a matching defect and IS corrected
  (MH-139's b·Δ). Here the strata are properties of the REAL side.

THE STANDING CAVEATS — wired in, not left to memory:
  * `n_fam == 1` is the INTERNAL NULL: beta is uniform there by construction, so the gap MUST be ~0. A
    non-zero value in that row is a design artifact, never a finding.
  * `apriori_class` is reported for continuity ONLY. ⚠ MH-147 RETIRED the rule as stated: the evidence axis
    (`w_max > median`) adds NOTHING given width (MWU p=0.293 within n_fam>=3) and the complement is NOT zero
    (−0.0114, p=6.3e-03). **The rule reduces to `n_fam >= 3`.** Do not re-headline the "27% domain".
  * CEILING is SWEPT, never reported at one threshold — MH-144: `ceiling <= 0` is a threshold count sitting
    exactly where the mass piles up (39.8% of genes within +/-0.01 of it); the robust point is `<= 0.02`.
  * `retention = gap_deconv/gap_core` is a RATIO and is GATED at |gap_core| >= RETENTION_GATE (axiom 5);
    below that it prints `n/a` rather than a large meaningless number.
  * Continuous descriptors are cut on QUARTILES of their own distribution as well as on the canonical floor,
    so no conclusion rests on one arbitrary cut.
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
BENCH = ROOT / "mirna_hallmark/output/learned/decoy_bench.tsv"
MAP = ROOT / "mirna_hallmark/output/learned/gene_competence_map.tsv"
OUT = ROOT / "mirna_hallmark/output/learned/decoy_by_gene_strata.tsv"
PROFILE = ROOT / "mirna_hallmark/output/learned/decoy_design_profile.tsv"

_FLOOR = np.log2(11)       # RPM>=10 — the canonical arm-expression floor used across the subproject
RETENTION_GATE = 0.005     # |gap_core| below this -> retention is a coin-flip with decimals (axiom 5)
MIN_N = 25                 # below this, no dose-adjusted intercept is fitted


def design_profile() -> pd.DataFrame:
    """Per-gene REAL-side design composition. No decoy, no Gibbs — a priori structure only."""
    from mirna_hallmark.eval import decoy_bench as DB

    ctx = DB._ctx()
    dose, fam = ctx["dose"], ctx["fam"]
    rows = []
    for g in sorted(ctx["he"].gene.unique()):
        arms = DB.real_regulators(g)
        if not arms:
            continue
        d = dose.reindex(arms).to_numpy(float)
        f = pd.Series([fam.get(a) for a in arms], index=arms)
        fc = f.value_counts()
        ab = [a for a, dv in zip(arms, d) if np.isfinite(dv) and dv >= _FLOOR]
        rows.append({
            "gene": g, "n_arm": len(arms), "n_fam": int(f.nunique(dropna=False)),
            "n_fam_multi": int((fc >= 2).sum()),
            "n_abund": len(ab), "frac_abund": len(ab) / len(arms),
            "n_fam_abund": int(pd.Series([fam.get(a) for a in ab]).nunique()) if ab else 0,
            "dose_max": np.nanmax(d) if np.isfinite(d).any() else np.nan,
            "dose_med": np.nanmedian(d) if np.isfinite(d).any() else np.nan,
        })
    P = pd.DataFrame(rows)
    P["collapse"] = P.n_arm / P.n_fam
    P.to_csv(PROFILE, sep="\t", index=False)
    print(f"[strata] design profile: {len(P):,} genes -> {PROFILE.name}")
    return P


def _adj(s: pd.DataFrame, col: str) -> tuple:
    """Gap at Δdose=0, b RE-DERIVED inside this stratum (never imported — MH-136b's retracted error)."""
    d = s.dropna(subset=[col, "d_dose"])
    if len(d) < MIN_N:
        return np.nan, np.nan
    X = np.c_[np.ones(len(d)), d.d_dose.to_numpy(float)]
    y = d[col].to_numpy(float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    sig2 = ((y - X @ b) ** 2).sum() / (len(d) - 2)
    se = np.sqrt(np.diag(sig2 * np.linalg.inv(X.T @ X)))
    return b[0], 1.96 * se[0]


def _row(lab: str, s: pd.DataFrame) -> dict:
    o = {"stratum": lab, "n": len(s)}
    if not len(s):
        return o
    o["gap_core"] = s.gap_core.mean()
    o["p_core"] = (stats.wilcoxon(s.real_core, s.dec_core, alternative="less").pvalue
                   if len(s) >= 10 else np.nan)
    o["win_pct"] = (s.gap_core < 0).mean()
    o["gap_core_adj"], o["adj_ci"] = _adj(s, "gap_core")
    if "gap_deconv" in s:
        d = s.dropna(subset=["gap_deconv"])
        o["gap_deconv"] = d.gap_deconv.mean() if len(d) else np.nan
        o["retention"] = (o["gap_deconv"] / o["gap_core"]
                          if abs(o["gap_core"]) >= RETENTION_GATE else np.nan)
    for c in ("d_dose", "d_fam", "real_core", "dec_core"):
        if c in s:
            o[c] = s[c].mean()
    return o


def _qcut(G: pd.DataFrame, col: str, rows: list, label: str) -> None:
    """Quartile strata on a continuous descriptor — avoids resting a conclusion on one arbitrary cut."""
    v = G[col].dropna()
    if v.nunique() < 4:
        return
    qs = v.quantile([0.25, 0.5, 0.75]).to_numpy()
    rows.append({"stratum": f"--- {label} (quartiles: {qs[0]:.2f} / {qs[1]:.2f} / {qs[2]:.2f}) ---"})
    edges = [-np.inf, *qs, np.inf]
    for i in range(4):
        sel = (G[col] > edges[i]) & (G[col] <= edges[i + 1])
        rows.append(_row(f"{col} Q{i+1}", G[sel]))


def main() -> pd.DataFrame:
    G = pd.read_csv(BENCH, sep="\t")
    if "n_fam_fake" in G:
        G["d_fam"] = G.n_fam_fake - G.n_fam
    P = design_profile()
    G = G.merge(P.drop(columns=[c for c in ("n_arm", "n_fam") if c in G]), on="gene", how="left")
    print(f"[strata] joined design profile: {G.collapse.notna().sum():,}/{len(G):,} genes profiled")
    if MAP.exists():
        M = pd.read_csv(MAP, sep="\t")
        keep = [c for c in ("gene", "apriori_class", "w_max", "ceiling_R2_oof_fam_core")
                if c in M.columns]
        G = G.merge(M[keep], on="gene", how="left", suffixes=("", "_map"))
        n = G.apriori_class.notna().sum() if "apriori_class" in G else 0
        print(f"[strata] joined competence map: {n:,}/{len(G):,} genes classed")
    else:
        print("  ⚠ no gene_competence_map.tsv — run `learned.analyses.gene_atlas` first")

    rows = [_row("ALL", G)]

    rows.append({"stratum": "--- DESIGN WIDTH: n_fam  (n_fam==1 is the INTERNAL NULL) ---"})
    for lab, sel in [("n_fam 1  [null]", G.n_fam == 1), ("n_fam 2", G.n_fam == 2),
                     ("n_fam 3-4", G.n_fam.between(3, 4)), ("n_fam 5+", G.n_fam >= 5)]:
        rows.append(_row(lab, G[sel]))
    rows.append(_row("n_fam >=3  [the MH-147 rule]", G[G.n_fam >= 3]))
    rows.append(_row("n_fam <3   [complement — NOT zero]", G[G.n_fam < 3]))

    rows.append({"stratum": "--- REGULATOR COUNT: n_arm ---"})
    for lab, sel in [("n_arm 1", G.n_arm == 1), ("n_arm 2", G.n_arm == 2),
                     ("n_arm 3-4", G.n_arm.between(3, 4)), ("n_arm 5+", G.n_arm >= 5)]:
        rows.append(_row(lab, G[sel]))

    rows.append({"stratum": "--- FAMILY BUNDLING: collapse = n_arm/n_fam ---"})
    rows.append(_row("collapse == 1  [no same-seed mates]", G[G.collapse == 1]))
    rows.append(_row("collapse > 1   [>=1 multi-arm family]", G[G.collapse > 1]))
    rows.append(_row("n_fam_multi 0", G[G.n_fam_multi == 0]))
    rows.append(_row("n_fam_multi 1", G[G.n_fam_multi == 1]))
    rows.append(_row("n_fam_multi 2+", G[G.n_fam_multi >= 2]))

    rows.append({"stratum": f"--- ABUNDANT REGULATORS (arm dose >= log2(11) = {_FLOOR:.2f}, RPM>=10) ---"})
    for lab, sel in [("n_abund 0  [all sparse]", G.n_abund == 0), ("n_abund 1", G.n_abund == 1),
                     ("n_abund 2", G.n_abund == 2), ("n_abund 3-4", G.n_abund.between(3, 4)),
                     ("n_abund 5+", G.n_abund >= 5)]:
        rows.append(_row(lab, G[sel]))
    rows.append({"stratum": "--- FAMILIES CARRYING AN ABUNDANT ARM ---"})
    for lab, sel in [("n_fam_abund 0", G.n_fam_abund == 0), ("n_fam_abund 1", G.n_fam_abund == 1),
                     ("n_fam_abund 2", G.n_fam_abund == 2), ("n_fam_abund 3+", G.n_fam_abund >= 3)]:
        rows.append(_row(lab, G[sel]))
    _qcut(G, "frac_abund", rows, "ABUNDANT FRACTION")
    _qcut(G, "dose_max", rows, "TOP REGULATOR DOSE")
    _qcut(G, "dose_med", rows, "MEDIAN REGULATOR DOSE")

    # ⭐ THE CROSS the user asked for: width x abundance. Is the width effect just abundance, or separable?
    rows.append({"stratum": "--- CROSS: WIDTH x ABUNDANCE (is width just abundance in disguise?) ---"})
    for wlab, wsel in [("n_fam<3", G.n_fam < 3), ("n_fam>=3", G.n_fam >= 3)]:
        for alab, asel in [("n_abund=0", G.n_abund == 0), ("n_abund 1-2", G.n_abund.between(1, 2)),
                           ("n_abund 3+", G.n_abund >= 3)]:
            rows.append(_row(f"{wlab} & {alab}", G[wsel & asel]))

    if "apriori_class" in G and G.apriori_class.notna().any():
        rows.append({"stratum": "--- APRIORI CLASS (⚠ MH-147: reduces to n_fam>=3) ---"})
        for k in ("A_COMPETENT", "B_PARTIAL", "C_WEAK", "D_UNDEFINED"):
            s = G[G.apriori_class == k]
            if len(s):
                rows.append(_row(k, s))

    cc = "ceiling_R2_oof_fam_core"
    if cc in G and G[cc].notna().any():
        rows.append({"stratum": "--- CEILING SWEEP (⚠ MH-144: <=0 is FRAGILE — mass piles up there) ---"})
        near = ((G[cc] > -0.01) & (G[cc] < 0.01)).mean()
        rows.append({"stratum": f"    [{near:.1%} of genes sit within +/-0.01 of zero ceiling]"})
        for t in (0.0, 0.01, 0.02, 0.05):
            s = G[G[cc] <= t]
            rows.append(_row(f"ceiling <= {t:.2f} [not measurable] ({len(s)/len(G):.0%})", s))
            rows.append(_row(f"ceiling >  {t:.2f} [measurable]", G[G[cc] > t]))

    R = pd.DataFrame(rows)
    R.to_csv(OUT, sep="\t", index=False)

    print(f"\n=== THE DECOY GAP BY GENE-DESIGN STRATUM — {len(G):,} genes ===")
    print("    [gap = raw; adj = Δdose residual removed WITHIN the stratum, b re-derived there]")
    print(f"    [retention gated at |gap_core| >= {RETENTION_GATE} (axiom 5); 'n/a' = denominator too small]\n")
    print(f"  {'stratum':38s} {'n':>5s} {'gap':>9s} {'adj':>9s} {'p':>9s} {'win':>5s} "
          f"{'deconv':>9s} {'ret':>6s} {'Δdose':>7s}")
    for _, r in R.iterrows():
        if pd.isna(r.get("n")):
            print(f"\n  {r.stratum}")
            continue
        def f(k, w, d="{:+.4f}"):
            v = r.get(k)
            return (" " * (w - 3) + "n/a") if v is None or pd.isna(v) else d.format(v).rjust(w)
        print(f"  {r.stratum:38s} {int(r.n):5d} {f('gap_core',9)} {f('gap_core_adj',9)} "
              f"{f('p_core',9,'{:.2e}')} {f('win_pct',5,'{:.0%}')} {f('gap_deconv',9)} "
              f"{f('retention',6,'{:.2f}')} {f('d_dose',7,'{:+.3f}')}")
    print(f"\n-> {OUT}")
    return R


if __name__ == "__main__":
    main()
