"""WHAT DOES THE LEARNED WEIGHTING BUY, PER GENE? — the β-weighted vs sum-abundance aggregate gap.

    .venv/bin/python3 -m mirna_hallmark.learned.weight_gain

MH-175 measured that the cross-cohort RANK transfer of the gene aggregate is very largely an ABUNDANCE
phenomenon (β-weighted +0.205 vs unweighted +0.189 in the independent cohort). This module asks the
per-gene version: **which genes does the weighting actually help, and why?**

    weight_gain(g) = rho_abund(g) − rho_agg(g)          POSITIVE = the learned weighting HELPS
                                                        (coupling is negative, so "more negative" = better)

⭐ THE INTERNAL NULL IS EXACT, NOT APPROXIMATE — and it is the reason this is worth measuring.
If a gene's β is CONSTANT across its arms, then `agg = Σ β·X = β·Σ X = β · abund`, a positive monotone
rescaling ⇒ **Spearman is IDENTICAL and `weight_gain` is EXACTLY 0**, not merely small. Since the card
broadcasts β equally within a seed family, this covers every ONE-FAMILY gene by construction. So:
  * a non-zero `weight_gain` in the constant-β stratum is a BUG, not a finding (a numerical check);
  * `beta_concentration` — how far β is from uniform — is the *mechanistic* predictor, and any
    characterisation that does not condition on it is describing design width in disguise.

PRE-REGISTERED PREDICTIONS (written before running — axiom 1a):
  P1  weight_gain == 0 (to numerical tolerance) for every constant-β gene. EXACT, an implementation check.
  P2  weight_gain increases with `beta_concentration` — this is the mechanism.
  P3  weight_gain increases with n_fam, but NOT after conditioning on beta_concentration (width would be
      acting only through "more families ⇒ more room for β to be uneven").
  P4  weight_gain is larger in TCGA (n≈1077) than in CPTAC (n≈100) — it is a signal, and CPTAC is noisy.

⚠ `weight_gain > 0` IS NOT EVIDENCE THE WEIGHTS ARE RIGHT. It is FITTED β against an UNFITTED sum;
MH-115/127 measured that a fitted matched DECOY also beats unfitted abundance, and decoy βs transport at
~80% of real. This module characterises WHERE the weighting acts, not whether it is correct — the
benchmark for that is a fitted fake (`eval/decoy_bench`).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "weight_gain_profile.tsv"
CONST_TOL = 1e-9        # |max β − min β| below this ⇒ β is constant ⇒ the exact internal null applies
HUB = ["PTEN", "ESR1", "ZEB1", "BCL2", "PDCD4", "CDKN1B", "GATA3"]


def _beta_profile() -> pd.DataFrame:
    """Per gene: how UNEVEN are the arm βs? Uniform β ⇒ the aggregate IS the abundance sum."""
    card = pd.read_csv(OUT / "edge_card_base.tsv", sep="\t",
                       usecols=["gene", "arm", "beta"], low_memory=False).dropna(subset=["beta"])
    rows = []
    for g, s in card.groupby("gene"):
        b = s.beta.to_numpy(float)
        b = b[np.isfinite(b)]
        if not len(b):
            continue
        n, tot = len(b), b.sum()
        p = b / tot if tot > 0 else np.full(n, 1.0 / n)
        H = float((p ** 2).sum())                       # Herfindahl; 1/n = uniform, 1 = one arm owns all
        conc = (H - 1.0 / n) / (1.0 - 1.0 / n) if n > 1 else 0.0   # 0 = uniform, 1 = fully concentrated
        rows.append({"gene": g, "n_arm_beta": n, "beta_conc": conc, "beta_eff_n": 1.0 / H,
                     "beta_const": bool(n == 1 or (b.max() - b.min()) < CONST_TOL)})
    return pd.DataFrame(rows)


def build() -> pd.DataFrame:
    gc = pd.read_csv(OUT / "gene_cptac_card.tsv", sep="\t")
    prof = _beta_profile()
    d = gc.merge(prof, on="gene", how="left")
    ctx = pd.read_csv(OUT / "gene_design_context.tsv", sep="\t")
    d = d.merge(ctx, on="gene", how="left", suffixes=("", "_ctx"))

    for pref, lay in (("tcga", "rna"), ("cptac_prosp", "rna"), ("cptac_prosp", "prot"),
                      ("cptac_t105", "rna"), ("cptac_t105", "prot")):
        a, b = f"{pref}_agg_rho_{lay}", f"{pref}_abund_rho_{lay}"
        if a in d and b in d:
            d[f"gain_{pref}_{lay}"] = d[b] - d[a]          # + = weighting helps
    d.to_csv(DEST, sep="\t", index=False)
    return d


def report(d: pd.DataFrame) -> None:
    print(f"=== WEIGHT GAIN (rho_abund − rho_agg; + = the learned β HELPS) — {len(d):,} genes ===\n")

    # ---------- P1: the EXACT internal null ----------
    g = "gain_tcga_rna"
    const = d[d.beta_const.fillna(False) & d[g].notna()]
    var = d[~d.beta_const.fillna(True) & d[g].notna()]
    print("--- P1  INTERNAL NULL: constant-β genes must give EXACTLY 0 (agg = β·abund, monotone) ---")
    print(f"  constant-β genes n={len(const):4d}   max |gain| = {const[g].abs().max():.3e}   "
          f"mean {const[g].mean():+.3e}   {'✅ EXACT' if const[g].abs().max() < 1e-9 else '⛔ NON-ZERO — BUG'}")
    print(f"  variable-β genes n={len(var):4d}   mean gain {var[g].mean():+.4f}  "
          f"median {var[g].median():+.4f}  helps in {(var[g] > 0).mean():.0%}")

    # ---------- P2/P3: mechanism ----------
    print("\n--- P2/P3  is the gain driven by β NON-UNIFORMITY, or just by design width? ---")
    v = var.dropna(subset=["beta_conc", "n_fam"])
    for x in ("beta_conc", "n_fam", "beta_eff_n", "n_arm_beta"):
        if x in v:
            s = stats.spearmanr(v[x], v[g])
            print(f"  spearman(gain, {x:12s}) = {s.correlation:+.4f}  p={s.pvalue:.2g}")
    # partial: does width survive conditioning on concentration, and vice versa?
    def _partial(x, z):
        s = v.dropna(subset=[x, z, g])
        rx = s[x] - np.poly1d(np.polyfit(s[z], s[x], 1))(s[z])
        ry = s[g] - np.poly1d(np.polyfit(s[z], s[g], 1))(s[z])
        r = stats.spearmanr(rx, ry)
        return r.correlation, r.pvalue, len(s)
    for x, z in (("beta_conc", "n_fam"), ("n_fam", "beta_conc")):
        r, p, n = _partial(x, z)
        print(f"  PARTIAL  gain ~ {x:10s} | {z:10s} = {r:+.4f}  p={p:.2g}  n={n}")

    # ---------- strata ----------
    print("\n--- gain by β-concentration quartile (variable-β genes only) ---")
    q = v.beta_conc.quantile([.25, .5, .75]).to_numpy()
    edges = [-np.inf, *q, np.inf]
    for i in range(4):
        s = v[(v.beta_conc > edges[i]) & (v.beta_conc <= edges[i + 1])]
        if len(s) < 10:
            continue
        print(f"  Q{i+1} (conc<= {edges[i+1]:.3f})  n={len(s):4d}  mean gain {s[g].mean():+.4f}  "
              f"helps {(s[g] > 0).mean():.0%}  mean n_fam {s.n_fam.mean():.1f}")

    print("\n--- gain by design width (variable-β genes only) ---")
    for lab, sel in [("n_fam 2", v.n_fam == 2), ("n_fam 3-4", v.n_fam.between(3, 4)),
                     ("n_fam 5+", v.n_fam >= 5)]:
        s = v[sel]
        if len(s) < 10:
            continue
        print(f"  {lab:10s} n={len(s):4d}  mean gain {s[g].mean():+.4f}  helps {(s[g] > 0).mean():.0%}  "
              f"mean conc {s.beta_conc.mean():.3f}")

    # ---------- P4: cohort ----------
    print("\n--- P4  cohort comparison (variable-β genes) ---")
    for c in ("gain_tcga_rna", "gain_cptac_prosp_rna", "gain_cptac_prosp_prot",
              "gain_cptac_t105_rna", "gain_cptac_t105_prot"):
        if c not in d:
            continue
        s = d[~d.beta_const.fillna(True)][c].dropna()
        if len(s) < 20:
            continue
        t = stats.wilcoxon(s, alternative="greater") if len(s) > 10 else None
        print(f"  {c:24s} n={len(s):4d}  mean {s.mean():+.4f}  median {s.median():+.4f}  "
              f"helps {(s > 0).mean():.0%}  p={t.pvalue:.2g}")


def hub_gibbs(d: pd.DataFrame) -> None:
    """⭐ MH-171 RE-RUN ON THE CANONICAL ESTIMATOR. Its retention used `ood_protein.fit_M` =
    `regression.fit_gene` = the ADAPTIVE LASSO, which MH-115 RETIRED. These columns use the card's
    dense Gibbs β instead — same genes, same cohorts, same composition block, canonical weights."""
    print("\n\n=== MH-171 HUB GENES, RE-RUN WITH THE GIBBS β (was: retired adaptive lasso) ===")
    print("    aggregate Σβ·X vs CPTAC protein · composition-adjusted / RAW / retention\n")
    lasso = {("PTEN", "prospective"): 0.461, ("ESR1", "prospective"): 0.797, ("ZEB1", "prospective"): 0.246,
             ("BCL2", "prospective"): 1.072, ("PDCD4", "prospective"): 0.659,
             ("CDKN1B", "prospective"): 1.000, ("GATA3", "prospective"): 0.802,
             ("PTEN", "tcga105"): 0.656, ("ESR1", "tcga105"): 0.886, ("ZEB1", "tcga105"): 0.308,
             ("BCL2", "tcga105"): 0.802, ("PDCD4", "tcga105"): 0.648,
             ("CDKN1B", "tcga105"): 0.938, ("GATA3", "tcga105"): 0.395}
    for cohort, pref in (("prospective", "cptac_prosp"), ("tcga105", "cptac_t105")):
        print(f"[{cohort}]")
        print(f"  {'gene':8s} {'rho_RAW':>9s} {'rho_adj':>9s} {'retention':>10s} {'class':>22s} "
              f"{'MH-171(lasso)':>14s} {'Δ':>8s}")
        rows = []
        for gname in HUB:
            r = d[d.gene == gname]
            if not len(r):
                print(f"  {gname:8s} {'—':>9s}  (not scored)")
                continue
            r = r.iloc[0]
            adj, raw = r.get(f"{pref}_agg_rho_prot"), r.get(f"{pref}_agg_rho_prot_raw")
            ret = r.get(f"{pref}_agg_ret_prot")
            cls = ("cell_intrinsic" if ret >= 0.7 else "composition_explained" if ret < 0.4 else "partial") \
                if pd.notna(ret) else "ungated"
            old = lasso.get((gname, cohort), np.nan)
            dlt = ret - old if pd.notna(ret) and pd.notna(old) else np.nan
            print(f"  {gname:8s} {raw:+9.3f} {adj:+9.3f} {ret:10.3f} {cls:>22s} "
                  f"{old:14.3f} {dlt:+8.3f}" if pd.notna(ret) else
                  f"  {gname:8s} {raw:+9.3f} {adj:+9.3f} {'n/a':>10s} {cls:>22s} {old:14.3f}")
            if pd.notna(ret):
                rows.append(ret)
        if rows:
            print(f"  → median retention (GIBBS) {np.median(rows):.3f}   "
                  f"[MH-171 lasso: {0.797 if cohort == 'prospective' else 0.656:.3f}]")


if __name__ == "__main__":
    df = build()
    report(df)
    hub_gibbs(df)
    print(f"\n-> {DEST}")
