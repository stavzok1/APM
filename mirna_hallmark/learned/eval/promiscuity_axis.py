"""Does REGULATOR PROMISCUITY predict where the learned β beats an unweighted abundance sum? (MH-208)

This is the aggregate-force design's one surviving idea, re-posed in the learned frame. That design
(`method_dev/aggregate_pressure/AGGREGATE_FORCE_VS_ABUNDANCE_DESIGN.md`, DROPPED 2026-08-03) proposed
`W(m,g) = w_eff(m,g) / D(m)` — discount a regulator by its **promiscuity**, the miRNA's output capacity
split across its whole targetome. Its estimator (`w_eff`) is the MH-115-retired heuristic and its
comparator was a bare abundance baseline, so the design itself is dead; but `D(m)` is a regulator
property the ensemble axes cannot see, because it is about the arm's behaviour **outside this gene**.

⛔ **IT MUST BE THE SEQUENCE TARGETOME, NOT THE CURATED ONE.** The pre-existing `he_degree*` annotation
is a FAME axis (ρ=+0.736 with an arm's distinct-PMID count; +0.556 with abundance; only +0.124 with the
sequence targetome, and the two top-10 lists share nothing). See `analyses/misc/genomewide_promiscuity`.
Both are scored here, precisely so the fame channel is visible rather than assumed away.

**PROVENANCE NOTE (why this module exists at all).** MH-201's `ood_cohort_regulator_features.tsv` and
`ood_cohort_modifier_scan.tsv` — which carry its headline `reg_dose_hhi` (q=2.1e-05) — have **NO
PRODUCER** in the repo (verified 2026-08-03: `ood_cohort.py` writes `ood_cohort_genes.tsv`, the
`edge_leg_*` tables and the manifest, and nothing anywhere writes those two). That is the MH-196 shape
that rotted the literature sets. This module therefore rebuilds its outcome from PRODUCED artifacts
only: `delta = rho_buffa_metagene − rho_abund_metagene` from `ood_cohort_genes.tsv` (verified to
reproduce the recorded `delta` at max|diff| = 1.0e-16) and the design from `gene_family_card.tsv`.

Run: `.venv/bin/python3 -m mirna_hallmark.learned.eval.promiscuity_axis`
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "ood_cohort"
FAMILY_CARD = OUT / "gene_family_card.tsv"
GENES = DEST / "ood_cohort_genes.tsv"
PROMISC_AXES = ("reg_promisc_med", "reg_promisc_max", "reg_promisc_min",
                "reg_promisc_sd", "reg_promisc_hhi")


def gene_promiscuity(col: str = "promisc_seq_strong") -> pd.DataFrame:
    """Per-gene regulator-promiscuity axes, over the arms of that gene's design families."""
    from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
    from mirna_hallmark.learned.gene_axes import hhi

    P = GP.load_promiscuity(col=col, fill="nan")
    fam = pd.read_csv(FAMILY_CARD, sep="\t", usecols=["gene", "family", "arms"])
    rows = []
    for g, sub in fam.groupby("gene"):
        arms = sorted({a for s in sub["arms"].dropna() for a in str(s).split(";") if a})
        if not arms:
            continue
        p = np.array([P.get(a, np.nan) for a in arms], float)
        ok = p[np.isfinite(p)]
        r = {"gene": g, "n_arms_design": len(arms),
             "reg_promisc_cov": float(np.mean(np.isfinite(p)))}
        if len(ok) >= 2:
            r.update({"reg_promisc_med": float(np.median(ok)), "reg_promisc_max": float(ok.max()),
                      "reg_promisc_min": float(ok.min()), "reg_promisc_sd": float(ok.std()),
                      "reg_promisc_hhi": hhi(np.expm1(ok))})
        rows.append(r)
    return pd.DataFrame(rows)


def outcome() -> pd.DataFrame:
    """β-vs-abundance margin per gene, rebuilt from the PRODUCED ood_cohort table (see provenance note)."""
    g = pd.read_csv(GENES, sep="\t")
    g["delta"] = g["rho_buffa_metagene"] - g["rho_abund_metagene"]
    return g[["gene", "delta", "rho_buffa_metagene", "rho_abund_metagene",
              "n_fam_buffa", "fam_cov"]]


def run(min_cov: float = 0.5) -> pd.DataFrame:
    seq = gene_promiscuity("promisc_seq_strong")
    cur = gene_promiscuity("promisc_he_expr").rename(
        columns={c: c.replace("reg_promisc", "cur_promisc") for c in
                 [*PROMISC_AXES, "reg_promisc_cov"]})
    M = outcome().merge(seq, on="gene", how="inner").merge(
        cur[["gene", *[c.replace("reg_promisc", "cur_promisc") for c in PROMISC_AXES]]],
        on="gene", how="left")
    print(f"[promisc] {len(M):,} genes with an outcome and a design")
    print(f"[promisc] K_D coverage of design arms: median {M.reg_promisc_cov.median():.2f}; "
          f"{(M.reg_promisc_cov >= min_cov).sum():,} genes at >= {min_cov}")

    # ⛔ axiom 8: where β is MATHEMATICALLY INERT the margin cannot be about weighting. Split, never pool.
    M["degenerate"] = M["n_fam_buffa"].fillna(0) <= 1
    print(f"[promisc] single-family (β inert) genes: {M.degenerate.sum():,} "
          f"({100*M.degenerate.mean():.1f}%) — reported separately, never pooled")

    ok = M[(M.reg_promisc_cov >= min_cov)]
    rows = []
    for lab, sub in (("ALL (cov-gated)", ok), ("MULTI-FAMILY (the live stratum)", ok[~ok.degenerate]),
                     ("single-family (inert control)", ok[ok.degenerate])):
        for ax in [*PROMISC_AXES, *[c.replace("reg_promisc", "cur_promisc") for c in PROMISC_AXES]]:
            d = sub[[ax, "delta"]].dropna()
            if len(d) < 60:
                continue
            r = stats.spearmanr(d[ax], d["delta"])
            rows.append({"stratum": lab, "axis": ax, "n": len(d),
                         "rho_vs_delta": r.statistic, "p": r.pvalue,
                         "defn": "sequence" if ax.startswith("reg_") else "curated"})
    R = pd.DataFrame(rows)
    # BH within each stratum × definition family (6 axes) — NOT against MH-201's producer-less 82-scan
    R["q"] = np.nan
    for (s, dfn), idx in R.groupby(["stratum", "defn"]).groups.items():
        p = R.loc[idx, "p"].to_numpy()
        order = np.argsort(p)
        q = np.empty_like(p)
        q[order] = np.minimum.accumulate((p[order] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
        R.loc[idx, "q"] = np.clip(q, 0, 1)

    print("\n=== promiscuity vs the β-over-abundance margin (delta; NEGATIVE = β wins more) ===")
    for s in R.stratum.unique():
        print(f"\n  {s}")
        for _, r in R[R.stratum == s].sort_values("p").iterrows():
            star = "  <-" if r.q < 0.05 else ""
            print(f"    {r.axis:20s} [{r.defn:8s}] n={int(r.n):4d}  rho={r.rho_vs_delta:+.4f}  "
                  f"p={r.p:.3g}  q={r.q:.3g}{star}")

    DEST.mkdir(parents=True, exist_ok=True)
    M.to_csv(DEST / "promiscuity_axis_genes.tsv", sep="\t", index=False)
    R.to_csv(DEST / "promiscuity_axis_scan.tsv", sep="\t", index=False)
    print(f"\n[write] {DEST/'promiscuity_axis_genes.tsv'} · {DEST/'promiscuity_axis_scan.tsv'}")
    return R


if __name__ == "__main__":
    run()
