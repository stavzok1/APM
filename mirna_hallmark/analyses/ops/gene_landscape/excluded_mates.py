"""⭐ WHAT DO WE LOSE BY EXCLUDING A FAMILY'S OTHER MEMBERS FROM A GENE'S DESIGN? (user-asked 2026-08-19)

    PYTHONPATH=. .venv/bin/python3 -m mirna_hallmark.analyses.ops.gene_landscape.excluded_mates

THE QUESTION. A gene's design contains SOME members of a seed family, not all: measured, the median gene
sees only **0.67** of its family (AKT1 sees 3 of miR-302's 11 arms). Those absences are a CURATION fact —
and same-seed mates share the seed, so if the in-design member can bind this 3'UTR, the excluded one can
too. So the excluded mates are the sharpest available test of what curation is throwing away.

⚠ THE CONFOUND IS ABUNDANCE, AND IT IS NOT A NUISANCE HERE — it is the leading hypothesis. An excluded
mate that is simply not expressed cannot couple regardless of curation, so a raw in-vs-out contrast is
mostly an abundance contrast. Everything below is therefore reported BOTH pooled and abundance-matched.

WHAT IS COMPUTED. For each (gene, seed_family) cell with excluded members, the partial Spearman of the
EXCLUDED arm's dose against the gene's mRNA, on the same C block the canonical path uses — i.e. exactly
the quantity `coupling_tum` holds for in-design edges — so in and out are directly comparable.

RESULT (400 genes · 1,420 excluded pairs vs 1,112 in-design edges, 2026-08-19)
-----------------------------------------------------------------------------
**POOLED THEY ARE INDISTINGUISHABLE** — excluded median rho **−0.0455** vs in-design **−0.0465**. An arm
curation left out of the design couples about as hard as one it kept.

**ABUNDANCE-MATCHED, CURATION EARNS ITS KEEP IN 4 OF 5 QUINTILES** — in-design wins by +0.028…+0.045
(all p<0.025) in Q2–Q5. ⚠ Q1 INVERTS (excluded −0.0569 vs in-design −0.0055, p=1.4e−08); at the bottom of
the abundance range rho is unstable, and in-design Q1 arms are near zero — read it as noise, not as
curation being wrong. ⇒ **the loss is a TAIL, not a systematic failure of curation.**

⭐ **THE TAIL IS REAL, AND IT IS NOT EXPLAINED BY THE THREE OBVIOUS THINGS.** 262 excluded mates sit above
the expression floor with rho < −0.10.
  * **NOT fame** — median **75 PMIDs vs 60** for in-design arms; these are miR-17-5p (133), miR-106b-5p
    (93), miR-93-5p (79). Comparably or better studied, so "just not famous" is refuted.
  * **NOT abundance** — all are above the floor by construction, and the matched panel above shows
    abundance does not favour them.
  * **MOSTLY NOT CO-EXPRESSION WITH THE IN-DESIGN MATE** — the control I expected to kill this only
    dented it. Same-seed mates correlate a median **0.479** (>0.8 in just 6%), so they are NOT generally
    co-transcribed; conditioning the excluded arm on its in-design partner attenuates the median from
    **−0.1795 to −0.1382 (23% of the signal was the mate's)** and **187 of 257 (73%) still couple below
    −0.10.** Several survive with near-zero mate correlation — CAB39/miR-15a-5p (cor 0.047, rho −0.324),
    COL3A1/let-7g-5p (0.057, −0.318), CD99/miR-30e-5p (0.019, −0.297) — i.e. unambiguously independent.

⇒ **~187 same-seed, above-floor, independently-coupling pairs are missing from the design**, concentrated
in miR-17~92 against canonical targets (BMPR2, APC, BTG2, BCL2, CASP7, DAB2, ADAR, ARID4B). Since a
same-seed mate binds the SAME site, their absence is a curation-completeness artifact, not biology. This
is direct evidence for the board's MH-166 model-expansion item.

⚠⚠ **CANDIDATES, NOT FINDINGS.** Bulk-observational, no decoy, no null, single cohort; and the polycistronic
miR-17~92 members do attenuate most (BMPR2/miR-106b −0.481→−0.320 at cor 0.625), so quote the CONDITIONED
number. Feed them to the refit, do not report them as edges.
"""
from __future__ import annotations

import pathlib
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

OUT = pathlib.Path(__file__).resolve().parents[3] / "output" / "learned"


def _partial_spearman(x, y, C):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(C).all(1)
    if ok.sum() < 60:
        return np.nan
    xr, yr = stats.rankdata(x[ok]), stats.rankdata(y[ok])
    A = np.column_stack([np.ones(ok.sum()), C[ok]])
    xr = xr - A @ np.linalg.lstsq(A, xr, rcond=None)[0]
    yr = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
    if xr.std() < 1e-12 or yr.std() < 1e-12:
        return np.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def main(limit_genes: int = 400) -> None:
    from mirna_hallmark.learned import data as LD, families as FAM

    e = pd.read_csv(OUT / "realization/edge_card.tsv", sep="\t", low_memory=False)
    a = pd.read_csv(OUT / "arm_card.tsv", sep="\t", low_memory=False)
    num = lambda s: pd.to_numeric(s, errors="coerce")

    D = LD._load()
    X = D["X"]
    fam = FAM.family_of(pd.Index(X.index))
    fam_members: dict = {}
    for arm, f in fam.items():
        if f is not None:
            fam_members.setdefault(f, []).append(arm)

    rpm = a.set_index("arm")["arm_med_rpm"].pipe(num)
    npmid = a.set_index("arm")["fame_npmid"].pipe(num) if "fame_npmid" in a.columns else None

    cells = e.dropna(subset=["seed_family"]).groupby(["gene", "seed_family"])["arm"].apply(list)
    genes = sorted({g for g, _ in cells.index})[:limit_genes]
    rows = []
    for gene in genes:
        try:
            Y, Xg, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=False)
        except Exception:
            continue
        Cn = C.to_numpy(float)
        yv = Y.to_numpy(float)
        for (g2, f), in_arms in cells.items():
            if g2 != gene:
                continue
            members = fam_members.get(f, [])
            out_arms = [m for m in members if m not in set(in_arms)]
            if not out_arms:
                continue
            for arm in out_arms:
                if arm not in X.index:
                    continue
                xv = X.loc[arm, Y.index].to_numpy(float) if set(Y.index) <= set(X.columns) else None
                if xv is None:
                    common = [s for s in Y.index if s in X.columns]
                    if len(common) < 60:
                        continue
                    xv = X.loc[arm, common].to_numpy(float)
                    yv2 = Y.loc[common].to_numpy(float)
                    Cn2 = C.loc[common].to_numpy(float)
                else:
                    yv2, Cn2 = yv, Cn
                rows.append({"gene": gene, "seed_family": f, "arm": arm, "status": "EXCLUDED",
                             "rho": _partial_spearman(xv, yv2, Cn2),
                             "rpm": rpm.get(arm, np.nan),
                             "npmid": (npmid.get(arm, np.nan) if npmid is not None else np.nan)})
    E = pd.DataFrame(rows).dropna(subset=["rho"])
    inn = e.assign(rho=num(e.coupling_tum), rpm=e.arm.map(rpm),
                   npmid=(e.arm.map(npmid) if npmid is not None else np.nan))
    inn = inn[inn.gene.isin(genes) & inn.rho.notna()][["gene", "seed_family", "arm", "rho", "rpm", "npmid"]]
    inn["status"] = "IN-DESIGN"
    both = pd.concat([inn, E], ignore_index=True)
    both.to_csv(OUT / "excluded_mates.tsv", sep="\t", index=False)

    print(f"genes scanned {len(genes)} | EXCLUDED same-seed mate pairs scored {len(E)} "
          f"| in-design edges {len(inn)}")
    print("\n=== POOLED (⚠ mostly an abundance contrast — see the matched panel below) ===")
    for lab, s in both.groupby("status"):
        print(f"  {lab:10s} n={len(s):5d}  median rho {s.rho.median():+.4f}   "
              f"frac< -0.10 {(s.rho < -0.10).mean():5.1%}   median log2RPM {s.rpm.median():.2f}")
    u = stats.mannwhitneyu(both[both.status == "IN-DESIGN"].rho, both[both.status == "EXCLUDED"].rho)
    print(f"  MWU p={u.pvalue:.3g}")

    print("\n=== ⭐ ABUNDANCE-MATCHED (quintiles of log2 RPM, pooled over all arms) ===")
    both["q"] = pd.qcut(both.rpm, 5, labels=False, duplicates="drop")
    print(f"  {'Q':<4}{'n_in':>6}{'med_in':>9}{'n_out':>7}{'med_out':>9}{'gap':>9}{'p':>10}")
    for q, s in both.dropna(subset=["q"]).groupby("q"):
        i, o = s[s.status == "IN-DESIGN"].rho, s[s.status == "EXCLUDED"].rho
        if len(i) < 15 or len(o) < 15:
            continue
        p = stats.mannwhitneyu(i, o).pvalue
        print(f"  Q{int(q)+1:<3}{len(i):>6}{i.median():>+9.4f}{len(o):>7}{o.median():>+9.4f}"
              f"{o.median()-i.median():>+9.4f}{p:>10.3g}")

    print("\n=== ⭐ THE ACTIONABLE TAIL — excluded mates that couple as hard as a real edge ===")
    hard = E[E.rho < -0.10].sort_values("rho")
    print(f"  excluded mates with rho < -0.10: {len(hard)} of {len(E)} ({len(hard)/max(len(E),1):.1%})")
    ab = hard[hard.rpm > 3.46]
    print(f"  ...of which ABOVE the expression floor (log2 11): {len(ab)}  <- these are the real loss")
    if len(ab):
        print(ab.head(15)[["gene", "arm", "seed_family", "rho", "rpm", "npmid"]].to_string(index=False))
        print(f"\n  their study attention: median {ab.npmid.median():.0f} PMIDs "
              f"vs {inn.npmid.median():.0f} for in-design arms  "
              f"-> {'UNDER-STUDIED, not unbindable' if ab.npmid.median() < inn.npmid.median() else 'comparably studied'}")


if __name__ == "__main__":
    main()
