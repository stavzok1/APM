"""⭐ DID CURATION NEVER LOOK, OR LOOK AND DOWNGRADE? — and what does ADDING the excluded mates do to the
family-level aggregate coupling? (user-asked 2026-08-19)

    PYTHONPATH=. .venv/bin/python3 -m mirna_hallmark.analyses.ops.gene_landscape.excluded_mates_evidence

PART A — THE EVIDENCE BASIS OF EACH EXCLUDED MATE. `excluded_mates.py` showed ~187 same-seed, above-floor
pairs coupling independently of their in-design partner. This asks WHY each one is absent, by joining every
excluded (arm, gene) pair to the full PMID-deduped ledger (1.1M rows; miRTarBase x TarBase-v9) and sorting
it into four states that mean different things:

    NEVER_LOOKED    no ledger row at all              -> a coverage hole; curation never assayed this pair
    WEAK_ONLY       rows exist, every one `weak`      -> curation LOOKED and DOWNGRADED
    HT_ONLY         non-weak, but no low-throughput   -> supported by ago_clip / qpcr_rna / chimeric only,
                    functional assay                     i.e. BINDING evidence without a repression assay
    LT_FUNC         reporter / western / proteomics   -> ⛔ properly evidenced AND STILL not in the design;
                                                         if this class is non-empty it is a DESIGN-CONSTRUCTION
                                                         defect, not a curation gap

⚠ The four are ordered by how much they would move the design if acted on, NOT by how likely they are.

PART B — THE AGGREGATE COST. The family rung is what the model actually fits, so the operational question
is not "does this arm couple" but "does the FAMILY's aggregate coupling to this gene improve when the
excluded mates are added back". Computed as the partial Spearman of the summed linear dose against the
gene's mRNA, in-design arms only vs all family members, on the same C block.

⚠ ADDING AN ARM CANNOT BE ASSUMED TO HELP. A mate that is abundant but uncoupled DILUTES the aggregate, so
the delta is signed and both tails are reported. That is the honest form of the question.
"""
from __future__ import annotations

import pathlib
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

OUT = pathlib.Path(__file__).resolve().parents[3] / "output" / "learned"
LT_FUNC = {"reporter", "western", "proteomics"}


def _pspear(x, y, C):
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


def main() -> None:
    from mirna_hallmark.learned import data as LD, families as FAM
    from mirna_hallmark.learned.evidence import ledger as LG

    E = pd.read_csv(OUT / "excluded_mates.tsv", sep="\t")
    exc = E[E.status == "EXCLUDED"].copy()
    L = LG.build_ledger()

    # ---------------- PART A ----------------
    key = L.assign(_k=L.arm + "||" + L.gene)
    g = key.groupby("_k")
    agg = pd.DataFrame({
        "n_rows": g.size(),
        "n_nonweak": g["weak"].apply(lambda s: int((~s.astype(bool)).sum())),
        "n_ltfunc": g.apply(lambda d: int(((~d.weak.astype(bool)) & d.assay_class.isin(LT_FUNC)).sum())),
        "classes": g["assay_class"].apply(lambda s: ",".join(sorted(set(s)))),
        "n_pmid": g["pmid"].nunique(),
    })
    exc["_k"] = exc.arm + "||" + exc.gene
    exc = exc.join(agg, on="_k")
    exc["n_rows"] = exc.n_rows.fillna(0)

    def state(r):
        if r.n_rows == 0:
            return "NEVER_LOOKED"
        if r.n_ltfunc > 0:
            return "LT_FUNC"
        if r.n_nonweak > 0:
            return "HT_ONLY"
        return "WEAK_ONLY"

    exc["ev_state"] = exc.apply(state, axis=1)
    exc.to_csv(OUT / "excluded_mates_evidence.tsv", sep="\t", index=False)

    print("=" * 96)
    print("PART A — WHY IS EACH EXCLUDED MATE ABSENT?  (n=%d excluded pairs)" % len(exc))
    print("=" * 96)
    order = ["NEVER_LOOKED", "WEAK_ONLY", "HT_ONLY", "LT_FUNC"]
    for st in order:
        s = exc[exc.ev_state == st]
        if not len(s):
            print(f"  {st:14s} n=    0")
            continue
        print(f"  {st:14s} n={len(s):5d} ({len(s)/len(exc):5.1%})  median rho {s.rho.median():+.4f}  "
              f"frac< -0.10 {(s.rho < -0.10).mean():5.1%}  median log2RPM {s.rpm.median():5.2f}")
    print("\n  ⭐ restricted to the ACTIONABLE tail (above floor AND rho < -0.10):")
    tail = exc[(exc.rpm > 3.46) & (exc.rho < -0.10)]
    for st in order:
        s = tail[tail.ev_state == st]
        print(f"     {st:14s} n={len(s):4d} ({len(s)/max(len(tail),1):5.1%})"
              + (f"   e.g. {', '.join((s.gene + '/' + s.arm).head(3))}" if len(s) else ""))
    lt = tail[tail.ev_state == "LT_FUNC"]
    if len(lt):
        print(f"\n  ⛔⛔ {len(lt)} pair(s) carry LOW-THROUGHPUT FUNCTIONAL evidence and are STILL not in the")
        print("      design — that is a DESIGN-CONSTRUCTION defect, not a curation gap. Inspect these first:")
        print(lt.head(12)[["gene", "arm", "seed_family", "rho", "rpm", "n_pmid", "classes"]].to_string(index=False))

    # ---------------- PART B ----------------
    print("\n" + "=" * 96)
    print("PART B — WHAT DOES ADDING THE EXCLUDED MATES DO TO THE FAMILY AGGREGATE?")
    print("=" * 96)
    e = pd.read_csv(OUT / "realization/edge_card.tsv", sep="\t", low_memory=False)
    X = LD._load()["X"]
    fam = FAM.family_of(pd.Index(X.index))
    members: dict = {}
    for arm, f in fam.items():
        if f is not None:
            members.setdefault(f, []).append(arm)
    cells = e.dropna(subset=["seed_family"]).groupby(["gene", "seed_family"])["arm"].apply(list)
    genes = sorted(exc.gene.unique())
    rows = []
    for gene in genes:
        try:
            Y, _Xg, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", deconv=False)
        except Exception:
            continue
        common = [s for s in Y.index if s in X.columns]
        if len(common) < 60:
            continue
        yv = Y.loc[common].to_numpy(float)
        Cn = C.loc[common].to_numpy(float)
        for (g2, f), in_arms in cells.items():
            if g2 != gene:
                continue
            ins = [m for m in in_arms if m in X.index]
            allm = [m for m in members.get(f, []) if m in X.index]
            if len(ins) < 1 or len(allm) <= len(ins):
                continue
            lin = lambda arms: np.nansum(np.power(2.0, X.loc[arms, common].to_numpy(float)) - 1.0, axis=0)
            r_in, r_all = _pspear(lin(ins), yv, Cn), _pspear(lin(allm), yv, Cn)
            if np.isfinite(r_in) and np.isfinite(r_all):
                rows.append({"gene": gene, "seed_family": f, "n_in": len(ins), "n_added": len(allm) - len(ins),
                             "rho_in": r_in, "rho_all": r_all, "delta": r_all - r_in})
    A = pd.DataFrame(rows)
    if A.empty:
        print("  (no cells scored)")
        return
    A.to_csv(OUT / "family_aggregate_with_mates.tsv", sep="\t", index=False)
    print(f"  (gene, family) cells scored: {len(A)}")
    print(f"  delta = rho(all members) - rho(in-design).  NEGATIVE = adding the mates STRENGTHENS repression.")
    print(f"    median delta {A.delta.median():+.4f}   improved {(A.delta < 0).mean():.1%}   "
          f"diluted {(A.delta > 0).mean():.1%}")
    w = stats.wilcoxon(A.delta)
    print(f"    Wilcoxon vs 0: p={w.pvalue:.3g}")
    print(f"\n  ⭐ THE MOST EXTREME GAINS — adding the mates takes the family aggregate DOWN the most:")
    print(A.nsmallest(15, "delta")[["gene", "seed_family", "n_in", "n_added", "rho_in", "rho_all", "delta"]]
          .to_string(index=False))
    print(f"\n  ⚠ AND THE OTHER TAIL — where adding them DILUTES the aggregate (the honest counterweight):")
    print(A.nlargest(8, "delta")[["gene", "seed_family", "n_in", "n_added", "rho_in", "rho_all", "delta"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
