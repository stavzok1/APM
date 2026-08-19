"""THE COMPARTMENT BLOCK — WHICH cell compartment drives a gene's apparent miRNA coupling, not just whether.

    .venv/bin/python3 -m mirna_hallmark.learned.compartment_card [--workers 8]

WHY (user-asked, 2026-08-01). The cards already carry the CONSEQUENCE of composition confounding —
`retention`, `composition_class` — so you can see that a gene's coupling collapses under adjustment. They
did **not** carry the CAUSE, so you could not see WHICH compartment did it. Diagnosing ZEB1 (MH-171:
protein↔CAFs +0.768, miR-200 budget↔CAFs −0.484, product −0.372 vs observed −0.413) required computing
those correlations by hand. This block puts the whole compartment vector on the record — **for every
compartment in the block, not only CAFs**, because the driver is not always CAFs and assuming it is would
repeat the error the block exists to catch.

THE MECHANISM IT EXPOSES, per gene, is a PRODUCT:
    r(target, compartment) × r(miRNA budget, compartment)  ≈  the spurious part of the raw coupling
An epithelial miRNA family falls as stroma rises; a stromal target rises as stroma rises; their product is
a negative "coupling" with no repression in it.

WHAT IS EMITTED
  long  `compartment_profile.tsv`  — gene × compartment × cohort: `r_target`, `r_budget`, and
        `retention_alone` = the coupling that SURVIVES residualising on **that one compartment only**.
  card  summary columns (joined by `card_context`):
        comp_driver          the compartment whose removal ALONE collapses the coupling most
        comp_driver_ret      what fraction of the raw coupling survives removing just that one
        comp_driver_share    1 − comp_driver_ret, i.e. the share that one compartment accounts for
        comp_top_target      the compartment the TARGET loads on most strongly (|r|)
        comp_top_target_r    that correlation
        comp_top_budget      the compartment the miRNA BUDGET loads on most strongly

⚠ READ THESE AS A DIAGNOSIS, NOT A VERDICT (axiom 2a: FLAG, don't delete). A miRNA acting through a
cell-STATE shift PRODUCES a composition change — the composition is then partly the MECHANISM, not a
nuisance. A high `comp_driver_share` says "this gene's coupling lives on that compartment axis"; it does
NOT by itself say the edge is fake. Adjudicate with the well-powered cohort, as MH-107 did for miR-200→ZEB1
mRNA (ρ=−0.209, p=8.7e-12 AFTER composition ⇒ the edge is real, only its bulk-PROTEIN magnitude inflated).

⚠ `purity` and `prolif` sit in the same block but are NOT lineages — treat a `comp_driver` of `purity` as
"tumour content", and of `prolif` as the MH-100 proliferation axis, not as a cell type.
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
LONG = OUT / "compartment_profile.tsv"
SUMM = OUT / "compartment_summary.tsv"
from mirna_hallmark.config import RHO_GATE   # ⭐ ONE home (MH-257); was a local 0.05 here


def _sp(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return spearmanr(a[m], b[m]).correlation if m.sum() > 10 else np.nan


def _resid1(v, x):
    """Residualise v on [1, x] — the leave-ONE-compartment-in design."""
    m = np.isfinite(v) & np.isfinite(x)
    out = np.full(v.shape, np.nan)
    if m.sum() > 4:
        D = np.column_stack([np.ones(m.sum()), x[m]])
        out[m] = v[m] - D @ np.linalg.lstsq(D, v[m], rcond=None)[0]
    return out


def _budget(gene, card, Xmat, fam, rung: str = "family"):
    """β-weighted budget, at an EXPLICIT rung (MH-195). `rung` is now an argument rather than an
    assumption, because the two answer different questions:

      rung="family"  pool arms into seed families, z-score, weight by the FAMILY-rung β. The canonical
                     choice: a family-pooled design column must carry a family-rung β (MH-193).
      rung="arm"     keep each arm its own column, z-score, weight by the ARM-rung β from the card.
                     The right choice when the question is about a SPECIFIC ARM rather than a seed.

    ⚠ These are NOT interchangeable and neither is "more correct" in general — they are different
    estimands (the standing both-rungs directive). What was wrong before was that the rung was IMPLIED
    by which matrix happened to be passed in, so nobody could tell which question was being answered.
    """
    s = card[card.gene == gene]
    arms = [a for a in s.arm if a in Xmat.index]
    if not arms:
        return None
    bet = s.set_index("arm").loc[arms, "beta"].to_numpy(float)
    lin = np.nan_to_num(np.power(2.0, Xmat.loc[arms].to_numpy(float)) - 1.0, nan=0.0)
    # ⛔ SAME DEFECT AS cptac_card, FIXED THE SAME WAY (MH-193): `bfam[f] = b` took the FIRST arm's β
    # as the family's, which is arbitrary because β is fit PER ARM and differs across arms in 99.8% of
    # multi-arm cells (MH-191). A family-POOLED design column must carry the FAMILY-rung β.
    from mirna_hallmark.learned.cptac_card import _family_beta
    acc, bfam, keys, seen = {}, {}, [], {}
    for a, b, row in zip(arms, bet, lin):
        f = str(fam.get(a))
        if f not in acc:
            acc[f] = np.zeros(lin.shape[1]); keys.append(f); seen[f] = []
        acc[f] = acc[f] + row
        seen[f].append(b)
    if rung == "arm":
        # ARM rung: no pooling. Each arm is its own z-scored column weighted by its OWN β.
        A = Xmat.loc[arms].to_numpy(float)
        Za = (A - np.nanmean(A, 1, keepdims=True)) / (np.nanstd(A, 1, keepdims=True) + 1e-9)
        return np.nansum(bet[:, None] * Za, axis=0)
    fb = _family_beta()
    for f in keys:
        v = fb.get((gene, f))
        bfam[f] = float(v) if v is not None and v == v else float(np.mean(seen[f]))
    F = np.log2(1.0 + np.vstack([acc[f] for f in keys]))
    Z = (F - np.nanmean(F, 1, keepdims=True)) / (np.nanstd(F, 1, keepdims=True) + 1e-9)
    return np.nansum(np.array([bfam[f] for f in keys])[:, None] * Z, axis=0)


def build() -> pd.DataFrame:
    from mirna_hallmark.eval import cptac_validation as CV
    from mirna_hallmark.eval import decoy_bench as DB
    from mirna_hallmark.learned import confounders as CF
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.eval import ood_protein as OP

    card = pd.read_csv(OUT / "edge_card_base.tsv", sep="\t",
                       usecols=["gene", "arm", "beta"], low_memory=False).dropna(subset=["beta"])
    fam = DB._ctx()["fam"]
    rows = []

    specs = []
    dd = LD._load()
    st = [s for s in dd["X"].columns if s in dd["Y"].columns]
    specs.append(("tcga_mrna", dd["X"][st], dd["Y"].reindex(columns=st), CF.build_C("tcga", st), st))
    d = OP._cptac("prospective"); L = CV.load_cptac_layers("prospective")
    sp = [s for s in d["arms"].columns if s in L["protein_z"].columns]
    specs.append(("cptac_prot", d["arms"][sp], L["protein_z"].reindex(columns=sp),
                  CF.build_C("cptac", sp), sp))

    for tag, Xm, Ylay, Cb, smp in specs:
        Cb = Cb.reindex(smp).apply(pd.to_numeric, errors="coerce")
        Cb = Cb.fillna(Cb.median(numeric_only=True))
        comps = list(Cb.columns)
        Cn = {c: Cb[c].to_numpy(float) for c in comps}
        genes = [g for g in card.gene.unique() if g in Ylay.index]
        print(f"[compartment] {tag}: {len(genes):,} genes × {len(comps)} compartments (n={len(smp)})")
        for g in genes:
            bud = _budget(g, card, Xm, fam)
            if bud is None:
                continue
            y = Ylay.loc[g].to_numpy(float)
            raw = _sp(bud, y)
            if not np.isfinite(raw):
                continue
            for c in comps:
                rt, rb = _sp(y, Cn[c]), _sp(bud, Cn[c])
                ra = _sp(_resid1(bud, Cn[c]), _resid1(y, Cn[c]))
                rows.append({"cohort": tag, "gene": g, "compartment": c,
                             "r_target": rt, "r_budget": rb, "rho_raw": raw, "rho_alone": ra,
                             "retention_alone": ra / raw if abs(raw) >= RHO_GATE else np.nan,
                             "product": rt * rb})
    P = pd.DataFrame(rows)
    P.to_csv(LONG, sep="\t", index=False)
    print(f"-> {LONG}  ({len(P):,} rows)")
    return P


def summarise(P: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (coh, g), s in P.groupby(["cohort", "gene"]):
        v = s.dropna(subset=["retention_alone"])
        rec = {"cohort": coh, "gene": g, "rho_raw": s.rho_raw.iloc[0]}
        if len(v):
            i = v.retention_alone.idxmin()
            rec.update(comp_driver=v.loc[i, "compartment"],
                       comp_driver_ret=v.loc[i, "retention_alone"],
                       comp_driver_share=1 - v.loc[i, "retention_alone"])
        t = s.dropna(subset=["r_target"])
        if len(t):
            j = t.r_target.abs().idxmax()
            rec.update(comp_top_target=t.loc[j, "compartment"], comp_top_target_r=t.loc[j, "r_target"])
        b = s.dropna(subset=["r_budget"])
        if len(b):
            k = b.r_budget.abs().idxmax()
            rec.update(comp_top_budget=b.loc[k, "compartment"], comp_top_budget_r=b.loc[k, "r_budget"])
        out.append(rec)
    S = pd.DataFrame(out)
    wide = None
    for coh, s in S.groupby("cohort"):
        s = s.drop(columns=["cohort"]).rename(
            columns={c: f"comp_{coh}_{c.replace('comp_', '')}" for c in s.columns if c != "gene"})
        wide = s if wide is None else wide.merge(s, on="gene", how="outer")
    wide.to_csv(SUMM, sep="\t", index=False)
    print(f"-> {SUMM}  ({len(wide):,} genes)")
    return S


def report(P: pd.DataFrame, S: pd.DataFrame) -> None:
    for coh, s in S.groupby("cohort"):
        print(f"\n=== [{coh}] which compartment DRIVES the coupling?  (n={len(s):,} genes) ===")
        d = s.dropna(subset=["comp_driver"])
        vc = d.comp_driver.value_counts()
        for k, n in vc.items():
            sub = d[d.comp_driver == k]
            print(f"  {k:20s} driver for {n:5d} genes ({n/len(d):5.1%})   "
                  f"median share {sub.comp_driver_share.median():.2f}")
        print(f"  --- what the TARGET loads on most ---")
        for k, n in s.dropna(subset=["comp_top_target"]).comp_top_target.value_counts().head(5).items():
            print(f"  {k:20s} {n:5d}")
        print(f"  --- what the miRNA BUDGET loads on most ---")
        for k, n in s.dropna(subset=["comp_top_budget"]).comp_top_budget.value_counts().head(5).items():
            print(f"  {k:20s} {n:5d}")


if __name__ == "__main__":
    argparse.ArgumentParser().parse_args()
    Pp = build()
    report(Pp, summarise(Pp))
