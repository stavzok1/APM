"""IS THE LEARNED β BUDGET PROGNOSTIC? — the question the pressure arc's negative left open.

    .venv/bin/python3 -m mirna_hallmark.learned.eval.beta_prognostic [--report]

⭐ WHY. The outcome arc (`analyses/outcome/`, 6 modules) is built on the **§6b-RETIRED pressure heuristic**,
and its verdict is emphatic: **pressure MAGNITUDE is exhaustively non-prognostic — 0 recurrence FDR across
3,368 features** (memory `outcome-prognostic-arc-verdict`). **But the LEARNED β has never been tested
against outcome.** That is board §C item 2, and it became nearly free once MH-201 built the Buffa transfer:
the C blocks, the family-dose matrix, the old-symbol resolver and the per-gene budget all already exist, and
**Buffa carries DRFS with 77 events over 10 years** — the only outcome-bearing independent cohort here.

⭐ IT IS ALSO THE RIGHT PLACE TO ASK IT. MH-201 measured that β's **level** transports into Buffa
(β −0.0186 vs unweighted abundance −0.0022, p=6.2e-13) while its **rank** does not separate from abundance.
So β demonstrably carries gene-level signal in this cohort; whether that signal is **prognostic** is a
different question with its own answer.

✅ **POSITIVE CONTROL, RUN FIRST AND NON-NEGOTIABLE: miR-210 → DRFS.** Buffa 2011's own published finding,
and the repo's standing outcome control. **Measured here: HR = 2.019 per SD, p = 1.0e-07, n=207/77 events**
— the harness recovers a known hypoxia-prognostic miRNA at the right effect size and direction. **Without
this a null is uninterpretable**; the pressure arc's own nulls are credible only because this control fires.

THE ESTIMANDS
-------------
  **P1 PER-GENE** — for each gene, Cox(DRFS ~ M_g) where `M_g(s) = Σ_f β_f·Z_f(s)` is the repressive budget
      on that gene, β FROZEN from TCGA. BH-FDR across genes. Reference: the **unweighted abundance budget**
      on the identical columns, so β is the only difference (MH-201's R1).
  **P2 AGGREGATE** — the per-sample total budget, i.e. "how much miRNA repressive pressure is this tumour
      under", the direct analogue of the retired heuristic's per-sample score.

⚠ **BOTH ARE RUN UNADJUSTED *AND* CLINICALLY ADJUSTED (ER, age, grade, nodes, size — 183 complete rows).**
An unadjusted prognostic hit in breast is nearly uninterpretable: grade and nodes dominate DRFS, and any
proliferation-linked axis will ride them. ⚠ And note the budget's own confound — `prolif` is IN the C block
used to build the coupling, but the BUDGET itself is a miRNA-abundance sum and carries its own
proliferation tone.

⚠ **DEGENERACY (MH-201, axiom 8):** with ONE measurable family the budget is `β·Z` and the abundance
reference is `Z` — for a Cox model the coefficient differs but the *ranking* of samples is identical up to
a positive scale, so β and abundance give the SAME hazard ratio. Multi-family genes are the only place the
weights can act; `n_fam` is reported and the split is applied.
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned/ood_cohort"
MIN_EVENTS = 20

#: module-level memo the forked workers inherit (warm-then-fork; never a per-gene read on NFS, axiom 3a).
_MEM: dict = {}


def _outcome():
    """DRFS time/event + the full clinical block, via the pressure arc's own loader (reuse, don't re-parse)."""
    from mirna_hallmark.analyses.pressure_dev.pressure_prognostic_gene_centric import _buffa_full_clin
    from mirna_hallmark.learned.eval import ood_cohort as OC
    mi = OC.buffa()["mirna"]
    T, E, full = _buffa_full_clin(mi)
    idx = list(mi.columns)
    return T.reindex(idx), E.reindex(idx), full.reindex(idx)


def _cox(x, T, E, cov=None):
    """Cox on a z-scored predictor; returns (HR per SD, p, n, events). `cov` adds clinical adjustment."""
    from lifelines import CoxPHFitter
    d = pd.DataFrame({"t": T, "e": E, "x": (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)})
    if cov is not None:
        d = pd.concat([d, cov], axis=1)
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    d = d[d["t"] > 0]
    if len(d) < 50 or d["e"].sum() < MIN_EVENTS or d["x"].std() < 1e-9:
        return (np.nan,) * 2 + (len(d), int(d["e"].sum()) if len(d) else 0)
    try:
        f = CoxPHFitter(penalizer=0.01).fit(d, "t", "e")
        return float(np.exp(f.params_["x"])), float(f.summary.loc["x", "p"]), len(d), int(d["e"].sum())
    except Exception:
        return (np.nan,) * 2 + (len(d), int(d["e"].sum()))


def gate(T, E) -> dict:
    """⛔ THE POSITIVE CONTROL. miR-210 → DRFS, Buffa 2011's own result. If this does not fire, STOP."""
    from mirna_hallmark.learned.eval import ood_cohort as OC
    mi = OC.buffa()["mirna"].apply(pd.to_numeric, errors="coerce")
    hr, p, n, ev = _cox(mi.loc["hsa-miR-210"].to_numpy(float), T, E)
    ok = bool(np.isfinite(hr) and hr > 1 and p < 1e-3)
    print(f"  ⭐ GATE  miR-210 → DRFS: HR/SD = {hr:.3f}  p = {p:.3g}  (n={n}, events={ev})   "
          f"{'✅ PASS — Buffa 2011 reproduced' if ok else '⛔ FAIL — harness broken, a null below means nothing'}")
    return {"hr": hr, "p": p, "n": n, "events": ev, "passed": ok}


def budgets():
    """(gene → per-sample β budget, abundance budget, n_fam), reusing MH-201's validated construction."""
    from mirna_hallmark.learned.eval import ood_cohort as OC
    FM, beta = OC.buffa_family_matrix(), OC.family_beta(False)
    out = {}
    for g in sorted({gg for gg, _ in FM.index}):
        sub = FM.xs(g, level="gene")
        fams = [f for f in sub.index if (g, f) in beta]
        if not fams:
            continue
        Z = OC._z(sub.loc[fams].to_numpy(float).T)
        b = np.array([beta[(g, f)] for f in fams], float)
        out[g] = (Z @ b, Z.sum(1), len(fams))
    return out


def run(*, report: bool = True) -> pd.DataFrame:
    from mirna_hallmark.learned.eval import ood_cohort as OC
    from mirna_hallmark.stats import bh_fdr

    T, E, full = _outcome()
    print(f"\n{'='*96}\nIS THE LEARNED β BUDGET PROGNOSTIC?  (Buffa DRFS, n={int(T.notna().sum())}, "
          f"{int(E.sum())} events)\n{'='*96}")
    g = gate(T, E)
    if not g["passed"]:
        raise SystemExit("[beta_prognostic] positive control failed — refusing to report a null")

    B = budgets()
    rows = []
    for gene, (mb, ab, nf) in B.items():
        y, _ = OC.buffa_row(gene)
        if y is None:
            continue
        r = {"gene": gene, "n_fam": nf}
        for lab, v in (("beta", mb), ("abund", ab)):
            for adj, cov in (("raw", None), ("clin", full)):
                hr, p, n, ev = _cox(v, T, E, cov)
                r[f"hr_{lab}_{adj}"] = hr; r[f"p_{lab}_{adj}"] = p
        r["n"], r["events"] = n, ev
        rows.append(r)
    R = pd.DataFrame(rows)
    for c in [c for c in R.columns if c.startswith("p_")]:
        R["q" + c[1:]] = bh_fdr(R[c].fillna(1.0).to_numpy())
    OUT.mkdir(parents=True, exist_ok=True)
    R.to_csv(OUT / "beta_prognostic_genes.tsv", sep="\t", index=False)

    if report:
        multi = R[R.n_fam >= 2]
        print(f"\n  --- P1 PER-GENE (n={len(R)} genes; {len(multi)} multi-family where β can act) ---")
        print(f"  {'predictor':22s} {'q<0.05':>7s} {'q<0.10':>7s} {'min q':>9s} {'median |log HR|':>16s}")
        for lab in ("beta", "abund"):
            for adj in ("raw", "clin"):
                q, hr = R[f"q_{lab}_{adj}"], R[f"hr_{lab}_{adj}"]
                print(f"  {lab+' / '+adj:22s} {int((q<0.05).sum()):7d} {int((q<0.10).sum()):7d} "
                      f"{q.min():9.3g} {np.abs(np.log(hr)).median():16.3f}")
        print(f"\n  --- multi-family only (β is mathematically inert at n_fam=1) ---")
        for lab in ("beta", "abund"):
            q = multi[f"q_{lab}_clin"]
            print(f"    {lab:6s} clin-adjusted:  q<0.05 {int((q<0.05).sum()):3d} / {len(multi)}   min q {q.min():.3g}")
        hits = R.nsmallest(10, "q_beta_clin")[["gene", "n_fam", "hr_beta_clin", "q_beta_clin",
                                               "hr_abund_clin", "q_abund_clin"]]
        print(f"\n  --- strongest β genes (clinically adjusted) ---\n{hits.to_string(index=False)}")

        # P2 aggregate
        print(f"\n  --- P2 AGGREGATE (per-sample total budget) ---")
        agg_b = np.nansum(np.array([v[0] for v in B.values()]), axis=0)
        agg_a = np.nansum(np.array([v[1] for v in B.values()]), axis=0)
        for lab, v in (("β total budget", agg_b), ("abundance total", agg_a)):
            for adj, cov in (("raw", None), ("clin", full)):
                hr, p, n, ev = _cox(v, T, E, cov)
                print(f"    {lab:18s} {adj:5s}: HR/SD {hr:.3f}  p={p:.3g}  (n={n}, ev={ev})")
        print(f"\n  ⚠ PRE-REGISTERED: the pressure arc found magnitude non-prognostic (0 FDR / 3,368 "
              f"features); I predicted β would be non-prognostic per-gene too.")
    (OUT / "beta_prognostic_manifest.json").write_text(json.dumps(
        {"module": "mirna_hallmark.learned.eval.beta_prognostic",
         "generated_utc": datetime.now(timezone.utc).isoformat(),
         "cohort": "Buffa GSE22216/22219, n=207 paired, DRFS only, 77 events, 10y",
         "positive_control_miR210": g,
         "estimand": "Cox(DRFS ~ z(budget)), beta FROZEN from TCGA; reference = unweighted abundance "
                     "budget on the IDENTICAL columns",
         "clinical_block": ["er", "age", "grade", "nodes", "size"],
         "n_genes": int(len(R))}, indent=2, default=str) + "\n")
    print(f"\n-> {OUT / 'beta_prognostic_genes.tsv'}")
    return R


#: TCGA arm — ~2x Buffa's events (OS 152 / DFI 99 vs 77) and a richer covariate block.
N_FOLDS = 5
GIBBS_ITER, GIBBS_BURN = 200, 80


def _tcga_one(gene: str):
    """Per-gene TCGA budgets: IN-SAMPLE (β from the shipped card) and OUT-OF-FOLD (β refit per fold).

    ⚠⚠ **WHY BOTH, AND WHY THE OOF ONE IS THE REAL TEST.** In Buffa β is genuinely frozen and external. In
    TCGA **β was fit on these very samples' target mRNA**, so an in-sample budget is optimised to track gene
    g here — and any gene-expression→outcome path can leak into it. The leak is bounded by the coupling
    (|ρ|≈0.1) rather than fatal, but it is real and it is avoidable: refit β within training folds only and
    score held-out patients. **The IN-SAMPLE minus OOF difference IS the leak, measured rather than argued.**
    Per-fold budgets are standardised before pooling (MH-181) — `Z·β` is an arbitrary-scale index.
    """
    import numpy as _np
    from sklearn.linear_model import LinearRegression

    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned import spike_slab as SS
    try:
        Y, X, Cm, w = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return None
    fam = FAM.family_of(pd.Index(X.columns))
    Xf, wf, _ = FAM.collapse_by_family(X, w, fam)
    if Xf.shape[1] < 1 or len(Y) < 200:
        return None
    v, y, Cmat = Xf.to_numpy(float), Y.to_numpy(float), Cm.to_numpy(float)
    n = len(y)
    mu, sd = v.mean(0), v.std(0)
    Z = (v - mu) / (sd + 1e-9); Z[:, sd < 0.1] = 0.0

    beta = _MEM.setdefault("fam_beta", None)
    if beta is None:
        f = pd.read_csv(C.REPO_ROOT / "mirna_hallmark/output/learned/gene_family_card.tsv", sep="\t",
                        usecols=["gene", "family", "beta"])
        beta = dict(zip(zip(f.gene, f.family.astype(str)), f.beta.astype(float)))
        _MEM["fam_beta"] = beta
    b_in = _np.array([beta.get((gene, c), 0.0) for c in Xf.columns], float)
    if not _np.isfinite(b_in).any() or _np.abs(b_in).sum() == 0:
        return None
    ins = Z @ b_in
    abu = Z.sum(1)

    fold = _np.random.default_rng(0).permutation(_np.arange(n) % N_FOLDS)
    oof = _np.full(n, _np.nan)
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        if te.sum() < 20 or tr.sum() < 100:
            continue
        lc = LinearRegression().fit(Cmat[tr], y[tr])
        m_, s_ = v[tr].mean(0), v[tr].std(0)
        Ztr = (v[tr] - m_) / (s_ + 1e-9); Ztr[:, s_ < 0.1] = 0.0
        Zte = (v[te] - m_) / (s_ + 1e-9); Zte[:, s_ < 0.1] = 0.0
        try:
            bb, _, _ = SS._gibbs_posterior(Ztr, -(y[tr] - lc.predict(Cmat[tr])), _np.ones(Ztr.shape[1]),
                                           n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
        except Exception:
            return None
        m2 = Zte @ bb
        s2 = _np.std(m2)
        oof[te] = (m2 - _np.mean(m2)) / s2 if s2 > 1e-12 else 0.0
    return {"gene": gene, "parts": list(Y.index), "insample": ins, "oof": oof,
            "abund": abu, "n_fam": Xf.shape[1]}


def tcga_arm(genes=None, workers: int = 6, endpoints=("os", "dfi")) -> pd.DataFrame:
    """The TCGA arm of the prognostic test — more events, but β is IN-SAMPLE, so OOF is the honest column."""
    from multiprocessing import get_context

    from mirna_hallmark.analyses.outcome.outcome_survival import load_tcga_outcome
    from mirna_hallmark.stats import bh_fdr

    cl = load_tcga_outcome()
    cov = pd.concat([cl["base"], cl["comp"]], axis=1)
    if genes is None:
        f = pd.read_csv(C.REPO_ROOT / "mirna_hallmark/output/learned/gene_family_card.tsv", sep="\t", usecols=["gene"])
        genes = sorted(f.gene.unique())
    # warm every static read in the PARENT before forking (axiom 3a): the family beta map, the
    # learned data matrices, and the CNV cache for the whole universe in ONE batched call - a per-gene
    # cache miss re-scans the ASCAT3 source and cost 2 h earlier in this arc.
    from mirna_hallmark import data_loaders as _D
    from mirna_hallmark.learned import data as _LD
    f0 = pd.read_csv(C.REPO_ROOT / "mirna_hallmark/output/learned/gene_family_card.tsv", sep="\t",
                     usecols=["gene", "family", "beta"])
    _MEM["fam_beta"] = dict(zip(zip(f0.gene, f0.family.astype(str)), f0.beta.astype(float)))
    _D.load_cnv_target_genes(sorted(set(genes) & set(_LD._load()["Y"].index)))
    print(f"[tcga_arm] {len(genes):,} genes · {workers} workers · endpoints {endpoints}", flush=True)
    with get_context("fork").Pool(workers) as p:
        res = [r for r in p.imap_unordered(_tcga_one, genes, chunksize=8) if r]
    print(f"[tcga_arm] budgets built for {len(res):,} genes")
    rows = []
    for r in res:
        row = {"gene": r["gene"], "n_fam": r["n_fam"]}
        for ep in endpoints:
            T, E = cl[f"{ep}_t"].reindex(r["parts"]), cl[f"{ep}_e"].reindex(r["parts"])
            cv = cov.reindex(r["parts"])
            for lab in ("insample", "oof", "abund"):
                hr, pv, n, ev = _cox(r[lab], T, E, cv)
                row[f"hr_{lab}_{ep}"] = hr; row[f"p_{lab}_{ep}"] = pv
            row[f"n_{ep}"], row[f"ev_{ep}"] = n, ev
        rows.append(row)
    R = pd.DataFrame(rows)
    for c in [c for c in R.columns if c.startswith("p_")]:
        R["q" + c[1:]] = bh_fdr(R[c].fillna(1.0).to_numpy())
    OUT.mkdir(parents=True, exist_ok=True)
    R.to_csv(OUT / "beta_prognostic_tcga.tsv", sep="\t", index=False)

    print(f"\n{'='*96}\nTCGA ARM — clinically adjusted (age, T/N/M stage, epi/imm/str/prolif)\n{'='*96}")
    for ep in endpoints:
        print(f"\n  --- {ep.upper()}  (n={int(R[f'n_{ep}'].median())}, events={int(R[f'ev_{ep}'].median())}) ---")
        print(f"    {'predictor':22s} {'q<0.05':>7s} {'q<0.10':>7s} {'min q':>10s}")
        for lab in ("insample", "oof", "abund"):
            q = R[f"q_{lab}_{ep}"]
            tag = ("  ⚠ β fit on these samples" if lab == "insample" else
                   "  ⭐ the honest column" if lab == "oof" else "")
            print(f"    {lab:22s} {int((q<0.05).sum()):7d} {int((q<0.10).sum()):7d} {q.min():10.3g}{tag}")
        ins, oo = int((R[f"q_insample_{ep}"] < 0.05).sum()), int((R[f"q_oof_{ep}"] < 0.05).sum())
        print(f"    ⇒ IN-SAMPLE − OOF = {ins - oo:+d} hits — that difference IS the leak from fitting β here.")
    print(f"\n-> {OUT / 'beta_prognostic_tcga.tsv'}")
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--tcga", action="store_true", help="run the TCGA arm (OS + DFI)")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.tcga:
        tcga_arm(workers=a.workers)
    else:
        run()
