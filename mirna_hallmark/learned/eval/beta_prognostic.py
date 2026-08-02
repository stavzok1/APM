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


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    ap.parse_args()
    run()
