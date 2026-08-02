"""RETIRED-ESTIMATOR AUDIT — every live call of `regression.fit_gene` (the adaptive lasso), CLASSIFIED.

    .venv/bin/python3 -m mirna_hallmark.eval._retired_lasso_audit

WHY (2026-08-01, user-demanded after the THIRD independent rediscovery). `regression.fit_gene` is the
**adaptive lasso** (`Lasso(positive=True)` on prior-scaled columns) that **MH-115 RETIRED** in favour of
the dense Gibbs. It is still called in ~30 places, and three separate sessions have each "discovered"
one of them and raised it as a bug: **MH-145** (`grid_full`/`grid_oof`), **MH-175** (`ood_protein.fit_M`),
and the `dossier`/`discovery` `he_agg` term. Each rediscovery costs a fresh investigation.

⭐ **THE POINT IS THE CLASSIFICATION, NOT THE COUNT.** A call to a retired estimator is only a defect if
the estimator is the ANSWER. Three classes:

  ESTIMAND    the lasso M *is* the reported quantity ⇒ **a real defect**; the number is not the
              canonical model's.
  NUISANCE    the lasso M is a covariate that absorbs something (`he_agg`), never reported ⇒ **inert
              unless measured otherwise.** MEASURED 2026-08-01 on 183 orphan (gene,arm) pairs: swapping
              `he_agg` lasso→Gibbs moves coupling by −0.0012 (Wilcoxon p=0.255, corr 0.959) ⇒ INERT.
  BASELINE    the lasso is deliberately run *as the comparator* (parity/shuffle/support tests) ⇒
              **correct as written — "fixing" it would delete the experiment.** This is MH-142's trap.

⚠ **DO NOT BULK-REPLACE.** Two of the three classes are fine, and MH-142 records what happens when a
load-bearing term is "fixed" without checking its role first.

HOW TO USE. Run it; anything in ESTIMAND that backs a LIVE claim needs either a re-run on the Gibbs or an
explicit "this is the lasso" label on the claim. The classification below is maintained BY HAND because
it encodes intent, which no grep can infer — if you add a `fit_gene` call, add it here.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

from mirna_hallmark import config as C

SUB = C.REPO_ROOT / "mirna_hallmark"

# path:line -> (class, why). Maintained by hand: the CLASS is intent, not syntax.
CLASSIFIED = {
    # ── ESTIMAND — the lasso M is the reported quantity. Real defects.
    "eval/grid_full.py":               ("ESTIMAND", "coupling grid M — flagged MH-145; relabel or re-run"),
    "eval/grid_oof.py":                ("ESTIMAND", "OOF coupling grid M — flagged MH-145"),
    "eval/within_family_showdown.py":  ("ESTIMAND", "within-family M comparison"),
    "learned/eval/ood_protein.py":     ("ESTIMAND", "fit_M for the CPTAC protein aggregate — found MH-175; "
                                                    "MH-171's retention re-run on Gibbs gave the same verdict"),
    "learned/mvp.py":                  ("ESTIMAND", "MVP aggregate M"),
    "learned/analyses/eb_shrink.py":   ("ESTIMAND", "EB shrinkage target"),
    "learned/states.py":               ("ESTIMAND", "`_aggregate_coupling` M for the NAT/TUM state aggregates "
                                                    "(states.py:288-289). ⚠ verify before citing a state number"),
    "learned/subtype.py":              ("ESTIMAND", "`subtype_coupling` M per subtype. ✅ **MH-165 is NOT "
                                                    "affected** — `genomewide_edge_heterogeneity` reads the "
                                                    "GIBBS card (`readouts_arm_edges.tsv`), verified 2026-08-01"),
    "analyses/misc/_oof_arm_vs_family_eiv.py":  ("ESTIMAND", "`_`-prefixed one-off probe"),
    "analyses/misc/_oof_arm_vs_family_eiv2.py": ("ESTIMAND", "`_`-prefixed one-off probe"),
    "analyses/misc/_oof_arm_vs_family_eiv3.py": ("ESTIMAND", "`_`-prefixed one-off probe"),
    # ── NUISANCE — a covariate, never reported. MEASURED INERT.
    "learned/discovery.py":            ("NUISANCE", "`he_agg` — partials an orphan on what curation captures "
                                                    "(MH-142: load-bearing). Swap measured INERT"),
    "learned/eval/dossier.py":         ("NUISANCE", "`he_agg`, same as discovery. Swap measured INERT "
                                                    "(Δ −0.0012, p=0.255, n=183)"),
    # ── BASELINE — the lasso IS the comparator. Correct as written.
    "learned/spike_slab.py":           ("BASELINE", "lasso-vs-Bayes selection comparison"),
    "learned/eval/bayes_parity.py":    ("BASELINE", "explicit parity harness — the whole point"),
    "learned/eval/shuffled_null.py":   ("BASELINE", "permutation null; only consistency real-vs-shuffled matters"),
    "learned/attribution_eb.py":       ("BASELINE", "lasso SUPPORT set (exact zeros) as the comparator"),
}
ORDER = ["ESTIMAND", "NUISANCE", "BASELINE", "UNCLASSIFIED"]


def main() -> None:
    out = subprocess.run(["grep", "-rn", r"fit_gene\s*(", str(SUB), "--include=*.py"],
                         capture_output=True, text=True).stdout
    sites = {}
    for line in out.splitlines():
        path, ln, _ = line.split(":", 2)
        rel = str(Path(path).relative_to(SUB))
        if "def fit_gene" in line or "output/learned/mh127" in rel:
            continue
        sites.setdefault(rel, []).append(ln)

    print(f"[retired_lasso_audit] {sum(len(v) for v in sites.values())} live calls in {len(sites)} modules\n")
    for cls in ORDER:
        rows = [(f, ls) for f, ls in sorted(sites.items())
                if CLASSIFIED.get(f, ("UNCLASSIFIED", ""))[0] == cls]
        if not rows:
            continue
        n = sum(len(ls) for _, ls in rows)
        head = {"ESTIMAND": "⛔ ESTIMAND — the lasso IS the reported quantity (REAL DEFECTS)",
                "NUISANCE": "⚠ NUISANCE — a covariate, never reported (MEASURED INERT)",
                "BASELINE": "✅ BASELINE — the lasso is the comparator (CORRECT AS WRITTEN)",
                "UNCLASSIFIED": "❓ UNCLASSIFIED — NEW call site, classify it in this file"}[cls]
        print(f"{head}   [{n} calls]")
        for f, ls in rows:
            print(f"    {f}:{','.join(ls)}\n        {CLASSIFIED.get(f, ('', 'no entry'))[1]}")
        print()
    unk = [f for f in sites if f not in CLASSIFIED]
    print("=" * 96)
    print("  ⚠ A call to a retired estimator is a defect ONLY in the ESTIMAND class. Do NOT bulk-replace:")
    print("    MH-142 records what happens when a load-bearing term is 'fixed' without checking its role.")
    if unk:
        print(f"  ❗ {len(unk)} UNCLASSIFIED site(s) — add them to CLASSIFIED with their intent.")
    print("=" * 96)


if __name__ == "__main__":
    main()
