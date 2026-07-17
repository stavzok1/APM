"""Phase-1 BOOTSTRAP for the canonical architecture map (ARCHITECTURE_PLAN.md): propose an `axis` tag for every
DISCOVERY_REGISTRY MH-row by keyword classification. **A BOOTSTRAP, NOT THE SOURCE OF TRUTH** — the join key
must be RIGHT (a wrong axis scatters a finding to the wrong place), so the emitted `axis_assignment_draft.tsv`
is DRAFT and needs a human review pass (ambiguous + unassigned + silent over/under-matches). The ongoing
convention is a manual `axis` column per row (ARCHITECTURE_PLAN Phase 4), not this classifier.

Run: .venv/bin/python3 -m mirna_hallmark.analyses.ops.axis_bootstrap
"""
from __future__ import annotations
import re
from pathlib import Path
from mirna_hallmark import config as C

REG = C.SUBPROJECT_DIR / "docs" / "DISCOVERY_REGISTRY.md"
OUT = C.SUBPROJECT_DIR / "docs" / "derived" / "axis_assignment_draft.tsv"

# (pattern, weight) per axis; dict order = tie-break priority (specific → general)
AX = {
 "cn-causal":[("copy.number",2),("\\bcn instrument",3),("2sls",3),("exclusion restriction",3),("hansen",3),("exogenous",2),("germline",2),("locus.cnv",2),("first.stage",2)],
 "protein":[("cptac",3),("\\bprotein",2),("tmt",2),("proteomic",2),("mediation",1)],
 "dcis-ev":[("dcis",3),("extracellular vesicle",3),("pre-malignant",3),("gse59",2),("gse270",2)],
 "external":[("buffa",3),("metabric",3),("scan-b",3),("external cohort",2),("ccle",2)],
 "outcome":[("survival",3),("prognostic",3),("hazard",2),("overall survival",3),("recurrence",2)],
 "subtype":[("pam50",3),("basal",2),("luminal",2),("her2",2),("intrinsic subtype",3),("subtype-strat",3),("subtype-specific",3)],
 "progression":[("gtex",2),("\\bnat\\b",2),("trajectory",3),("progression",3),("healthy.*anchor",2),("normal-like",2)],
 "decoy":[("decoy",3),("matched fake",3),("scrambled",3),("specificity",2),("\\bfake\\b",2),("positive control",1)],
 "discovery":[("discovery",3),("orphan",2),("scan_all",2),("convergent evidence",3),("site.free null",3),("novel",1),("a1.chimeric",3)],
 "attribution":[("shapley",3),("identity vs magnitude",3),("attribution",2),("\\blmg\\b",2),("who owns",2),("canonical regulator",1)],
 "model":[("gibbs",3),("posterior",2),("spike.slab",3),("readout",2),("gauge",3),("half-normal",2),("sampler",2),("niter",2),("shrinkage",2)],
 "edge-existence":[("coupling",2),("anti.correl",2),("partial spearman",2),("permutation null",2),("\\bfdr\\b",1),("repress",1)],
}


def classify(body: str):
    b = body.lower()[:1500]
    scores = {ax: sum(w * len(re.findall(p, b)) for p, w in pats) for ax, pats in AX.items()}
    mx = max(scores.values())
    if mx == 0:
        return "?", 0, False
    top = next(a for a in AX if scores[a] == mx)              # priority tie-break
    second = sorted(scores.values())[-2]
    return top, mx, second >= mx * 0.75


def run() -> None:
    reg = REG.read_text(encoding="utf-8")
    rows = re.findall(r'^\| (MH-[\w/]+) \|(.+?)(?=\n\| MH-|\n\Z)', reg, re.S | re.M)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        f.write("mh\taxis_draft\tscore\tneeds_review\n")
        for mh, body in rows:
            ax, s, amb = classify(body)
            f.write(f"{mh}\t{ax}\t{s}\t{'REVIEW' if (amb or ax=='?') else ''}\n")
    print(f"[axis-bootstrap] wrote {OUT} ({len(rows)} rows) — DRAFT, needs a human review pass")


if __name__ == "__main__":
    run()
