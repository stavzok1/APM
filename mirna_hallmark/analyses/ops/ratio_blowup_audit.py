"""⭐ FIND EVERY RATIO EVALUATED WHERE ITS DENOMINATOR VANISHES — axiom 5's class, swept in one pass.

Sibling of `nan_bool_audit.py`. That one finds booleans that cannot say "unknown"; this one finds the
**quantitative** version of the same disease: a ratio whose denominator is allowed near zero, so a handful
of rows carry values hundreds of times the bulk of the distribution.

⛔ **WHY IT MATTERS MORE THAN IT LOOKS.** The recorded cost in this subproject is not "an ugly outlier" —
it is inverted headlines. `share` read **999%** because a gene's βs cancelled to Σβ≈0 (MH-119); a ratio of
two SDs printed **112.86×** where the honest gated value was **~1.4×** (MH-146); `top_identity` reached
**+740** on the delivered gene card and made an attribution ranking meaningless (MH-249). ⇒ *"a statistic
evaluated where its denominator vanishes is not a finding; it is a coin-flip with decimals."*

**THE TEST — `TAIL = max|x| / p99|x|`, and it needs no knowledge of the producer.** A well-behaved quantity
sits at 1–3: its largest value is close to its 99th percentile. A vanishing-denominator ratio has a max
detached from its own bulk. ⭐ **CALIBRATED ON COLUMNS WHOSE STATUS WAS ALREADY KNOWN, and the separation
is clean with no overlap:**

    KNOWN BAD   identity 88.0 · identity_deconv 61.8 · identity_allocated 94.9 · top_identity 134.6 ·
                comp_tcga_mrna_driver_ret 8.6
    KNOWN GOOD  beta 2.0 · coupling_tum 1.3 · z 2.3 · cal_z 2.3 · ctx_ceiling 1.5 ·
                identity_abs 1.0 (clipped) · concentration 1.0 (bounded)

⇒ `TAIL_THR = 5.0` sits in the empty gap between 2.3 and 8.6. ⚠ It is a HEURISTIC, and a genuinely
heavy-tailed measurement (a read count, a target count) will trip it honestly — hence `--counts` to include
the count-like columns that are excluded by default, and hence the output is a REVIEW QUEUE, not a verdict.

A second, independent signal: **BOUND** — for a column whose NAME promises a range (`*_share`, `*_frac*`,
`*_pct*`, `*_ret*`, `*coherence*`, `identity*`), how many rows fall outside it. A name that promises [0,1]
and delivers −4.3 is a defect regardless of its tail shape.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.ratio_blowup_audit
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.ratio_blowup_audit --counts   # include count-like cols
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.ratio_blowup_audit --validate # re-run the calibration
"""
from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd

from mirna_hallmark.learned.card_rungs import OUT

CARDS = {"edge": "realization/edge_card.tsv", "gene": "realization/gene_card.tsv",
         "gene_family": "gene_family_card.tsv", "arm": "arm_card.tsv",
         "seed_family": "seed_family_card.tsv"}
TAIL_THR = 5.0
# counts, reads, distances and PMID tallies are honestly heavy-tailed — excluded unless --counts
COUNTY = re.compile(r"(^n_|_n_|_n$|_reads$|_count$|_studies$|_pmid|_degree|_size$|_len$|_dist$|_bp$"
                    r"|_rpm$|_targetome|_w$)")
# ⛔ NOTE WHAT IS *NOT* HERE: `*_ret` / `*_retention`. The first version of this module bounded them to
# [0,1] and produced 14 false positives at the top of the queue. **RETENTION IS NOT BOUNDED** — a sign flip
# under adjustment gives a negative value and amplification gives >1, both real (MH-253: 21.8% flip, 21.5%
# amplify on CPTAC protein), and **41.1% of `realized_retention` still sits outside [0,1] AFTER a correct
# denominator gate.** Its pathology shows up in TAIL, not in a bound. A bound test is only honest where the
# name promises an actual fraction.
BOUNDED = {  # name pattern -> (lo, hi) the NAME promises
    re.compile(r"(_share$|_frac_|_frac$|_pct$|coherence)"): (0.0, 1.0),
    re.compile(r"^identity"): (-1.0, 1.0),
}

# ⭐ ACCEPTED — a queued column with a RECORDED reason, so the report converges on what is actually open.
# ⚠ An entry here is a claim that someone looked, not that the number is clean. Each says which.
ACCEPTED = {
    # (a) FIXED AT SOURCE, awaiting the queued canonical rebuild. Their denominator is NOT on the delivered
    #     card, so they cannot be repaired in place (unlike `realized_retention`, which was).
    ("edge", "retention_rho"): "gated at source (MH-257); pending the canonical rebuild",
    ("edge", "own_specific_frac"): "gated at source (MH-257); pending the canonical rebuild",
    ("arm", "wshift_own_specific_frac"): "derived from own_specific_frac; gated at source, pending rebuild",
    # (b) DELIBERATELY UNGATED, because the RAW value is the primary record and a GATED COMPANION exists.
    #     Gating `identity` in place would silently change every attribution analysis; the doctrine is to
    #     ship the gated readout beside it (MH-249/255) and document the raw one loudly.
    ("edge", "identity"): "raw Shapley share; use identity_reliable / identity_abs (MH-255)",
    ("gene_family", "identity"): "raw Shapley share at family rung; use identity_reliable",
    ("edge", "identity_deconv"): "raw; no gated companion yet — OPEN, see MH-255",
    ("gene_family", "identity_deconv"): "raw; no gated companion yet — OPEN, see MH-255",
    ("edge", "identity_allocated"): "raw; no gated companion yet — OPEN, see MH-255",
    ("gene", "top_identity"): "raw; use top_identity_gated (MH-249)",
    # (c) KNOWN, with the fix recorded as a threshold change rather than a code defect.
    ("gene", "comp_tcga_mrna_driver_ret"): "RHO_GATE 0.05 too low here; quote the 0.10-gated value (MH-254)",
    ("gene_family", "comp_tcga_mrna_driver_ret"): "as the gene twin (MH-254)",
    ("gene", "comp_tcga_mrna_driver_share"): "= 1 - driver_ret exactly; same gate question (MH-254)",
    ("gene_family", "comp_tcga_mrna_driver_share"): "as the gene twin (MH-254)",
    ("gene_family", "comp_cptac_prot_driver_share"): "as above, protein layer (MH-254)",
    ("gene", "comp_cptac_prot_driver_share"): "as the gene_family twin, protein layer (MH-254)",
    # (d) NOT A RATIO AT ALL — an honest unbounded quantity whose tail is real and worth reading.
    ("edge", "sd_arm"): ("posterior SD, no upper bound. Tail is 6 rows (max 1.302 vs median 0.0036, 360x) "
                         "on miR-302 arms — LEFTY2 beta is UNIDENTIFIED, not mis-divided. Read, don't gate."),
}

KNOWN = [("edge", "identity", "BAD"), ("edge", "identity_deconv", "BAD"),
         ("edge", "identity_allocated", "BAD"), ("gene", "top_identity", "BAD"),
         ("gene", "comp_tcga_mrna_driver_ret", "BAD"), ("edge", "beta", "good"),
         ("edge", "coupling_tum", "good"), ("edge", "z", "good"), ("edge", "cal_z", "good"),
         ("gene", "ctx_ceiling", "good"), ("edge", "identity_abs", "good"),
         ("gene", "concentration", "good")]


def _tail(x: pd.Series) -> float:
    a = x.abs()
    p99 = float(a.quantile(0.99))
    return float(a.max()) / p99 if p99 > 0 else np.inf


def audit(counts: bool = False) -> pd.DataFrame:
    rows = []
    for card, rel in CARDS.items():
        p = OUT / rel
        if not p.exists():
            continue
        d = pd.read_csv(p, sep="\t", low_memory=False)
        for c in d.columns:
            if COUNTY.search(c) and not counts:
                continue
            x = pd.to_numeric(d[c], errors="coerce").dropna()
            if len(x) < 30 or x.nunique() < 10:
                continue
            t = _tail(x)
            lo = hi = None
            for pat, (a, b) in BOUNDED.items():
                if pat.search(c):
                    lo, hi = a, b
                    break
            out = int(((x < lo) | (x > hi)).sum()) if lo is not None else 0
            if t > TAIL_THR or out:
                rows.append({"card": card, "column": c, "n": len(x), "tail": round(t, 1),
                             "p99": round(float(x.abs().quantile(0.99)), 3),
                             "max": round(float(x.abs().max()), 3),
                             "bound": f"[{lo:g},{hi:g}]" if lo is not None else "—",
                             "n_out": out})
    q = pd.DataFrame(rows)
    return q if q.empty else q.sort_values(["tail", "n_out"], ascending=False)


def validate() -> int:
    """Re-run the calibration. The thresholds are only trustworthy while this still separates."""
    bad, good = [], []
    for card, col, lab in KNOWN:
        p = OUT / CARDS[card]
        if not p.exists():
            continue
        d = pd.read_csv(p, sep="\t", low_memory=False)
        if col not in d.columns:
            print(f"  ⚠ {card}.{col} absent — calibration incomplete")
            continue
        x = pd.to_numeric(d[col], errors="coerce").dropna()
        t = _tail(x)
        (bad if lab == "BAD" else good).append((col, t))
        print(f"  {col:<30}{lab:<6}tail {t:>8.1f}  {'⛔ fires' if t > TAIL_THR else '✅ clear'}")
    ok = all(t > TAIL_THR for _, t in bad) and all(t <= TAIL_THR for _, t in good)
    lo_bad = min((t for _, t in bad), default=np.inf)
    hi_good = max((t for _, t in good), default=0.0)
    print(f"\n  worst BAD {lo_bad:.1f}  vs  worst good {hi_good:.1f}   threshold {TAIL_THR}")
    print("  ✅ CALIBRATION HOLDS — the threshold sits in an empty gap." if ok else
          "  ⛔ CALIBRATION BROKEN — a known case is on the wrong side. Re-tune before trusting a run.")
    return 0 if ok else 1


def main() -> int:
    if "--validate" in sys.argv:
        return validate()
    q = audit("--counts" in sys.argv)
    if q.empty:
        print("✅ no column trips the tail or bound test.")
        return 0
    q["accepted"] = [ACCEPTED.get((r.card, r.column), "") for r in q.itertuples()]
    if "--all" not in sys.argv:
        acc_n = int((q.accepted != "").sum())
        q_open = q[q.accepted == ""]
        print(f"({acc_n} queued column(s) accepted with a recorded reason — `--all` to list them)\n")
        q = q_open
    if q.empty:
        print("✅ nothing OPEN — every queued column has a recorded reason.")
        return 0
    print(f"⚠ {len(q)} column(s) queued — heavy tail (max/p99 > {TAIL_THR}) or outside a name-promised "
          f"bound:\n")
    print(f"  {'card':<12}{'column':<34}{'n':>6}{'tail':>8}{'p99':>10}{'max':>12}{'bound':>10}{'out':>7}")
    for r in q.itertuples():
        print(f"  {r.card:<12}{r.column:<34}{r.n:>6}{r.tail:>8.1f}{r.p99:>10.3f}{r.max:>12.3f}"
              f"{r.bound:>10}{r.n_out:>7d}")
    print("\n  ⇒ gate the denominator and report the GATED statistic (axiom 5b), or sweep the threshold\n"
          "    and quote a value from where the arms agree (5c). A heavy tail on an honest measurement\n"
          "    (counts, reads) is not a defect — those are excluded by default; `--counts` shows them.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
