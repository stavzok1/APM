"""⭐ THE ONE RETENTION IMPLEMENTATION — arithmetic, gate, and coverage report in one place (MH-258).

> **User-directed, 2026-08-19:** *"why is it calculated in so many different places if it's essentially the
> same? logic should be canon."*

⛔⛔ **BUT THE FIRST THING TO SAY IS THAT IT IS *NOT* ONE QUANTITY — AND THAT IS THE DEEPER PROBLEM.**
Auditing every call site (MH-257 → MH-258) found **ten** computations under the name *retention*, spanning
**FOUR DIFFERENT ESTIMANDS** and **FIVE DIFFERENT GATE VALUES** (none · 1e-9/1e-6 · 0.005 · 0.02 · 0.05):

| estimand | ratio | what it asks | scale of the denominator |
|---|---|---|---|
| **BETA**  | `β_deconv / β_core`   | does the FITTED coefficient survive the C block? | β > 0, median ~0.02 |
| **RHO**   | `ρ_adj / ρ_raw`       | does the observed COUPLING survive it?           | correlation, |ρ| ≲ 0.7 |
| **GAP**   | `gap_deconv/gap_core` | does the β-vs-decoy MARGIN survive it?           | a difference, ~1e-2 |
| **STATE** | `ρ_H / ρ_Tsub`        | do HEALTHY weights transfer to tumour?           | ⚠ **not adjusted-vs-raw at all** |

⇒ **`CARD_RUNG_DOCTRINE.md` §5 says *"`retention` names two unrelated quantities"*. It names four.**
Collapsing them into one number would be the axiom-6 error (a naming collision "fixed" by deleting the
distinction). ⇒ **this module makes the LOGIC canon while keeping the ESTIMANDS distinct**: one gated
division, one place to fix it, and a name per estimand that a reader cannot confuse.

⚠⚠ **AND A SINGLE GATE CONSTANT WOULD BE WRONG.** A gate is a statement about the DENOMINATOR'S SCALE.
`RHO_GATE = 0.05` is right for a correlation and catastrophic for β (whose median is ~0.02 — it would drop
most of the card) and far too coarse for a gap (~1e-2). That is why `gate_q` exists: for a denominator with
no natural scale, gate at a QUANTILE of |denominator| and report what that resolved to.

    from mirna_hallmark.learned import retention as RET
    r = RET.ratio(num, den, gate=C.RHO_GATE, name="cptac_prosp_agg_ret_prot")   # absolute gate
    r = RET.ratio(num, den, gate_q=0.10,     name="retention")                  # quantile gate (β scale)
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from mirna_hallmark.config import RHO_GATE

# per-scale gates. ⚠ ONE constant per DENOMINATOR SCALE, not one for everything — see the module docstring.
GAP_GATE = 0.005        # β-vs-decoy margins live at ~1e-2; card_context has used this since MH-144
PROT_GATE = RHO_GATE    # ⚠ `eval/rppa_protein.py` used 0.02; unified to RHO_GATE (both are correlations)

_REPORT: list[dict] = []


def ratio(num, den, *, gate: Optional[float] = None, gate_q: Optional[float] = None,
          name: str = "", quiet: bool = False):
    """`num / den`, NaN wherever |den| falls below the gate. Records what the gate cost.

    Exactly one of `gate` (absolute) or `gate_q` (a quantile of |den|) must be given — a division with NO
    gate is the defect this module exists to prevent, so it is not an available option.

    ⭐ The point is not the arithmetic; it is that the gate is **declared, scale-checked and COUNTED**.
    `report()` prints how many rows each call dropped, so a gate that silently removes half the data is
    visible instead of being discovered later as a suspiciously small n.
    """
    if (gate is None) == (gate_q is None):
        raise ValueError(f"{name or 'retention.ratio'}: give exactly one of gate= or gate_q= "
                         "(an ungated division is the bug this module prevents)")
    n = pd.to_numeric(pd.Series(num), errors="coerce")
    d = pd.to_numeric(pd.Series(den), errors="coerce").reindex(n.index)
    a = d.abs()
    thr = float(gate) if gate is not None else float(a.quantile(gate_q))
    keep = a >= thr
    out = (n / d).where(keep & d.notna() & n.notna())
    defined_in = int((n.notna() & d.notna()).sum())
    _REPORT.append({"name": name, "gate": round(thr, 6),
                    "kind": "abs" if gate is not None else f"q{gate_q:g}",
                    "defined": defined_in, "kept": int(out.notna().sum()),
                    "dropped": defined_in - int(out.notna().sum())})
    if not quiet and defined_in and (defined_in - int(out.notna().sum())) / defined_in > 0.25:
        print(f"  ⚠ retention[{name}]: the |den| >= {thr:.4g} gate drops "
              f"{defined_in - int(out.notna().sum())}/{defined_in} rows "
              f"({1 - out.notna().sum()/defined_in:.1%}) — check the gate is on the right SCALE.")
    return out


def scalar(num: float, den: float, *, gate: float, name: str = "") -> float:
    """Scalar form of `ratio` — several call sites compute one number per gene/arm inside a loop.

    ⚠ Deliberately has NO `gate_q`: a quantile is meaningless for a single value, and silently accepting
    one would be a gate that is not a gate.
    """
    try:
        n, d = float(num), float(den)
    except (TypeError, ValueError):
        return float("nan")
    if not (np.isfinite(n) and np.isfinite(d)) or abs(d) < float(gate):
        return float("nan")
    return n / d


def report() -> pd.DataFrame:
    """What every gated division in this process cost. Empty until `ratio()` has been called."""
    return pd.DataFrame(_REPORT)


def reset() -> None:
    _REPORT.clear()
