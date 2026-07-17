"""Canonical per-arm expression-DETECTABILITY tiers + the "expressed-arm" floor.

**Why (see `method_dev/arm_expression/ARM_EXPRESSION_FLOOR.md`):** ~half the HE-edge arms are virtually unexpressed in
TCGA-BRCA (cohort-median RPM<1) and 234 never reach RPM>=10 in *any* tumor. Such arms cannot occupy RISC
(median-silent arms are <=2% of the miRNA pool even in the most extreme tumor) and cannot show coupling
(near-flat across samples -> ~zero cross-sample variance). They are noise in every cohort-level analysis.

**The floor is DETECTABILITY, not cohort-median** — a median cut would wrongly discard context-specifically
induced arms. We **KEEP any arm detected at >= ARM_EXPRESSED_MIN_RPM (10) in at least ONE tumor** (max RPM),
because induction in a tumor subset can be real/functional there; we **REMOVE only arms that NEVER reach the
floor in any tumor** (truly silent). This is a **noise filter, not a functional verdict** — surviving low
arms are still naturally down-weighted by the abundance-weighted aggregate, and function is decided by
coupling. Absolute miRNA:target stoichiometry is *not* used (miRNA-seq RPM and mRNA-seq TPM are separate
libraries with no copy-number calibration). `frac_expressed` (fraction of tumors at the floor) is reported
so the weakest conditional arms (detected in very few tumors) are visible.

Tiers: **robust** (cohort-median RPM>=MIN_RPM) > **conditional** (max RPM>=MIN_RPM but median<MIN_RPM —
context-specific/inducible) > **silent** (max RPM<MIN_RPM — never reaches the floor; removed).
`expressed` = robust | conditional.

API: `arm_expression_tiers()` (per-arm table, cached) · `expressed_arms()` (the filter set) ·
`filter_expressed_edges(edges)`. Run `-m mirna_hallmark.arm_expression` to (re)write the canonical table.
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
import mirna_hallmark.data_loaders as D

_OUT = C.OUTPUT_ROOT / "matrices" / "arm_expression_tiers.tsv"
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))


@lru_cache(maxsize=1)
def arm_expression_tiers(min_rpm: float = None, min_frac: float = None) -> pd.DataFrame:
    """Per-arm detectability table (index=arm): median_log2rpm, p90, max, frac_expressed, tier, expressed."""
    min_rpm = C.ARM_EXPRESSED_MIN_RPM if min_rpm is None else min_rpm
    min_frac = C.ARM_EXPRESSED_MIN_FRAC if min_frac is None else min_frac
    thr = np.log2(min_rpm + 1.0)                                   # log2(RPM+1) threshold
    M = D.load_mirna_arms().astype(float).clip(lower=0.0)
    df = pd.DataFrame({
        "median_log2rpm": M.median(axis=1),
        "p90_log2rpm": M.quantile(0.9, axis=1),
        "max_log2rpm": M.max(axis=1),
        "frac_expressed": (M >= thr).mean(axis=1),
    })
    # KEEP anything detected at >= MIN_RPM in >=1 tumor (context-specific induction is real);
    # REMOVE only arms that NEVER reach the floor in any tumor (max < thr) — truly silent.
    df["tier"] = np.where(df["median_log2rpm"] >= thr, "robust",
                          np.where(df["max_log2rpm"] >= thr, "conditional", "silent"))
    df["expressed"] = df["tier"] != "silent"
    _ = min_frac  # reported as frac_expressed; not the keep/remove criterion (that is max-based)
    df.index.name = "arm"
    return df


@lru_cache(maxsize=1)
def expressed_arms(min_rpm: float = None, min_frac: float = None) -> frozenset:
    """Set of arms passing the detectability floor (robust | conditional)."""
    t = arm_expression_tiers(min_rpm, min_frac)
    return frozenset(t.index[t["expressed"]])


def filter_expressed_edges(edges: pd.DataFrame, arm_col: str = "miRNA") -> pd.DataFrame:
    """Drop edges whose arm is silent across the cohort (locus-suffix tolerant)."""
    ok = expressed_arms()
    keep = edges[arm_col].map(lambda a: a in ok or _strip(a) in ok)
    return edges.loc[keep].copy()


def build(out: Path = _OUT) -> pd.DataFrame:
    df = arm_expression_tiers()
    out.parent.mkdir(parents=True, exist_ok=True)
    df.reset_index().to_csv(out, sep="\t", index=False)
    vc = df["tier"].value_counts()
    print(f"[arm-expr] {len(df):,} arms in TCGA matrix | floor: keep if RPM>={C.ARM_EXPRESSED_MIN_RPM:g} in "
          f">=1 tumor; remove if never reached in any tumor")
    print(f"  tiers: robust {vc.get('robust',0)} | conditional {vc.get('conditional',0)} | silent {vc.get('silent',0)}"
          f"  -> expressed (kept) {int(df['expressed'].sum())}")
    print(f"[arm-expr] wrote {out}")
    return df


if __name__ == "__main__":
    build()
