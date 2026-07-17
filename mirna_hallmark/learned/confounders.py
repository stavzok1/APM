"""ONE confounder block for FOUR cohorts — the shared C both cross-cohort channels depend on.

Why this exists
---------------
C was being constructed in three incompatible places: `learned/data.assemble_gene` (TCGA:
CPE + target_cn + mal_prolif + optional deconv), `learned/states.{_cibersortx_state_cov,
_state_metagene_cov}` (NAT from the *pooled* file; GTEx a **metagene fallback**), and
`eval/cptac_validation._covariates` (CPTAC: ESTIMATE purity + CIN, **no composition at all**).

Both channels being built are CROSS-COHORT — the state channel maps GTEx→TCGA, the protein channel
maps CPTAC-prospective→TCGA — and a channel's `pihat` is only commensurate with the model's β if the
two cohorts' β were defined after removing *the same kind* of variance. Three different C builders
made that impossible. This module is the single builder.

What the measurements say (2026-07-12, see LEARNED_MODEL_STATE_CHANNEL_PLAN.md §1)
----------------------------------------------------------------------------------
* **Composition is the load-bearing term.** β(deconv8) vs β(no C): ρ=0.83, and |β| drops 38%.
* **Purity is nearly redundant GIVEN composition.** β(deconv8) vs β(deconv8+CPE): ρ=0.984; CPE alone
  ≈ no-C (ρ=0.943). ⇒ the TCGA-CPE / CPTAC-ESTIMATE asymmetry is **benign** — neither cohort drops its
  purity term, and no commensurability surgery is needed. (CPE and the malignant fraction are NOT
  substitutes: r=0.554, means 0.732 vs 0.235.)
* **The panel / reference barely matters.** β is ρ=0.94 concordant across *entirely different atlases*
  and ρ=0.986 across Wu reference subsamples. So the one residual mismatch — TCGA/GTEx on the 300-cell
  Wu reference (REF-A), CPTAC on the 150-cell one (REF-B) — costs ~1.4% in rank and is absorbed by the
  channel's gauge constant + s² floor. **No re-runs.**

The availability matrix (`availability()`) is the CONTRACT: it says which terms exist in which cohort.
C blocks CANNOT be made column-identical — `CPE`/`target_cn` do not exist in healthy tissue — so the
policy is: **shared core everywhere, cohort-specific extras kept, and the asymmetry TESTED, not assumed.**

⚠ ADDITIVE — nothing is rewired yet. `data.assemble_gene` / `states.assemble_state` /
`cptac_validation` are untouched; migrating them is a separate, ripple-traced step (every persisted
output was computed against the current C).

CLI: `python -m mirna_hallmark.learned.confounders`
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

COHORTS = ("tcga", "nat", "gtex", "cptac")

# --- the deconvolution runs, one per cohort. All Wu-major/9-type: SAME cell-type semantics. -------
# `tcga` deliberately points at the re-exported wired file (byte-identical to
# `tcga_tumor_nat_pooled`, verified max|diff|=0.0) so this builder is bit-compatible with
# `data._deconv()` and introduces ZERO ripple.
_CX = Path("data/external/cibersortx")
_WIRED = Path("mirna_hallmark/output/brca_deconvolution/tcga_cibersortx_fractions.tsv")
# `tcga` and `nat` BOTH come from the pooled run (the wired file carries tumour rows AND '<part>-NAT'
# rows) — matching what `states._cibersortx_state_cov` already reads.
#
# ⚠ NAT deliberately uses the POOLED rows, NOT the standalone `tcga_nat_alone` run. Measured
# (2026-07-12): the two give β_NAT at ρ=0.957 — a second-order difference — but the standalone run has a
# materially WORSE decomposition. With S-mode fit on NAT samples alone there are no cancer cells to
# anchor the malignant signature, so it misassigns normal luminal cells into it: NAT `Cancer Epithelial`
# 0.033 (pooled, the MH-72 negative control) → 0.101 (alone, 3×), eaten straight out of `Normal
# Epithelial` (0.245 → 0.175) — which IS one of the 8 C columns. NAT's Cancer-Epi is misassignment, not
# occult tumour: r=−0.157 with the patient's tumour CPE, but +0.414 with NAT's own Normal Epithelial.
_DECONV_SRC = {
    "tcga":  _WIRED,                                                  # tumour rows of the pooled run
    "nat":   _WIRED,                                                  # '-NAT' rows of the SAME pooled run
    "gtex":  _CX / "gtex_wu_major" / "CIBERSORTx_Adjusted.txt",
    "cptac": _CX / "cptac_wu_alt" / "CIBERSORTx_Adjusted.txt",        # REF-A (300-cell) — matches the other three
}
# Wu reference actually used, by md5 of each run's `*_inferred_refsample.txt`. ⚠ THE FILENAMES LIE — only
# TWO Wu references exist: REF-A `7810db8e…` (`scref_celltype_labeled` ≡ `wu_major_300cells` ≡ `wu_major_alt`)
# and REF-B `d250e9e1…` (`scref_wu_major` ≡ `wu_major_150cells` ≡ `wu_major_C`). ALL FOUR cohorts now sit on
# REF-A (the user added `cptac_wu_alt` 2026-07-12) ⇒ the last reference mismatch is CLOSED. The superseded
# `cptac_wu_major` (REF-B) is kept as a sensitivity run: on the same 133 samples the two agree at r=0.88–0.996
# per cell type with identical fit R (0.866 vs 0.864), so the swap is immaterial — but matched is free.
_REF = {"tcga": "REF-A/300cell", "nat": "REF-A/300cell", "gtex": "REF-A/300cell", "cptac": "REF-A/300cell"}

MALIGNANT_COL = "Cancer Epithelial"                              # held OUT of C: the target is expressed there


def _composition_cols() -> list[str]:
    """The 8 non-malignant columns, in `data._DECONV_COLS` ORDER — so this builder is a drop-in for
    `data._deconv()` (verified: max|diff| = 0.0 over the 1059 shared participants)."""
    from mirna_hallmark.learned import data as LD
    return list(LD._DECONV_COLS)


COMPOSITION_COLS = ["CAFs", "T-cells", "Myeloid", "B-cells", "Endothelial", "PVL",
                    "Normal Epithelial", "Plasmablasts"]         # == data._DECONV_COLS (same order)
WU9 = COMPOSITION_COLS + [MALIGNANT_COL]

# No reference (Wu-9, HBCA-11 or LM22) has an adipocyte class, and GTEx breast is adipose-dominant.
# A marker metagene is the only way to close it — same pattern as data._mycaf / _mal_emt.
# ⚠ DEFAULT OFF, and it must stay off until the over-control test (plan P0b) rules. MEASURED on the 514
# GTEx samples: the 8 Wu fractions already explain **69.6%** of this axis (only 30.4% is new), and the
# metagene is r=−0.71 with `Normal Epithelial` / −0.54 with `Cancer Epithelial` — i.e. it is largely an
# INVERSE-EPITHELIAL axis. Residualising on it therefore risks removing the very compartment the target
# is expressed in — the exact failure mode `data._latent`'s global PCs hit (they killed miR-21→PTEN,
# −0.17→+0.03). Adding it is a coin-flip between de-confounding and over-controlling; gate it, don't assume.
ADIPOSE_MARKERS = ["ADIPOQ", "LEP", "PLIN1", "FABP4", "CFD", "CIDEC", "LIPE", "PPARG"]

_CACHE: dict = {}


# --------------------------------------------------------------------------- participant keying
def _to_participant(cohort: str, index: pd.Index) -> pd.Index:
    """Each cohort's deconvolution keys → the participant key its expression matrices use."""
    s = pd.Index(index.astype(str))
    if cohort == "nat":                                   # 'TCGA-A7-A0CE-NAT' → 'TCGA-A7-A0CE'
        return s.str.replace("-NAT$", "", regex=True)
    if cohort == "gtex":                                  # 'GTEX-1117F-2826-SM-5GZXL' → 'GTEX-1117F'
        return pd.Index(["-".join(x.split("-")[:2]) for x in s])
    return s                                              # tcga: 12-char barcode · cptac: 'X01BR001'


def deconv(cohort: str) -> pd.DataFrame:
    """Wu-9 cell-type fractions for a cohort, participant-keyed (all 9 columns, incl. malignant)."""
    key = f"dec_{cohort}"
    if key not in _CACHE:
        p = _DECONV_SRC[cohort]
        if not p.exists():
            _CACHE[key] = None
        else:
            f = pd.read_csv(p, sep="\t").drop_duplicates("Mixture").set_index("Mixture")
            is_nat = f.index.astype(str).str.contains("-NAT")
            if cohort == "tcga":                          # the pooled file carries both — take the tumours
                f = f[~is_nat]
            elif cohort == "nat":                         # ...and the NAT rows of that SAME run
                f = f[is_nat]
            f.index = _to_participant(cohort, f.index)
            f = f[~f.index.duplicated(keep="first")]
            _CACHE[key] = f[[c for c in WU9 if c in f.columns]].apply(pd.to_numeric, errors="coerce")
    return _CACHE[key]


# --------------------------------------------------------------------------- expression (metagenes)
def expression(cohort: str) -> Optional[pd.DataFrame]:
    """gene × participant log2 expression for a cohort (the substrate for marker metagenes)."""
    key = f"expr_{cohort}"
    if key in _CACHE:
        return _CACHE[key]
    E = None
    try:
        if cohort == "tcga":
            from mirna_hallmark.learned import data as LD
            E = LD._load()["Y"]
        elif cohort == "nat":
            from mirna_hallmark.learned import states as STS
            E = STS.state_matrices("11")[1]
        elif cohort == "gtex":
            from mirna_hallmark.learned import state as ST
            E = ST._gtex_mrna()
            E.columns = _to_participant("gtex", E.columns)
        elif cohort == "cptac":
            from mirna_hallmark import config as C
            p = getattr(C.PATHS, "cptac_rna", None) or Path("data/CPTAC/HS_CPTAC_BRCA_2018_RNA_GENE.cct")
            E = pd.read_csv(p, sep="\t", index_col=0)
            E = E[~E.index.isin(["IDX"])].apply(pd.to_numeric, errors="coerce")
    except Exception:
        E = None
    if E is not None:
        E = E[~E.index.duplicated(keep="first")]
    _CACHE[key] = E
    return E


def _metagene(E: pd.DataFrame, markers: Sequence[str]) -> Optional[pd.Series]:
    """Mean z-scored expression of a marker set (row-z across participants)."""
    mk = [g for g in markers if g in E.index]
    if len(mk) < 2:
        return None
    sub = E.loc[mk]
    z = sub.sub(sub.mean(1), axis=0).div(sub.std(1) + 1e-9, axis=0)
    return z.mean(0)


def _prolif(cohort: str) -> Optional[pd.Series]:
    """Proliferation. In cohorts WITH a malignant compartment (tcga/cptac) the E2F+G2M metagene is
    residualised on the `Cancer Epithelial` fraction → per-malignant-cell proliferation, matching
    `data._malignant_prolif`. In healthy tissue (gtex/nat) 'malignant proliferation' is undefined, so
    the plain metagene is used — an *intentional* asymmetry, recorded in `availability()`."""
    E = expression(cohort)
    if E is None:
        return None
    from mirna_hallmark.hallmark_sets import HallmarkSets
    hs = HallmarkSets.load()
    pg = sorted((set(hs.sets.get("HALLMARK_E2F_TARGETS", [])) |
                 set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))) & set(E.index))
    s = _metagene(E, pg)
    if s is None:
        return None
    d = deconv(cohort)
    if cohort in ("tcga", "cptac") and d is not None and MALIGNANT_COL in d.columns:
        ce = d[MALIGNANT_COL].dropna()
        common = s.index.intersection(ce.index)
        if len(common) > 25:                              # residualise on the malignant fraction
            lr = LinearRegression().fit(ce.loc[common].to_numpy().reshape(-1, 1), s.loc[common])
            s = s.copy()
            s.loc[common] = s.loc[common] - lr.predict(ce.loc[common].to_numpy().reshape(-1, 1))
    return s


def _purity(cohort: str) -> Optional[pd.Series]:
    """Each cohort's OWN validated purity estimate — deliberately NOT harmonised. Measured to be
    nearly redundant given composition (β ρ=0.984 with vs without), so the TCGA-CPE / CPTAC-ESTIMATE
    difference is benign; kept because dropping a validated confounder buys nothing."""
    try:
        if cohort == "tcga":
            from mirna_hallmark.learned import data as LD
            return pd.to_numeric(LD._load()["C"]["CPE"], errors="coerce")
        if cohort == "cptac":
            from mirna_hallmark.eval import cptac_validation as CV
            cl = CV.load_prospective_clinical().drop_duplicates("participant").set_index("participant")
            return pd.to_numeric(cl["purity"], errors="coerce")
    except Exception:
        return None
    return None                                            # nat / gtex: no purity estimate exists


# --------------------------------------------------------------------------- the builder
def build_C(cohort: str, participants: Sequence[str], *, composition: bool = True,
            purity: bool = True, prolif: bool = True, adipose: bool = False) -> pd.DataFrame:
    """The cohort's confounder block, participant-indexed and aligned to `participants`.

    Gene-INDEPENDENT terms only. `target_cn` is per-gene and stays in the caller (`assemble_gene`).
    Missing values are median-imputed (preserves n); a term absent for the cohort is simply omitted —
    see `availability()`.
    """
    if cohort not in COHORTS:
        raise ValueError(f"unknown cohort {cohort!r}; expected one of {COHORTS}")
    idx = pd.Index(list(dict.fromkeys(map(str, participants))))
    out = pd.DataFrame(index=idx)

    if composition:
        d = deconv(cohort)
        if d is not None:
            out = out.join(d.reindex(idx)[[c for c in COMPOSITION_COLS if c in d.columns]])
    if purity:
        p = _purity(cohort)
        if p is not None:
            out = out.join(p.reindex(idx).rename("purity"))
    if prolif:
        s = _prolif(cohort)
        if s is not None:
            out = out.join(s.reindex(idx).rename("prolif"))
    if adipose and cohort == "gtex":                       # GTEx breast is adipose-dominant; no ref has adipocytes
        E = expression(cohort)
        a = _metagene(E, ADIPOSE_MARKERS) if E is not None else None
        if a is not None:
            out = out.join(a.reindex(idx).rename("adipose"))

    out = out.apply(pd.to_numeric, errors="coerce")
    return out.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))


def availability() -> pd.DataFrame:
    """THE CONTRACT — which C terms exist in which cohort, and the coverage each achieves.
    The C blocks are NOT column-identical and cannot be made so; this states the asymmetry explicitly
    so the channels' gauges can test it rather than assume it away."""
    rows = []
    for c in COHORTS:
        d, e, p, s = deconv(c), expression(c), _purity(c), _prolif(c)
        rows.append({
            "cohort": c,
            "wu_reference": _REF[c],
            "n_deconv": 0 if d is None else len(d),
            "composition(8)": "yes" if d is not None else "NO",
            "malignant_frac": "yes" if (d is not None and MALIGNANT_COL in d.columns) else "NO",
            "purity": {"tcga": "CPE", "cptac": "ESTIMATE"}.get(c, "—  (none exists)") if p is not None else "—  (none exists)",
            "prolif": ("malignant-adj" if c in ("tcga", "cptac") else "plain metagene") if s is not None else "NO",
            "adipose": "available, DEFAULT-OFF (over-control risk)" if c == "gtex" else "n/a",
            "target_cn": {"tcga": "yes (ASCAT3)", "cptac": "yes (CNA.cct)"}.get(c, "—  (none)"),
            "n_expr": 0 if e is None else e.shape[1],
        })
    return pd.DataFrame(rows).set_index("cohort")


if __name__ == "__main__":
    with pd.option_context("display.width", 190, "display.max_columns", 20):
        print("\n=== C AVAILABILITY MATRIX (the contract) ===")
        print(availability().to_string())
        print("\n=== build_C() smoke — first 3 participants per cohort ===")
        for c in COHORTS:
            d = deconv(c)
            if d is None:
                print(f"  {c:<6} deconvolution MISSING"); continue
            C_ = build_C(c, list(d.index[:200]))
            print(f"  {c:<6} C shape={str(C_.shape):<10} cols={list(C_.columns)}")
