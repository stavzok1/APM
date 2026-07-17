"""In-situ -> invasive TIMING of the healthy->tumour miRNA-abundance signature (DCIS).

This is the *correct* use of the GSE59247 DCIS/IBC cohort -- NOT a linear
``healthy->NAT->DCIS->tumour`` index (that order is biologically wrong: NAT is a
field-effect normal, not a DCIS precursor; DCIS->IDC is non-obligate and grade/subtype
heterogeneous; the cohorts are cross-sectional + cross-platform). Instead we ask, per
arm, **when** along the in-situ -> invasive boundary an acquired healthy->tumour change
happens:

  EARLY  (pre-invasive)        already settled by the in-situ stage -> no further DCIS->IBC
                               movement (q >= alpha for DCIS-vs-IBC).
  LATE   (invasion-coupled)    still moving DCIS->IBC in the *acquired direction*
                               (q < alpha, concordant sign).
  DIVERGENT (non-monotonic)    moves DCIS->IBC *against* the acquired direction.

Design choices that keep this honest:
- **Same modality.** The DCIS array measures miRNA *abundance*, so the acquired
  reference is the abundance trajectory (``gtex_to_tumor_pct`` from the GTEx-anchored
  DE landscape), NOT the pressure dHT (which folds in edge weights the array can't see).
  The pressure-axis ``dHT``/``acquired_gainer`` are carried along as annotation only.
- **Same platform for the contrast.** DCIS-vs-IBC is computed *within* GSE59247
  (Agilent GPL15019), so the magnitude/significance never cross platforms.
- **Cross-platform only for DIRECTION.** TCGA enters solely as the acquired *sign*
  (rank-based), exactly as the loader documents.
- **No healthy anchor on GPL15019**, so "EARLY" = "no residual in-situ->invasive
  movement" under a monotone-progression assumption; we corroborate it with the
  intra-DCIS *grade* gradient (does the arm already track severity inside in-situ
  disease?), which needs no healthy anchor.
- **v16 -> modern arm bridge.** GSE59247 uses miRBase-v16 names (``hsa-let-7a``,
  ``hsa-miR-29a*``, ``_v16.0`` tags); TCGA uses arm-suffixed names. We bridge with an
  auditable rule (exact arm-suffix match first; else stem grouping with abundance
  resolution: v16 non-star = guide = most-abundant TCGA arm, ``*`` = passenger) and
  record ``bridge_method`` per arm; ambiguous/unmapped arms are dropped, not guessed.
- **Outcome is underpowered** (3 recurrence / 2 DMFS events) -> reported, not modelled.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.dcis_ev.dcis_timing
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

from analysis.utils.common.loaders import partial_spearman
from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_coupling import (
    EPI_MARKERS, IMMUNE_MARKERS, STROMA_MARKERS, _metagene, _prolif_metagene,
)
from mirna_hallmark.analyses.dcis_ev.dcis_geo_loader import load_gse59246, load_gse59247
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr

_MARKER_SET = set(EPI_MARKERS) | set(IMMUNE_MARKERS) | set(STROMA_MARKERS)


def composition_covariates(mrna: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Per-sample composition (epi/immune/stroma) + proliferation + BATCH(site) covariates.

    Metagenes from the paired mRNA (same markers as ``cross_state_coupling``); batch =
    collection ``source`` site dummies (reference dropped). Shared by the DCIS coupling and
    timing modules so a stromal/myoepithelial arm's signal is not a cell-fraction/site
    artifact (the MH-39 confound).
    """
    cov = pd.DataFrame({
        "epi": _metagene(mrna, EPI_MARKERS),
        "immune": _metagene(mrna, IMMUNE_MARKERS),
        "stroma": _metagene(mrna, STROMA_MARKERS),
        "prolif": _prolif_metagene(mrna, HallmarkSets.load()),
    })
    if "source" in meta.columns:
        site = meta["source"].reindex(cov.index).astype("object").fillna("NA")
        dummies = pd.get_dummies(site, prefix="site", drop_first=True).astype(float)
        cov = cov.join(dummies)
    return cov

OUT_DIR = C.OUTPUT_ROOT / "dcis_timing"
DE_LANDSCAPE = (C.TISSUE_REFERENCE_DIR / "cross_state_landscape"
                / "mirna_arm_de_landscape.tsv")
STATE_CLASS = (C.TISSUE_REFERENCE_DIR / "mirna_state_class" / "mirna_state_class.tsv")

ALPHA = C.FDR_ALPHA
MIN_PER_GROUP = 4          # min samples per group for a DCIS-vs-IBC test
_VER_RE = re.compile(r"_v\d+(?:\.\d+)?$", re.I)        # _v16.0 platform tag
_ARM_SUFFIX_RE = re.compile(r"-(?:5p|3p)(?:\.\d+)?$")  # -5p / -3p / -3p.1
_DOTN_RE = re.compile(r"\.\d+$")


# --------------------------------------------------------------------------- #
# arm-name bridge: miRBase v16 (GSE) -> modern arm (TCGA)
# --------------------------------------------------------------------------- #
def _stem(arm: str) -> str:
    """Hairpin stem: drop arm suffix (-5p/-3p[.N]) and trailing .N and star."""
    a = _ARM_SUFFIX_RE.sub("", arm)
    a = _DOTN_RE.sub("", a)
    return a.rstrip("*")


def bridge_v16_to_modern(
    gse_arms, tcga_arms, abundance: Dict[str, float]
) -> pd.DataFrame:
    """Map GSE59247 v16 probe names -> modern TCGA arm names.

    Returns a frame [gse_arm, modern_arm, bridge_method, is_star]. ``abundance``
    is a {modern_arm: tumour median} dict used to pick the guide (most-abundant)
    arm of a stem. Methods: exact_arm | abundance_guide | abundance_star;
    unmapped GSE arms are omitted (never guessed).
    """
    tcga = list(tcga_arms)
    tcga_set = set(tcga)
    # group modern arms by hairpin stem, ordered by abundance (guide first)
    stem_arms: Dict[str, list] = {}
    for a in tcga:
        stem_arms.setdefault(_stem(a), []).append(a)
    for s in stem_arms:
        stem_arms[s].sort(key=lambda a: abundance.get(a, -np.inf), reverse=True)

    rows = []
    for g in gse_arms:
        core = _VER_RE.sub("", g)             # strip _v16.0
        is_star = core.endswith("*")
        core = core.rstrip("*")
        # 1) already arm-suffixed and present in TCGA -> exact
        if _ARM_SUFFIX_RE.search(core) and core in tcga_set:
            rows.append((g, core, "exact_arm", is_star))
            continue
        # 2) stem grouping + abundance resolution
        members = stem_arms.get(_stem(core), [])
        if not members:
            continue                          # no modern arm for this hairpin -> drop
        if not is_star:
            rows.append((g, members[0], "abundance_guide", False))
        elif len(members) > 1:
            rows.append((g, members[1], "abundance_star", True))
        # star with only the guide measured in TCGA -> passenger absent -> drop
    out = pd.DataFrame(rows, columns=["gse_arm", "modern_arm", "bridge_method", "is_star"])
    # collapse GSE arms colliding on one modern arm: keep best method
    prio = {"exact_arm": 0, "abundance_guide": 1, "abundance_star": 2}
    out["_p"] = out["bridge_method"].map(prio)
    out = out.sort_values("_p").drop_duplicates("modern_arm", keep="first").drop(columns="_p")
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- #
# acquired healthy->tumour ABUNDANCE direction (GTEx-anchored)
# --------------------------------------------------------------------------- #
def acquired_direction(de: pd.DataFrame, state: Optional[pd.DataFrame]) -> pd.DataFrame:
    """Per modern arm: acquired abundance direction healthy->tumour (+1/-1/0).

    GTEx-detected -> sign(gtex_to_tumor_pct). GTEx-absent but DE class tumour_induced
    -> +1 (acquired from undetectable healthy, per the floor-0 healthy-anchor policy);
    tumour_suppressed -> -1. Otherwise 0 (not on the acquired axis). Pressure-axis
    ``dHT``/``acquired_gainer`` are merged as annotation.
    """
    d = de.rename(columns={"arm": "modern_arm"}).copy()
    det = d["gtex_detected"].fillna(False).astype(bool)
    dir_det = np.sign(d["gtex_to_tumor_pct"]).fillna(0.0)
    induced = d["de_class"].eq("tumor_induced")
    suppressed = d["de_class"].eq("tumor_suppressed")
    acq = np.where(det, dir_det,
                   np.where(induced, 1.0, np.where(suppressed, -1.0, 0.0)))
    d["acq_dir"] = acq.astype(int)
    d["acq_basis"] = np.where(det, "gtex_to_tumor_pct",
                              np.where(induced | suppressed, "tumor_de_class_floor0", "none"))
    keep = ["modern_arm", "acq_dir", "acq_basis", "gtex_to_tumor_pct",
            "log2fc_tumor_vs_nat", "trajectory", "de_class", "nat_field_effect",
            "tumor_median"]
    d = d[[c for c in keep if c in d.columns]]
    if state is not None:
        s = state.rename(columns={"miRNA": "modern_arm"})[
            ["modern_arm", "dHT", "acquired_gainer", "primary_class"]]
        d = d.merge(s, on="modern_arm", how="left")
    return d


# --------------------------------------------------------------------------- #
# in-situ -> invasive timing (DCIS vs IBC, same platform)
# --------------------------------------------------------------------------- #
def _mwu_rank_biserial(ibc: np.ndarray, dcis: np.ndarray) -> Tuple[float, float, float]:
    """Return (p, rank_biserial, median_delta). +rb / +delta => IBC > DCIS."""
    ibc = ibc[~np.isnan(ibc)]; dcis = dcis[~np.isnan(dcis)]
    if len(ibc) < MIN_PER_GROUP or len(dcis) < MIN_PER_GROUP:
        return np.nan, np.nan, np.nan
    u, p = mannwhitneyu(ibc, dcis, alternative="two-sided")
    auc = u / (len(ibc) * len(dcis))          # P(IBC > DCIS)
    return p, 2.0 * auc - 1.0, float(np.median(ibc) - np.median(dcis))


def insitu_invasive_timing(
    mirna: pd.DataFrame, meta: pd.DataFrame, bridge: pd.DataFrame, acq: pd.DataFrame,
    *, subset: Optional[pd.Index] = None, label: str = "all",
    cov: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """DCIS-vs-IBC per bridged arm, classified into EARLY/LATE/DIVERGENT timing.

    When ``cov`` (composition+prolif+site, from the paired mRNA) is supplied, also
    reports ``rank_biserial_adj`` = composition-adjusted partial Spearman of the arm vs
    the invasive-state indicator (the confound-controlled invasion movement).
    """
    cols = meta.index if subset is None else subset
    dcis_s = meta.index[(meta["state"] == "DCIS")].intersection(cols)
    ibc_s = meta.index[(meta["state"] == "IBC")].intersection(cols)
    state_bin = pd.Series(0.0, index=dcis_s).reindex(dcis_s.union(ibc_s)).fillna(1.0) \
        if cov is not None else None
    rows = []
    for _, br in bridge.iterrows():
        g = br["gse_arm"]
        if g not in mirna.index:
            continue
        p, rb, delta = _mwu_rank_biserial(
            mirna.loc[g, ibc_s].to_numpy(float), mirna.loc[g, dcis_s].to_numpy(float))
        row = {"gse_arm": g, "modern_arm": br["modern_arm"],
               "bridge_method": br["bridge_method"], "mwu_p": p,
               "rank_biserial": rb, "median_delta_ibc_minus_dcis": delta}
        if cov is not None:
            samp = state_bin.index
            rho_adj, _, _ = partial_spearman(mirna.loc[g, samp], state_bin, cov.reindex(samp))
            row["rank_biserial_adj"] = rho_adj   # +ve => higher in IBC, composition-adjusted
        rows.append(row)
    res = pd.DataFrame(rows)
    res = res.merge(acq, on="modern_arm", how="left")
    res["n_dcis"], res["n_ibc"] = len(dcis_s), len(ibc_s)
    # multiple-testing family = the acquired-axis arms (the only arms we time-classify)
    ax = res["acq_dir"].abs() > 0
    ok = res["mwu_p"].notna() & ax
    res["mwu_q"] = np.nan
    if ok.any():
        res.loc[ok, "mwu_q"] = bh_fdr(res.loc[ok, "mwu_p"].to_numpy())

    def _bh_class(r):
        if r.get("acq_dir", 0) == 0:
            return "not_on_acquired_axis"
        if not np.isfinite(r.get("mwu_q", np.nan)):
            return "untested"
        if r["mwu_q"] >= ALPHA:
            # NB: at n_ibc=14 this is "no detectable residual movement", i.e. EARLY
            # OR underpowered -- absence is NOT confirmed per-arm; see set-level test.
            return "EARLY_or_underpowered"
        return ("LATE_invasion_coupled"
                if np.sign(r["rank_biserial"]) == np.sign(r["acq_dir"])
                else "DIVERGENT_nonmonotonic")

    def _nominal(r):  # exploratory: uncorrected p<0.05, surfaces the weak movers
        if r.get("acq_dir", 0) == 0 or not np.isfinite(r.get("mwu_p", np.nan)):
            return ""
        if r["mwu_p"] >= 0.05:
            return "nominal_flat"
        return ("nominal_LATE" if np.sign(r["rank_biserial"]) == np.sign(r["acq_dir"])
                else "nominal_DIVERGENT")

    res["timing_class"] = res.apply(_bh_class, axis=1)
    res["nominal_timing"] = res.apply(_nominal, axis=1)
    res.insert(0, "stratum", label)
    return res.sort_values(["mwu_p"], na_position="last").reset_index(drop=True)


def set_level_aggregate(timing: pd.DataFrame, rhocol: str = "rank_biserial") -> dict:
    """Powered set-level test: does the acquired signature *as a whole* keep moving
    DCIS->IBC, or is it pre-invasively complete? Per-arm is underpowered at n_ibc=14;
    aggregating the signed effects is not.

    For acquired arms we orient each arm's effect (``rhocol``) by its acquired direction
    (``oriented = effect * acq_dir``); >0 means "still moving toward the tumour state at
    invasion". H0: oriented effects centred at 0 (Wilcoxon signed-rank). ``rhocol`` selects
    raw (``rank_biserial``) or composition-adjusted (``rank_biserial_adj``).
    """
    if rhocol not in timing.columns:
        return {}
    ax = timing[(timing["acq_dir"].abs() > 0) & timing[rhocol].notna()].copy()
    out: dict = {}
    for name, sub in (("gainers", ax[ax["acq_dir"] > 0]),
                      ("losers", ax[ax["acq_dir"] < 0])):
        if len(sub) < 6:
            out[name] = {"n": int(len(sub)), "note": "too few arms"}
            continue
        oriented = (sub[rhocol] * np.sign(sub["acq_dir"])).to_numpy()
        try:
            _, p = wilcoxon(oriented)
        except ValueError:
            p = np.nan
        med = float(np.median(oriented))
        out[name] = {
            "n": int(len(sub)),
            "median_oriented_rank_biserial": round(med, 4),
            "frac_still_moving_toward_tumour": round(float((oriented > 0).mean()), 3),
            "signed_rank_p": (float(f"{p:.3g}") if np.isfinite(p) else None),
            "interpretation": ("net further movement to invasion (loss continues)"
                               if np.isfinite(p) and p < ALPHA and med > 0
                               else "no net further movement -> consistent with "
                                    "pre-invasive establishment"),
        }
    return out


# --------------------------------------------------------------------------- #
# intra-DCIS grade gradient (severity inside in-situ disease; no healthy anchor needed)
# --------------------------------------------------------------------------- #
def _jonckheere(x: np.ndarray, g: np.ndarray) -> Tuple[float, float]:
    """Jonckheere-Terpstra ordered-trend test (normal approx). Returns (z, two-sided p)."""
    from scipy.stats import norm
    levels = np.unique(g)
    if levels.size < 2:
        return np.nan, np.nan
    J = 0.0
    for i in range(len(levels)):
        xi = x[g == levels[i]]
        for j in range(i + 1, len(levels)):
            xj = x[g == levels[j]]
            J += float(np.sum(xi[:, None] < xj[None, :]) + 0.5 * np.sum(xi[:, None] == xj[None, :]))
    n = len(x); ns = np.array([(g == l).sum() for l in levels], float)
    mu = (n ** 2 - np.sum(ns ** 2)) / 4.0
    var = (n ** 2 * (2 * n + 3) - np.sum(ns ** 2 * (2 * ns + 3))) / 72.0
    if var <= 0:
        return np.nan, np.nan
    z = (J - mu) / np.sqrt(var)
    return float(z), float(2 * norm.sf(abs(z)))


def grade_gradient(mirna: pd.DataFrame, meta: pd.DataFrame, bridge: pd.DataFrame,
                   acq: pd.DataFrame, *, cov: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Spearman(arm, DCIS grade 1/2/3) within DCIS + Jonckheere ordered test (+ adjusted)."""
    dcis = meta.index[meta["state"] == "DCIS"]
    grade = pd.to_numeric(meta.loc[dcis, "grade"], errors="coerce")
    keep = dcis[grade.notna()]
    gv = grade.loc[keep].to_numpy(float)
    rows = []
    if len(keep) >= MIN_PER_GROUP + 2:
        for _, br in bridge.iterrows():
            g = br["gse_arm"]
            if g not in mirna.index:
                continue
            x = mirna.loc[g, keep].to_numpy(float)
            m = ~np.isnan(x)
            if m.sum() < MIN_PER_GROUP + 2 or np.unique(gv[m]).size < 2:
                continue
            rho, p = spearmanr(x[m], gv[m])
            jz, jp = _jonckheere(x[m], gv[m])
            row = {"gse_arm": g, "modern_arm": br["modern_arm"], "grade_rho": rho,
                   "grade_p": p, "grade_jt_z": jz, "grade_jt_p": jp,
                   "n_dcis_graded": int(m.sum())}
            if cov is not None:
                kk = keep[m]
                rho_adj, p_adj, _ = partial_spearman(
                    mirna.loc[g, kk], grade.loc[kk], cov.reindex(kk))
                row["grade_rho_adj"] = rho_adj
                row["grade_p_adj"] = p_adj
            rows.append(row)
    res = pd.DataFrame(rows)
    if not res.empty:
        res["grade_q"] = bh_fdr(res["grade_p"].to_numpy())
        if "grade_p_adj" in res and res["grade_p_adj"].notna().any():
            ok = res["grade_p_adj"].notna()
            res.loc[ok, "grade_q_adj"] = bh_fdr(res.loc[ok, "grade_p_adj"].to_numpy())
        res = res.merge(acq[["modern_arm", "acq_dir", "trajectory"]], on="modern_arm", how="left")
    return res.sort_values("grade_p").reset_index(drop=True) if not res.empty else res


def grade_subtype_table(meta: pd.DataFrame) -> pd.DataFrame:
    """DCIS grade x PAM50 crosstab (is severity confounded with subtype?)."""
    dcis = meta[meta["state"] == "DCIS"]
    if "grade" not in dcis or "pam50" not in dcis:
        return pd.DataFrame()
    return pd.crosstab(dcis["grade"], dcis["pam50"]).reset_index()


def _summary(timing: pd.DataFrame, grade: pd.DataFrame, bridge: pd.DataFrame,
             meta: pd.DataFrame) -> dict:
    acq_axis = timing[timing["acq_dir"].abs() > 0]
    vc = acq_axis["timing_class"].value_counts().to_dict()
    nom = acq_axis["nominal_timing"].value_counts().to_dict()
    grade_sig = (int((grade["grade_q"] < ALPHA).sum()) if not grade.empty
                 and "grade_q" in grade else 0)
    grade_top = []
    if not grade.empty:
        for r in grade.nsmallest(5, "grade_p").itertuples():
            grade_top.append({"arm": r.modern_arm, "grade_rho": round(r.grade_rho, 3),
                              "grade_p": round(r.grade_p, 4),
                              "grade_jt_p": round(getattr(r, "grade_jt_p", np.nan), 4),
                              "grade_rho_adj": (round(r.grade_rho_adj, 3)
                                                if "grade_rho_adj" in grade
                                                and np.isfinite(getattr(r, "grade_rho_adj", np.nan)) else None)})
    # strongest nominal movers (the invasion-step signal that BH can't confirm at n=14)
    movers = []
    for r in acq_axis[acq_axis["nominal_timing"].isin(["nominal_LATE", "nominal_DIVERGENT"])
                      ].nsmallest(8, "mwu_p").itertuples():
        movers.append({"arm": r.modern_arm, "mwu_p": round(r.mwu_p, 4),
                       "rank_biserial": round(r.rank_biserial, 3),
                       "acq_dir": int(r.acq_dir), "nominal": r.nominal_timing,
                       "trajectory": getattr(r, "trajectory", None)})
    return {
        "n_gse_arms_bridged": int(len(bridge)),
        "bridge_methods": bridge["bridge_method"].value_counts().to_dict(),
        "n_on_acquired_axis_tested": int(acq_axis["mwu_q"].notna().sum()),
        "min_mwu_q": (round(float(acq_axis["mwu_q"].min()), 3)
                      if acq_axis["mwu_q"].notna().any() else None),
        "timing_counts_BH": vc,
        "timing_counts_nominal_p05": nom,
        "set_level": set_level_aggregate(timing),                       # raw
        "set_level_composition_adjusted": set_level_aggregate(timing, "rank_biserial_adj"),
        "top_nominal_invasion_movers": movers,
        "intra_dcis_grade_tracking_arms_q<alpha": grade_sig,
        "top_grade_gradient": grade_top,
        "n_dcis": int((meta["state"] == "DCIS").sum()),
        "n_ibc": int((meta["state"] == "IBC").sum()),
        "outcome_note": "recurrence/DMFS underpowered (3/2 events) -> not modelled",
        "caveats": [
            "n_ibc=14 underpowers per-arm BH; the powered statement is set_level (raw + composition-adjusted)",
            "set_level_composition_adjusted partials epi/immune/stroma+prolif+site from the paired GSE59246 mRNA",
            "cross-platform: TCGA used for acquired DIRECTION only (rank); magnitude is same-platform DCIS-vs-IBC",
            "no GPL15019 healthy anchor: EARLY = no residual in-situ->invasive movement under a monotone assumption, corroborated by the intra-DCIS grade gradient",
            "no reliable matched synchronous DCIS+IBC pairs in this cohort -> population-level contrasts only",
            "v16->modern arm bridge is best-effort; bridge_method recorded per arm",
        ],
        "alpha": ALPHA,
    }


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    print("[dcis_timing] loading GSE59247 (DCIS+IBC, GPL15019) + acquired signature ...")
    mirna, meta = load_gse59247()
    de = pd.read_csv(DE_LANDSCAPE, sep="\t")
    state = pd.read_csv(STATE_CLASS, sep="\t") if STATE_CLASS.exists() else None

    abundance = dict(zip(de["arm"], de["tumor_median"]))
    bridge = bridge_v16_to_modern(mirna.index, de["arm"], abundance)
    acq = acquired_direction(de, state)
    print(f"[dcis_timing] bridged {len(bridge)} GSE arms -> modern; "
          f"methods={bridge['bridge_method'].value_counts().to_dict()}")

    # composition+prolif+site covariates from the paired GSE59246 mRNA (MH-39 confound control)
    cov = None
    try:
        mrna, _ = load_gse59246()
        common = mirna.columns.intersection(mrna.columns)
        cov = composition_covariates(mrna[common], meta.loc[common])
        print(f"[dcis_timing] composition covariates from {len(common)} paired mRNA samples: "
              f"{list(cov.columns)}")
    except Exception as e:  # paired mRNA optional; raw still computed
        print(f"[dcis_timing] WARN: no paired mRNA covariates ({e}); raw only")

    timing = insitu_invasive_timing(mirna, meta, bridge, acq, cov=cov)
    timing.to_csv(out_dir / "dcis_insitu_invasive_timing.tsv", sep="\t", index=False)

    # PAM50-stratified (luminal vs non-luminal where n permits)
    grp = {"luminal": meta.index[meta["pam50"].isin(["LumA", "LumB"])],
           "nonluminal": meta.index[meta["pam50"].isin(["Her2", "Basal"])]}
    strat = []
    for lab, idx in grp.items():
        sub = insitu_invasive_timing(mirna, meta, bridge, acq, subset=idx, label=lab, cov=cov)
        if sub["mwu_q"].notna().any():
            strat.append(sub)
    if strat:
        pd.concat(strat, ignore_index=True).to_csv(
            out_dir / "dcis_timing_by_pam50.tsv", sep="\t", index=False)

    grade = grade_gradient(mirna, meta, bridge, acq, cov=cov)
    if not grade.empty:
        grade.to_csv(out_dir / "dcis_grade_gradient.tsv", sep="\t", index=False)
    gst = grade_subtype_table(meta)
    if not gst.empty:
        gst.to_csv(out_dir / "dcis_grade_subtype.tsv", sep="\t", index=False)

    summary = _summary(timing, grade, bridge, meta)
    manifest = {"module": "mirna_hallmark.analyses.dcis_ev.dcis_timing",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "cohort": "GSE59247 (DCIS+IBC, Agilent GPL15019)",
                "acquired_reference": str(DE_LANDSCAPE),
                **summary}
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"[dcis_timing] BH timing (acquired axis): {summary['timing_counts_BH']} "
          f"(min q={summary['min_mwu_q']})")
    print(f"[dcis_timing] SET-LEVEL raw      gainers: {summary['set_level'].get('gainers')}")
    print(f"[dcis_timing] SET-LEVEL raw      losers : {summary['set_level'].get('losers')}")
    if summary.get("set_level_composition_adjusted"):
        print(f"[dcis_timing] SET-LEVEL adjusted gainers: {summary['set_level_composition_adjusted'].get('gainers')}")
        print(f"[dcis_timing] SET-LEVEL adjusted losers : {summary['set_level_composition_adjusted'].get('losers')}")
    print(f"[dcis_timing] top nominal invasion movers: "
          f"{[m['arm'] for m in summary['top_nominal_invasion_movers']]}")
    print(f"[dcis_timing] intra-DCIS grade-tracking arms (q<{ALPHA}): "
          f"{summary['intra_dcis_grade_tracking_arms_q<alpha']}")
    print(f"[dcis_timing] wrote outputs under {out_dir}")
    return timing


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
