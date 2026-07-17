"""Target specificity of the CPTAC protein-residual couplings: cognate vs unrelated proteins.

Follow-up to MH-35 (`cptac_resid_survivors`). The survivor drivers establish that each focal
arm's *own target* protein sits below its mRNA-predicted level as the arm rises
(``protein_resid`` coupling **negative**). The open worry is **specificity**: if an arm were
merely a proxy for a global state (proliferation, tumor purity, a cell-composition axis) it
would anti-correlate with *every* protein, not just its target. A genuine miRNA→target edge
must instead be **negative for the cognate target but ~0 for unrelated proteins**.

This module runs that decoy / placebo control, per driver arm m:

  1. **Cognate**: partial Spearman(arm m, ``protein_resid`` of its target g), cov-residualized
     (prospective: purity+CIN) — re-derives the survivor coupling (expected negative).
  2. **Decoy pool**: the same partial Spearman of arm m against every Hallmark-universe protein
     measured in CPTAC that m does **not** target by *any* prior — miRTarBase (all evidence),
     TargetScan (seed, ts_weight≥MIN), or ENCORI (CLIP). These are the "unrelated proteins".
  3. **Specificity stats**: where the cognate rho sits in the decoy distribution
     (empirical left-tail p = frac(decoy ≤ cognate); z = (cognate − decoy_mean)/decoy_sd), and
     whether the decoy distribution itself is **broadly negative** (Wilcoxon of decoy rhos vs 0;
     a centered-on-0 decoy median is the "not broadly negative" result that rules out a global
     confound).

Cohort: **prospective** (CPTAC-2), the cohort where the protein-residual signal exists; loaders
reused from ``cptac_validation``. Decoys are drawn from the Hallmark universe (the edge/pressure
universe of the subproject), keeping the test inside the same gene space.

Outputs (``output/cptac_validation/prospective/target_specificity/``):
  - ``arm_specificity.tsv``       one row per focal (arm, target): cognate rho + decoy summary +
                                  empirical p / z + decoy-vs-0 Wilcoxon (the headline table).
  - ``decoy_couplings.tsv.gz``    every (arm, decoy protein) partial Spearman (the raw null).
  - ``method_manifest.json``

Run::

    .venv/bin/python3 -m mirna_hallmark.cptac_target_specificity
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from mirna_hallmark import config as C
from mirna_hallmark import cptac_batch as B
from mirna_hallmark import data_loaders as D
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.analyses.cptac.cptac_orphan_discovery import (
    MIN_N,
    _edge_partial_spearman,
    _rank_matrix,
    _residualize,
    _spearman_p,
    build_orphan_edges,
)
from mirna_hallmark.eval.cptac_validation import (
    OUT_DIR,
    load_cptac_layers,
    load_prospective_clinical,
    load_prospective_mirna_arms,
)

SPEC_DIR = OUT_DIR / "prospective" / "target_specificity"
SURVIVOR_DRIVERS = OUT_DIR / "prospective" / "resid_survivors" / "survivor_drivers.tsv"
MIN_DECOYS = 50  # need a populated null to call specificity

# "Batch" for CPTAC is the MS plex (literal, where the design is known) or the enrolling-site
# proxy; both live in mirna_hallmark.cptac_batch. For the prospective cohort the literal TMT
# plex (PDC000120) confounds the protein side; site is the accrual-level fallback shared by both
# assays. The batch pass adds those dummies to the purity+CIN residualization.


# --------------------------------------------------------------------------- #
# Focal edges + per-arm target (exclusion) sets
# --------------------------------------------------------------------------- #
def focal_edges() -> pd.DataFrame:
    """Survivor driver edges = (driver arm, cognate target) with established resid coupling."""
    if not SURVIVOR_DRIVERS.exists():
        raise FileNotFoundError(
            f"Missing {SURVIVOR_DRIVERS}. Run `mirna_hallmark.cptac_resid_survivors` first."
        )
    d = pd.read_csv(SURVIVOR_DRIVERS, sep="\t")
    return d.rename(columns={"gene_name": "target_gene", "driver_miRNA": "miRNA"})


def arm_target_sets(universe: Sequence[str]) -> Dict[str, Set[str]]:
    """For each arm: the set of genes it targets by ANY prior (the decoy-exclusion set).

    Union of miRTarBase (all-evidence hallmark edges) + TargetScan seed + ENCORI CLIP orphan
    priors. Anything outside this set (within the universe) is an *unrelated* protein for the arm.
    """
    targets: Dict[str, Set[str]] = {}

    mirtar = D.load_hallmark_edges()  # miRTarBase-derived (all evidence), arm-resolved
    for arm, sub in mirtar.groupby(mirtar["miRNA"].astype(str)):
        targets.setdefault(arm, set()).update(sub["gene"].astype(str))

    orphan = build_orphan_edges(list(universe))  # TargetScan seed ∪ ENCORI CLIP priors
    for arm, sub in orphan.groupby(orphan["miRNA"].astype(str)):
        targets.setdefault(arm, set()).update(sub["gene"].astype(str))

    return targets


# --------------------------------------------------------------------------- #
# Specificity engine
# --------------------------------------------------------------------------- #
def _arm_vs_pool(
    arm_rank: pd.Series, resid_rank: pd.DataFrame, genes: Sequence[str]
) -> pd.DataFrame:
    """Partial Spearman of one arm rank-row vs each gene's protein_resid rank-row."""
    rows: List[dict] = []
    for g in genes:
        if g not in resid_rank.index:
            continue
        rho, n = _edge_partial_spearman(arm_rank, resid_rank.loc[g])
        if n < MIN_N or not np.isfinite(rho):
            continue
        rows.append({"gene": g, "rho": rho, "n": n, "p": _spearman_p(rho, n)})
    return pd.DataFrame(rows)


def run_specificity(
    *,
    focal: pd.DataFrame,
    arm_rank: pd.DataFrame,
    resid_rank: pd.DataFrame,
    targets: Dict[str, Set[str]],
    universe: Set[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (per-focal-edge specificity table, long decoy-coupling table)."""
    measured = set(resid_rank.index) & universe
    spec_rows: List[dict] = []
    decoy_long: List[pd.DataFrame] = []

    for _, fr in focal.iterrows():
        arm, tgt = str(fr["miRNA"]), str(fr["target_gene"])
        if arm not in arm_rank.index or tgt not in resid_rank.index:
            continue

        a = arm_rank.loc[arm]
        cog_rho, cog_n = _edge_partial_spearman(a, resid_rank.loc[tgt])
        if not np.isfinite(cog_rho) or cog_n < MIN_N:
            continue  # no usable cognate estimate (sparse protein_resid) -> can't assess specificity

        # Decoys = measured universe proteins the arm does NOT target by any prior.
        excl = targets.get(arm, set()) | {tgt}
        decoy_genes = sorted(measured - excl)
        dec = _arm_vs_pool(a, resid_rank, decoy_genes)
        if len(dec) < MIN_DECOYS:
            continue
        dec.insert(0, "miRNA", arm)
        dec.insert(1, "cognate_target", tgt)
        decoy_long.append(dec)

        rhos = dec["rho"].to_numpy(dtype=float)
        sd = float(np.std(rhos, ddof=1))
        mean = float(np.mean(rhos))
        median = float(np.median(rhos))
        # Empirical left-tail p: how often an unrelated protein is at least as negative.
        emp_p = float((rhos <= cog_rho).mean())
        z = (cog_rho - mean) / sd if sd > 0 else np.nan
        # Is the decoy distribution itself broadly negative? (centered-on-0 => specific)
        try:
            w_stat, w_p = wilcoxon(rhos)
        except ValueError:
            w_stat, w_p = np.nan, np.nan
        # Mann-Whitney: is the single cognate more negative than the decoy bulk? (descriptive)
        spec_rows.append({
            "miRNA": arm,
            "cognate_target": tgt,
            "cognate_resid_rho": _r(cog_rho),
            "cognate_n": int(cog_n),
            "driver_fdr_sig": bool(fr.get("driver_fdr_sig", False)),
            "literature_status": fr.get("literature_status", ""),
            "n_decoys": int(len(dec)),
            "decoy_mean_rho": _r(mean),
            "decoy_median_rho": _r(median),
            "decoy_sd_rho": _r(sd),
            "decoy_frac_negative": _r(float((rhos < 0).mean())),
            "decoy_frac_neg_p05": _r(float(((rhos < 0) & (dec["p"] < 0.05)).mean())),
            "decoy_vs0_wilcoxon_p": _r(float(w_p)),
            "empirical_left_p": _r(emp_p),     # frac decoys ≤ cognate (small = specific)
            "specificity_z": _r(z),            # SDs the cognate sits below the decoy mean
            "specific": bool(np.isfinite(emp_p) and emp_p < 0.05 and cog_rho < 0),
        })

    spec = pd.DataFrame(spec_rows).sort_values("empirical_left_p")
    decoys = pd.concat(decoy_long, ignore_index=True) if decoy_long else pd.DataFrame()
    return spec, decoys


def _r(v) -> float:
    return round(float(v), 4) if v is not None and np.isfinite(v) else np.nan


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def _pass(focal, arms, resid_layer, cov, targets, universe):
    """Residualize+rank on a cov set, then run the specificity engine."""
    arm_rank = _rank_matrix(_residualize(arms, cov))
    resid_rank = _rank_matrix(_residualize(resid_layer, cov))
    return run_specificity(
        focal=focal, arm_rank=arm_rank, resid_rank=resid_rank,
        targets=targets, universe=universe,
    )


def run(*, out_dir: Path = SPEC_DIR, batch_kind: str = "plex") -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    universe = set(HallmarkSets.load().universe)

    focal = focal_edges()
    print(f"[specificity] {len(focal)} focal survivor-driver edges")

    arms = load_prospective_mirna_arms()
    layers = load_cptac_layers("prospective")
    clin = load_prospective_clinical().drop_duplicates("participant").set_index("participant")
    base_cov = clin[["purity", "cin"]]
    targets = arm_target_sets(universe)

    # Base pass: purity+CIN (matching the survivor/orphan screens).
    spec, decoys = _pass(focal, arms, layers["protein_resid"], base_cov, targets, universe)
    if spec.empty:
        print("[specificity] no focal edge had enough decoys / measured target; nothing written.")
        return

    # Batch pass: + MS-plex (or site) dummies on top of purity+CIN.
    batch_cov, bcols = B.augment_cov(base_cov, "prospective", batch_kind)
    batch_info: Dict[str, object] = {"enabled": bool(bcols), "kind": batch_kind}
    if bcols:
        bspec, _ = _pass(focal, arms, layers["protein_resid"], batch_cov, targets, universe)
        keep = ["miRNA", "cognate_target", "cognate_resid_rho", "n_decoys",
                "decoy_median_rho", "empirical_left_p", "specificity_z", "specific"]
        bspec = bspec[keep].rename(columns={c: f"{c}_batch" for c in keep[2:]})
        spec = spec.merge(bspec, on=["miRNA", "cognate_target"], how="left")
        spec["batch_robust"] = spec["specific"] & spec["specific_batch"].fillna(False)
        batch_info.update({
            "definition": f"{batch_kind} dummies (prospective): literal TMT plex (PDC000120) when "
                          "kind=plex/auto, else enrolling-site/TSS proxy; reference level dropped",
            "n_batch_dummies": len(bcols),
            "added_to_cov": ["purity", "cin"] + bcols,
            "rationale": "TMT plex is the canonical proteomics batch (confounds protein_resid); "
                         "site is the accrual-level fallback shared by both assays",
            "n_specific_base": int(spec["specific"].sum()),
            "n_specific_batch": int(spec["specific_batch"].fillna(False).sum()),
            "n_batch_robust": int(spec["batch_robust"].sum()),
            "n_focal_fdr_batch_robust": int((spec["batch_robust"] & spec["driver_fdr_sig"]).sum()),
        })

    spec.to_csv(out_dir / "arm_specificity.tsv", sep="\t", index=False)
    decoys.to_csv(out_dir / "decoy_couplings.tsv.gz", sep="\t", index=False, compression="gzip")

    # Pooled read: cognate bulk vs decoy bulk, and the grand decoy distribution vs 0.
    cog_rhos = spec["cognate_resid_rho"].dropna().to_numpy(dtype=float)
    all_decoy = decoys["rho"].to_numpy(dtype=float)
    try:
        mw_u, mw_p = mannwhitneyu(cog_rhos, all_decoy, alternative="less")
    except ValueError:
        mw_u, mw_p = np.nan, np.nan
    grand_decoy_median = float(np.median(all_decoy))
    cognate_median = float(np.median(cog_rhos)) if len(cog_rhos) else np.nan
    try:
        _, grand_decoy_p = wilcoxon(all_decoy)
    except ValueError:
        grand_decoy_p = np.nan

    n_specific = int(spec["specific"].sum())
    n_sig_focal = int(spec["driver_fdr_sig"].sum())
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "module": "mirna_hallmark.cptac_target_specificity",
        "cohort": "prospective (CPTAC-2)",
        "question": "cognate protein_resid coupling negative AND unrelated proteins not broadly negative",
        "focal_set": "cptac_resid_survivors survivor_drivers (driver arm -> cognate target)",
        "decoy_definition": "Hallmark-universe proteins measured in CPTAC NOT targeted by the arm "
                            "(miRTarBase any-evidence ∪ TargetScan seed ∪ ENCORI CLIP)",
        "cov": "purity + CIN (partial Spearman, cov-residualized ranks)",
        "min_decoys": MIN_DECOYS,
        "n_focal_edges_tested": int(len(spec)),
        "n_focal_fdr_driver": n_sig_focal,
        "n_specific_emp_p_lt_0p05": n_specific,
        "cognate_vs_decoy_mannwhitney_less_p": _r(float(mw_p)),
        "cognate_median_rho": _r(cognate_median),
        "grand_decoy_median_rho": _r(grand_decoy_median),
        "grand_decoy_vs0_wilcoxon_p": _r(float(grand_decoy_p)),
        "wilcoxon_note": "decoy-vs-0 Wilcoxon is over-powered (thousands of decoys); the "
                         "meaningful read is the MAGNITUDE gap: cognate_median_rho vs "
                         "grand_decoy_median_rho (decoys sit ~0, ~an order of magnitude weaker).",
        "interpretation": (
            "specific = cognate sits in the decoy left tail (emp p<0.05) while the decoy "
            "median is ~0 (grand_decoy_median_rho); a strongly-negative decoy median would "
            "instead indicate a global/composition confound rather than a target-specific edge."
        ),
        "caveat": "Seed-family sharing only partly removed via TargetScan exclusion; protein "
                  "anti-correlation != direct edge.",
        "batch_pilot": batch_info,
        "outputs": ["arm_specificity.tsv", "decoy_couplings.tsv.gz"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"[specificity] wrote {out_dir}")
    print(f"  focal edges tested: {len(spec)} ({n_sig_focal} FDR drivers) | "
          f"specific (emp p<0.05): {n_specific}")
    print(f"  cognate vs decoy (MWU less) p={_r(float(mw_p))} | "
          f"grand decoy median ρ={_r(grand_decoy_median)} (vs0 p={_r(float(grand_decoy_p))})")
    if batch_info.get("enabled"):
        print(f"  [batch:{batch_info['kind']}] +{batch_info['n_batch_dummies']} dummies: "
              f"specific base {batch_info['n_specific_base']} -> +batch "
              f"{batch_info['n_specific_batch']} | batch-robust {batch_info['n_batch_robust']} "
              f"({batch_info['n_focal_fdr_batch_robust']} of FDR drivers)")
    cols = ["miRNA", "cognate_target", "cognate_resid_rho", "empirical_left_p", "specific"]
    if "specific_batch" in spec.columns:
        cols += ["cognate_resid_rho_batch", "empirical_left_p_batch", "specific_batch", "batch_robust"]
    print(spec[cols].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=SPEC_DIR)
    ap.add_argument("--batch-kind", default="plex", choices=["none", "site", "plex", "auto"],
                    help="Batch confounder added on top of purity+CIN (default: literal TMT plex)")
    ap.add_argument("--min-purity", type=float, default=None, help="Accepted for runner compatibility (no-op)")
    args = ap.parse_args()
    run(out_dir=args.out_dir, batch_kind=args.batch_kind)


if __name__ == "__main__":
    main()
