"""Cell-composition confounder PILOT for the prediction-only orphan candidates.

Follow-up to MH-36 / `cptac_orphan_discovery`. The strongest *discovery* tier there is the
**prediction-only** set: 24 edges that are TargetScan-sequence-only (no CLIP, no miRTarBase)
yet show a translational `protein_resid` signal in CPTAC-prospective AND replicate in
TCGA-105. Many of the targets (VCL, VCAN, TPM4, ...) are cytoskeletal / ECM-adjacent, so a
real worry is that the protein anti-correlation is driven by **tumor cell-composition**
(CAF / endothelial / immune content) rather than a miRNA→target edge.

This pilot re-tests those edges after adding composition confounders to the residualization:
PAM50 subtype + CAF / endothelial / immune **marker composite scores** (mean of cohort
z-scored RNA over marker genes; a lightweight deconvolution proxy), on top of the existing
base confounders (prospective: purity+CIN; tcga105: CPE+HRD). Run in **both** cohorts.

Survival = `protein_resid` stays negative & significant after adjustment. Drops identify
edges whose protein coupling is a composition artifact.

Marker genes overlapping a tested target gene are dropped from the marker score (anti-circular).

Outputs (`output/cptac_validation/orphan_discovery/confound_pilot/`):
  - `pilot_edge_confound_comparison.tsv`  per-edge baseline vs +composition rho/p/q, both cohorts
  - `pilot_marker_loadings.tsv`           target protein ~ each composition score (why an edge drops)
  - `method_manifest.json`

Run::

    .venv/bin/python3 -m mirna_hallmark.analyses.cptac.cptac_orphan_confound_pilot
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.analyses.cptac.cptac_orphan_discovery import (
    DISCOVERY_DIR,
    MIN_N,
    _edge_partial_spearman,
    _rank_matrix,
    _residualize,
    _spearman_p,
)
from mirna_hallmark.eval.cptac_validation import (
    get_cohort_config,
    load_cct,
    load_cptac_layers,
    load_mapped,
    load_prospective_clinical,
    load_prospective_mirna_arms,
    sid_participant_maps,
)
from mirna_hallmark.estimate_scores import compute_estimate

PILOT_DIR = DISCOVERY_DIR / "confound_pilot"

# Composition marker panels (HGNC symbols). Lightweight signatures, not full deconvolution.
CAF_MARKERS = ["FAP", "ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "COL3A1",
               "THY1", "S100A4", "PDPN", "POSTN", "TAGLN", "DCN", "LUM", "FN1", "VCAN"]
ENDO_MARKERS = ["PECAM1", "VWF", "CDH5", "KDR", "ENG", "CLDN5", "TEK", "FLT1",
                "EGFL7", "ESAM", "ROBO4", "TIE1"]
IMMUNE_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD8A", "CD2", "CD19", "MS4A1", "NKG7",
                  "GZMB", "GZMA", "CD68", "LYZ", "ITGAM", "CD14", "CD27", "CXCL9"]


# --------------------------------------------------------------------------- #
def _cohort_raw_rna(cohort: str) -> pd.DataFrame:
    """Raw (pre-z) RNA matrix gene×sample, columns keyed like the cov index, for ssGSEA."""
    if cohort == "tcga105":
        cfg = get_cohort_config("tcga105")
        mapped = load_mapped("tcga105")
        sid_to_part, _ = sid_participant_maps(mapped)
        rna = load_cct(cfg.rna_cct)
        rna = rna[[c for c in mapped["cptac_sample_id"] if c in rna.columns]]
        rna = rna.rename(columns=sid_to_part)
    else:
        cfg = get_cohort_config("pancan122")
        rna = load_cct(cfg.rna_cct)
    return rna.loc[:, ~rna.columns.duplicated()]


def _estimate_fractions(cohort: str, out_dir: Path, force: bool = False) -> pd.DataFrame:
    """ESTIMATE ssGSEA ImmuneScore + StromalScore per sample (z-scored, participant-keyed)."""
    rna = _cohort_raw_rna(cohort)
    scores = compute_estimate(f"cptac_{cohort}", rna, out_dir=out_dir, force=force)
    keep = scores[["ImmuneScore", "StromalScore"]].apply(pd.to_numeric, errors="coerce")
    return (keep - keep.mean()) / keep.std(ddof=0)  # z-scale so comparable to other covariates


def _marker_score(rna_z: pd.DataFrame, markers: Sequence[str], drop: set) -> pd.Series:
    """Mean of (already-z-scored) RNA over present marker genes, anti-circular drop."""
    use = [m for m in markers if m in rna_z.index and m not in drop]
    if not use:
        return pd.Series(dtype=float)
    return rna_z.loc[use].mean(axis=0)


def _pam50_dummies(clin: pd.DataFrame, samples: Sequence[str]) -> pd.DataFrame:
    """One-hot PAM50 (NaN -> 'NA' category), drop the most common level as reference."""
    s = clin["PAM50_final"].reindex(samples).fillna("NA").astype(str)
    ref = s.value_counts().index[0]
    d = pd.get_dummies(s, prefix="pam50").drop(columns=[f"pam50_{ref}"], errors="ignore")
    return d.astype(float)


def _build_cov(base_cov: pd.DataFrame, rna_z: pd.DataFrame, clin: pd.DataFrame,
               target_genes: set, fractions: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Augment base confounders with PAM50 dummies + composition controls.

    ``fractions`` (ESTIMATE Immune/Stromal, z-scored) replaces the marker immune+CAF
    composites when given (true ssGSEA cell estimations); the endothelial marker score
    is always kept since ESTIMATE folds endothelium into StromalScore.
    """
    samples = list(base_cov.index)
    comp = pd.DataFrame(index=samples)
    comp["endo_score"] = _marker_score(rna_z, ENDO_MARKERS, target_genes).reindex(samples)
    if fractions is not None:
        comp["immune_estimate"] = fractions["ImmuneScore"].reindex(samples)
        comp["stromal_estimate"] = fractions["StromalScore"].reindex(samples)
    else:
        comp["caf_score"] = _marker_score(rna_z, CAF_MARKERS, target_genes).reindex(samples)
        comp["immune_score"] = _marker_score(rna_z, IMMUNE_MARKERS, target_genes).reindex(samples)
    pam = _pam50_dummies(clin, samples)
    return pd.concat([base_cov, comp, pam], axis=1)


def _edge_stats(edges: pd.DataFrame, arms: pd.DataFrame, layers: Dict[str, pd.DataFrame],
                cov: pd.DataFrame, tag: str) -> pd.DataFrame:
    """Partial-Spearman (protein_z & protein_resid) for the pilot edges under one cov set."""
    arm_rank = _rank_matrix(_residualize(arms, cov))
    prot_rank = _rank_matrix(_residualize(layers["protein_z"], cov))
    resid_rank = _rank_matrix(_residualize(layers["protein_resid"], cov))
    rows = []
    for mir, gene in zip(edges["miRNA"], edges["gene"]):
        if mir not in arm_rank.index or gene not in prot_rank.index:
            rows.append({"miRNA": mir, "gene": gene, f"{tag}_n": np.nan})
            continue
        pr, n = _edge_partial_spearman(arm_rank.loc[mir], prot_rank.loc[gene])
        rr, _ = _edge_partial_spearman(arm_rank.loc[mir], resid_rank.loc[gene])
        rows.append({"miRNA": mir, "gene": gene, f"{tag}_n": n,
                     f"{tag}_protein_rho": pr, f"{tag}_protein_p": _spearman_p(pr, n),
                     f"{tag}_resid_rho": rr, f"{tag}_resid_p": _spearman_p(rr, n)})
    return pd.DataFrame(rows)


def _marker_loadings(edges: pd.DataFrame, layers: Dict[str, pd.DataFrame],
                     cov_full: pd.DataFrame, cohort: str) -> pd.DataFrame:
    """For each target gene: Spearman(protein_z, composition score) — explains drops."""
    prot = layers["protein_z"]
    comp_cols = [c for c in ("caf_score", "endo_score", "immune_score",
                             "immune_estimate", "stromal_estimate") if c in cov_full.columns]
    out = []
    for gene in sorted(set(edges["gene"])):
        if gene not in prot.index:
            continue
        g = prot.loc[gene]
        rec = {"cohort": cohort, "gene": gene}
        for c in comp_cols:
            df = pd.concat([g.rename("g"), cov_full[c].rename("c")], axis=1).dropna()
            if len(df) >= MIN_N:
                rec[c] = round(float(np.corrcoef(df["g"].rank(), df["c"].rank())[0, 1]), 3)
        out.append(rec)
    return pd.DataFrame(out)


def _run_cohort(cohort: str, edges: pd.DataFrame, arms: pd.DataFrame,
                layers: Dict[str, pd.DataFrame], clin: pd.DataFrame,
                base_cov: pd.DataFrame,
                fractions: Optional[pd.DataFrame] = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_genes = set(edges["gene"])
    base = _edge_stats(edges, arms, layers, base_cov, f"{cohort}_base")
    cov_full = _build_cov(base_cov, layers["rna_z"], clin, target_genes, fractions=fractions)
    aug = _edge_stats(edges, arms, layers, cov_full, f"{cohort}_comp")
    merged = base.merge(aug, on=["miRNA", "gene"], how="outer")
    # BH-q within the pilot set, on the resid p-values (the translational layer).
    for tag in (f"{cohort}_base", f"{cohort}_comp"):
        p = merged[f"{tag}_resid_p"]
        merged[f"{tag}_resid_q"] = S.bh_fdr(p.fillna(1.0).values)
    loadings = _marker_loadings(edges, layers, cov_full, cohort)
    return merged, loadings


# --------------------------------------------------------------------------- #
def run(*, out_dir: Path = PILOT_DIR, fdr: float = 0.1, use_estimate: bool = False,
        force_estimate: bool = False, batch_kind: str = "none") -> None:
    if batch_kind != "none":
        out_dir = out_dir / f"batch_{batch_kind}"
    out_dir.mkdir(parents=True, exist_ok=True)
    from mirna_hallmark import cptac_batch as B
    sfx = "_estimate" if use_estimate else ""
    mode = "ESTIMATE ssGSEA fractions (Immune+Stromal) + endo marker + PAM50" if use_estimate \
        else "CAF/endo/immune marker composites + PAM50"
    print(f"[pilot] composition mode: {mode}")

    cand = pd.read_csv(DISCOVERY_DIR / "orphan_candidates.tsv", sep="\t")
    edges = cand[(cand["prior"] == "targetscan_only")
                 & cand["translational_candidate"].astype(bool)
                 & cand["tcga_replicated"].astype(bool)][["miRNA", "gene", "ts_weight"]].copy()
    print(f"[pilot] prediction-only replicated edges: {len(edges)}")

    pros_frac = _estimate_fractions("prospective", out_dir, force=force_estimate) if use_estimate else None
    tcga_frac = _estimate_fractions("tcga105", out_dir, force=force_estimate) if use_estimate else None

    # Prospective (primary).
    pros_arms = load_prospective_mirna_arms()
    pros_layers = load_cptac_layers("prospective")
    pros_clin = load_prospective_clinical().drop_duplicates("participant").set_index("participant")
    pros_base, pbc = B.augment_cov(pros_clin[["purity", "cin"]], "prospective", batch_kind)
    if pbc:
        print(f"[pilot] +batch ({batch_kind}) prospective base: {len(pbc)} dummies")
    pros, pros_load = _run_cohort("pros", edges, pros_arms, pros_layers, pros_clin, pros_base, pros_frac)

    # TCGA-105 (replication).
    tcga_arms = D.load_mirna_arms()
    tcga_layers = load_cptac_layers("tcga105")
    tcga_clin = D.load_clinical_strata().drop_duplicates("participant").set_index("participant")
    tcga_base, tbc = B.augment_cov(
        tcga_clin[[c for c in C.CONFOUNDER_NUMERIC if c in tcga_clin.columns]], "tcga105", batch_kind)
    if tbc:
        print(f"[pilot] +batch ({batch_kind}) tcga105 base: {len(tbc)} dummies")
    tcga, tcga_load = _run_cohort("tcga", edges, tcga_arms, tcga_layers, tcga_clin, tcga_base, tcga_frac)

    comp = edges.merge(pros, on=["miRNA", "gene"], how="left").merge(tcga, on=["miRNA", "gene"], how="left")

    # Survival flags: resid neg & q<fdr, base vs composition-adjusted.
    for coh in ("pros", "tcga"):
        comp[f"{coh}_base_surv"] = (comp[f"{coh}_base_resid_rho"] < 0) & (comp[f"{coh}_base_resid_q"] < fdr)
        comp[f"{coh}_comp_surv"] = (comp[f"{coh}_comp_resid_rho"] < 0) & (comp[f"{coh}_comp_resid_q"] < fdr)
    comp["pros_dropped"] = comp["pros_base_surv"] & ~comp["pros_comp_surv"]
    comp = comp.sort_values("pros_comp_resid_rho")

    comp.to_csv(out_dir / f"pilot_edge_confound_comparison{sfx}.tsv", sep="\t", index=False)
    pd.concat([pros_load, tcga_load], ignore_index=True).to_csv(
        out_dir / f"pilot_marker_loadings{sfx}.tsv", sep="\t", index=False)

    n_pros_base = int(comp["pros_base_surv"].sum())
    n_pros_comp = int(comp["pros_comp_surv"].sum())
    n_tcga_neg_base = int((comp["tcga_base_resid_rho"] < 0).sum())
    n_tcga_neg_comp = int((comp["tcga_comp_resid_rho"] < 0).sum())
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "module": "mirna_hallmark.analyses.cptac.cptac_orphan_confound_pilot",
        "pilot_set": "prediction-only (targetscan_only) translational + TCGA-replicated orphan edges",
        "composition_mode": "estimate_fractions" if use_estimate else "marker_composites",
        "batch_kind": batch_kind,
        "n_edges": int(len(edges)),
        "added_confounders": (["PAM50 subtype (dummies)", "ESTIMATE ImmuneScore (ssGSEA)",
                               "ESTIMATE StromalScore (ssGSEA; folds CAF+endothelium)", "endothelial marker score"]
                              if use_estimate else
                              ["PAM50 subtype (dummies)", "CAF marker score", "endothelial marker score", "immune marker score"]),
        "base_confounders": {"prospective": ["purity", "cin"], "tcga105": list(C.CONFOUNDER_NUMERIC)},
        "marker_panels": {"CAF": CAF_MARKERS, "endothelial": ENDO_MARKERS, "immune": IMMUNE_MARKERS},
        "fdr": fdr,
        "prospective_resid_survivors_base": n_pros_base,
        "prospective_resid_survivors_with_composition": n_pros_comp,
        "tcga105_resid_negative_base": n_tcga_neg_base,
        "tcga105_resid_negative_with_composition": n_tcga_neg_comp,
        "caveat": "Marker composite scores are a deconvolution proxy; PAM50 NaN folded as 'NA' level.",
        "outputs": ["pilot_edge_confound_comparison.tsv", "pilot_marker_loadings.tsv"],
    }
    (out_dir / f"method_manifest{sfx}.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"[pilot] prospective resid survivors: base {n_pros_base} -> +composition {n_pros_comp}")
    print(f"[pilot] tcga105 resid negative: base {n_tcga_neg_base} -> +composition {n_tcga_neg_comp}")
    cols = ["miRNA", "gene", "pros_base_resid_rho", "pros_base_resid_q",
            "pros_comp_resid_rho", "pros_comp_resid_q", "pros_comp_surv",
            "tcga_comp_resid_rho"]
    print(comp[cols].to_string(index=False))
    print(f"[pilot] wrote {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=PILOT_DIR)
    ap.add_argument("--fdr", type=float, default=0.1)
    ap.add_argument("--use-estimate", action="store_true",
                    help="Use ESTIMATE ssGSEA Immune+Stromal fractions instead of marker composites")
    ap.add_argument("--force-estimate", action="store_true", help="Recompute the ESTIMATE caches")
    ap.add_argument("--batch-kind", default="none", choices=["none", "site", "plex", "auto"],
                    help="Add MS-plex/site batch dummies to the base cov (writes to batch_<kind>/ subdir)")
    args = ap.parse_args()
    run(out_dir=args.out_dir, fdr=args.fdr, use_estimate=args.use_estimate,
        force_estimate=args.force_estimate, batch_kind=args.batch_kind)


if __name__ == "__main__":
    main()
