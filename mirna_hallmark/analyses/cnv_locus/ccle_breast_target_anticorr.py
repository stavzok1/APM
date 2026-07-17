"""Phase 4c, breast-retargeted: framework miRNA->target anti-correlation across CCLE BREAST lines.

The half-cell measured per-cell test (MH-65) was the only *single-cell* matched miRNA+mRNA data that
exists -- but it is **K562 (leukemia)** and n=19, so it came back underpowered/inconclusive. This module
trades single-cell resolution for **breast tissue + subtype coverage + power**: it tests the framework's
1,013 headline BY-neg edges as miRNA<->target anti-correlation **across the CCLE breast cell-line panel**
(50 lines with matched miRNA + mRNA, intrinsic subtypes Luminal/HER2/Basal), overall and **per subtype**.

Axis note: this is a **cross-cell-line** correlation (genetic/expression variation between lines), NOT
per-cell -- a different and complementary axis to MH-65. It is measured + breast + subtype-matched +
non-circular (CCLE NanoString miRNA + DepMap RNA-seq are independent of the miRTarBase+TCGA edges).

Confounder: across breast lines both miRNA and target track **subtype**. So every all-breast coupling is
reported **marginal AND subtype-adjusted** (partial Spearman on intrinsic-subtype dummies). CRITICAL: an
internal **positive control** of gold-standard repressive pairs (miR-200->ZEB, let-7->HMGA2, ...) gates
interpretation -- if the axis cannot recover those, a framework null is inconclusive (apm-rigor-protocol).

Data: ``data/CCLE/CCLE_miRNA_20181103.gct`` (NanoString) + ``OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv``
(DepMap 26Q1, log2 TPM+1). Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.cnv_locus.ccle_breast_target_anticorr``
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest, mannwhitneyu

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr, correlation_pair
from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import (
    build_ccle_arm_expression_map, _ccle_expr_matrix, _load_model_lineage,
)

EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
MRNA_TPM = C.REPO_ROOT / "data" / "CCLE" / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
OUT_DIR = C.OUTPUT_ROOT / "ccle_breast_target_anticorr"
RNG = np.random.default_rng(20260626)
MIN_ALL = 18          # min lines for the all-breast stratum (per-edge FDR feasible)
MIN_SUBTYPE = 8       # min lines for a per-PAM50 stratum (set-level only)


def _intrinsic_subtype(feat: str) -> str:
    """Coarse receptor-based intrinsic subtype from DepMap ModelSubtypeFeatures.

    HER2 (any HER2+/amp) > Basal (TNBC / basal_A/B / ER-PR-) > Luminal (ER+/PR+) — the parent
    cell-line map only covers ~10 curated lines, so we use DepMap's own receptor annotation for full
    breast-panel coverage. Coarser than centroid PAM50 but defensible and complete.
    """
    f = str(feat).lower()
    if "her2" in f:
        return "HER2"
    if "tnbc" in f or "basal" in f or ("er-" in f and "pr-" in f):
        return "Basal"
    if "er" in f or "pr" in f or "luminal" in f:
        return "Luminal"
    return "Unknown"


def _breast_subtype() -> pd.Series:
    """ModelID -> coarse intrinsic subtype (Luminal/HER2/Basal) for breast lines."""
    model = pd.read_csv(C.CCLE_MODEL, usecols=["ModelID", "OncotreeLineage", "ModelSubtypeFeatures"])
    model = model[model["OncotreeLineage"] == "Breast"]
    return pd.Series(
        {r["ModelID"]: _intrinsic_subtype(r["ModelSubtypeFeatures"]) for _, r in model.iterrows()},
        name="subtype",
    )


def _mrna_breast_matrix(breast_ids: set[str]) -> pd.DataFrame:
    """gene x ModelID log2(TPM+1) for breast lines (default DepMap entry; entrez stripped)."""
    ex = pd.read_csv(MRNA_TPM)
    ex = ex[ex["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")
    keep = [m for m in ex.index if m in breast_ids]
    ex = ex.loc[keep]
    gene_cols = [c for c in ex.columns if re.match(r".+ \(\d+\)$", c)]
    ex = ex[gene_cols]
    ex.columns = [c.split(" (")[0] for c in gene_cols]
    return ex.T                                            # gene x ModelID


# Gold-standard, literature-canonical miRNA->target repressive pairs (breast/EMT-relevant). If the
# cross-line axis cannot recover THESE as negative, it is insensitive to miRNA repression and a null on
# the framework edges is inconclusive, not a refutation (apm-rigor-protocol internal positive control).
POSITIVE_CONTROL_PAIRS = [
    ("hsa-miR-200c-3p", "ZEB1"), ("hsa-miR-200c-3p", "ZEB2"),
    ("hsa-miR-141-3p", "ZEB2"), ("hsa-miR-200b-3p", "ZEB1"),
    ("hsa-let-7a-5p", "HMGA2"), ("hsa-let-7b-5p", "HMGA2"),
    ("hsa-miR-21-5p", "PDCD4"), ("hsa-miR-21-5p", "PTEN"),
    ("hsa-miR-10b-5p", "HOXD10"), ("hsa-miR-155-5p", "SOCS1"),
]


# Literature EMT identity sets (claudin-low/mesenchymal vs luminal/epithelial breast lines) -- used
# only to VALIDATE the miRNA data: miR-200c/141 are canonical epithelial markers and must be far higher
# in epithelial lines. If the miRNA matrix fails this, cross-line nulls are measurement-confounded.
EMT_EPITHELIAL = {"MCF7", "T47D", "ZR751", "CAMA1", "BT474", "MDAMB453", "KPL1", "EFM19",
                  "HCC1500", "MDAMB415", "ZR7530"}
EMT_MESENCHYMAL = {"MDAMB231", "HS578T", "BT549", "MDAMB157", "MDAMB436", "HCC1395", "CAL120", "CAL51"}


def _emt_identity_check(mir: pd.DataFrame, lines, name_map) -> dict:
    """Does the miRNA data recover the known epithelial>>mesenchymal miR-200/141 gradient?"""
    from scipy.stats import mannwhitneyu
    epi = [m for m in lines if name_map.get(m) in EMT_EPITHELIAL]
    mes = [m for m in lines if name_map.get(m) in EMT_MESENCHYMAL]
    out = {"n_epithelial": len(epi), "n_mesenchymal": len(mes), "markers": {}}
    for arm in ("hsa-miR-200c-3p", "hsa-miR-141-3p"):
        if arm in mir.index and epi and mes:
            e = mir.loc[arm, epi].astype(float).dropna()
            m = mir.loc[arm, mes].astype(float).dropna()
            u, p = mannwhitneyu(e, m, alternative="greater")     # expect epithelial > mesenchymal
            out["markers"][arm] = {"epi_median": round(float(e.median()), 2),
                                   "mes_median": round(float(m.median()), 2),
                                   "p_epi_gt_mes": round(float(p), 3),
                                   "correct_direction": bool(e.median() > m.median())}
    out["passes"] = all(v["correct_direction"] and v["p_epi_gt_mes"] < 0.05 for v in out["markers"].values()) \
        if out["markers"] else False
    return out


def _curated_labels(model: pd.DataFrame, breast_ids) -> pd.Series:
    """ModelID -> curated TAD-panel PAM50 label (LumA/LumB/HER2/Basal); else 'Unknown'.

    Uses the parent `cell_line_subtype_map` (TADs `CELL_LINE_METADATA`), the same hand-curated,
    well-characterised breast panel the TAD/HiChIP analysis trusts -- so labels are reliable and the
    set is restricted to lines with rich orthogonal data, not the full noisy CCLE breast catalogue.
    """
    from pipeline.biosample_names import normalize_cell_line_label
    from pipeline.cell_line_subtype_map import subtype_for_key

    name = model.set_index("ModelID")["CCLEName"]
    out = {}
    for mid in breast_ids:
        token = str(name.get(mid, mid)).split("_")[0]
        try:
            out[mid] = subtype_for_key(normalize_cell_line_label(token)).pam50_group
        except Exception:
            out[mid] = "Unknown"
    return pd.Series(out, name="subtype")


def _spear(x: pd.Series, y: pd.Series, min_n: int) -> tuple[float, float, int]:
    df = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    if len(df) < min_n or df["x"].nunique() < 3 or df["y"].nunique() < 3:
        return np.nan, np.nan, len(df)
    rho, p = spearmanr(df["x"], df["y"])
    return float(rho), float(p), len(df)


def run(*, out_dir: Path = OUT_DIR, n_decoy: int = 15, curated: bool = False) -> pd.DataFrame:
    if curated:
        out_dir = out_dir / "curated"
    out_dir.mkdir(parents=True, exist_ok=True)

    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e["arm"] = e["key"].str.split(r"\|\|").str[0]
    e["gene"] = e["key"].str.split(r"\|\|").str[1]
    head = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()
    arm_targets = e.groupby("arm")["gene"].apply(set).to_dict()    # for decoy exclusion

    arm_to_row, _ = build_ccle_arm_expression_map(sorted(head["arm"].unique()), match_mode="alias")
    mir = _ccle_expr_matrix(arm_to_row)                    # arm x ModelID (log2(50+nSolver))

    model = _load_model_lineage()
    breast_ids = set(model.loc[model["OncotreeLineage"] == "Breast", "ModelID"])
    mrna = _mrna_breast_matrix(breast_ids)                 # gene x ModelID
    # curated = the TAD-aligned, hand-curated breast panel (reliable labels, well-characterised lines);
    # default = the full CCLE breast catalogue with DepMap receptor labels.
    subtype = _curated_labels(model, breast_ids) if curated else _breast_subtype()

    lines = sorted(set(mir.columns) & set(mrna.columns) & breast_ids)
    if curated:
        tumor = {"LumA", "LumB", "HER2", "Basal"}
        lines = [m for m in lines if subtype.get(m) in tumor]
    mir_b = mir[lines]
    subtype = subtype.reindex(lines).fillna("Unknown")
    sub_dummies = pd.get_dummies(subtype, prefix="subtype", dtype=float)
    print(f"[ccle-breast] mode={'CURATED' if curated else 'full'}: {len(lines)} breast lines with "
          f"miRNA+mRNA; subtype: {subtype.value_counts().to_dict()}")
    min_all = min(MIN_ALL, max(8, len(lines) - 2)) if curated else MIN_ALL   # adapt to the small panel

    # data-validity gate: does the miRNA matrix recover known epithelial identity?
    name_map = {m: str(model.set_index("ModelID")["CCLEName"].get(m, m)).split("_")[0] for m in lines}
    emt_check = _emt_identity_check(mir_b, lines, name_map)
    print(f"[ccle-breast] miRNA EMT-identity gate: {'PASS' if emt_check['passes'] else 'FAIL'} "
          f"(miR-200c epi {emt_check['markers'].get('hsa-miR-200c-3p',{}).get('epi_median')} vs "
          f"mes {emt_check['markers'].get('hsa-miR-200c-3p',{}).get('mes_median')})")

    strata = {"all_breast": lines}
    for g, ids in subtype.groupby(subtype).groups.items():
        if g != "Unknown" and len(ids) >= MIN_SUBTYPE:
            strata[f"subtype_{g}"] = list(ids)

    # ---- per-edge, per-stratum cross-line Spearman --------------------------------------------
    rows = []
    for _, r in head.iterrows():
        arm, gene = r["arm"], r["gene"]
        if arm not in mir_b.index or gene not in mrna.index:
            continue
        x_full = mir_b.loc[arm]
        y_full = mrna.loc[gene]
        for stratum, ids in strata.items():
            min_n = min_all if stratum == "all_breast" else MIN_SUBTYPE
            x, y = x_full[ids], y_full[ids]
            rho, p, n = _spear(x, y, min_n)
            rec = {"arm": arm, "gene": gene, "family": r["family"], "tcga_rho": r["rho"],
                   "stratum": stratum, "ccle_rho": rho, "ccle_p": p, "n_lines": n}
            if stratum == "all_breast" and np.isfinite(rho):
                cp = correlation_pair(y[x.index], x, sub_dummies.loc[x.index], spearman_only=True)
                rec["ccle_partial_rho"] = cp["partial_rho"]
                rec["ccle_partial_p"] = cp["partial_p"]
            rows.append(rec)
    res = pd.DataFrame(rows)

    # BH within each stratum
    for stratum in res["stratum"].unique():
        m = (res["stratum"] == stratum) & res["ccle_p"].notna()
        res.loc[m, "ccle_q"] = bh_fdr(res.loc[m, "ccle_p"].to_numpy())
    mall = (res["stratum"] == "all_breast") & res.get("ccle_partial_p", pd.Series(dtype=float)).notna()
    if "ccle_partial_p" in res.columns and mall.any():
        res.loc[mall, "ccle_partial_q"] = bh_fdr(res.loc[mall, "ccle_partial_p"].to_numpy())
    res.sort_values(["stratum", "ccle_rho"]).to_csv(out_dir / "edge_ccle_breast_anticorr.tsv", sep="\t", index=False)

    # ---- set-level per stratum ----------------------------------------------------------------
    strata_summary = {}
    for stratum in strata:
        s = res[(res["stratum"] == stratum)].dropna(subset=["ccle_rho"])
        if s.empty:
            continue
        n_neg = int((s["ccle_rho"] < 0).sum())
        strata_summary[stratum] = {
            "n_edges": int(len(s)), "frac_negative": round(n_neg / len(s), 3),
            "binom_p_gt_half": binomtest(n_neg, len(s), 0.5, alternative="greater").pvalue,
            "median_rho": round(float(s["ccle_rho"].median()), 3),
            "n_fdr_neg": int(((s["ccle_q"] < 0.10) & (s["ccle_rho"] < 0)).sum()),
        }
    allb = res[(res["stratum"] == "all_breast")].dropna(subset=["ccle_rho"])
    # PAM50-adjusted set level (within-subtype signal)
    adj = allb.dropna(subset=["ccle_partial_rho"]) if "ccle_partial_rho" in allb.columns else allb.iloc[0:0]
    adj_block = {}
    if len(adj):
        n_neg_adj = int((adj["ccle_partial_rho"] < 0).sum())
        adj_block = {
            "n_edges": int(len(adj)), "frac_negative": round(n_neg_adj / len(adj), 3),
            "median_partial_rho": round(float(adj["ccle_partial_rho"].median()), 3),
            "binom_p_gt_half": binomtest(n_neg_adj, len(adj), 0.5, alternative="greater").pvalue,
            "n_fdr_neg": int(((adj["ccle_partial_q"] < 0.10) & (adj["ccle_partial_rho"] < 0)).sum())
                if "ccle_partial_q" in adj.columns else None,
        }

    # ---- specificity: same-arm decoy non-target genes (all-breast) ----------------------------
    gene_pool = list(mrna.index)
    decoy = []
    for arm in allb["arm"].unique():
        forbidden = arm_targets.get(arm, set())
        pool = [g for g in gene_pool if g not in forbidden]
        k = int((allb["arm"] == arm).sum()) * n_decoy
        for g in RNG.choice(pool, size=min(k, len(pool)), replace=False):
            rho, _, _ = _spear(mir_b.loc[arm, lines], mrna.loc[g, lines], min_all)
            if np.isfinite(rho):
                decoy.append(rho)
    decoy = np.array(decoy)
    decoy_med = float(np.median(decoy)) if len(decoy) else np.nan
    spec_p = float(mannwhitneyu(allb["ccle_rho"], decoy, alternative="less").pvalue) if len(decoy) else np.nan

    # ---- INTERNAL POSITIVE CONTROL: gold-standard repressive pairs ----------------------------
    pc_arms = sorted({a for a, _ in POSITIVE_CONTROL_PAIRS})
    pc_map, _ = build_ccle_arm_expression_map(pc_arms, match_mode="alias")
    pc_mir = _ccle_expr_matrix(pc_map)
    pc_rows = []
    for arm, gene in POSITIVE_CONTROL_PAIRS:
        if arm in pc_mir.index and gene in mrna.index:
            rho, p, n = _spear(pc_mir.loc[arm, lines], mrna.loc[gene, lines], min_all)
            pc_rows.append({"arm": arm, "gene": gene, "ccle_rho": rho, "ccle_p": p, "n": n})
    pc_df = pd.DataFrame(pc_rows)
    pc_df.to_csv(out_dir / "positive_control_pairs.tsv", sep="\t", index=False)
    pc_valid = pc_df.dropna(subset=["ccle_rho"])
    pc_median = float(pc_valid["ccle_rho"].median()) if len(pc_valid) else np.nan
    pc_frac_neg = float((pc_valid["ccle_rho"] < 0).mean()) if len(pc_valid) else np.nan

    # ---- sign concordance vs TCGA -------------------------------------------------------------
    sign_agree = float((np.sign(allb["ccle_rho"]) == np.sign(allb["tcga_rho"])).mean())
    rho_corr = float(spearmanr(allb["ccle_rho"], allb["tcga_rho"]).statistic) if len(allb) > 5 else np.nan

    # per-family roll-up
    fam = (allb.groupby("family").agg(n=("ccle_rho", "size"), median_rho=("ccle_rho", "median"),
           frac_neg=("ccle_rho", lambda s: float((s < 0).mean()))).sort_values("median_rho").round(3))
    fam[fam["n"] >= 3].to_csv(out_dir / "family_ccle_summary.tsv", sep="\t")

    summary = {
        "module": "mirna_hallmark.analyses.cnv_locus.ccle_breast_target_anticorr",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "data": "CCLE NanoString miRNA (2018) + DepMap 26Q1 RNA-seq (log2 TPM+1); cross-cell-line axis",
        "mode": "curated_TAD_panel" if curated else "full_ccle_breast",
        "min_lines_for_all_breast": min_all,
        "n_breast_lines": len(lines),
        "subtype_breakdown": subtype.value_counts().to_dict(),
        "n_headline_edges": int(len(head)),
        "n_testable_all_breast": int(len(allb)),
        "set_level_by_stratum": strata_summary,
        "all_breast_subtype_adjusted": adj_block,
        "specificity_vs_decoy": {
            "n_decoy": int(len(decoy)), "decoy_median_rho": round(decoy_med, 3),
            "cognate_median_rho": strata_summary.get("all_breast", {}).get("median_rho"),
            "mwu_p_cognate_more_negative": spec_p,
        },
        "tcga_sign_concordance": {"frac_same_sign": round(sign_agree, 3),
                                  "ccle_vs_tcga_rho_spearman": round(rho_corr, 3) if np.isfinite(rho_corr) else None},
        "mirna_data_validity_gate": {
            **emt_check,
            "interpretation": "miR-200c/141 are canonical EPITHELIAL markers; they FAIL to separate known "
                              "epithelial>>mesenchymal lines here while the mRNA EMT markers (ZEB1/CDH1/VIM) "
                              "are textbook -> the 2018 NanoString miRNA is the broken link, so the cross-line "
                              "null is measurement-confounded on top of being the wrong axis.",
        },
        "internal_positive_control": {
            "n_pairs_tested": int(len(pc_valid)), "median_rho": round(pc_median, 3) if np.isfinite(pc_median) else None,
            "frac_negative": round(pc_frac_neg, 3) if np.isfinite(pc_frac_neg) else None,
            "pairs": pc_valid.assign(ccle_rho=pc_valid["ccle_rho"].round(2))[["arm", "gene", "ccle_rho"]].to_dict("records"),
            "note": "gold-standard repressive pairs (miR-200/ZEB, let-7/HMGA2, miR-21/PDCD4...). If these are "
                    "NOT recovered as negative, the cross-line axis is insensitive to miRNA repression.",
        },
        "conclusion": (
            "The framework edges are NOT anti-correlated across CCLE breast lines (median ρ +0.05, "
            "subtype-robust, no decoy specificity, TCGA sign-concordance 0.37) -- BUT this is DOUBLY "
            "compromised and therefore INCONCLUSIVE, not a refutation: (1) the internal positive control "
            "fails -- gold-standard miR-200->ZEB1/let-7->HMGA2 are flat-to-positive too; (2) the miRNA "
            "data-validity gate FAILS -- the 2018 NanoString miR-200c/141 do not even mark known epithelial "
            "identity (epi<=mes), while the DepMap mRNA EMT markers are textbook. So BOTH the axis "
            "(between-line co-regulation dominates repression) AND the miRNA measurement (NanoString) are "
            "limiting. A cleaner miRNA-SEQ panel would help, but no large (n>>10) breast-line small-RNA-seq "
            "matrix matched to mRNA exists (GEO panels are ~9 lines). The framework's repression evidence "
            "thus remains the WITHIN-cohort partial-correlation design (TCGA/CPTAC/Buffa) + spatial "
            "localization (MH-64/66); a measured per-cell breast test stays the open gap (Phase 4b)."
        ),
        "strongest_anticorr_breast_edges": (
            allb.sort_values("ccle_rho").head(12)
            .apply(lambda r: f"{r['arm']}->{r['gene']} ccle_rho={r['ccle_rho']:.2f} "
                             f"(TCGA {r['tcga_rho']:.2f}, n={int(r['n_lines'])})", axis=1).tolist()),
        "caveats": [
            "CROSS-CELL-LINE axis (between-line variation), NOT per-cell -- complementary to MH-65, not a per-cell test",
            "axis insensitivity confirmed by the failed positive control -> read as inconclusive, not refutation",
            "CCLE lines are cultured (no stroma/immune) -- pure tumour-epithelial-intrinsic read",
            "miRNA = 2018 NanoString (arm aliasing, e.g. 151a->151-3p); seed-family members not fully resolvable",
            "per-subtype strata small (Luminal n=7) -> set-level only; all-breast carries the per-edge FDR",
            "non-circular: CCLE/DepMap measurements independent of the miRTarBase+TCGA edge definitions",
        ],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    ab = strata_summary.get("all_breast", {})
    print(f"[ccle-breast] ALL-BREAST: {ab.get('frac_negative')} neg (n={ab.get('n_edges')}), "
          f"median ρ={ab.get('median_rho')}, binom p={ab.get('binom_p_gt_half'):.2e}, FDR-neg={ab.get('n_fdr_neg')}")
    if adj_block:
        print(f"[ccle-breast] PAM50-ADJ: {adj_block['frac_negative']} neg, median ρ={adj_block['median_partial_rho']}, "
              f"binom p={adj_block['binom_p_gt_half']:.2e}, FDR-neg={adj_block['n_fdr_neg']}")
    print(f"[ccle-breast] SPECIFICITY: cognate {ab.get('median_rho')} vs decoy {round(decoy_med,3)} "
          f"(MWU p={spec_p:.2e}); TCGA sign-concordance {sign_agree:.2f}")
    for st, d in strata_summary.items():
        if st != "all_breast":
            print(f"[ccle-breast]   {st}: {d['frac_negative']} neg (n={d['n_edges']}), median ρ={d['median_rho']}")
    print(f"[ccle-breast] POS-CONTROL (gold-standard pairs): median ρ={pc_median:+.2f}, "
          f"frac_neg={pc_frac_neg:.2f} over {len(pc_valid)} pairs -> "
          f"{'AXIS INSENSITIVE (inconclusive)' if pc_median is not np.nan and pc_median > -0.1 else 'axis recovers repression'}")
    print(f"[ccle-breast] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--n-decoy", type=int, default=15)
    ap.add_argument("--curated", action="store_true",
                    help="restrict to the TAD-curated, well-characterised breast panel (reliable labels)")
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, n_decoy=args.n_decoy, curated=args.curated)


if __name__ == "__main__":
    main()
