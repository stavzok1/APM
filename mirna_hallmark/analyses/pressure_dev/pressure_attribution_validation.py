"""Does miRNA PRESSURE matter over naive abundance-sum? — attribution + boosted-decoupling validation.

Joins existing outputs (no heavy compute) to make the case that the constructed pressure earns its
complexity, on concrete genes, two ways:

  (1) ATTRIBUTION — the structural softmax share re-ranks which arm dominates a gene's regulation
      (`dominant_flipped`); when the *surfaced* arm is a literature-backed regulator (total functional
      PMIDs, pan-context — NOT breast-only) that the abundance-dominant arm masks, pressure is revealing
      real regulation that a naive abundance-sum buries. (e.g. ABCA1: abundance→miR-93-5p but
      structural→miR-33b-5p, the textbook cholesterol-efflux ABCA1 regulator — note this canonical edge
      has LOWER miRTarBase evidence than the masked miR-93, so "canonical" is not the evidence_score.)
  (2) COUPLING BOOST — genes whose net-repression partial-ρ(P_agg, y_g | CPE+HRD+target_cn) is
      MORE negative under pressure than under abundance-sum (the MH-79 +96 net-new).

The headline `canonical_revealed` set = net-repressed AND pressure-boosted AND dominant_flipped AND the
surfaced arm carries ≥1 total PMID — i.e. pressure both *detects* repression abundance-sum misses and
*attributes* it to a literature-backed regulator abundance-sum masks. The `--robustness` mode re-tests
this shortlist with ESTIMATE composition + proliferation covariates (drops the lineage/stromal artifacts,
MH-51 miR-29→ECM class).

CAVEATS: (a) FORMULAS §5a — the share is the ATTRIBUTION axis and is ⊥ realized coupling; a high
structural share is NOT itself evidence of repression. Each case needs BOTH the coupling (that the gene
is repressed) AND the literature (that the surfaced arm is real); the share only says *which* arm. (b)
"Canonical" is a NOISY automated proxy (PMID count / evidence_score don't capture textbook-ness); the
shortlist is candidate-grade and the strongest claims need per-case curation.

Run:  .venv/bin/python3 -m mirna_hallmark.pressure_attribution_validation [--robustness]
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import gene_roles

OUT_DIR = C.OUTPUT_ROOT / "pressure_attribution_validation"
CMP = C.OUTPUT_ROOT / "coupling_predictor_comparison" / "gene_construction_grid_comparison.tsv.gz"
STRUCT = C.TISSUE_REFERENCE_DIR / "mirna_state_class" / "struct_vs_abundance_resolution.tsv"
# per-edge literature: evidence_score + total PMIDs (pan-context canonical, e.g. miR-33→ABCA1 cholesterol)
# + breast-specific PMIDs (secondary). Broader than breast-only.
RERANK = C.TISSUE_REFERENCE_DIR / "edge_breast_context" / "per_target_arm_reranking.tsv"
COREP = C.TISSUE_REFERENCE_DIR / "mirna_comovement" / "gene_corepression.tsv"

ABUND = "abundance_sum"
SPINE = "softmax_z_logrpm"          # the framework's pressure
NATSH = "shr|softmax|nat|p0"        # the MH-79 winner (stronger variant)
# Positive controls from the registry (MH-31/32): canonical breast regulators abundance-sum masks.
POSITIVE_CONTROLS = {
    "ABCA1": "hsa-miR-33b-5p", "PTEN": "hsa-miR-21-5p",
    "ERBB2": "hsa-miR-125b-5p", "EGFR": "hsa-miR-146a-5p", "ESR1": "hsa-miR-18a-5p",
}


def _gene_rho(cmp: pd.DataFrame) -> pd.DataFrame:
    """Per-gene partial-ρ + neg-sig for the 3 predictors at the (all, ungated) baseline config."""
    s = cmp[(cmp["prune"] == "all") & (cmp["gating"] == "ungated") &
            (cmp["predictor"].isin([ABUND, SPINE, NATSH]))]
    rho = s.pivot_table(index="gene", columns="predictor", values="partial_rho_cn")
    sig = s.pivot_table(index="gene", columns="predictor", values="neg_sig_cn").fillna(False).astype(bool)
    rho.columns = [f"rho_{c}" for c in rho.columns]
    sig.columns = [f"negsig_{c}" for c in sig.columns]
    g = rho.join(sig)
    g["boost_spine"] = g.get(f"rho_{SPINE}") - g.get(f"rho_{ABUND}")     # <0 = pressure more repressive
    g["boost_nat"] = g.get(f"rho_{NATSH}") - g.get(f"rho_{ABUND}")
    g["pressure_netnew_spine"] = g.get(f"negsig_{SPINE}", False) & ~g.get(f"negsig_{ABUND}", False)
    g["pressure_netnew_nat"] = g.get(f"negsig_{NATSH}", False) & ~g.get(f"negsig_{ABUND}", False)
    return g


def build(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    cmp = pd.read_csv(CMP, sep="\t")
    g = _gene_rho(cmp)

    struct = pd.read_csv(STRUCT, sep="\t").set_index("gene")
    struct = struct.rename(columns={"struct_dominant_mirna": "surfaced_arm",
                                    "abund_dominant_mirna": "masked_arm"})
    g = g.join(struct[["n_regulators", "surfaced_arm", "masked_arm", "dominant_flipped",
                       "abund_dominant_struct_rank", "within_gene_spearman"]], how="left")

    # canonical literature (non-circular, lit-only): per-edge total PMIDs + evidence_score + breast PMIDs.
    # Join for BOTH the surfaced (struct-dominant) and masked (abundant) arm, to ask: does the share
    # surface a BETTER-EVIDENCED regulator that abundance buried?
    rr = pd.read_csv(RERANK, sep="\t")[["miRNA", "gene", "evidence_score", "n_pmids", "n_breast_pmids"]]
    g = g.reset_index()
    surf = rr.rename(columns={"miRNA": "surfaced_arm", "evidence_score": "surfaced_evidence",
                              "n_pmids": "surfaced_n_pmids", "n_breast_pmids": "surfaced_n_breast_pmids"})
    g = g.merge(surf, on=["gene", "surfaced_arm"], how="left")
    msk = rr.rename(columns={"miRNA": "masked_arm", "evidence_score": "masked_evidence", "n_pmids": "masked_n_pmids"})
    g = g.merge(msk[["gene", "masked_arm", "masked_evidence", "masked_n_pmids"]], on=["gene", "masked_arm"], how="left")
    for c in ("surfaced_n_pmids", "surfaced_n_breast_pmids", "masked_n_pmids"):
        g[c] = g[c].fillna(0).astype(int)
    g["surfaced_stronger_evidence"] = g["surfaced_evidence"].fillna(0) >= g["masked_evidence"].fillna(0)

    corep = pd.read_csv(COREP, sep="\t").set_index("gene")
    g = g.merge(corep[["rho_gene_pressure_tumor", "gene_repression_class", "gene_net_repressed_tumor"]],
                left_on="gene", right_index=True, how="left")

    roles = gene_roles.load_gene_roles(g["gene"].tolist())
    if "gene" in roles.columns:
        roles = roles.set_index("gene")
    g = g.merge(roles[[c for c in ("role", "malignancy_sign") if c in roles.columns]],
                left_on="gene", right_index=True, how="left")

    g["is_positive_control"] = g.apply(
        lambda r: POSITIVE_CONTROLS.get(r["gene"]) == r["surfaced_arm"], axis=1)
    # headline: pressure DETECTS repression abundance-sum misses AND ATTRIBUTES it to a canonical arm
    g["net_repressed"] = g["gene_net_repressed_tumor"].fillna(False).astype(bool) | g["pressure_netnew_spine"]
    g["pressure_boosted"] = (g["boost_spine"] < 0) | (g["boost_nat"] < 0)
    # canonical = surfaced arm has literature support (total PMIDs ≥1, pan-context — not breast-only)
    g["canonical_revealed"] = (g["net_repressed"] & g["pressure_boosted"] &
                               g["dominant_flipped"].fillna(False).astype(bool) &
                               (g["surfaced_n_pmids"] >= 1))

    out_dir.mkdir(parents=True, exist_ok=True)
    g = g.sort_values(["canonical_revealed", "surfaced_n_pmids", "boost_nat"],
                      ascending=[False, False, True]).reset_index(drop=True)
    g.to_csv(out_dir / "validated_cases.tsv", sep="\t", index=False)
    rev = g[g["canonical_revealed"]].copy()
    rev.to_csv(out_dir / "canonical_revealed.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.pressure_attribution_validation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "question": "does pressure matter over abundance_sum: attribution (canonical arm masked by abundance) + coupling boost",
        "pressure_primary": SPINE, "pressure_strong": NATSH, "baseline": ABUND,
        "config": "prune=all, gating=ungated; partial-rho CPE+HRD+target_cn",
        "n_genes": int(len(g)), "n_dominant_flipped": int(g["dominant_flipped"].fillna(False).sum()),
        "n_pressure_boosted": int(g["pressure_boosted"].sum()),
        "n_canonical_revealed": int(g["canonical_revealed"].sum()),
        "caveat": "share is ATTRIBUTION (⊥ coupling, FORMULAS §5a); cases require coupling + canonical literature, not share alone",
        "inputs": [str(CMP), str(STRUCT), str(RERANK), str(COREP)],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return g


def robustness(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """Composition/proliferation robustness re-test on the canonical_revealed shortlist.

    Re-couples the SPINE pressure (softmax_z_logrpm) to gene expression adding ESTIMATE
    immune+stromal+epithelial+proliferation covariates on top of CPE+HRD+target_cn. A net-repression
    that is killed/flipped by composition is a lineage/composition artifact (MH-51 miR-29→ECM class),
    not regulation, and is dropped from the biological headline. Non-circular: same covariates regardless
    of edge."""
    from analysis.utils.common.loaders import partial_spearman
    from mirna_hallmark.hallmark_sets import HallmarkSets
    from mirna_hallmark.analyses.cross_state.cross_state_coupling import _state_covariates
    from mirna_hallmark.coupling_predictor_comparison import (
        _load_grid_inputs, _engine_gene_pressure, _cov_for_gene, _row, TARGET_NORM)

    rev = pd.read_csv(out_dir / "canonical_revealed.tsv", sep="\t")
    genes = rev["gene"].tolist()
    hs, edges, mirna, rna, cov_base, _conf, _pg, target_cn = _load_grid_inputs(genes)
    P = _engine_gene_pressure(edges, mirna, genes=genes, expr_mode=SPINE,
                              target_norm=TARGET_NORM, aggregate="sum", contrib_transform="signed")
    comp = _state_covariates("tumor", rna, hs, "estimate_epi")   # immune+stromal+epi+prolif (+batch)

    rows = []
    for _, r in rev.iterrows():
        gene = r["gene"]
        if gene not in P.index or gene not in rna.index:
            continue
        y = _row(rna, gene)
        x = _row(P, gene)
        cov_full = _cov_for_gene(cov_base, target_cn, gene).join(comp, how="left")
        shared = sorted(set(y.dropna().index) & set(x.dropna().index) & set(cov_full.dropna(how="all").index))
        if len(shared) < 30:
            continue
        adj_rho, adj_p, n = partial_spearman(y.loc[shared], x.loc[shared], cov_full.reindex(shared))
        base = r.get("rho_softmax_z_logrpm", np.nan)
        survives = bool((adj_rho < 0) and np.isfinite(adj_p) and (adj_p < C.FDR_ALPHA))
        stromal_class = bool(("miR-29" in str(r["masked_arm"])) or ("miR-29" in str(r["surfaced_arm"]))
                             or str(gene).startswith(("COL", "MMP", "FN1", "SPARC", "ELN")))
        rows.append({"gene": gene, "surfaced_arm": r["surfaced_arm"], "masked_arm": r["masked_arm"],
                     "role": r.get("role"), "rho_base_cn": base, "rho_composition_adj": round(float(adj_rho), 4),
                     "p_composition_adj": adj_p, "n": int(n),
                     "survives_composition": survives,
                     "composition_sensitive": bool((base < 0) and not survives),
                     "stromal_artifact_class": stromal_class,
                     "surfaced_n_pmids": int(r["surfaced_n_pmids"])})
    rob = pd.DataFrame(rows).sort_values(["survives_composition", "rho_composition_adj"],
                                         ascending=[False, True]).reset_index(drop=True)
    rob.to_csv(out_dir / "canonical_revealed_robustness.tsv", sep="\t", index=False)
    n_surv = int(rob["survives_composition"].sum())
    n_strom = int(rob["stromal_artifact_class"].sum())
    print(f"[pav-robustness] {len(rob)} canonical_revealed re-tested | survives_composition={n_surv} "
          f"| composition_sensitive={int(rob['composition_sensitive'].sum())} | stromal_class={n_strom}")
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print("\n=== survivors (composition+proliferation-robust net-repression) ===")
        print(rob[rob["survives_composition"]].head(25).to_string(index=False))
    return rob


def build_budget(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """Per-(gene, arm) REGULATORY BUDGET through three lenses, for all multi-regulator genes:
      • static weight    = evidence_score share within the gene (what curation says is a regulator)
      • abundance        = linear-RPM share within the gene (how MUCH miRNA — the abundance_sum view)
      • attribution      = realized pressure share `global_abs_share` (WHICH miRNA — our share).
    """
    from mirna_hallmark.pressure_engine import compute_gene_pressure_contributions
    from mirna_hallmark.coupling_predictor_comparison import _load_grid_inputs, TARGET_NORM

    struct = pd.read_csv(STRUCT, sep="\t")
    multi = struct.loc[struct["n_regulators"] >= 2, "gene"].tolist()
    _hs, edges, mirna, _rna, _cb, _cf, _pg, _tcn = _load_grid_inputs(multi)
    # attribution lens = the STRUCTURAL share (softmax of abundance-deviation, promiscuity-discounted; no
    # z) — the "which arm" attribution that defines surfaced/masked (matches struct_vs_abundance). The
    # z-weighted realized share is downstream (tracks coupling), not the attribution we reorder by.
    con = compute_gene_pressure_contributions(edges, mirna, genes=multi, expr_mode="softmax", target_norm=TARGET_NORM)
    con = con[["gene", "miRNA", "evidence_score", "median_log2rpm", "global_abs_share"]].copy()
    # 4th lens: realized magnitude share (softmax_z_logrpm) — where the actual FORCE lands (vs identity)
    mag = compute_gene_pressure_contributions(edges, mirna, genes=multi, expr_mode="softmax_z_logrpm",
                                              target_norm=TARGET_NORM).set_index(["gene", "miRNA"])["global_abs_share"]
    con["w_magnitude"] = con.set_index(["gene", "miRNA"]).index.map(mag).fillna(0.0)
    con["rpm"] = np.power(2.0, con["median_log2rpm"]) - 1.0
    gb = con.groupby("gene")
    con["w_evidence"] = (con["evidence_score"] / gb["evidence_score"].transform("sum")).fillna(0.0)
    con["w_abundance"] = (con["rpm"] / gb["rpm"].transform("sum")).fillna(0.0)
    con["w_attribution"] = con["global_abs_share"].fillna(0.0)
    # within-gene ranks (1 = top) under each lens
    for lens in ("w_evidence", "w_abundance", "w_attribution", "w_magnitude"):
        con[lens.replace("w_", "rank_")] = gb[lens].rank(ascending=False, method="min")
    out_dir.mkdir(parents=True, exist_ok=True)
    con.to_csv(out_dir / "regulatory_budget.tsv", sep="\t", index=False)
    return con


def identity_magnitude(out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """The identity vs magnitude diagnostic (ATTRIBUTION_IDENTITY_VS_MAGNITUDE §7, hyps 1 & 4).

    Per (gene, arm), three lenses from `compute_gene_pressure_contributions`:
      identity  = global_abs_share under expr_mode='softmax'        (who; abundance removed)
      level     = global_abs_share under expr_mode='softmax_logrpm' (identity × absolute level)
      magnitude = global_abs_share under expr_mode='softmax_z_logrpm' (realized responsive pressure)
    Plane = identity (x) vs magnitude (y); the off-diagonal quadrants are the biology:
      high-identity / low-magnitude = silenced specialist (the TSG-miRNA-loss signature)
      low-identity  / high-magnitude = loud promiscuous hub.
    Drift = per canonical gene, the top arm under softmax -> softmax_logrpm -> softmax_z_logrpm
    (hypothesis 1: adding magnitude drifts the surfaced arm specialist -> abundant)."""
    from mirna_hallmark.pressure_engine import compute_gene_pressure_contributions
    from mirna_hallmark.coupling_predictor_comparison import _load_grid_inputs, TARGET_NORM

    struct = pd.read_csv(STRUCT, sep="\t")
    multi = struct.loc[struct["n_regulators"] >= 2, "gene"].tolist()
    _hs, edges, mirna, _rna, _cb, _cf, _pg, _tcn = _load_grid_inputs(multi)
    modes = {"identity": "softmax", "level": "softmax_logrpm", "magnitude": "softmax_z_logrpm"}
    cols = {}
    for name, mode in modes.items():
        c = compute_gene_pressure_contributions(edges, mirna, genes=multi, expr_mode=mode, target_norm=TARGET_NORM)
        cols[name] = c.set_index(["gene", "miRNA"])["global_abs_share"].rename(name)
    pm = pd.concat(cols.values(), axis=1).reset_index()
    # global-median quadrant split on identity vs magnitude
    mi, mm = pm["identity"].median(), pm["magnitude"].median()
    pm["quadrant"] = np.where(pm["identity"] >= mi,
                              np.where(pm["magnitude"] >= mm, "loud_specialist", "silenced_specialist"),
                              np.where(pm["magnitude"] >= mm, "loud_hub", "minor"))
    out_dir.mkdir(parents=True, exist_ok=True)
    pm.to_csv(out_dir / "identity_magnitude_plane.tsv", sep="\t", index=False)

    # drift: top arm per gene under each mode, on the canonical genes
    canon = pd.read_csv(out_dir / "canonical_revealed.tsv", sep="\t")["gene"].tolist()
    drift = []
    for g in canon:
        sub = pm[pm["gene"] == g]
        if sub.empty:
            continue
        top = {m: sub.loc[sub[m].idxmax(), "miRNA"] for m in modes}
        drift.append({"gene": g, "top_identity": top["identity"], "top_level": top["level"],
                      "top_magnitude": top["magnitude"],
                      "drifted_id_to_mag": top["identity"] != top["magnitude"]})
    drift = pd.DataFrame(drift)
    drift.to_csv(out_dir / "magnitude_drift.tsv", sep="\t", index=False)
    n_drift = int(drift["drifted_id_to_mag"].sum()) if len(drift) else 0
    print(f"[pav-idmag] {len(pm)} (gene,arm) pairs | quadrants: "
          f"{pm['quadrant'].value_counts().to_dict()}")
    print(f"[pav-idmag] magnitude drift (identity-top != magnitude-top): {n_drift}/{len(drift)} canonical genes")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fdir = out_dir / "figures"; fdir.mkdir(parents=True, exist_ok=True)
    qcol = {"silenced_specialist": "#d1495b", "loud_hub": "#3a7ca5", "loud_specialist": "#6a4c93", "minor": "#cccccc"}
    fig, ax = plt.subplots(1, 2, figsize=(15, 6.2))
    for q, c in qcol.items():
        s = pm[pm["quadrant"] == q]
        ax[0].scatter(s["identity"], s["magnitude"], s=7, alpha=0.3, color=c, label=q.replace("_", " "))
    ax[0].axvline(mi, color="k", lw=0.7, ls="--"); ax[0].axhline(mm, color="k", lw=0.7, ls="--")
    mark = {"NCOA3": "hsa-miR-137", "MAP3K8": "hsa-miR-370-3p", "FOXO1": "hsa-miR-9-5p", "BTG2": "hsa-miR-21-5p"}
    for g, a in mark.items():
        r = pm[(pm["gene"] == g) & (pm["miRNA"] == a)]
        if not r.empty:
            ax[0].annotate(f"{a.replace('hsa-miR-','miR-')}→{g}", (r["identity"].iloc[0], r["magnitude"].iloc[0]),
                           fontsize=7, color="k")
            ax[0].scatter(r["identity"], r["magnitude"], s=40, edgecolor="k", facecolor="none", linewidth=1.2)
    ax[0].set_xlabel("identity share (softmax — who; abundance removed)")
    ax[0].set_ylabel("magnitude share (softmax_z_logrpm — realized force)")
    ax[0].set_title("Identity vs magnitude plane\n(off-diagonal = the biology)")
    ax[0].legend(fontsize=8, frameon=False, markerscale=2)
    # drift panel: how often top arm changes across the ladder
    lad = pd.DataFrame({"identity→level": [(drift.top_identity != drift.top_level).mean()],
                        "identity→magnitude": [(drift.top_identity != drift.top_magnitude).mean()],
                        "level→magnitude": [(drift.top_level != drift.top_magnitude).mean()]}).T[0]
    ax[1].bar(range(len(lad)), lad.values, color="#2a9d8f")
    ax[1].set_xticks(range(len(lad))); ax[1].set_xticklabels(lad.index, fontsize=9)
    ax[1].set_ylabel("fraction of canonical genes whose TOP arm changes")
    ax[1].set_title("Magnitude drift along the construction ladder\n(softmax → +level → +responsiveness)")
    ax[1].set_ylim(0, 1)
    fig.suptitle("Identity (who) vs Magnitude (how much) — they decouple", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(fdir / "fig4_identity_magnitude.png", dpi=140); plt.close(fig)
    print(f"[pav-idmag] wrote identity_magnitude_plane.tsv, magnitude_drift.tsv, fig4_identity_magnitude.png")
    return pm


def make_figures(out_dir: Path = OUT_DIR) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fdir = out_dir / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    bud = build_budget(out_dir)
    rob = pd.read_csv(out_dir / "canonical_revealed_robustness.tsv", sep="\t")
    LENS = ["w_evidence", "w_abundance", "w_attribution", "w_magnitude"]
    LAB = ["static\nevidence", "abundance\n(RPM share)", "attribution\n(identity)", "realized\n(magnitude)"]

    def arm_short(a: str) -> str:
        return str(a).replace("hsa-miR-", "miR-").replace("hsa-let-", "let-")

    # ---- Figure 1: the key reorderings — surfaced (canonical) vs masked (abundant) across 3 lenses ----
    surv = rob[rob["survives_composition"]].sort_values("rho_composition_adj").head(9).reset_index(drop=True)
    fig, axes = plt.subplots(3, 3, figsize=(14, 11)); axes = axes.ravel()
    x = np.arange(len(LENS))
    for i, r in surv.iterrows():
        ax = axes[i]; gene = r["gene"]; sa, ma = r["surfaced_arm"], r["masked_arm"]
        gsub = bud[bud["gene"] == gene].set_index("miRNA")
        sv = [gsub.loc[sa, l] if sa in gsub.index else 0 for l in LENS]
        mv = [gsub.loc[ma, l] if ma in gsub.index else 0 for l in LENS]
        ax.bar(x - 0.2, mv, 0.4, label=f"masked (abundant): {arm_short(ma)}", color="#b0b0b0")
        ax.bar(x + 0.2, sv, 0.4, label=f"surfaced (canonical): {arm_short(sa)}", color="#d1495b")
        ax.set_title(f"{gene}  ({r['role']}, ρ={r['rho_composition_adj']:.2f})", fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels(LAB, fontsize=7); ax.set_ylim(0, 1)
        ax.legend(fontsize=6.5, loc="upper left", frameon=False)
        if i % 3 == 0:
            ax.set_ylabel("within-gene share")
    fig.suptitle("Attribution flips dominance to the canonical regulator the abundant arm masks\n"
                 "(composition-robust net-repressed genes; share sums to 1 within each lens)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(fdir / "fig1_canonical_flips.png", dpi=140); plt.close(fig)

    # ---- Figure 2: full regulatory budget of exemplar genes (every regulator, 3 lenses stacked to 1) ----
    examples = [g for g in ["ABCA1", "FOXO1", "BTG2", "CASP8"] if g in set(bud["gene"])]
    fig, axes = plt.subplots(1, len(examples), figsize=(4.2 * len(examples), 5.2))
    if len(examples) == 1:
        axes = [axes]
    cmap = plt.get_cmap("tab20")
    for ax, gene in zip(axes, examples):
        gsub = bud[bud["gene"] == gene].sort_values("w_abundance", ascending=False).reset_index(drop=True)
        arms = gsub["miRNA"].tolist()
        colors = {a: cmap(j % 20) for j, a in enumerate(arms)}
        for li, lens in enumerate(LENS):
            bottom = 0.0
            for a in arms:
                v = float(gsub.loc[gsub["miRNA"] == a, lens].iloc[0])
                ax.bar(li, v, bottom=bottom, color=colors[a], edgecolor="white", linewidth=0.4)
                if v > 0.10:
                    ax.text(li, bottom + v / 2, arm_short(a), ha="center", va="center", fontsize=6.5)
                bottom += v
        ax.set_xticks(range(len(LENS))); ax.set_xticklabels(LAB, fontsize=7.5); ax.set_ylim(0, 1.001)
        ax.set_title(f"{gene}  ({len(arms)} regulators)", fontsize=11)
        ax.set_ylabel("within-gene share")
    fig.suptitle("How a gene's regulatory budget looks under each lens — the dominant regulator changes", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(fdir / "fig2_budget_composition.png", dpi=140); plt.close(fig)

    # ---- Figure 3: broad overview across all multi-regulator genes ----
    top1 = bud.groupby("gene").agg(top1_abund=("w_abundance", "max"), top1_attr=("w_attribution", "max"))
    struct = pd.read_csv(STRUCT, sep="\t")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].hist(top1_abund := top1["top1_abund"].dropna(), bins=25, alpha=0.6, label="abundance", color="#3a7ca5")
    ax[0].hist(top1["top1_attr"].dropna(), bins=25, alpha=0.6, label="attribution", color="#d1495b")
    ax[0].set_xlabel("top-1 regulator share (budget concentration)"); ax[0].set_ylabel("genes")
    ax[0].legend(frameon=False); ax[0].set_title("How concentrated is the budget?")
    wgs = pd.to_numeric(struct["within_gene_spearman"], errors="coerce").dropna()
    ax[1].hist(wgs, bins=30, color="#6a4c93")
    ax[1].axvline(0, color="k", lw=0.8, ls="--")
    ax[1].set_xlabel("within-gene Spearman(abundance rank, attribution rank)")
    ax[1].set_ylabel("genes")
    ax[1].set_title(f"Abundance vs attribution agreement\n(median {wgs.median():.2f}; "
                    f"{int((struct['dominant_flipped']==True).sum())}/{len(struct)} flip the #1)")
    s = bud.sample(min(6000, len(bud)), random_state=0)
    ax[2].scatter(s["w_abundance"], s["w_attribution"], s=5, alpha=0.25, color="#2a9d8f")
    ax[2].plot([0, 1], [0, 1], "k--", lw=0.8)
    ax[2].set_xlabel("abundance share"); ax[2].set_ylabel("attribution share")
    ax[2].set_title("Per-edge reordering\n(off-diagonal = abundance ≠ attribution)")
    fig.suptitle("Regulatory budget across all multi-regulator genes — three lenses disagree systematically", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(fdir / "fig3_budget_overview.png", dpi=140); plt.close(fig)
    print(f"[pav-figures] wrote 3 figures + regulatory_budget.tsv to {fdir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--robustness", action="store_true",
                    help="re-test the canonical_revealed shortlist with composition+proliferation covariates")
    ap.add_argument("--figures", action="store_true",
                    help="build the three-lens regulatory-budget table + figures")
    ap.add_argument("--identity-magnitude", action="store_true",
                    help="build the identity-vs-magnitude plane + magnitude-drift diagnostic")
    args = ap.parse_args()
    if args.identity_magnitude:
        identity_magnitude(args.out_dir)
        return
    if args.figures:
        make_figures(args.out_dir)
        return
    if args.robustness:
        robustness(args.out_dir)
        return
    g = build(args.out_dir)
    rev = g[g["canonical_revealed"]]
    print(f"[pav] genes={len(g)}  dominant_flipped={int(g['dominant_flipped'].fillna(False).sum())}  "
          f"pressure_boosted={int(g['pressure_boosted'].sum())}  canonical_revealed={len(rev)}")
    print("\n=== positive-control recovery ===")
    pc = g[g["is_positive_control"]]
    cols = ["gene", "surfaced_arm", "masked_arm", "surfaced_n_pmids", "surfaced_evidence", "masked_evidence",
            "surfaced_stronger_evidence", "boost_spine", "boost_nat", "net_repressed", "canonical_revealed", "role"]
    print((pc[cols].to_string(index=False)) if len(pc) else "(no positive controls matched — check naming)")
    print("\n=== top canonical_revealed cases ===")
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(rev.head(20)[cols + ["gene_repression_class"]].to_string(index=False))
    print(f"\n[done] wrote {args.out_dir}")


if __name__ == "__main__":
    main()
