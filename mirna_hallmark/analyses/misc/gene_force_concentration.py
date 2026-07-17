"""Aggregate "force" vs abundance-sum at the gene rung — concentration moderator (step 1)
and force-vs-base coupling on TCGA tumor (step 2).

Construction (per gene *g*, sample *s*):

  base  : c_base(m,g)  = a(m,s)                                   (Σ a — abundance dose)
  force : c_force(m,g) = a(m,s) · W(m,g)        W(m,g) = w_eff(m,g) / D(m)^α

``a`` = miRNA abundance (``log`` = max(log2RPM,0), the abundance_sum currency; ``lin`` =
raw-linear RPM, loading-dominance). ``w_eff`` = per-edge binding weight; ``D(m)`` = the
arm's **genome-wide** binding BUDGET (evidence-mass, NOT a degree) from
`genome_wide_promiscuity.genome_wide_budget` — so ``W`` is the *fraction* of m's budget
spent on g (a true budget-split; within-universe degree overstated it — fixed 2026-06-29).

Two binding currencies, each matched numerator↔denominator:
  - ``validated``  : w = log1p(evidence_score)  / D = Σ_genomewide log1p(evidence)  (HE edges)
  - ``targetscan`` : w = |weighted context++|   / D = Σ_genomewide |context++|
α sweeps promiscuity strength: α=0 ⇒ no penalty (raw binding weight), α=1 ⇒ full budget-share.

Per gene we report the **Hill effective-number spectrum** of the contribution shares
``p_m = c_m/Σc`` (H0 richness ≥ H1 exp-entropy ≥ H2 invSimpson ≥ H∞ 1/top1). H2 is the
primary moderator + falsification cut (`n_eff≥2`): n_eff≈1 ⇒ ~one edge ⇒ edge lemma ⇒
force cannot beat abundance; n_eff high ⇒ construction can matter. See skill
``apm-gene-question`` and ``[[aggregate-force-vs-abundance-design]]``.

Run: ``.venv/bin/python3 -m mirna_hallmark.gene_force_concentration`` (``--concentration-only``
to skip step 2). Outputs: ``gene_force_concentration.tsv`` + ``gene_force_coupling.tsv``.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

import mirna_hallmark.data_loaders as D
from mirna_hallmark import config as C
from mirna_hallmark import stats as S
from mirna_hallmark.coupling_predictor_comparison import OUT_DIR
from mirna_hallmark.analyses.misc.genome_wide_promiscuity import genome_wide_budget
from mirna_hallmark.pressure_build import load_mirtar_edges
from mirna_hallmark.hallmark_sets import HallmarkSets

# (w_src, D_src) matched binding currencies
_CURRENCY = {"val": ("validated", "validated"), "ts": ("targetscan", "targetscan")}
_ALPHAS = (0.0, 0.5, 1.0)


def _edge_weight(edges: pd.DataFrame, w_src: str) -> pd.Series:
    """Per-edge binding weight w_eff aligned to ``edges.index`` (missing → 0)."""
    if w_src == "validated":
        return np.log1p(pd.to_numeric(edges["evidence_score"], errors="coerce").fillna(0.0))
    if w_src == "targetscan":
        ts = _ts_weights_cached(tuple(sorted(edges["gene"].unique())))   # (miRNA,gene,ts_weight)
        m = edges.merge(ts, on=["miRNA", "gene"], how="left")
        return pd.to_numeric(m["ts_weight"], errors="coerce").fillna(0.0).set_axis(edges.index)
    raise ValueError(w_src)


@lru_cache(maxsize=2)
def _ts_weights_cached(genes_key: tuple) -> pd.DataFrame:
    from mirna_hallmark.robustness_checks import _load_targetscan_weights
    return _load_targetscan_weights(list(genes_key))


def _budget_share(edges: pd.DataFrame, w_src: str, d_src: str, alpha: float) -> pd.Series:
    """W(m,g) = w_eff / D(m)^α, aligned to edges.index. α=0 ⇒ raw w_eff (D irrelevant).
    Edges whose arm lacks a genome-wide budget at α>0 get W=0 (excluded from the aggregate)."""
    w = _edge_weight(edges, w_src)
    if alpha == 0.0:
        return w
    D_m = genome_wide_budget(d_src)
    Dvals = edges["miRNA"].map(D_m).astype(float)
    return (w / np.power(Dvals, alpha)).replace([np.inf, -np.inf], np.nan).fillna(0.0)


def _weight_dict(edges: pd.DataFrame, w_src: str, d_src: str, alpha: float) -> dict:
    W = _budget_share(edges, w_src, d_src, alpha)
    return {(m, g): float(v) for m, g, v in zip(edges["miRNA"], edges["gene"], W)}


# ---------------------------------------------------------------------------- #
# Step 1 — per-gene concentration / Hill effective-N spectrum
# ---------------------------------------------------------------------------- #
def _concentration(edges: pd.DataFrame, mass_col: str) -> pd.DataFrame:
    g = edges.groupby("gene")
    tot = g[mass_col].sum()
    sumsq = g[mass_col].apply(lambda s: float(np.square(s.to_numpy(dtype=float)).sum()))
    mx = g[mass_col].max()
    idx = edges.loc[edges.groupby("gene")[mass_col].idxmax(), ["gene", "miRNA"]].set_index("gene")["miRNA"]

    def _h1(s: pd.Series) -> float:
        v = s.to_numpy(dtype=float)
        v = v[np.isfinite(v) & (v > 0)]
        if v.sum() <= 0:
            return np.nan
        p = v / v.sum()
        return float(np.exp(-(p * np.log(p)).sum()))

    return pd.DataFrame({
        "n_regulators": g.size(),
        f"n_eff_h1_{mass_col}": g[mass_col].apply(_h1),
        f"n_eff_{mass_col}": (tot ** 2) / sumsq.replace(0, np.nan),
        f"n_eff_hinf_{mass_col}": tot / mx.replace(0, np.nan),
        f"top1_share_{mass_col}": mx / tot.replace(0, np.nan),
        f"top_arm_{mass_col}": idx,
    })


def _force_edges(genes: Sequence[str] | None):
    """edges + per-edge cohort-median abundance (log/lin) + matched-currency budget-share masses."""
    hs = HallmarkSets.load()
    edges = load_mirtar_edges(list(genes or hs.universe), resolve_arms=True)
    mirna = D.load_mirna_arms()
    a_log = mirna.astype(float).clip(lower=0.0).median(axis=1)
    a_lin = (2.0 ** a_log - 1.0).clip(lower=0.0)
    e = edges[edges["miRNA"].isin(a_log.index)].copy()
    e["a_log"], e["a_lin"] = e["miRNA"].map(a_log), e["miRNA"].map(a_lin)
    e["base_log"], e["base_lin"] = e["a_log"], e["a_lin"]
    for tag, (w_src, d_src) in _CURRENCY.items():
        W = _budget_share(e, w_src, d_src, 1.0)
        e[f"force_{tag}_log"] = e["a_log"] * W
        e[f"force_{tag}_lin"] = e["a_lin"] * W
    return e, edges, mirna


def build(genes: Sequence[str] | None = None, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    e, _edges, _mirna = _force_edges(genes)
    mass_cols = ["base_log", "base_lin", "force_val_log", "force_val_lin", "force_ts_log", "force_ts_lin"]
    e = e.dropna(subset=["base_log"])
    parts = [_concentration(e, c) for c in mass_cols]
    res = parts[0][["n_regulators"]].copy()
    for p in parts:
        res = res.join(p.drop(columns=["n_regulators"]))
    res["reorders_val_log"] = res["top_arm_base_log"] != res["top_arm_force_val_log"]
    res["reorders_ts_log"] = res["top_arm_base_log"] != res["top_arm_force_ts_log"]
    res = res.reset_index().rename(columns={"index": "gene"})

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "gene_force_concentration.tsv"
    res.to_csv(path, sep="\t", index=False)

    multi = res[res["n_regulators"] >= 2]
    print(f"[force-conc] {len(res):,} genes ({len(multi):,} multi-regulator), {len(e):,} edges")
    print("  Hill spectrum (median over multi-reg genes): H1 exp-entropy / H2 invSimpson / H∞ 1/top1")
    for base in mass_cols:
        h1, h2, hinf = (multi[f"n_eff_h1_{base}"].median(), multi[f"n_eff_{base}"].median(),
                        multi[f"n_eff_hinf_{base}"].median())
        print(f"  {base:14s} H1={h1:5.2f} H2={h2:5.2f} H∞={hinf:5.2f}  frac H2>=2: {(multi[f'n_eff_{base}']>=2).mean():.1%}")
    print(f"  budget-split reorders top regulator vs abundance — val:{res['reorders_val_log'].mean():.1%} "
          f"ts:{res['reorders_ts_log'].mean():.1%}")
    print(f"[force-conc] wrote {path}")
    return res


# ---------------------------------------------------------------------------- #
# Step 2 — force-vs-base coupling on TCGA tumor, stratified by n_eff
# ---------------------------------------------------------------------------- #
def _score_predictor(name: str, mat: pd.DataFrame, rna: pd.DataFrame, cov_base: pd.DataFrame,
                     target_cn: pd.DataFrame) -> list:
    from mirna_hallmark.coupling_predictor_comparison import _row, _cov_for_gene, MIN_N
    rna_genes = set(rna.index)
    rows = []
    for gene in mat.index:
        if gene not in rna_genes:
            continue
        x = pd.to_numeric(mat.loc[gene], errors="coerce")
        if x.notna().sum() < MIN_N or x.std(skipna=True) == 0:
            continue
        st = S.correlation_pair(_row(rna, gene), x, _cov_for_gene(cov_base, target_cn, gene))
        if st["n"] < MIN_N:
            continue
        rows.append({"predictor": name, "gene": gene, "n": st["n"],
                     "partial_rho_cn": st["partial_rho"], "partial_p_cn": st["partial_p"]})
    return rows


def force_vs_base_coupling(genes: Sequence[str] | None = None, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    from mirna_hallmark.coupling_predictor_comparison import _abundance_sum_pressure, _load_grid_inputs
    hs, edges, mirna, rna, cov_base, _conf, pressure_genes, target_cn = _load_grid_inputs(genes)
    mats_in = {"log": mirna, "lin": (2.0 ** mirna.astype(float) - 1.0).clip(lower=0.0)}

    # predictor set: base + force(currency × α) × abundance
    specs = {f"base_{ab}": (ab, None) for ab in mats_in}
    for tag, (w_src, d_src) in _CURRENCY.items():
        for a in _ALPHAS:
            wd = _weight_dict(edges, w_src, d_src, a)
            for ab in mats_in:
                specs[f"force_{tag}_a{a:g}_{ab}"] = (ab, wd)

    rows = []
    for name, (ab, wd) in specs.items():
        mat = _abundance_sum_pressure(edges, mats_in[ab], pressure_genes, weights=wd)
        rows += _score_predictor(name, mat, rna, cov_base, target_cn)
    df = pd.DataFrame(rows)
    for name in specs:
        m = df["predictor"] == name
        if m.any():
            df.loc[m, "partial_q_cn"] = S.bh_fdr(df.loc[m, "partial_p_cn"].fillna(1.0).values)
    df["neg_fdr_cn"] = (df["partial_rho_cn"] < 0) & (df["partial_q_cn"] < 0.05)

    conc = build(genes=genes, out_dir=out_dir)[["gene", "n_regulators", "n_eff_base_log"]]
    df = df.merge(conc, on="gene", how="left")
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "gene_force_coupling.tsv"
    df.to_csv(path, sep="\t", index=False)

    # summary: net-repression (neg-FDR) per predictor, overall + multi-driver tail (base n_eff>=2)
    print("\n===== FORCE vs BASE coupling (TCGA tumor; net repression = neg-FDR partial-ρ) =====")
    print(f"{'predictor':22s} {'n':>5s} {'negFDR(all)':>11s}  {'tail n':>6s} {'negFDR(neff>=2)':>15s}")
    tail_genes = set(conc.loc[conc["n_eff_base_log"] >= 2, "gene"])
    for name in specs:
        d = df[df.predictor == name]
        if d.empty:
            continue
        t = d[d.gene.isin(tail_genes)]
        print(f"{name:22s} {len(d):>5d} {int(d['neg_fdr_cn'].sum()):>11d}  "
              f"{len(t):>6d} {int(t['neg_fdr_cn'].sum()):>15d}")
    print(f"[force-coupling] wrote {path}")
    return df


# ---------------------------------------------------------------------------- #
# Step 3 — tumor↔GTEx Δρ (the change). Carries {base, validated-α0} only.
# ---------------------------------------------------------------------------- #
def _aggregate_pressure(edges, mir, genes, weights):
    from mirna_hallmark.coupling_predictor_comparison import _abundance_sum_pressure
    return _abundance_sum_pressure(edges, mir, genes, weights=weights)


def tumor_gtex_delta_rho(genes: Sequence[str] | None = None, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    """Per-gene Δρ = ρ_tumor − ρ_GTEx for base (Σa) vs validated-α0 force (Σ a·log1p(ev)).
    Symmetric composition+proliferation confounds in BOTH cohorts (the comparable model —
    tumor-only purity/HRD/CN have no GTEx analog). Tests the curated ± gene sets."""
    from analysis.utils.common.loaders import partial_spearman
    from mirna_hallmark.analyses.cross_state.cross_state_deep_dive import _state_bundles
    from mirna_hallmark.analyses.cross_state.cross_state_coupling import (
        _state_covariates, _states_mirna_for_pressure, _mirna_for_edge_pressure)

    hs = HallmarkSets.load()
    edges = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    prolif = set(hs.sets.get("HALLMARK_E2F_TARGETS", [])) | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))
    bundles = _state_bundles(sorted(set(edges["gene"]) | prolif))
    states = [s for s in ("tumor", "gtex") if s in bundles]
    cov = {s: _state_covariates(s, bundles[s][1], hs, "metagene") for s in states}
    mir_for = _mirna_for_edge_pressure(bundles, _states_mirna_for_pressure(impute_gtex=True))

    pgenes = sorted(set(edges["gene"]).intersection(*[set(bundles[s][1].index) for s in states]))
    Wval = _weight_dict(edges, "validated", "validated", 0.0)  # = log1p(evidence), no promiscuity

    # rectified ACQUIRED abundance: max(a − h_healthy, 0), h = QN healthy baseline on TCGA scale
    from mirna_hallmark.healthy_anchor import gtex_qn_baseline
    h = gtex_qn_baseline()
    def _acq(mir):
        return mir.sub(h.reindex(mir.index), axis=0).clip(lower=0.0)

    # predictor = (matrix transform, edge weights)
    preds = {"base": (lambda m: m, None),
             "force_val_a0": (lambda m: m, Wval),
             "acq": (_acq, None)}

    rows = []
    for s in states:
        rna = bundles[s][1]
        for pname, (tf, wd) in preds.items():
            mat = _aggregate_pressure(edges, tf(mir_for[s]), pgenes, wd)
            for g in mat.index:
                if g not in rna.index:
                    continue
                rho, p, n = partial_spearman(rna.loc[g], pd.to_numeric(mat.loc[g], errors="coerce"), cov[s])
                rows.append({"gene": g, "predictor": pname, "state": s, "rho": rho, "n": n})
    long = pd.DataFrame(rows)
    wide = long.pivot_table(index=["gene", "predictor"], columns="state", values="rho").reset_index()
    wide["delta"] = wide["tumor"] - wide.get("gtex", np.nan)   # ρ_tumor − ρ_GTEx (acquired repression < 0)

    # join ± labels
    pos = set(pd.read_csv(C.REPO_ROOT / "mirna_hallmark/annotations/breast_oncomir_tsg_positive_set.tsv",
                          sep="\t")["gene"])
    conc = build(genes=genes, out_dir=out_dir)[["gene", "n_regulators"]]
    neg = set(conc.loc[conc["n_regulators"] <= 1, "gene"]) - pos
    wide["label"] = np.where(wide["gene"].isin(pos), "positive",
                             np.where(wide["gene"].isin(neg), "negative", "other"))
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "gene_force_delta_rho.tsv"
    wide.to_csv(path, sep="\t", index=False)

    from scipy.stats import mannwhitneyu
    print(f"\n===== tumor↔GTEx Δρ (ρ_tumor − ρ_GTEx; acquired repression ⇒ Δρ<0) =====")
    print(f"states={states}  genes={wide['gene'].nunique()}  positives={len(pos & set(wide.gene))}  negatives={len(neg & set(wide.gene))}")
    print(f"{'predictor':14s} {'set':9s} {'n':>4s} {'med ρ_tum':>9s} {'med ρ_gtex':>10s} {'med Δρ':>8s}")
    for pname in preds:
        d = wide[wide.predictor == pname]
        for lab in ("positive", "negative"):
            dl = d[d.label == lab].dropna(subset=["delta"])
            print(f"{pname:14s} {lab:9s} {len(dl):>4d} {dl['tumor'].median():>9.3f} "
                  f"{dl.get('gtex').median():>10.3f} {dl['delta'].median():>8.3f}")
        dp = d[d.label == "positive"].dropna(subset=["delta"])["delta"]
        dn = d[d.label == "negative"].dropna(subset=["delta"])["delta"]
        if len(dp) >= 5 and len(dn) >= 5:
            u, pu = mannwhitneyu(dp, dn, alternative="less")  # positives more negative Δρ?
            print(f"  -> {pname}: Δρ(pos) < Δρ(neg) Mann-Whitney p={pu:.3g} (sep = {dp.median()-dn.median():+.3f})")
    print(f"[delta-rho] wrote {path}")
    return wide


def tumor_acquired_coupling(genes: Sequence[str] | None = None, out_dir: Path = OUT_DIR,
                            add_prolif: bool = True, expressed_only: bool = False) -> pd.DataFrame:
    """Step 3b — score by ρ_tumor (acquired folds in the healthy contrast, so no Δρ). Predictors
    {base=Σa, force_val_a0=Σ a·log1p(ev), acq=Σ max(a−h,0), acq_ev=Σ max(a−h,0)·log1p(ev)}; canonical
    confounds CPE+HRD+batch+target-CN **+ proliferation** (the MH-78 gate). Four independent reads that
    `acq` beats `base`: (1) ±-set MWU, (2) PAIRED per-gene Wilcoxon, (3) genome-wide net-repression count,
    (4) bootstrap CI of the paired margin + an orthogonal gene_roles TSG label."""
    from mirna_hallmark.coupling_predictor_comparison import _abundance_sum_pressure, _load_grid_inputs
    from mirna_hallmark.healthy_anchor import gtex_qn_baseline
    from mirna_hallmark.analyses.cross_state.cross_state_coupling import _prolif_metagene
    from mirna_hallmark import gene_roles
    from scipy.stats import mannwhitneyu, wilcoxon

    hs, edges, mirna, rna, cov_base, _conf, pgenes, target_cn = _load_grid_inputs(genes)
    if expressed_only:                                    # canonical detectability floor (arm_expression)
        from mirna_hallmark.arm_expression import filter_expressed_edges
        n0 = len(edges); edges = filter_expressed_edges(edges)
        print(f"[acq-tumor] expressed-arm floor: {len(edges)}/{n0} edges kept")
    if add_prolif:
        cov_base = cov_base.join(_prolif_metagene(rna, hs).rename("prolif"), how="left")
    h = gtex_qn_baseline()
    acq_mir = mirna.sub(h.reindex(mirna.index), axis=0).clip(lower=0.0)
    _id = lambda m: m
    _acq = lambda m: acq_mir
    Wv0 = _weight_dict(edges, "validated", "validated", 0.0)
    Wv1 = _weight_dict(edges, "validated", "validated", 1.0)
    Wt0 = _weight_dict(edges, "targetscan", "targetscan", 0.0)
    Wt1 = _weight_dict(edges, "targetscan", "targetscan", 1.0)
    # full head-to-head: abundance-form × weight currency × promiscuity (α=0 none / α=1 budget-split)
    specs = {
        "base":     (_id, None),    # Σ a
        "val_a0":   (_id, Wv0),     # Σ a·log1p(ev)                  validated weight, no promiscuity
        "val_a1":   (_id, Wv1),     # Σ a·log1p(ev)/D_val            + budget-split promiscuity
        "ts_a0":    (_id, Wt0),     # Σ a·|context++|                TargetScan weight, no promiscuity
        "ts_a1":    (_id, Wt1),     # Σ a·|context++|/D_ts           + budget-split promiscuity
        "acq":      (_acq, None),   # Σ max(a−h,0)                   acquired abundance
        "acq_ev":   (_acq, Wv0),    # Σ max(a−h,0)·log1p(ev)         acquired × validated weight
        "acq_ts":   (_acq, Wt0),    # Σ max(a−h,0)·|context++|       acquired × TargetScan weight
    }
    mats = {name: _abundance_sum_pressure(edges, tf(mirna), pgenes, weights=wd) for name, (tf, wd) in specs.items()}
    rows = []
    for name, mat in mats.items():
        rows += _score_predictor(name, mat, rna, cov_base, target_cn)
    df = pd.DataFrame(rows)
    for name in mats:
        m = df["predictor"] == name
        if m.any():
            df.loc[m, "partial_q_cn"] = S.bh_fdr(df.loc[m, "partial_p_cn"].fillna(1.0).values)
    df["neg_fdr_cn"] = (df["partial_rho_cn"] < 0) & (df["partial_q_cn"] < 0.05)

    pos = set(pd.read_csv(C.REPO_ROOT / "mirna_hallmark/annotations/breast_oncomir_tsg_positive_set.tsv",
                          sep="\t")["gene"])
    conc = build(genes=genes, out_dir=out_dir)[["gene", "n_regulators"]]
    neg = set(conc.loc[conc["n_regulators"] <= 1, "gene"]) - pos
    roles = gene_roles.load_gene_roles(pgenes)
    tsg = set(roles.loc[roles["malignancy_sign"] == -1, "gene"]) - pos
    wide = df.pivot_table(index="gene", columns="predictor", values="partial_rho_cn")
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "gene_acquired_tumor_coupling.tsv", sep="\t", index=False)

    cmodel = "CPE+HRD+batch+target-CN" + ("+prolif" if add_prolif else "")
    print(f"\n===== ACQUIRED tumor coupling (ρ_tumor; confounds={cmodel}) =====")
    print(f"{'predictor':14s} {'gwNegFDR':>8s} {'med ρ pos':>9s} {'pos negFDR':>10s} {'med ρ neg':>9s} "
          f"{'pos<neg MWU':>11s} {'sep':>7s}")
    for name in mats:
        d = df[df.predictor == name]
        dp = d[d.gene.isin(pos)].dropna(subset=["partial_rho_cn"])
        dn = d[d.gene.isin(neg)].dropna(subset=["partial_rho_cn"])
        mwu = mannwhitneyu(dp["partial_rho_cn"], dn["partial_rho_cn"], alternative="less")[1] if min(len(dp), len(dn)) >= 5 else float("nan")
        print(f"{name:14s} {int(d['neg_fdr_cn'].sum()):>8d} {dp['partial_rho_cn'].median():>9.3f} "
              f"{int(dp['neg_fdr_cn'].sum()):>10d} {dn['partial_rho_cn'].median():>9.3f} "
              f"{mwu:>11.4g} {dp['partial_rho_cn'].median()-dn['partial_rho_cn'].median():>7.3f}")

    # PAIRED per-gene Δρ = ρ_pred − ρ_base for EVERY construction (the clean head-to-head;
    # * = bootstrap 95%CI of the median margin excludes 0 ⇒ beats abundance)
    print("\n  paired per-gene Δρ = ρ_pred − ρ_base (negative ⇒ beats abundance):  * CI excludes 0")
    rng = np.random.default_rng(0)
    for gname, gset in (("positives(18)", pos), (f"TSG(gene_roles)", tsg), ("all genes", set(wide.index))):
        print(f"   [{gname}]")
        for name in mats:
            if name == "base" or name not in wide.columns:
                continue
            w = wide.loc[wide.index.isin(gset), ["base", name]].dropna()
            if len(w) < 6:
                continue
            dpair = (w[name] - w["base"]).to_numpy()
            wp = wilcoxon(dpair, alternative="less").pvalue if np.any(dpair != 0) else float("nan")
            boot = np.array([np.median(rng.choice(dpair, len(dpair), replace=True)) for _ in range(2000)])
            star = "*" if np.percentile(boot, 97.5) < 0 else " "
            print(f"     {name:8s} n={len(w):4d} medΔ={np.median(dpair):+.4f} "
                  f"95%CI[{np.percentile(boot,2.5):+.4f},{np.percentile(boot,97.5):+.4f}] Wilcoxon p={wp:.2g} {star}")
    print(f"[acq-tumor] wrote {out_dir / 'gene_acquired_tumor_coupling.tsv'}")
    return df


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--concentration-only", action="store_true", help="only step 1 concentration")
    ap.add_argument("--delta-rho", action="store_true", help="step 3: tumor↔GTEx Δρ")
    ap.add_argument("--acq-tumor", action="store_true", help="step 3b: acquired coupling by ρ_tumor")
    args = ap.parse_args()
    if args.concentration_only:
        build()
    elif args.delta_rho:
        tumor_gtex_delta_rho()
    elif args.acq_tumor:
        tumor_acquired_coupling()
    else:
        force_vs_base_coupling()


if __name__ == "__main__":
    main()
