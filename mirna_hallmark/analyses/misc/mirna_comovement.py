"""Stage 2 (extension 2a): miRNA / target co-movement + pressure programs.

Builds on the Stage-1 ``mirna_state_class`` outputs (per-edge cross-state
contributions ``c_{gtex,nat,tumor}`` + the per-miRNA trajectory table). Three
complementary views of how shifts are *shared* across regulators and targets:

1. **Co-moving miRNAs** (``comovement_modules.tsv``): arm x arm Spearman
   co-expression across tumor samples -> hierarchical modules; each arm's top
   co-expressed partners + the module's dominant trajectory / subtype class. Plus a
   clustering of the 3-state trajectory vectors ``(dHN, dNT)`` so arms that follow
   the same healthy->NAT->tumor path are grouped.

2. **Target-set co-repression vs a neutral reference** (``target_corepression.tsv``):
   per miRNA, a composite repression score = the mean within-state z of its target
   set (low = the whole module is down). Its composition+prolif-adjusted partial
   Spearman against the miRNA's abundance, per state -> tests whether the *entire*
   target module moves down together (realised module pressure) beyond any single
   gene. The cross-state delta (tumor - gtex) is the acquired module repression.

3. **NMF pressure programs** (``nmf_*``): NMF on the non-negative cohort-aggregate
   miRNA x gene pressure-contribution matrix ``c_tumor`` -> latent programs (a
   co-regulating miRNA module + its shared *target* module). Two arms co-load mostly
   because they share targets; there is no patient axis. Factor activity is recomputed
   per state to show which programs the tumor amplifies. Optional second NMF on the
   repressive partial-coupling magnitude ``max(-adj_rho, 0)``.

4. **Within-tumor (per-sample) NMF** (``nmf_within_tumor_*``): complementary to view 3.
   NMF on the arm x tumor-sample expression matrix (log2 RPM+1, non-negative) -> programs
   of miRNAs that actually **co-vary across tumor patients** (heterogeneity the aggregate
   matrix cannot see). Each factor gets a per-patient activity vector summarised per PAM50
   subtype (mean activity, fraction, Kruskal-Wallis) -> subtype-specific programs, and is
   labelled with its top arms (+ state class) and dominant target module (sum of
   ``c_tumor`` over the factor's top arms, to keep the link to *pressure*).

6. **Gene-pressure per-sample NMF** (``nmf_gene_pressure_*``): NMF on the gene x tumor-sample
   matrix of **total miRNA pressure received by each gene per patient** (sum over all arms via
   ``compute_gene_pressure``). Factors = co-pressured gene programs with per-PAM50 activity
   (mean/frac + Kruskal-Wallis), each labelled with top genes, dominant Hallmark set, and the
   top driver arms (most ``c_tumor`` onto the factor's genes). Gene-side analog of view 4. The
   same decomposition is also run **within each PAM50 subtype** (columns sliced;
   ``nmf_gene_pressure_within_subtype``) -> subtype-private co-pressured gene programs.

5. **Per-subtype NMF** (``nmf_subtype_*``): NMF run *inside a single PAM50 subtype*, in two
   flavours. (a) ``nmf_subtype_persample`` = per-sample NMF on only that subtype's tumor
   columns -> programs co-varying within e.g. Basal alone. (b) ``nmf_subtype_static`` =
   NMF on the aggregate miRNA x gene contribution matrix whose **cell values c(m,g) are
   recomputed from only that subtype's samples** (the softmax abundance multiplier is
   averaged over the subtype) - there is *no* sample axis; the subtype enters through how
   each cell is computed, exactly as asked. The two answer different questions: (a)
   within-subtype heterogeneity, (b) the subtype's target-sharing pressure structure.

Run:
  .venv/bin/python3 -m mirna_hallmark.analyses.misc.mirna_comovement
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_coupling import _METAGENE_GENES, _state_covariates
from mirna_hallmark.analyses.cross_state.cross_state_deep_dive import _state_bundles
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.mirna_state_class import OUT_DIR as STATE_DIR

OUT_DIR = C.TISSUE_REFERENCE_DIR / "mirna_comovement"
MIN_EXPR = 1.0          # tumor median log2(RPM+1) to enter co-expression clustering
MIN_TARGETS = 5         # min target genes for a composite repression score
N_MODULES = 25          # agglomerative co-expression modules
NMF_K = 12              # latent pressure programs
TOP_LOAD = 15           # top miRNAs / genes reported per factor


def _stable_factor_label(row: dict) -> str:
    """Content-derived, run-stable handle for an NMF factor.

    NMF factor *indices* (F1..Fk / P.. / S.. / D..) are permutation-arbitrary: sklearn
    returns components in an order that depends on the input matrix, so they reshuffle
    whenever the arm/gene universe changes (e.g. the MIMAT expansion). This label keys a
    factor on its *content* — dominant PAM50 subtype + dominant Hallmark (or lead arm for
    arm-space NMFs) — so the same program is identifiable across rebuilds regardless of its
    F-index. The index is kept on the `factor` column for joins to the loadings tables.
    """
    sub = row.get("dominant_subtype") or "NA"
    hm = row.get("dominant_hallmark")
    if hm:
        return f"{sub}:{str(hm).replace('HALLMARK_', '')}"
    lead = ""
    for key in ("top_mirnas", "top_gain_genes", "top_genes"):
        v = row.get(key)
        if v:
            lead = str(v).split(";")[0].replace("hsa-miR-", "miR-").replace("hsa-let-", "let-")
            break
    return f"{sub}:{lead}" if lead else sub


# --------------------------------------------------------------------------- #
# Inputs: Stage-1 tables
# --------------------------------------------------------------------------- #
def _load_state_class(state_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mir = pd.read_csv(state_dir / "mirna_state_class.tsv", sep="\t")
    edge = pd.read_csv(state_dir / "mirna_gene_edge_class.tsv", sep="\t")
    return mir, edge


# --------------------------------------------------------------------------- #
# 1. Co-moving miRNAs
# --------------------------------------------------------------------------- #
def _comovement_modules(
    tumor_mir: pd.DataFrame, mir_class: pd.DataFrame, *, n_modules: int = N_MODULES,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    arms = [a for a in mir_class["miRNA"] if a in tumor_mir.index]
    expr = tumor_mir.loc[arms].apply(pd.to_numeric, errors="coerce")
    keep = expr.index[expr.median(axis=1) >= MIN_EXPR]
    expr = expr.loc[keep]
    if expr.shape[0] < 5:
        return pd.DataFrame(), pd.DataFrame()
    corr = expr.T.corr(method="spearman")
    dist = 1.0 - corr.fillna(0.0).values
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2.0
    Z = linkage(squareform(dist, checks=False), method="average")
    k = min(n_modules, max(2, expr.shape[0] // 4))
    labels = fcluster(Z, t=k, criterion="maxclust")
    mod = pd.Series(labels, index=expr.index, name="comove_module")

    cls = mir_class.set_index("miRNA")
    rows = []
    for arm in expr.index:
        c = corr[arm].drop(arm).sort_values(ascending=False)
        partners = [f"{p}({c[p]:.2f})" for p in c.head(5).index]
        info = cls.loc[arm] if arm in cls.index else {}
        rows.append({
            "miRNA": arm,
            "comove_module": int(mod[arm]),
            "primary_class": info.get("primary_class") if hasattr(info, "get") else None,
            "secondary_class": info.get("secondary_class") if hasattr(info, "get") else None,
            "dHT": info.get("dHT") if hasattr(info, "get") else np.nan,
            "top_coexpressed": ";".join(partners),
        })
    modules = pd.DataFrame(rows)

    summ = (modules.groupby("comove_module")
            .agg(n_arms=("miRNA", "size"),
                 mean_dHT=("dHT", "mean"),
                 dominant_primary=("primary_class", lambda s: s.value_counts().idxmax() if s.notna().any() else None),
                 dominant_secondary=("secondary_class", lambda s: s.value_counts().idxmax() if s.notna().any() else None),
                 members=("miRNA", lambda s: ";".join(sorted(s)[:12])))
            .reset_index().sort_values("n_arms", ascending=False))
    return modules, summ


def _trajectory_clusters(mir_class: pd.DataFrame, *, n_clusters: int = 8) -> pd.Series:
    from sklearn.cluster import AgglomerativeClustering

    v = mir_class.dropna(subset=["dHN", "dNT"]).set_index("miRNA")[["dHN", "dNT"]]
    if v.shape[0] < n_clusters:
        return pd.Series(dtype=int)
    lab = AgglomerativeClustering(n_clusters=n_clusters).fit_predict(v.values)
    return pd.Series(lab, index=v.index, name="trajectory_cluster")


# --------------------------------------------------------------------------- #
# 2. Target-set co-repression vs neutral reference
# --------------------------------------------------------------------------- #
def _composite_repression(rna: pd.DataFrame, genes: Sequence[str],
                          weights: Optional[Dict[str, float]] = None) -> pd.Series:
    """Per-sample composite of the target set = (optionally edge-weighted) mean of each target
    gene's WITHIN-STATE z-score. The z is per-GENE across samples (z(g,s)=(x−mean_g)/sd_g), so it
    standardises every target onto a comparable scale before averaging; low composite = the module
    is coordinately down. `weights` (per-gene w(m,g), e.g. evidence_score) → a weighted mean so a
    strong edge counts more than a weak one (the equal-weight mean is the default/legacy)."""
    present = [g for g in genes if g in rna.index]
    if len(present) < MIN_TARGETS:
        return pd.Series(dtype=float)
    sub = rna.loc[present].apply(pd.to_numeric, errors="coerce")
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    if weights:
        w = pd.Series({g: float(weights.get(g, 0.0)) for g in present}).reindex(z.index).fillna(0.0)
        if w.sum() > 0:
            return z.mul(w, axis=0).sum(axis=0) / w.sum()
    return z.mean(axis=0)  # equal-weight per-sample composite (low = target module repressed)


def _target_corepression(
    edge: pd.DataFrame, hs: HallmarkSets, *, cov_source: str = "estimate_epi",
) -> pd.DataFrame:
    from analysis.utils.common.loaders import partial_spearman

    targets = edge.groupby("miRNA")["gene"].apply(lambda s: sorted(set(s)))
    arms = [a for a, gs in targets.items() if len(gs) >= MIN_TARGETS]
    # per-(arm,gene) edge weight w(m,g) for the WEIGHTED composite (evidence_score = static
    # mirTarBase/TargetScan evidence mass; fall back to c_tumor; else equal weight)
    wcol = "evidence_score" if "evidence_score" in edge.columns else (
        "c_tumor" if "c_tumor" in edge.columns else None)
    edge_w: Dict[str, Dict[str, float]] = {}
    if wcol is not None:
        for a, grp in edge.groupby("miRNA"):
            edge_w[a] = dict(zip(grp["gene"], pd.to_numeric(grp[wcol], errors="coerce").fillna(0.0)))
    # gene-role map (oncogene / TSG) for the ROLE-STRATIFIED composite. We deliberately do NOT
    # fold role into one weighted composite (that would conflate "is the module repressed" with
    # "is the repression pro/anti-tumor"); instead we compute the composite SEPARATELY over the
    # arm's TSG targets vs its oncogene targets, keeping the two questions cleanly apart.
    try:
        from mirna_hallmark.gene_roles import load_gene_roles
        _roles = load_gene_roles().set_index("gene")["malignancy_sign"].to_dict()
    except Exception:
        _roles = {}
    genes_needed = sorted(set().union(*[set(targets[a]) for a in arms])
                          | set(_METAGENE_GENES)
                          | set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
                          | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", [])))
    bundles = _state_bundles(genes_needed)
    states = list(bundles)
    cov: Dict[str, pd.DataFrame] = {}
    for s, (_mir, rna) in bundles.items():
        try:
            cov[s] = _state_covariates(s, rna, hs, cov_source)
        except Exception:
            cov[s] = _state_covariates(s, rna, hs, "metagene")

    rows: List[dict] = []
    ROLE_MIN = 4  # need at least this many role-annotated targets to form a stratified composite
    for i, arm in enumerate(arms):
        row = {"miRNA": arm, "n_targets": len(targets[arm])}
        tsg_targets = [g for g in targets[arm] if _roles.get(g, 0) < 0]
        onco_targets = [g for g in targets[arm] if _roles.get(g, 0) > 0]
        row["n_tsg_targets"] = len(tsg_targets)
        row["n_onco_targets"] = len(onco_targets)
        for s in states:
            mir, rna = bundles[s]
            if arm not in mir.index:
                row[f"rho_composite_{s}"] = np.nan
                row[f"q_composite_{s}"] = np.nan
                row[f"n_{s}"] = 0
                continue
            comp = _composite_repression(rna, targets[arm])
            comp_w = _composite_repression(rna, targets[arm], weights=edge_w.get(arm))
            x = pd.to_numeric(mir.loc[arm], errors="coerce")
            shared = sorted(set(comp.dropna().index) & set(x.dropna().index))
            if len(shared) < 25:
                row[f"rho_composite_{s}"] = np.nan
                row[f"rho_composite_w_{s}"] = np.nan
                row[f"q_composite_{s}"] = np.nan
                row[f"n_{s}"] = len(shared)
                continue
            adj_rho, adj_p, n = partial_spearman(comp.loc[shared], x.loc[shared], cov[s].reindex(shared))
            row[f"rho_composite_{s}"] = adj_rho
            row[f"q_composite_{s}"] = adj_p
            row[f"n_{s}"] = n
            # edge-WEIGHTED composite (strong edges count more) — parallel column
            if not comp_w.empty:
                sh_w = sorted(set(comp_w.dropna().index) & set(x.dropna().index))
                if len(sh_w) >= 25:
                    wr, _wp, _wn = partial_spearman(comp_w.loc[sh_w], x.loc[sh_w], cov[s].reindex(sh_w))
                    row[f"rho_composite_w_{s}"] = wr
                else:
                    row[f"rho_composite_w_{s}"] = np.nan
            else:
                row[f"rho_composite_w_{s}"] = np.nan
            # ROLE-STRATIFIED composites (clean separation, not a reweight): is the coordinated
            # repression landing on tumor-SUPPRESSOR targets (anti-tumor realized) or on
            # ONCOGENE targets (pro-tumor realized)?
            for tag, tgts in (("tsg", tsg_targets), ("onco", onco_targets)):
                rc = np.nan
                if len(tgts) >= ROLE_MIN:
                    cmp_r = _composite_repression(rna, tgts)
                    if not cmp_r.empty:
                        sh_r = sorted(set(cmp_r.dropna().index) & set(x.dropna().index))
                        if len(sh_r) >= 25:
                            rc, _rp, _rn = partial_spearman(cmp_r.loc[sh_r], x.loc[sh_r],
                                                            cov[s].reindex(sh_r))
                        row[f"rho_composite_{tag}_{s}"] = rc
                        continue
                row[f"rho_composite_{tag}_{s}"] = rc
        rows.append(row)
        if (i + 1) % 200 == 0:
            print(f"[comovement] target-corepression {i + 1}/{len(arms)}")
    df = pd.DataFrame(rows)
    from mirna_hallmark.stats import bh_fdr
    for s in states:
        if f"q_composite_{s}" in df.columns:
            df[f"q_composite_{s}"] = bh_fdr(df[f"q_composite_{s}"].fillna(1.0).values)
    if {"rho_composite_tumor", "rho_composite_gtex"}.issubset(df.columns):
        df["delta_tumor_gtex"] = df["rho_composite_tumor"] - df["rho_composite_gtex"]
    return df.sort_values("rho_composite_tumor").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# 2b. Gene-level realized net-repression (gene-side mirror of target_corepression)
# --------------------------------------------------------------------------- #
def gene_corepression(
    hs: Optional[HallmarkSets] = None, *, cov_source: str = "estimate_epi",
    out_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Is a gene GENUINELY under net miRNA repression? For each Hallmark target gene, build its
    aggregate incoming miRNA pressure per sample (`compute_gene_pressure`, sum over ALL its
    regulators) and take the composition+proliferation-adjusted partial Spearman of the gene's
    OWN expression vs that aggregate pressure, per state. ``rho_gene_pressure_{state} < 0`` = the
    gene's total miRNA pressure realizedly represses it (more pressure ↔ lower expression).

    Trajectory class lets a (de)coupled edge be read against its gene's whole-regulator state:
    `constitutive_repressed` (gene net-repressed in both reference and tumor → a single lost edge
    is *edge-local escape*, other regulators carry the gene), `lost_repression` (gene net-repressed
    in NAT/healthy but ≈0 in tumor → *gene-wide de-repression*), `gained_repression`, `never`."""
    from analysis.utils.common.loaders import partial_spearman

    from mirna_hallmark.analyses.cross_state.cross_state_landscape import ATTR_MODE
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    from mirna_hallmark.stats import bh_fdr

    if hs is None:
        hs = HallmarkSets.load()
    out_dir = out_dir or (C.TISSUE_REFERENCE_DIR / "mirna_comovement")
    out_dir.mkdir(parents=True, exist_ok=True)
    genes_all = sorted(hs.universe)
    edges = load_mirtar_edges(genes_all, resolve_arms=True)
    reg_genes = sorted(set(edges["gene"]) & set(genes_all))
    n_reg = edges.groupby("gene")["miRNA"].nunique()
    genes_needed = sorted(set(reg_genes) | set(_METAGENE_GENES)
                          | set(hs.sets.get("HALLMARK_E2F_TARGETS", []))
                          | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", [])))
    bundles = _state_bundles(genes_needed)
    states = list(bundles)
    cov: Dict[str, pd.DataFrame] = {}
    for s, (_mir, rna) in bundles.items():
        try:
            cov[s] = _state_covariates(s, rna, hs, cov_source)
        except Exception:
            cov[s] = _state_covariates(s, rna, hs, "metagene")
    # aggregate incoming pressure per state (gene x sample)
    Pcache: Dict[str, pd.DataFrame] = {}
    for s in states:
        mir, _rna = bundles[s]
        try:
            Pcache[s] = compute_gene_pressure(reg_genes, edges=edges, mirna=mir,  # type: ignore[arg-type]
                                              expr_mode=ATTR_MODE, resolve_arms=False)
        except Exception as exc:
            print(f"[comovement] gene_corepression: pressure failed for {s}: {exc}")
            Pcache[s] = pd.DataFrame()

    rows: List[dict] = []
    for gi, g in enumerate(reg_genes):
        row = {"gene": g, "n_regulators": int(n_reg.get(g, 0))}
        for s in states:
            _mir, rna = bundles[s]
            P = Pcache[s]
            rho = p = np.nan
            n = 0
            if (not P.empty) and g in rna.index and g in P.index:
                y = pd.to_numeric(rna.loc[g], errors="coerce")
                x = pd.to_numeric(P.loc[g], errors="coerce")
                shared = sorted(set(y.dropna().index) & set(x.dropna().index))
                if len(shared) >= 25 and x.loc[shared].std() > 0:
                    rho, p, n = partial_spearman(y.loc[shared], x.loc[shared], cov[s].reindex(shared))
            row[f"rho_gene_pressure_{s}"] = rho
            row[f"q_gene_pressure_{s}"] = p
            row[f"n_{s}"] = n
        rows.append(row)
        if (gi + 1) % 400 == 0:
            print(f"[comovement] gene-corepression {gi + 1}/{len(reg_genes)}")
    df = pd.DataFrame(rows)
    for s in states:
        if f"q_gene_pressure_{s}" in df.columns:
            df[f"q_gene_pressure_{s}"] = bh_fdr(df[f"q_gene_pressure_{s}"].fillna(1.0).values)
    if {"rho_gene_pressure_tumor", "rho_gene_pressure_gtex"}.issubset(df.columns):
        df["delta_tumor_gtex"] = df["rho_gene_pressure_tumor"] - df["rho_gene_pressure_gtex"]
    if {"rho_gene_pressure_tumor", "rho_gene_pressure_nat"}.issubset(df.columns):
        df["delta_tumor_nat"] = df["rho_gene_pressure_tumor"] - df["rho_gene_pressure_nat"]

    def _cls(r: pd.Series) -> str:
        t = r.get("rho_gene_pressure_tumor")
        ref = r.get("rho_gene_pressure_nat")
        if pd.isna(ref):
            ref = r.get("rho_gene_pressure_gtex")
        if pd.isna(t):
            return "na"
        rep_t = t < -0.10
        rep_ref = pd.notna(ref) and ref < -0.10
        if rep_t and rep_ref:
            return "constitutive_repressed"
        if rep_t and not rep_ref:
            return "gained_repression"
        if (not rep_t) and rep_ref:
            return "lost_repression"
        return "never_repressed"

    df["gene_repression_class"] = df.apply(_cls, axis=1)
    df["gene_net_repressed_tumor"] = ((df.get("rho_gene_pressure_tumor", np.nan) < -0.10)
                                      & (df.get("q_gene_pressure_tumor", 1.0) < 0.10))
    df = df.sort_values("rho_gene_pressure_tumor").reset_index(drop=True)
    df.to_csv(out_dir / "gene_corepression.tsv", sep="\t", index=False)
    vc = df["gene_repression_class"].value_counts().to_dict()
    print(f"[comovement] gene_corepression: {len(df)} genes; classes={vc}; "
          f"net-repressed in tumor={int(df['gene_net_repressed_tumor'].sum())}")
    try:
        edge_gene_repression_scope(df, out_dir=out_dir)
    except Exception as exc:
        print(f"[comovement] WARN edge_gene_repression_scope skipped ({exc})")
    return df


def edge_gene_repression_scope(
    gcr: Optional[pd.DataFrame] = None, *, out_dir: Optional[Path] = None,
    state_dir: Path = STATE_DIR,
) -> pd.DataFrame:
    """Direction-aware edge-level scope: read each edge's `joint_edge_class` against its target's
    gene-level realized-repression state (`gene_corepression`). Covers BOTH directions in one place
    (the de-coupling/escape view AND the symmetric acquired/gained view), so a gained coupling can be
    read as gene-wide acquired repression vs an edge-specific acquired coupling."""
    out_dir = out_dir or (C.TISSUE_REFERENCE_DIR / "mirna_comovement")
    if gcr is None:
        gcr = pd.read_csv(out_dir / "gene_corepression.tsv", sep="\t")
    edge = pd.read_csv(state_dir / "mirna_gene_edge_class.tsv", sep="\t")
    g = gcr.set_index("gene")
    keep = ["miRNA", "gene", "joint_edge_class"]
    keep += [c for c in ("gene_share_tumor", "gene_share_rank_tumor", "n_regulators_gene")
             if c in edge.columns]
    out = edge[keep].copy()
    out["gene_repression_class"] = out["gene"].map(g["gene_repression_class"]).fillna("na")
    out["gene_net_repressed_tumor"] = out["gene"].map(
        g["gene_net_repressed_tumor"]).fillna(False).astype(bool)

    def _scope(r: pd.Series) -> str:
        jec = str(r["joint_edge_class"])
        cls = str(r["gene_repression_class"])
        rep = bool(r["gene_net_repressed_tumor"])
        if jec in ("lost", "nat_decoupled"):
            if rep:
                return "edge_local_escape"        # gene still net-repressed by other regulators
            if cls == "lost_repression":
                return "gene_wide_derepression"    # gene lost net miRNA control
            return "gene_not_repressed"
        if jec in ("acquired_realized", "field_established_realized", "acquired_unrealized"):
            if cls == "gained_repression":
                return "gene_wide_acquired_repression"   # gene newly under net repression
            if cls == "constitutive_repressed":
                return "reinforces_constitutive"          # adds to pre-existing gene repression
            return "edge_specific_acquired"               # gene not broadly (newly) repressed
        if jec == "constitutive":
            return "constitutive_concordant" if cls == "constitutive_repressed" else "constitutive_edge_only"
        return f"other:{cls}"

    out["gene_repression_scope"] = out.apply(_scope, axis=1)
    out.to_csv(out_dir / "edge_gene_repression_scope.tsv", sep="\t", index=False)
    # report the two directional breakdowns
    lost = out[out["joint_edge_class"].isin(["lost", "nat_decoupled"])]
    acq = out[out["joint_edge_class"].isin(
        ["acquired_realized", "field_established_realized", "acquired_unrealized"]
    )]
    print(f"[comovement] edge_gene_repression_scope: lost={lost['gene_repression_scope'].value_counts().to_dict()}")
    print(f"[comovement] edge_gene_repression_scope: acquired={acq['gene_repression_scope'].value_counts().to_dict()}")
    return out


# --------------------------------------------------------------------------- #
# 3. NMF pressure programs
# --------------------------------------------------------------------------- #
def _nmf_share(edge: pd.DataFrame, *, k: int = NMF_K) -> Dict[str, pd.DataFrame]:
    """NMF on the *regulatory-share* matrices rather than raw pressure. `c_tumor` NMF (`_nmf_programs`)
    is dominated by the abundance backbone; normalising removes that and exposes a complementary
    structure:
    - **dominant-regulator programs** — NMF on `gene_share_tumor` (= c / Σ_regulators c, column-
      normalised so every gene's incoming pressure sums to 1): surfaces arms that *own* their targets'
      budgets even when globally low-abundance (the backbone top factor stays, but secondaries diverge
      from the c_tumor NMF — miR-125b/24/27~29 dominant-regulator modules).
    - **grip-change (acquired-share) program** — signed NMF on `gene_share_tumor − gene_share_gtex`
      split into gain/loss channels: which arms *co-gain* vs *co-lose* regulatory grip across genes,
      independent of overall pressure inflation (the share analogue of the signed gene-pressure NMF)."""
    from sklearn.decomposition import NMF
    out: Dict[str, pd.DataFrame] = {}
    GS = _pivot(edge, "gene_share_tumor")
    if GS.empty or min(GS.shape) < k:
        return out
    m = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=500)
    W = pd.DataFrame(m.fit_transform(GS.values), index=GS.index, columns=[f"D{j+1}" for j in range(k)])
    H = pd.DataFrame(m.components_, index=[f"D{j+1}" for j in range(k)], columns=GS.columns)
    rows = []
    for f in W.columns:
        rows.append({"factor": f,
                     "top_mirnas": ";".join(W[f].sort_values(ascending=False).head(8).index),
                     "top_genes": ";".join(H.loc[f].sort_values(ascending=False).head(8).index)})
    out["nmf_share_dominant"] = pd.DataFrame(rows)

    if {"gene_share_tumor", "gene_share_gtex"}.issubset(edge.columns):
        d = edge.copy()
        d["d_share"] = (pd.to_numeric(d["gene_share_tumor"], errors="coerce")
                        - pd.to_numeric(d["gene_share_gtex"], errors="coerce"))
        d["gain"] = d["d_share"].clip(lower=0.0)
        d["loss"] = (-d["d_share"]).clip(lower=0.0)
        Gn = d.pivot_table(index="miRNA", columns="gene", values="gain", aggfunc="mean").fillna(0.0)
        Ln = d.pivot_table(index="miRNA", columns="gene", values="loss", aggfunc="mean").fillna(0.0)
        Ln.columns = [f"{c}|loss" for c in Ln.columns]
        Gn.columns = [f"{c}|gain" for c in Gn.columns]
        M = pd.concat([Gn, Ln], axis=1).fillna(0.0)
        if min(M.shape) >= k:
            ms = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=500)
            Ws = pd.DataFrame(ms.fit_transform(M.values), index=M.index, columns=[f"G{j+1}" for j in range(k)])
            Hs = pd.DataFrame(ms.components_, index=[f"G{j+1}" for j in range(k)], columns=M.columns)
            srows = []
            for f in Ws.columns:
                gmass = Hs.loc[f, [c for c in M.columns if c.endswith("|gain")]].sum()
                lmass = Hs.loc[f, [c for c in M.columns if c.endswith("|loss")]].sum()
                direction = "gain" if gmass >= lmass else "loss"
                top = Hs.loc[f].sort_values(ascending=False).head(8).index
                srows.append({"factor": f, "direction": direction,
                              "gain_mass_frac": float(gmass / (gmass + lmass + 1e-9)),
                              "top_mirnas": ";".join(Ws[f].sort_values(ascending=False).head(8).index),
                              "top_genes": ";".join(g.rsplit("|", 1)[0] for g in top)})
            out["nmf_share_grip_change"] = pd.DataFrame(srows)
    return out


def _pivot(edge: pd.DataFrame, value: str) -> pd.DataFrame:
    if value not in edge.columns:
        return pd.DataFrame()
    m = edge.pivot_table(index="miRNA", columns="gene", values=value, aggfunc="mean")
    return m.fillna(0.0).clip(lower=0.0)


def _nmf_programs(edge: pd.DataFrame, *, k: int = NMF_K) -> Dict[str, pd.DataFrame]:
    from sklearn.decomposition import NMF

    C_t = _pivot(edge, "c_tumor")
    if C_t.empty or min(C_t.shape) < k:
        return {}
    model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=500)
    W = model.fit_transform(C_t.values)   # miRNA x k
    H = model.components_                  # k x gene
    Wn = pd.DataFrame(W, index=C_t.index, columns=[f"F{j + 1}" for j in range(k)])
    Hn = pd.DataFrame(H, index=[f"F{j + 1}" for j in range(k)], columns=C_t.columns)

    mir_load, gene_load = [], []
    for f in Wn.columns:
        for arm in Wn[f].sort_values(ascending=False).head(TOP_LOAD).index:
            mir_load.append({"factor": f, "miRNA": arm, "loading": float(Wn.loc[arm, f])})
        for g in Hn.loc[f].sort_values(ascending=False).head(TOP_LOAD).index:
            gene_load.append({"factor": f, "gene": g, "loading": float(Hn.loc[f, g])})

    # per-state factor activity: weight each edge's state contribution by factor membership
    states = [s for s in ("gtex", "nat", "tumor") if f"c_{s}" in edge.columns]
    e = edge.set_index(["miRNA", "gene"])
    act_rows = []
    for f in Wn.columns:
        memb = {}
        for (arm, g) in e.index:
            if arm in Wn.index and g in Hn.columns:
                memb[(arm, g)] = Wn.loc[arm, f] * Hn.loc[f, g]
        memb_s = pd.Series(memb)
        row = {"factor": f}
        for s in states:
            cvals = e[f"c_{s}"].reindex(memb_s.index).fillna(0.0)
            row[f"activity_{s}"] = float((memb_s * cvals).sum())
        act_rows.append(row)
    act = pd.DataFrame(act_rows)
    for s in states:
        tot = act[f"activity_{s}"].sum()
        if tot > 0:
            act[f"frac_{s}"] = act[f"activity_{s}"] / tot
    if {"frac_tumor", "frac_gtex"}.issubset(act.columns):
        act["shift_tumor_gtex"] = act["frac_tumor"] - act["frac_gtex"]

    out = {"nmf_factor_loadings_mirna": pd.DataFrame(mir_load),
           "nmf_factor_loadings_gene": pd.DataFrame(gene_load),
           "nmf_factor_state_activity": act.sort_values("shift_tumor_gtex", ascending=False)
           if "shift_tumor_gtex" in act.columns else act}

    # optional: repressive partial-coupling NMF
    if "adj_rho_tumor" in edge.columns:
        edge_rep = edge.copy()
        edge_rep["rep"] = (-edge_rep["adj_rho_tumor"]).clip(lower=0.0)
        R = _pivot(edge_rep, "rep")
        if not R.empty and min(R.shape) >= k:
            rmod = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=500)
            RW = pd.DataFrame(rmod.fit_transform(R.values), index=R.index,
                              columns=[f"R{j + 1}" for j in range(k)])
            RH = pd.DataFrame(rmod.components_, index=[f"R{j + 1}" for j in range(k)], columns=R.columns)
            rep_rows = []
            for f in RW.columns:
                top_m = ";".join(RW[f].sort_values(ascending=False).head(8).index)
                top_g = ";".join(RH.loc[f].sort_values(ascending=False).head(8).index)
                rep_rows.append({"factor": f, "top_mirnas": top_m, "top_genes": top_g})
            out["nmf_repressive_coupling_programs"] = pd.DataFrame(rep_rows)
    return out


# --------------------------------------------------------------------------- #
# 4. Within-tumor (per-sample) NMF, PAM50-resolved
# --------------------------------------------------------------------------- #
def _nmf_within_tumor(
    tumor_mat: pd.DataFrame, edge: pd.DataFrame, mir_class: pd.DataFrame, *, k: int = NMF_K,
) -> Dict[str, pd.DataFrame]:
    """NMF on the arm x tumor-sample expression matrix -> programs of miRNAs that
    co-vary across patients, with per-PAM50 activity (Q3 subtype specifics)."""
    from sklearn.decomposition import NMF
    from scipy.stats import kruskal

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _PAM50
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.stats import bh_fdr

    expr = tumor_mat.apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[expr.median(axis=1) >= MIN_EXPR].dropna(how="any")
    if expr.shape[0] < k or expr.shape[1] < k:
        return {}
    X = expr.clip(lower=0.0)
    model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=1000)
    W = pd.DataFrame(model.fit_transform(X.values), index=X.index,
                     columns=[f"T{j + 1}" for j in range(k)])   # arm loadings
    H = pd.DataFrame(model.components_, index=[f"T{j + 1}" for j in range(k)],
                     columns=X.columns)                          # per-sample activity

    clin = D.load_clinical_strata()
    part_pam = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
                if "PAM50_final" in clin.columns else {})
    col_pam = pd.Series({c: part_pam.get(_participant(c)) for c in X.columns})

    cls = mir_class.set_index("miRNA")
    mir_rows: List[dict] = []
    for f in W.columns:
        for arm in W[f].sort_values(ascending=False).head(TOP_LOAD).index:
            info = cls.loc[arm] if arm in cls.index else {}
            mir_rows.append({
                "factor": f, "miRNA": arm, "loading": float(W.loc[arm, f]),
                "primary_class": info.get("primary_class") if hasattr(info, "get") else None,
                "secondary_class": info.get("secondary_class") if hasattr(info, "get") else None,
            })
    mir_load = pd.DataFrame(mir_rows)

    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None
    act_rows: List[dict] = []
    for f in H.index:
        a = H.loc[f].astype(float)
        row: dict = {"factor": f}
        means, groups = [], []
        for sub in _PAM50:
            cols = [c for c in col_pam.index[col_pam == sub] if c in a.index]
            vals = a.loc[cols].dropna()
            m = float(vals.mean()) if len(vals) else np.nan
            row[f"mean_act_{sub}"] = m
            means.append(m)
            if len(vals) >= 10:
                groups.append(vals.values)
        means_arr = np.array(means, float)
        if np.nansum(means_arr) > 0:
            for sub, m in zip(_PAM50, means_arr):
                row[f"frac_{sub}"] = (m / np.nansum(means_arr)) if np.isfinite(m) else np.nan
            row["dominant_subtype"] = _PAM50[int(np.nanargmax(means_arr))]
        else:
            row["dominant_subtype"] = None
        try:
            row["kw_p"] = float(kruskal(*groups)[1]) if len(groups) >= 2 else np.nan
        except ValueError:
            row["kw_p"] = np.nan
        top_arms = list(W[f].sort_values(ascending=False).head(8).index)
        row["top_mirnas"] = ";".join(top_arms)
        if ec is not None:
            sub_e = ec[ec.index.get_level_values("miRNA").isin(top_arms)]
            top_g = sub_e.groupby(level="gene").sum().sort_values(ascending=False).head(8)
            row["top_target_module"] = ";".join(top_g.index)
        row["factor_label"] = _stable_factor_label(row)
        act_rows.append(row)
    act = pd.DataFrame(act_rows)
    act["kw_q"] = bh_fdr(act["kw_p"].fillna(1.0).values)
    return {"nmf_within_tumor_loadings_mirna": mir_load,
            "nmf_within_tumor_factor_subtype": act.sort_values("kw_q")}


# --------------------------------------------------------------------------- #
# 6. NMF on the per-sample GENE-PRESSURE matrix (gene x tumor-sample)
# --------------------------------------------------------------------------- #
def _gene_factor_label(top_genes: Sequence[str], hs: HallmarkSets, ec) -> dict:
    """Dominant Hallmark set + top driver arms for a gene-program factor."""
    tg = set(top_genes)
    best, best_frac = None, 0.0
    for hset, members in hs.sets.items():
        frac = len(tg & set(members)) / max(1, len(tg))
        if frac > best_frac:
            best, best_frac = hset, frac
    row = {"top_genes": ";".join(list(top_genes)[:8]),
           "dominant_hallmark": best, "dominant_hallmark_frac": round(best_frac, 3)}
    if ec is not None:
        sub_e = ec[ec.index.get_level_values("gene").isin(top_genes)]
        row["top_driver_mirnas"] = ";".join(
            sub_e.groupby(level="miRNA").sum().sort_values(ascending=False).head(8).index)
    return row


def _nmf_gene_pressure(
    tumor_mat: pd.DataFrame, edge: pd.DataFrame, hs: HallmarkSets, *, k: int = NMF_K,
    expr_mode: Optional[str] = None, baseline: Optional[pd.Series] = None, tag: str = "",
) -> Dict[str, pd.DataFrame]:
    """NMF on the gene x tumor-sample TOTAL-pressure matrix (sum of every miRNA's
    contribution per gene per patient) -> programs of genes that are co-pressured across
    patients, with per-PAM50 activity. Gene-side analog of the within-tumor arm NMF.

    Also runs the SAME decomposition *within each PAM50 subtype* (columns sliced to that
    subtype; pressure is per-(gene,sample) so slicing is exact) -> subtype-private
    co-pressured gene programs.

    `expr_mode`/`baseline`/`tag`: default = abundance pressure (`ATTR_MODE`, softmax_logrpm).
    When `baseline` is given with the dev-anchored mode, the per-(gene,sample) pressure is the
    ACQUIRED deviation vs the healthy baseline (clipped at 0 below, since NMF needs non-negative
    -> 'acquired-only gain'); `tag='_acquired'` namespaces the outputs. This isolates programs of
    genes that co-GAIN pressure vs healthy across patients (vs the abundance co-variation view)."""
    from sklearn.decomposition import NMF
    from scipy.stats import kruskal

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import ATTR_MODE, _PAM50
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    from mirna_hallmark.stats import bh_fdr

    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    mode = expr_mode or ATTR_MODE
    kw = {"healthy_baseline": baseline} if baseline is not None else {}
    P = compute_gene_pressure(genes, edges=edges, mirna=tumor_mat, expr_mode=mode,  # type: ignore[arg-type]
                              resolve_arms=False, **kw)
    if P.empty:
        return {}
    P = P.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)
    P = P.loc[P.sum(axis=1) > 0]
    if min(P.shape) < k:
        return {}
    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None

    def _fit(mat: pd.DataFrame, prefix: str):
        m = mat.loc[mat.sum(axis=1) > 0]
        if min(m.shape) < k:
            return None, None
        model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=1000)
        W = pd.DataFrame(model.fit_transform(m.values), index=m.index,
                         columns=[f"{prefix}{j + 1}" for j in range(k)])
        H = pd.DataFrame(model.components_, index=W.columns, columns=m.columns)
        return W, H

    W, H = _fit(P, "P")
    if W is None:
        return {}

    clin = D.load_clinical_strata()
    part_pam = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
                if "PAM50_final" in clin.columns else {})
    col_pam = pd.Series({c: part_pam.get(_participant(c)) for c in P.columns})

    gene_rows: List[dict] = []
    for f in W.columns:
        for g in W[f].sort_values(ascending=False).head(TOP_LOAD).index:
            gene_rows.append({"factor": f, "gene": g, "loading": float(W.loc[g, f])})
    gene_load = pd.DataFrame(gene_rows)

    act_rows: List[dict] = []
    for f in H.index:
        a = H.loc[f].astype(float)
        row: dict = {"factor": f}
        means, groups = [], []
        for sub in _PAM50:
            cols = [c for c in col_pam.index[col_pam == sub] if c in a.index]
            vals = a.loc[cols].dropna()
            m = float(vals.mean()) if len(vals) else np.nan
            row[f"mean_act_{sub}"] = m
            means.append(m)
            if len(vals) >= 10:
                groups.append(vals.values)
        means_arr = np.array(means, float)
        if np.nansum(means_arr) > 0:
            for sub, m in zip(_PAM50, means_arr):
                row[f"frac_{sub}"] = (m / np.nansum(means_arr)) if np.isfinite(m) else np.nan
            row["dominant_subtype"] = _PAM50[int(np.nanargmax(means_arr))]
        else:
            row["dominant_subtype"] = None
        try:
            row["kw_p"] = float(kruskal(*groups)[1]) if len(groups) >= 2 else np.nan
        except ValueError:
            row["kw_p"] = np.nan
        top_genes = list(W[f].sort_values(ascending=False).head(30).index)
        row.update(_gene_factor_label(top_genes, hs, ec))
        row["factor_label"] = _stable_factor_label(row)
        act_rows.append(row)
    act = pd.DataFrame(act_rows)
    act["kw_q"] = bh_fdr(act["kw_p"].fillna(1.0).values)

    # ---- within-subtype gene-program NMF (columns sliced per PAM50) ----------- #
    sub_rows: List[dict] = []
    sub_load: List[dict] = []
    for sub in _PAM50:
        cols = [c for c in col_pam.index[col_pam == sub] if c in P.columns]
        if len(cols) < 30:
            continue
        Ws, _ = _fit(P[cols], f"{sub}_G")
        if Ws is None:
            continue
        for f in Ws.columns:
            top_genes = list(Ws[f].sort_values(ascending=False).head(30).index)
            _sr = {"subtype": sub, "factor": f, **_gene_factor_label(top_genes, hs, ec)}
            _sr["factor_label"] = _stable_factor_label({**_sr, "dominant_subtype": sub})
            sub_rows.append(_sr)
            for g in Ws[f].sort_values(ascending=False).head(TOP_LOAD).index:
                sub_load.append({"subtype": sub, "factor": f, "gene": g, "loading": float(Ws.loc[g, f])})

    out = {f"nmf_gene_pressure{tag}_loadings_gene": gene_load,
           f"nmf_gene_pressure{tag}_factor_subtype": act.sort_values("kw_q")}
    if sub_rows:
        out[f"nmf_gene_pressure{tag}_within_subtype"] = pd.DataFrame(sub_rows)
        out[f"nmf_gene_pressure{tag}_within_subtype_loadings"] = pd.DataFrame(sub_load)
    return out


# --------------------------------------------------------------------------- #
# 6b. SIGNED gene-pressure NMF (gain AND loss in one decomposition)
# --------------------------------------------------------------------------- #
def _nmf_signed_gene_pressure(
    tumor_mat: pd.DataFrame, edge: pd.DataFrame, hs: HallmarkSets, baseline: pd.Series, *,
    k: int = NMF_K,
) -> Dict[str, pd.DataFrame]:
    """Handle gain AND loss in ONE NMF (vs the acquired-only clip). The per-(gene,sample) dev
    pressure (`softmax_devhealthy_logrpm` vs the healthy baseline) is SIGNED — positive = more
    pressure than healthy, negative = relieved/less. NMF needs non-negative, so we split each gene
    into two channels `<gene>|gain = max(dev,0)` and `<gene>|loss = max(-dev,0)`, stack them
    (2·genes × samples), and decompose once. A factor can then load on a gene's gain channel,
    its loss channel, or both — so a single program captures co-ordinated gain AND relief, with
    direction preserved (unlike NMF on |dev|). Per-factor we report the top gain-genes and top
    loss-genes separately, the dominant Hallmark, and per-PAM50 activity (KW)."""
    from sklearn.decomposition import NMF
    from scipy.stats import kruskal

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _PAM50
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges
    from mirna_hallmark.config import PRESSURE_HEALTHY_ANCHOR_MODE
    from mirna_hallmark.stats import bh_fdr

    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    P = compute_gene_pressure(genes, edges=edges, mirna=tumor_mat,  # type: ignore[arg-type]
                              expr_mode=PRESSURE_HEALTHY_ANCHOR_MODE, resolve_arms=False,
                              healthy_baseline=baseline)
    if P.empty:
        return {}
    P = P.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    P = P.loc[P.abs().sum(axis=1) > 0]
    pos = P.clip(lower=0.0); pos.index = [f"{g}|gain" for g in P.index]
    neg = (-P).clip(lower=0.0); neg.index = [f"{g}|loss" for g in P.index]
    M = pd.concat([pos, neg]); M = M.loc[M.sum(axis=1) > 0]
    if min(M.shape) < k:
        return {}
    model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=1000)
    W = pd.DataFrame(model.fit_transform(M.values), index=M.index,
                     columns=[f"S{j + 1}" for j in range(k)])
    H = pd.DataFrame(model.components_, index=W.columns, columns=M.columns)
    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None

    clin = D.load_clinical_strata()
    part_pam = (clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
                if "PAM50_final" in clin.columns else {})
    col_pam = pd.Series({c: part_pam.get(_participant(c)) for c in M.columns})

    rows: List[dict] = []
    for f in H.index:
        load = W[f].sort_values(ascending=False)
        gain = [i.replace("|gain", "") for i in load.index if i.endswith("|gain")][:8]
        loss = [i.replace("|loss", "") for i in load.index if i.endswith("|loss")][:8]
        # dominant direction = which channel carries more of the factor's loading mass
        gain_mass = float(W.loc[[i for i in W.index if i.endswith("|gain")], f].sum())
        loss_mass = float(W.loc[[i for i in W.index if i.endswith("|loss")], f].sum())
        a = H.loc[f].astype(float)
        means, groups = [], []
        row: dict = {"factor": f, "direction": "gain" if gain_mass >= loss_mass else "loss",
                     "gain_mass_frac": gain_mass / (gain_mass + loss_mass + 1e-9)}
        for sub in _PAM50:
            cols = [c for c in col_pam.index[col_pam == sub] if c in a.index]
            vals = a.loc[cols].dropna()
            m = float(vals.mean()) if len(vals) else np.nan
            row[f"frac_{sub}"] = m
            means.append(m)
            if len(vals) >= 10:
                groups.append(vals.values)
        arr = np.array(means, float)
        if np.nansum(arr) > 0:
            for sub in _PAM50:
                row[f"frac_{sub}"] = row[f"frac_{sub}"] / np.nansum(arr) if np.isfinite(row[f"frac_{sub}"]) else np.nan
            row["dominant_subtype"] = _PAM50[int(np.nanargmax(arr))]
        try:
            row["kw_p"] = float(kruskal(*groups)[1]) if len(groups) >= 2 else np.nan
        except ValueError:
            row["kw_p"] = np.nan
        row["top_gain_genes"] = ";".join(gain)
        row["top_loss_genes"] = ";".join(loss)
        lab = _gene_factor_label(gain + loss, hs, ec)
        row["dominant_hallmark"] = lab.get("dominant_hallmark")
        row["top_driver_mirnas"] = lab.get("top_driver_mirnas")
        row["factor_label"] = _stable_factor_label(row)
        rows.append(row)
    act = pd.DataFrame(rows)
    act["kw_q"] = bh_fdr(act["kw_p"].fillna(1.0).values)
    load_rows = [{"factor": f, "gene_channel": g, "loading": float(W.loc[g, f])}
                 for f in W.columns for g in W[f].sort_values(ascending=False).head(TOP_LOAD).index]
    # Persist the FULL W (gene_channel x factor) and H (factor x sample) so this fit is the
    # single source of truth for the signed NMF. nmf_sample_signatures reuses these instead of
    # re-fitting its own copy — otherwise the two independent fits return permutation-arbitrary
    # factor indices (e.g. the LIFR/ABCE1/IGFBP5 backbone lands on S9 here but S5 there), which
    # makes a bare-index join across the loadings/factor_subtype and per-sample tables silently wrong.
    return {"nmf_gene_pressure_signed_factor_subtype": act.sort_values("kw_q"),
            "nmf_gene_pressure_signed_loadings": pd.DataFrame(load_rows),
            "nmf_gene_pressure_signed_W": W.rename_axis("gene_channel").reset_index(),
            "nmf_gene_pressure_signed_H": H.rename_axis("factor").reset_index()}


# --------------------------------------------------------------------------- #
# 5. Per-subtype NMF (two flavours: per-sample subset + static recomputed-c)
# --------------------------------------------------------------------------- #
def _nmf_per_subtype(
    tumor_mat: pd.DataFrame, edge: pd.DataFrame, hs: HallmarkSets, *, k: int = 8,
) -> Dict[str, pd.DataFrame]:
    """NMF inside each PAM50 subtype. (a) per-sample on the subtype's columns; (b) static
    on the contribution matrix whose cells are recomputed from only the subtype's samples."""
    from sklearn.decomposition import NMF

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import ATTR_MODE, _PAM50
    from mirna_hallmark.family_normal_reference import _participant
    from mirna_hallmark.pressure_build import compute_gene_pressure_contributions, load_mirtar_edges

    clin = D.load_clinical_strata()
    if "PAM50_final" not in clin.columns:
        return {}
    part_pam = clin.dropna(subset=["PAM50_final"]).set_index("participant")["PAM50_final"].to_dict()
    col_pam = {c: part_pam.get(_participant(c)) for c in tumor_mat.columns}

    genes = sorted(hs.universe)
    edges = load_mirtar_edges(genes, resolve_arms=True)
    ec = edge.groupby(["miRNA", "gene"])["c_tumor"].mean() if "c_tumor" in edge.columns else None

    persample_rows: List[dict] = []
    static_rows: List[dict] = []
    for sub in _PAM50:
        cols = [c for c in tumor_mat.columns if col_pam.get(c) == sub]
        if len(cols) < 30:
            continue
        # (a) per-sample NMF within this subtype
        expr = tumor_mat[cols].apply(pd.to_numeric, errors="coerce")
        expr = expr.loc[expr.median(axis=1) >= MIN_EXPR].dropna(how="any").clip(lower=0.0)
        if min(expr.shape) >= k:
            W = pd.DataFrame(
                NMF(n_components=k, init="nndsvda", random_state=0, max_iter=1000).fit_transform(expr.values),
                index=expr.index, columns=[f"{sub}_S{j + 1}" for j in range(k)])
            for f in W.columns:
                top_arms = list(W[f].sort_values(ascending=False).head(8).index)
                row = {"subtype": sub, "flavor": "per_sample", "factor": f, "top_mirnas": ";".join(top_arms)}
                if ec is not None:
                    sub_e = ec[ec.index.get_level_values("miRNA").isin(top_arms)]
                    row["top_target_module"] = ";".join(
                        sub_e.groupby(level="gene").sum().sort_values(ascending=False).head(8).index)
                persample_rows.append(row)
        # (b) static NMF on c recomputed from this subtype's samples (no sample axis)
        contrib = compute_gene_pressure_contributions(
            genes, edges=edges, mirna=tumor_mat[cols], expr_mode=ATTR_MODE, resolve_arms=False)  # type: ignore[arg-type]
        if contrib.empty:
            continue
        Cm = (contrib.pivot_table(index="miRNA", columns="gene", values="mean_abs_contribution", aggfunc="mean")
              .fillna(0.0).clip(lower=0.0))
        if min(Cm.shape) >= k:
            model = NMF(n_components=k, init="nndsvda", random_state=0, max_iter=500)
            Wc = pd.DataFrame(model.fit_transform(Cm.values), index=Cm.index,
                              columns=[f"{sub}_C{j + 1}" for j in range(k)])
            Hc = pd.DataFrame(model.components_, index=Wc.columns, columns=Cm.columns)
            for f in Wc.columns:
                static_rows.append({
                    "subtype": sub, "flavor": "static", "factor": f,
                    "top_mirnas": ";".join(Wc[f].sort_values(ascending=False).head(8).index),
                    "top_genes": ";".join(Hc.loc[f].sort_values(ascending=False).head(8).index)})

    out: Dict[str, pd.DataFrame] = {}
    if persample_rows:
        out["nmf_subtype_persample"] = pd.DataFrame(persample_rows)
    if static_rows:
        out["nmf_subtype_static"] = pd.DataFrame(static_rows)
    return out


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def run(*, out_dir: Path = OUT_DIR, state_dir: Path = STATE_DIR,
        cov_source: str = "estimate_epi", do_corepression: bool = True,
        do_subtype_nmf: bool = True, do_gene_pressure_nmf: bool = True,
        do_acquired_nmf: bool = True) -> Dict[str, pd.DataFrame]:
    from mirna_hallmark.analyses.cross_state.cross_state_landscape import _state_matrices

    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    mir_class, edge = _load_state_class(state_dir)
    print(f"[comovement] loaded {len(mir_class)} arms, {len(edge)} edges from {state_dir.name}")

    states = _state_matrices()
    modules, mod_summary = _comovement_modules(states["tumor"], mir_class)
    traj = _trajectory_clusters(mir_class)
    if not modules.empty and not traj.empty:
        modules = modules.merge(traj.rename("trajectory_cluster").reset_index().rename(columns={"index": "miRNA"}),
                                on="miRNA", how="left")

    nmf = _nmf_programs(edge)
    nmf_share = _nmf_share(edge)
    nmf_wt = _nmf_within_tumor(states["tumor"], edge, mir_class)
    nmf_sub = _nmf_per_subtype(states["tumor"], edge, hs) if do_subtype_nmf else {}
    nmf_gp = _nmf_gene_pressure(states["tumor"], edge, hs) if do_gene_pressure_nmf else {}
    # ACQUIRED (dev-anchored) gene-pressure NMF: programs of genes co-GAINING pressure vs the
    # healthy baseline (complements the abundance view above). Needs the QN'd healthy baseline.
    nmf_gp_acq: Dict[str, pd.DataFrame] = {}
    nmf_gp_signed: Dict[str, pd.DataFrame] = {}
    if do_gene_pressure_nmf and do_acquired_nmf:
        from mirna_hallmark.config import PRESSURE_HEALTHY_ANCHOR_MODE
        from mirna_hallmark.healthy_anchor import load_baseline
        try:
            base = load_baseline()
            nmf_gp_acq = _nmf_gene_pressure(states["tumor"], edge, hs,
                                            expr_mode=PRESSURE_HEALTHY_ANCHOR_MODE,
                                            baseline=base, tag="_acquired")
            # SIGNED: gain AND loss in one decomposition (±-channel stack)
            nmf_gp_signed = _nmf_signed_gene_pressure(states["tumor"], edge, hs, base)
        except Exception as exc:
            print(f"[comovement] WARN acquired/signed gene-pressure NMF skipped ({exc})")

    corep = pd.DataFrame()
    if do_corepression:
        corep = _target_corepression(edge, hs, cov_source=cov_source)
        try:
            gene_corepression(hs, cov_source=cov_source, out_dir=out_dir)
        except Exception as exc:
            print(f"[comovement] WARN gene_corepression skipped ({exc})")

    if not modules.empty:
        modules.to_csv(out_dir / "comovement_modules.tsv", sep="\t", index=False)
        mod_summary.to_csv(out_dir / "comovement_module_summary.tsv", sep="\t", index=False)
    if not corep.empty:
        corep.to_csv(out_dir / "target_corepression.tsv", sep="\t", index=False)
    for name, df in {**nmf, **nmf_share, **nmf_wt, **nmf_sub, **nmf_gp,
                     **nmf_gp_acq, **nmf_gp_signed}.items():
        if not df.empty:
            df.to_csv(out_dir / f"{name}.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.analyses.misc.mirna_comovement",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "state_class_dir": str(state_dir),
        "n_arms": int(len(mir_class)), "n_edges": int(len(edge)),
        "n_comove_modules": int(modules["comove_module"].nunique()) if not modules.empty else 0,
        "nmf_k": NMF_K, "nmf_within_tumor": bool(nmf_wt), "nmf_per_subtype": bool(nmf_sub),
        "nmf_gene_pressure": bool(nmf_gp), "cov_source": cov_source if do_corepression else None,
        "params": {"min_expr": MIN_EXPR, "min_targets": MIN_TARGETS, "n_modules": N_MODULES,
                   "top_load": TOP_LOAD},
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    _print_summary(modules, mod_summary, corep, nmf, nmf_wt, nmf_sub, nmf_gp)
    if nmf_gp_acq:
        acq = nmf_gp_acq.get("nmf_gene_pressure_acquired_factor_subtype")
        if acq is not None and not acq.empty:
            print("\n[comovement] ACQUIRED (dev) gene-pressure NMF - most PAM50-divergent programs:")
            for _, r in acq.head(5).iterrows():
                print(f"  {r.factor}: dom={r.get('dominant_subtype')} kw_q={r.get('kw_q', np.nan):.1e} "
                      f"[{r.get('dominant_hallmark')}] genes={r.get('top_genes', '')}")
    if nmf_gp_signed:
        sg = nmf_gp_signed.get("nmf_gene_pressure_signed_factor_subtype")
        if sg is not None and not sg.empty:
            print("\n[comovement] SIGNED (gain+loss) gene-pressure NMF programs:")
            for _, r in sg.head(6).iterrows():
                print(f"  {r.factor}: dir={r.get('direction')} (gain_frac={r.get('gain_mass_frac', np.nan):.2f}) "
                      f"dom={r.get('dominant_subtype')} kw_q={r.get('kw_q', np.nan):.1e} "
                      f"gain=[{r.get('top_gain_genes', '')}] loss=[{r.get('top_loss_genes', '')}]")
    return {"modules": modules, "module_summary": mod_summary, "corepression": corep,
            **nmf, **nmf_share, **nmf_wt, **nmf_sub, **nmf_gp, **nmf_gp_acq, **nmf_gp_signed}


def _print_summary(modules, mod_summary, corep, nmf, nmf_wt=None, nmf_sub=None, nmf_gp=None) -> None:
    if not mod_summary.empty:
        print("\n[comovement] top co-expression modules (n_arms, mean dHT, dominant class):")
        for _, r in mod_summary.head(6).iterrows():
            print(f"  M{r.comove_module}: n={r.n_arms} dHT={r.mean_dHT:+.2f} "
                  f"{r.dominant_primary}/{r.dominant_secondary} :: {r.members}")
    if not corep.empty:
        print("\n[comovement] strongest acquired target-module repression (delta tumor-gtex):")
        sub = corep.dropna(subset=["delta_tumor_gtex"]).sort_values("delta_tumor_gtex").head(8) \
            if "delta_tumor_gtex" in corep.columns else corep.head(8)
        for _, r in sub.iterrows():
            g = r.get("rho_composite_gtex", np.nan)
            t = r.get("rho_composite_tumor", np.nan)
            print(f"  {r.miRNA}: composite rho gtex {g:+.2f} -> tumor {t:+.2f} (n_targets={int(r.n_targets)})")
    act = nmf.get("nmf_factor_state_activity")
    if act is not None and not act.empty and "shift_tumor_gtex" in act.columns:
        print("\n[comovement] NMF programs the tumor amplifies most (frac shift tumor-gtex):")
        ml = nmf.get("nmf_factor_loadings_mirna")
        for _, r in act.head(4).iterrows():
            tops = ""
            if ml is not None and not ml.empty:
                tops = ";".join(ml.loc[ml.factor.eq(r.factor), "miRNA"].head(5))
            print(f"  {r.factor}: shift {r.shift_tumor_gtex:+.3f} :: {tops}")
    if nmf_wt:
        wt = nmf_wt.get("nmf_within_tumor_factor_subtype")
        if wt is not None and not wt.empty:
            print("\n[comovement] within-tumor (per-sample) NMF programs - most PAM50-divergent:")
            for _, r in wt.head(6).iterrows():
                q = r.get("kw_q", np.nan)
            print(f"  {r.factor}: dom={r.get('dominant_subtype')} kw_q={q:.1e} :: "
                  f"{r.get('top_mirnas', '')}  -> {r.get('top_target_module', '')}")
    if nmf_sub:
        for key, label in (("nmf_subtype_static", "static"), ("nmf_subtype_persample", "per-sample")):
            df = nmf_sub.get(key)
            if df is not None and not df.empty:
                print(f"\n[comovement] per-subtype NMF ({label}) - first program per subtype:")
                for sub, grp in df.groupby("subtype"):
                    r = grp.iloc[0]
                    tail = r.get("top_genes", r.get("top_target_module", ""))
                    print(f"  {sub} {r.factor}: {r.get('top_mirnas', '')} -> {tail}")
    if nmf_gp:
        gp = nmf_gp.get("nmf_gene_pressure_factor_subtype")
        if gp is not None and not gp.empty:
            print("\n[comovement] gene-pressure (per-sample) NMF - most PAM50-divergent programs:")
            for _, r in gp.head(6).iterrows():
                q = r.get("kw_q", np.nan)
                print(f"  {r.factor}: dom={r.get('dominant_subtype')} kw_q={q:.1e} "
                      f"[{r.get('dominant_hallmark')}] genes={r.get('top_genes', '')} "
                      f":: drivers={r.get('top_driver_mirnas', '')}")
        ws = nmf_gp.get("nmf_gene_pressure_within_subtype")
        if ws is not None and not ws.empty:
            print("\n[comovement] within-subtype gene-pressure programs (first per subtype):")
            for sub, grp in ws.groupby("subtype"):
                r = grp.iloc[0]
                print(f"  {sub} {r.factor}: [{r.get('dominant_hallmark')}] {r.get('top_genes', '')} "
                      f":: drivers={r.get('top_driver_mirnas', '')}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state-dir", type=Path, default=STATE_DIR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--cov-source", choices=["estimate", "estimate_epi", "metagene"],
                    default="estimate_epi",
                    help="composition covariates for target-set corepression (default: estimate_epi)")
    ap.add_argument("--no-corepression", action="store_true",
                    help="skip the (heavier) target-set composite partial-rho step")
    ap.add_argument("--no-subtype-nmf", action="store_true",
                    help="skip the per-subtype NMF (per-sample + static-c flavours)")
    ap.add_argument("--no-gene-pressure-nmf", action="store_true",
                    help="skip the per-sample gene-pressure NMF (gene x sample)")
    ap.add_argument("--no-acquired-nmf", action="store_true",
                    help="skip the acquired (dev-anchored) gene-pressure NMF")
    args = ap.parse_args()
    run(out_dir=args.out_dir, state_dir=args.state_dir, cov_source=args.cov_source,
        do_corepression=not args.no_corepression, do_subtype_nmf=not args.no_subtype_nmf,
        do_gene_pressure_nmf=not args.no_gene_pressure_nmf,
        do_acquired_nmf=not args.no_acquired_nmf)


if __name__ == "__main__":
    main()
