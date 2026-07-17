"""Seed-family pooling (Design §Decision F).

Same-seed arms are near-collinear ⇒ plain lasso keeps one member arbitrarily, unstably across folds.
The identified estimand is **family→gene**; the per-arm split is a nomination, not a coefficient.

MVP variant = **collapse-to-family predictor** (the design's cheap, valid variant; the repo's family rung):
pool member arms into one predictor (abundance summed in linear RPM, re-logged), fit the gene-focused
NN adaptive-lasso on family predictors → family→gene weights. Uses the real TargetScan seed family
(`miRNA_family` in the edge table), which correctly splits miR-200a/141 from miR-200b/c/429.

(The full hierarchical partial-pooling — family weight + shrunk per-member δ with posterior width — is
Phase 3b; collapse is its infinite-shrinkage limit.)
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import data_loaders as D


_FAM_MAP: dict = {}


def _arm_family_map() -> pd.Series:
    """arm → miRNA_family Series, CACHED (deterministic from the edge table). Was rebuilt on EVERY `family_of` call —
    a full `load_hallmark_edges` (+ defensive copy) + `drop_duplicates` — which the exclusion batch pays ~1000×/run
    (profiled = 45% of runtime). Cache it once."""
    if "m" not in _FAM_MAP:
        e = D.load_hallmark_edges()
        _FAM_MAP["m"] = e[["miRNA", "miRNA_family"]].drop_duplicates("miRNA").set_index("miRNA")["miRNA_family"]
    return _FAM_MAP["m"]


def family_of(arms, edges: Optional[pd.DataFrame] = None) -> pd.Series:
    """arm → TargetScan seed family (from the edge table's `miRNA_family`; 100% coverage). Cached mapping for the
    default edge table (`_arm_family_map`); an explicit `edges` overrides it (uncached, for custom tables)."""
    if edges is not None:
        m = edges[["miRNA", "miRNA_family"]].drop_duplicates("miRNA").set_index("miRNA")["miRNA_family"]
    else:
        m = _arm_family_map()
    fam = m.reindex(list(arms))
    return fam.fillna(pd.Series(list(arms), index=list(arms)))  # singletons keep their own name


_COMP_CACHE: dict = {}
_SEED_COMP_TSV = "mirna_hallmark/output/learned/seed_composition.tsv.gz"


def _seed_comp() -> pd.DataFrame:
    """Cached per-(participant, arm, target_family) seed-composition `c_{s,m,f}` (`within_family.seed_composition`, Phase 1)
    — the fraction of arm m's reads in participant s that carry family f's seed (canonical + 5'-isomiR JUMP legs). Persisted
    once (parses ~1078 isoform files). Columns: participant, arm, target_family, frac, n_reads, is_canonical."""
    if "c" not in _COMP_CACHE:
        import os
        if os.path.exists(_SEED_COMP_TSV):
            _COMP_CACHE["c"] = pd.read_csv(_SEED_COMP_TSV, sep="\t")
        else:
            from mirna_hallmark.learned import within_family as WF
            c = WF.seed_composition()
            os.makedirs(os.path.dirname(_SEED_COMP_TSV), exist_ok=True)
            c.to_csv(_SEED_COMP_TSV, sep="\t", index=False)
            _COMP_CACHE["c"] = c
    return _COMP_CACHE["c"]


def collapse_by_family(
    X: pd.DataFrame,
    w_prior: pd.Series,
    fam: pd.Series,
    *,
    isomir: bool = False,
    comp: pd.DataFrame | None = None,
) -> Tuple[pd.DataFrame, pd.Series, dict]:
    """Pool member arms into one family predictor. Returns (X_fam, w_fam, members).

    X_fam : participant × family   pooled abundance = log2(1 + Σ_members (2^x − 1))  (true RPM pool)
    w_fam : family                 prior weight = max member weight (family ≈ strongest member)
    members : family → [arms]

    **`isomir=True` (Phase 2, `ISOMIR_AWARE_MODELING.md`, DEFAULT OFF — the 22 callers are unaffected):** the
    seed-heterogeneity edge re-definition — `X_fam,f[s] = log2(1 + Σ_m c_{s,m,f}·RPM_{s,m})`, so a family GAINS an arm's
    5'-isomiR reads whose shifted seed jumps INTO f (e.g. miR-29b +1 → miR-767) and LOSES its own reads that jump away
    (orphan/other). Arms below the isomiR read-floor (no `c` row) keep the **canonical assumption** (full RPM to their own
    family), so families with no isomiR coverage reduce EXACTLY to the plain pool. Gated (Phase 4): keep OFF unless a
    family's held-out coupling improves — the estimand is "reads carrying family f's seed," a cleaner instrument (§8-enforcing)."""
    fam = fam.reindex(X.columns)
    members: dict = {}
    w = {}
    lin = np.power(2.0, X.to_numpy(dtype=float)) - 1.0          # log2(RPM+1) → RPM
    lin = pd.DataFrame(lin, index=X.index, columns=X.columns)
    for f, arms in fam.groupby(fam).groups.items():
        members[f] = list(arms)
        w[f] = float(w_prior.reindex(list(arms)).max())
    if not isomir:
        cols = {f: np.log2(1.0 + lin[members[f]].sum(axis=1)) for f in members}
        return pd.DataFrame(cols, index=X.index), pd.Series(w, name="w_prior"), members

    # ---- Phase 2: isomiR-corrected pool (opt-in, gated) ----
    c = comp if comp is not None else _seed_comp()
    fams = set(members)
    c = c[(c["arm"].isin(X.columns)) & (c["target_family"].isin(fams)) & (c["frac"] > 0.0)]
    covered = set(c["arm"].unique())                                   # arms with isomiR coverage (enough reads)
    # part B — covered arms' c-weighted contribution to each family (canonical share + any jump-ins)
    ll = lin.copy(); ll.index.name = "participant"
    ll = ll.reset_index().melt(id_vars="participant", var_name="arm", value_name="rpm")
    mg = c.merge(ll, on=["participant", "arm"], how="inner")
    mg["wr"] = mg["frac"] * mg["rpm"]
    partB = (mg.groupby(["participant", "target_family"])["wr"].sum()
             .unstack(fill_value=0.0).reindex(index=X.index, columns=sorted(fams), fill_value=0.0))
    # part A — UNCOVERED canonical members keep their full RPM (no isomiR data ⇒ canonical assumption)
    partA = pd.DataFrame({f: lin[[m for m in members[f] if m not in covered]].sum(axis=1) if
                          any(m not in covered for m in members[f]) else pd.Series(0.0, index=X.index)
                          for f in sorted(fams)}, index=X.index)
    X_fam = np.log2(1.0 + partA.add(partB, fill_value=0.0))[sorted(fams)]
    return X_fam, pd.Series(w, name="w_prior"), members


def resolution_report(gene: str, *, alpha: float = 0.005) -> pd.DataFrame:
    """Phase-3b: for each model-selected multi-member family, which member (if any) is the DRIVER?

    The identified estimand is family→gene (Design §F). A member is the driver only if its OWN variation,
    independent of its family-mates, still tracks the target — so nomination is **target-driven**: residualise
    the member arm AND the target on (the other family members + C), and test whether the anti-correlation
    SURVIVES (conditioned partial-Spearman < 0). This unifies identifiability and attribution — a collinear
    member has ~no residual variance after conditioning ⇒ automatically unidentified (no 0.7 threshold needed);
    a member that merely rides the family aggregate drops to ~0; the driver is the one whose conditioned
    anti-corr stays negative. within_corr (abundance collinearity) is kept only as the identifiability
    diagnostic. (Replaces the old abundance-divergence heuristic, which never looked at the target.)
    """
    from scipy.stats import spearmanr
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import states as ST                       # lazy: avoid families↔states import cycle
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = family_of(X.columns)
    Xf, wf, members = collapse_by_family(X, w, fam)
    Mf = ST.canonical_M(gene, "01", arm_level=False)                      # canonical family weight (gates which to examine)
    Cm = C.to_numpy(float); yv = Y.to_numpy(float)

    def _resid(v, Z):
        return v - LinearRegression().fit(Z, v).predict(Z)

    rows = []
    for f, wt in Mf[Mf > 0].sort_values(ascending=False).items():
        mem = [m for m in members.get(f, [f]) if m in X.columns]
        if len(mem) < 2:
            continue
        cc = X[mem].corr(method="spearman").abs()
        maxc = float(cc.where(~np.eye(len(mem), dtype=bool)).max().max())  # abundance collinearity (diagnostic)
        part = {}
        for arm in mem:                                                   # conditioned partial anti-corr per member
            others = [m for m in mem if m != arm]
            Z = np.column_stack([Cm, X[others].to_numpy(float)]) if others else Cm
            xr = _resid(X[arm].to_numpy(float), Z)
            if float(np.std(xr)) < 1e-6:                                  # collinear → no independent variation
                part[arm] = np.nan; continue
            part[arm] = float(spearmanr(xr, _resid(yv, Z)).correlation)
        ps = pd.Series(part)
        surv = ps[ps < -0.1].sort_values()                              # survivors (real conditioned repression), most-negative first
        n_surv = len(surv)
        # the nomination reflects HOW the surviving anti-corr distributes, not a forced single winner:
        klass = "family-only" if n_surv == 0 else ("single-driver" if n_surv == 1 else "co-drivers")
        rows.append({"gene": gene, "family": f, "family_weight": round(float(wt), 3),
                     "n_members": len(mem), "within_corr": round(maxc, 2),
                     "resolution": klass, "n_drivers": n_surv,
                     "driver": surv.index[0].replace("hsa-", "") if n_surv else "(family-only)",   # strongest survivor (compat)
                     "drivers": ", ".join(m.replace("hsa-", "") for m in surv.index) if n_surv else "(family-only)",
                     "driver_rho": round(float(surv.min()), 3) if n_surv else np.nan,
                     "member_partials": ", ".join(f"{m.replace('hsa-','')}={v:+.2f}" if v == v else f"{m.replace('hsa-','')}=collinear"
                                                   for m, v in ps.items()),
                     "members": ", ".join(mem)})
    return pd.DataFrame(rows)


def cluster_resolution(gene: str, *, rho_thr: float = 0.7) -> pd.DataFrame:
    """Kind-B (genomic-cluster) collinearity resolver — see COLLINEARITY_AND_IDENTIFIABILITY.md §2–3.

    Among a gene's regulator FAMILIES (already seed-collapsed, so this is CROSS-seed), find abundance-collinear
    clusters (|ρ_abund| ≥ rho_thr = co-transcribed riders) and resolve who actually regulates g by two levers
    abundance cannot separate but which ARE arm-specific here (different seeds ⇒ different targetomes):
      (1) SEQUENCE — biochemical potential on g (scanMiR/TargetScan `strength`); a co-transcribed mate with NO
          site is a rider, not a regulator (the cheap sequence-first prune);
      (2) conditioned COUPLING — partial-Spearman of the family's abundance vs g, residualising on the OTHER
          cluster-mates + C (generalises `resolution_report`/`instrument.cluster_attribution` to the
          observational cross-seed case).
    role: owner (site ∧ couples) · rider (no site ∧ no residual coupling) · designed-not-coupling (site, no
    residual coupling — silenced or a mate owns it) · couples-no-site (reference-blind: SNV/editing/APA/uncurated).
    """
    import itertools
    from scipy.stats import spearmanr
    from sklearn.linear_model import LinearRegression
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import structural_identity as SI
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = family_of(X.columns)
    Xf, _wf, members = collapse_by_family(X, w, fam)
    fams = list(Xf.columns)
    if len(fams) < 2:
        print(f"{gene}: <2 families"); return pd.DataFrame()
    corr = Xf.corr(method="spearman").abs()
    adj = {f: set() for f in fams}
    for a, b in itertools.combinations(fams, 2):
        if corr.loc[a, b] >= rho_thr:
            adj[a].add(b); adj[b].add(a)
    seen: set = set(); clusters = []
    for f in fams:                                                        # connected components of the collinearity graph
        if f in seen or not adj[f]:
            continue
        comp: set = set(); stack = [f]
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u); stack += [v for v in adj[u] if v not in comp]
        seen |= comp
        if len(comp) > 1:
            clusters.append(sorted(comp))
    if not clusters:
        print(f"{gene}: no cross-seed collinear clusters (|ρ_abund|≥{rho_thr})"); return pd.DataFrame()

    strength, _, _ = SI._potential_matrix("fused")

    def _pot(f):
        mem = [m for m in members.get(f, [f]) if m in strength.index]
        return float(strength.loc[mem, gene].max()) if (mem and gene in strength.columns) else 0.0

    Cm = C.to_numpy(float); yv = Y.to_numpy(float)

    def _resid(v, Z):
        return v - LinearRegression().fit(Z, v).predict(Z)

    try:                                                                  # CN-overlay (B2): locus CN as instrument
        from mirna_hallmark.learned import instrument as INST
        cn = INST.arm_cn()
    except Exception:
        cn = None

    def _cluster_cn(comp):
        """Per-family locus-CN reduced-form for a cluster (B2). Different loci ⇒ CN separable even where
        expression is collinear: partial-Spearman(CN_f, g | C + other-loci CN) = f's EXOGENOUS effect on g;
        first-stage = CN_f→abundance. Returns {family: (cn_reduced, cn_fs)}; empty if <2 loci CN separable."""
        if cn is None:
            return {}
        fam_arm = {}
        for f in comp:                                                   # a member arm carrying locus CN
            for a in members.get(f, [f]):
                if a in cn.columns and cn[a].notna().any():
                    fam_arm[f] = a; break
        if len(fam_arm) < 2:
            return {}
        parts = [p for p in Xf.index if p in cn.index]
        Zd = pd.DataFrame({f: pd.to_numeric(cn.loc[parts, a], errors="coerce") for f, a in fam_arm.items()},
                          index=parts).dropna()
        if Zd.shape[0] < 40:
            return {}
        p2 = Zd.index; Cm2 = C.loc[p2].to_numpy(float); yv2 = Y.loc[p2].to_numpy(float)
        out = {}
        for f in Zd.columns:
            oth = [o for o in Zd.columns if o != f]
            Zbase = np.column_stack([Cm2, Zd[oth].to_numpy(float)]) if oth else Cm2
            cnf = _resid(Zd[f].to_numpy(float), Zbase)
            cnf_c = _resid(Zd[f].to_numpy(float), Cm2)                   # CN_f | C only
            sd_c = float(np.std(cnf_c))
            if float(np.std(cnf)) < 1e-6 or (sd_c > 0 and float(np.std(cnf)) / sd_c < 0.4):
                continue                                                 # CN_f mostly shared w/ a mate ⇒ SAME locus (B1) → CN cannot separate
            red = float(spearmanr(cnf, _resid(yv2, Zbase)).correlation)  # exogenous CN→g | other loci + C
            fs = float(spearmanr(cnf_c, _resid(Xf.loc[p2, f].to_numpy(float), Cm2)).correlation)
            out[f] = (round(red, 3), round(fs, 2))
        return out

    rows = []
    for ci, comp in enumerate(clusters):
        cnmap = _cluster_cn(comp)
        for f in comp:
            others = [o for o in comp if o != f]
            Z = np.column_stack([Cm, Xf[others].to_numpy(float)]) if others else Cm
            xr = _resid(Xf[f].to_numpy(float), Z)
            coup = np.nan if float(np.std(xr)) < 1e-6 else float(spearmanr(xr, _resid(yv, Z)).correlation)
            p = _pot(f); couples = bool(coup == coup and coup < -0.1)
            role = ("owner" if (p > 0 and couples) else "rider" if (p <= 0 and not couples)
                    else "designed-not-coupling" if (p > 0 and not couples) else "couples-no-site")
            cnr, cnfs = cnmap.get(f, (np.nan, np.nan))
            # per-family CAUSAL verdict — CN can only speak where its first stage is adequate (fs≥0.15, CN↑⇒arm↑)
            cn_causal = ("CN-repressor" if (cnfs == cnfs and cnfs >= 0.15 and cnr == cnr and cnr <= -0.1)
                         else "CN-null" if (cnfs == cnfs and cnfs >= 0.15) else "CN-blind")
            rows.append({"cluster": ci, "family": f, "n_mates": len(comp),
                         "max_rho_abund": round(float(corr.loc[f, others].max()), 2),
                         "potential": round(p, 2), "has_site": bool(p > 0),
                         "cond_coupling": round(coup, 3) if coup == coup else np.nan, "couples": couples,
                         "cn_reduced": cnr, "cn_fs": cnfs, "cn_causal": cn_causal, "role": role})
    df = pd.DataFrame(rows).sort_values(["cluster", "potential"], ascending=[True, False])
    with pd.option_context("display.width", 165):
        print(f"\n=== {gene} — Kind-B cluster resolution (cross-seed abundance-collinear ≥{rho_thr}) ===")
        print(df.to_string(index=False))
    for ci, comp in enumerate(clusters):
        sub = df[df["cluster"] == ci]
        own = list(sub.loc[sub["role"] == "owner", "family"]); rid = list(sub.loc[sub["role"] == "rider", "family"])
        # CN is a per-family CAUSAL verdict where its instrument is strong — NOT a forced single owner. Report
        # which mates CN can adjudicate and whether it causally supports repression (agreement w/ coupling is
        # only meaningful for CN-adjudicable families).
        rep = list(sub.loc[sub["cn_causal"] == "CN-repressor", "family"])
        adj = list(sub.loc[sub["cn_causal"].isin(["CN-repressor", "CN-null"]), "family"])
        cn_tag = (f"  CN-causal repressor: {rep}" if rep
                  else f"  CN-adjudicable {adj} but no causal repressor (mild/null)" if adj
                  else "  CN-blind (weak locus instruments)")
        print(f"  cluster {ci} {list(comp)}: coupling-owner={own or '—'}  rider={rid or '—'}{cn_tag}")
    return df
