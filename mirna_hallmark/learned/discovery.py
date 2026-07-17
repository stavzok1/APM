"""Discovery lane (Design §Decision D/E soft-discovery) — where the prior actually pays off.

Bar 3 showed the prior's value is INCLUSION (which arms), not magnitude. Curation's inclusion is
incomplete: the heavily-studied edges are often buffered in bulk breast while less-studied breast-active
edges carry the cross-sample variance. So use the *sequence* inclusion prior (TargetScan-predicted sites,
`assemble_gene(orphans=True)`) to NOMINATE arms outside the curated HE set, and keep only those that add
**repressive coupling beyond the curated aggregate** with **biochemical support** (scanMiR K_D).

Criterion for a discovery (all three):
  1. orphan (TargetScan site, NOT in curated HE for this gene);
  2. **partial** repressive coupling — corr(arm, Y | C, HE-aggregate) < threshold  (adds signal curation missed);
  3. scanMiR predicted repression < 0 (a biochemical duplex, not a spurious correlation).
Weak sub-HE ledger evidence (if any) is reported as corroboration, not required.

CLI: `python -m mirna_hallmark.learned.discovery`  (scan a panel, rank cohort-wide discoveries)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.linear_model import LinearRegression

from mirna_hallmark import config as CFG
from mirna_hallmark.learned import attribution_eb as AE
from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import kd as KD
from mirna_hallmark.learned import regression as LR
from mirna_hallmark.learned.evidence import ledger as LG

PANEL = ["PTEN", "ESR1", "ZEB1", "CDKN1A", "GATA3", "BCL2", "MYC", "CDKN1B", "RB1", "SMAD4"]
C_CARD = CFG.OUTPUT_ROOT / "learned" / "arm_attribution_card.tsv"   # canonical Gibbs β per (gene, arm)


def _resid(v, Cmat):
    return v - LinearRegression().fit(Cmat, v).predict(Cmat)


def discover_gene(gene: str, *, min_partial: float = -0.12, alpha: float = 0.01,
                  min_scanmir: float = 0.0, permute: int | None = None,
                  deconv_check: bool = True) -> pd.DataFrame:
    """Orphan edges that add repressive coupling to `gene` beyond the curated HE aggregate + scanMiR support.
    `permute` shuffles the target (breaks the miRNA→target relation) → the FDR null. `min_scanmir` requires
    scanMiR predicted repression below it (stringency)."""
    try:
        Y, Xo, C, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True)
        _, Xhe, _, whe = LD.assemble_gene(gene, w_prior_source="ledger")
    except Exception:
        return pd.DataFrame()
    if permute is not None:
        Y = pd.Series(np.random.default_rng(permute).permutation(Y.to_numpy()), index=Y.index)
    he = set(Xhe.columns)
    orphans = [a for a in Xo.columns if a not in he]
    if not orphans:
        return pd.DataFrame()
    Cm = C.to_numpy(float)
    # curated HE aggregate (fit on HE only), added to C so orphan coupling is PARTIAL on what curation captures
    M_he = LR.fit_gene(Y, Xhe, C, whe, alpha=alpha)
    he_agg = Xhe.to_numpy(float) @ M_he.reindex(Xhe.columns).fillna(0).to_numpy()
    Cext = np.column_stack([Cm, he_agg])
    yr = _resid(Y.to_numpy(float), Cext)
    aff = KD.genome_affinity()      # genome-wide (unbiased) biochemical support — was HE-restricted KD.affinity()
    affg = aff.loc[aff["gene"] == gene].set_index("arm")["repression"]
    lw = LG.edge_weights(); lwg = lw.loc[lw["gene"] == gene].set_index("arm")["ledger_weight"]
    rows = []
    for a in orphans:
        xr = _resid(Xo[a].to_numpy(float), Cext)
        rho = spearmanr(xr, yr).correlation
        rows.append({"gene": gene, "arm": a, "partial_coupling": round(float(rho), 3),
                     "scanmir_rep": round(float(affg.get(a, np.nan)), 2) if a in affg.index else np.nan,
                     "sub_he_evidence": round(float(lwg.get(a, 0.0)), 2)})
    df = pd.DataFrame(rows)
    keep = df[(df["partial_coupling"] < min_partial) & (df["scanmir_rep"] < min_scanmir)].copy()
    if keep.empty or not deconv_check:
        return keep
    # composition control: recompute partial coupling with CIBERSORTx non-malignant fractions in C.
    # A real edge survives; a composition artefact (e.g. luminal-fraction-driven) collapses.
    try:
        Yd, Xod, Cd, _ = LD.assemble_gene(gene, w_prior_source="ledger", orphans=True, deconv=True)
        if permute is not None:
            Yd = pd.Series(np.random.default_rng(permute).permutation(Yd.to_numpy()), index=Yd.index)
        _, Xhed, _, whed = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
        Cmd = Cd.to_numpy(float)
        aggd = Xhed.to_numpy(float) @ LR.fit_gene(Yd, Xhed, Cd, whed, alpha=alpha).reindex(Xhed.columns).fillna(0).to_numpy()
        Cextd = np.column_stack([Cmd, aggd])
        yrd = _resid(Yd.to_numpy(float), Cextd)
        keep["partial_deconv"] = [round(float(spearmanr(_resid(Xod[a].to_numpy(float), Cextd), yrd).correlation), 3)
                                  if a in Xod.columns else np.nan for a in keep["arm"]]
    except Exception:
        keep["partial_deconv"] = np.nan
    keep["robust"] = keep["partial_deconv"] < min_partial * 0.6           # retains >60% of signal under deconv
    return keep.sort_values("partial_coupling")


_BATCH: dict = {}


def _batch_np(parts) -> np.ndarray:
    """plate_both batch dummies (post-fit conditioning only) for the given participants, cached full-cohort."""
    if "bd" not in _BATCH:
        from mirna_hallmark import tcga_batch as TB
        _BATCH["bd"] = TB.tcga_batch_dummies(list(LD._load()["X"].columns), kind="plate_both")
    bd = _BATCH["bd"].reindex(list(parts)).fillna(0.0)
    return bd.loc[:, bd.sum(axis=0) > 0].to_numpy(float)


def _he_arms_map() -> dict:
    """gene → set of curated POOLED-HE regulator arms present in the arm matrix (one build, avoids a 2nd
    assemble). POOLED-HE (migration) = the 'known' baseline, so orphans are edges beyond it."""
    he = LG.pooled_he_edges()
    xidx = set(LD._load()["X"].index)
    return {str(g): set(sub["miRNA"]) & xidx for g, sub in he.groupby("gene")}


def _resid_mat(V: np.ndarray, Cmat: np.ndarray) -> np.ndarray:
    """Residualise a BLOCK of columns on Cmat in one lstsq (vs one LinearRegression fit per arm)."""
    A = np.column_stack([np.ones(len(Cmat)), Cmat])
    beta, *_ = np.linalg.lstsq(A, V, rcond=None)
    return V - A @ beta


def _rho_block(Vr: np.ndarray, yr: np.ndarray) -> np.ndarray:
    """Pearson on pre-residualised RANKS ≡ partial Spearman, vectorised over columns of Vr."""
    a = Vr - Vr.mean(0)
    b = yr - yr.mean()
    d = np.sqrt((a * a).sum(0) * (b * b).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(d > 0, (a * b[:, None]).sum(0) / d, np.nan)


# ---- shared read-only state for the worker pool (module-level so forked workers inherit it) ----
_SCAN: dict = {}


def _scan_one_gene(g: str):
    """One gene's candidates + site-free null draws, for the worker pool. Returns (rows, nulls).

    THE SITE-FREE NULL for THIS gene, through the SAME Cext / yr / sample support as the candidates — the
    IDENTICAL estimator (`spearmanr(_resid(x, Cext), yr)`, verified bit-equal to the vectorised block form).
    RNG is seeded by (null_seed, gene name) so the draw is independent of pool ordering.
    """
    hemap, aff, lw, Xall, DETECTED, ab_all = (_SCAN["hemap"], _SCAN["aff"], _SCAN["lw"],
                                              _SCAN["Xall"], _SCAN["DETECTED"], _SCAN["ab_all"])
    min_partial, min_scanmir, alpha = _SCAN["min_partial"], _SCAN["min_scanmir"], _SCAN["alpha"]
    batch, n_null_arms, with_null, permute = _SCAN["batch"], _SCAN["n_null_arms"], _SCAN["with_null"], _SCAN["permute"]
    rows, nulls = [], []
    try:
        Y, Xo, C, w = LD.assemble_gene(g, w_prior_source="ledger", orphans=True)
    except Exception:
        return rows, nulls
    if permute is not None:
        Y = pd.Series(np.random.default_rng(abs(hash(("perm", g))) % (2**32)).permutation(Y.to_numpy()), index=Y.index)
    he_arms = [a for a in hemap.get(g, set()) if a in Xo.columns]
    orphans = [a for a in Xo.columns if a not in set(he_arms)]
    if not orphans or len(he_arms) < 1:
        return rows, nulls
    affg = aff.loc[aff["gene"] == g].set_index("arm")["repression"]
    cand = [a for a in orphans if a in affg.index and affg[a] < min_scanmir]   # scanMiR prefilter (cheap)
    if not cand and not with_null:
        return rows, nulls
    Cm = C.to_numpy(float)
    Xhe = Xo[he_arms]
    # CURATED-HE AGGREGATE = the CANONICAL model's predicted curated repression (MH-155, user-caught).
    # he_agg is a single "curated pressure" covariate added to Cext so the orphan's coupling is residual on
    # what curation already explains. It must be the CANONICAL Gibbs estimate, NOT `LR.fit_gene` — that is the
    # §6b-RETIRED adaptive lasso (memory: Lasso RETIRED; MH-145 flagged this exact pattern) and was also the
    # ~12 s/gene bottleneck.
    #   `_SCAN["card"]` = dense posterior mean β per (gene, arm) from arm_attribution_card.tsv (100% HE coverage
    #   on the discovery universe). ⚠ SCALE: that β is the coefficient on the model's `_prep` design — X
    #   **z-scored** and C-residualised (attribution_eb._prep) — so β must be applied to a z-scored Xhe, not raw
    #   log2RPM. We z-score Xhe here on the SAME convention; the curated β was fit WITHOUT the orphans present,
    #   so there is no leakage into the orphan estimate. Fallback to the lasso only if a gene is absent from the
    #   card (0% in practice).
    cardg = _SCAN["card"].get(g)
    if cardg is not None:
        _, Xhe_z, hz_cols = AE._prep(Y, Xhe, C)                          # z-scored, C-residualised (β's scale)
        beta_he = np.array([cardg.get(a, 0.0) for a in hz_cols], float)
        he_agg = Xhe_z @ beta_he                                         # canonical predicted curated repression
    else:
        M_he = LR.fit_gene(Y, Xhe, C, w.reindex(he_arms), alpha=alpha)
        he_agg = Xhe.to_numpy(float) @ M_he.reindex(Xhe.columns).fillna(0).to_numpy()
    # BATCH DELIBERATELY OFF (MH-155, user-asked twice). The canonical C = [CPE, target_cn, mal_prolif] and
    # the canonical model / MH-123's calibration-standard null carry NO batch — so he_agg's β was fit in a
    # no-batch world and adding batch here scores against a world it never saw. And batch is SHARED population
    # structure ⇒ the site-free arms feel the same batch as real edges ⇒ it is already IN the null;
    # residualising it narrows both classes equally (σ₀ 0.094→0.087, ~7%) without improving calibration.
    Cext = np.column_stack([Cm, he_agg, _batch_np(Y.index)]) if batch else np.column_stack([Cm, he_agg])
    # ⚡ ONE residualiser per gene (MH-155 optim). The projection onto Cext is the SAME for every arm, so
    # build Q = QR([1, Cext]) once and residualise any block as V − Q(QᵀV) — a single matmul — instead of the
    # per-arm sklearn `_resid` + `spearmanr` (which re-fit OLS and re-rank on every candidate). Estimator is
    # UNCHANGED: spearmanr(_resid(x,Cext), yr) ≡ pearson(rank(Q-resid x), rank(Q-resid Y)); verified bit-equal.
    Qc, _ = np.linalg.qr(np.column_stack([np.ones(len(Cext)), Cext]))
    def _rblock(V):                                            # residualise columns of V on Cext (with intercept)
        return V - Qc @ (Qc.T @ V)
    yr = _rblock(Y.to_numpy(float)[:, None]).ravel()
    yr_rank = rankdata(yr)
    lwg = lw.loc[lw["gene"] == g].set_index("arm")["ledger_weight"]
    fam = _SCAN["fam"]
    he_fams = {fam.get(h) for h in he_arms} - {None}
    # FLAGS (MH-155, user-asked). `he_max_corr` = max |Spearman| of the orphan vs the gene's curated HE arms,
    # RAW (pre-residualisation) — he_agg residualises out the curated *aggregate*, but this flags a candidate
    # that is a near-duplicate of one SPECIFIC known regulator. `same_seed_he` = the gene already has a curated
    # HE edge to another arm of this candidate's seed family (a paralogue of a known regulator, not novel).
    Rhe = rankdata(Xhe.to_numpy(float), axis=0) if he_arms else None       # samples × n_he
    if Rhe is not None:
        Rhe = Rhe - Rhe.mean(0)
        Rhe_norm = np.sqrt((Rhe * Rhe).sum(0))
    if cand:
        # VECTORISED candidate coupling: residualise all candidate arms at once, rank, correlate with yr-rank.
        Vc = Xo[cand].to_numpy(float)
        Vcr = rankdata(_rblock(Vc), axis=0)
        rho_c = _rho_block(Vcr, yr_rank)
        Rc = rankdata(Vc, axis=0); Rc = Rc - Rc.mean(0)                    # raw ranks for he_max_corr
        for j, a in enumerate(cand):
            rho = rho_c[j]
            if np.isfinite(rho) and rho < min_partial:
                hmc = np.nan
                if Rhe is not None:
                    ra = Rc[:, j]; den = np.sqrt((ra * ra).sum()) * Rhe_norm
                    with np.errstate(invalid="ignore", divide="ignore"):
                        cc = np.where(den > 0, (Rhe * ra[:, None]).sum(0) / den, 0.0)
                    hmc = round(float(np.abs(cc).max()), 3) if len(cc) else np.nan
                rows.append({"gene": g, "arm": a, "partial_coupling": round(float(rho), 3),
                             "scanmir_rep": round(float(affg[a]), 2),
                             "sub_he_evidence": round(float(lwg.get(a, 0.0)), 2),
                             "arm_abundance": round(float(ab_all.get(a, np.nan)), 3),
                             "seed_family": fam.get(a), "same_seed_he": bool(fam.get(a) in he_fams),
                             "he_max_corr": hmc, "_is_null": False})
    if with_null and permute is None:
        sited = set(affg.index[affg < min_scanmir]) | set(Xo.columns)   # duplex OR TargetScan site OR curated
        pool = [a for a in DETECTED if a not in sited]
        if len(pool) >= 5:
            rng = np.random.default_rng((_SCAN["null_seed"], abs(hash(g)) % (2**32)))
            take = [pool[j] for j in rng.choice(len(pool), size=min(n_null_arms, len(pool)), replace=False)]
            V = Xall.loc[take, Y.index].T.astype(float).fillna(0.0).to_numpy()      # samples × arms
            ok = V.std(0) > 1e-9
            if ok.any():
                V, take = V[:, ok], [t for t, k in zip(take, ok) if k]
                Vr = rankdata(_rblock(V), axis=0)                          # same residualiser as the candidates
                rho_n = _rho_block(Vr, yr_rank)
                for a, r in zip(take, rho_n):
                    if np.isfinite(r):
                        nulls.append({"gene": g, "arm": a, "partial_coupling": round(float(r), 4),
                                      "arm_abundance": round(float(ab_all.get(a, np.nan)), 3), "_is_null": True})
    return rows, nulls


def scan_all(genes=None, *, permute: int | None = None, min_partial: float = -0.15,
             min_scanmir: float = -0.5, alpha: float = 0.01, progress: int = 200,
             batch: bool = False, n_null_arms: int = 40, null_seed: int = 0,
             with_null: bool = True, workers: int | None = None) -> pd.DataFrame:
    """ONE assemble per gene: orphan partial coupling (beyond curated HE aggregate) + scanMiR prefilter.
    No deconv (that's applied to survivors afterwards). ⚡ PARALLEL over genes (MH-155) — each gene
    re-assembles its own design, so the loop is embarrassingly parallel (~2.2 s/gene × 1,571 ≈ 58 min serial).

    **MH-155 — THE NULL IS FITTED IN THIS LOOP, ON EACH GENE'S OWN `Cext`.** For each gene we also score
    `n_null_arms` SITE-FREE arms — arms that are (a) not curated for this gene, (b) carry no TargetScan
    context++ site (⇒ absent from `assemble_gene(orphans=True)`'s X), and (c) carry no genome-wide scanMiR
    duplex (`repression > min_scanmir`). They **cannot repress by any modelled mechanism** yet pass through the
    IDENTICAL estimator, covariates (C + HE-aggregate + batch), gene, and sample support. Their ρ distribution
    IS this lane's null.

    Why in-loop rather than calling `eval.site_free_null.fit()`: that module residualises on **C only**, while
    this lane residualises on **C + he_agg + batch**. Scoring one against the other mis-scales the null — the
    exact defect MH-154 documented. Same estimand or no null.

    `permute` (raw-Y shuffle) is **DEPRECATED** — a single draw, and per MH-154 permutation reproduces the
    theoretical null because it destroys the unmodelled structure that widens the real one. `with_null=True`
    is the calibrated path.

    Returns candidate rows; site-free null rows are carried alongside under `_is_null=True`.
    """
    genes = genes or sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))  # POOLED-HE gene universe (migration)
    Xall = LD._load()["X"]
    # The site-free null POOL must mirror the candidate arms' detection profile. assemble_gene(orphans=True)
    # admits only curated/TargetScan-site arms — a DETECTED set — and merely fills their occasional gaps
    # (data.py:447 `.fillna(0.0)`). A site-free arm drawn from the full matrix is often near-all-zeros
    # (only 9.1% of arms are all-finite; the rest are sparsely detected); such arms carry almost no rank
    # variance, so their partial correlations collapse toward 0 and DEFLATE the null (robust sd 0.045 vs the
    # honest 0.083–0.13). Restrict to arms detected in >50% of samples — the SAME population the eval-module
    # null uses (`eval.site_free_null`) — then impute the residual gaps as the regression does.
    DETECTED = [a for a in Xall.index if np.isfinite(Xall.loc[a].to_numpy(float)).mean() > 0.5]
    # Canonical Gibbs β per (gene, arm) for the HE aggregate control — see _scan_one_gene. gene -> {arm: β}.
    try:
        _card = pd.read_csv(C_CARD, sep="\t")
        card = {g: dict(zip(grp.arm, grp.beta)) for g, grp in _card.groupby("gene")}
    except FileNotFoundError:
        print(f"[scan] ⚠ {C_CARD} absent — HE aggregate falls back to the RETIRED lasso (slow, non-canonical)")
        card = {}
    from mirna_hallmark.coupling_inference import seed_family_map
    _SCAN.clear()
    _SCAN.update(dict(
        hemap=_he_arms_map(), aff=KD.genome_affinity(), lw=LG.edge_weights(), Xall=Xall, card=card,
        fam=seed_family_map(list(Xall.index)),                 # arm -> seed family (for same-seed-HE flag)
        DETECTED=DETECTED, ab_all=Xall.loc[DETECTED].fillna(0.0).median(axis=1),
        min_partial=min_partial, min_scanmir=min_scanmir, alpha=alpha, batch=batch,
        n_null_arms=n_null_arms, with_null=with_null, permute=permute, null_seed=null_seed))
    rows, nulls = [], []
    if workers is None:                       # leave 2 cores for the parent + OS (8-core box → 6 workers);
        import os                             # oversubscribing all cores starved the parent (~2× not 8×).
        workers = max(1, min(16, (os.cpu_count() or 4) - 2))
    if workers and workers > 1 and len(genes) > 1:
        import multiprocessing as mp
        with mp.Pool(workers) as pool:
            for i, (r, n) in enumerate(pool.imap_unordered(_scan_one_gene, genes, chunksize=4)):
                if progress and i % progress == 0:
                    print(f"[scan{'·permnull' if permute else ''}] {i}/{len(genes)} genes, "
                          f"{len(rows)} hits, {len(nulls)} null draws", flush=True)
                rows.extend(r); nulls.extend(n)
    else:
        for i, g in enumerate(genes):
            if progress and i % progress == 0:
                print(f"[scan] {i}/{len(genes)} genes, {len(rows)} hits, {len(nulls)} null draws", flush=True)
            r, n = _scan_one_gene(g)
            rows.extend(r); nulls.extend(n)
    out = pd.DataFrame(rows + nulls)
    return out if len(out) else pd.DataFrame(columns=["gene", "arm", "partial_coupling", "_is_null"])


def deconv_validate(cand: pd.DataFrame, *, min_partial: float = -0.15, alpha: float = 0.01,
                    progress: int = 50, batch: bool = False) -> pd.DataFrame:
    """Composition-robustness for candidate (gene, arm) edges — **ONE deconv assemble per gene**, validating
    all of that gene's candidates in the same pass (not a re-`discover_gene` per gene). Adds `partial_deconv`
    + `robust` (retains >60% of the coupling under CIBERSORTx non-malignant fractions)."""
    # canonical DECONV he_agg = card `beta_deconv` (NOT the retired lasso, matching scan_all); QR residualiser.
    try:
        _card = pd.read_csv(C_CARD, sep="\t")
        cardd = {g: dict(zip(grp.arm, grp.beta_deconv)) for g, grp in _card.groupby("gene")}
    except FileNotFoundError:
        cardd = {}
    hemap = _he_arms_map()
    out = []
    for i, (g, grp) in enumerate(cand.groupby("gene")):
        if progress and i % progress == 0:
            print(f"[deconv] {i}/{cand['gene'].nunique()} genes", flush=True)
        try:
            Yd, Xod, Cd, wd = LD.assemble_gene(g, w_prior_source="ledger", orphans=True, deconv=True)
            _, Xhed, _, whed = LD.assemble_gene(g, w_prior_source="ledger", deconv=True)
        except Exception:
            continue
        Cmd = Cd.to_numpy(float)
        cg = cardd.get(g)
        if cg is not None:
            _, Xhez, hz = AE._prep(Yd, Xhed, Cd)                 # deconv β's scale (z-scored, deconv-C-residualised)
            aggd = Xhez @ np.array([cg.get(a, 0.0) for a in hz], float)
        else:
            aggd = Xhed.to_numpy(float) @ LR.fit_gene(Yd, Xhed, Cd, whed, alpha=alpha).reindex(Xhed.columns).fillna(0).to_numpy()
        Cextd = np.column_stack([Cmd, aggd, _batch_np(Yd.index)]) if batch else np.column_stack([Cmd, aggd])
        Qd, _ = np.linalg.qr(np.column_stack([np.ones(len(Cextd)), Cextd]))
        yrd = rankdata((Yd.to_numpy(float) - Qd @ (Qd.T @ Yd.to_numpy(float))))
        arms = [a for a in grp["arm"] if a in Xod.columns]
        if arms:
            Vd = Xod[arms].to_numpy(float)
            Vdr = rankdata(Vd - Qd @ (Qd.T @ Vd), axis=0)
            rd = dict(zip(arms, _rho_block(Vdr, yrd)))
        else:
            rd = {}
        for _, r in grp.iterrows():
            a = r["arm"]
            pdc = rd.get(a, np.nan)
            row = r.to_dict()
            row["partial_deconv"] = round(float(pdc), 3) if pdc == pdc else np.nan
            ret = (pdc / r["partial_coupling"]) if r["partial_coupling"] else np.nan
            row["retention"] = round(float(ret), 2) if ret == ret else np.nan
            # composition-robust = retains ≥60% of coupling under deconv AND still couples (not a stroma edge)
            row["robust"] = bool(ret >= 0.6 and pdc < -0.1)
            out.append(row)
    return pd.DataFrame(out).sort_values("partial_coupling") if out else pd.DataFrame()


def calibrate(scan: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score candidates against the in-loop SITE-FREE null (MH-155). Efron form: FIT the null's location+scale
    per arm-abundance quintile (robust median / 1.4826·MAD), score each candidate against N(μ₀,σ₀) → continuous
    p. Exceedance counting is resolution-limited (min p = 1/N) and can never fire BH — see MH-123.

    **TWO senses of "family", both reported (MH-155, matching MH-45):**
      * `q_by_arm`  — BY across ALL arm-edges. **BY not BH** because arm-edges are dependent (shared Y/C within
        a gene, shared population structure across genes) and BH is only valid under positive dependence.
      * `q_seedfamily` — the DE-DUPLICATED hypothesis. Paralogues sharing a seed are non-independent, so each
        (gene, seed-family) is collapsed to ONE **Simes** p, then BY across those. Simes is NOT itself an FDR —
        it reduces the multiplicity (arm-edges → gene×family); BY is the FDR on top. This is the honest count.

    ⚠ **Any q-COUNT is a THRESHOLD statistic and is FRAGILE** (axiom 5): mass piles at the boundary, so a few-%
    shift in σ₀ moves the count severalfold. `null_z` and the set-level shift are the robust readouts.
    """
    from mirna_hallmark import stats as S
    from mirna_hallmark.coupling_inference import benjamini_yekutieli, simes_p
    from mirna_hallmark.eval.site_free_null import _fit_bins, N_BINS
    from scipy.stats import norm

    nul = scan[scan["_is_null"]].rename(columns={"partial_coupling": "rho", "arm_abundance": "ab"})
    real = scan[~scan["_is_null"]].copy()
    nul = nul.dropna(subset=["rho", "ab"])
    if len(nul) < 250 or not len(real):
        print(f"[calibrate] site-free null too small (n={len(nul)}) — refusing to calibrate")
        return real, pd.DataFrame()
    bins, params = _fit_bins(nul)
    bi = np.clip(np.digitize(real["arm_abundance"].to_numpy(float), bins) - 1, 0, N_BINS - 1)
    mu = np.array([params[b][0] for b in bi]); sd = np.array([params[b][1] for b in bi])
    real["null_mu0"], real["null_sd0"] = mu.round(4), sd.round(4)
    real["null_z"] = ((real["partial_coupling"].to_numpy(float) - mu) / sd).round(3)
    p = norm.cdf(real["null_z"].to_numpy(float))
    real["p_sitefree"] = p
    real["q_by_arm"] = benjamini_yekutieli(np.nan_to_num(p, nan=1.0)).round(4)      # arm-edge, dependence-robust
    # seed-family collapse: one Simes p per (gene, seed-family), BY across families, mapped back to members
    key = list(zip(real["gene"], real["seed_family"].fillna("NA_" + real["arm"].astype(str))))
    fam_p = {}
    for k in set(key):
        ps = p[[kk == k for kk in key]]
        fam_p[k] = simes_p(ps)
    fk = sorted(fam_p)
    fq = dict(zip(fk, benjamini_yekutieli([fam_p[k] for k in fk])))
    real["p_seedfamily"] = [round(float(fam_p[k]), 6) for k in key]
    real["q_seedfamily"] = [round(float(fq[k]), 4) for k in key]
    if "same_seed_he" in real.columns:
        real["same_seed_he"] = real["same_seed_he"].astype(bool)   # object→bool so ~ counts, not bitwise-negates
    par = pd.DataFrame([{"bin": b, "ab_lo": bins[b], "ab_hi": bins[b + 1],
                         "mu0": params[b][0], "sd0": params[b][1],
                         "n_null": int(((nul.ab >= bins[b]) & (nul.ab < bins[b + 1])).sum())}
                        for b in range(N_BINS)])
    return real.sort_values("null_z"), par


def attach_evidence(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate each discovery with the WEAKER evidence it carries (MH-155, user-asked) — beyond the scanMiR /
    TargetScan PRECONDITIONS (every candidate already has both a context++ site and a K_D duplex; `scanmir_rep`
    is present, `ts_mag` added here). Attached, NOT gated:
      * `ts_mag`        — TargetScan context++ magnitude (larger = stronger predicted repression).
      * `chimeric_wt` / `chimeric_src` — DIRECT miRNA↔target duplex from `chimeric_evidence` (Manakov eCLIP +
        TarBase CLASH/qCLASH), the one source type that RESOLVES the mature arm. ⚠ NO breast chimeric exists in
        any source — cross-tissue corroboration, not breast-specific.
      * `ev_classes` / `ev_npmid` — the ledger's per-edge assay classes and #DISTINCT PMIDs, **deduplicated
        across mirTarBase+TarBase by (edge × PMID × class)** (that is the ledger's whole job). This surfaces the
        higher-throughput / lower-tier experiments (ago_clip, qpcr_rna, degradome, …) that sit below HE.
    """
    from mirna_hallmark.learned import chimeric_evidence as CE
    from mirna_hallmark.learned.evidence import ledger as _LG
    df = df.copy()
    # ts_mag
    ts = LD._targetscan_context()
    if ts is not None:
        tsm = ts.set_index(["arm", "gene"])["ts_mag"]
        df["ts_mag"] = [round(float(tsm.get((a, g), 0.0)), 3) for a, g in zip(df["arm"], df["gene"])]
    # ledger per-class (deduped) — one long table keyed by (arm, gene, pmid, assay_class)
    L = _LG.build_ledger()
    Lg = {k: v for k, v in L.groupby(["arm", "gene"])}
    cls, npm = [], []
    for a, g in zip(df["arm"], df["gene"]):
        sub = Lg.get((a, g))
        if sub is None or not len(sub):
            cls.append(""); npm.append(0)
        else:
            cls.append(",".join(sorted(sub["assay_class"].dropna().unique())))
            npm.append(int(sub["pmid"].dropna().nunique()))
    df["ev_classes"], df["ev_npmid"] = cls, npm
    # chimeric direct-duplex weight (per gene, over that gene's candidate arms)
    cw, cs = {}, {}
    for g, grp in df.groupby("gene"):
        try:
            ev = CE.evidence(g, list(grp["arm"].unique()))            # per-arm pooled chimeric weight
            mat = CE.evidence_matrix(g, list(grp["arm"].unique()))    # per-source, to name the provenance
        except Exception:
            ev, mat = None, None
        for a in grp["arm"].unique():
            cw[(g, a)] = float(ev.get(a, 0.0)) if ev is not None else 0.0
            if mat is not None and a in getattr(mat, "index", []):
                srcs = [c for c in mat.columns if float(mat.loc[a, c]) > 0]
                cs[(g, a)] = ",".join(srcs)
            else:
                cs[(g, a)] = ""
    df["chimeric_wt"] = [round(cw.get((g, a), 0.0), 3) for a, g in zip(df["arm"], df["gene"])]
    df["chimeric_src"] = [cs.get((g, a), "") for a, g in zip(df["arm"], df["gene"])]
    return df


def run_all(*, min_partial: float = -0.15, min_scanmir: float = -0.5, top: int | None = None,
            n_null_arms: int = 40, out="mirna_hallmark/output/learned/discoveries.tsv") -> pd.DataFrame:
    """Genome-wide discovery, calibrated against the in-loop site-free null (MH-155).

    ⚠ **RESERVE "DISCOVERY" LANGUAGE FOR THE AGGREGATE LANE** (MH-123). Per-edge miRNA→target discovery in
    TCGA bulk at n≈1041 is essentially impossible: the honest null is ~2.5–3.7× wider than the theoretical
    one, so per-edge q-values here are a calibrated *ranking*, not a licence for a per-edge claim.
    """
    from pathlib import Path
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    scan = scan_all(min_partial=min_partial, min_scanmir=min_scanmir, n_null_arms=n_null_arms)
    real, par = calibrate(scan)
    real = attach_evidence(real)
    scan[scan["_is_null"]].to_csv(Path(out).parent / "discovery_sitefree_null.tsv", sep="\t", index=False)
    real.to_csv(Path(out).parent / "discovery_candidates.tsv", sep="\t", index=False)   # persist BEFORE deconv
    if len(par):
        print(f"[run_all] site-free null fitted on {int(par.n_null.sum()):,} arm-gene draws "
              f"(σ₀ {par.sd0.min():.4f}–{par.sd0.max():.4f}, μ₀ {par.mu0.min():+.4f}…{par.mu0.max():+.4f})")
        print(par.to_string(index=False))
        n_arm = int((real["q_by_arm"] < 0.05).sum())
        fam_sig = real[real["q_seedfamily"] < 0.05].groupby(["gene", "seed_family"]).ngroups
        n_novel = int((~real["same_seed_he"]).sum())
        print(f"[run_all] {len(real)} arm-candidates past ρ<{min_partial} | {n_arm} q_by_arm<0.05 | "
              f"{fam_sig} distinct (gene,seed-family) at q_seedfamily<0.05 | median null_z {real['null_z'].median():+.2f}")
        print(f"[run_all] {n_novel} candidates are NOT paralogues of a curated regulator (same_seed_he=False)")
        print("[run_all] ⚠ q-counts are threshold-FRAGILE (axiom 5); read null_z / the seed-family collapse.")
    cand = real.nsmallest(top, "null_z") if top else real
    dv = deconv_validate(cand, min_partial=min_partial)                                  # single-pass
    if len(dv):
        dv.to_csv(out, sep="\t", index=False)
        rob = dv[dv["robust"]] if "robust" in dv else dv
        print(f"[run_all] {len(dv)} validated → {len(rob)} deconv-robust across {rob['gene'].nunique()} genes "
              f"| {int((rob['sub_he_evidence']==0).sum())} fully novel | wrote {out}")
        print(rob.sort_values('partial_deconv').head(15).to_string(index=False))
    return dv


def run(genes=None, *, min_partial: float = -0.12) -> pd.DataFrame:
    genes = genes or PANEL
    out = pd.concat([discover_gene(g, min_partial=min_partial) for g in genes], ignore_index=True) \
        if genes else pd.DataFrame()
    if out.empty:
        print("no discoveries at this threshold")
        return out
    out = out.sort_values("partial_coupling")
    with pd.option_context("display.width", 160):
        print(f"=== Discovered orphan edges (partial coupling < {min_partial} beyond curated HE + scanMiR-supported) ===")
        print(out.to_string(index=False))
    nr = int(out["robust"].sum()) if "robust" in out else 0
    print(f"\n{len(out)} candidates across {out['gene'].nunique()} genes | {nr} SURVIVE composition control "
          f"(deconv) | {int((out['sub_he_evidence'] == 0).sum())} fully novel (sequence+expression+K_D only)")
    if nr:
        print("robust discoveries:", ", ".join(f"{r.arm}→{r.gene}" for r in out[out['robust']].itertuples()))
    return out


if __name__ == "__main__":
    run()
