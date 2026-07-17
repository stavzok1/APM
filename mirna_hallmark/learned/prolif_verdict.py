"""Per-gene proliferation confound-vs-mechanism VERDICT (MH-100).

The confounder block C's `mal_prolif` (E2F+G2M metagene, composition-adjusted) is materially INCOMPLETE:
residual proliferation R = resid(broad-cell-cycle metagene | C) still explains |rho|~0.3-0.8 of most targets'
residual variance. But "control proliferation more completely" is WRONG-headed globally, because proliferation
is BOTH confounder AND mechanism (miR-15/16->CCND1, let-7->CDK6 act THROUGH proliferation). The right question
is per-gene appropriateness, and it is decided by a NON-CIRCULAR held-out OOF 2x2:

  train {C, C+R} x evaluate {C-space, C+R-space}  ->  does controlling R help held-out coupling?
    control HELPS both spaces      -> R is a genuine CONFOUND        -> control is appropriate  (e.g. CDKN1B)
    control HURTS both spaces      -> R is the MECHANISM (mediator)  -> controlling = over-control (e.g. BCL2/ZEB1)
    winner FLIPS with eval space   -> estimand-dependent             -> FRAGILE, OOF cannot settle it

This is a DIAGNOSTIC/annotation, not a change to the production estimator: C stays as-is (global control rejected;
the design's clean fix -- AL2 intronic nascent-transcription latent -- is data-blocked, no intronic reads).
Emits per-gene {rho_resid, delta_Cspace, delta_CRspace, class, mean_rel_dbeta} so downstream readouts know which
couplings are proliferation-fragile. Sibling of the coding-pleiotropy gate's "flag, don't subtract" philosophy.

CLI: `.venv/bin/python3 -m mirna_hallmark.learned.prolif_verdict [GENE ...]`  (default = representative panel)
"""
from __future__ import annotations
import sys, os, numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import KFold
from mirna_hallmark.learned import data as LD, spike_slab as SS, families as FAM, regression as LR
from mirna_hallmark.hallmark_sets import HallmarkSets

OUT = "mirna_hallmark/output/learned/prolif_verdict.tsv"               # 28-gene panel, 5-fold (high-quality ref)
EDGE_OUT = "mirna_hallmark/output/learned/prolif_verdict_edges.tsv"
GW_OUT = "mirna_hallmark/output/learned/prolif_verdict_genomewide.tsv"  # 1507-gene screen (4-fold)
GW_EDGE_OUT = "mirna_hallmark/output/learned/prolif_verdict_edges_genomewide.tsv"
TOL = 0.02                                                    # Spearman points; below this = "neutral"/no-verdict
TOL_E = 0.15                                                  # per-edge rel-Δβ magnitude for a directional edge call

# broad cell-cycle metagene (spans MKI67/PCNA/mitotic beyond the E2F+G2M mean already in mal_prolif)
_PROLIF = ["MKI67","PCNA","TOP2A","CCNB1","CCNB2","AURKA","AURKB","BUB1","CDK1","CDC20",
           "CCNA2","FOXM1","PLK1","RRM2","TYMS","UBE2C","BIRC5","CENPF","KIF11","MCM2"]
_EFFECTORS = {"CCND1","CCNE1","CCNE2","CDK4","CDK6","CDC25A","MYB","FOXM1","E2F1","E2F2","E2F3",
              "CDK1","CCNB1","CCNA2","MCM2","MCM7","PCNA","BUB1","AURKA","AURKB","TOP2A","MKI67"}
# canonical anti-proliferative miRNA families (matched by miR-number substring on the seed-family column)
_ANTIPROLIF = ["miR-15","miR-16","miR-195","miR-424","miR-497","let-7","miR-34","miR-449",
               "miR-192","miR-215","miR-193","miR-26","miR-503","miR-98","miR-100","miR-99"]

def _prolif_regulator(col):
    return any(tok in col for tok in _ANTIPROLIF)
PANEL = ["CCND1","CCNE1","CCNE2","CDK6","CDK4","E2F1","E2F3","FOXM1","MYB","MYC","CCNB1","CDC25A","BUB1","AURKA",
         "PTEN","ZEB1","ZEB2","ESR1","GATA3","BCL2","VEGFA","THBS1","CDKN1A","CDKN1B","SPARC","SMAD7","PDCD4",
         "VIM","CDH1","TGFBR2","RB1","BCL2L11"]

def _resid(V, C):
    C1 = np.column_stack([np.ones(len(C)), C]); b, *_ = np.linalg.lstsq(C1, V, rcond=None); return V - C1 @ b

def _prolif_set():
    hs = HallmarkSets.load()
    return (set(hs.sets.get("HALLMARK_E2F_TARGETS", [])) | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", [])) |
            set(hs.sets.get("HALLMARK_MYC_TARGETS_V1", [])) | _EFFECTORS)

def _residual_prolif(gene, C, Yfull, parts):
    """R = residual-proliferation axis orthogonal to C (target gene EXCLUDED from the metagene)."""
    pg = [g for g in _PROLIF if g in Yfull.index and g != gene]
    P = Yfull.loc[pg, parts].astype(float)
    P = P.sub(P.mean(1), axis=0).div(P.std(1) + 1e-9, axis=0).mean(0)
    return _resid(P.to_numpy(float), C.to_numpy(float))

def _oof_pred(Yg, X, C, w, folds=5, n_iter=900, seed=0):
    n = len(Yg); pred = np.full(n, np.nan)
    for tr, te in KFold(folds, shuffle=True, random_state=seed).split(X):
        M, _ = SS.fit_gene_ss(Yg.iloc[tr], X.iloc[tr], C.iloc[tr], w, seed=seed, n_iter=n_iter, burn=n_iter // 3)
        pred[te] = LR.aggregate(X.iloc[te], M)
    return pred

def _classify(dC, dCR):
    if (dC > TOL and dCR < -TOL) or (dC < -TOL and dCR > TOL):
        return "fragile"                                     # winner flips with eval space -> undetermined
    davg = 0.5 * (dC + dCR)
    if davg > TOL:  return "confound"                        # control helps -> genuine confound
    if davg < -TOL: return "over_control"                   # control hurts -> proliferation is the mechanism
    return "neutral"

def _edge_verdict(gene_class, rel):
    """Per-edge call: gene OOF class = causal context; edge Δβ direction = the refinement within it."""
    if rel > TOL_E:                                          # edge GAINS under control -> prolif was HIDING it
        return "unmasked_real"
    if rel < -TOL_E:                                         # edge SHRINKS under control
        if gene_class == "confound":     return "deconfound"     # spurious prolif-driven edge, correctly removed
        if gene_class == "over_control": return "mechanism"      # coupling runs THROUGH prolif -> don't strip
        return "ambiguous"                                       # fragile/neutral gene -> can't tell
    return "prolif_robust"                                   # stable -> insensitive to prolif control (cleanest)

def gene_verdict(gene, prolif_set=None, Yfull=None, *, folds=5, n_iter=900):
    prolif_set = prolif_set if prolif_set is not None else _prolif_set()
    Yfull = Yfull if Yfull is not None else LD._load()["Y"]
    Yg, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    R = _residual_prolif(gene, C, Yfull, Yg.index)
    Caug = C.copy(); Caug["resid_prolif"] = R
    yv = Yg.to_numpy(float); Cc, Ca = C.to_numpy(float), Caug.to_numpy(float)
    rho_resid = spearmanr(R, _resid(yv, Cc)).correlation      # residual prolif still in the target
    # OOF 2x2 — REPRESSION-DIRECTED coupling = −ρ (repression ⇒ pred anti-correlates target ⇒ ρ<0 ⇒ −ρ>0). NOT abs():
    # abs conflated real repression-strengthening (−0.35→−0.47) with sign-FLIPS to anti-repression (−0.04→+0.16). Fixed
    # 2026-07-11 (user-caught, via the host-gene lens); persisted prolif_verdict{,_genomewide,_edges}.tsv REGENERATED
    # under −ρ same day (headline robust: GW over_control 389 vs confound 170 still 2.3:1, net Δ_Cspace −0.017; panel
    # confound 5→3). The regen also picked up the live _arm_name_map universe correction (GW edges 1191→1242).
    pC  = _oof_pred(Yg, X, C, w, folds, n_iter)
    pCR = _oof_pred(Yg, X, Caug, w, folds, n_iter)
    def rho(pred, S): return -spearmanr(_resid(pred, S), _resid(yv, S)).correlation
    MC_C, MCR_C   = rho(pC, Cc), rho(pCR, Cc)
    MC_CR, MCR_CR = rho(pC, Ca), rho(pCR, Ca)
    dC, dCR = MCR_C - MC_C, MCR_CR - MC_CR                     # >0 = controlling R HELPS held-out coupling
    gclass = _classify(dC, dCR)
    # full-data beta movement (how much each coupling shifts under R) -> the per-EDGE layer
    M_base, _ = SS.fit_gene_ss(Yg, X, C, w, gene=gene, n_iter=n_iter, burn=n_iter // 3)
    M_aug, _  = SS.fit_gene_ss(Yg, X, Caug, w, gene=gene, n_iter=n_iter, burn=n_iter // 3)
    coupled = M_base[M_base > 0.005]
    rel = ((M_aug.reindex(coupled.index) - coupled) / (coupled + 1e-6))
    edges = [dict(gene=gene, family=f, g_prolif=gene in prolif_set, prolif_reg=_prolif_regulator(f),
                  beta_base=round(float(coupled[f]), 4), beta_aug=round(float(M_aug[f]), 4),
                  rel=round(float(rel[f]), 3), rho_resid=round(float(rho_resid), 3),
                  gene_verdict=gclass, edge_verdict=_edge_verdict(gclass, float(rel[f]))) for f in coupled.index]
    grow = dict(gene=gene, n=len(Yg), n_edge=int(len(coupled)), g_prolif=gene in prolif_set,
                rho_resid=round(float(rho_resid), 3), coupling_C=round(MC_C, 3), coupling_CR_space=round(MC_CR, 3),
                delta_Cspace=round(float(dC), 3), delta_CRspace=round(float(dCR), 3),
                mean_rel_dbeta=round(float(rel.mean()) if len(rel) else 0.0, 3), verdict=gclass)
    return grow, edges

def run(genes=None, *, folds=5, n_iter=900, quiet=False, out=OUT, edge_out=EDGE_OUT):
    genes = genes or PANEL
    prolif_set, Yfull, rows, erows = _prolif_set(), LD._load()["Y"], [], []
    for i, g in enumerate(genes):
        try:
            grow, edges = gene_verdict(g, prolif_set, Yfull, folds=folds, n_iter=n_iter)
            rows.append(grow); erows.extend(edges)
            if not quiet:
                print(f"  {g:9s} rho_resid={grow['rho_resid']:+.3f} dC={grow['delta_Cspace']:+.3f} "
                      f"dCR={grow['delta_CRspace']:+.3f} -> {grow['verdict']}", flush=True)
            elif i % 100 == 0:
                print(f"  [{i}/{len(genes)}] {g}", flush=True)
        except Exception as e:
            print(f"  [skip {g}] {e}", flush=True)
    df = pd.DataFrame(rows); ed = pd.DataFrame(erows)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    df.to_csv(out, sep="\t", index=False); ed.to_csv(edge_out, sep="\t", index=False)
    print(f"\nwrote {out} ({len(df)} genes) + {edge_out} ({len(ed)} edges)")
    print("gene verdicts:", df["verdict"].value_counts().to_dict())
    if len(ed): print("edge verdicts:", ed["edge_verdict"].value_counts().to_dict())
    print(f"mean |rho_resid|={df.rho_resid.abs().mean():.3f}  mean delta_Cspace={df.delta_Cspace.mean():+.3f}")
    return df, ed

def _coupled_universe():
    from mirna_hallmark.learned.evidence import ledger as LG
    ed = LG.pooled_he_edges(); X = LD._load()["X"]; Y = LD._load()["Y"]
    return sorted(g for g in set(ed["gene"]) & set(Y.index)
                  if ed.loc[ed["gene"] == g, "miRNA"].isin(X.index).any())

# --- parallel genome-wide (fork COW: parent loads data once, workers inherit; BLAS pinned to 1 thread/worker) ---
_GW = {"folds": 3, "n_iter": 700}
def _gw_worker(gene):
    try:
        return ("ok", gene_verdict(gene, _GW["ps"], _GW["Y"], folds=_GW["folds"], n_iter=_GW["n_iter"]))
    except Exception as e:
        return ("err", (gene, str(e)))

def run_genomewide(*, workers=8, folds=3, n_iter=700):
    """Full coupled universe, PARALLEL over genes. Writes SEPARATE files (`GW_OUT`) so the 28-gene 5-fold
    panel (`OUT`) stays the high-quality reference; borderline genes (|delta|<~0.03) are low-confidence at
    screen resolution. Load data in the PARENT before the pool so forked workers share it copy-on-write."""
    import time
    from concurrent.futures import ProcessPoolExecutor
    _GW["folds"], _GW["n_iter"] = folds, n_iter
    _GW["ps"], _GW["Y"] = _prolif_set(), LD._load()["Y"]; LD._load()["X"]     # warm parent caches pre-fork
    genes = _coupled_universe()
    print(f"genome-wide PARALLEL: {len(genes)} genes · {workers} workers · {folds}-fold n_iter={n_iter}", flush=True)
    rows, erows, t0 = [], [], time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (tag, payload) in enumerate(ex.map(_gw_worker, genes, chunksize=4)):
            if tag == "ok":
                grow, edges = payload; rows.append(grow); erows.extend(edges)
            elif i < 5000:
                print(f"  [skip {payload[0]}] {payload[1]}", flush=True)
            if i % 200 == 0:
                print(f"  [{i}/{len(genes)}] {(time.time()-t0)/60:.1f}min", flush=True)
    df, ed = pd.DataFrame(rows), pd.DataFrame(erows)
    os.makedirs(os.path.dirname(GW_OUT), exist_ok=True)
    df.to_csv(GW_OUT, sep="\t", index=False); ed.to_csv(GW_EDGE_OUT, sep="\t", index=False)
    print(f"\nwrote {GW_OUT} ({len(df)} genes) + {GW_EDGE_OUT} ({len(ed)} edges) in {(time.time()-t0)/60:.1f}min")
    print("gene verdicts:", df["verdict"].value_counts().to_dict())
    print("edge verdicts:", ed["edge_verdict"].value_counts().to_dict())
    print(f"mean |rho_resid|={df.rho_resid.abs().mean():.3f}  mean delta_Cspace={df.delta_Cspace.mean():+.3f}")
    return df, ed

# --- #2: WITHIN-FAMILY arm-level proliferation-SOURCE decomposition (dose-delivery, NOT the verdict) --------------
def arm_prolif_source(gene, family, Yfull=None):
    """For a (gene, family) edge, decompose the family's proliferation correlation to its member ARMS:
    which arm's abundance carries the proliferation signal that confounds the family dose. This grades the
    SOURCE (dose-delivery, §8) — it does NOT subdivide the confound-vs-mechanism VERDICT (that's family-level).
    Per arm: abundance share of the family dose + corr(arm abundance, residual-proliferation R)."""
    Yfull = Yfull if Yfull is not None else LD._load()["Y"]
    Yg, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")            # X = ARM-level (pre-collapse)
    fam_of = FAM.family_of(X.columns)
    members = [a for a in X.columns if fam_of.get(a) == family]
    if not members:
        return pd.DataFrame()
    R = _residual_prolif(gene, C, Yfull, Yg.index)
    lin = np.power(2.0, X[members].to_numpy(float)) - 1.0                    # de-log to RPM for the dose share
    share = lin.mean(0) / (lin.mean(0).sum() + 1e-9)
    rows = []
    for j, a in enumerate(members):
        xa = X[a].to_numpy(float)
        rho = spearmanr(xa, R).correlation if np.std(xa) > 1e-9 else np.nan
        rows.append(dict(gene=gene, family=family, arm=a, abund_share=round(float(share[j]), 3),
                         rho_arm_prolif=round(float(rho), 3) if rho == rho else np.nan,
                         contrib=round(float(share[j]) * abs(float(rho)) if rho == rho else 0.0, 3)))
    return pd.DataFrame(rows).sort_values("contrib", ascending=False)

# intronic miRNA -> its host coding gene (proliferation-relevant hosts; the co-transcription confound source)
_INTRONIC_HOST = {"hsa-miR-93-5p": "MCM7", "hsa-miR-106b-5p": "MCM7", "hsa-miR-25-3p": "MCM7",
                  "hsa-miR-301a-3p": "SKA2", "hsa-miR-338-3p": "AATK", "hsa-miR-33a-5p": "SREBF2",
                  "hsa-miR-128-1": "R3HDM1", "hsa-miR-152-3p": "COPZ2"}

def confound_mechanism(arms=None):
    """Decompose a source arm's proliferation confound into CN vs HOST-GENE transcription. TESTS (not assumes)
    whether the confound is CN-instrument-fixable: for miR-93 the locus CN explains ~5% but the host gene MCM7
    (miR-93/106b~25 is intronic in MCM7, a replication gene) explains ~72% ⇒ the confound is HOST-GENE
    TRANSCRIPTIONAL co-regulation (design axis BP intragenic-host pleiotropy), NOT copy number. So the handle is
    the host mRNA, not the CN instrument."""
    from mirna_hallmark.learned import instrument as IN
    X, Y, L = LD._load()["X"], LD._load()["Y"], IN.locus_cn()
    pg = [g for g in _PROLIF if g in Y.index]
    P = Y.loc[pg]; P = P.sub(P.mean(1), axis=0).div(P.std(1) + 1e-9, axis=0).mean(0)
    arms = arms or list(_INTRONIC_HOST)
    rows = []
    for arm in arms:
        if arm not in X.index:
            continue
        parts = X.columns.intersection(P.index).intersection(Y.columns).intersection(L.index)
        a = X.loc[arm, parts].astype(float); p = P.reindex(parts).astype(float)
        r_ap = spearmanr(a, p).correlation
        # CN leg
        foc = IN._arm_focal_locus(arm); loc = foc.get("locus_id")
        cn_frac = np.nan
        if loc in L.columns and L[loc].std() > 1e-9:
            cn = L.loc[parts, loc].astype(float)
            m = a.notna() & cn.notna() & p.notna()
            cn_frac = 1 - spearmanr(a[m], _resid(p[m].to_numpy(), cn[m].to_numpy().reshape(-1, 1))).correlation / r_ap
        # host-gene leg
        host = _INTRONIC_HOST.get(arm); host_frac = np.nan; r_ah = np.nan
        if host in Y.index:
            h = Y.loc[host, parts].astype(float)
            m = a.notna() & h.notna() & p.notna()
            r_ah = spearmanr(a[m], h[m]).correlation
            host_frac = 1 - spearmanr(a[m], _resid(p[m].to_numpy(), h[m].to_numpy().reshape(-1, 1))).correlation / r_ap
        rows.append(dict(arm=arm, host=host, rho_abund_prolif=round(float(r_ap), 3),
                         rho_abund_host=round(float(r_ah), 3) if r_ah == r_ah else np.nan,
                         cn_mediated_frac=round(float(cn_frac), 2) if cn_frac == cn_frac else np.nan,
                         host_mediated_frac=round(float(host_frac), 2) if host_frac == host_frac else np.nan))
    df = pd.DataFrame(rows)
    out = "mirna_hallmark/output/learned/prolif_confound_mechanism.tsv"
    os.makedirs(os.path.dirname(out), exist_ok=True); df.to_csv(out, sep="\t", index=False)
    print(f"wrote {out} ({len(df)} arms)")
    print(df.to_string(index=False))
    if len(df):
        print(f"\nmedian CN-mediated={df.cn_mediated_frac.median():.2f}  median host-mediated={df.host_mediated_frac.median():.2f}"
              "  => the confound is HOST-GENE transcriptional (axis BP), NOT CN.")
    return df

def arm_source_scan(verdict_edges=None, classes=("confound", "fragile", "mechanism")):
    """Run `arm_prolif_source` over every edge whose gene/edge verdict is proliferation-implicated, to see
    whether the family's confound is ONE arm (fixable, e.g. a proliferation-amplicon arm) or diffuse."""
    ed = verdict_edges if verdict_edges is not None else pd.read_csv(EDGE_OUT, sep="\t")
    hit = ed[ed["edge_verdict"].isin(classes) | ed["gene_verdict"].isin(classes)][["gene", "family"]].drop_duplicates()
    Yfull, out = LD._load()["Y"], []
    for _, r in hit.iterrows():
        d = arm_prolif_source(r["gene"], r["family"], Yfull)
        if len(d) > 1:                                                       # only multi-arm families are interesting
            top = d.iloc[0]; out.append(dict(gene=r["gene"], family=r["family"], n_arm=len(d),
                                             top_arm=top["arm"], top_share=top["abund_share"],
                                             top_rho=top["rho_arm_prolif"], top_contrib=top["contrib"],
                                             concentration=round(float(top["contrib"] / (d["contrib"].sum() + 1e-9)), 2)))
    res = pd.DataFrame(out).sort_values("top_contrib", ascending=False) if out else pd.DataFrame()
    if len(res):
        res.to_csv("mirna_hallmark/output/learned/prolif_arm_source.tsv", sep="\t", index=False)
        print(f"wrote prolif_arm_source.tsv ({len(res)} multi-arm families); median concentration={res.concentration.median():.2f}")
    return res

if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--genomewide":
        run_genomewide()
    elif a and a[0] == "--arm-scan":
        arm_source_scan()
    else:
        run(a or None)
