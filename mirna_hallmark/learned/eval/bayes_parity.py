"""Does the learned-τ² Bayes estimator inherit the lasso's credibility? Parity on the two checks that
GATE credibility, not just OOF coupling: the CPTAC protein OOD (anti-circularity crux) and the
shuffled-evidence null (does the prior do real work). Same harnesses, estimator swapped.

`fit_gene_bayes` mirrors `regression.fit_gene` EXACTLY (same residualisation, adaptive w-scaling, sign) —
only Lasso → learned-τ² dense-ridge posterior mean (`spike_slab._gibbs_posterior`, π≡1). So any difference
is the estimator, nothing else.

    python -m mirna_hallmark.learned.eval.bayes_parity --cptac      # CPTAC protein parity (lasso vs Bayes)
    python -m mirna_hallmark.learned.eval.bayes_parity --shuffle    # shuffled-evidence null parity
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import regression as LR
from mirna_hallmark.learned import spike_slab as SS


def _resid(V, Cmat):
    return V - LinearRegression().fit(Cmat, V).predict(Cmat)


def fit_gene_bayes(Y, X, C, w_prior, *, n_iter: int = 1200, burn: int = 400, **_) -> pd.Series:
    """Learned-τ² dense-ridge posterior-mean M — drop-in for regression.fit_gene (same signature/scale)."""
    Cmat = C.to_numpy(float)
    yr = _resid(Y.to_numpy(float), Cmat)
    Xr = _resid(X.to_numpy(float), Cmat)
    w = np.nan_to_num(w_prior.reindex(X.columns).to_numpy(float), nan=0.0)
    w = np.clip(w, 1e-3, None)
    coef, _, _ = SS._gibbs_posterior(Xr * w, -yr, np.ones(Xr.shape[1]), n_iter=n_iter, burn=burn, seed=0)
    return pd.Series(coef * w, index=X.columns, name="M")


# ─────────────────────────── #3b  CPTAC protein parity ───────────────────────────
def cptac_parity(genes=None, *, cohort: str = "prospective") -> pd.DataFrame:
    """Fit M on TCGA mRNA with BOTH estimators, score BOTH on CPTAC protein via the SAME ood_protein.score.
    Prospective = independent patients + protein layer = the non-circular OOD the credibility claim rests on."""
    from mirna_hallmark.learned.eval import ood_protein as OP
    genes = genes or OP.HUB
    rows = []
    for g in genes:
        try:
            Y, X, C, w = LD.assemble_gene(g, w_prior_source="ledger", deconv=True)  # same setup as OP.fit_M
            M_lasso = LR.fit_gene(Y, X, C, w, alpha=0.005)
            M_bayes = fit_gene_bayes(Y, X, C, w)
            for tag, M in [("lasso", M_lasso), ("bayes", M_bayes)]:
                s = OP.score(g, M, cohort=cohort)
                rows.append({"gene": g, "est": tag, "n": s["n"], "n_reg": s["n_reg"],
                             "rho_protein": s["rho_protein"], "rho_discord": s["rho_discord"],
                             "beats_abund": s["model_beats_abund"]})
        except Exception as e:
            rows.append({"gene": g, "est": "-", "error": repr(e)[:60]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 170):
        print(df.round(3).to_string(index=False))
    ok = df[df.get("rho_protein").notna()] if "rho_protein" in df else df.iloc[0:0]
    print(f"\n=== CPTAC {cohort} ({'INDEPENDENT OOD' if cohort=='prospective' else 'layer-transfer'}) — lasso vs Bayes ===")
    for tag in ["lasso", "bayes"]:
        t = ok[ok["est"] == tag]
        if len(t):
            print(f"  {tag:6s}: protein-coupled {int((t['rho_protein']<0).sum())}/{len(t)} | "
                  f"discordance-coupled {int((t['rho_discord']<0).sum())}/{len(t)} | "
                  f"beats abundance {int(t['beats_abund'].sum())}/{len(t)} | mean ρ_protein {t['rho_protein'].mean():+.3f}")
    # paired: does Bayes match/beat lasso per gene?
    piv = ok.pivot_table(index="gene", columns="est", values="rho_protein")
    if {"lasso", "bayes"}.issubset(piv.columns):
        pv = piv.dropna()
        print(f"  paired: Bayes ρ_protein ≤ lasso on {int((pv['bayes']<=pv['lasso']+1e-6).sum())}/{len(pv)} genes "
              f"| mean Δ(bayes−lasso) {(pv['bayes']-pv['lasso']).mean():+.3f}")
    return df


# ─────────────────────────── #3a  shuffled-evidence null parity ───────────────────────────
def _safe_rho(a, b):
    if np.nanstd(a) == 0 or np.nanstd(b) == 0:
        return np.nan
    return spearmanr(a, b).correlation


def shuffle_null_bayes(gene: str, *, n_shuffle: int = 30, folds: int = 5, seed: int = 0,
                       n_iter: int = 800, burn: int = 300) -> dict:
    """Bar-3 with the Bayes estimator: does the REAL evidence→arm assignment couple better than shuffled?"""
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    Cmat = C.to_numpy(float)
    yr = _resid(Y.to_numpy(float), Cmat)
    splits = list(KFold(n_splits=folds, shuffle=True, random_state=seed).split(X))

    def oof_rho(wv):
        oof = np.full(len(Y), np.nan)
        for tr, te in splits:
            M = fit_gene_bayes(Y.iloc[tr], X.iloc[tr], C.iloc[tr], wv, n_iter=n_iter, burn=burn)
            oof[te] = X.iloc[te].to_numpy(float) @ M.reindex(X.columns).fillna(0).to_numpy()
        return _safe_rho(_resid(oof, Cmat), yr)

    real = oof_rho(w)
    rng = np.random.default_rng(seed)
    nulls = np.array([oof_rho(pd.Series(rng.permutation(w.to_numpy()), index=w.index)) for _ in range(n_shuffle)])
    nulls = nulls[np.isfinite(nulls)]
    p = (1 + int((nulls <= real).sum())) / (len(nulls) + 1) if len(nulls) else np.nan
    z = (real - nulls.mean()) / (nulls.std() + 1e-9) if len(nulls) else np.nan
    return {"gene": gene, "real_rho": round(real, 3), "null_mean": round(float(nulls.mean()), 3),
            "p_perm": round(p, 3), "z": round(float(z), 2), "n_null": len(nulls)}


def shuffle_parity(genes=None, **kw) -> pd.DataFrame:
    genes = genes or ["PTEN", "ESR1", "ZEB1", "GATA3", "CDKN1A"]
    rows = [shuffle_null_bayes(g, **kw) for g in genes]
    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    d = df[df["real_rho"].notna()]
    print(f"\n=== shuffled-evidence null (Bayes) === real beats null-mean on "
          f"{int((d['real_rho']<d['null_mean']).sum())}/{len(d)} | prior-does-work (p<0.1) "
          f"{int((d['p_perm']<0.1).sum())}/{len(d)} | mean z {d['z'].mean():+.2f}")
    return df


# ─────────────────────────── #1  spike-slab-inclusion discovery vs partial-Spearman ───────────────────────────
def discovery_parity(genes=None, *, min_partial: float = -0.12, pip_thr: float = 0.5) -> pd.DataFrame:
    """Can the INCLUSION mode discover orphan edges (the job the dense winner can't do)? Fit a spike-and-slab
    on [HE ∪ orphan] arms with evidence-π — orphans have ~0 ledger weight ⇒ LOW π, so the DATA must overcome
    the prior to pull them in (PIP > thr). Compare the spike-slab orphans to the existing partial-Spearman
    discovery (`discovery.discover_gene`). Overlap ⇒ the inclusion mode recovers the same discoveries."""
    from mirna_hallmark.learned import discovery as DISC
    from mirna_hallmark.learned import priors as PR
    genes = genes or DISC.PANEL
    rows = []
    for g in genes:
        try:
            Y, Xo, C, w = LD.assemble_gene(g, w_prior_source="ledger", orphans=True)
            _, Xhe, _, _ = LD.assemble_gene(g, w_prior_source="ledger")
            he = set(Xhe.columns)
            orphans = [a for a in Xo.columns if a not in he]
            if not orphans:
                continue
            Cm = C.to_numpy(float)
            yr = _resid(Y.to_numpy(float), Cm)
            Xr = _resid(Xo.to_numpy(float), Cm)
            sd = Xr.std(0)
            Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
            wv = np.nan_to_num(w.reindex(Xo.columns).to_numpy(float), nan=0.0)
            pi = np.where(sd > 1e-9, PR.inclusion_prior(wv), 0.0)             # orphans → low π (near-0 ledger weight)
            _, _, pip = SS._gibbs_posterior(Xz, -yr, np.clip(pi, 0, 1), n_iter=1500, burn=500, seed=0)
            pip = pd.Series(pip, index=Xo.columns)
            ss_orph = set(a for a in orphans if pip[a] > pip_thr)             # spike-slab discoveries
            ref = DISC.discover_gene(g, min_partial=min_partial, deconv_check=False)  # partial-Spearman reference
            ps_orph = set(ref["arm"]) if len(ref) else set()
            rows.append({"gene": g, "n_orphan": len(orphans),
                         "ss_pip": len(ss_orph), "partial": len(ps_orph),
                         "overlap": len(ss_orph & ps_orph),
                         "ss_only": len(ss_orph - ps_orph), "partial_only": len(ps_orph - ss_orph),
                         "top_ss": ", ".join(sorted(a.replace("hsa-", "") for a in ss_orph)[:4])})
        except Exception as e:
            rows.append({"gene": g, "error": repr(e)[:50]})
    df = pd.DataFrame(rows)
    with pd.option_context("display.width", 170, "display.max_colwidth", 44):
        print(df.to_string(index=False))
    d = df[df.get("overlap").notna()] if "overlap" in df else df.iloc[0:0]
    if len(d):
        tp = d["overlap"].sum(); ss = d["ss_pip"].sum(); ps = d["partial"].sum()
        print(f"\n=== spike-slab-inclusion discovery vs partial-Spearman ({len(d)} genes) ===")
        print(f"  spike-slab PIP>{pip_thr}: {ss} orphans | partial-Spearman: {ps} | overlap {tp} "
              f"(precision {tp/ss:.2f}, recall {tp/ps:.2f})" if ss and ps else "  (sparse)")
        print(f"  ⇒ the INCLUSION mode {'recovers' if tp else 'does NOT recover'} the partial-Spearman discoveries; "
              f"ss-only {d['ss_only'].sum()} (new nominations), partial-only {d['partial_only'].sum()} (missed)")
    return df


def _ss_orphan_pip(g, *, deconv=False, permute=None, prior_mode="learned", orphan_source="targetscan"):
    """Per gene: orphan PIP from a spike-slab on [HE∪orphan] arms. `prior_mode`: 'learned' (Beta–Bernoulli
    base rate, no hand-set prior — POLICY DEFAULT), 'uniform' (π=0.5, high-recall), 'evidence' (skeptical,
    deprecated). deconv=True puts CIBERSORTx fractions in C; permute shuffles Y (FDR null)."""
    from mirna_hallmark.learned import priors as PR
    Y, Xo, C, w = LD.assemble_gene(g, w_prior_source="ledger", orphans=True, deconv=deconv, orphan_source=orphan_source)
    _, Xhe, _, _ = LD.assemble_gene(g, w_prior_source="ledger", deconv=deconv)
    orphans = [a for a in Xo.columns if a not in set(Xhe.columns)]
    if not orphans:
        return None, []
    if permute is not None:
        Y = pd.Series(np.random.default_rng(permute).permutation(Y.to_numpy()), index=Y.index)
    Cm = C.to_numpy(float); yr = _resid(Y.to_numpy(float), Cm); Xr = _resid(Xo.to_numpy(float), Cm)
    sd = Xr.std(0); Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    wv = np.nan_to_num(w.reindex(Xo.columns).to_numpy(float), nan=0.0)
    mask = sd > 1e-9
    learn = False
    if prior_mode == "densez":                                    # dense ridge (π≡1) → return |z| as the score
        m, s, _ = SS._gibbs_posterior(Xz, -yr, mask.astype(float), n_iter=1200, burn=400, seed=0)
        return pd.Series(np.abs(m / (s + 1e-9)), index=Xo.columns), orphans   # score = |z|, threshold 2 in caller
    if prior_mode == "evidence":
        pi = np.where(mask, PR.inclusion_prior(wv), 0.0)
    elif prior_mode == "uniform":
        pi = np.where(mask, 0.5, 0.0)
    else:                                                          # 'learned' Beta–Bernoulli base rate
        pi = mask.astype(float); learn = True
    _, _, pip = SS._gibbs_posterior(Xz, -yr, np.clip(pi, 0, 1), n_iter=1200, burn=400, seed=0, learn_pi0=learn)
    return pd.Series(pip, index=Xo.columns), orphans


def discovery_full(*, out: str | None = None, prior_mode: str = "learned", orphan_source: str = "targetscan",
                   pip_thr: float = 0.5, min_scanmir: float = 0.0, limit: int | None = None,
                   progress: int = 50, fdr_permute: int = 1) -> pd.DataFrame:
    """GENOME-WIDE gated spike-slab discovery (POLICY: `prior_mode='learned'` base rate). PIP>thr + scanMiR
    duplex, then a **composition (deconv) gate** + a **permutation FDR**. Cross-tab vs `discoveries.tsv`."""
    from pathlib import Path
    from mirna_hallmark.learned import kd as KD
    from mirna_hallmark.learned.evidence import ledger as LG
    if out is None:
        sfx = "" if orphan_source == "targetscan" else f"_{orphan_source}"
        out = f"mirna_hallmark/output/learned/discovery_bayes_{prior_mode}{sfx}.tsv"
    pip_thr = 2.0 if prior_mode == "densez" else pip_thr           # densez score is |z| (threshold 2), else PIP (0.5)
    genes = sorted(set(LG.pooled_he_edges()["gene"].dropna().astype(str)))
    if limit:
        genes = genes[:limit]
    aff = KD.affinity()
    rows, null_hits = [], 0
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[discovery_bayes] {i}/{len(genes)} ({len(rows)} hits)", flush=True)
        affg = aff.loc[aff["gene"] == g].set_index("arm")["repression"]
        try:                                                                   # FDR null: same scan on shuffled Y
            pipn, orphn = _ss_orphan_pip(g, permute=1000 + i, prior_mode=prior_mode, orphan_source=orphan_source)
            if pipn is not None:
                null_hits += sum(1 for a in orphn if pipn[a] > pip_thr and affg.get(a, 0) < min_scanmir)
        except Exception:
            pass
        try:
            pip, orphans = _ss_orphan_pip(g, prior_mode=prior_mode, orphan_source=orphan_source)
            if pip is None:
                continue
            cand = [a for a in orphans if pip[a] > pip_thr and (affg.get(a, np.nan) == affg.get(a, np.nan))
                    and affg.get(a, np.nan) < min_scanmir]
            if not cand:
                continue
            pipd, _ = _ss_orphan_pip(g, deconv=True, prior_mode=prior_mode, orphan_source=orphan_source)                           # composition gate
            for a in cand:
                pd_ = float(pipd[a]) if (pipd is not None and a in pipd.index) else np.nan
                rows.append({"gene": g, "arm": a, "pip": round(float(pip[a]), 2),
                             "scanmir_rep": round(float(affg.get(a, np.nan)), 2),
                             "pip_deconv": round(pd_, 2) if pd_ == pd_ else np.nan,
                             "robust": bool(pd_ == pd_ and pd_ > pip_thr)})     # survives composition control
        except Exception:
            pass
    df = pd.DataFrame(rows)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    rob = df[df["robust"]] if "robust" in df else df
    fdr = null_hits / max(len(df), 1)
    print(f"\n=== GENOME-WIDE GATED spike-slab discovery: {len(df)} raw / {len(rob)} composition-robust / "
          f"{df['gene'].nunique()} genes ===")
    print(f"  permutation FDR: {null_hits} null hits vs {len(df)} real ⇒ empirical FDR ≈ {fdr:.3f}")
    try:
        ps = pd.read_csv("mirna_hallmark/output/learned/discoveries.tsv", sep="\t")
        ps = ps[ps["robust"]]
        se = set(zip(rob["gene"], rob["arm"])); pe = set(zip(ps["gene"], ps["arm"]))
        ov = len(se & pe)
        print(f"  robust-vs-robust: overlap {ov} | ss-only {len(se-pe)} | partial-only {len(pe-se)} "
              f"(recall {ov/max(len(pe),1):.2f}, precision {ov/max(len(se),1):.2f})")
    except Exception as e:
        print(f"  (cross-tab skipped: {e})")
    return df


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cptac", action="store_true")
    ap.add_argument("--shuffle", action="store_true")
    ap.add_argument("--discovery", action="store_true")
    ap.add_argument("--discovery-full", action="store_true")
    ap.add_argument("--prior-mode", default="learned", choices=["learned","uniform","evidence","densez"])
    ap.add_argument("--orphan-source", default="targetscan", choices=["targetscan","kd"])
    ap.add_argument("--cohort", default="prospective")
    a = ap.parse_args()
    if a.cptac:
        cptac_parity(cohort=a.cohort)
    elif a.shuffle:
        shuffle_parity()
    elif a.discovery_full:
        discovery_full(prior_mode=a.prior_mode, orphan_source=a.orphan_source)
    elif a.discovery:
        discovery_parity()
    else:
        cptac_parity()
