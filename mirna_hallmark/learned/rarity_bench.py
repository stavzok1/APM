"""Converge the specificity/rarity term — WHERE does seed-rarity enter the final Bayes model, and at what gain?

Two loci, benched in ONE pass per gene (shared prep) on both readouts of the one posterior:

  * INCLUSION π   (spike-and-slab / SELECTION+DISCOVERY readout)  — logit(π) += g_r·z_rarity.
                  Inert in the dense ridge (π≡1), so it can ONLY move the SS-selected aggregate.
  * SLAB SCALE    (dense learned-τ² ridge = model 'A', AND the SS slab) — scale_m = 1 + μ·τ + rarity-term.
                  This is the channel that lets specificity act on the MAIN (dense) model. Two flavours:
                    gated (·τ)   : trust rarity only where there is evidence   → 1 + μ·τ + s_r·rarity·τ
                    free (stand.): pure-sequence, no evidence needed           → 1 + μ·τ + s_r·rarity

Readouts per (gene, variant):
  * ρ_dense : 5-fold OOF coupling of the DENSE ridge (π≡1, half-normal slab, learned ν²) with the variant's slab.
  * ρ_ss    : 5-fold OOF coupling of the spike-and-slab SELECTED aggregate (evidence±rarity π, variant slab).
  * top arm + its seed-rarity (dense: max |M| attribution; ss: max PIP selection) — does rarity re-nominate the
    collinear-soup specialist? + n_selected + selection stability (mean Jaccard across folds).

Panel: collinear SOUPS {PTEN, CCND1, VEGFA} (where rarity should bite) vs CLEAN {ESR1, CDH1} (where it must be
inert — the harm test). Baseline = no rarity anywhere; π-rarity affects only ρ_ss (ρ_dense≡baseline by design).

CLI: .venv/bin/python3 -m mirna_hallmark.learned.rarity_bench            # panel, default sweep
     .venv/bin/python3 -m mirna_hallmark.learned.rarity_bench --gains 0 0.8 1.6
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
from mirna_hallmark.learned import priors as PR
from mirna_hallmark.learned import seed_rarity as SR
from mirna_hallmark.learned import spike_slab as SS

SOUP = ["PTEN", "CCND1", "VEGFA"]
CLEAN = ["ESR1", "CDH1"]
MU_GAIN = 1.0
SLAB_FLOOR = 0.25


def _resid(V, Cmat):
    return V - LinearRegression().fit(Cmat, V).predict(Cmat)


def _prep(gene: str):
    """Assemble family predictors once; per-edge evidence w, biochem μ, confidence τ, seed rarity."""
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
    cols = X.columns
    wv = w.reindex(cols).fillna(0.0).to_numpy(dtype=float)
    mu = PR._biochem_magnitude(gene, cols, by="family").to_numpy(dtype=float)
    tau = PR._confidence(wv)
    rar = SR.rarity(cols).reindex(cols).to_numpy(dtype=float)                 # NaN where uncovered
    return {"Y": Y, "X": X, "C": C, "cols": cols, "w": wv, "mu": mu, "tau": tau, "rarity": rar}


def _slab(p, *, locus: str, gain: float) -> np.ndarray:
    """Per-edge slab-scale multiplier. base = 1 + μ_gain·μ·τ (the shipped biochem×evidence widening).
    locus 'slabG' adds s_r·rarity·τ (evidence-gated); 'slabF' adds s_r·rarity (standalone). rarity NaN → 0."""
    base = 1.0 + MU_GAIN * p["mu"] * p["tau"]
    r = np.where(np.isnan(p["rarity"]), 0.0, p["rarity"])
    if locus == "slabG":
        base = base + gain * r * p["tau"]
    elif locus == "slabF":
        base = base + gain * r
    return np.clip(base, SLAB_FLOOR, None)


def _pi(p, *, rarity_on: bool, gain: float) -> np.ndarray:
    """Evidence-graded inclusion π, optionally with the rarity ORDERING channel at strength `gain`."""
    rr = p["rarity"] if rarity_on else None
    return PR.inclusion_prior(p["w"], rarity=rr, rarity_gain=gain)


def _oof(p, *, dense: bool, pi_full, slab, folds=5, seed=0, n_iter=900, burn=300):
    """5-fold OOF aggregate coupling for one variant. dense=True forces π≡1 (ridge, model A);
    else uses `pi_full` (spike-and-slab selection). Returns (rho, top_arm, top_rarity, n_active, stability)."""
    Y = p["Y"].to_numpy(dtype=float); X = p["X"].to_numpy(dtype=float); Cmat = p["C"].to_numpy(dtype=float)
    cols = p["cols"]; n, q = X.shape
    oof = np.full(n, np.nan)
    sets = []
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    for fi, (tr, te) in enumerate(kf.split(X)):
        Ctr = Cmat[tr]
        r = -_resid(Y[tr], Ctr)
        Xr = _resid(X[tr], Ctr)
        sd = Xr.std(0)
        Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
        pv = np.ones(q) if dense else np.where(sd > 1e-9, pi_full, 0.0)
        beta, pip = SS._gibbs_ss(Xz, r, pv, slab_scale=slab, n_iter=n_iter, burn=burn, seed=seed + fi)
        M = np.where(sd > 1e-9, beta / (sd + 1e-9), 0.0)
        oof[te] = X[te] @ M
        sel = pip if not dense else (beta > 0).astype(float)
        sets.append(frozenset(cols[np.asarray(sel) > 0.5]))
    yr = _resid(Y, Cmat)
    rho = spearmanr(_resid(oof, Cmat), yr).correlation
    # full-data fit → identity of the top edge
    r = -_resid(Y, Cmat); Xr = _resid(X, Cmat); sd = Xr.std(0)
    Xz = np.where(sd > 1e-9, (Xr - Xr.mean(0)) / (sd + 1e-9), 0.0)
    pv = np.ones(q) if dense else np.where(sd > 1e-9, pi_full, 0.0)
    beta, pip = SS._gibbs_ss(Xz, r, pv, slab_scale=slab, n_iter=n_iter, burn=burn, seed=seed)
    M = np.where(sd > 1e-9, beta / (sd + 1e-9), 0.0)
    score = M if dense else pip * M                            # dense: attribution magnitude; ss: PIP-weighted
    top_i = int(np.argmax(score)) if np.any(score > 0) else -1
    top_arm = cols[top_i] if top_i >= 0 else "-"
    top_rar = round(float(p["rarity"][top_i]), 2) if top_i >= 0 and not np.isnan(p["rarity"][top_i]) else np.nan
    n_active = int((pip > 0.5).sum()) if not dense else int((M > 0).sum())
    # SPECIALIST metrics (the thing rarity is FOR): PIP of the rarest covered edge, and the PIP↔rarity tilt.
    rar = p["rarity"]; cov = ~np.isnan(rar)
    sel = pip if not dense else M / (np.max(M) + 1e-12)         # dense: normalised magnitude as pseudo-PIP
    if cov.sum() >= 1:
        spec_i = int(np.nanargmax(np.where(cov, rar, -np.inf)))
        spec_pip = round(float(sel[spec_i]), 3)
        spec_arm = str(cols[spec_i])
    else:
        spec_pip, spec_arm = np.nan, "-"
    tilt = (round(float(spearmanr(sel[cov], rar[cov]).correlation), 3)
            if cov.sum() >= 3 and np.ptp(sel[cov]) > 0 else np.nan)   # >0 ⇒ selection tilts to rare seeds
    return (round(float(rho), 3), str(top_arm), top_rar, n_active, round(_mean_jacc(sets), 2),
            spec_arm, spec_pip, tilt)


def _mean_jacc(sets):
    js = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            js.append(len(sets[i] & sets[j]) / len(u) if u else 1.0)
    return float(np.mean(js)) if js else float("nan")


def variants(gains):
    """(name, locus, gain, rarity_in_pi) tuples. base once; π-rarity and slab-rarity(gated/free) per gain."""
    v = [("base", "base", 0.0, False)]
    for g in gains:
        if g == 0:
            continue
        v.append((f"pi@{g}", "base", g, True))
        v.append((f"slabG@{g}", "slabG", g, False))
        v.append((f"slabF@{g}", "slabF", g, False))
    return v


def bench(genes=None, *, gains=(0.4, 0.8, 1.2, 1.6), folds=5, out=None) -> pd.DataFrame:
    genes = genes or (SOUP + CLEAN)
    rows = []
    for g in genes:
        klass = "soup" if g in SOUP else "clean"
        try:
            p = _prep(g)
        except Exception as e:
            rows.append({"gene": g, "error": repr(e)[:80]}); continue
        for name, locus, gain, rar_pi in variants(gains):
            slab = _slab(p, locus=locus, gain=gain)
            pi_full = _pi(p, rarity_on=rar_pi, gain=gain)
            # DENSE (model A): π ignored → only slab-locus variants differ from base
            if locus != "base" or name == "base":
                rd, ad, rrd, nd, sd_, sa, sp, tl = _oof(p, dense=True, pi_full=None, slab=slab, folds=folds)
                rows.append({"gene": g, "class": klass, "readout": "dense", "variant": name,
                             "rho": rd, "top": ad, "top_rar": rrd, "n_active": nd, "stab": sd_,
                             "spec": sa, "spec_pip": sp, "tilt": tl})
            # SS (selection): all variants differ
            rs, as_, rrs, ns, ss_, sa, sp, tl = _oof(p, dense=False, pi_full=pi_full, slab=slab, folds=folds)
            rows.append({"gene": g, "class": klass, "readout": "ss", "variant": name,
                         "rho": rs, "top": as_, "top_rar": rrs, "n_active": ns, "stab": ss_,
                         "spec": sa, "spec_pip": sp, "tilt": tl})
        print(f"[rarity_bench] {g} done ({klass})", flush=True)
    df = pd.DataFrame(rows)
    if out:
        from pathlib import Path
        Path(out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, sep="\t", index=False)
        print(f"wrote {out}")
    _summ(df)
    return df


def _summ(df: pd.DataFrame):
    if "rho" not in df:
        return
    with pd.option_context("display.width", 200, "display.max_rows", 300, "display.max_colwidth", 30):
        print("\n=== per (gene, readout, variant) ===")
        print(df.to_string(index=False))
    # Δρ vs base, by readout × class × variant (more-negative rho = better coupling ⇒ Δ<0 = rarity helps)
    print("\n=== Δρ vs base (negative = rarity improves coupling) ===")
    for readout in ["dense", "ss"]:
        sub = df[df["readout"] == readout]
        base = sub[sub["variant"] == "base"].set_index("gene")["rho"]
        for klass in ["soup", "clean"]:
            ss = sub[(sub["class"] == klass) & (sub["variant"] != "base")]
            if ss.empty:
                continue
            agg = (ss.assign(d=ss.apply(lambda r: r["rho"] - base.get(r["gene"], np.nan), axis=1))
                     .groupby("variant")["d"].mean())
            print(f"  [{readout} · {klass}]  " + "  ".join(f"{v}:{d:+.3f}" for v, d in agg.items()))
    # SPECIALIST nomination (ss only): Δ spec_pip vs base + mean PIP↔rarity tilt — the discovery payoff of rarity.
    print("\n=== specialist nomination (ss): Δ spec_pip vs base | tilt=Spearman(PIP,rarity) ===")
    sub = df[df["readout"] == "ss"]
    base_sp = sub[sub["variant"] == "base"].set_index("gene")["spec_pip"]
    for klass in ["soup", "clean"]:
        ss = sub[(sub["class"] == klass) & (sub["variant"] != "base")]
        if ss.empty:
            continue
        dsp = (ss.assign(d=ss.apply(lambda r: r["spec_pip"] - base_sp.get(r["gene"], np.nan), axis=1))
                 .groupby("variant")["d"].mean())
        tl = ss.groupby("variant")["tilt"].mean()
        print(f"  [{klass}] Δspec_pip  " + "  ".join(f"{v}:{d:+.3f}" for v, d in dsp.items()))
        print(f"  [{klass}] tilt       " + "  ".join(f"{v}:{t:+.3f}" for v, t in tl.items()))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--gains", nargs="*", type=float, default=[0.4, 0.8, 1.2, 1.6])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--out", default="mirna_hallmark/output/learned/rarity_bench.tsv")
    a = ap.parse_args()
    bench(a.genes, gains=a.gains, folds=a.folds, out=a.out)
