"""TASK 4 — PINPOINT the composition confound to a SPECIFIC Wu-major lineage, per EDGE and per GENE.

Estimator = the canonical §6b one, unchanged: `attribution_eb._prep` (C-residualised, z-scored,
variance-floored design) -> `spike_slab._gibbs_posterior(pi=1)` (dense, n_iter=2000, burn=700, seed=0).
Verified to reproduce `readouts_edges.tsv` beta exactly (PTEN, max|diff| ~1e-6).

Per gene, per level ({family, arm}), 18 posteriors:
  core                  C = core (CPE + target_cn + mal_prolif)                 -> beta_core
  full                  C = core + all 8 non-malignant Wu lineages              -> beta_deconv (== readouts)
  +L (x8)               C = core + ONE lineage                                  -> beta_L
  +shuf(L) (x8)         C = core + ONE lineage PERMUTED ACROSS PATIENTS         -> NULL (df-matched placebo)

lineage_loss_L  = 1 - beta_L / beta_core      (fraction of the coupling that lineage L ALONE explains)
full_loss       = 1 - beta_full / beta_core   (= 1 - retention)
pinpointed      = argmax_L lineage_loss_L, but ONLY if it beats the max over the 8 SHUFFLED lineages
                  (the null takes the SAME max-over-8, so the argmax's own selection bias is controlled)
                  AND explains >= 50% of full_loss. Otherwise -> mixed_no_single_lineage.

ORIENTATION (MH-112): sign(spearman(miRNA, L)) x sign(spearman(target, L)).
  OPPOSITE -> bulk-mixing arithmetic.  SAME -> candidate compartment-INTRINSIC edge.

Outputs (mirna_hallmark/output/learned/):
  composition_pinpoint_edges_family.tsv · composition_pinpoint_edges_arm.tsv   (per-(gene,unit) CARD)
  composition_pinpoint_genes.tsv                                               (per-GENE roll-up)

CLI:  .venv/bin/python3 -m mirna_hallmark.analyses.cptac.composition_pinpoint_lineage [--limit N] [--workers 8]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

OUT = Path("mirna_hallmark/output/learned")
N_ITER, BURN, SEED = 2000, 700, 0
LINEAGES = ["CAFs", "T-cells", "Myeloid", "B-cells", "Endothelial", "PVL",
            "Normal Epithelial", "Plasmablasts"]
MIN_BETA = 1e-9


def _fit(Y, Xf, C):
    from mirna_hallmark.learned import attribution_eb as AE, spike_slab as SS
    yr, Xz, cols = AE._prep(Y, Xf, C)
    b, sd, _ = SS._gibbs_posterior(Xz, yr, np.ones(Xz.shape[1]), n_iter=N_ITER, burn=BURN, seed=SEED)
    return pd.Series(b, index=cols), pd.Series(sd, index=cols)


def gene_pinpoint(gene: str):
    """-> (edges_family, edges_arm) frames with the per-lineage betas / losses / nulls / orientation."""
    from mirna_hallmark.learned import data as LD, families as FAM

    Y, X, C0, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=False)
    dv = LD._deconv()
    if dv is None:
        return None
    dvj = dv.reindex(C0.index)
    dvj = dvj.fillna(dvj.median())                      # EXACTLY as assemble_gene(deconv=True) does
    dvj = dvj[[c for c in LINEAGES if c in dvj.columns]]

    rng = np.random.default_rng(abs(hash(gene)) % (2**31))
    shuf = pd.DataFrame({c: rng.permutation(dvj[c].to_numpy()) for c in dvj.columns}, index=dvj.index)

    # target's own compartment loading (raw spearman, MH-112 convention)
    tgt_rho = {L: spearmanr(Y.to_numpy(float), dvj[L].to_numpy(float)).statistic for L in dvj.columns}

    out = {}
    fam_map = FAM.family_of(pd.Index(X.columns))
    for level in ("family", "arm"):
        if level == "family":
            Xf, wf, members = FAM.collapse_by_family(X, w, fam_map)
        else:
            Xf, members = X, {a: [a] for a in X.columns}
        if Xf.shape[1] < 1:
            continue

        b_core, sd_core = _fit(Y, Xf, C0)
        b_full, sd_full = _fit(Y, Xf, C0.join(dvj))

        rows = pd.DataFrame({"gene": gene, "unit": b_core.index, "level": level,
                             "n": len(Y), "p_unit": Xf.shape[1],
                             "beta_core": b_core.to_numpy(), "beta_sd_core": sd_core.to_numpy(),
                             "beta_full": b_full.reindex(b_core.index).to_numpy()})
        rows["z_core"] = np.where(rows.beta_sd_core > 1e-12, rows.beta_core / rows.beta_sd_core, 0.0)
        rows["identified"] = rows.z_core.abs() > 2.0
        with np.errstate(divide="ignore", invalid="ignore"):
            rows["retention"] = np.where(rows.beta_core.abs() > MIN_BETA,
                                         rows.beta_full / rows.beta_core, np.nan)
        rows["full_loss"] = 1.0 - rows.retention
        rows["composition_class"] = pd.cut(rows.retention, [-np.inf, 0.4, 0.7, np.inf],
                                           labels=["composition_explained", "partial", "cell_intrinsic"]
                                           ).astype(object)

        for L in dvj.columns:
            bL, _ = _fit(Y, Xf, C0.join(dvj[[L]]))
            bN, _ = _fit(Y, Xf, C0.join(shuf[[L]].rename(columns={L: L})))
            bL = bL.reindex(b_core.index).to_numpy()
            bN = bN.reindex(b_core.index).to_numpy()
            with np.errstate(divide="ignore", invalid="ignore"):
                rows[f"loss_{L}"] = np.where(rows.beta_core.abs() > MIN_BETA, 1.0 - bL / rows.beta_core, np.nan)
                rows[f"null_{L}"] = np.where(rows.beta_core.abs() > MIN_BETA, 1.0 - bN / rows.beta_core, np.nan)
            rows[f"beta_{L}"] = bL
            # miRNA-side compartment loading of THIS unit (family-collapsed predictor or arm)
            rows[f"mirrho_{L}"] = [spearmanr(Xf[u].to_numpy(float), dvj[L].to_numpy(float)).statistic
                                   for u in b_core.index]
            rows[f"tgtrho_{L}"] = tgt_rho[L]

        rows["arms"] = [";".join(sorted(members.get(u, []))) for u in b_core.index]
        rows["n_arms"] = [len(members.get(u, [])) for u in b_core.index]
        if level == "arm":
            rows["seed_family"] = rows.unit.map(fam_map)
        out[level] = rows
    return out


def _one(g):
    try:
        return gene_pinpoint(g)
    except Exception as e:                                        # noqa: BLE001
        sys.stderr.write(f"[skip] {g}: {type(e).__name__}: {e}\n")
        return None


def classify(d: pd.DataFrame) -> pd.DataFrame:
    """argmax-over-8 pinpoint + max-over-8 SHUFFLED null + orientation -> final tag."""
    lc = [f"loss_{L}" for L in LINEAGES if f"loss_{L}" in d.columns]
    nc = [f"null_{L}" for L in LINEAGES if f"null_{L}" in d.columns]
    Lnames = [c[5:] for c in lc]
    Lo = d[lc].to_numpy(float)
    Nu = d[nc].to_numpy(float)

    with np.errstate(invalid="ignore"):
        best = np.nanargmax(np.where(np.isfinite(Lo), Lo, -np.inf), axis=1)
    d["pinpointed_lineage"] = [Lnames[i] for i in best]
    d["lineage_loss_fraction"] = Lo[np.arange(len(d)), best]
    srt = np.sort(np.where(np.isfinite(Lo), Lo, -np.inf), axis=1)
    d["second_loss"] = srt[:, -2]
    d["loss_margin"] = d.lineage_loss_fraction - d.second_loss          # argmax-tie guard (rule d)
    d["max_null_loss"] = np.nanmax(np.where(np.isfinite(Nu), Nu, -np.inf), axis=1)
    d["mean_null_loss"] = np.nanmean(Nu, axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        d["explained_share"] = np.where(d.full_loss.abs() > 1e-6,
                                        d.lineage_loss_fraction / d.full_loss, np.nan)

    # ⚠ itertuples() MANGLES column names with spaces/hyphens ('T-cells' -> _7): index by POSITION instead.
    M = d[[f"mirrho_{L}" for L in Lnames]].to_numpy(float)
    T = d[[f"tgtrho_{L}" for L in Lnames]].to_numpy(float)
    d["mir_rho"] = M[np.arange(len(d)), best]
    d["tgt_rho"] = T[np.arange(len(d)), best]
    d["orientation"] = np.where(np.sign(d.mir_rho) == np.sign(d.tgt_rho), "same", "opposite")

    beats_null = d.lineage_loss_fraction > d.max_null_loss
    dominant = d.explained_share >= 0.5
    d["single_lineage"] = beats_null & dominant

    # ⚠ MEASURED (this run, 5117 family couplings): at the PINPOINTED lineage the SAME/OPPOSITE tag is
    # 98% DETERMINED BY sign(beta_core) — P(same | beta<0) = 0.981, P(same | beta>0) = 0.159. That is
    # MECHANICAL, not biological: the lineage that most REDUCES beta is by construction the one whose
    # confounding path rho(miR,L)·rho(tgt,L) carries the SAME SIGN as the coupling it is eating, and the
    # design's target is entered as -resid(Y) (so beta>0 == repression). ⇒ an ANTI-repressive edge
    # (beta<0) is AUTOMATICALLY "same-compartment" and its orientation carries NO information.
    # Gate the compartment-intrinsic call on beta_core > 0 (a real repressive coupling); anti-repressive
    # couplings get their own class instead of being silently promoted to a stromal-network "discovery".
    conf = d.composition_class.isin(["composition_explained", "partial"])
    tag = np.where(~conf, "cell_intrinsic",
          np.where(d.beta_core < 0, "anti_repressive_unresolved",
          np.where(~d.single_lineage, "mixed_no_single_lineage",
          np.where(d.orientation.eq("same"),
                   "compartment_intrinsic:" + d.pinpointed_lineage,
                   "artifact_opposite_compartment:" + d.pinpointed_lineage))))
    d["final_tag"] = tag
    d["orientation_informative"] = d.beta_core > 0          # else orientation == sign(beta), tautological
    # ...and a sign() of a ~0 loading is NOISE. n=1041 ⇒ |rho|>0.061 is p<0.05; require BOTH sides to carry
    # a non-trivial loading before an "orientation" is allowed to mean anything. MEASURED: 0 of the 9 (family)
    # / 11 (arm) compartment-intrinsic candidates survive |rho|>0.10 on both sides.
    d["orientation_loading_ok"] = (d.mir_rho.abs() > 0.10) & (d.tgt_rho.abs() > 0.10)
    d["retention_reliable"] = d.identified                              # MH-119: gate the ratio readout
    return d


def rollup(d: pd.DataFrame) -> pd.DataFrame:
    """Per-GENE: Shapley of a LINEAR aggregate is additive => gene_retention = sum(beta_deconv)/sum(beta).
    ⚠ MH-119: the denominator can CANCEL. `denom_coherence = |sum b| / sum|b|` gates the ratio."""
    rows = []
    for (g, lev), s in d.groupby(["gene", "level"]):
        tot = s.beta_core.sum()
        coh = abs(tot) / (s.beta_core.abs().sum() + 1e-12)
        r = {"gene": g, "level": lev, "n_units": len(s),
             "sum_beta_core": tot, "sum_beta_full": s.beta_full.sum(),
             "denom_coherence": coh,
             "gene_retention": s.beta_full.sum() / tot if abs(tot) > 1e-9 else np.nan}
        r["gene_reliable"] = bool(abs(tot) > 1e-6 and coh > 0.5 and s.identified.any())
        losses = {}
        for L in LINEAGES:
            bl = f"beta_{L}"
            if bl in s.columns:
                losses[L] = 1.0 - s[bl].sum() / tot if abs(tot) > 1e-9 else np.nan
                r[f"gene_loss_{L}"] = losses[L]
        nulls = []
        for L in LINEAGES:
            nc = f"null_{L}"
            if nc in s.columns:                                  # gene-level null: sum of the per-unit null betas
                bnull = (1.0 - s[nc]) * s.beta_core              # reconstruct beta_null from the ratio
                nulls.append(1.0 - bnull.sum() / tot if abs(tot) > 1e-9 else np.nan)
        r["gene_max_null_loss"] = np.nanmax(nulls) if nulls else np.nan
        if losses:
            ser = pd.Series(losses)
            r["gene_pinpointed_lineage"] = ser.idxmax()
            r["gene_lineage_loss"] = ser.max()
            r["gene_loss_margin"] = ser.max() - ser.sort_values().iloc[-2]
            r["gene_full_loss"] = 1.0 - r["gene_retention"] if pd.notna(r["gene_retention"]) else np.nan
            r["gene_explained_share"] = (r["gene_lineage_loss"] / r["gene_full_loss"]
                                         if pd.notna(r["gene_full_loss"]) and abs(r["gene_full_loss"]) > 1e-6
                                         else np.nan)
            r["gene_single_lineage"] = bool(r["gene_lineage_loss"] > r["gene_max_null_loss"]
                                            and pd.notna(r["gene_explained_share"])
                                            and r["gene_explained_share"] >= 0.5)
        r["gene_composition_class"] = ("cell_intrinsic" if pd.notna(r["gene_retention"]) and r["gene_retention"] >= 0.7
                                       else "partial" if pd.notna(r["gene_retention"]) and r["gene_retention"] >= 0.4
                                       else "composition_explained" if pd.notna(r["gene_retention"]) else "na")
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--genes", type=str, default=None)
    a = ap.parse_args()

    from mirna_hallmark.learned import data as LD
    r = pd.read_csv(OUT / "readouts_edges.tsv", sep="\t")
    genes = sorted(r.gene.unique())
    if a.genes:
        genes = a.genes.split(",")
    if a.limit:
        genes = genes[:a.limit]
    print(f"[pinpoint] {len(genes)} genes · workers={a.workers} · 18 posteriors/gene/level", flush=True)

    LD._load()                                                   # warm the cache BEFORE forking
    LD._deconv()
    t0 = time.time()
    fam, arm = [], []
    if a.workers > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.workers) as pool:
            for i, res in enumerate(pool.imap_unordered(_one, genes, chunksize=4), 1):
                if res:
                    if "family" in res:
                        fam.append(res["family"])
                    if "arm" in res:
                        arm.append(res["arm"])
                if i % 100 == 0:
                    print(f"  {i}/{len(genes)}  {time.time()-t0:.0f}s", flush=True)
    else:
        for g in genes:
            res = _one(g)
            if res:
                fam.append(res.get("family"))
                arm.append(res.get("arm"))

    done = {}
    for name, chunks in (("family", fam), ("arm", arm)):
        d = classify(pd.concat([c for c in chunks if c is not None], ignore_index=True))
        d.to_csv(OUT / f"composition_pinpoint_edges_{name}.tsv", sep="\t", index=False)
        done[name] = d
        print(f"[write] composition_pinpoint_edges_{name}.tsv  {d.shape}")
    gd = rollup(pd.concat(done.values(), ignore_index=True))
    gd.to_csv(OUT / "composition_pinpoint_genes.tsv", sep="\t", index=False)
    print(f"[write] composition_pinpoint_genes.tsv  {gd.shape}   {time.time()-t0:.0f}s total")


if __name__ == "__main__":
    main()
