"""§Y — DOES RISC/AGO AVAILABILITY MODIFY THE DE-CONFOUNDED SLOPE β? (MH-164)

    .venv/bin/python3 -m mirna_hallmark.learned.ago_interaction            # observed + shuffled-ago null

THE TEST (board §Y / METHODS §4; never run before — the pressure arc applied the AGO gate as a
multiplier and reported "gated AND ungated", calling it "a sensitivity layer, NOT a causal model").
Per HE gene g (~1000 shared TCGA patients):
    r   = -resid(Y | C)                     repression readout (AE._prep; C = core CPE/target_cn/mal_prolif)
    M   = Xhe_z @ beta_card                  canonical de-confounded budget (== discovery.py he_agg, verbatim)
    ago = ago_capacity_z                     per-sample RISC capacity (co-limited AGO1-4 ∧ TNRC6, ago_gate.py)
    fit  r ~ M + ago + prolif + M:ago + M:prolif   → test the M:ago coefficient, aggregate across genes.

RIGOR (why the naive answer is wrong):
  * ago_perp = ago ⊥ confounder-block (composition+prolif+purity; R²(ago~C)=0.395 ⇒ 40% of ago IS the
    confounders) ⇒ the M:ago_perp interaction is ago-SPECIFIC, not a proliferation-modifier leak. ⚠ but the
    OVER-CONTROL check shows this barely matters: raw-ago +0.00671 vs ago⊥C +0.00656 (2% shift) — the null is
    NOT manufactured by residualisation (rigor-auditor cleared this).
  * The per-gene t's are NOT independent (shared patients + shared ago) ⇒ a naive Stouffer Z=+7.8 is INVALID
    pseudo-replication. The ONLY honest arbiter is the SHUFFLED-AGO null: permute ago across patients (holds
    every gene's M and the full gene-gene correlation fixed, breaks only the ago–patient link).
  * M:prolif (Stouffer −19) = a working POSITIVE CONTROL for the interaction machinery (a real slope-modifier
    IS detected). ⚠ M:purity / M:composition interactions are NOT carried (conservative-safe for a null).

VERDICT (MH-164, tag P): the M:ago_perp interaction is UNDETECTED at n≈1000, with the observed trend in the
PREDICTED (positive) direction — perm p≈0.10–0.12, ~40% power, MDE ≈ 1.7× the observed effect. NOT "AGO does
not modify β"; it retires the gate to "biological rationale, not demonstrated coupling-modification", no stronger.
"""
from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd

OUT = None                                                  # set in run() to avoid import-time config load


def run(n_null: int = 300, seed: int = 0):
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    warnings.filterwarnings("ignore")
    from scipy import stats
    from mirna_hallmark import config as CFG
    from mirna_hallmark import ago_gate as AG
    from mirna_hallmark.learned import attribution_eb as AE
    from mirna_hallmark.learned import confounders as CF
    from mirna_hallmark.learned import data as LD
    out = CFG.OUTPUT_ROOT / "learned"

    card = pd.read_csv(out / "readouts_arm_edges.tsv", sep="\t")
    CARD = {g: dict(zip(sub["arm"], sub["beta"])) for g, sub in card.groupby("gene")}
    ago = AG.compute_ago_gate()["ago_capacity_z"]

    # ago ⊥ confounder block (gene-independent composition + prolif + purity)
    Cg = CF.build_C("tcga", ago.index, composition=True).apply(pd.to_numeric, errors="coerce")
    Cg = Cg.loc[:, Cg.std() > 1e-9].fillna(Cg.median(numeric_only=True))
    common = ago.index.intersection(Cg.index)
    Cm = np.column_stack([np.ones(len(common)), Cg.loc[common].to_numpy(float)])
    av = ago.loc[common].to_numpy(float)
    ago_perp = pd.Series(av - Cm @ np.linalg.lstsq(Cm, av, rcond=None)[0], index=common)
    r2 = 1 - np.var(ago_perp) / np.var(av)
    print(f"[ago⊥C] R²(ago~C)={r2:.3f}  ({r2:.0%} of ago is the confounder block)", flush=True)
    pro = Cg["mal_prolif"] if "mal_prolif" in Cg else Cg.iloc[:, 0]
    apos = {p: i for i, p in enumerate(ago_perp.index)}
    ago_raw = ago.reindex(ago_perp.index)

    # per-gene standardized arrays + positions into ago
    G, per_gene = [], []
    for g in sorted(CARD):
        try:
            Y, Xhe, C, _ = LD.assemble_gene(g, w_prior_source="ledger")
        except Exception:
            continue
        if Xhe is None or Xhe.shape[1] == 0:
            continue
        yr, Xz, cols = AE._prep(Y, Xhe, C)
        beta = np.array([CARD[g].get(a, 0.0) for a in cols], float)
        if not np.any(beta):
            continue
        M = Xz @ beta
        idx = list(Y.index)
        keep = [(k, apos[p]) for k, p in enumerate(idx) if p in apos]
        if len(keep) < 200:
            continue
        kk = np.array([k for k, _ in keep]); pp = np.array([j for _, j in keep])
        Mv, rv = M[kk], yr[kk]
        pv = pro.reindex([idx[k] for k in kk]).to_numpy(float)
        ok = np.isfinite(Mv) & np.isfinite(rv) & np.isfinite(pv)
        if ok.sum() < 200 or np.std(Mv[ok]) < 1e-9:
            continue
        zz = lambda v: (v - v.mean()) / (v.std() + 1e-12)
        G.append((zz(Mv[ok]), zz(rv[ok]), zz(pv[ok]), pp[ok]))
        per_gene.append({"gene": g, "n": int(ok.sum()),
                         "base_rho": stats.spearmanr(Mv[ok], rv[ok]).correlation})
    print(f"[genes] {len(G)} usable", flush=True)

    def aggregate(ago_vec):
        bs = []
        for Mz, rz, pz, pos in G:
            az = ago_vec[pos]
            if az.std() < 1e-9:
                bs.append(np.nan); continue
            az = (az - az.mean()) / (az.std() + 1e-12)
            D = np.column_stack([np.ones_like(Mz), Mz, az, pz, Mz * az, Mz * pz])
            bs.append(np.linalg.lstsq(D, rz, rcond=None)[0][4])
        bs = np.array(bs)
        return np.nanmean(bs), np.nanmean(bs > 0), bs

    perp = ago_perp.to_numpy(float)
    raw = ago_raw.to_numpy(float)
    obs_mean, obs_frac, obs_bs = aggregate(perp)
    raw_mean, _, _ = aggregate(raw)                          # over-control check
    print(f"\nOBSERVED  ago⊥C mean b_Mago={obs_mean:+.5f}  frac+={obs_frac:.1%}   |   raw-ago {raw_mean:+.5f} "
          f"(over-control shift {abs(obs_mean-raw_mean)/abs(raw_mean):.0%})", flush=True)

    rng = np.random.default_rng(seed)
    nm, nf = [], []
    for j in range(n_null):
        m, f, _ = aggregate(perp[rng.permutation(len(perp))])
        nm.append(m); nf.append(f)
    nm, nf = np.array(nm), np.array(nf)
    p_mean = (1 + (np.abs(nm - nm.mean()) >= abs(obs_mean - nm.mean())).sum()) / (n_null + 1)
    z_mean = (obs_mean - nm.mean()) / (nm.std() + 1e-12)
    z_frac = (obs_frac - nf.mean()) / (nf.std() + 1e-12)
    mde80 = 2.80 * nm.std()                                  # ~80%-power minimum detectable effect
    print(f"\nSHUFFLED-AGO NULL (N={n_null}) — the calibrated arbiter")
    print(f"  mean b : obs {obs_mean:+.5f}  null {nm.mean():+.5f}±{nm.std():.5f}  z={z_mean:+.2f}  perm p={p_mean:.3f}")
    print(f"  frac+  : obs {obs_frac:.1%}  null {nf.mean():.1%}±{nf.std():.1%}  z={z_frac:+.2f}")
    print(f"  power  : MDE@80% ≈ {mde80:.5f} = {mde80/abs(obs_mean):.1f}× the observed effect  ⇒ UNDERPOWERED")
    print(f"  VERDICT: {'DETECTED' if p_mean<0.05 else 'UNDETECTED (trend in predicted direction; retire the gate to rationale-only)'}")

    pd.DataFrame(per_gene).assign(b_Mago_perp=obs_bs).to_csv(out / "ago_interaction.tsv", sep="\t", index=False)
    pd.DataFrame({"null_mean_b": nm, "null_frac_pos": nf}).to_csv(out / "ago_interaction_null.tsv", sep="\t", index=False)
    print(f"\n[wrote] {out/'ago_interaction.tsv'} , {out/'ago_interaction_null.tsv'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-null", type=int, default=300)
    a = ap.parse_args()
    run(n_null=a.n_null)
