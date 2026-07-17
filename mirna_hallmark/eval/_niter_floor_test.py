"""Is readouts' N_ITER/BURN=2000/700 reducible? The decoy_bench precedent (200/80, verified against a
Monte-Carlo noise floor 7x larger) applies to a POINT ESTIMATE. readouts also emits `beta_sd`, and a
posterior SD is far more chain-length-sensitive than a mean -- and `identified` (|z|>2) divides by it.

DESIGN: the arbiter is the sampler's OWN seed-to-seed jitter at the reference config. If a shorter chain
moves beta/beta_sd LESS than reference-vs-reference-different-seed does, the cut is free.
"""
import sys, time, numpy as np, pandas as pd
from mirna_hallmark.learned import readouts as R, spike_slab as SS

genes = sys.argv[1].split(",")
CONFIGS = [(2000, 700), (1000, 350), (500, 200), (200, 80)]
REF = CONFIGS[0]

def fit(core, n_iter, burn, seed):
    b, sd, _, _ = SS._gibbs_posterior(core["Xz"], core["yr"], np.ones(core["p"]),
                                      n_iter=n_iter, burn=burn, seed=seed, return_samples=True)
    return b, sd

rows = []
for g in genes:
    core = R._posteriors(g, False, 2000, 700, 0, "family", discovery=False)
    if core is None or core["p"] < 2:
        continue
    b_ref, sd_ref = fit(core, *REF, 0)
    # THE NOISE FLOOR: same config, different seeds
    for s in (1, 2, 3):
        b, sd = fit(core, *REF, s)
        rows.append(dict(gene=g, p=core["p"], cfg="FLOOR(2000/700,seed%d)" % s, n_iter=REF[0],
                         d_beta=np.abs(b-b_ref).max(), d_sd=np.abs(sd-sd_ref).max(),
                         flip=int((( np.abs(b/np.where(sd>1e-12,sd,np.nan))>2) != (np.abs(b_ref/np.where(sd_ref>1e-12,sd_ref,np.nan))>2)).sum())))
    for (ni, bu) in CONFIGS[1:]:
        t=time.time(); b, sd = fit(core, ni, bu, 0); dt=time.time()-t
        rows.append(dict(gene=g, p=core["p"], cfg=f"{ni}/{bu}", n_iter=ni, secs=dt,
                         d_beta=np.abs(b-b_ref).max(), d_sd=np.abs(sd-sd_ref).max(),
                         flip=int((( np.abs(b/np.where(sd>1e-12,sd,np.nan))>2) != (np.abs(b_ref/np.where(sd_ref>1e-12,sd_ref,np.nan))>2)).sum())))
pd.DataFrame(rows).to_csv("/tmp/claude-207054/-sci-labs-michall-stavzok-APM/95cdd8f8-0eb4-4f39-aee5-5f5cfea5d600/scratchpad/niter_sweep.tsv", sep="\t", index=False)
print("  wrote niter_sweep.tsv", len(rows), "rows")
