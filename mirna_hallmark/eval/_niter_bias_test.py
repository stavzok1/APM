"""Is the 1000/350 deviation NOISE or BIAS? The noise floor is random reseeding (signed mean should be ~0).
A short burn-in is systematic (the chain retains its init) -- so a signed mean != 0 is a BIAS, and a bias
the size of the noise is NOT the same thing as noise."""
import sys, numpy as np, pandas as pd
from mirna_hallmark.learned import readouts as R, spike_slab as SS
def fit(core, ni, bu, seed):
    b, sd, _, _ = SS._gibbs_posterior(core["Xz"], core["yr"], np.ones(core["p"]),
                                      n_iter=ni, burn=bu, seed=seed, return_samples=True)
    return b, sd
dB_floor, dB_cand, dS_floor, dS_cand = [], [], [], []
for g in sys.argv[1].split(","):
    core = R._posteriors(g, False, 2000, 700, 0, "family", discovery=False)
    if core is None or core["p"] < 2: continue
    b_ref, sd_ref = fit(core, 2000, 700, 0)
    for s in (1, 2, 3):
        b, sd = fit(core, 2000, 700, s); dB_floor += list(b - b_ref); dS_floor += list(sd - sd_ref)
    b, sd = fit(core, 1000, 350, 0);     dB_cand  += list(b - b_ref); dS_cand  += list(sd - sd_ref)
from scipy import stats
for name, fl, ca in [("beta", dB_floor, dB_cand), ("beta_sd", dS_floor, dS_cand)]:
    fl, ca = np.array(fl), np.array(ca)
    tf, tc = stats.ttest_1samp(fl, 0), stats.ttest_1samp(ca, 0)
    print(f"  {name}:")
    print(f"    FLOOR (reseed)   signed mean {fl.mean():+.3e}  p={tf.pvalue:.3f}  {'~0 (noise, as it must be)' if tf.pvalue>0.05 else 'BIASED?!'}")
    print(f"    1000/350         signed mean {ca.mean():+.3e}  p={tc.pvalue:.3f}  {'~0 => NOISE, not bias' if tc.pvalue>0.05 else '*** SYSTEMATIC BIAS ***'}")
    print(f"    |mean_cand| / sd(floor) = {abs(ca.mean())/fl.std():.3f}   (effect size vs the sampler's own jitter)")
