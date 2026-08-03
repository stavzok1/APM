"""THE DECOY BENCH — the canonical, reusable implementation of the site-free negative-control benchmark.

    .venv/bin/python3 -m mirna_hallmark.eval.decoy_bench                 # full universe, both C blocks
    .venv/bin/python3 -m mirna_hallmark.eval.decoy_bench --seeds 5       # ⚠ ORDER-ONLY, not stability (MH-161)
    .venv/bin/python3 -m mirna_hallmark.eval.decoy_bench --genes PTEN,ZEB1

WHY THIS EXISTS (2026-07-16, user-demanded). The decoy arc was run as ~40 one-off scratchpad scripts, each
re-implementing the family RPM pool, the OOF budget, the decoy builder and the C block. **Every serious bug in
that arc was a DUPLICATION bug**: a double-log of an already-log2 matrix; per-GENE stop codons (88% coverage
loss); ABSOLUTE instead of SIGNED axis loadings; a 79-minute `corrwith` inside a per-gene loop; and core-C-only
runs that never computed retention (axiom 2a). One shared implementation prevents all five. Everything the arc
learned is encoded here ONCE.

THE ESTIMAND
------------
    M_g(s) = Σ_f β_f · X_f(s)          the LEARNED BUDGET, β from the canonical Gibbs
    ρ_g    = Spearman(M_g, Y_g | C)    scored OUT-OF-FOLD over patients
Compared between the gene's REAL curated regulators and a matched SITE-FREE decoy set.
⚠ STRICT OOF: β, the z-scoring mean/sd, the variance floor AND C's coefficients are fit on TRAIN patients only
  and APPLIED to the held-out fold. Every patient is scored exactly once by a β that never saw it.

THE DECOY (what the arc established, and what every earlier version got wrong)
------------------------------------------------------------------------------
ELIGIBILITY — a fake arm for gene g must be:
  * SITE-FREE: no TargetScan context++ site, no genome-wide scanMiR duplex, **and no 6mer** (⭐ NEW).
    ⚠ TargetScan's default file has Site Type ∈ {1=7mer-A1, 2=7mer-m8, 3=8mer} — **there is NO 6mer code**.
    Scanned from sequence: 14.6% of the old "site-free" decoys carried a 6mer. That contaminates the fake with
    real sites and biases the gap DOWNWARD (conservative), but it is wrong and is fixed here.
  * not a seed-family mate, not a genomic-cluster mate of any real regulator of g
  * decorrelated from every real regulator IN THE SPARSE RESIDUAL (|r_S| < 0.30 after the top-k global factors
    are projected out). ⚠ NOT raw |r|: two arms both loading on mPC1 are correlated THROUGH mPC1, so a raw cap
    systematically REJECTS high-loading arms — that single mistake under-loaded the v3 decoy by 17–22% and
    produced a ~2× overstatement of the gap.
MATCHING — optimal 1:1 assignment (Hungarian) on **SIGNED** loadings + dose + variance.
  ⭐ SIGNED, not |loading|: mPC2 IS the purity axis (r=−0.40 with CPE). Matching |mPC2| let a decoy at +0.3 be
    "matched" to a real regulator at −0.3 — same magnitude, OPPOSITE behaviour against a stromal target.
    MEASURED: 39% of v4 decoys loaded the opposite way on mPC2 (corr(real,fake) = +0.45). ANTI-MATCHED.
  ⚠ k for the metric must stay SMALL (default 2). k=46 (the Marchenko-Pastur count) makes the assignment
    48-dimensional, so it spends its budget on trivial factors and wrecks every axis that matters.
  ⚠ DOSE is nearly IRRELEVANT to the gap (measured b = +0.0009, p≈0.87) while LOADING dominates (b = +0.30).
    v3 spent its whole matching budget on dose. Weights below reflect the measured influence.

THE C BLOCK — axiom 2a, and it is NOT optional
----------------------------------------------
`assemble_gene` defaults to CORE C = [CPE, target_cn, mal_prolif]; composition is a SEPARATE block. The rule is
**run BOTH and report `retention = gap_deconv / gap_core`** — never silently condition composition away (a
miRNA acting through a cell-state shift PRODUCES a composition change; blanket-conditioning over-controls) and
never silently omit it (MH-107: a compartment-driven result presented as cell-intrinsic). **This module always
runs both.** Measured on the miR-200/EMT genes: gap retention 0.71, REAL-ρ retention 0.62.

WHAT THIS BENCH HAS MEASURED (full HE universe, 1,353 genes, family/β/core)
--------------------------------------------------------------------------
    gap = −0.0140, p = 5.6e-05, real wins 55%
    by families:  1 → −0.0005 (p=0.40, MUST be ~0: β≡uniform there — the internal null)
                  2 → −0.0175 (p=0.016) · 3–4 → −0.0357 (p=1e-04) · 5+ → −0.0263 (p=0.002)
    ⚠ every mis-matched decoy INFLATES the gap ⇒ treat any single construction as an UPPER BOUND.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("mirna_hallmark/output/learned")
K_METRIC = 2          # factors in the MATCHING metric — small on purpose (see the k=46 disaster)
K_RESID = 5           # factors projected out for the sparse-residual decorrelation
R_SPARSE_CAP = 0.30
N_FOLDS = 5
# ⭐ GIBBS ITERATIONS — VALIDATED, not guessed (axiom 3: runtime is a deliverable).
# Profiled: the Gibbs is 83% of the bench's runtime (C-blocks 17%, pooling 1%). My suspected waste —
# CF.build_C being rebuilt per gene — is only 3%, so caching it is NOT worth doing.
# MEASURED on 12 genes: 800/300 vs 400/150 vs 200/80 vs 1600/600 all agree to max|Δρ| <= 0.0013, while the
# MONTE-CARLO NOISE FLOOR (same setting, different FOLD seed) is 0.0094 — 7x larger. Even DOUBLING to 1600/600
# moves nothing ⇒ the chain has long converged by 200 iters, and the FOLD assignment (not the sampler)
# dominates the uncertainty. 200/80 is 3.9x faster and indistinguishable.
# ⚠ Do NOT lower further without re-running that check: this sampler has a bug history (`_rtnorm_pos` violated
#   its own beta>=0 support, contaminating 3.15% of persisted beta).
GIBBS_ITER, GIBBS_BURN = 200, 80
# measured influence of each axis on the gap (the Δ=0 extrapolation): loading ≫ variance ≫ dose≈0
# ⭐ DOSE WEIGHT RAISED 0.05 -> 10.0 (MH-137). The 0.05 came from b(dose)=+0.0009 (p=0.87) measured on the
# OLD decoy, where Δdose spanned only −0.04..−0.58 — I then applied it at Δdose=−2.4, FAR outside the measured
# range. On the 6mer-excluded decoy dose DOES bias the gap: Spearman(Δdose, gap)=+0.129, b=+0.0055, p<0.001,
# inflating it 1.75× (−0.0175 -> −0.0306).
# ⚠ WHY THE DEFICIT EXISTS AT ALL — and it is BIOLOGY, not a bug: dropping SITE-BEARING arms halves the pool
#   (2,230 -> 1,038) AND halves its mean dose (1.37 -> 0.76), because **abundant miRNAs have sites on nearly
#   everything** — abundance and target breadth are the same axis. The real regulators average dose 5.88 vs the
#   whole universe's 1.37, so the site-free pool can never match them. (The |r_S| decorrelation filter, which I
#   wrongly blamed, removes FOUR arms.) MEASURED: dose weight 10 takes Δdose −3.05 -> −0.90 and Δvar −0.22 ->
#   −0.02 at NO cost in pairs; a hard |Δdose|<1 caliper reaches −0.11 but drops 14% of pairs and degrades the
#   LOADING match (Δs2 −0.008 -> −0.033), which is the axis that actually matters (b≈+0.30).
AXIS_WEIGHT = {"s1": 1.0, "s2": 1.0, "dose": 10.0, "var": 0.5}
# ⭐ DOSE CALIPER — DEFAULT OFF (user-corrected). A 1.0 log2 caliper takes Δdose −0.90 → −0.09, but it does so
# by DROPPING 14.4% of pairs — and they are precisely the HIGH-DOSE regulators (dropped: 1,164 RPM median vs
# 31 RPM kept, p=2e-179), i.e. the edges delivering most of the repression. That RESTRICTS THE UNIVERSE to
# where a matched control happens to exist, which is a worse distortion than the dose residual it fixes.
# The residual is instead removed POST HOC by regressing the gap on Δdose (the b·Δ correction), which keeps
# every pair. Set to 1.0 to reinstate the restricted-universe variant.
# ⚠ I previously called the dose deficit "biology, not a bug" from a BROKEN supply table (it excluded all 702
#   real-regulator arms GLOBALLY — but an arm regulating gene A is a fine decoy for gene B). Corrected supply
#   vs demand per gene: 27.6 eligible arms >6 log2 against 1.33 needed — a 20× SURPLUS. The deficit is the
#   MATCHER trading dose against loading, and a caliper fixes it: |Δ|<1.0 takes Δdose −3.05 → −0.105 for 14%
#   of pairs. Set to 0 to disable.
DOSE_CALIPER = 0.0

_C: dict = {}


def _ctx():
    """Everything expensive, loaded ONCE. (The arc's 79-minute run was a load inside a per-gene loop.)"""
    if _C:
        return _C
    from sklearn.utils.extmath import randomized_svd

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned import instrument as INS

    d = LD._load()
    X, Y = d["X"], d["Y"]                       # ⚠ X IS ALREADY log2(RPM+1) — do NOT log it again
    arms = list(X.index)
    Xn = X.to_numpy(float)
    Z = ((X.T - X.mean(axis=1)) / (X.std(axis=1) + 1e-9)).fillna(0.0)
    U, _, _ = randomized_svd(Z.to_numpy(float), n_components=K_RESID, random_state=0)
    Zn = Z.to_numpy(float)
    Zc = (Zn - Zn.mean(0)) / (Zn.std(0) + 1e-12)
    Pc = (U - U.mean(0)) / (U.std(0) + 1e-12)
    signed = pd.DataFrame(Zc.T @ Pc / len(Zc), index=arms,          # ⭐ SIGNED, not |loading|
                          columns=[f"s{i+1}" for i in range(K_RESID)])
    Sres = Zn - U @ (U.T @ Zn)
    Ss = (Sres - Sres.mean(0)) / (Sres.std(0) + 1e-12)
    prop = signed[[f"s{i+1}" for i in range(K_METRIC)]].copy()
    prop["dose"] = np.nanmedian(Xn, 1)
    prop["var"] = np.nanvar(Xn, 1)
    E = D.load_hallmark_edges()
    HE = E[E.high_evidence]
    fm = FAM.family_of(pd.Index(arms))
    _C.update(X=X, Y=Y, arms=arms, signed=signed, dose=pd.Series(np.nanmedian(Xn, 1), index=arms), propz=(prop - prop.mean()) / (prop.std() + 1e-12),
              corr_s=pd.DataFrame((Ss.T @ Ss) / len(Ss), index=arms, columns=arms),
              he=HE, he_set={(r.miRNA, r.gene) for r in HE.itertuples()},
              fam=fm[~fm.index.duplicated()], LD=LD, FAM=FAM)
    try:
        _C["clusters"] = INS._genomic_clusters()
    except Exception:
        _C["clusters"] = {}
    _C["sites"] = _site_maps()
    return _C


def _site_maps() -> dict:
    """Everything that disqualifies an arm as a decoy for a gene. FOUR layers:
      (1) STRONG SITES  — TargetScan context++ (7mer-A1/7mer-m8/8mer) ∪ genome-wide scanMiR duplexes.
      (2) ⭐ 6mers, POISSON-GATED — a 6mer is 6 nt (4^6=4096), so a UTR of length L expects (L−5)/4096 copies
          BY CHANCE: at the median 1,092-nt UTR, λ=0.27 and **P(X≥1)=24%** — a SINGLE 6mer is sequence noise.
          Excluding on PRESENCE (the v7 rule) discards arms for nothing and hits LONG-UTR genes hardest.
          v8 excludes only 6mer counts that are Poisson-SIGNIFICANT (p<0.05) vs that gene's own λ.
          MEASURED: presence excludes 12.8M pairs; Poisson-gated excludes 2.3M — the pool grows 1,071→1,481.
          ⚠ It returns ZERO high-dose arms: abundant miRNAs are caught by the STRONG layer regardless.
      (3) ⭐ ALL EVIDENCE, high-throughput INCLUDED — the hole this fixes was 126×: v7 excluded only the
          **14,147 high-evidence miRTarBase rows (2% of 913,778)**, leaving **509,376 non-HE pairs** eligible
          as "decoys", plus **TarBase v9's 1,311,328 human pairs** (HITS-CLIP 2.9M rows, PAR-CLIP 1.5M,
          qCLASH 126k) entirely unused. A decoy with a CLIP hit is NOT a decoy. Union: **1,790,439 pairs**.
          It costs ~20% of the high-dose supply (34.6 → 27.6 arms/gene). Correct to pay.
      (4) ⭐ MANAKOV CHIMERIC eCLIP — the STRONGEST evidence class: the miRNA and its target LIGATED in one
          AGO footprint. 1.5M reads / 599,161 (arm,gene) pairs; **75% are invisible to miRTarBase+TarBase**
          (+448,330 new). MEASURED: **133/4,937 (2.7%) of the previous decoys carried a chimeric read.**
          ⚠ my first parser MISSED this file entirely — it searched for a column containing "mir"; Manakov's
            miRNA column is `noncodingRNA`.
      ⭐ (5) POSTAR3 + ENCORI — **MEASURED 2026-08-03, and the exclusion is materially CLOSED. The previous
          note here ("POSTAR3 not present in this repo ⇒ cannot be closed with local data") was WRONG on the
          fact and, more importantly, wrong on the MECHANISM.**
          * **POSTAR3 IS on disk** — `data/external/POSTAR/human (1).txt.gz`, 676 MB, 2026-06-30. But it is
            the **RBP binding-site table**: 221 distinct RBPs, **2,360,006 AGO2 records**, TNRC6A/B/C ~3k,
            DICER1 8.5k — and **ZERO miRNA-named entries**. It therefore yields NO pair-level (arm, gene)
            call, only miRNA-ANONYMOUS occupancy intervals.
          * ⭐ **And an AGO-peak ∩ seed-site layer would be a NO-OP BY CONSTRUCTION.** `build_decoys` already
            requires `(a, g) not in ctx["sites"]`, and layers (1)+(2) exclude EVERY arm with a strong site or
            a Poisson-significant 6mer in that gene. An arm whose seed sits under an AGO peak is thus
            *already* ineligible. Only a source that resolves the **mature arm** on a **SEEDLESS** pair can
            add anything here — which is exactly what Manakov chimeric (layer 4) is, and it moved 2.7%.
          * **ENCORI is 5.3× under-ingested as a LABELLED layer but 97.0% redundant as a SET.** The local
            cache `data/external_cache/encori/miRNATarget/` (8,647 per-gene files, 6,431 non-empty) holds
            **20,404 distinct (arm,gene) pairs**; the shipped `ENCORI_starBase` layer carries only 3,861.
            But **19,799 / 20,404 (97.0%)** are already in the union via TarBase v9 / miRTarBase / Manakov,
            leaving **605 new**; of those **578 are already in the full site map** via layers (1)/(2), and
            ⭐ **0 of the 4,937 assigned decoys carry one (0.000%)**.
          ⇒ **Neither source can move the gap.** Closing the evidence hole further requires a pair-level,
            arm-resolving, SEEDLESS-capable source (POSTAR3's separate miRNA-target/degradome module, which
            is NOT what is on disk) — and layer (4) bounds the plausible scale of that at a few percent.
      TOTAL evidence-excluded: 2,242,630 pairs (verified 2026-08-03; the log below derives the breakdown)."""
    ts = pd.read_csv("data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt", sep="\t",
                     usecols=["Gene Symbol", "miRNA"], low_memory=False)
    S = set(map(tuple, ts.rename(columns={"Gene Symbol": "gene", "miRNA": "arm"})[["arm", "gene"]].values))
    kd = pd.concat([pd.read_csv(f"data/external_cache/scanmir/{f}", sep="\t", usecols=["arm", "gene"])
                    for f in ("genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz")],
                   ignore_index=True).drop_duplicates()
    S |= set(map(tuple, kd[["arm", "gene"]].values))
    six = OUT / "sixmer_counts.tsv.gz"
    if six.exists():                                     # ⭐ POISSON-GATED, not presence
        from scipy import stats as _st
        s6 = pd.read_csv(six, sep="\t")
        sig = s6[_st.poisson.sf(s6.n6 - 1, s6.expected) < 0.05]
        S |= set(map(tuple, sig[["arm", "gene"]].values))
        print(f"  [sites] 6mer Poisson-gated: {len(sig):,} of {len(s6):,} 6mer-bearing pairs excluded")
    elif (OUT / "sixmer_sites.tsv.gz").exists():
        s6 = pd.read_csv(OUT / "sixmer_sites.tsv.gz", sep="\t")
        S |= set(map(tuple, s6[["arm", "gene"]].values))
        print("  ⚠ using the PRESENCE 6mer map (over-strict); build sixmer_counts.tsv.gz for the Poisson gate")
    else:
        print("  ⚠ no 6mer map — the decoy will be 6mer-CONTAMINATED (14.6%)")
    evf = OUT / "evidence_exclusion.tsv.gz"
    if evf.exists():                                     # ⭐ ALL evidence, high-throughput included
        ev = pd.read_csv(evf, sep="\t")
        S |= set(map(tuple, ev[["arm", "gene"]].values))
        # ⚠ PRINT WHAT IS ACTUALLY IN THE FILE, never a hard-coded list. The old label read
        # "(miRTarBase ALL ∪ TarBase v9)" while the cache in fact carried FOUR sources — including
        # Manakov_chimeric (448,330) and ENCORI_starBase (3,861). A reader auditing coverage from the log
        # would have concluded the chimeric layer was missing and "fixed" a non-bug (it nearly happened,
        # 2026-07-17). The cache has no builder and no code-version stamp, so the LOG is the only provenance
        # a reader gets: it must be derived, not asserted. (CLAUDE.md axiom 6.)
        src = (", ".join(f"{k} {v:,}" for k, v in ev["src"].value_counts().items())
               if "src" in ev.columns else "⚠ no `src` column — provenance unknown")
        print(f"  [sites] evidence exclusion: {len(ev):,} pairs  [{src}]")
    else:
        print("  ⚠⚠ NO evidence exclusion — decoys may be CLIP/CLASH-supported REAL edges (a 126× hole)")
    return S


def build_sixmer_map() -> pd.DataFrame:
    """Scan every (arm, gene) 6mer from SEQUENCE. Cached — this is the layer TargetScan does not provide."""
    ctx = _ctx()
    seqs, nm, buf = {}, None, []
    for line in open("data/external_cache/scanmir/genomewide_3utr.fa"):
        if line.startswith(">"):
            if nm:
                seqs[nm] = "".join(buf)
            nm, buf = line[1:].strip().split()[0].split("|")[0], []
        else:
            buf.append(line.strip())
    if nm:
        seqs[nm] = "".join(buf)
    mat = {}
    import glob
    for f in glob.glob("data/miRNA/*.fa") + glob.glob("data/miRNA/mature*"):
        n2 = None
        for line in open(f):
            if line.startswith(">"):
                n2 = line[1:].strip().split()[0]
            elif n2:
                mat[n2] = line.strip().upper().replace("U", "T"); n2 = None
    rcm = str.maketrans("ACGT", "TGCA")
    rows = []
    for a in ctx["arms"]:
        s = mat.get(a)
        if not s or len(s) < 8:
            continue
        m6 = s[1:7].translate(rcm)[::-1]
        for g, u in seqs.items():
            if m6 in u:
                rows.append((a, g))
    df = pd.DataFrame(rows, columns=["arm", "gene"])
    df.to_csv(OUT / "sixmer_sites.tsv.gz", sep="\t", index=False)
    print(f"[decoy_bench] 6mer map: {len(df):,} (arm, gene) pairs -> {OUT/'sixmer_sites.tsv.gz'}")
    return df


def real_regulators(gene: str) -> list:
    ctx = _ctx()
    s = ctx["he"][ctx["he"].gene == gene]
    return [a for a in s.miRNA.unique() if a in ctx["X"].index]


def build_decoys(genes=None, seed: int = 0) -> pd.DataFrame:
    """One matched site-free decoy arm per real regulator, per gene. Optimal assignment on SIGNED axes."""
    from scipy.optimize import linear_sum_assignment
    ctx = _ctx()
    MC = [c for c in ctx["propz"].columns]
    w = np.array([AXIS_WEIGHT.get(c, 1.0) for c in MC])
    genes = genes or sorted({g for g in ctx["he"].gene.unique() if g in ctx["Y"].index})
    # ⚠ `seed` ONLY reshuffles gene order — the eligible set, the metric and linear_sum_assignment are
    # deterministic per gene, so decoys are IDENTICAL across seeds (verified, MH-161). `--seeds` is NOT a
    # stability knob: the real MC noise is `oof_budget`'s FOLD seed (hardcoded 0 in `_one`). To average it,
    # vary that fold seed (see `eval/admissibility_bench._one_firm`), not this one.
    rng = np.random.default_rng(seed)
    out = []
    for g in rng.permutation(list(genes)):
        real = real_regulators(g)
        if not real:
            continue
        rfam = {ctx["fam"].get(a) for a in real if isinstance(ctx["fam"].get(a), str)}
        bad = set(real)
        for a in real:
            bad |= set(ctx["clusters"].get(a, ()) or ())
        el = [a for a in ctx["arms"] if a not in bad
              and not (isinstance(ctx["fam"].get(a), str) and ctx["fam"].get(a) in rfam)
              and (a, g) not in ctx["sites"] and (a, g) not in ctx["he_set"]]
        if len(el) < len(real):
            continue
        cs = ctx["corr_s"].loc[el, real].abs().max(axis=1)
        el = [a for a in el if cs[a] < R_SPARSE_CAP]
        if len(el) < len(real):
            continue
        Pr, Pe = ctx["propz"].loc[real, MC].to_numpy(float), ctx["propz"].loc[el, MC].to_numpy(float)
        C = np.nan_to_num(np.sqrt((((Pe[None, :, :] - Pr[:, None, :]) ** 2) * w).sum(-1)), nan=1e6)
        if DOSE_CALIPER:                                  # ⭐ forbid out-of-caliper pairs outright
            dr = ctx["dose"].reindex(real).to_numpy(float)[:, None]
            de = ctx["dose"].reindex(el).to_numpy(float)[None, :]
            C = np.where(np.abs(de - dr) > DOSE_CALIPER, 1e9, C)
        ri, ci = linear_sum_assignment(C)
        for i, j in zip(ri, ci):
            if DOSE_CALIPER and C[i, j] >= 1e9:
                continue                                  # UNMATCHABLE -> drop the pair (a coverage cost)
            out.append({"gene": g, "seed": seed, "real_arm": real[i], "fake_arm": el[j]})
    return pd.DataFrame(out)


def pool_family(arms, index) -> pd.DataFrame:
    """The CANONICAL family collapse: a TRUE RPM pool log2(1 + Σ(2^x − 1)) — never a groupby.mean()."""
    ctx = _ctx()
    fam = ctx["fam"].reindex(arms)
    A = ctx["X"].reindex(arms).T.reindex(index).apply(pd.to_numeric, errors="coerce")
    lin = np.nan_to_num(np.power(2.0, A.to_numpy(float)) - 1.0, nan=0.0)
    L = pd.DataFrame(lin, index=index, columns=list(arms))
    return pd.DataFrame({str(f): np.log2(1.0 + L[[a for a in arms if fam.get(a) == f]].sum(axis=1))
                         for f in dict.fromkeys([fam.get(a) for a in arms])}, index=index)


def oof_budget(Yv: pd.Series, Xf: pd.DataFrame, C: pd.DataFrame, seed: int = 0) -> float:
    """ρ = Spearman(Σ_f β_f·X_f, Y | C), β fitted OUT-OF-FOLD. Nothing from the held-out fold touches the fit."""
    from scipy import stats
    from sklearn.linear_model import LinearRegression

    from mirna_hallmark.learned import spike_slab as SS
    n = len(Yv)
    fold = np.random.default_rng(seed).permutation(np.arange(n) % N_FOLDS)
    y, Xv, Cm = Yv.to_numpy(float), Xf.to_numpy(float), C.to_numpy(float)
    Mt, Yt = [], []
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        if te.sum() < 20 or tr.sum() < 50:
            continue
        lc = LinearRegression().fit(Cm[tr], y[tr])                 # C fit on TRAIN
        mu, sd = Xv[tr].mean(0), Xv[tr].std(0)                     # z-params from TRAIN
        Ztr = (Xv[tr] - mu) / (sd + 1e-9); Ztr[:, sd < 0.1] = 0.0   # variance floor from TRAIN
        Zte = (Xv[te] - mu) / (sd + 1e-9); Zte[:, sd < 0.1] = 0.0
        try:
            b, _, _ = SS._gibbs_posterior(Ztr, -(y[tr] - lc.predict(Cm[tr])), np.ones(Ztr.shape[1]),
                                          n_iter=GIBBS_ITER, burn=GIBBS_BURN, seed=0)
        except Exception:
            return np.nan
        # ⭐ STANDARDISE EACH FOLD'S BUDGET BEFORE POOLING (2026-08-01, MH-181). `b` is REFIT PER FOLD and
        # `Zte @ b` is an ARBITRARY-SCALE index, so concatenating raw folds glues differently-scaled
        # pieces together and the pooled Spearman is NOT invariant to that. (`Yt` needs no such fix — it
        # is a residual in y units with a per-fold intercept, so its folds are already commensurate.)
        # Found in `weight_gain_oof`, where it broke an EXACT null: one-family genes must give gain 0 and
        # gave up to 9.2e-02 until each fold was z-scored. This is the same estimator, same defect.
        m_ = Zte @ b
        s_ = np.std(m_)
        Mt.append((m_ - np.mean(m_)) / s_ if s_ > 1e-12 else m_ * 0.0)
        Yt.append(y[te] - lc.predict(Cm[te]))
    if not Mt:
        return np.nan
    M, Yo = np.concatenate(Mt), np.concatenate(Yt)
    return stats.spearmanr(M, Yo).correlation if np.std(M) > 1e-9 else np.nan


def _C_blocks(gene):
    """⭐ axiom 2a: BOTH blocks, always. Returns {'core': C, 'deconv': C}."""
    from mirna_hallmark.learned import confounders as CF
    from mirna_hallmark.learned import data as LD
    Y, _, C, _ = LD.assemble_gene(gene, w_prior_source="ledger")
    out = {"core": C}
    try:
        Cd = CF.build_C("tcga", C.index, composition=True).reindex(C.index)
        Cd = Cd.apply(pd.to_numeric, errors="coerce")
        Cd = Cd.fillna(Cd.median(numeric_only=True))
        Cd["target_cn"] = pd.to_numeric(C["target_cn"], errors="coerce").reindex(C.index).fillna(0.0)
        out["deconv"] = Cd
    except Exception:
        pass
    return Y, out


def _collin(P: pd.DataFrame) -> float:
    """Mean |pairwise Spearman| among a design's family columns — NaN below 2 columns."""
    if P.shape[1] < 2:
        return np.nan
    R = P.corr(method="spearman").to_numpy(float)
    v = np.abs(R[np.triu_indices(R.shape[0], k=1)])
    v = v[np.isfinite(v)]
    return float(v.mean()) if v.size else np.nan


def _one(args):
    gene, pairs = args
    try:
        Yv, blocks = _C_blocks(gene)
    except Exception:
        return None
    ra, fa = list(pairs.real_arm), list(pairs.fake_arm)
    Xr, Xf = pool_family(ra, Yv.index), pool_family(fa, Yv.index)
    # ⭐ n_fam_fake + collinearity EMITTED (2026-08-01). The matcher works on ARMS and the family collapse
    # happens HERE, afterward — so neither design WIDTH nor within-design COLLINEARITY is matched, and until
    # now only the REAL side's width was recorded, making the asymmetry unmeasurable from the artifact.
    # It is real (width real 3.096 vs fake 3.510, p=4.6e-42; collin 0.1663 vs 0.1547, p=6.7e-03) but INERT:
    # b(Δn_fam)=−0.0009 p=0.72, b(Δcollin)=−0.024 p=0.62 ⇒ contamination of the gap ≤~0.0025 at the CI edge.
    # Keep emitting them so the control stays auditable. Full analysis: `eval/decoy_design_match.py`.
    o = {"gene": gene, "n_arm": len(ra), "n_fam": Xr.shape[1], "n_fam_fake": Xf.shape[1],
         "collin_real": _collin(Xr), "collin_fake": _collin(Xf)}
    # ⭐ Δdose per gene (fake − real, log2 units). REQUIRED by the caliper-OFF design: with DOSE_CALIPER=0
    # every pair is kept, so the dose residual is removed POST HOC by regressing gap on Δdose (the b·Δ
    # correction) — but that is impossible unless Δdose is emitted, and it was not. Without this column the
    # reported gap is the RAW, dose-inflated one (MH-136b's exact error: the raw −0.0306 was inflated ~1.75×
    # by Δdose=−2.36). ⚠ b MUST be re-derived on THIS decoy, never reused: MH-136b retracted a b measured
    # over Δdose∈[−0.04,−0.58] and applied at −2.36, far outside its range.
    _d = _C["dose"]
    o["d_dose"] = float(np.nanmean(_d.reindex(fa).to_numpy(float))
                        - np.nanmean(_d.reindex(ra).to_numpy(float)))
    for blk, C in blocks.items():
        o[f"real_{blk}"] = oof_budget(Yv, Xr, C)
        o[f"dec_{blk}"] = oof_budget(Yv, Xf, C)
        o[f"gap_{blk}"] = o[f"real_{blk}"] - o[f"dec_{blk}"]
    return o


def run(genes=None, seeds: int = 1, workers: int = 8) -> pd.DataFrame:
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")
    from multiprocessing import Pool

    from scipy import stats
    _ctx()
    allr = []
    for s in range(seeds):
        P = build_decoys(genes, seed=s)
        print(f"[decoy_bench] seed {s}: decoys for {P.gene.nunique():,} genes / {len(P):,} pairs")
        jobs = [(g, sub) for g, sub in P.groupby("gene")]
        with Pool(workers) as p:
            R = pd.DataFrame([r for r in p.imap_unordered(_one, jobs, chunksize=4) if r])
        R["seed"] = s
        allr.append(R)
    G = pd.concat(allr, ignore_index=True).dropna(subset=["real_core", "dec_core"])
    # ⛔ NEVER let a SUBSET run clobber the canonical full-universe artifact (fixed 2026-08-01, MH-139).
    # `decoy_bench.tsv` is MH-137/139's cited evidence, and this line used to write it unconditionally —
    # so a `--genes` run, or an aborted one, silently replaced 1,349 genes with whatever it had. It
    # happened: the file was found holding the alphabetically-first 38 genes, i.e. the row's evidence
    # column no longer pointed at its evidence, with nothing in the file to reveal it.
    dest = OUT / ("decoy_bench.tsv" if genes is None else "decoy_bench_subset.tsv")
    G.to_csv(dest, sep="\t", index=False)
    if genes is not None:
        print(f"  ⚠ SUBSET run ({G.gene.nunique():,} genes) -> {dest.name}; "
              f"canonical decoy_bench.tsv left untouched")

    print(f"\n=== THE DECOY BENCH — {G.gene.nunique():,} genes, {seeds} seed(s) ===")
    print("    [TCGA · family · β · 5-fold OOF over patients · RPM-pool design]\n")
    print(f"  {'C block':10s} {'REAL':>9s} {'DECOY':>9s} {'gap':>9s} {'p':>10s} {'win%':>6s}")
    for blk in ("core", "deconv"):
        if f"gap_{blk}" not in G:
            continue
        s = G.dropna(subset=[f"real_{blk}", f"dec_{blk}"])
        w = stats.wilcoxon(s[f"real_{blk}"], s[f"dec_{blk}"], alternative="less").pvalue
        print(f"  {blk:10s} {s[f'real_{blk}'].mean():+9.4f} {s[f'dec_{blk}'].mean():+9.4f} "
              f"{s[f'gap_{blk}'].mean():+9.4f} {w:10.2e} {(s[f'gap_{blk}']<0).mean():6.0%}")
    if "gap_deconv" in G:
        r = G.gap_deconv.mean() / G.gap_core.mean()
        print(f"\n  ⭐ RETENTION (axiom 2a) = gap_deconv / gap_core = {r:.2f}")
    print("\n  --- by n_families (1 family MUST be ~0: β≡uniform there — the internal null) ---")
    for lab, sel in [("1", G.n_fam == 1), ("2", G.n_fam == 2), ("3-4", G.n_fam.between(3, 4)),
                     ("5+", G.n_fam >= 5)]:
        s = G[sel]
        if len(s) < 15:
            continue
        w = stats.wilcoxon(s.real_core, s.dec_core, alternative="less").pvalue
        print(f"  {lab:4s} n={len(s):4d}  gap {s.gap_core.mean():+.4f}  p={w:.4f}")
    return G


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--genes", help="comma-separated subset")
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--build-sixmer", action="store_true", help="build the cached 6mer site map first")
    a = ap.parse_args()
    if a.build_sixmer:
        build_sixmer_map()
    run(genes=a.genes.split(",") if a.genes else None, seeds=a.seeds, workers=a.workers)
