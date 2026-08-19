"""THE CPTAC EDGE BLOCK — the same per-edge coupling question as the TCGA card, in the CPTAC cohorts,
PLUS the two layers TCGA cannot supply: protein, and the protein-vs-mRNA translational residual.

    .venv/bin/python3 -m mirna_hallmark.learned.cptac_card

WHY (user-asked, 2026-08-01). The cards carry TCGA coupling per edge. CPTAC adds (a) an INDEPENDENT patient
set, and (b) a readout TCGA does not have at all — **protein**, which is what miRNAs actually control.

THREE LAYERS PER EDGE, per cohort (all partial Spearman, arm abundance vs the layer):
    rho_rna    arm vs target **mRNA**      — the direct analogue of the TCGA card's coupling column
    rho_prot   arm vs target **protein**   — ⭐ the layer TCGA cannot give
    rho_disc   arm vs the **protein-vs-mRNA residual** — translational repression beyond mRNA loss
Each in a composition-ADJUSTED and a RAW form, plus their gated retention (axiom 2a: FLAG, don't delete).

TWO COHORTS, AND THEY ARE NOT INTERCHANGEABLE — the prefixes are the warning:
    `cptac_prosp_*`  CPTAC-2 prospective, n≈101 shared samples. **INDEPENDENT patients AND its own
                     CPTAC-2 miRNA-seq** ⇒ a genuine out-of-cohort test.
    `cptac_t105_*`   CPTAC TCGA-105, n≈105. **CIRCULAR for the model: the SAME patients as the TCGA fit,
                     and its miRNA IS the TCGA miRNA.** Only its PROTEIN is new data ⇒ read it as a
                     LAYER transfer (mRNA→protein), never as cohort replication.

⚠⚠ WHAT THESE COLUMNS ARE NOT.
  * **NOT a coupling lever.** Protein carries **4–6%** of the mRNA channel's Fisher information about β,
    **≤7.6% at any a_g** (pre-registered ceiling, MH-103). Adding protein does not sharpen β and must not
    be sold as doing so. It is an independent READOUT, not a better measurement of the same thing.
  * **NOT per-edge significance.** n≈100–120 here against n≈1000 in TCGA, and MH-123/124 measured the
    per-edge null as 3–4× too narrow even at TCGA's n. Use these to STRATIFY and to CROSS-CHECK, never to
    declare a single edge validated.
  * **NOT decoy-controlled.** MH-130d: a site-free fitted fake reaches **99%** of the protein transfer.
    A negative `rho_prot` is real arithmetic; it is **not** evidence the curated edge is right.
  * **`βᵗ` IS DEAD (MH-103)** — `rho_disc` is a descriptive residual column, NOT a revival of the
    translational-repression latent, which was falsified at n=101. Do not rebuild it from this.

COMPOSITION. Both cohorts get the same `confounders.build_C` block `ood_protein._cov` uses (CPTAC for
prospective, TCGA for tcga105), and BOTH the adjusted and RAW value are emitted so retention is visible
rather than silently conditioned — MH-107's exact defect was a channel that only ever ran one block.
Measured context: on gold edges the median protein retention is 0.414 with 241/495 `composition_explained`
(MH-172); on the 7 canonical hub genes it is 0.797 (MH-171). Expect layer- and universe-dependence.
"""
from __future__ import annotations

import os
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "edge_cptac_card.tsv"
DEST_GENE = OUT / "gene_cptac_card.tsv"
# the BUILD INTERMEDIATE (see card_context) — used only to enumerate the HE (gene, arm) universe + beta
CARD = OUT / "edge_card_base.tsv"

COHORTS = {"prospective": "cptac_prosp", "tcga105": "cptac_t105"}
LAYERS = {"rna_z": "rna", "protein_z": "prot", "protein_resid": "disc"}
from mirna_hallmark.config import RHO_GATE   # ⭐ ONE home (MH-257); was a local 0.05 here


_FAMBETA: dict = {}


def _family_beta() -> dict:
    """(gene, family) -> the FAMILY-rung β, from `readouts.run(level="family")`. Memoised.

    ⚠ This is a DIFFERENT number from the card's `beta`, which is the same estimator at `level="arm"`
    (MH-191/192). A family-pooled design column must be weighted by the family-rung β; using an arm's β
    there is a unit mismatch, and picking the FIRST arm's β is additionally arbitrary."""
    if not _FAMBETA:
        for f in ("gene_family_card.tsv", "readouts_edges.tsv"):
            p = OUT / f
            if p.exists():
                d = pd.read_csv(p, sep="\t", low_memory=False)
                if {"gene", "family", "beta"} <= set(d.columns):
                    _FAMBETA.update({(g, str(fm)): b for g, fm, b in
                                     zip(d.gene, d.family, d.beta)})
                    break
    return _FAMBETA


def _resid_rows(M: np.ndarray, D: np.ndarray) -> np.ndarray:
    """Residualise every ROW of M on the design D, NaN-safe. Done ONCE per matrix, not per edge."""
    out = np.full(M.shape, np.nan)
    for i in range(M.shape[0]):
        v = M[i]
        m = np.isfinite(v)
        if m.sum() > D.shape[1] + 2:
            Dm = D[m]
            out[i, m] = v[m] - Dm @ np.linalg.lstsq(Dm, v[m], rcond=None)[0]
    return out


def _sp(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return spearmanr(a[m], b[m]).correlation if m.sum() > 10 else np.nan


def build(edges: pd.DataFrame | None = None) -> pd.DataFrame:
    from mirna_hallmark.eval import cptac_validation as V
    from mirna_hallmark.learned.eval import ood_protein as OP

    if edges is None:
        edges = pd.read_csv(CARD, sep="\t", usecols=["gene", "arm"]).drop_duplicates()
    res = edges.copy()

    for cohort, pref in COHORTS.items():
        d = OP._cptac(cohort)
        arms = d["arms"]
        L = V.load_cptac_layers(cohort)
        samples = [s for s in arms.columns if s in L["protein_z"].columns]
        Dcov = OP._cov(cohort, samples)                       # SAME block ood_protein uses (reuse contract)
        A = arms[samples].to_numpy(float)
        Ar = _resid_rows(A, Dcov)
        ai = {a: i for i, a in enumerate(arms.index)}

        Lr, Lraw, gi = {}, {}, {}
        for key, tag in LAYERS.items():
            M = L[key].reindex(columns=samples)
            Lraw[tag] = M.to_numpy(float)
            Lr[tag] = _resid_rows(Lraw[tag], Dcov)
            gi[tag] = {g: i for i, g in enumerate(M.index)}

        n_ok = 0
        cols = {f"{pref}_{k}": np.full(len(res), np.nan) for k in
                ("n", "rho_rna", "rho_rna_raw", "ret_rna", "rho_prot", "rho_prot_raw", "ret_prot",
                 "rho_disc", "rho_disc_raw")}
        cls = pd.Series(pd.NA, index=res.index, dtype="object")
        for r, (g, a) in enumerate(zip(res.gene, res.arm)):
            j = ai.get(a)
            if j is None:
                continue
            got = False
            for tag in ("rna", "prot", "disc"):
                i = gi[tag].get(g)
                if i is None:
                    continue
                adj, raw = _sp(Ar[j], Lr[tag][i]), _sp(A[j], Lraw[tag][i])
                cols[f"{pref}_rho_{tag}"][r] = adj
                cols[f"{pref}_rho_{tag}_raw"][r] = raw
                if tag in ("rna", "prot") and np.isfinite(raw) and abs(raw) >= RHO_GATE:
                    cols[f"{pref}_ret_{tag}"][r] = adj / raw
                got = True
            if got:
                cols[f"{pref}_n"][r] = len(samples)
                n_ok += 1
        rp = cols[f"{pref}_ret_prot"]
        cls[rp >= 0.7] = "cell_intrinsic"
        cls[(rp >= 0.4) & (rp < 0.7)] = "partial"
        cls[rp < 0.4] = "composition_explained"
        for k, v in cols.items():
            res[k] = v
        res[f"{pref}_comp_class_prot"] = cls
        print(f"[cptac_card] {cohort:12s} n={len(samples):4d} samples · {Dcov.shape[1]-1} covariates · "
              f"{n_ok:,}/{len(res):,} edges scored ({n_ok/len(res):.0%})")

    res.to_csv(DEST, sep="\t", index=False)
    print(f"-> {DEST}")
    return res


def report(res: pd.DataFrame) -> None:
    print("\n=== CPTAC EDGE BLOCK ===")
    for cohort, pref in COHORTS.items():
        tag = "INDEPENDENT patients + own miRNA" if cohort == "prospective" \
              else "CIRCULAR (same patients + TCGA miRNA; only PROTEIN is new)"
        print(f"\n[{cohort}] {tag}")
        for lay, lab in (("rna", "mRNA  (TCGA analogue)"), ("prot", "protein ⭐ new layer"),
                         ("disc", "prot-vs-mRNA resid")):
            a, r = res[f"{pref}_rho_{lay}"], res[f"{pref}_rho_{lay}_raw"]
            k = a.notna()
            if not k.any():
                continue
            print(f"  {lab:24s} n={int(k.sum()):5d}  mean rho adj {a[k].mean():+.4f} (raw {r[k].mean():+.4f})"
                  f"  neg {(a[k] < 0).mean():.0%}")
        ret = res[f"{pref}_ret_prot"].dropna()
        if len(ret):
            print(f"  protein composition retention: median {ret.median():.3f}  "
                  f"cell_intrinsic {int((ret>=0.7).sum())} · partial {int(((ret>=0.4)&(ret<0.7)).sum())} · "
                  f"composition_explained {int((ret<0.4).sum())}  (n={len(ret)})")


if __name__ == "__main__":
    report(build())


# --------------------------------------------------------------------------------------------------
# GENE-LEVEL AGGREGATE — "aggregate up" the edge block into one partial correlation per gene
# --------------------------------------------------------------------------------------------------
def build_gene() -> pd.DataFrame:
    """Per GENE: the AGGREGATE miRNA budget vs each CPTAC layer, plus a sum-abundance reference.

        agg_g(s)   = SUM_a beta_a * X_a(s)     beta = the DENSE GIBBS posterior mean, FIT ON TCGA mRNA
                                               at n~1000 — the canonical estimator (MH-115/§6b).
        abund_g(s) = SUM_a X_a(s)              the same arms, UNWEIGHTED — the reference baseline.

    Both are scored by partial Spearman against `rna_z`, `protein_z`, and `protein_resid` (protein
    residualised on its OWN mRNA — the translational layer), composition-ADJUSTED and RAW.
    A TCGA reference is computed with the IDENTICAL construction so the transfer is measurable rather
    than asserted.

    ⚠⚠ **`agg_beats_abund` IS A REFERENCE, NOT A CONTROL.** It compares a FITTED beta against an
    UNFITTED sum. MH-115/127 measured that a *fitted matched DECOY* also beats unfitted abundance, and
    decoy betas transport at ~80% of real. Beating abundance is therefore the FLOOR, not evidence the
    curated weights are right — the honest benchmark is a fitted fake, which is `eval/decoy_bench`.

    ⛔ **NOT the weights `ood_protein.fit_M` uses.** That calls `regression.fit_gene` = the ADAPTIVE
    LASSO, which MH-115 RETIRED to a baseline (verified in code: `Lasso(positive=True)` on
    prior-scaled columns). This function deliberately uses the Gibbs `beta` off the card instead. Same
    "two estimators, one label" trap as MH-145 — read the call, not the name.
    """
    from mirna_hallmark.eval import cptac_validation as V
    from mirna_hallmark.learned import confounders as CF
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.eval import ood_protein as OP

    card = pd.read_csv(CARD, sep="\t", usecols=["gene", "arm", "beta"], low_memory=False)
    card = card.dropna(subset=["beta"])
    genes = sorted(card.gene.unique())
    res = pd.DataFrame({"gene": genes})

    def _score(Xmat, layers, Dcov, prefix, cols):
        """Fill `cols` for one cohort: aggregate + abundance vs every layer, adjusted and raw.

        ⭐ RUNG-CORRECT (fixed 2026-08-01, MH-179). β is ESTIMATED on **z-scored SEED-FAMILY columns**
        (`readouts.py`), so it must be APPLIED to them. This previously did `Σ β_a · X_a` on RAW
        ARM-LEVEL log2-RPM — a units AND rung mismatch that halved the aggregate's measured coupling
        (OOF ρ_agg −0.145 on the correct construction vs −0.066 on the raw-arm one). Now: pool arms to
        families (a TRUE RPM pool, `decoy_bench.pool_family`), z-score, then apply the family β.
        The UNWEIGHTED reference is the sum of the SAME z-scored family columns, so the only difference
        between the two aggregates is β."""
        from mirna_hallmark.eval import decoy_bench as DB
        ai = {a: i for i, a in enumerate(Xmat.index)}
        Xr = _resid_rows(Xmat.to_numpy(float), Dcov)
        Xn = Xmat.to_numpy(float)
        Lr, Lraw, gi = {}, {}, {}
        for tag, M in layers.items():
            Lraw[tag] = M.to_numpy(float)
            Lr[tag] = _resid_rows(Lraw[tag], Dcov)
            gi[tag] = {g: i for i, g in enumerate(M.index)}
        for r, gene_name in enumerate(genes):
            g = gene_name
            sub = card[card.gene == g]
            idx = [(ai[a], b) for a, b in zip(sub.arm, sub.beta) if a in ai]
            if not idx:
                continue
            rows = np.array([i for i, _ in idx]); bet = np.array([b for _, b in idx], float)
            # ⛔ BUG FIXED 2026-08-01. `_resid_rows` returns an ALL-NaN row for any arm it cannot
            # residualise, and `bet @ Xr[rows]` then propagates that NaN across the WHOLE aggregate —
            # so ONE bad arm nulled the gene. It biased exactly the wrong way: the more arms a gene had,
            # the likelier it was dropped (PTEN's 88 arms → NaN, GATA3's 2 → fine), and gain RISES with
            # arm count, so the surviving subset under-represented the genes the weighting helps most.
            # Measured before the fix: 48.3% of t105 genes had NaN rho against 1.6% NaN n_arms.
            keep = np.isfinite(Xr[rows]).sum(1) > 10
            if not keep.any():
                continue
            rows, bet = rows[keep], bet[keep]
            # --- FAMILY RUNG: pool the gene's arms into seed families, then z-score, then apply β ---
            arms_g = [Xmat.index[i] for i in rows]
            fam = DB._ctx()["fam"].reindex(arms_g)
            # ⚠ pool from THIS COHORT's own arm matrix — `decoy_bench.pool_family` reads the TCGA X and
            # would return all-NaN for the prospective cohort, which carries its OWN CPTAC-2 miRNA-seq.
            # True RPM pool (log2(1 + Σ(2^x − 1))), never a mean — the canonical family collapse.
            lin = np.nan_to_num(np.power(2.0, Xmat.loc[arms_g].to_numpy(float)) - 1.0, nan=0.0)
            # ⛔ WAS: `bfam[f] = b` — the FIRST arm's β used as the family's. Arbitrary, because β is fit
            # PER ARM (`readouts.run(level="arm")`) and DIFFERS across arms in 99.8% of multi-arm cells
            # (MH-191). Picking one arm silently chose a winner among near-collinear same-seed arms.
            # ⭐ NOW: read the FAMILY-rung β from `gene_family_card.tsv` — the same estimator at
            # `level="family"`, i.e. the number that actually belongs on a family-pooled column
            # (MH-193). Falls back to the arms' MEAN, which is at least symmetric, if a cell is absent.
            fkeys, bfam, acc, seen = [], {}, {}, {}
            for a, b, row in zip(arms_g, bet, lin):
                f = str(fam.get(a))
                if f not in acc:
                    acc[f] = np.zeros(lin.shape[1]); fkeys.append(f); seen[f] = []
                acc[f] = acc[f] + row
                seen[f].append(b)
            fam_beta = _family_beta()
            for f in fkeys:
                bf = fam_beta.get((gene_name, f))
                bfam[f] = float(bf) if bf is not None and bf == bf else float(np.mean(seen[f]))
            fcols = fkeys
            bvec = np.array([bfam[f] for f in fcols], float)
            Fn = np.log2(1.0 + np.vstack([acc[f] for f in fcols]))   # family x sample
            Fz = (Fn - np.nanmean(Fn, 1, keepdims=True)) / (np.nanstd(Fn, 1, keepdims=True) + 1e-9)
            Fza = _resid_rows(Fz, Dcov)
            ok = np.isfinite(Fza).sum(1) > 10
            if not ok.any():
                continue
            Fza, Fzr, bvec = Fza[ok], Fz[ok], bvec[ok]
            agg_a = np.nansum(bvec[:, None] * Fza, axis=0); agg_r = np.nansum(bvec[:, None] * Fzr, axis=0)
            abu_a = np.nansum(Fza, axis=0); abu_r = np.nansum(Fzr, axis=0)   # UNWEIGHTED, same columns
            cols[f"{prefix}_n_arms"][r] = int(ok.sum())                      # FAMILIES actually used
            for tag in layers:
                i = gi[tag].get(g)
                if i is None:
                    continue
                cols[f"{prefix}_agg_rho_{tag}"][r] = _sp(agg_a, Lr[tag][i])
                cols[f"{prefix}_agg_rho_{tag}_raw"][r] = _sp(agg_r, Lraw[tag][i])
                cols[f"{prefix}_abund_rho_{tag}"][r] = _sp(abu_a, Lr[tag][i])

    def _blank(prefix, tags):
        k = [f"{prefix}_n_arms"] + [f"{prefix}_agg_rho_{t}" for t in tags] + \
            [f"{prefix}_agg_rho_{t}_raw" for t in tags] + [f"{prefix}_abund_rho_{t}" for t in tags]
        return {x: np.full(len(genes), np.nan) for x in k}

    for cohort, pref in COHORTS.items():
        d = OP._cptac(cohort)
        L = V.load_cptac_layers(cohort)
        samples = [s for s in d["arms"].columns if s in L["protein_z"].columns]
        Dcov = OP._cov(cohort, samples)
        layers = {t: L[k].reindex(columns=samples) for k, t in LAYERS.items()}
        cols = _blank(pref, LAYERS.values())
        _score(d["arms"][samples], layers, Dcov, pref, cols)
        for k, v in cols.items():
            res[k] = v
        rp, rr = res[f"{pref}_agg_rho_prot"], res[f"{pref}_agg_rho_prot_raw"]
        from mirna_hallmark.learned import retention as _RET
        res[f"{pref}_agg_ret_prot"] = _RET.ratio(rp, rr, gate=RHO_GATE,
                                                 name=f"{pref}_agg_ret_prot").to_numpy()
        # ⚠ a FLOOR, not a control — see the docstring
        # ⛔⛔ FIXED 2026-08-19 (column review unit 17): this was a BARE `<`, and `NaN < NaN` is **False**,
        # not NaN — so a gene whose protein coupling was never computed read as "beta does NOT beat
        # abundance". It was defined on all 1,420 genes while its two inputs cover 1,127 (prospective) and
        # 932 (tcga105); **ALL 293 / 488 unmeasurable rows read False, none True.** That silently deflated
        # the headline rate: prospective **25.1% -> 31.7%**, tcga105 **21.1% -> 32.1%** once the denominator
        # is the measurable set. Fourth instance of this failure class in one session (`echim_any`,
        # `fst_is_dominant_*`, `arb_max_identity`'s inert gate) — **a boolean built by comparing two
        # possibly-missing numbers MUST be masked, because `<` can never return "unknown".**
        # ✅ Verified: on rows where both inputs exist the logic itself was already correct (100% match).
        _ab = res[f"{pref}_abund_rho_prot"]
        res[f"{pref}_agg_beats_abund_prot"] = (rp < _ab).where(rp.notna() & _ab.notna())
        print(f"[cptac_gene] {cohort:12s} n={len(samples)} · "
              f"{int(res[f'{pref}_n_arms'].notna().sum()):,}/{len(genes):,} genes scored")

    # --- TCGA reference, IDENTICAL construction (n~1000) so the transfer is measured, not assumed ---
    dd = LD._load()
    Xt, Yt = dd["X"], dd["Y"]
    st = [s for s in Xt.columns if s in Yt.columns]
    Ct = CF.build_C("tcga", st).reindex(st).apply(pd.to_numeric, errors="coerce")
    Ct = Ct.fillna(Ct.median(numeric_only=True))
    Dt = np.column_stack([np.ones(len(st)), Ct.to_numpy(float)])
    cols = _blank("tcga", ["rna"])
    _score(Xt[st], {"rna": Yt.reindex(columns=st)}, Dt, "tcga", cols)
    for k, v in cols.items():
        res[k] = v
    print(f"[cptac_gene] tcga (ref)  n={len(st)} · {int(res['tcga_n_arms'].notna().sum()):,} genes scored")

    res.to_csv(DEST_GENE, sep="\t", index=False)
    print(f"-> {DEST_GENE}")
    return res
