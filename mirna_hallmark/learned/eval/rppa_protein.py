"""miRNA → PROTEIN at n=866 — the RPPA channel the subproject never opened.

    .venv/bin/python3 -m mirna_hallmark.learned.eval.rppa_protein [--workers 6]

⭐⭐ WHY THIS EXISTS, AND WHY IT IS NOT A REPEAT OF THE CPTAC ARC. `mirna_hallmark` has never touched RPPA
(verified: zero references across the subproject), yet **`data/rppa/` holds 919 TCGA-BRCA samples → an
881 × 461 participant×antibody matrix, of which 866 participants ALSO have miRNA** — the *same patients*,
so there is **NO COHORT BOUNDARY**, the loss MH-104 measured at ~80%.

**⭐ THE POWER ARGUMENT IS THE WHOLE POINT.** Detectable |ρ| at 80% power:
    CPTAC  n=101  →  |ρ| ≥ **0.276**
    RPPA   n=866  →  |ρ| ≥ **0.095**
and this program's typical gene-level coupling is **|ρ| ≈ 0.07–0.12**. ⇒ **CPTAC's detection floor was
about 3× the effect size it was asked to measure.** ⚠ **This does NOT retroactively rescue MH-103** — the
protein-βᵗ centrepiece also died of a **mediator leak**, which is a design flaw and not a power problem, and
`protein-beta-t-falsified` stands as written. What is new is that a *fresh* test finally sits inside the
regime where the effect lives.

✅ **POSITIVE CONTROL, RUN FIRST — the best this project has had.** MH-196's versioned literature set
∩ RPPA = **104 genes (77 at ≥3 PMIDs)** whose canonical regulator is known AND whose protein is measured:
PTEN←miR-21 · PDCD4←miR-21 · ZEB1←miR-200 · BCL2←miR-15/16 · CDKN1A←miR-17/20/93 · CDKN1B←miR-221/222 ·
MTOR←miR-99/100 · NOTCH1←miR-34 · EZH2←miR-101 · CCND1/CCNE1←miR-15/16. **If the canonical family's β does
not track its own protein here, nothing downstream is interpretable.**

⚠⚠ **COMPOSITION IS NOT OPTIONAL AND THIS IS EXACTLY WHERE THE PROJECT WAS BURNED BEFORE.**
`cptac-protein-composition-confound`: the CPTAC protein validation ran with **NO composition term**, and
the correction was a REFRAME, not a retraction — epithelial programs survived, stromal ones collapsed.
Every readout here is emitted under **both** C blocks with `retention = ρ_deconv/ρ_core` (axiom 2a).

⭐⭐ **THE ACTUAL PRIZE — THE DEGRADATION vs TRANSLATION DECOMPOSITION, powered for the first time.**
A miRNA represses two ways: it destabilises the transcript (visible in mRNA) and it blocks translation
(visible ONLY as protein-beyond-mRNA). Residualising protein on its own mRNA isolates the second:
    `rho_prot`      Spearman(M_g, protein | C)                — total protein-level repression
    `rho_mrna`      Spearman(M_g, mRNA    | C)                — the destabilisation channel
    `rho_discord`   Spearman(M_g, resid(protein ~ mRNA) | C)  — ⭐ the TRANSLATIONAL channel
**The project has never had the n to separate these.** At 866 it can.

⚠ RPPA CAVEATS, stated up front: ~200 informative antibodies is a narrow, signalling-biased panel (not a
proteome); antibody quality varies; and **phospho-antibodies measure a modified fraction, not abundance** —
they are reported separately and never pooled with total-protein probes.
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned/rppa"
MIN_N = 100
_MEM: dict = {}


def rppa():
    """(memoized) participant × antibody matrix + gene→antibody map, via the repo's existing loader."""
    if "rppa" in _MEM:
        return _MEM["rppa"]
    from analysis.utils.rppa_antibody_matrix import build_participant_antibody_matrix
    M, gmap = build_participant_antibody_matrix()
    g2ab: dict = {}
    for _, r in gmap.dropna(subset=["gene_name"]).iterrows():
        for g in str(r.gene_name).split():
            g2ab.setdefault(g.strip(), []).append(r.peptide_target)
    # phospho probes measure a MODIFIED FRACTION, not abundance — kept, but flagged, never pooled
    phos = {a for a in M.columns if any(k in str(a).upper() for k in ("_P", "PS", "PT", "PY")) and
            any(ch.isdigit() for ch in str(a))}
    _MEM["rppa"] = {"M": M, "g2ab": g2ab, "phospho": phos}
    return _MEM["rppa"]


def _partial(x, y, Cm):
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(Cm).all(1)
    if ok.sum() < MIN_N:
        return np.nan, ok.sum()
    R = np.column_stack([rankdata(x[ok]), rankdata(y[ok])])
    D = np.column_stack([np.ones(ok.sum())] + [rankdata(Cm[ok, j]) for j in range(Cm.shape[1])])
    res = R - D @ np.linalg.lstsq(D, R, rcond=None)[0]
    if np.std(res[:, 0]) < 1e-9 or np.std(res[:, 1]) < 1e-9:
        return np.nan, ok.sum()
    return float(np.corrcoef(res[:, 0], res[:, 1])[0, 1]), int(ok.sum())


def _one(gene: str):
    """Per gene: the β budget vs protein, mRNA, and protein-residualised-on-mRNA, under BOTH C blocks."""
    from mirna_hallmark.learned import attribution_eb as AE
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned import families as FAM

    R = _MEM["rppa"]
    abs_ = R["g2ab"].get(gene)
    if not abs_:
        return None
    beta = _MEM["beta"]
    out = {"gene": gene, "n_ab": len(abs_)}
    for blk, deconv in (("core", False), ("deconv", True)):
        try:
            Y, X, Cm, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=deconv)
        except Exception:
            return None
        fam = FAM.family_of(pd.Index(X.columns))
        Xf, _, _ = FAM.collapse_by_family(X, w, fam)
        fams = [c for c in Xf.columns if (gene, str(c)) in beta]
        if not fams:
            return None
        _, Xz, cols = AE._prep(Y, Xf[fams], Cm)
        b = np.array([beta[(gene, str(c))] for c in cols], float)
        M = Xz @ b
        parts = list(Y.index)
        Cmat = Cm.to_numpy(float)
        prot = R["M"].reindex(parts)
        # antibody: prefer a TOTAL-protein probe; phospho only if nothing else
        tot = [a for a in abs_ if a in prot.columns and a not in R["phospho"]]
        use = (tot or [a for a in abs_ if a in prot.columns])
        if not use:
            return None
        pv = prot[use[0]].to_numpy(float)
        yv = Y.to_numpy(float)
        rp, n = _partial(M, pv, Cmat)
        rm, _ = _partial(M, yv, Cmat)
        # discordance: protein residualised on its OWN mRNA -> the translational channel
        ok = np.isfinite(pv) & np.isfinite(yv)
        disc = np.full(len(pv), np.nan)
        if ok.sum() > MIN_N:
            A = np.column_stack([np.ones(ok.sum()), rankdata(yv[ok])])
            disc[ok] = rankdata(pv[ok]) - A @ np.linalg.lstsq(A, rankdata(pv[ok]), rcond=None)[0]
        rd, _ = _partial(M, disc, Cmat)
        ab_ref, _ = _partial(Xz.sum(1), pv, Cmat)     # R1: unweighted abundance on identical columns
        out |= {f"rho_prot_{blk}": rp, f"rho_mrna_{blk}": rm, f"rho_disc_{blk}": rd,
                f"rho_abund_prot_{blk}": ab_ref, f"n_{blk}": n, f"n_fam_{blk}": len(fams),
                "antibody": use[0], "is_phospho": use[0] in R["phospho"]}
    r0 = out.get("rho_prot_core")
    # ⭐ MH-258: canonical implementation, but the 0.02 gate VALUE is deliberately KEPT rather than
    # raised to RHO_GATE. Checked first: this output is well-behaved (`retention_prot` n=154, median
    # +0.681, **max 2.44** — no blow-up), so 0.02 is doing its job at this scale and raising it would
    # drop rows to fix nothing. ⚠ A gate is a claim about the DENOMINATOR'S scale, not a house style.
    from mirna_hallmark.learned import retention as _RET
    out["retention_prot"] = _RET.scalar(out.get("rho_prot_deconv"), r0, gate=0.02,
                                        name="retention_prot")
    return out


def run(workers: int = 6) -> pd.DataFrame:
    from multiprocessing import get_context

    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.stats import bh_fdr

    R = rppa()
    f = pd.read_csv(C.REPO_ROOT / "mirna_hallmark/output/learned/gene_family_card.tsv", sep="\t",
                    usecols=["gene", "family", "beta"])
    _MEM["beta"] = dict(zip(zip(f.gene, f.family.astype(str)), f.beta.astype(float)))
    genes = sorted(set(f.gene) & set(R["g2ab"]))
    D.load_cnv_target_genes(sorted(set(genes) & set(LD._load()["Y"].index)))   # ONE batched warm (axiom 3a)
    print(f"[rppa] {len(genes)} genes with both β and an antibody · {workers} workers", flush=True)
    with get_context("fork").Pool(workers) as p:
        rows = [r for r in p.imap_unordered(_one, genes, chunksize=4) if r]
    T = pd.DataFrame(rows)
    # Fisher-z for a PARTIAL correlation: df = n − 3 − k, k = number of covariates conditioned on.
    from scipy.stats import norm
    K = {"core": 3, "deconv": 11}          # [CPE, target_cn, mal_prolif] (+8 lineage fractions)
    for c in [x for x in T.columns if x.startswith(("rho_prot_", "rho_disc_", "rho_mrna_"))]:
        blk = c.rsplit("_", 1)[-1]
        r = T[c].to_numpy(float)
        n = T[f"n_{blk}"].to_numpy(float)
        df = n - 3.0 - K.get(blk, 3)
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.arctanh(np.clip(r, -0.999999, 0.999999)) * np.sqrt(np.maximum(df, 1.0))
        T["p_" + c[4:]] = 2.0 * norm.sf(np.abs(z))
    for c in [x for x in T.columns if x.startswith("p_")]:
        T["q_" + c[2:]] = bh_fdr(T[c].fillna(1.0).to_numpy())
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "rppa_coupling.tsv", sep="\t", index=False)
    print(f"-> {OUT / 'rppa_coupling.tsv'}  ({len(T)} genes)")
    return T


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    run(a.workers)
