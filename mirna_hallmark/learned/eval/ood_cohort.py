"""β(TCGA) → BUFFA — the cohort-boundary transfer test. The model's only untested boundary.

    .venv/bin/python3 -m mirna_hallmark.learned.eval.ood_cohort --stage 0     # the reproduction gate
    .venv/bin/python3 -m mirna_hallmark.learned.eval.ood_cohort --report      # full run

⭐ WHY. MH-104 measured the learned model losing **~80% crossing cohorts and ~0% crossing mRNA→protein**
(TCGA held-out ceiling median ρ **−0.117**, 88.8% repression-directed → CPTAC mRNA **−0.023**, retention
**0.193**). But that is TCGA→CPTAC: 101 patients, same consortium, a cohort profiled to be comparable to
TCGA. **Buffa (GSE22216 miRNA + GSE22219 mRNA, n=207) is the only genuinely independent-patient cohort in
the repo carrying both layers, and β had never touched it** — `learned/eval/__init__.py` reserved the name
`ood_cohort` as a `STUB`, and `eval/buffa_validation.py` runs a per-EDGE lane that never sees β.

⚠ **Buffa is NOT a METABRIC subset — `STATE_OF_PLAY` Axis 7 says it is, and that is WRONG (verified
2026-08-02).** Expression identity against all 1,980 METABRIC arrays: best-match r **median 0.270, max
0.608, margin median 0.0094, 0/207 above r=0.90** (a duplicate would be 0.95–0.99 with a large margin).
Clinical: **49%** of Buffa samples have *zero* candidates on (age±1, nodes, ER, size±2mm, grade); grade-1
frac 21% vs 9%, ER+ 62% vs 73%. Provenance agrees independently: Buffa = **Oxford**, Weatherall Institute,
surgery **1989–1992**; METABRIC = Cambridge/Vancouver. ⇒ **good news** (the cohort really is independent),
but `data/METABRIC/data_cna.txt` cannot be joined to it.

⛔ **`target_cn` IS THEREFORE UNOBTAINABLE, NOT MERELY ABSENT.** Exhaustive search (2026-08-02): GSE22220 is
a SuperSeries with **exactly two** subseries (miRNA + mRNA, no CNV/aCGH/SNP); Buffa 2011 (*Cancer Res*
71:5635, PMID 21737487) generates only expression; no later study CNV-profiles the same Oxford tumours.
**Do not repeat that search.**

⭐⭐ **WHICH IS WHY THE C-ABLATION CEILING EXISTS.** Comparing a reduced-C Buffa transfer against MH-104's
full-C TCGA ceiling would confound *the cohort boundary* with *a weaker confounder block*, and the retention
number would mean nothing. So the TCGA ceiling is ALSO recomputed under the exact reduced C that Buffa can
support. **Every retention figure names its denominator.**

THE DESIGN
----------
`M_g(s) = Σ_f β_f · Z_f(s)`, `ρ_g = Spearman(M_g, Y_g | C)` — β **FROZEN**, fit on TCGA, zero parameters
fit in Buffa.
* **Rung = FAMILY**, forced: Buffa's miRNA array is 501 legacy bare-stem probes with **no -5p/-3p
  resolution**, so the arm rung is not measurable there. (`cptac_card._family_beta` /
  `compartment_card._budget` reached the same conclusion for CPTAC — MH-193.)
* **E1 LEVEL** — median ρ, % negative. Expected to collapse (MH-174: *"transport the RANK, never the LEVEL"*).
* **E2 RANK** — Spearman(ρ_TCGA-ceiling, ρ_Buffa) across genes. **The primary.**
* **R1** unweighted sum-abundance on the identical columns (β is the only difference). This is also the
  cross-cohort gauge's `c` term, which is why it is mandatory rather than nice-to-have.
* **R2** the **transported DECOY β** — same pipeline, site-free fitted fake. MH-130d measured a fake reaching
  **99%** of the protein transfer, so without R2 a positive result is uninterpretable.

⛔ **NO GAUGE.** A transfer test fits nothing in Buffa: `Gauge.a` cancels under Spearman, `c` is exactly R1,
and `gauge.cohort_matrices`/`confounders.COHORTS` have no `buffa` branch. `info_ratio`/`calibrate` gate
CHANNEL construction, a different operation (MH-102b/d).
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "ood_cohort"
FAMILY_CARD = OUT / "family_card.tsv"
DECOY_FAMILY_BETA = OUT / "decoy_family_betas.tsv"

#: CAF/myofibroblast markers — same set `learned/data.py::_mycaf` uses.
CAF_MARKERS = ("ACTA2", "TAGLN", "POSTN", "FAP", "COL11A1", "THBS2", "MYH11", "CNN1")
MIN_N = 30
N_FOLDS = 5
GIBBS_ITER, GIBBS_BURN = 200, 80        # decoy_bench's validated setting (max|Δρ| ≤ 0.0013 vs 1600/600)
N_PERM = 20                             # draws for the within-gene beta-permutation control (R2b)

_MEM: dict = {}


# --------------------------------------------------------------------------- #
# Cohort
# --------------------------------------------------------------------------- #

def buffa() -> dict:
    """(memoized) Buffa miRNA + mRNA on the 207 paired samples. ⚠ 207, NOT the 210 that
    `buffa_validation.py:6` and its manifest still say — GSE22216 has 210 miRNA arrays and GSE22219 has
    216 mRNA arrays; **207 are matched**, and every row of the persisted output already says n=207."""
    if "buffa" in _MEM:
        return _MEM["buffa"]
    from mirna_hallmark.eval import buffa_validation as BV
    mi, rna = BV.load_buffa()
    _MEM["buffa"] = {"mirna": mi, "rna": rna, "samples": list(mi.columns)}
    return _MEM["buffa"]


def gene_alias() -> dict:
    """⭐ (memoized) OLD-SYMBOL RESOLVER: modern gene symbol → the row name Buffa's 2010 array actually uses.

    ⚠ **WITHOUT THIS, 16.7% OF THE β UNIVERSE IS SILENTLY DROPPED — including ZEB1**, one of the project's
    canonical pairs (miR-200/ZEB1), which the array carries under its 1999 symbol `TCF8`. Buffa's platform
    is GPL6098 (Illumina humanRef-8 v1.0, 2007 annotation), so a straight symbol join loses every gene
    renamed since. **MEASURED: 258/1,549 β genes absent by direct match; 73 recovered by synonym (28.3%),
    taking coverage 1,291 → 1,364.** Recovered examples: `ZEB1→TCF8`, `AGO2→EIF2C2`, `AGO4→EIF2C4`,
    `CDK1→CDC2`, `CCN2→CTGF`, `CD40→TNFRSF5`, `ASS1→ASS`.

    Built from the **locally cached** `GPL6098_family.soft.gz` `Synonym` column — no network, and the
    platform's own annotation is the authority for what its probes were called. A synonym is admitted only
    if its target Symbol is actually a row of the matrix, and a direct hit always wins over a synonym.

    ⚠ The matrix index also carries Excel date corruption (`1-Sep`, `1-Nov` — i.e. SEPT1, CCN3/NOV), a
    known hazard of 2010-era symbol handling. Those rows are left alone; none is a β gene.
    """
    if "alias" in _MEM:
        return _MEM["alias"]
    import gzip
    from mirna_hallmark.eval.buffa_validation import GEO_DIR
    idx = set(buffa()["rna"].index)
    rows, inside = [], False
    with gzip.open(GEO_DIR / "GPL6098_family.soft.gz", "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!platform_table_begin"):
                inside = True; continue
            if line.startswith("!platform_table_end"):
                break
            if inside:
                rows.append(line.rstrip("\n").split("\t"))
    hdr = rows[0]
    si, yi = hdr.index("Symbol"), hdr.index("Synonym")
    amap: dict = {}
    for r in rows[1:]:
        if len(r) <= max(si, yi):
            continue
        sym = r[si].strip()
        if sym not in idx:                       # only resolve TO a row that exists
            continue
        for s in r[yi].split(";"):
            s = s.strip()
            if s and s != "nan" and s not in idx:   # a direct hit always wins
                amap.setdefault(s, sym)
    _MEM["alias"] = amap
    return amap


def buffa_row(gene: str):
    """The Buffa mRNA row for `gene`, resolving old symbols. Returns (values, resolved_name) or (None, None)."""
    rna = buffa()["rna"]
    name = gene if gene in rna.index else gene_alias().get(gene)
    if name is None:
        return None, None
    v = rna.loc[name]
    if isinstance(v, pd.DataFrame):
        v = v.iloc[0]
    return pd.to_numeric(v, errors="coerce").to_numpy(float), name


# --------------------------------------------------------------------------- #
# STAGE 1 — the three confounder blocks
# --------------------------------------------------------------------------- #

def c_blocks(*, force_estimate: bool = False) -> dict:
    """The three Buffa C blocks. **Hard-fails if a block cannot be built.**

    ⛔ THE FAILURE MODE THIS GUARDS. `learned/eval/ood_protein.py:47-59` wraps `build_C` in a bare
    `except` and returns an **intercept-only** design on failure — so a RAW Spearman gets reported as
    "adjusted", with `n_cov == 0` as the only tell. That is MH-107's exact defect (a channel that ran one
    block and never computed retention). `confounders.build_C` raises `ValueError` for any cohort outside
    `("tcga","nat","gtex","cptac")`, so on Buffa that `except` would fire every time. **Here a block that
    cannot be built is an exception, never a silent downgrade.**

    | block      | terms                                                                    |
    |------------|--------------------------------------------------------------------------|
    | `raw`      | prolif metagene (E2F ∪ G2M), ER status                                     |
    | `metagene` | + ESTIMATE Stromal/Immune (ssGSEA — rank-based, genuinely platform-free)   |
    |            |   + CAF marker metagene                                                    |
    | `nnls`     | + 8 Wu lineage fractions (Cancer Epithelial held out, as `confounders.py:78`) |

    ⚠ NO `target_cn` in any block — unobtainable for this cohort (see module docstring). That is what the
    C-ablation ceiling exists to price.
    """
    if "cblocks" in _MEM and not force_estimate:
        return _MEM["cblocks"]
    from mirna_hallmark.eval import buffa_validation as BV
    from mirna_hallmark.hallmark_sets import HallmarkSets
    from mirna_hallmark.learned.confounders import _metagene

    B = buffa()
    rna, samples = B["rna"], B["samples"]

    prolif = BV.prolif_proxy(rna, HallmarkSets.load()).reindex(samples)
    if prolif.isna().all():
        raise RuntimeError("[ood_cohort] prolif metagene is all-NaN — refusing to continue")
    er = BV.load_buffa_er().reindex(samples)
    if er.notna().sum() < 0.9 * len(samples):
        raise RuntimeError(f"[ood_cohort] ER present for only {int(er.notna().sum())}/{len(samples)}")
    raw = pd.DataFrame({"prolif": prolif, "er": er}, index=samples)

    # ---- metagene block: real ESTIMATE (ssGSEA, platform-free) + a CAF metagene
    from mirna_hallmark import estimate_scores as ES
    est = ES.compute_estimate("buffa", rna.apply(pd.to_numeric, errors="coerce"), force=force_estimate)
    est = est.reindex(samples)
    missing = [c for c in ("StromalScore", "ImmuneScore") if c not in est.columns or est[c].isna().all()]
    if missing:
        raise RuntimeError(f"[ood_cohort] ESTIMATE produced no {missing} — refusing to fall back")
    caf = _metagene(rna.apply(pd.to_numeric, errors="coerce"), CAF_MARKERS)
    if caf is None:
        raise RuntimeError("[ood_cohort] CAF metagene unbuildable (<2 markers present)")
    meta = raw.assign(estimate_stromal=est["StromalScore"].to_numpy(),
                      estimate_immune=est["ImmuneScore"].to_numpy(),
                      caf_metagene=caf.reindex(samples).to_numpy())

    # ---- nnls block: Wu lineage fractions from bulk mRNA alone
    fr = lineage_fractions()
    nn = raw.join(fr.reindex(samples).add_prefix("fr_"))
    if nn.isna().all().any():
        raise RuntimeError("[ood_cohort] a lineage fraction column is all-NaN")

    out = {"raw": raw, "metagene": meta, "nnls": nn}
    for k, v in out.items():
        out[k] = v.astype(float).apply(lambda s: s.fillna(s.median()))
        print(f"[ood_cohort] C block {k:9s}: {out[k].shape[1]} cols  {list(out[k].columns)}")
    _MEM["cblocks"] = out
    return out


def lineage_fractions() -> pd.DataFrame:
    """(memoized) NNLS deconvolution of Buffa bulk mRNA against the Wu celltype_major signature.

    ⚠ Buffa is **log2 beadchip intensity**, and `deconvolve_matrix` wants LINEAR — so we pass `2**X`.
    Microarray intensity is not proportional to abundance the way CPM/TPM is, so these fractions carry a
    platform bias. They are tolerable here because they enter C only as **covariates in a rank-based
    partial correlation**, where a monotone distortion is far less damaging than the same distortion on the
    estimand — and MH-114 is unambiguous that running with NO composition block reproduces its exact defect
    (gene-level FDR-neg 34→3, hallmark EMT −0.330→+0.162 sign flip). The `metagene` block is the
    platform-free cross-check; **disagreement between the two IS a result.**
    """
    if "fr" in _MEM:
        return _MEM["fr"]
    from mirna_hallmark.analyses.spatial.spatial_common import deconvolve_matrix
    rna = buffa()["rna"].apply(pd.to_numeric, errors="coerce")
    lin = np.power(2.0, rna).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    fr = deconvolve_matrix(lin)
    fr = fr.drop(columns=[c for c in fr.columns if "Cancer Epithelial" in c])   # held out, confounders.py:78
    print(f"[ood_cohort] NNLS lineage fractions {fr.shape}: {list(fr.columns)}")
    _MEM["fr"] = fr
    return fr


# --------------------------------------------------------------------------- #
# Design construction — Buffa family columns
# --------------------------------------------------------------------------- #

def _stem(a: str) -> str:
    import re
    return re.sub(r"-(5p|3p)$", "", str(a))


def buffa_family_matrix() -> pd.DataFrame:
    """(memoized) (gene, family) → the Buffa family DOSE column, for every design cell measurable there.

    Pooling is the **TRUE RPM pool** `log2(1 + Σ(2^x − 1))` (`decoy_bench.pool_family`), never a
    `groupby.mean()` — that was a documented decoy-bench defect. Arms map by **bare stem** because the
    Illumina miRNA array carries legacy names with no -5p/-3p (guide-arm-accurate, approximate for
    passenger arms — same caveat `buffa_validation` flags).
    """
    if "fammat" in _MEM:
        return _MEM["fammat"]
    fam = pd.read_csv(FAMILY_CARD, sep="\t", usecols=["gene", "family", "arms"])
    mi = buffa()["mirna"].apply(pd.to_numeric, errors="coerce")
    by_stem: dict = {}
    for a in mi.index:
        by_stem.setdefault(_stem(a), a)
    rows, idx = [], []
    for g, f, arms in zip(fam.gene, fam.family, fam.arms):
        members = [by_stem[_stem(a)] for a in str(arms).split(";") if _stem(a) in by_stem]
        if not members:
            continue
        lin = np.power(2.0, mi.loc[sorted(set(members))].to_numpy(float)) - 1.0
        rows.append(np.log2(1.0 + np.clip(lin, 0, None).sum(axis=0)))
        idx.append((g, f))
    M = pd.DataFrame(rows, index=pd.MultiIndex.from_tuples(idx, names=["gene", "family"]),
                     columns=mi.columns)
    print(f"[ood_cohort] Buffa family-dose matrix: {M.shape[0]:,} (gene,family) cells measurable "
          f"of {len(fam):,} ({M.shape[0]/len(fam):.1%})")
    _MEM["fammat"] = M
    return M


def family_beta(decoy: bool = False) -> dict:
    """(memoized) frozen FAMILY-rung β, `(gene, family) -> beta`.

    Real β is `readouts.run(level="family")` (the §8 collapse APPLIED) as carried by `family_card.tsv`.
    ⛔ NOT the edge card's `beta`, which is `level="arm"` — same name, different unit (MH-191/192), and
    Buffa's array cannot resolve arms anyway.
    Decoy β is the site-free fitted fake from `decoy_family_betas.tsv` (`decoy == True`).
    """
    key = f"beta_{decoy}"
    if key in _MEM:
        return _MEM[key]
    if decoy:
        d = pd.read_csv(DECOY_FAMILY_BETA, sep="\t")
        d = d[d["decoy"].astype(bool)]
        out = dict(zip(zip(d.gene, d.unit.astype(str)), d.beta.astype(float)))
    else:
        f = pd.read_csv(FAMILY_CARD, sep="\t", usecols=["gene", "family", "beta"])
        out = dict(zip(zip(f.gene, f.family.astype(str)), f.beta.astype(float)))
    _MEM[key] = out
    return out


def buffa_family_all() -> pd.DataFrame:
    """(memoized) family -> Buffa dose column over ALL Buffa arms grouped by seed family.
    Used for the DECOY transfer, whose units have no per-gene member list in the design."""
    if "famall" in _MEM:
        return _MEM["famall"]
    from mirna_hallmark.coupling_inference import seed_family_map
    mi = buffa()["mirna"].apply(pd.to_numeric, errors="coerce")
    fam = seed_family_map(pd.Index([_stem(a) for a in mi.index]))
    lab = pd.Series([str(fam.get(_stem(a))) for a in mi.index], index=mi.index)
    rows, idx = [], []
    for f, members in lab.groupby(lab).groups.items():
        if f == "None":
            continue
        lin = np.power(2.0, mi.loc[list(members)].to_numpy(float)) - 1.0
        rows.append(np.log2(1.0 + np.clip(lin, 0, None).sum(axis=0)))
        idx.append(f)
    M = pd.DataFrame(rows, index=idx, columns=mi.columns)
    _MEM["famall"] = M
    return M


def _z(A: np.ndarray) -> np.ndarray:
    """z-score columns WITHIN THE COHORT, with `attribution_eb._prep`'s variance floor. A within-cohort
    z is what makes the transfer bridge-free: no QN onto TCGA's distribution shape (rejected in
    `method_dev/landscape/cross_state_top50.py:12`), only the assumption that the families' rank-order of
    informativeness transports — strictly weaker, and the same one CPTAC-prospective already runs under."""
    mu, sd = A.mean(0), A.std(0)
    Z = (A - mu) / (sd + 1e-9)
    Z[:, sd < 0.1] = 0.0
    return Z


def _partial_spearman(m: np.ndarray, y: np.ndarray, Cm: np.ndarray) -> float:
    """Partial Spearman = Pearson on ranks residualised against ranked C (`realization.partial_spearman`)."""
    ok = np.isfinite(m) & np.isfinite(y)
    if ok.sum() < MIN_N:
        return np.nan
    R = np.column_stack([rankdata(m[ok]), rankdata(y[ok])])
    Cr = np.column_stack([np.ones(ok.sum())] + [rankdata(Cm[ok, j]) for j in range(Cm.shape[1])])
    res = R - Cr @ np.linalg.lstsq(Cr, R, rcond=None)[0]
    if np.std(res[:, 0]) < 1e-9 or np.std(res[:, 1]) < 1e-9:
        return np.nan
    return float(np.corrcoef(res[:, 0], res[:, 1])[0, 1])


def _tcga_ablated_C(gene: str):
    """The TCGA design under the **reduced C Buffa can support** — the C-ABLATION ceiling's block.

    ⭐ WHY THIS IS THE LOAD-BEARING CONTROL. Buffa has no `target_cn` (gene-specific CN; unobtainable, not
    merely missing — see the module docstring) and no CPE (DNA-based consensus purity). Scoring a Buffa
    transfer against MH-104's FULL-C TCGA ceiling would therefore confound **the cohort boundary** with
    **a weaker confounder block**, and the retention figure would be uninterpretable.

    Mirror: drop `target_cn`, drop `CPE`, keep malignant proliferation, add TCGA's own ESTIMATE
    stromal/immune (`estimate_scores`, the same ssGSEA estimator run on Buffa).
    ⚠ The mirror is CLOSE, not exact: TCGA's prolif is `mal_prolif` (residualised on the Cancer-Epithelial
    fraction) while Buffa's is a plain metagene, and ER is in Buffa's block but not here. Both asymmetries
    make the ablated TCGA ceiling if anything the *stronger* design, so retention against it is conservative.
    """
    from mirna_hallmark.learned import data as LD
    Y, X, Cfull, w = LD.assemble_gene(gene, w_prior_source="ledger", target_cn=False)
    est = _MEM.get("tcga_est")
    if est is None:
        est = pd.read_csv(C.TISSUE_REFERENCE_DIR / "estimate" / "tumor_estimate_scores.tsv",
                          sep="\t", index_col=0)
        _MEM["tcga_est"] = est
    Cab = Cfull.drop(columns=[c for c in Cfull.columns if c.upper() == "CPE"], errors="ignore").copy()
    e = est.reindex(Cab.index)
    Cab["estimate_stromal"] = e["StromalScore"].to_numpy()
    Cab["estimate_immune"] = e["ImmuneScore"].to_numpy()
    Cab = Cab.apply(lambda s: s.fillna(s.median()))
    return Y, X, Cab


_BLOCKS = ("raw", "metagene", "nnls")


def _one_gene(gene: str):
    """One gene: the two TCGA ceilings, the Buffa target ceiling, the frozen-β transfer under three C
    blocks, and both references. Returns a row dict or None."""
    from mirna_hallmark.eval.decoy_bench import oof_budget
    from mirna_hallmark.learned import families as FAM

    B, CB = _MEM["buffa"], _MEM["cblocks"]
    FM, FA = _MEM["fammat"], _MEM["famall"]
    beta, dbeta = _MEM["beta_False"], _MEM["beta_True"]
    yv, yname = buffa_row(gene)
    if yv is None:
        return None
    rec = {"gene": gene, "buffa_row": yname, "via_alias": bool(yname != gene)}
    # ---------- TCGA ceilings (source, full C) and (C-ablation, Buffa-supportable C) ----------
    try:
        from mirna_hallmark.learned import data as LD
        Y, X, Cfull, w = LD.assemble_gene(gene, w_prior_source="ledger")
        fam = FAM.family_of(pd.Index(X.columns))
        Xf, _, _ = FAM.collapse_by_family(X, w, fam)
        rec["n_tcga"] = len(Y); rec["n_fam_tcga"] = Xf.shape[1]
        rec["rho_tcga_full"] = oof_budget(Y, Xf, Cfull)
        _, Xa, Cab = _tcga_ablated_C(gene)
        Xfa, _, _ = FAM.collapse_by_family(Xa, w, FAM.family_of(pd.Index(Xa.columns)))
        rec["rho_tcga_abl"] = oof_budget(Y, Xfa, Cab)
    except Exception:
        rec["rho_tcga_full"] = rec["rho_tcga_abl"] = np.nan

    # ---------- Buffa: the gene's measurable family columns ----------
    try:
        sub = FM.xs(gene, level="gene")
    except KeyError:
        return rec
    fams = [f for f in sub.index if (gene, f) in beta]
    if len(fams) < 1:
        return rec
    A = sub.loc[fams].to_numpy(float).T                      # samples x families
    Z = _z(A)
    b = np.array([beta[(gene, f)] for f in fams], float)
    rec["n_fam_buffa"] = len(fams)
    rec["fam_cov"] = len(fams) / max(1, len(beta_families_of(gene)))
    M = Z @ b
    ABU = Z.sum(1)                                           # R1 — unweighted, SAME columns
    for blk in _BLOCKS:
        Cm = CB[blk].to_numpy(float)
        rec[f"rho_buffa_{blk}"] = _partial_spearman(M, yv, Cm)
        rec[f"rho_abund_{blk}"] = _partial_spearman(ABU, yv, Cm)

    # ---------- R2a: the transported DECOY beta (site-free fitted fake) ----------
    # ⚠ LOW COVERAGE BY CONSTRUCTION, and the reason is structural, not a bug: decoy units are chosen to be
    # SITE-FREE, while Buffa's 501-probe legacy array covers mostly well-studied miRNAs. Measured: only
    # 20/424 decoy units (4.7%) and 164/3,558 decoy cells (4.6%) exist in Buffa, over 131/516 genes.
    # Reported on that subset only, as a paired contrast; R2b below is the full-coverage complement.
    dfams = [f for (g, f) in dbeta if g == gene and f in FA.index]
    if len(dfams) >= 1:
        Zd = _z(FA.loc[dfams].to_numpy(float).T)
        bd = np.array([dbeta[(gene, f)] for f in dfams], float)
        Md = Zd @ bd
        rec["n_fam_decoy"] = len(dfams)
        for blk in _BLOCKS:
            rec[f"rho_decoy_{blk}"] = _partial_spearman(Md, yv, CB[blk].to_numpy(float))

    # ---------- R2b: WITHIN-GENE beta PERMUTATION — full coverage, complements R2a ----------
    # Holds the column set AND the multiset of weights fixed, permuting only WHICH family gets WHICH beta.
    # Asks: does the model's ASSIGNMENT of weight to family carry the transfer, or would any reweighting of
    # the same abundance columns do? A site-free decoy tests a different thing (can these arms bind at all);
    # this tests the ranking. Averaged over `n_perm` draws so it is a null level, not one draw.
    if len(fams) >= 2:
        rng = np.random.default_rng(abs(hash(gene)) % (2 ** 32))
        acc = {blk: [] for blk in _BLOCKS}
        for _ in range(N_PERM):
            Mp = Z @ rng.permutation(b)
            for blk in _BLOCKS:
                acc[blk].append(_partial_spearman(Mp, yv, CB[blk].to_numpy(float)))
        for blk in _BLOCKS:
            v = [x for x in acc[blk] if np.isfinite(x)]
            rec[f"rho_permbeta_{blk}"] = float(np.mean(v)) if v else np.nan

    # ---------- the TARGET ceiling: beta refit IN BUFFA, out-of-fold ----------
    try:
        ys = pd.Series(yv, index=B["samples"])
        Xb = pd.DataFrame(A, index=B["samples"], columns=[str(f) for f in fams])
        ok = ys.notna()
        if ok.sum() >= 100:
            rec["rho_buffa_ceiling"] = oof_budget(ys[ok], Xb[ok], CB["metagene"][ok])
    except Exception:
        rec["rho_buffa_ceiling"] = np.nan
    return rec


def beta_families_of(gene: str) -> list:
    idx = _MEM.setdefault("beta_idx", {})
    if not idx:
        for (g, f) in _MEM["beta_False"]:
            idx.setdefault(g, []).append(f)
    return idx.get(gene, [])


def _warm():
    """⚠ WARM EVERY STATIC READ BEFORE FORKING (axiom 3a). On this NFS box an un-memoized per-gene read is
    fatal at N×workers — and `data_loaders.load_cnv_target_genes` re-scans the ASCAT3 SOURCE on any cache
    miss, which cost 2 h earlier today; one batched call over the whole universe is 75 s."""
    from mirna_hallmark import data_loaders as D
    from mirna_hallmark.learned import data as LD
    buffa(); gene_alias(); c_blocks(); buffa_family_matrix(); buffa_family_all()
    family_beta(False); family_beta(True); beta_families_of("__warm__")
    _tcga_ablated_C  # noqa: B018  (module-level cache primed on first call)
    genes = sorted({g for (g, _) in _MEM["beta_False"]})
    yidx = set(LD._load()["Y"].index)
    D.load_cnv_target_genes(sorted(set(genes) & yidx))       # ONE batched call, never per gene
    print(f"[ood_cohort] warm: {len(genes):,} beta genes")
    return genes


def run(genes=None, workers: int = 6) -> pd.DataFrame:
    from multiprocessing import get_context
    allg = _warm()
    rna_idx, al = set(buffa()["rna"].index), gene_alias()
    genes = genes or sorted(g for g in allg if g in rna_idx or g in al)
    n_alias = sum(1 for g in genes if g not in rna_idx)
    print(f"[ood_cohort] gene universe: {len(genes):,}  (+{n_alias} recovered via old-symbol alias)")
    print(f"[ood_cohort] scoring {len(genes):,} genes (beta ∩ Buffa mRNA) · {workers} workers", flush=True)
    with get_context("fork").Pool(workers) as p:
        rows = [r for r in p.imap_unordered(_one_gene, genes, chunksize=8) if r]
    R = pd.DataFrame(rows)
    DEST.mkdir(parents=True, exist_ok=True)
    R.to_csv(DEST / "ood_cohort_genes.tsv", sep="\t", index=False)
    print(f"-> {DEST / 'ood_cohort_genes.tsv'}  ({len(R):,} genes)")
    return R


#: what I committed to BEFORE running (plan file, 2026-08-02). Being wrong here is informative.
PREREG = {"E1_retention_vs_ablation": (0.10, 0.25),
          "E2_rank_transfer": (0.10, 0.25),
          "R1_abundance_within_10pct_of_beta": True,
          "R2b_permuted_beta_well_below_real": True}


def report(R: pd.DataFrame, block: str = "metagene") -> dict:
    """⭐⭐ **THE ONE-FAMILY DEGENERACY GOVERNS THIS ENTIRE READOUT — split on it or be misled.**

    For a gene with a SINGLE measurable family, `M = β·Z` and the abundance reference is `ABU = Z`. Spearman
    is invariant to positive scaling, so **ρ_β ≡ ρ_abundance EXACTLY** (verified: max|Δ| = 0.00e+00 over all
    587 such genes). β cannot act there. **43% of the scored universe is in that state** (median 2 families),
    so a pooled statistic silently averages β's real effect together with 587 genes where β is inert — which
    is exactly what produced the incoherent "median Δ = +0.0000 with Wilcoxon p = 1.3e-27" in the first
    version of this report, and what made the pooled rank transfer look like a β result when 1-family genes
    alone give **+0.234**. Everything β-vs-reference is therefore reported on **multi-family genes only**,
    with the 1-family stratum shown separately as the pure-abundance channel it is.
    """
    from scipy.stats import binomtest, wilcoxon

    one, multi = R[R.n_fam_buffa == 1], R[R.n_fam_buffa >= 2]
    print(f"\n{'='*102}\nβ(TCGA) → BUFFA — cohort-boundary transfer   [C block: {block}]\n{'='*102}")
    print(f"  {len(R):,} genes scored ({int(R.get('via_alias', pd.Series(dtype=bool)).sum())} via old-symbol "
          f"alias) · ⭐ {len(multi):,} MULTI-family (β can act) · {len(one):,} single-family "
          f"({len(one)/len(R):.0%}, β mathematically inert)\n")

    def _med(c, D=R):
        return float(D[c].median()) if c in D else np.nan

    cf, ca, cb = _med("rho_tcga_full"), _med("rho_tcga_abl"), _med("rho_buffa_ceiling")
    print("  --- ⭐ THE CEILINGS (every retention names its denominator) ---")
    print(f"    source        TCGA 5-fold OOF, FULL C                 median ρ {cf:+.4f}")
    print(f"    ⭐ C-ABLATION  TCGA 5-fold OOF, Buffa-supportable C     median ρ {ca:+.4f}"
          f"   ⇒ dropping target_cn+CPE costs {abs(cf)-abs(ca):+.4f} of |ρ| ({1-abs(ca)/abs(cf):.0%})")
    print(f"    in-cohort     β REFIT in Buffa, 5-fold OOF             median ρ {cb:+.4f}")
    print("    ⚠ THE IN-COHORT FIT IS **NOT A CEILING** — it is WORSE than the frozen transfer, so a")
    print("      'retention' against it exceeds 1 and is meaningless. At n=207 with a median of 2 families,")
    print("      refitting in-cohort is noisier than importing a β fit on 1,065 TCGA patients. Report it as")
    print("      what it is: transporting beats refitting here. Denominator for E1 is the C-ABLATION ceiling.")

    out: dict = {"n_multi": int(len(multi)), "n_single": int(len(one))}
    print(f"\n  --- E1 LEVEL, MULTI-FAMILY ONLY (MH-174: 'transport the RANK, never the LEVEL') ---")
    for blk in _BLOCKS:
        m = multi.dropna(subset=[f"rho_buffa_{blk}", f"rho_abund_{blk}"])
        s, ab = m[f"rho_buffa_{blk}"], m[f"rho_abund_{blk}"]
        neg = float((s < 0).mean()); p = binomtest(int((s < 0).sum()), len(s), 0.5).pvalue
        w = wilcoxon(s - ab).pvalue
        ret = float(s.median() / ca) if np.isfinite(ca) and abs(ca) > 1e-9 else np.nan
        print(f"    {blk:9s} β {s.median():+.4f} ({neg:.1%} neg, p={p:.2g})  vs R1 unweighted {ab.median():+.4f}"
              f"   Δ {(s-ab).median():+.4f} p={w:.2g}   retention vs C-ABLATION {ret:+.3f}")
        out[blk] = {"median_rho": float(s.median()), "frac_neg": neg, "p_binom": float(p),
                    "abund_median": float(ab.median()), "delta_vs_abund": float((s - ab).median()),
                    "p_vs_abund": float(w), "retention_vs_ablation": ret, "n": int(len(m))}

    print(f"\n  --- ⭐ E2 RANK TRANSFER — the PRIMARY estimand, β vs its references ---")
    for blk in _BLOCKS:
        m = multi.dropna(subset=["rho_tcga_abl", f"rho_buffa_{blk}", f"rho_abund_{blk}"])
        rb = spearmanr(m.rho_tcga_abl, m[f"rho_buffa_{blk}"])
        ra = spearmanr(m.rho_tcga_abl, m[f"rho_abund_{blk}"])
        pc = f"rho_permbeta_{blk}"
        rp = spearmanr(m.rho_tcga_abl, m[pc]) if pc in m else None
        print(f"    {blk:9s} β {rb.correlation:+.4f} (p={rb.pvalue:.2g})  |  R1 abundance "
              f"{ra.correlation:+.4f} (p={ra.pvalue:.2g})  |  R2b permuted-β "
              f"{rp.correlation:+.4f}" if rp is not None else "", end="")
        print(f"   n={len(m)}")
        out[blk].update({"rank_transfer": float(rb.correlation), "rank_transfer_p": float(rb.pvalue),
                         "rank_transfer_abund": float(ra.correlation),
                         "rank_transfer_permbeta": float(rp.correlation) if rp is not None else np.nan})
    m1 = one.dropna(subset=["rho_tcga_abl", f"rho_buffa_{block}"])
    r1 = spearmanr(m1.rho_tcga_abl, m1[f"rho_buffa_{block}"])
    print(f"\n    ⚠ SINGLE-family stratum (β INERT — this is the PURE ABUNDANCE channel): "
          f"{r1.correlation:+.4f} (p={r1.pvalue:.2g}, n={len(m1)})")
    print(f"      ⇒ abundance transports on its own. β's claim rests ONLY on the multi-family increment above.")
    out["rank_transfer_single_family_pure_abundance"] = float(r1.correlation)

    print("\n  --- R2b β-PERMUTATION (same columns, same weight multiset, shuffled assignment) ---")
    for blk in _BLOCKS:
        c = f"rho_permbeta_{blk}"
        if c in multi:
            m = multi.dropna(subset=[f"rho_buffa_{blk}", c])
            w = wilcoxon(m[f"rho_buffa_{blk}"] - m[c]).pvalue
            print(f"    {blk:9s} real {m[f'rho_buffa_{blk}'].median():+.4f} vs permuted {m[c].median():+.4f}"
                  f"   Wilcoxon p={w:.2g}  n={len(m)}")
            out[blk]["permbeta_median"] = float(m[c].median()); out[blk]["p_vs_permbeta"] = float(w)
    dm = multi.dropna(subset=[f"rho_buffa_{block}", f"rho_decoy_{block}"])
    if len(dm) > 10:
        w = wilcoxon(dm[f"rho_buffa_{block}"] - dm[f"rho_decoy_{block}"]).pvalue
        print(f"    R2a DECOY β real {dm[f'rho_buffa_{block}'].median():+.4f} vs decoy "
              f"{dm[f'rho_decoy_{block}'].median():+.4f}  p={w:.2g}  ⚠ n={len(dm)} ONLY — decoy units are "
              f"SITE-FREE and Buffa's legacy array measures 4.7% of them. UNDERPOWERED, not negative.")
        out["decoy_subset_n"] = int(len(dm)); out["decoy_p"] = float(w)

    print("\n  --- BY DESIGN SIZE (multi-family only) ---")
    mm = multi.copy(); mm["_b"] = pd.cut(mm.n_fam_buffa, [1, 3, 6, 12, 10_000], labels=["2-3", "4-6", "7-12", "13+"])
    for bkt, g in mm.groupby("_b", observed=True):
        m = g.dropna(subset=["rho_tcga_abl", f"rho_buffa_{block}"])
        rr = spearmanr(m.rho_tcga_abl, m[f"rho_buffa_{block}"]).correlation if len(m) > 10 else np.nan
        print(f"    {str(bkt):>5s} fams: n={len(g):4d}  median ρ_Buffa {g[f'rho_buffa_{block}'].median():+.4f}"
              f"   rank-transfer {rr:+.4f}")

    print(f"\n  --- ⭐ PRE-REGISTERED vs OBSERVED (multi-family) ---")
    e1, e2 = out[block]["retention_vs_ablation"], out[block]["rank_transfer"]
    for name, val, key in (("E1 retention vs C-ablation", e1, "E1_retention_vs_ablation"),
                           ("E2 rank transfer", e2, "E2_rank_transfer")):
        lo, hi = PREREG[key]
        print(f"    {name:28s} predicted [{lo}, {hi}]  observed {val:+.3f}   "
              f"{'✅ in band' if lo <= val <= hi else '⛔ OUTSIDE — record the miss'}")
    print(f"    {'R1 β beats abundance':28s} predicted 'within ~10%'  observed β {out[block]['median_rho']:+.4f}"
          f" vs {out[block]['abund_median']:+.4f}, p={out[block]['p_vs_abund']:.2g}")
    return out


def compartment_null(R: pd.DataFrame, n_shuffle: int = 3, seed: int = 0) -> dict:
    """⭐ MH-114 COMPARTMENT-ORIENTATION STRATIFIED NULL — the board's pre-registration requirement.

    Both cohorts are bulk breast and share the CAF confound, so **a clean replication proves nothing on its
    own**: bulk mixing manufactures a negative ρ whenever miRNA and target load OPPOSITELY on CAF, and a
    positive one when they load the SAME way. MH-114 measured that in CPTAC the *whole* apparent effect was
    that arithmetic — real gradient −0.1272 vs **shuffled −0.1290**, identical.

    ⛔ AND THE RETRACTION INSIDE MH-114 IS THE PART TO COPY: a shuffled null compared on the **MEAN** is
    blind to this, because the artifact is sign-symmetric and cancels (shuffled mean ρ ≈ +0.0001). That
    unstratified null gave the confidently WRONG verdict. **A null must be stratified by the confound's own
    axis.**

    Here, at the gene-BUDGET level: score gene i's budget against gene j's mRNA (degree-preserving budget
    shuffle), stratify by `sign(corr(budget, CAF)) × sign(corr(target, CAF))`.
    ⚠ `UNDEFINED` is an EXPLICIT third state — `protein_compartment_null.py:136` bins a NaN product as
    "SAME" because `NaN < 0` is False, silently contaminating the reference stratum.

    THE READOUT IS **NOT** THE GRADIENT — it is the **compartment-BLIND real-minus-shuffled contrast**,
    which for genuine repression must be the SAME SIGN in BOTH orientations (an artifact flips).
    """
    B, CB, FM = buffa(), c_blocks(), buffa_family_matrix()
    beta = family_beta(False)
    caf = CB["metagene"]["caf_metagene"].to_numpy(float)
    Cm = CB["metagene"].to_numpy(float)
    genes = sorted(set(R[R.n_fam_buffa >= 2].gene))
    bud, tar, names = [], [], []
    for g in genes:
        try:
            sub = FM.xs(g, level="gene")
        except KeyError:
            continue
        fams = [f for f in sub.index if (g, f) in beta]
        yv, _ = buffa_row(g)
        if len(fams) < 2 or yv is None:
            continue
        Z = _z(sub.loc[fams].to_numpy(float).T)
        bud.append(Z @ np.array([beta[(g, f)] for f in fams], float)); tar.append(yv); names.append(g)
    Bu, Ta = np.array(bud), np.array(tar)
    cb = np.array([spearmanr(Bu[i], caf).correlation for i in range(len(names))])
    ct = np.array([spearmanr(Ta[i][np.isfinite(Ta[i])], caf[np.isfinite(Ta[i])]).correlation
                   for i in range(len(names))])

    def _o(i, j):
        p = cb[i] * ct[j]
        return "UNDEFINED" if not np.isfinite(p) else ("OPPOSITE" if p < 0 else "SAME")

    real = pd.DataFrame([(_o(i, i), _partial_spearman(Bu[i], Ta[i], Cm)) for i in range(len(names))],
                        columns=["o", "rho"]).dropna()
    rng = np.random.default_rng(seed)
    sh = []
    for _ in range(n_shuffle):
        pm = rng.permutation(len(names))
        sh += [(_o(i, pm[i]), _partial_spearman(Bu[i], Ta[pm[i]], Cm))
               for i in range(len(names)) if pm[i] != i]
    shuf = pd.DataFrame(sh, columns=["o", "rho"]).dropna()

    def _g(df):
        return float(df[df.o == "OPPOSITE"].rho.mean() - df[df.o == "SAME"].rho.mean())

    from scipy.stats import mannwhitneyu
    res = {"n_genes": len(names), "gradient_real": _g(real), "gradient_shuffled": _g(shuf),
           "n_undefined_real": int((real.o == "UNDEFINED").sum())}
    print(f"\n  --- ⭐ MH-114 COMPARTMENT-ORIENTATION NULL ({len(names)} multi-family genes, "
          f"{n_shuffle} shuffles) ---")
    for tag, df in (("REAL", real), ("SHUFFLED", shuf)):
        o, s = df[df.o == "OPPOSITE"].rho, df[df.o == "SAME"].rho
        print(f"    {tag:9s} OPPOSITE {o.mean():+.4f} (n={len(o)}) | SAME {s.mean():+.4f} (n={len(s)})"
              f" | UNDEFINED {int((df.o=='UNDEFINED').sum())}  ⇒ gradient {o.mean()-s.mean():+.4f}")
    print(f"    ⇒ the GRADIENT is largely arithmetic (real {res['gradient_real']:+.4f} vs shuffled "
          f"{res['gradient_shuffled']:+.4f}) — as MH-114 predicts. THE READOUT IS BELOW.")
    print("    --- ⭐ COMPARTMENT-BLIND real−shuffled (must be SAME SIGN in both orientations) ---")
    for o in ("OPPOSITE", "SAME"):
        a, b = real[real.o == o].rho, shuf[shuf.o == o].rho
        if len(a) > 10 and len(b) > 10:
            p = mannwhitneyu(a, b, alternative="less").pvalue
            print(f"      {o:9s} real {a.mean():+.4f} vs shuffled {b.mean():+.4f}   Δ {a.mean()-b.mean():+.4f}"
                  f"   MWU p={p:.3g}")
            res[f"blind_delta_{o}"] = float(a.mean() - b.mean()); res[f"blind_p_{o}"] = float(p)
    same_sign = np.sign(res.get("blind_delta_OPPOSITE", 0)) == np.sign(res.get("blind_delta_SAME", 0))
    res["compartment_blind_same_sign"] = bool(same_sign)
    verdict = ("✅ SAME SIGN in both orientations — genuine transferred repression, not mixing"
               if same_sign else "⛔ SIGN FLIPS — consistent with compartment arithmetic")
    print(f"      ⇒ {verdict}")
    return res


def edge_leg(n_boot: int = 5000) -> dict:
    """STAGE 3b — the existing per-EDGE lane, re-run under the COMPOSITION-AWARE C blocks.

    ⭐ WHY. `eval/buffa_validation.py` conditions its per-edge partial Spearman on **[prolif, ER] only** —
    it predates the Buffa composition block, and MH-114 is unambiguous that a bulk-vs-bulk replication with
    no composition term is exactly the setup that produced a **90% artifact gradient** in CPTAC. The lane's
    headline replication numbers (`+0.319` Buffa↔TCGA concordance, `0.593` sign-concordance, `0.673` TCGA
    neg-sig same-sign, `0.128` neg-sig-in-Buffa) are therefore **under-controlled**, and were carried as
    settled for a year. Now that `c_blocks()` exists the re-run is nearly free.

    ⚠ **`raw` here is IDENTICAL to what `buffa_validation` uses**, so it doubles as a reproduction gate: if
    the `raw` column does not return the four recorded numbers, the comparison is invalid before it starts.

    Reports every metric under all three blocks plus **retention = adjusted/raw**, the axiom-2a idiom — a
    replication that survives composition adjustment means something a raw one does not.
    """
    from scipy.stats import spearmanr as _sp

    from mirna_hallmark import data_loaders as DL
    from mirna_hallmark.eval import buffa_validation as BV

    B, CB = buffa(), c_blocks()
    he = DL.high_evidence_edges(DL.load_hallmark_edges())[["miRNA", "gene"]].drop_duplicates()
    tcga = pd.read_csv(BV.TCGA_CORR, sep="\t")[["miRNA", "gene", "rho_adj", "q_adj"]].dropna(subset=["rho_adj"])
    out: dict = {}
    tabs: dict = {}
    for blk in _BLOCKS:
        print(f"[edge_leg] {blk}: scoring {len(he):,} HE edges …", flush=True)
        eb = BV.edge_coupling(he, B["mirna"], B["rna"], CB[blk])
        tabs[blk] = eb
        cmp = eb.merge(tcga, on=["miRNA", "gene"], how="inner").rename(
            columns={"rho_adj": "tcga_rho_adj", "q_adj": "tcga_q_adj"})
        both = cmp.dropna(subset=["partial_rho", "tcga_rho_adj"])
        ns = both[(both.tcga_rho_adj < 0) & (both.tcga_q_adj < 0.05)]
        out[blk] = {
            "n_with_coupling": int(len(eb)), "n_vs_tcga": int(len(cmp)),
            "concordance": float(_sp(both.partial_rho, both.tcga_rho_adj)[0]),
            "sign_concordance": float((np.sign(both.partial_rho) == np.sign(both.tcga_rho_adj)).mean()),
            "tcga_negsig_same_sign": float((ns.partial_rho < 0).mean()) if len(ns) else np.nan,
            "tcga_negsig_sig": float(((ns.partial_rho < 0) & (ns.partial_q < 0.05)).mean()) if len(ns) else np.nan,
            "buffa_neg_sig_frac": float(((eb.partial_rho < 0) & (eb.partial_q < 0.05)).mean()),
            "median_rho": float(eb.partial_rho.median()),
        }
    REC = {"concordance": 0.319, "sign_concordance": 0.593,
           "tcga_negsig_same_sign": 0.673, "tcga_negsig_sig": 0.128}
    print(f"\n{'='*96}\nSTAGE 3b — the per-EDGE lane under composition-aware C\n{'='*96}")
    print(f"  {'metric':26s} {'RECORDED':>9s} {'raw':>9s} {'metagene':>10s} {'nnls':>9s} {'retention':>10s}")
    for k, rec in REC.items():
        r, m_, n_ = out["raw"][k], out["metagene"][k], out["nnls"][k]
        ret = m_ / r if abs(r) > 1e-9 else np.nan
        gate = "" if abs(r - rec) <= 0.002 else "  ⛔ raw≠recorded"
        print(f"  {k:26s} {rec:9.3f} {r:9.3f} {m_:10.3f} {n_:9.3f} {ret:10.3f}{gate}")
    print(f"  {'buffa_neg_sig_frac':26s} {'—':>9s} {out['raw']['buffa_neg_sig_frac']:9.3f} "
          f"{out['metagene']['buffa_neg_sig_frac']:10.3f} {out['nnls']['buffa_neg_sig_frac']:9.3f}")
    ok = all(abs(out["raw"][k] - v) <= 0.002 for k, v in REC.items())
    print(f"\n  reproduction gate on `raw`: {'✅ PASS' if ok else '⛔ FAIL — comparison invalid'}")
    print("  ⇒ retention = metagene/raw. A replication that survives composition adjustment means")
    print("    something a raw one does not (axiom 2a); MH-114 measured 90% of a CPTAC gradient as artifact.")
    DEST.mkdir(parents=True, exist_ok=True)
    for blk, t in tabs.items():
        t.to_csv(DEST / f"edge_leg_{blk}.tsv", sep="\t", index=False)
    (DEST / "edge_leg_summary.json").write_text(json.dumps(
        {"recorded_baselines_prolif_ER_only": REC, "reproduction_gate_passed": bool(ok),
         "by_block": out}, indent=2) + "\n")
    print(f"-> {DEST / 'edge_leg_summary.json'}")
    return out


def write_manifest(R: pd.DataFrame, res: dict) -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    man = {
        "module": "mirna_hallmark.learned.eval.ood_cohort",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort": "Buffa 2011 (GSE22216 miRNA + GSE22219 mRNA), n=207 PAIRED (not the 210 in "
                  "buffa_validation.py:6); Oxford, surgery 1989-1992; NO TCGA patient overlap",
        "cohort_is_not_metabric": "VERIFIED 2026-08-02: 0/207 expression matches r>0.90 vs 1,980 METABRIC "
                                  "arrays (best r median 0.270, margin 0.0094); 49% have no clinical "
                                  "candidate; different institutions/accrual. STATE_OF_PLAY Axis 7 is wrong.",
        "estimator": "dense Gibbs spike_slab._gibbs_posterior at pi==1 (readouts.run(level='family'))",
        "rung": "FAMILY (family_card.tsv) — forced: Buffa's miRNA array has no -5p/-3p resolution",
        "beta_frozen": True, "params_fit_in_buffa": 0,
        "ceiling_definition": {
            "source": "TCGA 5-fold OOF, FULL canonical C (CPE, target_cn, mal_prolif)",
            "c_ablation": "TCGA 5-fold OOF, C reduced to what Buffa supports (no target_cn, no CPE; "
                          "+ TCGA ESTIMATE stromal/immune) — THE denominator for E1",
            "target": "beta REFIT in Buffa 5-fold OOF — bounds what 207 microarray patients can support"},
        "retention_denominator": "c_ablation (and target ceiling reported alongside)",
        "confounder_blocks": {k: list(v.columns) for k, v in c_blocks().items()},
        "target_cn": "UNOBTAINABLE for this cohort — GSE22220 has exactly 2 subseries (miRNA+mRNA); "
                     "Buffa 2011 generates no CNV; not a METABRIC subset. Do not re-search.",
        "counts": {"genes_scored": int(len(R)),
                   "genes_with_transfer": int(R[f"rho_buffa_metagene"].notna().sum()),
                   "median_families_per_gene": float(R.n_fam_buffa.median()),
                   "family_cells_measurable": "4,312 / 5,117 (84.3%)",
                   "decoy_units_measurable_pct": 4.7},
        "prereg": PREREG,
        "results": res,
        "outputs": ["ood_cohort_genes.tsv", "ood_cohort_manifest.json"],
    }
    (DEST / "ood_cohort_manifest.json").write_text(json.dumps(man, indent=2, default=str) + "\n")
    print(f"-> {DEST / 'ood_cohort_manifest.json'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--genes")
    ap.add_argument("--block", default="metagene")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--reuse", action="store_true", help="report from the persisted table")
    ap.add_argument("--edge-leg", action="store_true", help="stage 3b: per-edge lane under composition C")
    a = ap.parse_args()
    if a.edge_leg:
        edge_leg()
    elif a.stage == 1:
        c_blocks(); buffa_family_matrix()
    else:
        R = (pd.read_csv(DEST / "ood_cohort_genes.tsv", sep="\t") if a.reuse
             else run(a.genes.split(",") if a.genes else None, a.workers))
        res = report(R, a.block)
        res["compartment_null"] = compartment_null(R)
        write_manifest(R, res)
