"""IS THE PROTEIN NULL *BIOLOGY* OR *INSTRUMENT*? — the decoy-controlled layer test at proteome scale.

    .venv/bin/python3 -m mirna_hallmark.learned.eval.layer_boundary --stage rppa      # the matched-gene control
    .venv/bin/python3 -m mirna_hallmark.learned.eval.layer_boundary --stage tcga105   # layer-only, 879 genes
    .venv/bin/python3 -m mirna_hallmark.learned.eval.layer_boundary --stage prospective
    .venv/bin/python3 -m mirna_hallmark.learned.eval.layer_boundary --report

⭐⭐ WHY THIS EXISTS. MH-203/204 established that β beats a **fitted** site-free decoy on mRNA in TCGA
(−0.0189, p=6.7e-09) and in Buffa (−0.0374, p=3.5e-06) but **not on protein** (−0.0086, p=0.554, n=180
RPPA genes). That was recorded as a **layer split**. This module asks whether the protein cell is a real
property of the layer or an artifact of the instrument that measured it — because **three things were
confounded inside that one number**:

  1. **GENE SELECTION.** RPPA's ~200-antibody panel is signalling-biased; the 180 genes it can measure are
     not a random draw from the 1,354.
  2. **POWER.** n=180 genes, not n=866 patients, is what sets the test's resolution.
  3. **LAYER.** The thing we actually wanted to know.

⭐ **THE ANSWER TO (1) AND (2) NEEDED NO NEW DATA, AND IT REFRAMES THE RESULT** (`--stage rppa`):
**on the SAME 180 genes the *mRNA* gap is −0.0103, p=0.0565 — also not significant.** So protein is not
losing to mRNA; **those genes are weak on both layers.** Matched, protein retains **84%** of the mRNA gap
(−0.0086 vs −0.0103), against the **45%** the naive whole-set comparison implies. And the MDE at n=180 is
**|gap| ≥ 0.0320**, so *neither* number is resolvable: **−0.0086 and −0.0103 are both inside the noise
floor.** ⇒ the recorded "protein null" was **underpowered and gene-selected**, not a measured negative.
⚠ The paired within-gene difference (+0.0175, Wilcoxon p=0.027) does *not* survive a sign test
(protein worse in 101/180 = 56.1%, p=0.117), and the protein estimator is the *less* dispersed of the two
(sd ratio 0.75, r=+0.585) — a noisier estimator on identical truth reproduces exactly this pattern.

⭐⭐ **SO THE BINDING CONSTRAINT IS GENES, NOT PATIENTS** — and that is what CPTAC fixes. To detect the
mRNA-sized gap (−0.0189) at 80% power needs **~514 genes**; RPPA supplies 180, CPTAC supplies ~880–935.

WHICH CPTAC — THERE ARE TWO, AND THEY ARE NOT INTERCHANGEABLE
-------------------------------------------------------------
| cohort | patients | miRNA | proteome | decoy-set genes ≥80% cov | boundary |
|---|---|---|---|---|---|
| **`tcga105`** | 105, **THE SAME PATIENTS as the TCGA fit** | **IS the TCGA miRNA** | iTRAQ-4 | **879** | ⭐ **LAYER ONLY** |
| `prospective` | 122, independent | own CPTAC-2 miRNA-seq | TMT-11 | 935 | cohort **+** layer |

⭐ **`tcga105`'s patient-circularity — which `ood_protein.py` correctly flags as a WEAKNESS for validation —
is a FEATURE here**, because this question is about the layer, and sharing patients removes the cohort
boundary entirely. It is the proteome-scale analogue of the RPPA design. ⚠ Its proteome is the **weaker**
one (iTRAQ-4: shallower depth, ratio compression) — so `prospective` (TMT-11) is run too, and there the
**within-cohort mRNA-vs-protein contrast** is the readout, because that differences the cohort boundary out.

⚠⚠ **COMPOSITION IS NOT OPTIONAL AND THIS IS THE EXACT PLACE THE PROJECT WAS BURNED.**
`cptac-protein-composition-confound` / MH-107: the CPTAC protein validation ran with **NO composition
term**, and the correction was a REFRAME (epithelial programs survive, stromal collapse). MH-172
re-derived the cost on gold edges: **27.6% of raw protein-coupled edges SIGN-FLIP under composition and
mean |ρ| falls 39%.** Every cell here is emitted under **both** blocks with retention (axiom 2a).

⛔ **WHAT THIS IS NOT.** Not a revival of `βᵗ` (dead at MH-103, a mediator leak — a design flaw, not a
power problem). Not per-edge significance. Not a coupling lever: protein carries 4–6% of the mRNA
channel's Fisher information about β.
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
from scipy.stats import binomtest, norm, rankdata, wilcoxon

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned/layer_boundary"
LEARNED = C.REPO_ROOT / "mirna_hallmark/output/learned"
BUDGETS = LEARNED / "decoy_oof_budgets.parquet"
SCORED = LEARNED / "decoy_oof_scored.tsv"
FAMILY_CARD = LEARNED / "gene_family_card.tsv"
DECOY_FAMILY_BETA = LEARNED / "decoy_family_betas_oof.tsv"
PAIRS = LEARNED / "decoy_full_pairs.tsv"

MIN_N = 40          # patients per gene; CPTAC cohorts are ~105-122
POWER_Z = norm.ppf(0.975) + norm.ppf(0.80)
ARE = 0.955         # Wilcoxon signed-rank asymptotic relative efficiency vs the t-test

_MEM: dict = {}


# --------------------------------------------------------------------------- #
# Shared statistics
# --------------------------------------------------------------------------- #

def mde(sd: float, n: int) -> float:
    """Minimum detectable |gap| at 80% power, two-sided, signed-rank.

    ⭐ **REPORTED FOR EVERY CELL, BECAUSE THIS TEST'S HISTORY IS A POWER FAILURE MISREAD AS A NEGATIVE.**
    A null whose effect sits inside its own MDE is *unresolved*, not *absent*."""
    return POWER_Z * sd / np.sqrt(max(n, 1)) / ARE


def summarise(gap: pd.Series, label: str) -> dict:
    """Median gap + Wilcoxon + sign test + the MDE it must clear to mean anything."""
    g = pd.Series(gap).dropna()
    if len(g) < 8:
        return {"label": label, "n": len(g), "gap": np.nan, "p": np.nan, "mde": np.nan}
    n_worse = int((g > 0).sum())
    m = mde(float(g.std()), len(g))
    return {"label": label, "n": len(g), "gap": float(g.median()),
            "p": float(wilcoxon(g).pvalue), "p_sign": float(binomtest(n_worse, len(g)).pvalue),
            "frac_real_better": float((g < 0).mean()), "sd": float(g.std()), "mde": float(m),
            "resolvable": bool(abs(g.median()) >= m)}


def _partial(m: np.ndarray, y: np.ndarray, Cm: np.ndarray | None) -> tuple[float, int]:
    """Partial Spearman — Pearson on ranks residualised against ranked C (`realization.partial_spearman`)."""
    ok = np.isfinite(m) & np.isfinite(y)
    if Cm is not None:
        ok &= np.isfinite(Cm).all(1)
    if ok.sum() < MIN_N:
        return np.nan, int(ok.sum())
    R = np.column_stack([rankdata(m[ok]), rankdata(y[ok])])
    cols = [np.ones(ok.sum())]
    if Cm is not None:
        cols += [rankdata(Cm[ok, j]) for j in range(Cm.shape[1])]
    D = np.column_stack(cols)
    res = R - D @ np.linalg.lstsq(D, R, rcond=None)[0]
    if np.std(res[:, 0]) < 1e-9 or np.std(res[:, 1]) < 1e-9:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(res[:, 0], res[:, 1])[0, 1]), int(ok.sum())


# --------------------------------------------------------------------------- #
# Stage RPPA — the matched-gene control (no new data required)
# --------------------------------------------------------------------------- #

def rppa_leg() -> dict:
    """⭐ THE CONTROL THAT REFRAMES MH-203's PROTEIN NULL, computed from artifacts already on disk.

    Compares the decoy gap on **all** mRNA genes, on the **180 protein-measurable genes at the mRNA
    layer**, and on those same genes **at the protein layer**. The middle row is the one that matters:
    it holds the layer fixed and varies only the gene set."""
    T = pd.read_csv(SCORED, sep="\t")
    P = T.dropna(subset=["gap_prot", "gap"])
    rows = [summarise(T.gap, "mRNA · all decoy-paired genes"),
            summarise(P.gap, "mRNA · the 180 protein-measurable genes"),
            summarise(P.gap_prot, "PROTEIN · the same 180 genes")]
    paired = P.gap_prot - P.gap
    rows.append(summarise(paired, "PAIRED protein−mRNA (same genes)"))

    # ── the ANTIBODY audit: is the panel's phospho contamination doing the work? (measured: no)
    from mirna_hallmark.learned.eval import rppa_protein as RP
    R = RP.rppa()
    used = [next((a for a in R["g2ab"].get(g, []) if a in R["M"].columns), None) for g in P.gene]
    phos = pd.Series([u in R["phospho"] for u in used], index=P.index)
    rows.append(summarise(P.loc[~phos, "gap_prot"], "PROTEIN · total-protein probes only"))
    rows.append(summarise((P.gap_prot - P.gap)[~phos], "PAIRED · total-protein probes only"))

    out = {"rows": rows, "n_phospho": int(phos.sum()),
           "retention_matched": float(P.gap_prot.median() / P.gap.median()),
           "retention_naive": float(P.gap_prot.median() / T.gap.median()),
           "sd_ratio_prot_over_mrna": float(P.gap_prot.std() / P.gap.std()),
           "r_between_gaps": float(np.corrcoef(P.gap, P.gap_prot)[0, 1]),
           "n_genes_needed_for_mrna_effect": float((POWER_Z * T.gap.std() / (abs(T.gap.median()) * ARE)) ** 2)}
    return out


# --------------------------------------------------------------------------- #
# Cohort layers
# --------------------------------------------------------------------------- #

def layers(cohort: str) -> dict:
    """(memoized) CPTAC protein / mRNA / protein-residualised-on-mRNA for one cohort.

    ⚠ `rna_z` here is CPTAC's OWN mRNA on the SAME patients and the SAME cohort-z pipeline as the protein —
    that is what makes the mRNA-vs-protein contrast apples-to-apples. Scoring against TCGA's `Y` instead
    would reintroduce a platform difference into the very contrast the design is trying to isolate."""
    key = f"layers_{cohort}"
    if key not in _MEM:
        from mirna_hallmark.eval import cptac_validation as CV
        _MEM[key] = CV.load_cptac_layers(cohort)
    return _MEM[key]


def comp_C(cohort: str, samples) -> np.ndarray | None:
    """The composition block. ⚠ HARD-FAILS rather than silently falling back to intercept-only — an
    unadjusted protein number is exactly the MH-107 defect this module exists to avoid."""
    from mirna_hallmark.learned import confounders as CF
    key = f"C_{cohort}"
    if key not in _MEM:
        _MEM[key] = CF.build_C("cptac" if cohort == "prospective" else "tcga", list(samples))
    Cb = _MEM[key]
    if Cb is None:
        raise RuntimeError(f"composition block unavailable for {cohort} — refusing to score unadjusted")
    Cm = Cb.reindex(list(samples)).apply(pd.to_numeric, errors="coerce")
    return Cm.fillna(Cm.median()).to_numpy(float)


# --------------------------------------------------------------------------- #
# Stage tcga105 — LAYER ONLY, at proteome scale
# --------------------------------------------------------------------------- #

def tcga105_leg(min_cov: float = 0.8) -> pd.DataFrame:
    """⭐⭐ THE DESIGN-APPROPRIATE INSTRUMENT: same patients as the TCGA fit ⇒ **no cohort boundary**.

    Because `tcga105` reuses the TCGA miRNA, the existing **out-of-fold** budgets apply directly — the
    identical estimator, folds and C that produced the mRNA headline, so the layer is the only thing that
    changes. No re-transport, and therefore no transport artifact to argue about."""
    B = pd.read_parquet(BUDGETS)
    L = layers("tcga105")
    prot, rna, resid = L["protein_z"], L["rna_z"], L["protein_resid"]
    parts = sorted(set(B.participant) & set(prot.columns))
    genes = sorted(set(B.gene) & set(prot.index))
    cov = prot.loc[genes, parts].notna().mean(axis=1)
    genes = [g for g in genes if cov[g] >= min_cov]
    print(f"[layer_boundary/tcga105] {len(genes):,} genes × {len(parts)} patients "
          f"(≥{min_cov:.0%} protein coverage)", flush=True)
    Cm = comp_C("tcga105", parts)
    Bi = B[B.participant.isin(parts)].set_index(["gene", "participant"]).sort_index()
    rows = []
    for g in genes:
        try:
            sub = Bi.loc[g].reindex(parts)
        except KeyError:
            continue
        br, bd = sub.budget_real.to_numpy(float), sub.budget_decoy.to_numpy(float)
        rec = {"gene": g, "n_fam_real": float(sub.n_fam_real.iloc[0])}
        for lab, src in (("prot", prot), ("mrna", rna), ("disc", resid)):
            yv = pd.to_numeric(src.loc[g, parts], errors="coerce").to_numpy(float)
            for blk, Cx in (("core", None), ("comp", Cm)):
                rr, n = _partial(br, yv, Cx)
                rd, _ = _partial(bd, yv, Cx)
                rec |= {f"real_{lab}_{blk}": rr, f"decoy_{lab}_{blk}": rd,
                        f"gap_{lab}_{blk}": rr - rd, f"n_{lab}": n}
        rows.append(rec)
    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "tcga105_gaps.tsv", sep="\t", index=False)
    print(f"-> {OUT/'tcga105_gaps.tsv'}  ({len(T):,} genes)")
    return T


# --------------------------------------------------------------------------- #
# Stage prospective — independent patients, TMT-11; β must be TRANSPORTED
# --------------------------------------------------------------------------- #

def _pool(members, mat: pd.DataFrame) -> np.ndarray | None:
    """TRUE RPM family pool `log2(1 + Σ(2^x − 1))` — never a `groupby.mean()` (a documented decoy-bench defect)."""
    mem = [m for m in members if m in mat.index]
    if not mem:
        return None
    lin = np.nan_to_num(np.power(2.0, mat.loc[sorted(set(mem))].to_numpy(float)) - 1.0, nan=0.0)
    return np.log2(1.0 + np.clip(lin, 0, None).sum(axis=0))


def _z(A: np.ndarray) -> np.ndarray:
    """z WITHIN the cohort — bridge-free: no QN onto TCGA's distribution shape, only the (weaker)
    assumption that the families' rank-order of informativeness transports."""
    mu, sd = np.nanmean(A, 0), np.nanstd(A, 0)
    Z = (A - mu) / (sd + 1e-9)
    Z[:, sd < 0.1] = 0.0
    return np.nan_to_num(Z)


def prospective_leg(min_cov: float = 0.8) -> pd.DataFrame:
    """Independent patients + TMT-11 protein. Confounds cohort with layer — so the readout here is the
    **within-cohort mRNA-vs-protein contrast**, in which the cohort boundary is common to both terms and
    differences out."""
    from mirna_hallmark.coupling_inference import seed_family_map
    from mirna_hallmark.learned import cptac_data as CD

    L = layers("prospective")
    prot, rna = L["protein_z"], L["rna_z"]
    resid = L["protein_resid"]
    X = CD.arms().apply(pd.to_numeric, errors="coerce")
    parts = sorted(set(prot.columns) & set(rna.columns) & set(X.columns))
    if len(parts) < MIN_N:
        raise RuntimeError(f"only {len(parts)} shared prospective participants")
    X = X[parts]

    fam = pd.read_csv(FAMILY_CARD, sep="\t", usecols=["gene", "family", "arms", "beta"])
    beta_r = {(g, str(f)): b for g, f, b in zip(fam.gene, fam.family, fam.beta)}
    d = pd.read_csv(DECOY_FAMILY_BETA, sep="\t")
    d = d[d["decoy"].astype(bool)]
    beta_d = {(g, str(u)): b for g, u, b in zip(d.gene, d.unit, d.beta)}
    pairs = pd.read_csv(PAIRS, sep="\t")
    fake_by_gene = pairs.groupby("gene").fake_arm.apply(list).to_dict()
    fmap = seed_family_map(pd.Index(X.index))

    genes = sorted(set(fam.gene) & set(prot.index) & set(fake_by_gene))
    cov = prot.loc[genes, parts].notna().mean(axis=1)
    genes = [g for g in genes if cov[g] >= min_cov]
    print(f"[layer_boundary/prospective] {len(genes):,} genes × {len(parts)} patients", flush=True)
    Cm = comp_C("prospective", parts)
    fam_arms = fam.groupby("gene").apply(lambda s: list(zip(s.family, s.arms)), include_groups=False).to_dict()

    rows = []
    for g in genes:
        # real design: this gene's own families, pooled from CPTAC arms
        cols_r, br_ = [], []
        for f, arms in fam_arms.get(g, []):
            v = _pool([a for a in str(arms).split(";")], X)
            if v is not None and (g, str(f)) in beta_r:
                cols_r.append(v); br_.append(beta_r[(g, str(f))])
        # decoy design: this gene's matched fake arms, grouped by seed family
        by_f: dict = {}
        for a in fake_by_gene.get(g, []):
            by_f.setdefault(str(fmap.get(a)), []).append(a)
        cols_d, bd_ = [], []
        for f, mem in by_f.items():
            if f == "None" or (g, f) not in beta_d:
                continue
            v = _pool(mem, X)
            if v is not None:
                cols_d.append(v); bd_.append(beta_d[(g, f)])
        if not cols_r or not cols_d:
            continue
        mr = _z(np.column_stack(cols_r)) @ np.array(br_, float)
        md = _z(np.column_stack(cols_d)) @ np.array(bd_, float)
        rec = {"gene": g, "n_fam_real": len(cols_r), "n_fam_decoy": len(cols_d)}
        for lab, src in (("prot", prot), ("mrna", rna), ("disc", resid)):
            yv = pd.to_numeric(src.loc[g, parts], errors="coerce").to_numpy(float)
            for blk, Cx in (("core", None), ("comp", Cm)):
                rr, n = _partial(mr, yv, Cx)
                rd, _ = _partial(md, yv, Cx)
                rec |= {f"real_{lab}_{blk}": rr, f"decoy_{lab}_{blk}": rd,
                        f"gap_{lab}_{blk}": rr - rd, f"n_{lab}": n}
        rows.append(rec)
    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "prospective_gaps.tsv", sep="\t", index=False)
    print(f"-> {OUT/'prospective_gaps.tsv'}  ({len(T):,} genes)")
    return T


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #

def _fmt(r: dict) -> str:
    if not np.isfinite(r.get("gap", np.nan)):
        return f"  {r['label']:46s}  n={r['n']:4d}   (insufficient)"
    # ⚠ "<80% power" ≠ "null": the signed-rank significance threshold is 1.96/2.80 = 0.70×MDE, so a result
    # can be p<0.05 while the study had under 80% power for that effect. Read the flag as a winner's-curse
    # warning on the MAGNITUDE, never as a verdict on the SIGN.
    flag = "80% powered" if r["resolvable"] else "<80% power"
    return (f"  {r['label']:46s}  n={r['n']:4d}  gap={r['gap']:+.4f}  p={r['p']:<9.3g}"
            f"  MDE={r['mde']:.4f}  [{flag}]")


def report() -> dict:
    print("\n" + "=" * 100)
    print("IS THE PROTEIN NULL BIOLOGY OR INSTRUMENT?".center(100))
    print("=" * 100)

    R = rppa_leg()
    print("\n── STAGE 1 · RPPA matched-gene control (no new data) " + "─" * 46)
    for r in R["rows"]:
        print(_fmt(r))
    print(f"\n  retention on MATCHED genes : {R['retention_matched']:.0%}"
          f"      (naive, vs the full mRNA set: {R['retention_naive']:.0%})")
    print(f"  sd(protein gap)/sd(mRNA gap): {R['sd_ratio_prot_over_mrna']:.2f}"
          f"   r between the two gaps: {R['r_between_gaps']:+.3f}")
    print(f"  ⇒ genes needed to detect the mRNA-sized gap at 80% power: "
          f"{R['n_genes_needed_for_mrna_effect']:,.0f}   (RPPA supplies 180)")

    out = {"rppa": R, "cptac": {}}
    for coh, f in (("tcga105", OUT / "tcga105_gaps.tsv"), ("prospective", OUT / "prospective_gaps.tsv")):
        if not f.exists():
            continue
        T = pd.read_csv(f, sep="\t")
        tag = "LAYER ONLY — same patients as the fit" if coh == "tcga105" else "cohort + layer"
        print(f"\n── STAGE 2 · CPTAC {coh}  ({tag}) " + "─" * max(4, 60 - len(coh) - len(tag)))
        rows = []
        for blk, bl in (("core", "no composition"), ("comp", "COMPOSITION-ADJUSTED")):
            for lab, nice in (("mrna", "mRNA"), ("prot", "protein"), ("disc", "protein-beyond-mRNA")):
                c = f"gap_{lab}_{blk}"
                if c in T:
                    rows.append(summarise(T[c], f"{nice:20s} [{bl}]"))
            paired = T[f"gap_prot_{blk}"] - T[f"gap_mrna_{blk}"]
            rows.append(summarise(paired, f"{'PAIRED protein−mRNA':20s} [{bl}]"))
        for r in rows:
            print(_fmt(r))
        # ⚠ DENOMINATOR GATE (my own recorded trap). A retention is a ratio; when the denominator is a
        # near-zero gap the ratio is noise amplification, not a measurement. `prospective`'s raw protein
        # gap is −0.0003 and the ungated ratio printed **−3.53**, which means nothing at all.
        for lab, num, den in (("mRNA", "gap_mrna_comp", "gap_mrna_core"),
                              ("protein", "gap_prot_comp", "gap_prot_core")):
            d = float(T[den].median())
            msg = f"{T[num].median()/d:+.2f}" if abs(d) >= 0.005 else f"UNDEFINED (denominator {d:+.4f} ≈ 0)"
            print(f"  composition retention · {lab:8s}: {msg}")
        out["cptac"][coh] = {"n_genes": int(len(T)), "rows": rows}
    # ── THE DECISIVE STATISTIC: the paired protein−mRNA contrast, across every instrument.
    # Paired on GENE, so gene selection cancels; within a cohort, so the cohort boundary cancels too.
    # Whatever is left is the LAYER.
    T0 = pd.read_csv(SCORED, sep="\t").dropna(subset=["gap_prot", "gap"])
    inst = [("RPPA (866 pts, 180 genes)", (T0.gap_prot - T0.gap).dropna())]
    for coh, f in (("CPTAC tcga105 (105 pts)", OUT / "tcga105_gaps.tsv"),
                   ("CPTAC prospective (101 pts)", OUT / "prospective_gaps.tsv")):
        if f.exists():
            D = pd.read_csv(f, sep="\t")
            inst.append((coh, (D.gap_prot_comp - D.gap_mrna_comp).dropna()))
    print("\n── STAGE 3 · THE LAYER CONTRAST, PAIRED ON GENE " + "─" * 51)
    print("     (>0 ⇒ protein WORSE than mRNA on the SAME genes; gene selection & cohort both cancel)")
    zs = []
    for lab, d in inst:
        w = wilcoxon(d).pvalue
        zs.append(norm.isf(w / 2) * np.sign(d.median()))
        print(f"  {lab:32s} n={len(d):4d}  paired={d.median():+.4f}  p={w:<9.3g}  MDE={mde(float(d.std()), len(d)):.4f}")
    z = float(np.sum(zs) / np.sqrt(len(zs)))
    print(f"\n  Stouffer across instruments: z={z:+.2f}  p={2*norm.sf(abs(z)):.3g}")
    print("  ⚠ ANTI-CONSERVATIVE — the gene sets OVERLAP (151 of RPPA's 180 are also in tcga105),")
    print("    so these are not independent tests. Read it as a direction check, not a p-value.")
    out["layer_contrast"] = {"instruments": [{"label": l, "n": int(len(d)), "paired": float(d.median()),
                                              "p": float(wilcoxon(d).pvalue)} for l, d in inst],
                             "stouffer_z": z, "stouffer_p": float(2 * norm.sf(abs(z))),
                             "caveat": "gene sets overlap; not independent"}

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "layer_boundary_report.json").write_text(json.dumps(
        {"generated": datetime.now(timezone.utc).isoformat(), **out}, indent=2, default=float))
    print(f"\n-> {OUT/'layer_boundary_report.json'}")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["rppa", "tcga105", "prospective", "all"], default="all")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--min-cov", type=float, default=0.8)
    a = ap.parse_args()
    if a.stage in ("tcga105", "all"):
        tcga105_leg(a.min_cov)
    if a.stage in ("prospective", "all"):
        prospective_leg(a.min_cov)
    report()
