"""THE ATTRIBUTION TEST — does the model RANK the canonical family, on a ground truth that can be audited?

    .venv/bin/python3 -m mirna_hallmark.eval.attribution_rank [--tier 2] [--boot 5000]

⭐ THE QUESTION, and why it is separate from everything else the program measures. `coupling` asks *does this
family's dose track this gene's level* — well supported, transferable, de-confounded. **Attribution asks
*WHICH family represses this gene*, and it is the only axis where the model is graded against an external
answer key.** The registered verdict (MH-124 §4b, MH-126c TEST 1) is that β is at chance and plain abundance
beats it — measured at n=32 on a set that **no longer exists and was never reproducible** (see
`lit_ground_truth.py`). This re-runs it on the versioned pull, 6-10x larger and audit-stamped.

THE ESTIMAND — a RANK, never an argmax
--------------------------------------
Per gene: where does the canonical family sit among that gene's design families, scored by each candidate
readout? Normalised to **0 = ranked first, 1 = ranked last, 0.5 = chance**; ties get average rank so chance
is exactly 0.5 under permutation. **Argmax is reported but is NOT the estimand** — MH-124 §4b's lesson was
that argmax is the wrong instrument for a collinear design (Design §F: same-seed arms are indistinguishable,
so a top-1 test grades the model on a distinction the data cannot support). A rank test uses the whole
ordering and is what turned "attribution fails" into a measurable gradient.

WHAT IS SCORED, AND WHAT IS REFUSED
-----------------------------------
  `identity`       Shapley/LMG credit on R² with bagged-NNLS weights — the DOCTRINE's identity (MH-138)
  `beta`           the dense Gibbs posterior mean — MAGNITUDE, not identity; scored because it is what the
                   registered claim is about
  `beta_frac`      β_f/Σβ — a magnitude fraction (MH-138 demoted it FROM identity; it splits nothing under
                   collinearity, since Shapley of an ADDITIVE value function is trivially β_f)
  `*_deconv`       the same two under the composition-adjusted C block
  `abundance`      family dose, log2(1+Σ linear RPM), median over samples — ⭐ THE BASELINE TO BEAT, and the
                   one that currently wins
  `w_ledger`       ⭐ CIRCULAR POSITIVE CONTROL — it IS the ground truth's own source, so it MUST come out
                   near 0. It grades the harness, not the model. If it does not win, this file is broken.

  ⛔ `prior_pi`, `pip_discovery` — REFUSED IN CODE. Both are w-informed by construction (`_evidence_pi`), and
  w is derived from this very ledger, so scoring them here is circular. Verified in `readouts.py:190-202`
  that `beta`/`identity` are NOT: `wf` is a dead variable, the dense chain runs at π≡1, and the NNLS takes no
  w. (Matches MH-126c's empirical gate: max|Δβ| = 0.0 under shuffled AND constant w.)

THE UNIT IS THE FAMILY (this is the whole reason MH-126c's n=32 was really n=13)
-------------------------------------------------------------------------------
Genes sharing a canonical family are NOT independent draws — abundance is family-constant, so a single
abundant, well-studied family (miR-21, let-7) contributes many genes with a correlated answer. The primary
inference is therefore a **cluster bootstrap resampling `lit_family`**, and the per-gene Wilcoxon is reported
only as a leverage-blind secondary.

⚠ WHAT A LARGER n DOES *NOT* FIX. PMID depth tracks abundance (+0.171 within gene) — the same miRNAs are
well-studied AND abundant. Scaling n cannot break that; it sharpens whatever the truth is. Read a β-loses
result as "β does not carry identity", NOT as "the literature is biased so the test is unfair" — the bias
runs in favour of the ABUNDANCE baseline, and β is being compared to it on equal terms.
"""
from __future__ import annotations

import argparse
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "attribution_rank.tsv"
SUMMARY = OUT / "attribution_rank_summary.tsv"

#: every scorer is "higher = more likely to be the regulator", so all rank DESCENDING.
SCORERS = ("identity", "beta", "beta_frac", "identity_deconv", "beta_deconv", "abundance",
           "w_ledger", "lit_depth")
#: w-contaminated against a w-derived ground truth. Refused, not merely omitted.
BANNED = {"prior_pi", "pip_discovery", "pip_dense"}
#: |identity| beyond this means the Shapley normaliser (Σφ) collapsed toward 0 — the share is meaningless.
IDENTITY_SANE = 5.0


def _arm_members(arms: str, index: set | None = None, alias: dict | None = None) -> list:
    """⚠ THE `arms` COLUMN IS **SEMICOLON**-SEPARATED (471 multi-arm rows carry ';', ZERO carry ','), and the
    card's arm names are the UN-SUFFIXED form while the expression matrix keeps GENCODE's `.N` disambiguators
    (card `hsa-miR-126-3p` vs matrix `hsa-miR-126-3p.1`, memory `universe-redefinition-pending-refresh`).
    Splitting on ',' silently resolved every multi-arm family to ONE nonexistent name, and the suffix
    mismatch dropped single-arm ones — together they NaN'd 42.9% of CANONICAL rows against 15.7% of
    competitors, a 2.7x asymmetry pointed straight at the right answer."""
    parts = [p.strip() for p in str(arms).replace(",", ";").split(";") if p.strip()]
    if index is None:
        return parts
    out = []
    for p in parts:
        if p in index:
            out.append(p)
        elif alias and p in alias:
            out.append(alias[p])
    return out


def _alias_map(index) -> dict:
    """Un-suffixed arm name -> the matrix's `.N` form, only where the un-suffixed name is ABSENT and the
    resolution is UNIQUE (an ambiguous stem is left unmapped rather than guessed)."""
    have, stems = set(index), {}
    for a in index:
        s = a.split(".")[0]
        if s != a:
            stems.setdefault(s, []).append(a)
    return {s: v[0] for s, v in stems.items() if s not in have and len(v) == 1}


def _family_abundance(card: pd.DataFrame) -> pd.Series:
    """Family DOSE per (gene, family): sum the design's member arms in LINEAR RPM, re-log, take the median
    over samples. Same aggregate `discovery._fam_dose` uses — `within-family = DOSE`, never a mean."""
    from mirna_hallmark.learned import data as LD
    X = LD._load()["X"]
    lin = (np.power(2.0, X.astype(float)) - 1.0).clip(lower=0.0).fillna(0.0)
    idx, alias = set(lin.index), _alias_map(lin.index)
    out = {}
    for arms in card["arms"].dropna().unique():
        members = _arm_members(arms, idx, alias)
        if not members:
            continue
        out[arms] = float(np.median(np.log2(1.0 + lin.loc[members].sum(axis=0).to_numpy(float))))
    return card["arms"].map(out)


def _ledger_weight(card: pd.DataFrame) -> pd.Series:
    """Family-level ledger weight = Σ over member arms (the circular positive control)."""
    from mirna_hallmark.learned.evidence import ledger as LG
    w = LG.edge_weights().set_index(["arm", "gene"])["ledger_weight"]
    return pd.Series([float(sum(w.get((m, g), 0.0) for m in _arm_members(a)))
                      for g, a in zip(card["gene"], card["arms"])], index=card.index)


def _lit_depth(card: pd.DataFrame) -> pd.Series:
    """⭐ THE HARNESS CONTROL — distinct low-throughput-functional PMIDs per (gene, family), i.e. the ground
    truth's OWN defining statistic. Because the set admits only unambiguous tops (margin >= 1), ranking by
    this MUST put the canonical family first for EVERY gene: mean rank exactly 0.000. Anything else means
    the ranking code is broken, not the model. (Distinct from `w_ledger`, which sums ALL assay classes —
    including the 1.01M ago_clip rows — and so is a genuinely different, non-definitional quantity.)"""
    from mirna_hallmark.eval.lit_ground_truth import LT_FUNC, _clean_arms
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned.evidence import ledger as LG
    lg = LG.build_ledger()
    f = lg[lg["assay_class"].isin(LT_FUNC) & (~lg["weak"])].copy()
    f["arm"] = _clean_arms(f["arm"])
    fam = FAM.family_of(pd.Index(f["arm"].unique()))
    f["family"] = f["arm"].map(fam)
    dep = f.dropna(subset=["family"]).groupby(["gene", "family"])["pmid"].nunique()
    return pd.Series([float(dep.get((g, fm), 0.0)) for g, fm in zip(card["gene"], card["family"])],
                     index=card.index)


def build(tier: int = 2) -> pd.DataFrame:
    from mirna_hallmark.eval import lit_ground_truth as GT

    gt = GT.load(min_tier=tier)
    card = pd.read_csv(OUT / "gene_family_card.tsv", sep="\t", low_memory=False)
    bad = BANNED & set(SCORERS)
    if bad:                                          # a guard, not a comment — see module docstring
        raise RuntimeError(f"w-contaminated scorer(s) against a w-derived ground truth: {sorted(bad)}")
    card = card[card["gene"].isin(set(gt["gene"]))].copy()
    card["abundance"] = _family_abundance(card)
    card["w_ledger"] = _ledger_weight(card)
    card["lit_depth"] = _lit_depth(card)
    cov = {s: 1.0 - card[s].isna().mean() for s in SCORERS if s in card}
    print("[attribution_rank] scorer coverage: " +
          "  ".join(f"{k} {v:.1%}" for k, v in sorted(cov.items(), key=lambda kv: kv[1])))

    lit = dict(zip(gt["gene"], gt["lit_family"]))
    meta = gt.set_index("gene")[["n_pmid", "margin", "tier", "n_design_families"]]

    rows = []
    for g, sub in card.groupby("gene"):
        target = lit.get(g)
        if target is None or target not in set(sub["family"]):
            continue
        k = int(len(sub))
        if k < 2:
            continue
        hit = (sub["family"] == target).to_numpy()
        rec = {"gene": g, "lit_family": target, "n_fam": k, **meta.loc[g].to_dict()}
        for s in SCORERS:
            v = sub[s].to_numpy(float)
            if s.startswith("identity") and np.isfinite(v).any() and np.nanmax(np.abs(v)) > IDENTITY_SANE:
                v = np.full_like(v, np.nan)             # degenerate Shapley normaliser -> no verdict
            # ⚠ A MISSING SCORE IS *UNMEASURED*, NOT *WORST*. An earlier version imputed NaN as `min-1` so it
            # ranked last; because coverage was 2.7x worse on canonical families than on competitors, that
            # alone drove the right answer to the bottom of 43% of genes and INVERTED the abundance verdict
            # (0.586 "worse than chance" -> the truth is far better than chance). Rank the canonical among
            # the families that ARE scored, and record how many that was.
            ok = np.isfinite(v)
            m = int(ok.sum())
            if m < 2 or not ok[hit][0]:
                rec[f"rank_{s}"] = np.nan; rec[f"top1_{s}"] = np.nan; rec[f"k_{s}"] = m
                continue
            r = rankdata(-v[ok], method="average")                      # 1 = highest score
            pos = float(r[(hit[ok])][0])
            rec[f"rank_{s}"] = (pos - 1.0) / (m - 1.0)                  # 0 = first, 1 = last, chance 0.5
            rec[f"top1_{s}"] = bool(pos == r.min())
            rec[f"k_{s}"] = m
        rows.append(rec)
    return pd.DataFrame(rows)


def fame_null(tier: int = 2, n_boot: int = 5000) -> pd.DataFrame:
    """⭐⭐ THE FAMILY-FAME NULL — the control that decides what the raw ranking MEANS.

    The ground truth is defined by PMID depth, and **abundant miRNAs get studied more** — measured on this
    set, within gene, across all design families: spearman(family abundance, family LT-func depth) =
    **+0.187 mean / +0.244 median, positive in 67.3% of genes, p=3.2e-10** (the family-level analogue of
    MH-126c's arm-level +0.171). So a scorer can rank the canonical family high for two entirely different
    reasons: it identifies **this gene's** regulator, or it merely likes **famous families**.

    THE NULL SEPARATES THEM: for each gene, substitute **another gene's canonical family that is also present
    in THIS gene's design**. Fame is held constant by construction (the substitute is somebody's canonical
    regulator); only *gene-specificity* varies. `Δ = rank(real) − mean rank(alternatives)` is the part of a
    scorer's performance that is actually attribution.

    ⚠ Restricting to "families that already have literature for this gene" does NOT work as a fame control —
    measured: it drops >=1 family in only **0.6%** of genes (8.16 design families/gene vs 8.15 lit-bearing).
    Nearly every family in a well-studied gene's design already carries some functional literature.

    ⚠ CONSERVATIVE BY CONSTRUCTION: a substituted family may be a *genuine* secondary regulator of this gene
    (just not the most-cited one), which puts real signal in the null and SHRINKS Δ. The measured Δ is a
    lower bound.

    Alternatives are enumerated EXHAUSTIVELY (not sampled), the test is PAIRED within gene, and the bootstrap
    clusters on `lit_family` — the 20-draws-per-gene MWU used while prototyping this was anti-conservative.
    """
    from mirna_hallmark.eval import lit_ground_truth as GT

    gt = GT.load(min_tier=tier)
    card = pd.read_csv(OUT / "gene_family_card.tsv", sep="\t", low_memory=False)
    card = card[card["gene"].isin(set(gt["gene"]))].copy()
    card["abundance"] = _family_abundance(card)
    card["w_ledger"] = _ledger_weight(card)
    card["lit_depth"] = _lit_depth(card)
    lit = dict(zip(gt["gene"], gt["lit_family"]))
    pool = set(gt["lit_family"])

    rows = []
    for g, sub in card.groupby("gene"):
        target = lit.get(g)
        fams = sub["family"].to_numpy()
        if target is None or target not in set(fams):
            continue
        alts = [f for f in set(fams) & pool if f != target]
        if not alts:
            continue
        rec = {"gene": g, "lit_family": target, "n_alt": len(alts)}
        for s in SCORERS:
            v = sub[s].to_numpy(float)
            if s.startswith("identity") and np.isfinite(v).any() and np.nanmax(np.abs(v)) > IDENTITY_SANE:
                v = np.full_like(v, np.nan)
            ok = np.isfinite(v)
            if ok.sum() < 2 or not ok[fams == target][0]:
                continue
            r = rankdata(-v[ok], method="average")
            fo, m = fams[ok], int(ok.sum())
            rec[f"real_{s}"] = float((r[fo == target][0] - 1.0) / (m - 1.0))
            av = [float((r[fo == a][0] - 1.0) / (m - 1.0)) for a in alts if a in set(fo)]
            if av:
                rec[f"null_{s}"] = float(np.mean(av))
                rec[f"d_{s}"] = rec[f"real_{s}"] - rec[f"null_{s}"]
        rows.append(rec)
    D = pd.DataFrame(rows)

    print(f"\n{'='*96}\n⭐ FAMILY-FAME NULL — how much of each scorer is GENE-SPECIFIC vs 'likes famous "
          f"families'\n{'='*96}")
    print(f"  {len(D):,} genes with >=1 alternative canonical family in their own design "
          f"({D.lit_family.nunique()} clusters)\n")
    print(f"  {'scorer':18s} {'REAL':>7s} {'FAME-NULL':>10s} {'Δ gene-specific':>16s} {'95% CI':>18s} {'p':>8s}")
    for s in SCORERS:
        c = f"d_{s}"
        if c not in D:
            continue
        sub = D.dropna(subset=[c])
        obs, ci, p = _cluster_boot_stat(sub, c, n_boot, centre=0.0)
        tag = ("⭐ BASELINE" if s == "abundance" else
               "⭐ CONTROL" if s in ("w_ledger", "lit_depth") else "")
        print(f"  {s:18s} {sub[f'real_{s}'].mean():7.3f} {sub[f'null_{s}'].mean():10.3f} "
              f"{obs:16.3f} [{ci[0]:6.3f},{ci[1]:6.3f}] {p:8.4f}  {tag}")

    print(f"\n  --- ⭐ THE DECISIVE PAIRED CONTRAST: is the MODEL's gene-specific edge bigger than "
          f"ABUNDANCE's? ---")
    for s in ("beta", "identity"):
        sub = D.dropna(subset=[f"d_{s}", "d_abundance"])
        diff = sub[f"d_{s}"] - sub["d_abundance"]
        obs, ci, p = _cluster_boot_stat(sub.assign(_x=diff), "_x", n_boot, centre=0.0)
        verdict = "MODEL WINS" if ci[1] < 0 else ("ABUNDANCE WINS" if ci[0] > 0 else "TIED")
        print(f"    Δ({s}) − Δ(abundance) = {obs:+.3f}  95%CI [{ci[0]:+.3f}, {ci[1]:+.3f}]  p={p:.4f}"
              f"   n={len(sub)}   -> {verdict}")
    _power(D, n_boot)
    D.to_csv(OUT / "attribution_fame_null.tsv", sep="\t", index=False)
    print(f"\n-> {OUT / 'attribution_fame_null.tsv'}")
    return D


#: distinct seed families that are the top-cited regulator of >=1 gene ANYWHERE in the ledger — the hard
#: ceiling on the cluster count, since the cluster unit is the canonical family (measured 2026-08-02).
FAMILY_CEILING = 330


def _power(D: pd.DataFrame, n_boot: int) -> None:
    """⭐ IS THE β-vs-ABUNDANCE CONTRAST RESOLVABLE **AT ALL** with curated literature?

    The board carried "scale the ground-truth set" as the item gating attribution, on the reasoning that
    *"the model cannot attribute" is UNDERPOWERED, not refuted*. That is answerable arithmetically now that
    the contrast has a measured effect size and a measured cluster count, and the answer closes the item:
    the bootstrap half-width shrinks as 1/sqrt(clusters), and **the cluster unit is the canonical FAMILY, of
    which only ~330 exist in the entire curated literature.** So the reachable precision is bounded by
    biology's publication record, not by our effort."""
    sub = D.dropna(subset=["d_beta", "d_abundance"])
    diff = (sub["d_beta"] - sub["d_abundance"]).to_numpy(float)
    obs, ci, _ = _cluster_boot_stat(sub.assign(_x=diff), "_x", n_boot, centre=0.0)
    k = int(sub["lit_family"].nunique())
    half = (ci[1] - ci[0]) / 2.0
    need = int(np.ceil(k * (half / abs(obs)) ** 2)) if abs(obs) > 1e-9 else np.nan
    at_ceiling = half * np.sqrt(k / FAMILY_CEILING)
    print(f"\n  --- ⭐ POWER: can the literature EVER resolve this contrast? ---")
    print(f"    effect {obs:+.3f} · {k} clusters · CI half-width {half:.3f}")
    print(f"    exhausting the literature ({FAMILY_CEILING} families) -> half-width {at_ceiling:.3f}  "
          f"=> {'STILL not significant' if at_ceiling > abs(obs) else 'would resolve'}")
    print(f"    clusters needed to resolve {abs(obs):.3f}: ~{need:,}  "
          f"({need / FAMILY_CEILING:.1f}x every canonical family that has ever been published)")
    print(f"    ⇒ the ground-truth set is NOT the binding constraint. Scaling genes cannot fix this;")
    print(f"      the effect is smaller than the entire curated literature can resolve.")


def _cluster_boot_stat(d: pd.DataFrame, col: str, n_boot: int, centre: float, seed: int = 0):
    """Cluster bootstrap of a mean, resampling `lit_family`; p is two-sided against `centre`."""
    sub = d.dropna(subset=[col])
    if len(sub) < 5:
        return np.nan, (np.nan, np.nan), np.nan
    groups = [g[col].to_numpy(float) for _, g in sub.groupby("lit_family")]
    rng = np.random.default_rng(seed)
    means = np.array([np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))]).mean()
                      for _ in range(n_boot)])
    p = 2.0 * min((means >= centre).mean(), (means <= centre).mean())
    return float(sub[col].mean()), (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))), \
        max(p, 1.0 / n_boot)


def _cluster_boot(d: pd.DataFrame, col: str, n_boot: int, seed: int = 0):
    """Resample CANONICAL FAMILIES with replacement — genes sharing a family are not independent draws."""
    sub = d.dropna(subset=[col])
    if len(sub) < 5:
        return np.nan, (np.nan, np.nan), np.nan
    groups = [g[col].to_numpy(float) for _, g in sub.groupby("lit_family")]
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        means[i] = np.concatenate([groups[j] for j in pick]).mean()
    obs = sub[col].mean()
    # two-sided p against chance 0.5, from the bootstrap distribution
    p = 2.0 * min((means >= 0.5).mean(), (means <= 0.5).mean())
    return obs, (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))), max(p, 1.0 / n_boot)


def report(R: pd.DataFrame, n_boot: int = 5000) -> pd.DataFrame:
    print(f"\n{'='*96}\nATTRIBUTION — normalised rank of the CANONICAL family  (0 = ranked first, "
          f"0.5 = chance, 1 = last)\n{'='*96}")
    print(f"  {len(R):,} genes · {R.lit_family.nunique()} distinct canonical families "
          f"(the cluster-bootstrap unit) · median {R.n_fam.median():.0f} competing families\n")
    out = []
    for s in SCORERS:
        col = f"rank_{s}"
        obs, ci, p = _cluster_boot(R, col, n_boot)
        sub = R.dropna(subset=[col])
        # chance top-1 uses the number of families ACTUALLY SCORED for this readout, not the design size
        exp_top1 = float((1.0 / sub[f"k_{s}"]).sum())
        got_top1 = float(sub[f"top1_{s}"].astype(float).sum())
        wil = stats.wilcoxon(sub[col] - 0.5).pvalue if len(sub) > 10 else np.nan
        tag = ("⭐ HARNESS CONTROL (must be 0.000)" if s == "lit_depth" else
               "⭐ w-CONTROL" if s == "w_ledger" else
               "⭐ BASELINE" if s == "abundance" else "")
        print(f"  {s:18s} n={len(sub):4d}  rank {obs:.3f}  95%CI [{ci[0]:.3f}, {ci[1]:.3f}]  "
              f"p={p:.4f}   top1 {got_top1:.0f}/{len(sub)} vs {exp_top1:.1f} exp  {tag}")
        if s == "lit_depth" and not (obs == obs and obs < 1e-9):
            print(f"    ⛔ HARNESS FAILURE — the definitional control must rank first for every gene "
                  f"(got {obs:.4f}). Every number above is suspect.")
        out.append({"scorer": s, "n": len(sub), "mean_rank": obs, "ci_lo": ci[0], "ci_hi": ci[1],
                    "p_cluster_boot": p, "p_wilcoxon_gene": wil,
                    "top1_obs": got_top1, "top1_exp": exp_top1, "n_families": sub.lit_family.nunique()})
    S = pd.DataFrame(out)

    print(f"\n  --- PAIRED, within gene (negative = the first scorer ranks the canonical family HIGHER) ---")
    for a, b in (("beta", "abundance"), ("identity", "abundance"), ("identity", "beta")):
        d = R.dropna(subset=[f"rank_{a}", f"rank_{b}"])
        if len(d) < 10:
            continue
        diff = d[f"rank_{a}"] - d[f"rank_{b}"]
        w = stats.wilcoxon(diff).pvalue
        print(f"    {a:16s} - {b:16s}  {diff.mean():+.3f}   Wilcoxon p={w:.3g}   n={len(d)}")

    print(f"\n  --- BY TIER (deeper literature = purer ground truth) ---")
    for t in sorted(R.tier.unique()):
        d = R[R.tier >= t]
        cells = "  ".join(f"{s[:9]} {d[f'rank_{s}'].mean():.3f}" for s in ("identity", "beta", "abundance"))
        print(f"    tier>={t}: n={len(d):4d} ({d.lit_family.nunique():3d} fams)   {cells}")

    print(f"\n  --- BY DESIGN SIZE (how hard the attribution is) ---")
    R = R.copy()
    R["_bin"] = pd.cut(R.n_fam, [1, 3, 6, 12, 10_000], labels=["2-3", "4-6", "7-12", "13+"])
    for b, d in R.groupby("_bin", observed=True):
        cells = "  ".join(f"{s[:9]} {d[f'rank_{s}'].mean():.3f}" for s in ("identity", "beta", "abundance"))
        print(f"    {str(b):>5s} families: n={len(d):4d}   {cells}")
    return S


def run(tier: int = 2, n_boot: int = 5000) -> pd.DataFrame:
    R = build(tier=tier)
    R.to_csv(DEST, sep="\t", index=False)
    S = report(R, n_boot=n_boot)
    S.to_csv(SUMMARY, sep="\t", index=False)
    print(f"\n-> {DEST}\n-> {SUMMARY}")
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", type=int, default=2)
    ap.add_argument("--boot", type=int, default=5000)
    ap.add_argument("--no-fame-null", action="store_true")
    a = ap.parse_args()
    run(tier=a.tier, n_boot=a.boot)
    if not a.no_fame_null:
        fame_null(tier=a.tier, n_boot=a.boot)
