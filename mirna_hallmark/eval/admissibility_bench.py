"""CURATED-EDGE ADMISSIBILITY — can an HE edge repress IN BREAST AT ALL? (board DO-FIRST #6)

    .venv/bin/python3 -m mirna_hallmark.eval.admissibility_bench --census        # tag edges (cheap)
    .venv/bin/python3 -m mirna_hallmark.eval.admissibility_bench --bench 3        # firm-up bench
    .venv/bin/python3 -m mirna_hallmark.eval.admissibility_bench --analyze        # width-stratified tests

Two a-priori, Y-blind, model-blind gates per HE edge:
  SITE       — arm has a seed site in the target 3'UTR (TargetScan ctx++ ∪ scanMiR ∪ Poisson-6mer;
               the SITE layers of `decoy_bench._site_maps`, NOT the evidence layer — a foreign-tissue
               CLIP hit is exactly the "curation-transferred" case this is meant to catch).
  EXPRESSION — arm median > _FLOOR=log2(11) (RPM≥10) in TCGA tumour OR GTEx healthy breast.
  admissible = has_site ∧ expressed   (need a site to bind AND dose to be present).

⚠ REUSE, do NOT re-implement: the decoy pool / OOF budget / decoy builder / C block all come from
`decoy_bench` verbatim (its whole docstring is a warning against duplicating them). This module only
FILTERS the real-regulator set and inlines run()'s thin orchestration so the canonical decoy_bench.tsv
is not overwritten.

VERDICT (MH-160, firmed K=3 fold-seeds; see DISCOVERY_REGISTRY):
  * the pre-registered "admissible gap GROWS" is REFUTED (admissible vs all gap_core, pooled p=0.39,
    no width stratum <0.10) — admissibility is NOT a decoy-gap sharpener;
  * guardrail (c) FAILS — inadmissible edges still couple (−0.03 at n_fam 2 & 5+), not ~0;
  * the ONE firmed signal: admissible edges are more cell-intrinsic (composition-invariant) than
    inadmissible, but ONLY in single-family genes (n_fam=1: admis −0.0154 ret 1.07 vs inadmis −0.0027
    ret 0.38, MWU p=0.023; pooled p=0.016 is essentially this one stratum, axiom 5).
  ⇒ keep admissibility as a DIAGNOSTIC flag (axiom 2a); do NOT adopt it as a headline filter.

NB — `decoy_bench --seeds` is ORDER-ONLY (decoys are deterministic across bench seeds; verified
identical). The real MC noise is the FOLD seed (hardcoded 0 in `_one`); this bench averages it.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark.eval import decoy_bench as DB

OUT = DB.OUT                                    # mirna_hallmark/output/learned
_FLOOR = np.log2(11)
TAGS = OUT / "admissibility_edge_tags.tsv"
RES = OUT / "admissibility_bench.tsv"


def _site_only() -> set:
    """SITE layers only (TargetScan ∪ scanMiR ∪ Poisson-6mer) — NOT the evidence-exclusion layer."""
    ts = pd.read_csv("data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt",
                     sep="\t", usecols=["Gene Symbol", "miRNA"], low_memory=False)
    S = set(map(tuple, ts.rename(columns={"Gene Symbol": "gene", "miRNA": "arm"})[["arm", "gene"]].values))
    kd = pd.concat([pd.read_csv(f"data/external_cache/scanmir/{f}", sep="\t", usecols=["arm", "gene"])
                    for f in ("genomewide_kd.tsv.gz", "genomewide_kd_new.tsv.gz", "genomewide_kd_disc.tsv.gz")],
                   ignore_index=True).drop_duplicates()
    S |= set(map(tuple, kd[["arm", "gene"]].values))
    six = OUT / "sixmer_counts.tsv.gz"
    if six.exists():
        from scipy import stats as st
        s6 = pd.read_csv(six, sep="\t")
        S |= set(map(tuple, s6[st.poisson.sf(s6.n6 - 1, s6.expected) < 0.05][["arm", "gene"]].values))
    return S


def build_tags() -> pd.DataFrame:
    """Census: tag every HE edge has_site / expressed / admissible. Writes TAGS."""
    from mirna_hallmark.learned import state as STA
    ctx = DB._ctx()
    HE = ctx["he"][["miRNA", "gene"]].copy()
    HE = HE[HE.miRNA.isin(ctx["X"].index)].drop_duplicates()
    HE.columns = ["arm", "gene"]
    SITE = _site_only()
    tcga_med, gtex_med = ctx["dose"], STA._gtex_mirna().median(axis=1)
    HE["has_site"] = [(a, g) in SITE for a, g in zip(HE.arm, HE.gene)]
    HE["expr_tcga"] = HE.arm.map(lambda a: bool(tcga_med.get(a, -np.inf) > _FLOOR))
    HE["expr_gtex"] = HE.arm.map(lambda a: bool(gtex_med.get(a, -np.inf) > _FLOOR))
    HE["expressed"] = HE.expr_tcga | HE.expr_gtex
    HE["admissible"] = HE.has_site & HE.expressed
    HE.to_csv(TAGS, sep="\t", index=False)
    g = HE.groupby("gene").agg(n_arm=("arm", "size"), n_adm=("admissible", "sum"), n_site=("has_site", "sum"))
    print(f"[census] {len(HE):,} HE edges | has_site {HE.has_site.mean():.1%} | expressed {HE.expressed.mean():.1%} "
          f"| ADMISSIBLE {HE.admissible.mean():.1%} | inadmissible {(~HE.admissible).mean():.1%}")
    print(f"  genes all-seedless {(g.n_site==0).sum()} (MH-136 ~187) | 0-admissible {(g.n_adm==0).sum()} | "
          f"median n_arm {g.n_arm.median():.0f}->{g.n_adm.median():.0f}  ->  {TAGS}")
    return HE


def _rr(keep_admit, admit_set, ctx):
    def rr(gene):
        s = ctx["he"][ctx["he"].gene == gene]
        return [a for a in s.miRNA.unique()
                if a in ctx["X"].index and ((a, gene) in admit_set) == keep_admit]
    return rr


def _one_firm(args):
    """decoy/C/pool ONCE per gene; average real & dec ρ over K FOLD seeds (the real MC axis)."""
    gene, ra, fa, K = args
    try:
        Yv, blocks = DB._C_blocks(gene)
    except Exception:
        return None
    Xr, Xf = DB.pool_family(ra, Yv.index), DB.pool_family(fa, Yv.index)
    o = {"gene": gene, "n_arm": len(ra), "n_fam": Xr.shape[1]}
    for blk, C in blocks.items():
        rc = np.nanmean([DB.oof_budget(Yv, Xr, C, seed=s) for s in range(K)])
        dc = np.nanmean([DB.oof_budget(Yv, Xf, C, seed=s) for s in range(K)])
        o[f"real_{blk}"], o[f"dec_{blk}"], o[f"gap_{blk}"] = rc, dc, rc - dc
    return o


def run_bench(K: int = 3, workers: int = 8) -> pd.DataFrame:
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(v, "1")
    from multiprocessing import Pool
    if not TAGS.exists():
        build_tags()
    tags = pd.read_csv(TAGS, sep="\t")
    ADMIT = {(r.arm, r.gene) for r in tags.itertuples() if r.admissible}
    ctx = DB._ctx()
    orig = DB.real_regulators
    frames = []
    for mode in ("all", "admissible", "inadmissible"):
        DB.real_regulators = {"all": orig, "admissible": _rr(True, ADMIT, ctx),
                              "inadmissible": _rr(False, ADMIT, ctx)}[mode]
        P = DB.build_decoys(genes=None, seed=0)                    # deterministic — once per mode
        jobs = [(g, list(sub.real_arm), list(sub.fake_arm), K) for g, sub in P.groupby("gene")]
        with Pool(workers) as p:
            R = pd.DataFrame([r for r in p.imap_unordered(_one_firm, jobs, chunksize=4) if r])
        R["mode"] = mode
        frames.append(R.dropna(subset=["real_core", "dec_core"]))
        print(f"[{mode}] {len(frames[-1])} genes (K={K})", flush=True)
    DB.real_regulators = orig
    G = pd.concat(frames, ignore_index=True)
    G.to_csv(RES, sep="\t", index=False)
    print(f"  -> {RES}")
    return G


def analyze():
    from scipy import stats
    G = pd.read_csv(RES, sep="\t")
    G["b"] = pd.cut(G.n_fam.round(), [0, 1, 2, 4, 999], labels=["1", "2", "3-4", "5+"])
    print(f"{'nfam':5s} {'mode':13s} {'n':>4s} {'gap_core':>9s} {'gap_deconv':>10s} {'reten':>6s}")
    for nf in ["1", "2", "3-4", "5+"]:
        for m in ("all", "admissible", "inadmissible"):
            s = G[(G.b == nf) & (G["mode"] == m)]
            if len(s) < 10:
                continue
            ret = s.gap_deconv.mean() / s.gap_core.mean() if s.gap_core.mean() else np.nan
            print(f"{nf:5s} {m:13s} {len(s):4d} {s.gap_core.mean():+9.4f} {s.gap_deconv.mean():+10.4f} {ret:6.2f}")
        print()

    def mwu(nf, m1, m2, col):
        a = G[(G["mode"] == m1) & ((G.b == nf) | (nf == "ALL"))][col].dropna()
        b = G[(G["mode"] == m2) & ((G.b == nf) | (nf == "ALL"))][col].dropna()
        return (a.mean(), b.mean(), stats.mannwhitneyu(a, b, alternative="less").pvalue) if len(a) > 9 and len(b) > 9 else None
    for nf in ["1", "2", "3-4", "5+", "ALL"]:
        s, c = mwu(nf, "admissible", "all", "gap_core"), mwu(nf, "admissible", "inadmissible", "gap_deconv")
        if s:
            print(f"  nfam {nf:4s} | sharpen admis {s[0]:+.4f} vs all {s[1]:+.4f} p={s[2]:.3f}"
                  f"  | cell-intrinsic admis {c[0]:+.4f} vs inadmis {c[1]:+.4f} p={c[2]:.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--census", action="store_true")
    ap.add_argument("--bench", type=int, metavar="K", help="firm-up bench, K fold-seeds")
    ap.add_argument("--analyze", action="store_true")
    a = ap.parse_args()
    if a.census:
        build_tags()
    if a.bench:
        run_bench(K=a.bench)
    if a.analyze:
        analyze()
