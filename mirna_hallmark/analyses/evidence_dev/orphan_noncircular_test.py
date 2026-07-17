"""TASK 2 step 3-5: test the NON-CIRCULAR orphan universe vs curated HE vs the OLD circular orphan set,
against a SITE-FREE, arm-degree-matched (=abundance/detection-matched) NULL.

Per edge: partial Spearman(arm, gene | C)  ==  Spearman(resid(arm|C), resid(target|C))   [codebase convention]
  C core   = CPE + target_cn + mal_prolif                     (assemble_gene default)
  C deconv = core + 8 Wu-major non-malignant lineage fractions
p from the t-approximation with df = n - 2 - rank(C); BH-FDR **within each universe**; a "hit" = FDR<0.05
AND rho<0 (a repressive edge).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.stats import t as tdist
from statsmodels.stats.multitest import multipletests

from mirna_hallmark.learned import data as LD

from mirna_hallmark.analyses.evidence_dev.orphan_noncircular_universe import build_universes

OUT = Path("mirna_hallmark/output/learned/orphan_noncircular")


# ───────────────────────────── null universes ─────────────────────────────
def null_degree_preserving(seq: set, has_site: set, he_genes: list, *, seed: int = 0):
    """Site-free null preserving BOTH marginals: each arm keeps its edge COUNT (=> arm abundance/detection
    matched EXACTLY), and gene partners are drawn ∝ the gene's remaining degree quota (=> gene marginals
    approximately preserved). Returns (null_pairs, matched_seq_pairs) where matched_seq has the SAME per-arm
    degrees as the null (some arms' site-free pools are smaller than their degree, so the null is capped;
    matched_seq applies the same cap so the two sets are exactly arm-comparable)."""
    rng = np.random.default_rng(seed)
    by_arm: dict = {}
    gene_deg: dict = {}
    for g, a in seq:
        by_arm.setdefault(a, []).append(g)
        gene_deg[g] = gene_deg.get(g, 0) + 1
    genes = np.array(he_genes)
    quota = np.array([gene_deg.get(g, 0) for g in genes], float)
    null, matched = set(), set()
    for a in sorted(by_arm, key=lambda x: -len(by_arm[x])):        # high-degree arms first (tightest pools)
        gs = by_arm[a]
        free = np.array([i for i, g in enumerate(genes) if (g, a) not in has_site])
        if free.size == 0:
            continue
        k = min(len(gs), free.size)
        w = quota[free] + 1e-6                                     # gene-degree-preserving propensity
        p = w / w.sum()
        picks = rng.choice(free, size=k, replace=False, p=p)
        for i in picks:
            null.add((genes[i], a))
            quota[i] = max(quota[i] - 1, 0)
        for g in rng.permutation(np.array(gs, dtype=object))[:k]:  # same cap on the real set
            matched.add((g, a))
    return null, matched


# ───────────────────────────── the coupling engine ─────────────────────────────
def _resid(M: np.ndarray, C: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(C, M, rcond=None)
    return M - C @ beta


def _spearman_block(xr: np.ndarray, yr: np.ndarray, dof: int):
    """Vectorised Spearman of each column of xr against yr (both already residualised), + two-sided p."""
    rx = rankdata(xr, axis=0)
    ry = rankdata(yr)
    rx = rx - rx.mean(0)
    ry = ry - ry.mean()
    sx = np.sqrt((rx ** 2).sum(0))
    sy = np.sqrt((ry ** 2).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = (rx * ry[:, None]).sum(0) / (sx * sy)
    rho = np.clip(np.nan_to_num(rho, nan=np.nan), -0.999999, 0.999999)
    tstat = rho * np.sqrt(dof / np.maximum(1 - rho ** 2, 1e-12))
    p = 2 * tdist.sf(np.abs(tstat), dof)
    return rho, p


def run(seed: int = 0, progress: int = 200) -> pd.DataFrame:
    u = build_universes()
    X, Y = u["X"], u["Y"]
    seq, curated, has_site = u["seq"], u["curated"], u["has_site"]
    he_genes = u["he_genes"]

    old = pd.read_csv("mirna_hallmark/output/learned/discovery_dossier.tsv", sep="\t")
    old_pairs = set(zip(old["gene"].astype(str), old["arm"].astype(str)))

    null, seq_matched = null_degree_preserving(seq, has_site, he_genes, seed=seed)
    print(f"[universes] curated={len(curated)} old_circular={len(old_pairs)} seq={len(seq)} "
          f"seq_matched={len(seq_matched)} null={len(null)}", flush=True)

    UNIV = {"curated_HE": curated, "old_circular": old_pairs, "seq_noncircular": seq,
            "seq_matched": seq_matched, "sitefree_null": null}

    need: dict = {}
    for name, s in UNIV.items():
        for g, a in s:
            need.setdefault(g, set()).add(a)

    # ---- confounder blocks (participant-level parts, gene-specific target_cn appended per gene)
    d = LD._load()
    clin = d["C"]
    conf = [c for c in LD.confounder_columns() if c in clin.columns]      # ('CPE',)
    parts = clin.index.intersection(Y.columns).intersection(X.columns)
    base = clin.loc[parts, conf].apply(pd.to_numeric, errors="coerce")
    mp = LD._malignant_prolif(base.index)
    base = base.join(mp.rename("mal_prolif"))
    dv = LD._deconv()
    dvj = dv.reindex(base.index)
    dvj = dvj.fillna(dvj.median())
    base_dec = base.join(dvj)
    tcn_all = LD._target_cn(sorted(need))
    xarms = set(X.index)

    rows = []
    genes = sorted(need)
    for i, g in enumerate(genes):
        if progress and i % progress == 0:
            print(f"[couple] {i}/{len(genes)} genes, {len(rows)} edges", flush=True)
        arms = sorted(a for a in need[g] if a in xarms)
        if not arms or g not in Y.index:
            continue
        yv = pd.to_numeric(Y.loc[g, parts], errors="coerce")
        tcn = tcn_all.get(g)
        Ccore = base.copy()
        Cdec = base_dec.copy()
        if tcn is not None:
            tv = pd.to_numeric(tcn.reindex(base.index), errors="coerce")
            Ccore["target_cn"] = tv
            Cdec["target_cn"] = tv
        Ccore = Ccore.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))
        Cdec = Cdec.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))
        keep = yv.notna() & base[conf].notna().all(axis=1)
        if keep.sum() < 60:
            continue
        pk = base.index[keep.to_numpy()]
        xv = X.loc[arms, pk].T.astype(float).fillna(0.0).to_numpy()       # participant × arm
        yy = yv.loc[pk].to_numpy(float)
        out = {}
        for tag, Cdf in [("core", Ccore), ("dec", Cdec)]:
            Cm = np.column_stack([np.ones(len(pk)), Cdf.loc[pk].to_numpy(float)])
            dof = len(pk) - 2 - np.linalg.matrix_rank(Cm)
            rho, p = _spearman_block(_resid(xv, Cm), _resid(yy, Cm), dof)
            out[tag] = (rho, p)
        for j, a in enumerate(arms):
            rows.append({"gene": g, "arm": a, "n": len(pk),
                         "rho_core": out["core"][0][j], "p_core": out["core"][1][j],
                         "rho_dec": out["dec"][0][j], "p_dec": out["dec"][1][j]})
    res = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    res.to_parquet(OUT / "edge_coupling_all.parquet")

    # ---- per-universe BH-FDR
    res["edge"] = list(zip(res.gene, res.arm))
    idx = res.set_index("edge")
    summ = []
    for name, s in UNIV.items():
        sub = idx.loc[idx.index.intersection(list(s))].copy()
        sub = sub[sub["p_core"].notna() & sub["p_dec"].notna()]
        row = {"universe": name, "n_edges_defined": len(s), "n_tested": len(sub)}
        for tag in ["core", "dec"]:
            rej, q, *_ = multipletests(sub[f"p_{tag}"].to_numpy(), alpha=0.05, method="fdr_bh")
            neg = rej & (sub[f"rho_{tag}"].to_numpy() < 0)
            pos = rej & (sub[f"rho_{tag}"].to_numpy() > 0)
            row[f"sig_neg_{tag}"] = int(neg.sum())
            row[f"rate_neg_{tag}"] = round(float(neg.mean()), 4)
            row[f"sig_pos_{tag}"] = int(pos.sum())
            row[f"mean_rho_{tag}"] = round(float(sub[f"rho_{tag}"].mean()), 4)
            sub[f"q_{tag}"] = q
        summ.append(row)
        sub.reset_index().to_parquet(OUT / f"edges_{name}.parquet")
    S = pd.DataFrame(summ)
    S.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    print("\n=== FDR-SIGNIFICANT NEGATIVE (repressive) RATE, BH within universe ===")
    print(S.to_string(index=False))
    return S


if __name__ == "__main__":
    run(seed=int(sys.argv[1]) if len(sys.argv) > 1 else 0)
