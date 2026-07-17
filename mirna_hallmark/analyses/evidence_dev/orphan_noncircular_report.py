"""TASK 2 readout: FDR rates per universe from the persisted per-edge coupling table, + the controls.

⚠ The null is recovered from the pairs ACTUALLY TESTED (set-iteration order is hash-randomised per process,
so re-drawing it here would not reproduce the scan's null). site-free pool := tested edges with NO TargetScan
context++ site, NO scanMiR K_D entry (genome-wide), and not curated.

Two matched comparisons, because arm ABUNDANCE/DETECTION is the dominant confound:
  (A) ABUNDANCE-STRATIFIED (Cochran-Mantel-Haenszel over 10 arm-abundance deciles x detection): the
      abundance-adjusted OR of each universe vs the site-free null. Uses all edges.
  (B) EXACT ARM-MATCHED: for each universe, draw from the site-free pool the SAME number of edges PER ARM
      => arm abundance & detection marginals identical by construction. Deterministic (sorted).
DIRECTIONAL CONTROL: the sig-POSITIVE rate is reported next to the negative one. Real sequence-driven
repression must show its excess specifically in the NEGATIVE tail.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests

from mirna_hallmark.learned import data as LD
from mirna_hallmark.analyses.evidence_dev.orphan_noncircular_universe import build_universes

OUT = Path("mirna_hallmark/output/learned/orphan_noncircular")


def _hits(sub, tag):
    """BH within the given edge set; hit = FDR<0.05 & rho<0 (repressive). Also returns the sig-positive mask."""
    rej, _q, *_ = multipletests(sub[f"p_{tag}"].to_numpy(), alpha=0.05, method="fdr_bh")
    rho = sub[f"rho_{tag}"].to_numpy()
    return (rej & (rho < 0)), (rej & (rho > 0))


def _arm_matched(target: pd.DataFrame, pool: pd.DataFrame, *, seed: int = 0) -> pd.DataFrame:
    """Draw from `pool` the same number of edges PER ARM as `target` has (deterministic)."""
    rng = np.random.default_rng(seed)
    pool_by_arm = {a: sub.index.to_numpy() for a, sub in pool.groupby("arm", sort=True)}
    take = []
    for a, k in sorted(target["arm"].value_counts().items()):
        idx = pool_by_arm.get(a)
        if idx is None or len(idx) == 0:
            continue
        take.append(rng.choice(idx, size=min(k, len(idx)), replace=False))
    return pool.loc[np.concatenate(take)] if take else pool.iloc[:0]


def main(seed: int = 0):
    res = pd.read_parquet(OUT / "edge_coupling_all.parquet")
    res = res[res.p_core.notna() & res.p_dec.notna()].reset_index(drop=True)
    res["edge"] = list(zip(res.gene, res.arm))
    print(f"[table] {len(res)} tested (gene,arm) rows over {res.gene.nunique()} genes | n median {res.n.median():.0f}")

    u = build_universes()
    seq, curated, has_site = u["seq"], u["curated"], u["has_site"]
    old = pd.read_csv("mirna_hallmark/output/learned/discovery_dossier.tsv", sep="\t")
    old_pairs = set(zip(old.gene.astype(str), old.arm.astype(str)))

    d = LD._load(); Xn = d["X"].astype(float)
    abund = Xn.median(axis=1); det = (Xn.fillna(0) > 0).mean(axis=1)
    res["abund"] = abund.reindex(res.arm).to_numpy()
    res["det"] = det.reindex(res.arm).to_numpy()

    is_site = res.edge.isin(has_site)                       # TS site OR any scanMiR K_D entry OR curated
    POOL = res[~is_site].copy()                             # the SITE-FREE null pool, as actually tested
    U = {"curated_HE": res[res.edge.isin(curated)].copy(),
         "old_circular": res[res.edge.isin(old_pairs)].copy(),
         "seq_noncircular": res[res.edge.isin(seq)].copy()}
    print(f"[pool] site-free tested edges = {len(POOL)} over {POOL.arm.nunique()} arms / {POOL.gene.nunique()} genes")

    # ---------- raw rates ----------
    print("\n=== (1) RAW FDR RATES (BH within set, alpha=0.05; hit = FDR<0.05 & rho<0) ===")
    print(f"{'set':17s} {'n':>7s} {'ab_med':>7s} {'det_med':>7s} | {'core neg':>8s} {'rate':>6s} {'pos':>6s} "
          f"| {'dec neg':>8s} {'rate':>6s} {'pos':>6s}")
    rows = []
    for name, df in list(U.items()) + [("sitefree_null_ALL", POOL)]:
        line = f"{name:17s} {len(df):7d} {df.abund.median():7.3f} {df.det.median():7.3f} |"
        rec = {"set": name, "n": len(df), "arm_abund_med": round(float(df.abund.median()), 3),
               "arm_det_med": round(float(df.det.median()), 3)}
        for tag in ["core", "dec"]:
            neg, pos = _hits(df, tag)
            line += f" {int(neg.sum()):8d} {neg.mean():6.3f} {int(pos.sum()):6d} |"
            rec |= {f"sig_neg_{tag}": int(neg.sum()), f"rate_neg_{tag}": round(float(neg.mean()), 4),
                    f"sig_pos_{tag}": int(pos.sum()), f"mean_rho_{tag}": round(float(df[f"rho_{tag}"].mean()), 4)}
        rows.append(rec)
        print(line)

    # ---------- (A) abundance/detection-stratified CMH vs the site-free null ----------
    print("\n=== (2A) ABUNDANCE+DETECTION-STRATIFIED (CMH) vs SITE-FREE NULL — the essential control ===")
    qa = np.unique(np.quantile(POOL.abund.dropna(), np.linspace(0, 1, 11)))
    qd = np.unique(np.quantile(POOL.det.dropna(), np.linspace(0, 1, 4)))
    def _strat(df):
        return (pd.cut(df.abund, qa, include_lowest=True, duplicates="drop").astype(str) + "|" +
                pd.cut(df.det, qd, include_lowest=True, duplicates="drop").astype(str))
    for tag in ["core", "dec"]:
        pneg, ppos = _hits(POOL, tag)
        POOL[f"neg_{tag}"], POOL[f"pos_{tag}"] = pneg, ppos
    POOL["stratum"] = _strat(POOL)
    for name, df in U.items():
        df["stratum"] = _strat(df)
        for tag in ["core", "dec"]:
            neg, pos = _hits(df, tag)
            df[f"neg_{tag}"], df[f"pos_{tag}"] = neg, pos
        for tag, lab in [("core", "core"), ("dec", "deconv")]:
            for dirn in ["neg", "pos"]:
                tabs = []
                for s in sorted(set(df.stratum) & set(POOL.stratum)):
                    a = df[df.stratum == s]; b = POOL[POOL.stratum == s]
                    if len(a) < 20 or len(b) < 20:
                        continue
                    ah = int(a[f"{dirn}_{tag}"].sum()); bh = int(b[f"{dirn}_{tag}"].sum())
                    tabs.append([[ah, len(a) - ah], [bh, len(b) - bh]])
                if not tabs:
                    continue
                st = StratifiedTable(np.array(tabs).transpose(1, 2, 0))
                ptest = st.test_null_odds()
                mark = "  <- directional control" if dirn == "pos" else ""
                print(f"  {name:17s} [{lab:6s}] sig-{dirn.upper()}  CMH OR={st.oddsratio_pooled:.3f} "
                      f"(95% CI {st.oddsratio_pooled_confint()[0]:.3f}-{st.oddsratio_pooled_confint()[1]:.3f})  "
                      f"p={ptest.pvalue:.2e}  [{len(tabs)} strata]{mark}")
        print()

    # ---------- (B) exact arm-matched null ----------
    print("=== (2B) EXACT ARM-MATCHED site-free null (identical arm abundance/detection marginals) ===")
    out = []
    for name, df in U.items():
        nul = _arm_matched(df, POOL, seed=seed)
        chk = (f"abund med {df.abund.median():.3f} vs {nul.abund.median():.3f} | "
               f"det med {df.det.median():.3f} vs {nul.det.median():.3f}")
        print(f"  {name}: n={len(df)} vs matched-null n={len(nul)}  [MATCH CHECK: {chk}]")
        for tag in ["core", "dec"]:
            neg, pos = _hits(df, tag); nneg, npos = _hits(nul, tag)
            o, p = fisher_exact([[int(neg.sum()), len(df) - int(neg.sum())],
                                 [int(nneg.sum()), len(nul) - int(nneg.sum())]], alternative="greater")
            o2, p2 = fisher_exact([[int(pos.sum()), len(df) - int(pos.sum())],
                                   [int(npos.sum()), len(nul) - int(npos.sum())]], alternative="greater")
            print(f"    [{tag:4s}] sig-NEG {neg.mean():.4f} vs null {nneg.mean():.4f}  OR={o:.3f} p={p:.2e}  "
                  f"excess={neg.mean()-nneg.mean():+.4f}   || sig-POS {pos.mean():.4f} vs {npos.mean():.4f} OR={o2:.3f}")
            out.append({"set": name, "C": tag, "rate_neg": round(float(neg.mean()), 4),
                        "null_rate_neg": round(float(nneg.mean()), 4), "OR": round(float(o), 3), "p": p,
                        "rate_pos": round(float(pos.mean()), 4), "null_rate_pos": round(float(npos.mean()), 4)})
    pd.DataFrame(rows).to_csv(OUT / "summary_rates.tsv", sep="\t", index=False)
    pd.DataFrame(out).to_csv(OUT / "summary_vs_matched_null.tsv", sep="\t", index=False)

    # ---------- (3) sequence route + dose-response ----------
    print("\n=== (3) NON-CIRCULAR ORPHANS by SEQUENCE ROUTE, each vs its OWN arm-matched null (core) ===")
    df = U["seq_noncircular"]
    df = df.assign(_ts=df.edge.isin(u["ts_pairs"]), _kd=df.edge.isin(u["kd_pairs"]))
    for name, m in [("TS-only", df._ts & ~df._kd), ("KD-only", ~df._ts & df._kd), ("TS & KD both", df._ts & df._kd)]:
        sub = df[m].copy()
        nul = _arm_matched(sub, POOL, seed=seed)
        neg, pos = _hits(sub, "core"); nneg, npos = _hits(nul, "core")
        o, p = fisher_exact([[int(neg.sum()), len(sub) - int(neg.sum())],
                             [int(nneg.sum()), len(nul) - int(nneg.sum())]], alternative="greater")
        print(f"  {name:13s} n={len(sub):6d} (ab_med {sub.abund.median():.2f} vs null {nul.abund.median():.2f}) "
              f"sig-neg {neg.mean():.4f} vs null {nneg.mean():.4f}  OR={o:.3f} p={p:.2e}  meanRho {sub.rho_core.mean():+.4f}")

    print("\n=== (4) TargetScan context++ STRENGTH dose-response (sequence-only stratifier, core) ===")
    tsm = LD._targetscan_context().set_index(["arm", "gene"])["ts_mag"]
    sub = df[df._ts].copy()
    sub["ts_mag"] = tsm.reindex(pd.MultiIndex.from_arrays([sub.arm, sub.gene])).to_numpy()
    sub = sub[sub.ts_mag.notna()]
    qs = np.quantile(sub.ts_mag, np.linspace(0, 1, 5))
    print(f"{'ts_mag quartile':>20s} {'n':>7s} {'ab_med':>7s} {'sig-neg':>8s} {'mean rho':>9s}")
    for lo, hi in zip(qs[:-1], qs[1:]):
        s2 = sub[(sub.ts_mag >= lo) & (sub.ts_mag < hi)]
        neg, _ = _hits(s2, "core")
        print(f"  [{lo:6.3f},{hi:6.3f})   {len(s2):7d} {s2.abund.median():7.2f} {neg.mean():8.4f} {s2.rho_core.mean():+9.4f}")


if __name__ == "__main__":
    main()
