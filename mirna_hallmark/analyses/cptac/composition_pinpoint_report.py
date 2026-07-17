"""TASK 4 report — read the pinpoint tables and print every number, with the CONTROLS.

.venv/bin/python3 -m mirna_hallmark.analyses.cptac.composition_pinpoint_report
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, wilcoxon

OUT = Path("mirna_hallmark/output/learned")
LIN = ["CAFs", "T-cells", "Myeloid", "B-cells", "Endothelial", "PVL", "Normal Epithelial", "Plasmablasts"]
CONF = ["composition_explained", "partial"]


def load(level):
    return pd.read_csv(OUT / f"composition_pinpoint_edges_{level}.tsv", sep="\t")


def sec(t):
    print("\n" + "=" * 100 + f"\n{t}\n" + "=" * 100)


def main():
    for level in ("family", "arm"):
        d = load(level)
        sec(f"LEVEL = {level.upper()}   ({len(d)} couplings, {d.gene.nunique()} genes, n={int(d.n.iloc[0])} patients)")

        print("composition_class:")
        print("  " + d.composition_class.value_counts().to_string().replace("\n", "\n  "))

        # ---- (0) aggregate add-one-lineage attribution (reproduces MH-111's 10.7/6.3/5.6/4.1) ----
        c = d[d.composition_class.isin(CONF)].copy()
        rel = c[c.retention_reliable]
        print(f"\n--- ADD-ONE-LINEAGE median loss, ALL {len(d)} couplings (MH-111 reproduction) ---")
        agg = pd.Series({L: d[f"loss_{L}"].median() for L in LIN}).sort_values(ascending=False)
        aggn = pd.Series({L: d[f"null_{L}"].median() for L in LIN})
        print(pd.DataFrame({"median_loss": agg, "median_SHUFFLED_null": aggn.reindex(agg.index),
                            "joint(all 8)": d.full_loss.median()}).round(4).to_string())

        # ---- (1) PINPOINT distribution on the confounded edges ----
        print(f"\n--- (1) PINPOINT on {len(c)} confounded edges ({len(rel)} retention-reliable) ---")
        pin = c[c.single_lineage]
        mix = c[~c.single_lineage]
        print(f"  single-lineage pinpointed : {len(pin):5d}  ({100*len(pin)/len(c):.1f}%)")
        print(f"  MIX (no single lineage)   : {len(mix):5d}  ({100*len(mix)/len(c):.1f}%)")
        print("\n  pinpointed lineage distribution (single-lineage edges):")
        vc = pin.pinpointed_lineage.value_counts()
        for L in LIN:
            print(f"    {L:20s} {vc.get(L,0):5d}  ({100*vc.get(L,0)/max(len(pin),1):5.1f}%)")

        # CONTROL: same argmax on the SHUFFLED-lineage losses -> the null's own pinpoint distribution
        NU = c[[f"null_{L}" for L in LIN]].to_numpy(float)
        nb = np.nanargmax(np.where(np.isfinite(NU), NU, -np.inf), 1)
        nullpin = pd.Series([LIN[i] for i in nb]).value_counts()
        print("\n  CONTROL — argmax over the 8 SHUFFLED lineages (df placebo, same argmax):")
        for L in LIN:
            print(f"    {L:20s} {nullpin.get(L,0):5d}  ({100*nullpin.get(L,0)/len(c):5.1f}%)")
        rate = float((c.lineage_loss_fraction > c.max_null_loss).mean())
        print(f"\n  CONTROL — real max-loss > SHUFFLED max-loss in {100*rate:.1f}% of confounded edges")
        print(f"     real  max-over-8 loss: median {c.lineage_loss_fraction.median():.4f}")
        print(f"     NULL  max-over-8 loss: median {c.max_null_loss.median():.4f}")
        try:
            w = wilcoxon(c.lineage_loss_fraction.dropna(), c.max_null_loss.loc[c.lineage_loss_fraction.dropna().index])
            print(f"     Wilcoxon (real > null): p={w.pvalue:.3g}")
        except Exception as e:                                     # noqa: BLE001
            print("     Wilcoxon failed:", e)

        # ---- (2) ORIENTATION cross ----
        # ⚠ TAUTOLOGY GUARD (measured): at the pinpointed lineage, orientation is 98% determined by
        # sign(beta_core) — P(same|b<0)=0.981, P(same|b>0)=0.159. The lineage that most REDUCES beta is by
        # construction the one whose confound path carries beta's own sign. So the orientation test is only
        # meaningful on REPRESSIVE couplings (beta_core > 0).
        print("\n--- (2a) TAUTOLOGY CHECK: orientation vs sign(beta_core) ---")
        d["b_sign"] = np.where(d.beta_core > 0, "repressive(b>0)", "anti-repressive(b<0)")
        print("    " + pd.crosstab(d.b_sign, d.orientation).to_string().replace("\n", "\n    "))
        d = d[d.beta_core > 0].copy()                # <-- the tautology-free subset for everything below
        print(f"\n--- (2b) ORIENTATION x class, REPRESSIVE couplings only (n={len(d)}) ---")
        ct = pd.crosstab(d.composition_class, d.orientation)
        ct["%same"] = (100 * ct.get("same", 0) / ct.sum(1)).round(1)
        print(ct.to_string())
        ce = d[d.composition_class == "composition_explained"]
        ci = d[d.composition_class == "cell_intrinsic"]
        tab = [[int((ce.orientation == "same").sum()), int((ce.orientation == "opposite").sum())],
               [int((ci.orientation == "same").sum()), int((ci.orientation == "opposite").sum())]]
        orr, p = fisher_exact(tab)
        print(f"\n  Fisher (composition_explained SAME vs cell_intrinsic SAME): OR={orr:.3f}  p={p:.3g}   table={tab}")

        print("\n  FINAL TAG counts (confounded edges only):")
        ft = c.final_tag.value_counts()
        print("    " + ft.to_string().replace("\n", "\n    "))

        # ---- compartment-intrinsic candidates ----
        c = c[c.beta_core > 0]
        cand = c[c.final_tag.str.startswith("compartment_intrinsic")].sort_values("lineage_loss_fraction",
                                                                                  ascending=False)
        cand = cand[cand.retention_reliable]
        print(f"\n  COMPARTMENT-INTRINSIC candidates (beta>0, same-compartment, single-lineage, reliable): {len(cand)}")
        print("  LOADING-MAGNITUDE gate ladder (a sign() of a ~0 loading is noise; n=1041 => |rho|>0.061 is p<.05):")
        for t in (0.0, 0.05, 0.10, 0.15):
            print(f"    both |rho| > {t:.2f} : {int(((cand.mir_rho.abs()>t)&(cand.tgt_rho.abs()>t)).sum())} survive")
        cols = ["gene", "unit", "beta_core", "retention", "pinpointed_lineage", "lineage_loss_fraction",
                "explained_share", "mir_rho", "tgt_rho"]
        print(cand[cols].head(15).round(3).to_string(index=False))

        # ---- CONFOUND CHECKS: abundance & detection ----
        print("\n  CONFOUND CHECK (abundance / detection are the two dominant confounds here):")
        if level == "arm":
            from mirna_hallmark.learned import data as LD
            X = LD._load()["X"]
            det = pd.Series(np.isfinite(X.to_numpy(float)).mean(1), index=X.index)
            ab = pd.Series(np.nanmean(X.to_numpy(float), 1), index=X.index)
            c2 = c.copy()
            c2["det"] = c2.unit.map(det)
            c2["abund"] = c2.unit.map(ab)
            for v in ("det", "abund"):
                a_ = c2.loc[c2.single_lineage, v].dropna()
                b_ = c2.loc[~c2.single_lineage, v].dropna()
                if len(a_) > 10 and len(b_) > 10:
                    u = mannwhitneyu(a_, b_)
                    print(f"    {v:6s} pinpointed {a_.median():.3f} vs MIX {b_.median():.3f}  MWU p={u.pvalue:.3g}")

    # ---- (4) GENE ROLL-UP ----
    g = pd.read_csv(OUT / "composition_pinpoint_genes.tsv", sep="\t")
    for level in ("family", "arm"):
        gl = g[g.level == level]
        sec(f"(4) GENE ROLL-UP — {level.upper()}  ({len(gl)} genes)   gene_retention = sum(beta_deconv)/sum(beta)")
        print("  gene_composition_class:")
        print("    " + gl.gene_composition_class.value_counts().to_string().replace("\n", "\n    "))
        rel = gl[gl.gene_reliable]
        print(f"\n  reliable genes (|sum b|>1e-6, denom_coherence>0.5, >=1 identified unit): {len(rel)}/{len(gl)}")
        print(f"  median gene_retention (reliable): {rel.gene_retention.median():.3f}")
        conf = rel[rel.gene_composition_class.isin(CONF)]
        pin = conf[conf.gene_single_lineage]
        print(f"\n  confounded genes: {len(conf)}   pinpointed to ONE lineage: {len(pin)} ({100*len(pin)/max(len(conf),1):.1f}%)"
              f"   MIX: {len(conf)-len(pin)}")
        print("  gene pinpointed lineage:")
        vc = pin.gene_pinpointed_lineage.value_counts()
        for L in LIN:
            print(f"    {L:20s} {vc.get(L,0):5d}  ({100*vc.get(L,0)/max(len(pin),1):5.1f}%)")
        print("\n  top-12 pinpointed genes by gene_lineage_loss:")
        print(pin.nlargest(12, "gene_lineage_loss")[["gene", "n_units", "gene_retention",
                                                     "gene_pinpointed_lineage", "gene_lineage_loss",
                                                     "gene_explained_share", "gene_max_null_loss"]]
              .round(3).to_string(index=False))

    sec("FAMILY vs ARM concordance of the pinpoint")
    f = load("family")[["gene", "unit", "pinpointed_lineage", "single_lineage", "composition_class", "arms"]]
    a = load("arm")[["gene", "unit", "pinpointed_lineage", "single_lineage", "composition_class"]]
    f1 = f[f.arms.astype(str).str.count(";").eq(0)].copy()          # singleton families == the arm itself
    f1["arm"] = f1.arms
    m = f1.merge(a, left_on=["gene", "arm"], right_on=["gene", "unit"], suffixes=("_fam", "_arm"))
    cc = m[m.composition_class_fam.isin(CONF) | m.composition_class_arm.isin(CONF)]
    print(f"  singleton-family couplings matched to their arm: {len(m)} (confounded in either: {len(cc)})")
    print(f"  pinpointed lineage AGREES: {100*(cc.pinpointed_lineage_fam == cc.pinpointed_lineage_arm).mean():.1f}%")
    print(f"  single_lineage flag agrees: {100*(cc.single_lineage_fam == cc.single_lineage_arm).mean():.1f}%")

    mf = load("family"); ma = load("arm")
    multi = ma[ma.seed_family.map(ma.groupby("seed_family").size().to_dict()).gt(1)] if "seed_family" in ma else None
    print(f"\n  ARM-only rows (multi-member families, invisible at family level): "
          f"{len(ma) - len(mf)} extra couplings")


if __name__ == "__main__":
    main()
