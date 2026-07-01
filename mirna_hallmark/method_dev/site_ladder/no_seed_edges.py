"""Characterise the HE edges whose miRNA seed matches NOWHERE on the target's MANE 3'UTR — and check whether
a non-MANE **3'UTR isoform** rescues them.

The universal seed-scan finds a canonical 6/7/8mer for 4,143/5,219 HE edges; the remaining **1,076 have NO
canonical seed on the MANE 3'UTR**. Two questions: (1) are those real (experimental support)? (2) does an
alternative 3'UTR isoform (from GENCODE, all transcripts) carry the missing site?

Panels:
  A  experimental support (ENCORI/TarBase/Manakov/any) — no-seed vs sited edges.
  E  3'UTR-isoform rescue — of the no-seed edges, how many gain a canonical site in a NON-MANE isoform.
  F  experimental support (any) by isoform-rescue status — does an isoform site predict real support?
  B  target MANE-3'UTR length — no-seed on shorter UTRs.
  C / D  miRNAs / genes with the most no-seed edges.

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.no_seed_edges`  (evidence + isoform scan cached)
Out: `figures/no_seed_edges.png` (+ `no_seed_edges.tsv.gz`)
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mirna_hallmark.method_dev.site_ladder.utr_seed_scan import _ANNOT, _SUFFIX, _arm_seq, _fetch, _revcomp, _scan

HERE = Path(__file__).parent
FIG = Path(__file__).parent.parent / "figures"
_TARBASE = Path("data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz")
_MANAKOV = "data/external_cache/manakov_chimeric/*.tsv"
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))
_REG = ["iso3UTR", "CDS", "5UTR", "iso6", "none"]
_REG_C = {"iso3UTR": "#1b7837", "CDS": "#3182bd", "5UTR": "#8856a7", "iso6": "#f1a340", "none": "#969696"}
_REG_LAB = {"iso3UTR": "7mer+ in\nalt 3′UTR iso", "CDS": "7mer+ in\nCDS (ORF)", "5UTR": "7mer+ in\n5′UTR",
            "iso6": "6mer only\nin isoform", "none": "seedless\neverywhere"}
_REG_SHORT = {"iso3UTR": "alt 3′UTR", "CDS": "CDS", "5UTR": "5′UTR", "iso6": "iso 6mer", "none": "seedless"}


def _pairsets(genes):
    tb = pd.read_csv(_TARBASE, sep="\t", usecols=["mirna_name", "gene_name"], low_memory=False)
    tbp = {(_strip(m), g) for m, g in zip(tb["mirna_name"], tb["gene_name"])}
    man = set()
    for fp in glob.glob(_MANAKOV):
        m = pd.read_csv(fp, sep="\t", usecols=["gene_name", "noncodingRNA", "noncodingRNA_type"], low_memory=False)
        m = m[m["noncodingRNA_type"].astype(str).str.contains("miR", case=False, na=False)]
        man |= {(_strip(a), g) for a, g in zip(m["noncodingRNA"], m["gene_name"])}
    try:
        from mirna_hallmark.encori_edges import enc_depth_lookup
        enc = enc_depth_lookup(genes=sorted(genes))
        encp = {(_strip(m), g) for m, g in zip(enc["miRNA"], enc["gene"])}
    except Exception:
        encp = set()
    return tbp, man, encp


def _isoform_3utrs(genes: set) -> dict:
    """gene -> list of NON-MANE 3'UTR isoform sequences (per transcript, strand-aware)."""
    cols = ["feature", "seqname", "start", "end", "strand", "gene_name", "transcript_id", "is_mane_select"]
    g = pd.read_parquet(_ANNOT, columns=cols)
    g = g[g["gene_name"].isin(genes) & g["feature"].isin(["UTR", "stop_codon"])]
    out = {}
    for gene, gs in g.groupby("gene_name"):
        seqs = []
        for _tx, sub in gs.groupby("transcript_id"):
            st, utr = sub[sub["feature"] == "stop_codon"], sub[sub["feature"] == "UTR"]
            if st.empty or utr.empty or bool(sub["is_mane_select"].any()):     # skip non-coding + MANE (already scanned)
                continue
            strand, chrom = sub["strand"].iloc[0], sub["seqname"].iloc[0]
            if strand == "+":
                u3 = utr[utr["start"] >= int(st["start"].min())].sort_values("start")
                seq = "".join(_fetch(chrom, int(r.start), int(r.end)) for r in u3.itertuples())
            else:
                u3 = utr[utr["end"] <= int(st["end"].max())].sort_values("start", ascending=False)
                seq = "".join(_revcomp(_fetch(chrom, int(r.start), int(r.end))) for r in u3.itertuples())
            if seq:
                seqs.append(seq)
        if seqs:
            out[gene] = seqs
    return out


def _mane_cds_5utr(genes: set) -> dict:
    """gene -> (MANE CDS seq, MANE 5'UTR seq), strand-aware."""
    cols = ["feature", "seqname", "start", "end", "strand", "gene_name", "is_mane_select"]
    g = pd.read_parquet(_ANNOT, columns=cols)
    g = g[g["is_mane_select"] & g["gene_name"].isin(genes) & g["feature"].isin(["CDS", "UTR", "start_codon"])]
    out = {}
    for gene, sub in g.groupby("gene_name"):
        sc, cds, utr = sub[sub["feature"] == "start_codon"], sub[sub["feature"] == "CDS"], sub[sub["feature"] == "UTR"]
        if sc.empty or cds.empty:
            out[gene] = ("", ""); continue
        strand, chrom = sub["strand"].iloc[0], sub["seqname"].iloc[0]
        if strand == "+":
            cseq = "".join(_fetch(chrom, int(r.start), int(r.end)) for r in cds.sort_values("start").itertuples())
            u5 = utr[utr["end"] <= int(sc["start"].min())].sort_values("start")
            useq = "".join(_fetch(chrom, int(r.start), int(r.end)) for r in u5.itertuples())
        else:
            cseq = "".join(_revcomp(_fetch(chrom, int(r.start), int(r.end))) for r in cds.sort_values("start", ascending=False).itertuples())
            u5 = utr[utr["start"] >= int(sc["end"].max())].sort_values("start", ascending=False)
            useq = "".join(_revcomp(_fetch(chrom, int(r.start), int(r.end))) for r in u5.itertuples())
        out[gene] = (cseq, useq)
    return out


def _region_rescue(noc: pd.DataFrame, force=False) -> pd.DataFrame:
    """Where does the missing canonical seed sit, if anywhere? priority alt-3'UTR-iso > CDS > 5'UTR > iso-6mer > none."""
    cache = HERE / "_noseed_region.tsv.gz"
    if cache.exists() and not force:
        return pd.read_csv(cache, sep="\t")
    iso, cds5, arms = _isoform_3utrs(set(noc["gene"])), _mane_cds_5utr(set(noc["gene"])), _arm_seq()
    rec = []
    for e in noc.itertuples():
        m = arms.get(e.miRNA) or arms.get(_SUFFIX.sub("", e.miRNA))
        i7 = i6 = c7 = u7 = False
        if m:
            for seq in iso.get(e.gene, []):
                n6, n7, _n8, _p = _scan(seq, m); i7 = i7 or n7 >= 1; i6 = i6 or n6 >= 1
            cseq, useq = cds5.get(e.gene, ("", ""))
            c7, u7 = _scan(cseq, m)[1] >= 1, _scan(useq, m)[1] >= 1
        region = "iso3UTR" if i7 else ("CDS" if c7 else ("5UTR" if u7 else ("iso6" if i6 else "none")))
        rec.append((e.miRNA, e.gene, region))
    out = pd.DataFrame(rec, columns=["miRNA", "gene", "region"])
    out.to_csv(cache, sep="\t", index=False)
    return out


def build(force=False):
    ss = pd.read_csv(HERE / "utr_seed_sites.tsv.gz", sep="\t")
    noc = ss[(ss["n_7mer_plus"] == 0) & (ss["n_6mer"] == 0)]
    sited = ss[ss["n_7mer_plus"] > 0]
    cache = HERE / "_noseed_evidence.tsv.gz"
    if cache.exists() and not force:
        ev = pd.read_csv(cache, sep="\t")
    else:
        tbp, man, encp = _pairsets(set(ss["gene"]))
        def fl(df, grp):
            k = list(zip(df["miRNA"].map(_strip), df["gene"]))
            return pd.DataFrame({"miRNA": df["miRNA"].values, "gene": df["gene"].values, "utr_len": df["utr_len"].values,
                                 "grp": grp, "TarBase": [p in tbp for p in k], "ENCORI": [p in encp for p in k],
                                 "Manakov": [p in man for p in k]})
        ev = pd.concat([fl(noc, "noc"), fl(sited, "sited")], ignore_index=True)
        ev["any"] = ev[["TarBase", "ENCORI", "Manakov"]].any(axis=1)
        ev.to_csv(cache, sep="\t", index=False)
    fn, fs = ev[ev["grp"] == "noc"].copy(), ev[ev["grp"] == "sited"]
    fn = fn.merge(_region_rescue(noc, force=force), on=["miRNA", "gene"], how="left")

    fig, axes = plt.subplots(2, 3, figsize=(19, 10))
    (axA, axE, axF), (axB, axC, axD) = axes
    fig.subplots_adjust(hspace=0.55, wspace=0.28, top=0.9, bottom=0.08, left=0.05, right=0.99)

    # A — experimental support
    cats = ["ENCORI", "TarBase", "Manakov", "any"]; x = np.arange(len(cats)); w = 0.38
    nv, sv = [fn[c].mean() for c in cats], [fs[c].mean() for c in cats]
    axA.bar(x - w / 2, nv, w, color="#d7301f", label=f"no-seed ({len(fn):,})")
    axA.bar(x + w / 2, sv, w, color="#3182bd", label=f"sited 7mer+ ({len(fs):,})")
    for i, v in enumerate(nv): axA.text(i - w / 2, v + .01, f"{v:.0%}", ha="center", fontsize=8, color="#d7301f")
    for i, v in enumerate(sv): axA.text(i + w / 2, v + .01, f"{v:.0%}", ha="center", fontsize=8, color="#3182bd")
    axA.set_xticks(x); axA.set_xticklabels(["ENCORI\nAGO-CLIP", "TarBase", "Manakov", "ANY"]); axA.set_ylim(0, max(sv) + 0.12)
    axA.set_ylabel("fraction supported"); axA.legend(fontsize=8)
    axA.set_title(f"A · experimental support — no-seed far weaker\nENCORI AGO-CLIP {nv[0]:.0%} vs {sv[0]:.0%}; {fn['any'].mean():.0%} any support",
                  fontsize=11, loc="left")

    # E — where's the missing site? rescue by region (alt 3'UTR isoform / CDS / 5'UTR)
    rc = fn["region"].value_counts()
    vals = [int(rc.get(k, 0)) for k in _REG]
    axE.bar(range(len(_REG)), vals, color=[_REG_C[k] for k in _REG])
    for i, v in enumerate(vals): axE.text(i, v + 6, f"{v}\n{v/len(fn):.0%}", ha="center", fontsize=8)
    axE.set_xticks(range(len(_REG))); axE.set_xticklabels([_REG_LAB[k] for k in _REG], fontsize=7)
    axE.set_ylabel("# no-seed edges"); axE.set_ylim(0, max(vals) * 1.2)
    axE.set_title(f"E · where's the missing site? region rescue\n{1 - vals[-1]/len(fn):.0%} rescued (alt-isoform/CDS/5′UTR); {vals[-1]/len(fn):.0%} seedless everywhere",
                  fontsize=11, loc="left")

    # F — experimental support (any) by rescue region
    grp = fn.groupby("region")["any"].mean().reindex(_REG)
    axF.bar(range(len(_REG)), grp.values, color=[_REG_C[k] for k in _REG])
    for i, v in enumerate(grp.values):
        if pd.notna(v): axF.text(i, v + .01, f"{v:.0%}", ha="center", fontsize=8.5)
    axF.set_xticks(range(len(_REG))); axF.set_xticklabels([_REG_SHORT[k] for k in _REG], fontsize=8)
    axF.set_ylabel("fraction with ANY support"); axF.set_ylim(0, (grp.max() + 0.1) if grp.notna().any() else 0.6)
    axF.set_title(f"F · does the rescue find real sites?\nalt-3′UTR & CDS rescues better-supported ({grp['CDS']:.0%}) than seedless ({grp['none']:.0%})",
                  fontsize=11, loc="left")

    # B — MANE UTR length
    axB.boxplot([np.clip(fn["utr_len"], 1, None), np.clip(fs["utr_len"], 1, None)], showfliers=False, widths=0.5)
    axB.set_yscale("log"); axB.set_xticks([1, 2]); axB.set_xticklabels(["no-seed", "sited 7mer+"]); axB.set_ylabel("MANE 3'UTR length (nt, log)")
    axB.set_title(f"B · MANE 3′UTR length — no-seed on shorter UTRs\nmedian {fn['utr_len'].median():.0f} vs {fs['utr_len'].median():.0f} nt",
                  fontsize=11, loc="left")

    tm = fn["miRNA"].str.replace("hsa-", "", regex=False).value_counts().head(12)[::-1]
    axC.barh(range(len(tm)), tm.values, color="#756bb1"); axC.set_yticks(range(len(tm))); axC.set_yticklabels(tm.index, fontsize=8)
    axC.set_xlabel("# no-seed edges"); axC.set_title(f"C · miRNAs with most no-seed edges ({fn['miRNA'].nunique()})\ne.g. miR-375/137/320a — atypical seeds", fontsize=11, loc="left")

    tg = fn["gene"].value_counts().head(12)[::-1]
    axD.barh(range(len(tg)), tg.values, color="#31a354"); axD.set_yticks(range(len(tg))); axD.set_yticklabels(tg.index, fontsize=8)
    axD.set_xlabel("# no-seed edges"); axD.set_title(f"D · genes with most no-seed edges ({fn['gene'].nunique()})\nstudied hubs (HIF1A/MYC/DNMT1/PTEN)", fontsize=11, loc="left")

    resc = 1 - vals[-1] / len(fn)
    fig.suptitle(f"HE edges with NO canonical seed on the MANE 3′UTR — {len(fn):,} of {len(ss):,}; "
                 f"{resc:.0%} rescued in another region (alt 3′UTR isoform / CDS / 5′UTR), {vals[-1]/len(fn):.0%} seedless everywhere",
                 fontsize=13, fontweight="bold", y=0.965)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "no_seed_edges.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    fn.to_csv(HERE / "no_seed_edges.tsv.gz", sep="\t", index=False)
    print(f"[no-seed] {len(fn):,} no-MANE-3′UTR-seed edges; region rescue "
          f"{dict(zip(_REG, vals))}; {resc:.0%} rescued / {vals[-1]/len(fn):.0%} seedless everywhere; wrote figures/no_seed_edges.png")


if __name__ == "__main__":
    build()
