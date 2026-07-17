"""Experimental validation of the site filter ladder — does climbing the ladder predict experimental
miRNA-target support? Three independent experimental layers:

  ENCORI  AGO-CLIP interaction support     (in repo)
  TarBase v9  experimentally-validated interactions (4.7M rows)  + HE union/intersection
  Manakov chimeric eCLIP  direct miRNA-target duplexes (HEK293T)

For each, P(experimentally-supported | ladder feature) should rise with rung. Also reports TarBase×HE
union/intersection. Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.validate_ladder_experimental`
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import pandas as pd

from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges

OUT = Path(__file__).parent
_TARBASE = Path("data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz")
_MANAKOV = "data/external_cache/manakov_chimeric/*.tsv"
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))


def _tarbase_pairs():
    tb = pd.read_csv(_TARBASE, sep="\t", usecols=["mirna_name", "gene_name", "tissue"], low_memory=False)
    allp = {(_strip(m), g) for m, g in zip(tb["mirna_name"], tb["gene_name"])}
    br = tb[tb["tissue"].astype(str).str.contains("reast", na=False)]
    breast = {(_strip(m), g) for m, g in zip(br["mirna_name"], br["gene_name"])}
    return tb, allp, breast


def _manakov_pairs():
    pairs = set()
    for fp in glob.glob(_MANAKOV):
        m = pd.read_csv(fp, sep="\t", usecols=["gene_name", "noncodingRNA", "noncodingRNA_type"], low_memory=False)
        m = m[m["noncodingRNA_type"].astype(str).str.contains("miR", case=False, na=False)]
        pairs |= {(_strip(a), g) for a, g in zip(m["noncodingRNA"], m["gene_name"])}
    return pairs


def build() -> None:
    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    egenes = set(e["gene"])
    he = {(_strip(m), g) for m, g in zip(e["miRNA"], e["gene"])}

    tb, tb_all, tb_breast = _tarbase_pairs()
    tb_he = {p for p in tb_all if p[1] in egenes}                 # restrict to our genes for union
    man = _manakov_pairs()
    try:
        from mirna_hallmark.encori_edges import enc_depth_lookup
        enc = enc_depth_lookup(genes=sorted(egenes))
        encori = {(_strip(m), g) for m, g in zip(enc["miRNA"], enc["gene"])} if not enc.empty else set()
    except Exception:
        encori = set()

    print(f"TarBase v9: {len(tb):,} rows → {len(tb_all):,} distinct (miRNA,gene); breast {len(tb_breast):,} pairs")
    print(f"Manakov chimeric: {len(man):,} (miRNA,gene) duplex pairs")
    print(f"\nHE edges = {len(he):,}")
    print(f"  HE ∩ TarBase (HE edges experimentally validated): {len(he & tb_he):,} ({len(he & tb_he)/len(he):.0%})")
    print(f"  HE ∩ Manakov (direct duplex):                     {len(he & man):,} ({len(he & man)/len(he):.0%})")
    print(f"  HE ∩ ENCORI (AGO-CLIP):                           {len(he & encori):,} ({len(he & encori)/len(he):.0%})")
    print(f"  HE only (no TarBase): {len(he - tb_he):,}   |   HE ∪ TarBase(HE genes) = {len(he | tb_he):,}")

    d = pd.read_csv(OUT / "utr_site_ladder.tsv.gz", sep="\t")
    g = d.groupby(["miRNA", "gene"]).agg(has_8mer=("type", lambda s: (s == "8mer").any()),
                                         conserved=("conserved", "any"),
                                         has_3psupp=("has_3p_supp", "any")).reset_index()
    g["key"] = list(zip(g["miRNA"].map(_strip), g["gene"]))
    for name, ref in [("ENCORI", encori), ("TarBase", tb_he), ("Manakov", man)]:
        g[name] = g["key"].isin(ref)
    g.drop(columns="key").to_csv(OUT / "ladder_experimental_support.tsv.gz", sep="\t", index=False)

    print("\nP(experimentally supported | ladder feature)  — should rise up the ladder:")
    print(f"  {'feature':30s} {'n':>5s} {'ENCORI':>7s} {'TarBase':>8s} {'Manakov':>8s}")
    masks = [("weak (none)", ~g["has_8mer"] & ~g["conserved"] & ~g["has_3psupp"]),
             ("8mer (L1)", g["has_8mer"]), ("3'-supp (L3)", g["has_3psupp"]),
             ("conserved PCT≥0.5 (L4)", g["conserved"]),
             ("8mer+conserved+3'-supp (top)", g["has_8mer"] & g["conserved"] & g["has_3psupp"])]
    for lab, mk in masks:
        s = g[mk]
        print(f"  {lab:30s} {len(s):>5d} {s['ENCORI'].mean():>6.0%} {s['TarBase'].mean():>7.0%} {s['Manakov'].mean():>7.0%}")
    print(f"\n[validate] wrote ladder_experimental_support.tsv.gz")


if __name__ == "__main__":
    build()
