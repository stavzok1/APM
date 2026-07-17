"""L5 site-level occupancy — lift our MANE-3'UTR-relative predicted sites to GRCh38 genomic coordinates
and intersect with Manakov chimeric-eCLIP duplexes (experimentally observed miRNA-target sites).

Coordinate systems verified GRCh38 both sides: GENCODE v49 UTR exons (gencode 'chrN') ↔ Manakov chr.g
('N', no prefix) — Manakov's PTEN 3'UTR duplexes (chr10:87,965,469…) fall inside PTEN's GENCODE 3'UTR.
A predicted site is `site_clip` if it overlaps a Manakov duplex of the **same miRNA in the same gene**.

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.site_genomic_l5`
Out: adds `chrom,g_start,g_end,site_clip` to `utr_site_ladder_genomic.tsv.gz`.
"""
from __future__ import annotations

import bisect
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges

OUT = Path(__file__).parent
_ANNOT = Path("data/gencode.v49.annotation.parquet")
_MANAKOV = "data/external_cache/manakov_chimeric/*.tsv"
_TARBASE = Path("data/miRNA/Homo_sapiens_TarBase-v9.tsv.gz")
# POSTAR AGO peaks — the OPTIONAL, gene-agnostic (non-arm-resolving) secondary source. The arm-resolving
# sources (Manakov chimeric, TarBase CLIP) are the ones that matter. Prefer a persisted data location; the
# dead-session /tmp path is a legacy fallback (MH-149). If neither exists, POSTAR is skipped (module still runs).
_POSTAR_CANDIDATES = [
    Path("data/external_cache/postar_ago.bed"),
    Path("/tmp/claude-207054/-sci-labs-michall-stavzok-APM/e8a8c248-2ca1-467d-ae0d-4b23f52cc38e/scratchpad/postar_ago.bed"),
]
_POSTAR_AGO = next((p for p in _POSTAR_CANDIDATES if p.exists()), _POSTAR_CANDIDATES[0])
_TB_CLIP = {"HITS-CLIP", "PAR-CLIP", "CLASH", "qCLASH", "Chimeric fragments"}   # site-bearing CLIP methods
_strip = lambda a: re.sub(r"\.\d+$", "", str(a))


def _utr_blocks(genes: set) -> dict:
    """gene -> (chrom_nochr, strand, [(gstart,gend)...] in transcript 5'→3' order) for the MANE 3'UTR."""
    g = pd.read_parquet(_ANNOT, columns=["feature", "seqname", "start", "end", "strand",
                                         "gene_name", "is_mane_select"])
    g = g[(g["is_mane_select"]) & (g["gene_name"].isin(genes)) & (g["feature"].isin(["UTR", "stop_codon"]))]
    out = {}
    for gene, sub in g.groupby("gene_name"):
        st = sub[sub["feature"] == "stop_codon"]; utr = sub[sub["feature"] == "UTR"]
        if st.empty or utr.empty:
            continue
        strand = sub["strand"].iloc[0]; chrom = str(sub["seqname"].iloc[0]).replace("chr", "")
        if strand == "+":
            u3 = utr[utr["start"] >= int(st["start"].min())].sort_values("start")          # genomic asc = tx order
        else:
            u3 = utr[utr["end"] <= int(st["end"].max())].sort_values("start", ascending=False)  # high→low = tx order
        blocks = [(int(r.start), int(r.end)) for r in u3.itertuples()]
        if blocks:
            out[gene] = (chrom, strand, blocks)
    return out


def _lift(blocks_entry, pos: int, length: int = 7):
    """UTR-relative pos (0-based, mRNA 5'→3') -> (g0, g1) genomic, using exon blocks in transcript order."""
    chrom, strand, blocks = blocks_entry
    cum = 0
    for gs, ge in blocks:
        blen = ge - gs + 1
        if pos < cum + blen:
            off = pos - cum
            if strand == "+":
                g0 = gs + off
            else:
                g0 = ge - off - (length - 1)
            return chrom, g0, g0 + length - 1
        cum += blen
    return None, None, None


def _manakov_intervals(genes: set) -> dict:
    """(miRNA, gene) -> list of (start,end) Manakov 3'UTR duplex intervals (GRCh38)."""
    iv = {}
    for fp in glob.glob(_MANAKOV):
        m = pd.read_csv(fp, sep="\t", low_memory=False,
                        usecols=["chr.g", "start.g", "end.g", "gene_name", "noncodingRNA",
                                 "noncodingRNA_type", "feature"])
        m = m[m["gene_name"].isin(genes) & m["noncodingRNA_type"].astype(str).str.contains("miR", case=False, na=False)
              & m["feature"].astype(str).str.contains("three_prime_utr", na=False)]
        for a, gn, s, en in zip(m["noncodingRNA"], m["gene_name"], m["start.g"], m["end.g"]):
            iv.setdefault((_strip(a), gn), []).append((int(s), int(en)))
    return iv


def _tarbase_intervals(genes: set) -> dict:
    """(miRNA, gene) -> list of (start,end) TarBase v9 **site-level CLIP** intervals (GRCh38), our genes only."""
    iv = {}
    for ch in pd.read_csv(_TARBASE, sep="\t", low_memory=False, chunksize=1_000_000,
                          usecols=["mirna_name", "gene_name", "start", "end", "experimental_method"]):
        ch = ch[ch["gene_name"].isin(genes) & ch["experimental_method"].isin(_TB_CLIP)].dropna(subset=["start", "end"])
        for a, gn, s, en in zip(ch["mirna_name"], ch["gene_name"], ch["start"], ch["end"]):
            iv.setdefault((_strip(a), gn), []).append((int(s), int(en)))
    return iv


def _postar_index() -> dict:
    """chrom -> (starts ndarray sorted, ends ndarray) for POSTAR AGO peaks (gene-agnostic occupancy).
    OPTIONAL: returns {} if the bed is unavailable — the module runs on Manakov + TarBase alone."""
    if not _POSTAR_AGO.exists():
        print(f"[L5] POSTAR AGO bed absent ({_POSTAR_AGO}) — skipping the gene-agnostic occupancy source.")
        return {}
    bed = pd.read_csv(_POSTAR_AGO, sep="\t", header=None, names=["chrom", "start", "end"])
    idx = {}
    for c, sub in bed.groupby("chrom"):
        sub = sub.sort_values("start")
        idx[c] = (sub["start"].to_numpy(), sub["end"].to_numpy())
    return idx


def _postar_hit(idx, chrom, g0, g1) -> bool:
    if chrom not in idx:
        return False
    starts, ends = idx[chrom]
    lo = bisect.bisect_left(starts, g0 - 500); hi = bisect.bisect_right(starts, g1)
    for k in range(lo, hi):
        if ends[k] >= g0 and starts[k] <= g1:
            return True
    return False


def build() -> pd.DataFrame:
    d = pd.read_csv(OUT / "utr_site_ladder.tsv.gz", sep="\t")
    genes = set(d["gene"])
    blocks = _utr_blocks(genes)
    man = _manakov_intervals(genes)
    tb = _tarbase_intervals(genes)
    post = _postar_index()
    print(f"[L5] {len(blocks)}/{len(genes)} genes lifted; Manakov pairs {len(man):,}; "
          f"TarBase-CLIP pairs {len(tb):,}; POSTAR AGO peaks {sum(len(v[0]) for v in post.values()):,}")

    chrom, g0, g1, s_man, s_tb, s_post = [], [], [], [], [], []
    for r in d.itertuples():
        be = blocks.get(r.gene)
        if be is None:
            chrom.append(None); g0.append(None); g1.append(None)
            s_man.append(False); s_tb.append(False); s_post.append(False); continue
        c, a, b = _lift(be, int(r.utr_pos)); cc = "chr" + c if c else None
        chrom.append(cc); g0.append(a); g1.append(b)
        key = (_strip(r.miRNA), r.gene)
        s_man.append(any(a <= e and b >= s for s, e in man.get(key, [])))        # miRNA-specific duplex
        s_tb.append(any(a <= e and b >= s for s, e in tb.get(key, [])))          # miRNA-specific TarBase CLIP
        s_post.append(_postar_hit(post, cc, a, b) if a is not None else False)   # AGO occupancy (any miRNA)
    d["chrom"], d["g_start"], d["g_end"] = chrom, g0, g1
    d["site_manakov"], d["site_tarbase"], d["site_postar_ago"] = s_man, s_tb, s_post
    d["site_clip_any"] = d["site_manakov"] | d["site_tarbase"] | d["site_postar_ago"]
    path = OUT / "utr_site_ladder_genomic.tsv.gz"
    d.to_csv(path, sep="\t", index=False)

    print(f"[L5] lifted {int(d['g_start'].notna().sum()):,}/{len(d):,} sites; site-level support — "
          f"Manakov {d['site_manakov'].mean():.0%}, TarBase-CLIP {d['site_tarbase'].mean():.0%}, "
          f"POSTAR-AGO {d['site_postar_ago'].mean():.0%}, any {d['site_clip_any'].mean():.0%}")
    print("  site-level occupancy by rung (predicted site coincides with experimental binding):")
    print(f"    {'rung':24s} {'n':>5s} {'Manakov':>7s} {'TarBase':>7s} {'POSTAR':>7s} {'any':>5s}")
    for lab, mk in [("7mer-A1", d["type"] == "7mer-A1"), ("7mer-m8", d["type"] == "7mer-m8"),
                    ("8mer", d["type"] == "8mer"), ("conserved", d["conserved"]),
                    ("8mer+conserved+3'-supp", (d["type"] == "8mer") & d["conserved"] & d["has_3p_supp"])]:
        s = d[mk]
        print(f"    {lab:24s} {len(s):>5d} {s['site_manakov'].mean():>6.0%} {s['site_tarbase'].mean():>6.0%} "
              f"{s['site_postar_ago'].mean():>6.0%} {s['site_clip_any'].mean():>4.0%}")
    print(f"[L5] wrote {path}")
    return d


if __name__ == "__main__":
    build()
