"""Universal 3'UTR seed-scan — site counts + positions for EVERY (miRNA, gene) HE edge.

Fixes TargetScan's coverage gap (only 321 conserved-family arms): scans the MANE-Select 3'UTR of each
HE-edge gene for each miRNA's seed, so we get site **counts AND positions** for all arms. Positions also
enable the S0/S1/S2 site-pair classification (`site_interaction_schematic.py`).

Pipeline: GENCODE annotation (MANE 3'UTR exon coords) + GRCh38 genome FASTA → 3'UTR sequence;
mature.fa → per-arm seed; canonical seed-match scan (6mer / 7mer-A1 / 7mer-m8 / 8mer, Bartel/TargetScan).

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.utr_seed_scan`
Out: `utr_seed_sites.tsv.gz` (miRNA, gene, utr_len, n_6mer, n_7mer_plus, n_8mer, site_pos)
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import pandas as pd

from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges

OUT = Path(__file__).parent
_GENOME = Path("data/genome_assembly/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
_ANNOT = Path("data/gencode.v49.annotation.parquet")
_MATURE = Path("data/miRNA/mature.fa")
_SUFFIX = re.compile(r"\.\d+$")
_COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ---- genome access via .fai (no pysam) ------------------------------------- #
@lru_cache(maxsize=1)
def _fai() -> dict:
    fai = {}
    for line in open(str(_GENOME) + ".fai"):
        p = line.rstrip("\n").split("\t")
        fai[p[0]] = (int(p[2]), int(p[3]), int(p[4]))   # offset, linebases, linewidth
    return fai


_FH = None
def _fetch(chrom: str, start: int, end: int) -> str:
    """Genomic [start,end] 1-based inclusive, + strand. Maps gencode 'chrN' -> Ensembl 'N'."""
    global _FH
    if _FH is None:
        _FH = open(_GENOME, "rb")
    name = chrom[3:] if chrom.startswith("chr") else chrom
    if name == "M":
        name = "MT"
    fai = _fai()
    if name not in fai:
        return ""
    offset, lb, lw = fai[name]
    bp = lambda p: offset + (p // lb) * lw + (p % lb)
    _FH.seek(bp(start - 1))
    return _FH.read(bp(end) - bp(start - 1)).decode().replace("\n", "").replace("\r", "").upper()


# ---- 3'UTR sequences (MANE) ------------------------------------------------ #
def _gene_3utr(genes: set) -> dict:
    cols = ["feature", "seqname", "start", "end", "strand", "gene_name", "is_mane_select"]
    g = pd.read_parquet(_ANNOT, columns=cols)
    g = g[(g["is_mane_select"]) & (g["gene_name"].isin(genes)) & (g["feature"].isin(["UTR", "stop_codon"]))]
    out = {}
    for gene, sub in g.groupby("gene_name"):
        st = sub[sub["feature"] == "stop_codon"]
        utr = sub[sub["feature"] == "UTR"]
        if st.empty or utr.empty:
            continue
        strand, chrom = sub["strand"].iloc[0], sub["seqname"].iloc[0]
        if strand == "+":
            u3 = utr[utr["start"] >= int(st["start"].min())].sort_values("start")
            seq = "".join(_fetch(chrom, int(r.start), int(r.end)) for r in u3.itertuples())
        else:
            u3 = utr[utr["end"] <= int(st["end"].max())].sort_values("start", ascending=False)
            seq = "".join(_revcomp(_fetch(chrom, int(r.start), int(r.end))) for r in u3.itertuples())
        if seq:
            out[gene] = seq
    return out


@lru_cache(maxsize=1)
def _arm_seq() -> dict:
    out, name, buf = {}, None, []
    for line in open(_MATURE):
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf).upper().replace("U", "T")
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf).upper().replace("U", "T")
    return {k: v for k, v in out.items() if k.startswith("hsa-")}


def _scan(utr: str, m: str) -> tuple:
    """Return (n_6mer, n_7mer_plus, n_8mer, [functional-site UTR positions]) for mature DNA seq m."""
    if len(m) < 8 or not utr:
        return 0, 0, 0, []
    core6 = _revcomp(m[1:7])                 # mRNA 6mer site (miRNA nt2-7)
    c_m8 = _COMP_CHAR(m[7])                   # complement of miRNA nt8
    n6 = n7 = n8 = 0; pos = []
    i = utr.find(core6)
    while i != -1:
        m8 = i >= 1 and utr[i - 1] == c_m8                     # 5' of 6mer pairs miRNA nt8
        a1 = i + 6 < len(utr) and utr[i + 6] == "A"            # 3' of 6mer = A1 rule
        n6 += 1
        if m8 and a1:
            n8 += 1
        if m8 or a1:
            n7 += 1; pos.append(i)
        i = utr.find(core6, i + 1)
    return n6, n7, n8, pos


def _COMP_CHAR(b: str) -> str:
    return {"A": "T", "C": "G", "G": "C", "T": "A"}.get(b, "N")


def build(out_dir: Path = OUT) -> pd.DataFrame:
    hs = HallmarkSets.load()
    edges = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    genes = set(edges["gene"])
    utrs = _gene_3utr(genes)
    arms = _arm_seq()
    print(f"[utr-scan] 3'UTRs extracted for {len(utrs)}/{len(genes)} genes (MANE); {len(edges)} HE edges to scan")

    seqcache = {}
    rows = []
    for r in edges.itertuples():
        utr = utrs.get(r.gene)
        m = arms.get(r.miRNA) or arms.get(_SUFFIX.sub("", r.miRNA))
        if utr is None or m is None:
            rows.append((r.miRNA, r.gene, len(utr) if utr else 0, 0, 0, 0, "", utr is None, m is None))
            continue
        key = (r.miRNA, r.gene)
        if key not in seqcache:
            seqcache[key] = _scan(utr, m)
        n6, n7, n8, pos = seqcache[key]
        rows.append((r.miRNA, r.gene, len(utr), n6, n7, n8, ",".join(map(str, pos)), False, False))
    df = pd.DataFrame(rows, columns=["miRNA", "gene", "utr_len", "n_6mer", "n_7mer_plus", "n_8mer",
                                     "site_pos", "no_utr", "no_seq"])
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "utr_seed_sites.tsv.gz"
    df.to_csv(path, sep="\t", index=False)

    ok = df[~(df["no_utr"] | df["no_seq"])]
    print(f"[utr-scan] scanned {len(ok)}/{len(df)} edges (no UTR: {int(df['no_utr'].sum())}, no seq: {int(df['no_seq'].sum())})")
    print(f"           edges with ≥1 functional (7mer+) site: {int((ok['n_7mer_plus']>=1).sum())} "
          f"({(ok['n_7mer_plus']>=1).mean():.0%}); multi-site (≥2): {int((ok['n_7mer_plus']>=2).sum())}; "
          f"max sites {int(ok['n_7mer_plus'].max())}")
    print(f"           vs TargetScan coverage (321 arms): seed-scan covers {ok['miRNA'].nunique()} arms")
    print(f"[utr-scan] wrote {path}")
    return df


if __name__ == "__main__":
    build()
