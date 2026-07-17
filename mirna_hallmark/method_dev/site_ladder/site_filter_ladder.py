"""miRNA target-site confidence FILTER LADDER (plan: SITE_FILTER_LADDER_PLAN.md).

Per-SITE table over the universal seed-scan, with published filters layered as rungs of rising precision:
  L1 site type (drop 6mer; weight 8mer>7m-m8>7m-A1)        Bartel 2009 / Agarwal 2015
  L2 local AU content + position-in-UTR (drop <15nt of stop) Grimson 2007
  L3 3'-supplementary pairing (miRNA nt 13-16)              Grimson 2007
  L4 conservation PCT (score, not hard filter)              Friedman 2009   [added next]
  APA-aware UTR mask (parallel breast set)                  APAatlas        [added next]
  L5 AGO-CLIP / chimeric occupancy                          POSTAR3/Manakov [needs download]

This file builds L1–L3 (pure sequence); L4/APA/L5 extend it. Reuses 3'UTR + seed machinery from
`utr_seed_scan`. Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.site_filter_ladder`
Out: `utr_site_ladder.tsv.gz` (one row per functional site, 7mer+).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges
from mirna_hallmark.method_dev.site_ladder.utr_seed_scan import _arm_seq, _gene_3utr, _revcomp, _SUFFIX

OUT = Path(__file__).parent
_CONSERVED = Path("data/miRNA/Conserved_Family_Info.txt")
_APA = Path("data/external_cache/apaatlas/breast_pdui_per_gene.tsv")
_APA_DISTAL = 0.67       # site in distal third of UTR + gene APA-shortened in breast -> not retained
_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
# Agarwal 2015 (TargetScan7) mean context++ magnitude by site type (relative repression weight)
_TYPE_W = {"8mer": 0.31, "7mer-m8": 0.161, "7mer-A1": 0.099}
_STOP_BUFFER = 15        # Grimson 2007: sites within 15 nt of stop codon are poorly effective
_AU_WIN = 30             # nt up/down-stream for local-AU
_SUPP_WIN = 17           # upstream window to look for 3'-supplementary pairing


def _local_au(utr: str, s: int, e: int) -> float:
    """Distance-weighted AU fraction in the 30 nt up/down-stream of site [s,e) (Grimson 2007 style)."""
    L = len(utr); num = den = 0.0
    for off in range(1, _AU_WIN + 1):
        w = 1.0 / off
        for pos in (s - off, e - 1 + off):
            if 0 <= pos < L:
                den += w
                if utr[pos] in "AT":
                    num += w
    return num / den if den else 0.0


def _supp_pairs(utr: str, p: int, m: str) -> int:
    """Max contiguous 3'-supplementary match of miRNA nt13-16 (revcomp) in the window 5' of the seed."""
    if len(m) < 16:
        return 0
    motif = _revcomp(m[12:16])              # 4 nt, mRNA-sense, for miRNA nt13-16
    win = utr[max(0, p - _SUPP_WIN):max(0, p - 1)]
    best = 0
    for k in (4, 3):                        # require >=3 contiguous
        if motif[:k] in win or motif[-k:] in win:
            best = max(best, k)
    return best


def _sites(utr: str, m: str) -> list:
    """All 7mer+ sites for mature DNA seq m, each with L1–L3 features."""
    if len(m) < 8 or not utr:
        return []
    core6 = _revcomp(m[1:7]); c_m8 = _COMP.get(m[7], "N"); L = len(utr)
    out = []
    i = utr.find(core6)
    while i != -1:
        m8 = i >= 1 and utr[i - 1] == c_m8
        a1 = i + 6 < L and utr[i + 6] == "A"
        typ = "8mer" if (m8 and a1) else "7mer-m8" if m8 else "7mer-A1" if a1 else "6mer"
        if typ != "6mer":                   # L1: drop 6mers (kept as count only, see utr_seed_scan)
            out.append(dict(
                utr_pos=i, type=typ, type_weight=_TYPE_W[typ],
                local_AU=round(_local_au(utr, i, i + 6), 3),          # L2
                dist_to_stop=i, min_dist=min(i, L - (i + 6)),
                position_ok=bool(i >= _STOP_BUFFER),
                supp_pairs=_supp_pairs(utr, i, m), has_3p_supp=bool(_supp_pairs(utr, i, m) >= 3)))  # L3
        i = utr.find(core6, i + 1)
    return out


def _arm_family() -> dict:
    f = pd.read_csv(C.MIR_FAMILY_INFO, sep="\t")
    f = f[f["Species ID"].astype(str) == "9606"]
    return {r["MiRBase ID"]: r["miR family"] for _, r in f.iterrows()}


def _conserved_pct() -> dict:
    """(miR family, gene) -> max PCT over conserved human sites (Friedman 2009). NULL/missing -> absent."""
    pct = {}
    for ch in pd.read_csv(_CONSERVED, sep="\t", usecols=["miR Family", "Gene Symbol", "Species ID", "PCT"],
                          chunksize=2_000_000, low_memory=False):
        ch = ch[ch["Species ID"].astype(str) == "9606"]
        ch["PCT"] = pd.to_numeric(ch["PCT"], errors="coerce")
        ch = ch.dropna(subset=["PCT"])
        for fam, gene, v in zip(ch["miR Family"], ch["Gene Symbol"], ch["PCT"]):
            k = (fam, gene)
            if v > pct.get(k, -1):
                pct[k] = float(v)
    return pct


def _apa() -> dict:
    """gene -> (apa_shortened in breast, breast PDUI median). APAatlas (Hong 2020)."""
    a = pd.read_csv(_APA, sep="\t")
    return {r["gene_name"]: (bool(r["breast_apa_shortened"]), float(r["breast_pdui_median"]))
            for _, r in a.iterrows()}


def build(out_dir: Path = OUT) -> pd.DataFrame:
    hs = HallmarkSets.load()
    edges = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    utrs = _gene_3utr(set(edges["gene"]))
    arms = _arm_seq()
    fam = _arm_family()
    pct = _conserved_pct()
    apa = _apa()
    rows = []
    for r in edges.itertuples():
        utr = utrs.get(r.gene); m = arms.get(r.miRNA) or arms.get(_SUFFIX.sub("", r.miRNA))
        if utr is None or m is None:
            continue
        p = pct.get((fam.get(r.miRNA, fam.get(_SUFFIX.sub("", r.miRNA))), r.gene))   # L4 edge-level PCT
        sh, pdui = apa.get(r.gene, (False, None))                                     # APA modifier
        for s in _sites(utr, m):
            frac = s["utr_pos"] / len(utr)
            rows.append({"miRNA": r.miRNA, "gene": r.gene, "utr_len": len(utr), **s,
                         "PCT": p, "conserved": bool(p is not None and p >= 0.5),
                         "apa_shortened": sh, "breast_pdui": pdui,
                         "apa_retained": not (sh and frac > _APA_DISTAL)})
    df = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "utr_site_ladder.tsv.gz"
    df.to_csv(path, sep="\t", index=False)

    n_edges = df.groupby(["miRNA", "gene"]).ngroups
    print(f"[ladder] {len(df):,} functional (7mer+) sites over {n_edges:,} edges")
    print("  L1 site types:    " + ", ".join(f"{t}={int((df['type']==t).sum())}" for t in _TYPE_W))
    print(f"  L2 position_ok (≥15nt from stop): {int(df['position_ok'].sum())} ({df['position_ok'].mean():.0%}); "
          f"median local_AU {df['local_AU'].median():.2f}")
    print(f"  L3 3'-supp pairing (≥3nt): {int(df['has_3p_supp'].sum())} ({df['has_3p_supp'].mean():.0%})")
    has_pct = df["PCT"].notna()
    print(f"  L4 conservation: {int(has_pct.sum())} sites have PCT (conserved-family arms); "
          f"conserved (PCT≥0.5): {int(df['conserved'].sum())} ({df['conserved'].mean():.0%}); "
          f"median PCT(where present) {df.loc[has_pct,'PCT'].median():.2f}")
    n_dropped = int((~df['apa_retained']).sum())
    print(f"  APA (breast): {int(df['apa_shortened'].sum())} sites in APA-shortened genes; "
          f"{n_dropped} distal sites dropped from the breast-APA-aware set ({n_dropped/len(df):.0%})")
    print(f"[ladder] wrote {path}")
    return df


if __name__ == "__main__":
    build()
