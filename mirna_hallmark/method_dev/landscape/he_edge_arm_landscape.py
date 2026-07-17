"""Genome-wide landscape of the miRNA arm loci that carry HE (analysis) edges.

Each arm in the analysis edge set (`load_mirtar_edges`, confidence_logclass≥2 + ENCORI; edges are unique
per (miRNA,gene)) placed at its genomic locus, dot size ∝ HE-edge degree (# distinct genes it targets).
Multi-locus miRNAs appear once per locus. Loci from `mirna_mature_loci.csv` matched by name (full +
locus-suffix-stripped) with a `mature.fa` MIMAT-bridge fallback. Non-conserved miRNAs absent from the
(MirGeneDB-style) annotation stay unmapped — their coords live only in the external miRBase gff3.

**Polycistronic clusters** = loci within `CLUSTER_KB` on the same chromosome+strand (≥2 members), drawn
with a red ring. **Per-family** figures (one per seed family) via `--families`.

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.landscape.he_edge_arm_landscape [--families let-7,miR-17 | --families top]`
Out: `figures/he_edge_arm_landscape.png`, `figures/families/<family>.png`, `he_edge_arm_loci.tsv`
"""
from __future__ import annotations

import argparse
import re
from functools import lru_cache
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import load_mirtar_edges

OUT = Path(__file__).parent                  # landscape/ — data outputs (loci.tsv, gene_seed_diversity.tsv)
FIG = Path(__file__).parent.parent / "figures"   # shared figures dir at method_dev root
CLUSTER_KB = 10  # max gap between consecutive same-strand loci to call a polycistronic cluster
HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
    "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
    "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
_CHR_ORDER = list(HG38)
_YMAP = {c: len(_CHR_ORDER) - i for i, c in enumerate(_CHR_ORDER)}
_SUFFIX = re.compile(r"\.\d+$")
_strip = lambda a: _SUFFIX.sub("", str(a))


def _name2mimat() -> dict:
    out = {}
    for line in open(C.MIRNA_MATURE_FA if hasattr(C, "MIRNA_MATURE_FA") else "data/miRNA/mature.fa"):
        if line.startswith(">"):
            p = line[1:].split()
            if p[0].startswith("hsa-") and len(p) > 1 and p[1].startswith("MIMAT"):
                out[p[0]] = p[1]
    return out


def _he_degree() -> pd.Series:
    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    return e.groupby("miRNA")["gene"].nunique().rename("he_degree")


def _tform_arm(a: str) -> str:
    """arm name -> GENCODE miRNA-gene name (hsa-miR-1260a-5p -> MIR1260A; hsa-let-7a-5p -> MIRLET7A)."""
    a = _SUFFIX.sub("", str(a)); a = re.sub(r"^hsa-", "", a); a = re.sub(r"-[35]p$", "", a)
    a = a.upper().replace("MIR-", "MIR").replace("LET-", "LET")
    return "MIR" + a if a.startswith("LET") else a


@lru_cache(maxsize=1)
def _gencode_arm_loci() -> dict:
    """GENCODE-gene-name -> [(chrom,start,strand)] for miRNA genes (recovers non-MirGeneDB arms; in-repo)."""
    g = pd.read_parquet("data/gencode.v49.annotation.parquet",
                        columns=["feature", "seqname", "start", "strand", "gene_name", "gene_type"])
    g = g[(g["feature"] == "gene") & (g["gene_type"] == "miRNA") & (g["seqname"].isin(HG38))]
    out = {}
    for n, c, s, strand in zip(g["gene_name"].astype(str).str.upper(), g["seqname"], g["start"], g["strand"]):
        for key in {n, re.sub(r"-\d+$", "", n)}:
            out.setdefault(key, []).append((c, int(s), strand))
    return out


def mapped_loci() -> pd.DataFrame:
    """arm, he_degree, chrom, start, strand, source. MirGeneDB loci primary; GENCODE miRNA-gene fallback."""
    deg = _he_degree()
    loci = pd.read_csv(C.MIRNA_MATURE_LOCI, usecols=["chrom", "start", "strand", "mirbase_mature_id", "mature_accession"])
    loci = loci[loci["chrom"].isin(HG38)].copy()
    # long key table: name | stripped name | MIMAT  ->  coords
    keyed = pd.concat([
        loci.assign(key=loci["mirbase_mature_id"].astype(str)),
        loci.assign(key=loci["mirbase_mature_id"].map(_strip)),
        loci.assign(key=loci["mature_accession"].astype(str)),
    ])[["key", "chrom", "start", "strand"]].drop_duplicates()
    n2m = _name2mimat()
    am = pd.DataFrame({"arm": deg.index, "he_degree": deg.values})
    am_keyed = pd.concat([
        am.assign(key=am["arm"].astype(str)),
        am.assign(key=am["arm"].map(_strip)),
        am.assign(key=am["arm"].map(lambda a: n2m.get(_strip(a), n2m.get(a, "")))),
    ])
    m = (am_keyed.merge(keyed, on="key").drop_duplicates(["arm", "chrom", "start"]).drop(columns="key"))
    m["start"] = m["start"].astype(int); m["source"] = "mirgenedb"

    # GENCODE miRNA-gene fallback for arms MirGeneDB lacks (non-conserved miRNAs)
    mapped_arms = set(m["arm"]); gmap = _gencode_arm_loci(); extra = []
    for arm, dg in deg.items():
        if arm in mapped_arms:
            continue
        t = _tform_arm(arm)
        locs = gmap.get(t) or gmap.get(re.sub(r"-\d+$", "", t))
        if locs:
            extra += [{"arm": arm, "he_degree": dg, "chrom": c, "start": s, "strand": st, "source": "gencode"}
                      for c, s, st in locs]
    if extra:
        m = pd.concat([m, pd.DataFrame(extra)], ignore_index=True).drop_duplicates(["arm", "chrom", "start"])

    from mirna_hallmark.arm_expression import arm_expression_tiers          # canonical expression tier
    tdict = arm_expression_tiers()["tier"].to_dict()
    m["tier"] = m["arm"].map(lambda a: tdict.get(a, tdict.get(_strip(a), "silent")))

    n_arms, n_map = len(deg), m["arm"].nunique()
    n_gc = m[m["source"] == "gencode"]["arm"].nunique()
    print(f"[landscape] {n_arms} HE-edge arms; {n_map} mapped ({len(m)} loci) — MirGeneDB primary + "
          f"{n_gc} recovered via GENCODE miRNA genes; {n_arms-n_map} still unmapped. "
          f"mapped = {m.drop_duplicates('arm')['he_degree'].sum()/deg.sum():.0%} of total HE-edge degree.")
    return m


def assign_clusters(m: pd.DataFrame) -> pd.DataFrame:
    """Tag loci within CLUSTER_KB on the same chrom+strand as one polycistronic cluster (size≥2)."""
    m = m.sort_values(["chrom", "strand", "start"]).copy()
    cid, sizes = [], {}
    cur, last = -1, None
    for _, r in m.iterrows():
        key = (r["chrom"], r["strand"])
        if last is None or key != last[0] or r["start"] - last[1] > CLUSTER_KB * 1000:
            cur += 1
        cid.append(cur); sizes[cur] = sizes.get(cur, 0) + 1
        last = (key, r["start"])
    m["cluster"] = cid
    m["cluster_size"] = m["cluster"].map(sizes)
    m["clustered"] = m["cluster_size"] >= 2
    return m


# hg38 centromere midpoints (approx; for p/q-arm separation)
CEN = {
    "chr1": 123400000, "chr2": 93900000, "chr3": 90900000, "chr4": 50000000, "chr5": 48800000,
    "chr6": 59800000, "chr7": 60100000, "chr8": 45200000, "chr9": 43000000, "chr10": 39800000,
    "chr11": 53400000, "chr12": 35500000, "chr13": 17700000, "chr14": 17200000, "chr15": 19000000,
    "chr16": 36800000, "chr17": 25100000, "chr18": 18500000, "chr19": 26200000, "chr20": 28100000,
    "chr21": 12000000, "chr22": 15000000, "chrX": 61000000, "chrY": 10400000,
}


def _draw_chroms(ax, m) -> None:
    """p/q-arm-separated chromosome lines (gap + dot at centromere). Per-chromosome miRNA-arm and
    locus counts are written in a clearly-headed column at the right margin (y-ticks = chr number only)."""
    cnt = {c: (sub["arm"].nunique(), len(sub)) for c, sub in m.groupby("chrom")}
    xr = max(HG38.values())
    for c in _CHR_ORDER:
        y, L, cen, gap = _YMAP[c], HG38[c], CEN[c], HG38[c] * 0.004
        ax.plot([0, cen - gap], [y, y], color="0.86", lw=4, solid_capstyle="round", zorder=1)   # p-arm
        ax.plot([cen + gap, L], [y, y], color="0.86", lw=4, solid_capstyle="round", zorder=1)   # q-arm
        ax.plot([cen], [y], marker="o", ms=4.5, color="0.55", zorder=2)                          # centromere
        a, l = cnt.get(c, (0, 0))
        if a:
            ax.text(xr * 1.05, y, f"{a:>3d}     {l:>3d}", fontsize=6.5, va="center", ha="left",
                    color="0.4", fontfamily="monospace")
    ax.text(xr * 1.05, len(_CHR_ORDER) + 1.3, "miRNA  loci\narms", fontsize=6.8, va="center",
            ha="left", color="0.3", fontweight="bold", fontfamily="monospace")
    ax.set_yticks([_YMAP[c] for c in _CHR_ORDER])
    ax.set_yticklabels([c[3:] for c in _CHR_ORDER], fontsize=8)


def _place_labels(ax, items, fontsize=6.2) -> None:
    """Non-overlapping labels: stack vertically per chromosome with leader lines; no two collide in x at a level."""
    span = max(HG38.values())
    dx = span * 0.11
    placed = []  # (x0, x1, level_y)
    for x, chrom, text in sorted(items, key=lambda t: (_CHR_ORDER.index(t[1]), t[0])):
        y0, lvl = _YMAP[chrom], 0
        while lvl < 8:
            ly = round(y0 + 0.32 + lvl * 0.34, 2)
            if not any(abs(ly - py) < 0.01 and not (x + dx / 2 < px0 or x - dx / 2 > px1) for px0, px1, py in placed):
                break
            lvl += 1
        placed.append((x - dx / 2, x + dx / 2, ly))
        ax.annotate(text, (x, y0), xytext=(x, ly), fontsize=fontsize, ha="center", va="bottom",
                    color="0.12", arrowprops=dict(arrowstyle="-", lw=0.35, color="0.6"), zorder=5)


def _scatter(ax, m, vmin, vmax):
    s = 16 + 9 * m["he_degree"].clip(upper=40)
    return ax.scatter(m["start"], m["chrom"].map(_YMAP), s=s, c=m["he_degree"], cmap="viridis",
                      vmin=vmin, vmax=vmax, alpha=0.75, edgecolors="white", linewidths=0.3, zorder=3)


def _finish(fig, ax, sc, title):
    ax.set_xlabel("genomic position (bp, hg38)  ·  ● = centromere (p-arm left, q-arm right)")
    ax.set_ylabel("chromosome")
    ax.set_title(title)
    ax.set_xlim(-max(HG38.values()) * 0.015, max(HG38.values()) * 1.16)
    ax.set_ylim(0.2, len(_CHR_ORDER) + 2.6)
    if sc is not None:
        fig.colorbar(sc, ax=ax, label="HE-edge degree (# target genes)", shrink=0.6)
    fig.tight_layout()


def build_genome() -> pd.DataFrame:
    m = assign_clusters(mapped_loci())
    m.to_csv(OUT / "he_edge_arm_loci.tsv", sep="\t", index=False)
    exp = m[m["tier"] != "silent"]; sil = m[m["tier"] == "silent"]
    fig, ax = plt.subplots(figsize=(13, 9))
    _draw_chroms(ax, m)
    # silent arms (RPM<10 cohort-wide) drawn as small grey hollow rings; expressed coloured by degree
    ax.scatter(sil["start"], sil["chrom"].map(_YMAP), s=12, facecolors="none", edgecolors="0.6",
               linewidths=0.5, zorder=2, label=f"silent ({len(sil['arm'].unique())} arms, never ≥10 RPM in any tumor)")
    sc = _scatter(ax, exp, exp["he_degree"].min(), exp["he_degree"].max())
    top = exp.sort_values("he_degree", ascending=False).head(14)
    _place_labels(ax, [(r.start, r.chrom, f"{r.arm} ({int(r.he_degree)})") for r in top.itertuples()])
    _finish(fig, ax, sc, f"Genome-wide loci of HE-edge miRNA arms — {exp['arm'].nunique()} expressed (filled, "
                         f"colour∝degree) + {sil['arm'].nunique()} silent (grey rings)")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.savefig(FIG / "he_edge_arm_landscape.png", dpi=150)
    print(f"[landscape] {exp['arm'].nunique()} expressed + {sil['arm'].nunique()} silent of "
          f"{m['arm'].nunique()} MAPPED arms shown (unmapped HE arms excluded; canonical tier totals over all "
          f"789 HE arms = 482/307 in arm_expression_tiers.tsv); wrote figure + loci tsv")
    return m


def _arm_family() -> pd.Series:
    fam = pd.read_csv(C.MIR_FAMILY_INFO, sep="\t")
    h = fam[fam["Species ID"].astype(str) == "9606"]
    return h.set_index("MiRBase ID")["miR family"]


def build_families(which=None, top=10) -> None:
    fam_dir = FIG / "families"; fam_dir.mkdir(parents=True, exist_ok=True)
    m = mapped_loci()
    a2f = _arm_family()
    m = m.copy(); m["family"] = m["arm"].map(lambda a: a2f.get(a, a2f.get(_strip(a))))
    m = m.dropna(subset=["family"])
    if which and which != ["top"]:
        fams = which
    else:
        fams = (m.drop_duplicates("arm").groupby("family")["he_degree"].sum()
                .sort_values(ascending=False).head(top).index.tolist())
    for f in fams:
        sub = m[m["family"] == f]
        if sub.empty:
            print(f"  [family] {f}: no mapped loci, skip"); continue
        expd = sub[sub["tier"] != "silent"]; sild = sub[sub["tier"] == "silent"]
        fig, ax = plt.subplots(figsize=(12, 8))
        _draw_chroms(ax, sub)
        if len(sild):                                          # silent arms = grey hollow rings (as in build_genome)
            ax.scatter(sild["start"], sild["chrom"].map(_YMAP), s=15, facecolors="none", edgecolors="0.6",
                       linewidths=0.6, zorder=2, label=f"silent ({sild['arm'].nunique()}, never ≥10 RPM)")
        sc = _scatter(ax, expd, expd["he_degree"].min(), expd["he_degree"].max()) if len(expd) else None
        _place_labels(ax, [(r.start, r.chrom, f"{r.arm} ({int(r.he_degree)})") for r in sub.itertuples()], fontsize=6.5)
        _finish(fig, ax, sc, f"Seed family '{f}' — {expd['arm'].nunique()} expressed + {sild['arm'].nunique()} silent "
                             f"(grey rings) / {len(sub)} loci")
        if len(sild):
            ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
        fig.savefig(fam_dir / f"family_{f.replace('/', '_')}.png", dpi=140); plt.close(fig)
        print(f"  [family] {f}: {expd['arm'].nunique()} expressed + {sild['arm'].nunique()} silent, {len(sub)} loci")


def _gene_loci() -> pd.DataFrame:
    """gene_name -> (chrom, start) from GENCODE v49 (gene rows on canonical chromosomes)."""
    from pipeline.config import PATHS
    g = pd.read_parquet(PATHS.gencode_gtf_pq, columns=["feature", "seqname", "start", "gene_name"])
    g = g[(g["feature"] == "gene") & (g["seqname"].isin(HG38))].drop_duplicates("gene_name").copy()
    g["start"] = g["start"].astype(int)
    return g.set_index("gene_name")[["seqname", "start"]].rename(columns={"seqname": "chrom"})


def build_family_bipartite(family: str, max_edges: int = 600) -> None:
    """PILOT: bipartite genome 'mirror' for one seed family — miRNA arm loci on the RIGHT, their HE-edge
    target-gene loci on the LEFT (mirrored genome), edges connecting each arm to its targets (coloured by arm)."""
    m = mapped_loci(); a2f = _arm_family()
    m = m.copy(); m["family"] = m["arm"].map(lambda a: a2f.get(a, a2f.get(_strip(a))))
    arms_m = m[m["family"] == family]
    if arms_m.empty:
        print(f"  [bipartite] {family}: no mapped arm loci"); return
    fam_arms = sorted(arms_m["arm"].unique())
    arm_pos = arms_m.sort_values("start").groupby("arm").first()[["chrom", "start"]]  # representative locus/arm

    hs = HallmarkSets.load()
    edges = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    gl = _gene_loci()
    fe = edges[edges["miRNA"].isin(fam_arms)][["miRNA", "gene"]].drop_duplicates()
    fe = fe.join(gl, on="gene").dropna(subset=["chrom"])
    fe = fe.merge(arm_pos.rename(columns={"chrom": "m_chrom", "start": "m_start"}), left_on="miRNA", right_index=True)
    gdeg = fe.groupby("gene")["miRNA"].nunique()

    X = max(HG38.values())
    fig, ax = plt.subplots(figsize=(15, 9))
    for c in _CHR_ORDER:                                   # mirrored chromosome lines + centromeres
        y, L, cen, gap = _YMAP[c], HG38[c], CEN[c], HG38[c] * 0.004
        for sgn in (1, -1):
            ax.plot([0, sgn * (cen - gap)], [y, y], color="0.88", lw=3, solid_capstyle="round", zorder=1)
            ax.plot([sgn * (cen + gap), sgn * L], [y, y], color="0.88", lw=3, solid_capstyle="round", zorder=1)
            ax.plot([sgn * cen], [y], marker="o", ms=3.5, color="0.6", zorder=2)
    ax.axvline(0, color="0.6", lw=0.8, zorder=1)

    cmap = plt.cm.tab10
    arm_color = {a: cmap(i % 10) for i, a in enumerate(fam_arms)}
    draw = fe.sample(max_edges, random_state=0) if len(fe) > max_edges else fe
    for r in draw.itertuples():                            # edges: arm (right) -> target gene (left)
        ax.plot([r.m_start, -r.start], [_YMAP[r.m_chrom], _YMAP[r.chrom]],
                color=arm_color[r.miRNA], lw=0.3, alpha=0.22, zorder=2)
    for a, row in arm_pos.iterrows():                      # miRNA arms (right), labelled
        ax.scatter([row["start"]], [_YMAP[row["chrom"]]], s=55, color=arm_color[a],
                   edgecolors="k", linewidths=0.4, zorder=5)
    _place_labels(ax, [(row["start"], row["chrom"], a) for a, row in arm_pos.iterrows()], fontsize=7)
    g_in = gl.loc[gl.index.isin(gdeg.index)].copy(); g_in["indeg"] = gdeg.reindex(g_in.index)
    ax.scatter(-g_in["start"], g_in["chrom"].map(_YMAP), s=8 + 10 * g_in["indeg"], c="0.35",
               alpha=0.5, edgecolors="none", zorder=3)
    top_g = g_in.sort_values("indeg", ascending=False).head(12)
    _place_labels(ax, [(-r.start, r.chrom, f"{g}") for g, r in top_g.iterrows()], fontsize=6)

    ax.set_yticks([_YMAP[c] for c in _CHR_ORDER]); ax.set_yticklabels([c[3:] for c in _CHR_ORDER], fontsize=8)
    ax.set_xlim(-X * 1.05, X * 1.05); ax.set_ylim(0.2, len(_CHR_ORDER) + 2.6)
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{abs(x) / 1e8:.1f}"))  # both sides = true positive locus
    ax.set_xlabel("◄ target-gene locus              genomic position (×1e8 bp, hg38)              miRNA-arm locus ►")
    ax.set_ylabel("chromosome")
    ax.set_title(f"Bipartite mirror — family '{family}': {len(fam_arms)} arms → "
                 f"{g_in.shape[0]} target genes, {len(fe)} HE edges ({len(draw)} drawn)")
    handles = [plt.Line2D([], [], color=arm_color[a], lw=2, label=a) for a in fam_arms]
    ax.legend(handles=handles, fontsize=6.5, loc="lower left", ncol=2, framealpha=0.9, title="miRNA arm")
    fig.tight_layout()
    out = FIG / "families" / f"bipartite_{family.replace('/', '_')}.png"
    fig.savefig(out, dpi=150); plt.close(fig)
    print(f"  [bipartite] {family}: {len(fam_arms)} arms, {g_in.shape[0]} genes, {len(fe)} edges -> {out}")


def _seed_map() -> dict:
    f = pd.read_csv(C.MIR_FAMILY_INFO, sep="\t")
    f = f[f["Species ID"].astype(str) == "9606"]
    return {r["MiRBase ID"]: r["Seed+m8"] for _, r in f.iterrows()}


@lru_cache(maxsize=1)
def _site_counts() -> dict:
    """(arm, gene) -> # predicted 3'UTR sites (TargetScan hsa rows per pair). For edge weighting."""
    from mirna_hallmark.analyses.misc.genome_wide_promiscuity import TARGETSCAN_CONTEXT
    parts = []
    for ch in pd.read_csv(TARGETSCAN_CONTEXT, sep="\t", usecols=["Gene Symbol", "miRNA"], chunksize=2_000_000):
        ch = ch[ch["miRNA"].astype(str).str.startswith("hsa-")]
        if not ch.empty:
            parts.append(ch)
    s = pd.concat(parts).groupby(["miRNA", "Gene Symbol"]).size()
    return {(m, g): int(n) for (m, g), n in s.items()}


def build_seed_hubs(top: int = 30) -> None:
    """Per-gene 3'UTR seed diversity (# DISTINCT seeds targeting the gene, genome-wide) → the convergence hubs."""
    seed = _seed_map()
    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    e = e.assign(seed=e["miRNA"].map(seed))
    out = pd.concat([e.dropna(subset=["seed"]).groupby("gene")["seed"].nunique().rename("n_distinct_seeds"),
                     e.groupby("gene")["miRNA"].nunique().rename("n_miRNAs")], axis=1)
    out = out.sort_values("n_distinct_seeds", ascending=False).reset_index()
    out.to_csv(OUT / "gene_seed_diversity.tsv", sep="\t", index=False)
    print("[seed-hubs] multi-seed 3'UTR hubs -> gene_seed_diversity.tsv; top: "
          + ", ".join(f"{r.gene}({int(r.n_distinct_seeds)})" for r in out.head(top).itertuples()))


def build_cluster_coordination(anchor: str = "hsa-miR-17-5p") -> None:
    """Bipartite COORDINATION view of the polycistronic cluster containing ``anchor``: member arms (left,
    coloured by seed) → genes co-targeted by ≥2 **distinct seeds** (right, ordered by # converging seeds).
    Shows genuine combinatorial coordination (different seeds converging), not seed-redundancy."""
    m = assign_clusters(mapped_loci())
    hit = m[(m["arm"] == anchor) | (m["arm"].map(_strip) == _strip(anchor))]
    if hit.empty:
        print(f"  [coord] anchor {anchor} not mapped"); return
    cl = hit["cluster"].iloc[0]
    g = m[m["cluster"] == cl]
    arms = sorted(g["arm"].unique())
    chrom, lo, hi = g["chrom"].iloc[0], int(g["start"].min()), int(g["start"].max())
    seed = _seed_map()
    hs = HallmarkSets.load()
    edges = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    tgt = edges.groupby("miRNA")["gene"].apply(set).to_dict()

    seedsets = {}
    for a in arms:
        seedsets.setdefault(seed.get(a, a), set()).update(tgt.get(a, set()))
    nseed = {}
    for s, ts in seedsets.items():
        for gn in ts:
            nseed[gn] = nseed.get(gn, 0) + 1
    shared = sorted([gn for gn, c in nseed.items() if c >= 2], key=lambda gn: (-nseed[gn], gn))
    if not shared:
        print(f"  [coord] cluster {cl}: no genes co-targeted by ≥2 seeds"); return

    seeds_sorted = sorted(set(seed.get(a, a) for a in arms))
    scol = {s: plt.cm.tab10(i % 10) for i, s in enumerate(seeds_sorted)}
    arms = sorted(arms, key=lambda a: (seed.get(a, a), a))
    ay = {a: i for i, a in enumerate(arms)}
    gy = {gn: i for i, gn in enumerate(shared)}
    ng, na = len(shared), len(arms)
    fig, ax = plt.subplots(figsize=(11, max(6, 0.34 * max(ng, na) + 1)))
    # scale arm/gene columns to a common height
    def _ya(a): return ay[a] / max(na - 1, 1) * (ng - 1) if ng > 1 else ay[a]
    sites = _site_counts()                                 # edge width ∝ # 3'UTR sites (TargetScan)
    for a in arms:
        for gn in tgt.get(a, set()) & set(shared):
            ns = sites.get((a, gn), sites.get((_strip(a), gn), 1))
            ax.plot([0, 1], [_ya(a), gy[gn]], color=scol[seed.get(a, a)],
                    lw=0.5 + 0.5 * min(ns, 8), alpha=min(0.4 + 0.08 * ns, 0.9), zorder=2)
    for a in arms:
        ax.scatter([0], [_ya(a)], s=60, color=scol[seed.get(a, a)], edgecolors="k", linewidths=0.4, zorder=4)
        ax.text(-0.02, _ya(a), a.replace("hsa-", ""), ha="right", va="center", fontsize=7.5)
    for gn in shared:
        ax.scatter([1], [gy[gn]], s=30 + 20 * nseed[gn], color="0.3", zorder=4)
        ax.text(1.02, gy[gn], f"{gn} ({nseed[gn]})", ha="left", va="center", fontsize=7.5,
                fontweight="bold" if nseed[gn] >= 3 else "normal")
    ax.set_xlim(-0.28, 1.28); ax.set_ylim(-2.5, ng + 1.4); ax.axis("off")
    ax.text(0, ng + 0.7, "miRNA arms", ha="center", fontsize=9, fontweight="bold")
    ax.text(1, ng + 0.7, "co-targeted genes  (label = # converging seeds)", ha="center", fontsize=9, fontweight="bold")
    ax.set_title(f"Coordination — polycistronic cluster {chrom}:{lo:,}–{hi:,}:  {na} arms / {len(seeds_sorted)} seeds "
                 f"→ {len(shared)} genes co-targeted by ≥2 DISTINCT seeds  (edge width ∝ # 3'UTR sites)", fontsize=10.5, pad=12)
    handles = [plt.Line2D([], [], color=scol[s], lw=2.5, label=s) for s in seeds_sorted]
    ax.legend(handles=handles, fontsize=7, loc="upper center", ncol=min(len(seeds_sorted), 8),
              framealpha=0.9, title="seed (Seed+m8)", bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout()
    out = FIG / "families" / f"coordination_{anchor.replace('hsa-', '').replace('/', '_')}.png"
    fig.savefig(out, dpi=150); plt.close(fig)
    print(f"  [coord] cluster {chrom}:{lo:,}-{hi:,}: {na} arms / {len(seeds_sorted)} seeds → {len(shared)} co-targeted genes -> {out}")


def build_gene_hub(gene: str = "PTEN") -> None:
    """Convergence hub for one gene: the distinct seeds (and their arms) targeting it, ordered by
    # 3'UTR sites; left-node size ∝ # arms in the seed, edge width/colour ∝ # TargetScan sites."""
    from collections import defaultdict
    seed = _seed_map(); sites = _site_counts()
    hs = HallmarkSets.load()
    e = load_mirtar_edges(sorted(hs.universe), resolve_arms=True)
    sub = e[e["gene"] == gene]
    if sub.empty:
        print(f"  [hub] {gene}: no edges"); return
    from mirna_hallmark.arm_expression import arm_expression_tiers       # show only expressed regulators
    _tier = arm_expression_tiers()["tier"]
    n_all = sub["miRNA"].nunique()
    arms = sorted(a for a in sub["miRNA"].unique()
                  if _tier.get(a, _tier.get(_strip(a), "silent")) != "silent")
    byseed = defaultdict(list)
    for a in arms:
        byseed[seed.get(a, a)].append(a)
    rows = []
    for s, alist in byseed.items():
        ns = sum(sites.get((a, gene), sites.get((_strip(a), gene), 1)) for a in alist)
        rep = max(alist, key=lambda a: sites.get((a, gene), 1))
        rows.append((s, len(alist), ns, rep))
    rows.sort(key=lambda r: (-r[2], -r[1]))
    n = len(rows); gy = (n - 1) / 2
    vmax = max(r[2] for r in rows)
    norm = plt.Normalize(1, vmax); cmap = plt.cm.viridis
    fig, ax = plt.subplots(figsize=(10, max(7, 0.23 * n)))
    for i, (s, na, ns, rep) in enumerate(rows):
        y = n - 1 - i; col = cmap(norm(ns))
        ax.plot([0, 1], [y, gy], color=col, lw=0.6 + 0.5 * min(ns, 8), alpha=0.75, zorder=2)
        ax.scatter([0], [y], s=18 + 13 * na, color=col, edgecolors="0.4", linewidths=0.3, zorder=3)
        ax.text(-0.015, y, rep.replace("hsa-", "") + (f"  +{na-1}" if na > 1 else ""),
                ha="right", va="center", fontsize=5.6)
    ax.scatter([1], [gy], s=650, color="#b3261e", edgecolors="k", linewidths=0.6, zorder=5)
    ax.text(1.035, gy, gene, ha="left", va="center", fontsize=14, fontweight="bold")
    ax.set_xlim(-0.36, 1.12); ax.set_ylim(-1.5, n + 0.5); ax.axis("off")   # tight right margin: PTEN at frame edge, no dead fan-space
    ax.set_title(f"{gene} — converged on by {n} distinct seeds / {len(arms)} EXPRESSED miRNAs (of {n_all}; silent arms excluded)\n"
                 f"left-node size ∝ # arms in seed · edge width/colour ∝ # TargetScan 3'UTR sites", fontsize=11)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=ax, label="# 3'UTR sites on " + gene + " (TargetScan)", shrink=0.5)
    fig.tight_layout()
    out = FIG / f"hub_{gene}.png"
    fig.savefig(out, dpi=150); plt.close(fig)
    print(f"  [hub] {gene}: {n} seeds / {len(arms)} expressed miRNAs (of {n_all}; {n_all - len(arms)} silent dropped) -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--families", help="comma-list of seed families, or 'top' for top-10 by HE-degree")
    ap.add_argument("--bipartite", help="seed family name for the bipartite gene-mirror pilot")
    ap.add_argument("--cluster", help="anchor arm (e.g. hsa-miR-17-5p) for the cluster-coordination figure")
    ap.add_argument("--seed-hubs", action="store_true", help="write per-gene 3'UTR seed-diversity table")
    ap.add_argument("--gene-hub", help="gene name for the convergence-hub figure (e.g. PTEN)")
    args = ap.parse_args()
    build_genome()
    if args.families:
        build_families(["top"] if args.families == "top" else args.families.split(","))
    if args.bipartite:
        build_family_bipartite(args.bipartite)
    if args.cluster:
        build_cluster_coordination(args.cluster)
    if args.seed_hubs:
        build_seed_hubs()
    if args.gene_hub:
        build_gene_hub(args.gene_hub)


if __name__ == "__main__":
    main()
