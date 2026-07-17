"""Multi-resolution genome dispersion figures for miRNA hairpin loci.

Builds summary tables and PNG panels at chromosome, cytoband-arm, fixed-bin,
and megacluster zoom resolutions. Outputs under ``mirna_locus_cnv/maps/figures/``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cnv_locus.mirna_cnv_genome_maps import (
    CHR19_MEGACLUSTER_ID,
    CHR_ORDER,
    CHR_SIZE_MB,
    PAM50_ORDER,
    _chrom_offsets,
    _cytoband_arm,
    build_chromosome_arm_cn_summary,
    build_chromosome_cn_summary,
    build_dense_locus_clusters,
)
CHR19_ZOOM_MB = (53.4, 54.0)
BIN_SIZE_MB = 5.0


def build_dispersion_summary_tables(
    locus_map: pd.DataFrame,
    paralog_map: pd.DataFrame,
    *,
    bin_size_mb: float = BIN_SIZE_MB,
) -> Dict[str, pd.DataFrame]:
    """Tables at chromosome, arm, and fixed-bin resolution."""
    chrom = build_chromosome_cn_summary(locus_map, paralog_map)
    arm = build_chromosome_arm_cn_summary(locus_map, paralog_map)

    loc = locus_map.copy()
    loc["cytoband_arm"] = [_cytoband_arm(r.chrom, r.midpoint) for r in loc.itertuples(index=False)]
    loc["mid_mb"] = loc["midpoint"] / 1e6
    loc["bin_start_mb"] = (np.floor(loc["mid_mb"] / bin_size_mb) * bin_size_mb).round(3)
    loc["genome_bin"] = loc["chrom"] + ":" + loc["bin_start_mb"].astype(str)

    bins = (
        loc.groupby(["chrom", "genome_bin", "bin_start_mb"], as_index=False)
        .agg(
            n_hairpin_loci=("locus_id", "count"),
            mean_locus_delta_vs_diploid=("delta_vs_diploid", "mean"),
        )
        .sort_values(["chrom", "bin_start_mb"])
    )

    clusters = build_dense_locus_clusters(locus_map)
    return {"chromosome": chrom, "cytoband_arm": arm, "genome_bins": bins, "dense_clusters": clusters}


def build_dispersion_summary_tables_for_pam50(
    locus_map: pd.DataFrame,
    paralog_map: pd.DataFrame,
    pam50: str,
    *,
    bin_size_mb: float = BIN_SIZE_MB,
) -> Dict[str, pd.DataFrame]:
    """Subtype-stratified tables using ``delta_vs_diploid_{pam50}`` from the locus map."""
    delta_col = f"delta_vs_diploid_{pam50}"
    if delta_col not in locus_map.columns:
        return {}

    loc = locus_map.copy()
    loc["stratum_delta"] = pd.to_numeric(loc[delta_col], errors="coerce")
    loc["cytoband_arm"] = [_cytoband_arm(r.chrom, r.midpoint) for r in loc.itertuples(index=False)]
    loc["mid_mb"] = loc["midpoint"] / 1e6
    loc["bin_start_mb"] = (np.floor(loc["mid_mb"] / bin_size_mb) * bin_size_mb).round(3)
    loc["genome_bin"] = loc["chrom"] + ":" + loc["bin_start_mb"].astype(str)

    chrom_rows: List[dict] = []
    mimat_by_chrom = paralog_map.groupby("chrom")["mimat"].nunique()
    for chrom in CHR_ORDER:
        sub = loc.loc[loc["chrom"] == chrom]
        if sub.empty:
            continue
        chrom_rows.append(
            {
                "PAM50_final": pam50,
                "chrom": chrom,
                "chr_size_mb": round(CHR_SIZE_MB.get(chrom, np.nan), 2),
                "n_hairpin_loci": int(len(sub)),
                "n_mature_arms": int(mimat_by_chrom.get(chrom, 0)),
                "mean_locus_delta_vs_diploid": round(float(sub["stratum_delta"].mean()), 4),
                "median_locus_delta_vs_diploid": round(float(sub["stratum_delta"].median()), 4),
            }
        )
    chrom = pd.DataFrame(chrom_rows)

    arm_rows: List[dict] = []
    for arm_name, sub in loc.groupby("cytoband_arm"):
        if not arm_name:
            continue
        chrom_name = arm_name[:-1] if arm_name.endswith(("p", "q")) else arm_name
        arm_rows.append(
            {
                "PAM50_final": pam50,
                "cytoband_arm": arm_name,
                "chrom": chrom_name,
                "n_hairpin_loci": int(len(sub)),
                "mean_locus_delta_vs_diploid": round(float(sub["stratum_delta"].mean()), 4),
            }
        )
    arm = pd.DataFrame(arm_rows).sort_values("n_hairpin_loci", ascending=False)

    bins = (
        loc.groupby(["chrom", "genome_bin", "bin_start_mb"], as_index=False)
        .agg(
            PAM50_final=(delta_col, lambda _: pam50),
            n_hairpin_loci=("locus_id", "count"),
            mean_locus_delta_vs_diploid=("stratum_delta", "mean"),
        )
        .sort_values(["chrom", "bin_start_mb"])
    )
    bins["mean_locus_delta_vs_diploid"] = bins["mean_locus_delta_vs_diploid"].round(4)

    clusters = build_dense_locus_clusters(locus_map)
    if not clusters.empty:
        cl = clusters.copy()
        cl["PAM50_final"] = pam50
        cl["mean_locus_delta_vs_diploid"] = cl.apply(
            lambda r: round(
                float(
                    loc.loc[
                        (loc["chrom"] == r["chrom"])
                        & (loc["mid_mb"] >= r["cluster_start_mb"] - 0.001)
                        & (loc["mid_mb"] <= r["cluster_end_mb"] + 0.001)
                    ]["stratum_delta"].mean()
                ),
                4,
            ),
            axis=1,
        )
    else:
        cl = clusters

    return {"chromosome": chrom, "cytoband_arm": arm, "genome_bins": bins, "dense_clusters": cl}


def _attach_genome_x(locus_map: pd.DataFrame) -> pd.DataFrame:
    catalog = locus_map.rename(columns={"locus_id": "entity_id"})
    offsets = _chrom_offsets(catalog)
    out = locus_map.copy()
    out["genome_x"] = out.apply(lambda r: offsets.get(r.chrom, 0) + r.midpoint, axis=1)
    out["mid_mb"] = out["midpoint"] / 1e6
    out["cytoband_arm"] = [_cytoband_arm(r.chrom, r.midpoint) for r in out.itertuples(index=False)]
    return out


def plot_genome_dispersion_multires(
    locus_map: pd.DataFrame,
    summaries: Dict[str, pd.DataFrame],
    *,
    out_path: Path,
    chr19_cluster_id: str = CHR19_MEGACLUSTER_ID,
    delta_col: str = "delta_vs_diploid",
    title_suffix: str = "",
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    loc = _attach_genome_x(locus_map)
    if delta_col not in loc.columns:
        return
    loc["plot_delta"] = pd.to_numeric(loc[delta_col], errors="coerce")
    chrom_df = summaries["chromosome"]
    arm_df = summaries["cytoband_arm"]
    bin_df = summaries["genome_bins"]

    fig = plt.figure(figsize=(16, 14), dpi=150)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 1.0, 1.0], hspace=0.35, wspace=0.28)

    chrom_palette = plt.cm.tab20.colors
    chroms = [c for c in CHR_ORDER if c in set(loc["chrom"])]
    color_map = {c: chrom_palette[i % len(chrom_palette)] for i, c in enumerate(chroms)}

    # A — whole genome Manhattan (Δ CN)
    ax_a = fig.add_subplot(gs[0, :])
    for chrom in chroms:
        sub = loc.loc[loc["chrom"] == chrom]
        ax_a.scatter(
            sub["genome_x"] / 1e6,
            sub["plot_delta"],
            c=[color_map[chrom]],
            s=16,
            alpha=0.75,
            linewidths=0,
        )
    ax_a.axhline(0, color="0.35", lw=0.8, ls="--")
    ax_a.axhline(1, color="0.7", lw=0.6, ls=":")
    ax_a.set_ylabel("Locus Δ CN vs diploid")
    ax_a.set_xlabel("Pseudo-genome position (Gb)")
    ax_a.set_title(f"A — 506 hairpin loci across hg38 (colored by chromosome){title_suffix}")
    ax_a.set_xticks([])

    # B — loci per chromosome
    ax_b = fig.add_subplot(gs[1, 0])
    plot_ch = chrom_df.sort_values("chrom", key=lambda s: s.map({c: i for i, c in enumerate(CHR_ORDER)}))
    colors = [color_map.get(c, "0.5") for c in plot_ch["chrom"]]
    ax_b.bar(range(len(plot_ch)), plot_ch["n_hairpin_loci"], color=colors, edgecolor="none")
    ax_b.set_xticks(range(len(plot_ch)))
    ax_b.set_xticklabels([c.replace("chr", "") for c in plot_ch["chrom"]], rotation=90, fontsize=7)
    ax_b.set_ylabel("# hairpin loci")
    ax_b.set_title("B — loci per chromosome")

    # C — loci per cytoband arm (top 24 by count)
    ax_c = fig.add_subplot(gs[1, 1])
    top_arm = arm_df.sort_values("n_hairpin_loci", ascending=False).head(24)
    ax_c.barh(top_arm["cytoband_arm"][::-1], top_arm["n_hairpin_loci"][::-1], color="steelblue", edgecolor="none")
    ax_c.set_xlabel("# hairpin loci")
    ax_c.set_title("C — top cytoband arms (p/q)")

    # D — 5 Mb bin density (loci per bin)
    ax_d = fig.add_subplot(gs[2, 0])
    bin_plot = bin_df.copy()
    bin_plot["bin_idx"] = np.arange(len(bin_plot))
    ax_d.bar(
        bin_plot["bin_idx"],
        bin_plot["n_hairpin_loci"],
        width=1.0,
        color="darkorange",
        edgecolor="none",
        alpha=0.85,
    )
    ax_d.set_ylabel("# loci / 5 Mb bin")
    ax_d.set_xlabel("Genome bins (chr order)")
    ax_d.set_title(f"D — {BIN_SIZE_MB:g} Mb bin occupancy")
    ax_d.set_xticks([])

    # E — chr19 megacluster zoom
    ax_e = fig.add_subplot(gs[2, 1])
    z0, z1 = CHR19_ZOOM_MB
    chr19 = loc.loc[loc["chrom"] == "chr19"].copy()
    zoom = chr19.loc[(chr19["mid_mb"] >= z0) & (chr19["mid_mb"] <= z1)]
    cluster_loci = set()
    cl = summaries["dense_clusters"]
    if not cl.empty:
        row = cl.loc[cl.apply(lambda r: f"{r.chrom}:{r.cluster_start_mb:.3f}-{r.cluster_end_mb:.3f}", axis=1) == chr19_cluster_id]
        if not row.empty:
            r0 = row.iloc[0]
            cluster_loci = set(
                chr19.loc[
                    (chr19["mid_mb"] >= r0["cluster_start_mb"] - 0.001)
                    & (chr19["mid_mb"] <= r0["cluster_end_mb"] + 0.001)
                ]["locus_id"].astype(str)
            )
    for _, r in chr19.iterrows():
        in_cluster = str(r.locus_id) in cluster_loci
        pd_val = loc.loc[loc["locus_id"] == r.locus_id, "plot_delta"]
        y = float(pd_val.iloc[0]) if len(pd_val) else np.nan
        ax_e.scatter(
            r.mid_mb,
            y,
            s=80 if in_cluster else 24,
            c="crimson" if in_cluster else "0.65",
            alpha=0.9 if in_cluster else 0.5,
            edgecolors="k" if in_cluster else "none",
            linewidths=0.4,
        )
    ax_e.axhline(0, color="0.35", lw=0.8, ls="--")
    ax_e.set_xlim(z0, z1)
    ax_e.set_xlabel("chr19 position (Mb)")
    ax_e.set_ylabel("Δ CN vs diploid")
    ax_e.set_title(f"E — chr19 zoom ({z0}–{z1} Mb); megacluster highlighted")

    fig.suptitle(
        f"miRNA hairpin locus genome dispersion (multi-resolution){title_suffix}",
        fontsize=13,
        y=0.995,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_chromosome_facets(
    locus_map: pd.DataFrame,
    out_dir: Path,
    *,
    delta_col: str = "delta_vs_diploid",
    filename: str = "mirna_locus_dispersion_by_chromosome_facets.png",
    title_suffix: str = "",
) -> Path:
    """One panel per chromosome: position vs Δ CN."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    loc = _attach_genome_x(locus_map)
    if delta_col not in loc.columns:
        return out_dir / filename
    loc["plot_delta"] = pd.to_numeric(loc[delta_col], errors="coerce")
    chroms = [c for c in CHR_ORDER if c in set(loc["chrom"])]
    n = len(chroms)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 2.2 * nrows), dpi=130)
    axes = np.array(axes).reshape(-1)

    for i, chrom in enumerate(chroms):
        ax = axes[i]
        sub = loc.loc[loc["chrom"] == chrom]
        ax.scatter(sub["mid_mb"], sub["plot_delta"], s=12, alpha=0.8, c="tab:blue", linewidths=0)
        ax.axhline(0, color="0.4", lw=0.6, ls="--")
        ax.set_title(chrom, fontsize=9)
        ax.set_xlabel("Mb", fontsize=7)
        ax.set_ylabel("Δ CN", fontsize=7)
        ax.tick_params(labelsize=6)
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"miRNA loci by chromosome (within-chr Mb){title_suffix}", fontsize=12)
    fig.tight_layout()
    out = out_dir / filename
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return out


def write_dispersion_outputs(
    *,
    locus_map_path: Optional[Path] = None,
    paralog_map_path: Optional[Path] = None,
    out_dir: Optional[Path] = None,
    include_by_pam50: bool = True,
) -> Dict[str, Path]:
    locus_map_path = Path(locus_map_path or C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_locus_genome_map.tsv")
    paralog_map_path = Path(paralog_map_path or C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_mimat_paralog_map.tsv")
    out_dir = Path(out_dir or C.MIRNA_LOCUS_CNV_DIR / "maps" / "figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    locus_map = pd.read_csv(locus_map_path, sep="\t")
    paralog_map = pd.read_csv(paralog_map_path, sep="\t")
    summaries = build_dispersion_summary_tables(locus_map, paralog_map)

    paths: Dict[str, Path] = {}
    for name, df in summaries.items():
        p = out_dir.parent / f"mirna_locus_dispersion_{name}.tsv"
        df.to_csv(p, sep="\t", index=False)
        paths[f"table_{name}"] = p

    paths["multires_fig"] = out_dir / "mirna_locus_genome_dispersion_multires.png"
    plot_genome_dispersion_multires(locus_map, summaries, out_path=paths["multires_fig"])
    paths["chrom_facets_fig"] = plot_chromosome_facets(locus_map, out_dir)

    paths["chrom_facets_fig"] = plot_chromosome_facets(locus_map, out_dir)

    pam50_figures: Dict[str, Dict[str, str]] = {}
    if include_by_pam50:
        pam_dir = out_dir / "by_pam50"
        table_root = out_dir.parent / "by_pam50"
        for pam in PAM50_ORDER:
            delta_col = f"delta_vs_diploid_{pam}"
            if delta_col not in locus_map.columns:
                continue
            pam_summaries = build_dispersion_summary_tables_for_pam50(locus_map, paralog_map, pam)
            if not pam_summaries:
                continue
            pam_table_dir = table_root / pam
            pam_table_dir.mkdir(parents=True, exist_ok=True)
            pam_fig_dir = pam_dir / pam
            pam_fig_dir.mkdir(parents=True, exist_ok=True)
            for name, df in pam_summaries.items():
                p = pam_table_dir / f"mirna_locus_dispersion_{name}.tsv"
                df.to_csv(p, sep="\t", index=False)
                paths[f"table_{pam}_{name}"] = p
            suffix = f" — PAM50 {pam} (mean locus Δ within subtype)"
            multires = pam_fig_dir / "mirna_locus_genome_dispersion_multires.png"
            plot_genome_dispersion_multires(
                locus_map,
                pam_summaries,
                out_path=multires,
                delta_col=delta_col,
                title_suffix=suffix,
            )
            facets = plot_chromosome_facets(
                locus_map,
                pam_fig_dir,
                delta_col=delta_col,
                title_suffix=suffix,
            )
            paths[f"multires_fig_{pam}"] = multires
            paths[f"chrom_facets_fig_{pam}"] = facets
            pam50_figures[pam] = {"multires": str(multires), "facets": str(facets)}

    manifest = {
        "n_loci": int(len(locus_map)),
        "n_chromosomes": int(locus_map["chrom"].nunique()),
        "n_cytoband_arms": int(summaries["cytoband_arm"]["cytoband_arm"].nunique()),
        "n_genome_bins_5mb": int(len(summaries["genome_bins"])),
        "n_dense_clusters": int(len(summaries["dense_clusters"])),
        "chr19_megacluster_id": CHR19_MEGACLUSTER_ID,
        "figures": {k: str(v) for k, v in paths.items() if str(v).endswith(".png") and "multires_fig_" not in k and "chrom_facets_fig_" not in k},
        "figures_by_pam50": pam50_figures,
    }
    manifest_path = out_dir / "dispersion_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    paths["manifest"] = manifest_path

    print(
        f"[mirna_locus_genome_dispersion] loci={len(locus_map)}; "
        f"figures -> {out_dir}"
    )
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--locus-map", type=Path, default=None)
    ap.add_argument("--paralog-map", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument(
        "--no-by-pam50",
        action="store_true",
        help="Skip per-PAM50 dispersion tables and figures",
    )
    args = ap.parse_args()
    write_dispersion_outputs(
        locus_map_path=args.locus_map,
        paralog_map_path=args.paralog_map,
        out_dir=args.out_dir,
        include_by_pam50=not args.no_by_pam50,
    )


if __name__ == "__main__":
    main()
