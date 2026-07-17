"""chr19 megacluster participant gallery (no BAM required).

One schematic per primary participant over the polycistronic block
(``chr19:53.667-53.762``, 46 hairpins). Hairpins are colored by dosage layer;
DEL/DUP SVs and the ASCAT3 segment are drawn on shared tracks.

Outputs under ``sv_overlap/chr19_megacluster/gallery/``.
"""

from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cnv_locus.mirna_cnv_genome_maps import (
    CHR19_MEGACLUSTER_ID,
    _loci_in_cluster_block,
    build_dense_locus_clusters,
)
from mirna_hallmark.analyses.cnv_locus.mirna_locus_sv_case_plots import (
    _clip_interval,
    _parse_interval,
    _safe_slug,
    _ucsc_url,
    draw_segment_track,
    draw_sv_event_stack,
    format_span_bp,
    schematic_y_layout,
    summarize_interval,
    summarize_sv_events,
    sv_events_summary_text,
)

LAYER_COLORS = {
    "neither": "#bdc3c7",
    "cnv_only": "#e67e22",
    "both": "#8e44ad",
    "sv_only": "#16a085",
}
PAD_BP = 80_000


def _cluster_window(cluster_row: pd.Series, *, pad_bp: int = PAD_BP) -> Tuple[str, int, int]:
    chrom = str(cluster_row["chrom"])
    s = int(float(cluster_row["cluster_start_mb"]) * 1e6) - pad_bp
    e = int(float(cluster_row["cluster_end_mb"]) * 1e6) + pad_bp
    return chrom, max(0, s), e


def _load_cluster_context() -> Tuple[pd.DataFrame, pd.Series, List[str], pd.DataFrame]:
    locus_map = pd.read_csv(C.MIRNA_LOCUS_CNV_DIR / "maps" / "mirna_cnv_locus_genome_map.tsv", sep="\t")
    clusters = build_dense_locus_clusters(locus_map)
    cl = clusters.loc[
        clusters.apply(
            lambda r: f"{r['chrom']}:{r['cluster_start_mb']:.3f}-{r['cluster_end_mb']:.3f}",
            axis=1,
        )
        == CHR19_MEGACLUSTER_ID
    ].iloc[0]
    block = _loci_in_cluster_block(locus_map, cl)
    loci = block["locus_id"].astype(str).tolist()
    coords = block[["locus_id", "locus_name", "chrom", "start", "end", "midpoint"]].copy()
    return locus_map, cl, loci, coords


def plot_chr19_participant_schematic(
    participant: str,
    ploci: pd.DataFrame,
    coords: pd.DataFrame,
    sv_rows: pd.DataFrame,
    *,
    window: Tuple[str, int, int],
    out_path: Path,
) -> Dict[str, object]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    wchrom, w0, w1 = window
    span_kb = (w1 - w0) / 1e3
    events = summarize_sv_events(sv_rows)
    seg_y, sv_y_top, hairpin_y, y_max = schematic_y_layout(len(events))

    fig_h = 4.0 + max(0, len(events) - 1) * 0.35
    fig, ax = plt.subplots(figsize=(15, fig_h), dpi=140)
    ax.set_xlim(w0, w1)
    ax.set_ylim(-0.45, y_max)
    ax.set_xlabel(f"{wchrom} (bp)")
    ax.set_yticks([hairpin_y, (sv_y_top + 0.2) / 2 if events else 0.8, seg_y])
    ax.set_yticklabels(["46 hairpins (by layer)", f"SV events (n={len(events)})", "ASCAT3 segment"])
    ax.grid(axis="x", alpha=0.25, lw=0.5)

    seg_label = ploci["overlap_segment"].dropna().astype(str)
    seg_raw = seg_label.iloc[0] if len(seg_label) else ""
    cn = float(ploci["copy_number"].median())
    cstate = str(ploci["cn_state"].mode().iloc[0]) if not ploci["cn_state"].mode().empty else ""
    cpe_val = pd.to_numeric(ploci["CPE"], errors="coerce").dropna()
    cpe = float(cpe_val.iloc[0]) if len(cpe_val) else None
    seg_info = draw_segment_track(ax, seg_raw, window, seg_y, copy_number=cn, cn_state=cstate, cpe=cpe)

    if events:
        draw_sv_event_stack(ax, events, window, y_top=sv_y_top)
    else:
        ax.text(w0, sv_y_top, "  (no overlapping SV)", fontsize=7, va="center")

    merged = coords.merge(
        ploci[["locus_id", "dosage_layer", "copy_number"]],
        on="locus_id",
        how="left",
    )
    for row in merged.itertuples(index=False):
        layer = str(getattr(row, "dosage_layer", "neither") or "neither")
        color = LAYER_COLORS.get(layer, "#7f8c8d")
        mid = int(row.midpoint)
        ax.vlines(mid, hairpin_y - 0.12, hairpin_y + 0.12, colors=color, linewidth=1.8, zorder=4)

    pam = ploci["PAM50_final"].dropna().astype(str).iloc[0] if ploci["PAM50_final"].notna().any() else ""
    layer_counts = ploci["dosage_layer"].value_counts().to_dict()
    lc_str = ", ".join(f"{k}={v}" for k, v in sorted(layer_counts.items()))
    sv_types = {str(e["svtype"]) for e in events}
    med_vaf = pd.to_numeric(sv_rows.get("tumor_vaf"), errors="coerce").median() if not sv_rows.empty else np.nan
    med_ratio = pd.to_numeric(sv_rows.get("tumor_vaf_ratio_adj"), errors="coerce").median() if not sv_rows.empty else np.nan
    title = (
        f"chr19 megacluster | {participant} | {pam} | mean CN={cn:.2f} | {lc_str}"
        + (f" | SV: {','.join(sorted(sv_types))}" if sv_types else "")
        + (f" | med VAF={med_vaf:.3f}" if np.isfinite(med_vaf) else "")
    )
    ax.set_title(title, fontsize=9, loc="left")

    legend_elems = [Line2D([0], [0], color=c, lw=3, label=k) for k, c in LAYER_COLORS.items()]
    legend_elems += [
        Line2D([0], [0], color="#e74c3c", lw=3, label="DEL"),
        Line2D([0], [0], color="#2980b9", lw=3, label="DUP"),
    ]
    ax.legend(handles=legend_elems, loc="upper right", fontsize=6, framealpha=0.9, ncol=2)

    igv_window = f"{wchrom}:{w0}-{w1}"
    sv_txt = sv_events_summary_text(events)
    footer = f"Window {span_kb:.0f} kb · {CHR19_MEGACLUSTER_ID} · {igv_window}"
    if sv_txt:
        footer += f"\nSVs: {sv_txt}"
    fig.text(0.01, 0.02, footer, fontsize=6.5, color="0.35")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)

    seg_span = seg_info or summarize_interval(seg_raw) or {}
    n_del = sum(1 for e in events if e["svtype"] == "DEL")
    n_dup = sum(1 for e in events if e["svtype"] == "DUP")
    return {
        "segment_span_label": seg_span.get("span_label", ""),
        "segment_span_bp": seg_span.get("span_bp", np.nan),
        "overlap_segment": seg_raw,
        "sv_events_summary": sv_txt,
        "n_sv_del": n_del,
        "n_sv_dup": n_dup,
        "has_mixed_del_dup": bool(n_del and n_dup),
        "median_tumor_vaf": round(float(med_vaf), 4) if np.isfinite(med_vaf) else np.nan,
        "median_vaf_ratio_adj": round(float(med_ratio), 4) if np.isfinite(med_ratio) else np.nan,
        "CPE": cpe,
    }


def write_chr19_megacluster_gallery(
    *,
    sv_dir: Optional[Path] = None,
    layers_path: Optional[Path] = None,
    hits_path: Optional[Path] = None,
    out_dir: Optional[Path] = None,
    pad_bp: int = PAD_BP,
) -> Dict[str, Path]:
    sv_dir = Path(sv_dir or C.MIRNA_LOCUS_SV_DIR)
    layers_path = Path(layers_path or sv_dir / "mirna_locus_cnv_sv_layers.tsv.gz")
    hits_path = Path(hits_path or sv_dir / "mirna_locus_sv_hits.tsv.gz")
    out_dir = Path(out_dir or sv_dir / "chr19_megacluster" / "gallery")
    out_dir.mkdir(parents=True, exist_ok=True)

    _, cluster_row, loci, coords = _load_cluster_context()
    window = _cluster_window(cluster_row, pad_bp=pad_bp)
    wchrom, w0, w1 = window
    igv_window = f"{wchrom}:{w0}-{w1}"

    layers = pd.read_csv(layers_path, sep="\t", low_memory=False)
    layers["locus_id"] = layers["entity_id"].astype(str)
    cl_layers = layers.loc[layers["locus_id"].isin(loci)].copy()

    hits = pd.read_csv(hits_path, sep="\t", low_memory=False)
    cl_hits = hits.loc[hits["locus_id"].astype(str).isin(loci)].copy()

    png_dir = out_dir / "participants"
    png_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: List[dict] = []

    participants = sorted(cl_layers["participant"].astype(str).unique())
    for participant in participants:
        ploci = cl_layers.loc[cl_layers["participant"] == participant]
        sv_sub = cl_hits.loc[cl_hits["participant"] == participant]
        slug = _safe_slug(participant)
        png = png_dir / f"{slug}.png"
        meta = plot_chr19_participant_schematic(
            participant, ploci, coords, sv_sub, window=window, out_path=png
        )
        layer_counts = ploci["dosage_layer"].value_counts().to_dict()
        dominant = ploci["dosage_layer"].mode().iloc[0] if not ploci["dosage_layer"].mode().empty else ""
        manifest_rows.append(
            {
                "participant": participant,
                "PAM50_final": ploci["PAM50_final"].dropna().astype(str).iloc[0]
                if ploci["PAM50_final"].notna().any()
                else "",
                "sample_vial": ploci["sample_vial"].dropna().astype(str).iloc[0]
                if "sample_vial" in ploci.columns and ploci["sample_vial"].notna().any()
                else "",
                "mean_copy_number": round(float(ploci["copy_number"].mean()), 3),
                "dominant_dosage_layer": dominant,
                "n_cnv_only": int(layer_counts.get("cnv_only", 0)),
                "n_both": int(layer_counts.get("both", 0)),
                "n_sv_only": int(layer_counts.get("sv_only", 0)),
                "n_neither": int(layer_counts.get("neither", 0)),
                "n_sv_events": int(sv_sub["sv_id"].nunique()) if not sv_sub.empty else 0,
                "n_sv_del": int(meta.get("n_sv_del", 0)),
                "n_sv_dup": int(meta.get("n_sv_dup", 0)),
                "has_mixed_del_dup": bool(meta.get("has_mixed_del_dup", False)),
                "segment_span_label": meta.get("segment_span_label", ""),
                "segment_span_bp": meta.get("segment_span_bp", np.nan),
                "overlap_segment": meta.get("overlap_segment", ""),
                "sv_events_summary": meta.get("sv_events_summary", ""),
                "median_tumor_vaf": meta.get("median_tumor_vaf", np.nan),
                "median_vaf_ratio_adj": meta.get("median_vaf_ratio_adj", np.nan),
                "CPE": meta.get("CPE", np.nan),
                "has_sv": bool(not sv_sub.empty),
                "igv_window": igv_window,
                "ucsc_url": _ucsc_url(igv_window),
                "plot_png": str(png.relative_to(out_dir)),
            }
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_dir / "chr19_gallery_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)

    index_path = out_dir / "index.html"
    index_path.write_text(_render_chr19_html(manifest, out_dir, cluster_row), encoding="utf-8")

    meta = {
        "cluster_id": CHR19_MEGACLUSTER_ID,
        "n_participants": int(len(manifest)),
        "n_hairpin_loci": len(loci),
        "window": igv_window,
        "pad_bp": pad_bp,
        "layer_summary": cl_layers["dosage_layer"].value_counts().to_dict(),
    }
    meta_path = out_dir / "chr19_gallery.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[chr19_megacluster_gallery] {len(manifest)} participants -> {out_dir}")
    return {"manifest": manifest_path, "index": index_path, "meta": meta_path}


def _render_chr19_html(manifest: pd.DataFrame, out_dir: Path, cluster_row: pd.Series) -> str:
    if manifest.empty:
        return "<html><body><p>No chr19 gallery plots.</p></body></html>"

    def _cards(df: pd.DataFrame) -> str:
        parts = []
        for row in df.itertuples(index=False):
            parts.append(
                f"""
                <div class="card" data-layer="{html.escape(str(row.dominant_dosage_layer))}"
                     data-pam50="{html.escape(str(row.PAM50_final))}" data-sv="{1 if row.has_sv else 0}"
                     data-mixed="{1 if getattr(row, 'has_mixed_del_dup', False) else 0}">
                  <img src="{html.escape(row.plot_png)}" alt="{html.escape(str(row.participant))}" loading="lazy"/>
                  <div class="meta">
                    <strong>{html.escape(str(row.participant))}</strong>
                    · {html.escape(str(row.PAM50_final))}
                    · mean CN {row.mean_copy_number}<br/>
                    <span class="seg">ASCAT3 segment: {html.escape(str(getattr(row, 'segment_span_label', '') or '—'))}</span>
                    · CPE {getattr(row, 'CPE', '')}<br/>
                    cnv_only={row.n_cnv_only} both={row.n_both}
                    sv_only={row.n_sv_only} neither={row.n_neither}
                    · DEL={getattr(row, 'n_sv_del', 0)} DUP={getattr(row, 'n_sv_dup', 0)}
                    {' · <b>mixed DEL+DUP</b>' if getattr(row, 'has_mixed_del_dup', False) else ''}
                    · med VAF={getattr(row, 'median_tumor_vaf', '')}<br/>
                    <span class="sv">{html.escape(str(getattr(row, 'sv_events_summary', '') or 'no SV'))}</span><br/>
                    <a href="{html.escape(row.ucsc_url)}" target="_blank">UCSC</a>
                  </div>
                </div>
                """
            )
        return "".join(parts)

    priority = manifest.loc[manifest["has_sv"]].sort_values(["n_sv_only", "n_both"], ascending=False)
    rest = manifest.loc[~manifest["has_sv"]].sort_values(["n_cnv_only", "mean_copy_number"], ascending=False)

    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"/>
<title>chr19 miRNA megacluster — all participants</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 1.5rem; background: #fafafa; }}
h1 {{ font-size: 1.25rem; }}
.toolbar {{ margin: 1rem 0; display: flex; gap: 0.75rem; flex-wrap: wrap; align-items: center; }}
select, input {{ padding: 0.35rem 0.5rem; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(440px, 1fr)); gap: 1rem; }}
.card {{ background: #fff; border: 1px solid #ddd; border-radius: 6px; overflow: hidden; }}
.card.hidden {{ display: none; }}
.card img {{ width: 100%; display: block; }}
.meta {{ padding: 0.6rem 0.75rem; font-size: 0.82rem; line-height: 1.35; }}
.meta .seg {{ color: #555; }}
.meta .sv {{ font-size: 0.78rem; color: #333; }}
.stats {{ font-size: 0.9rem; color: #444; margin-bottom: 0.5rem; }}
</style>
<script>
function applyFilters() {{
  const layer = document.getElementById('layerFilter').value;
  const pam = document.getElementById('pamFilter').value;
  const sv = document.getElementById('svFilter').value;
  const mixed = document.getElementById('mixedFilter').value;
  const q = document.getElementById('searchBox').value.toLowerCase();
  document.querySelectorAll('.card').forEach(el => {{
    const okLayer = !layer || el.dataset.layer === layer;
    const okPam = !pam || el.dataset.pam50 === pam;
    const okSv = !sv || el.dataset.sv === sv;
    const okMixed = !mixed || el.dataset.mixed === mixed;
    const okQ = !q || el.innerText.toLowerCase().includes(q);
    el.classList.toggle('hidden', !(okLayer && okPam && okSv && okMixed && okQ));
  }});
}}
</script>
</head>
<body>
<h1>chr19 megacluster ({html.escape(CHR19_MEGACLUSTER_ID)})</h1>
<p class="stats">{int(cluster_row['n_hairpin_loci'])} hairpins in {float(cluster_row['cluster_span_kb']):.0f} kb
· {len(manifest)} primary participants · window {html.escape(manifest.iloc[0]['igv_window'])}</p>
<p>ASCAT3 segment size labeled above gray bar; each SV on its own row (red=DEL, blue=DUP) with span.
Serve via <code>python3 -m http.server</code> in this folder for images to load in Cursor Simple Browser.</p>
<div class="toolbar">
  <label>Layer <select id="layerFilter" onchange="applyFilters()">
    <option value="">All</option>
    <option value="cnv_only">cnv_only</option>
    <option value="both">both</option>
    <option value="sv_only">sv_only</option>
    <option value="neither">neither</option>
  </select></label>
  <label>PAM50 <select id="pamFilter" onchange="applyFilters()">
    <option value="">All</option>
    {' '.join(f'<option value="{html.escape(p)}">{html.escape(p)}</option>' for p in sorted(manifest['PAM50_final'].dropna().unique()))}
  </select></label>
  <label>SV <select id="svFilter" onchange="applyFilters()">
    <option value="">All</option>
    <option value="1">Has SV</option>
    <option value="0">No SV</option>
  </select></label>
  <label>Mixed SV <select id="mixedFilter" onchange="applyFilters()">
    <option value="">All</option>
    <option value="1">DEL+DUP both</option>
  </select></label>
  <input id="searchBox" placeholder="Search participant…" oninput="applyFilters()"/>
</div>

<h2>Participants with SV overlap ({len(priority)})</h2>
<div class="grid">{_cards(priority)}</div>

<h2>All other participants ({len(rest)})</h2>
<div class="grid">{_cards(rest)}</div>
</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sv-dir", type=Path, default=None)
    ap.add_argument("--pad-bp", type=int, default=PAD_BP, help="Padding around cluster (default 80 kb)")
    args = ap.parse_args()
    write_chr19_megacluster_gallery(sv_dir=args.sv_dir, pad_bp=args.pad_bp)


if __name__ == "__main__":
    main()
