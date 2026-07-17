"""Static genomic schematics for miRNA SV/CNV review cases (no BAM/CRAM required).

Reads review-queue tables and overlapping SV hits; draws ASCAT3 segment, hairpin
locus, and DEL/DUP intervals on a shared x-axis. Writes PNGs plus a browsable HTML
index under ``sv_overlap/review_queues/case_plots/``.
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

_INTERVAL_RE = re.compile(r"^(chr\w+):(\d+)-(\d+)$")
_SVTYPE_COLORS = {"DEL": "#e74c3c", "DUP": "#2980b9"}


def format_span_bp(span_bp: object) -> str:
    bp = abs(int(pd.to_numeric(span_bp, errors="coerce") or 0))
    if bp >= 1_000_000:
        return f"{bp / 1e6:.2f} Mb"
    if bp >= 1_000:
        return f"{bp / 1e3:.1f} kb"
    return f"{bp} bp"


def summarize_interval(label: object) -> Optional[Dict[str, object]]:
    parsed = _parse_interval(label)
    if not parsed:
        return None
    chrom, start, end = parsed
    s, e = min(start, end), max(start, end)
    span_bp = e - s
    return {
        "chrom": chrom,
        "start": s,
        "end": e,
        "span_bp": span_bp,
        "span_label": format_span_bp(span_bp),
        "raw_label": str(label).strip(),
    }


def summarize_sv_events(sv_rows: pd.DataFrame) -> List[Dict[str, object]]:
    if sv_rows.empty:
        return []
    events: List[Dict[str, object]] = []
    for row in sv_rows.drop_duplicates("sv_id").itertuples(index=False):
        start = int(pd.to_numeric(getattr(row, "sv_start", 0), errors="coerce") or 0)
        end_raw = getattr(row, "sv_end", np.nan)
        end = int(pd.to_numeric(end_raw, errors="coerce") or start)
        span_bp = pd.to_numeric(getattr(row, "sv_span_bp", np.nan), errors="coerce")
        if not np.isfinite(span_bp) or span_bp <= 0:
            span_bp = abs(end - start)
        svtype = str(getattr(row, "svtype", ""))
        events.append(
            {
                "sv_id": str(getattr(row, "sv_id", "")),
                "svtype": svtype,
                "chrom": str(getattr(row, "sv_chrom", "")),
                "start": min(start, end),
                "end": max(start, end),
                "span_bp": int(span_bp),
                "span_label": format_span_bp(span_bp),
                "tumor_vaf": pd.to_numeric(getattr(row, "tumor_vaf", np.nan), errors="coerce"),
                "tumor_vaf_adj": pd.to_numeric(getattr(row, "tumor_vaf_adj", np.nan), errors="coerce"),
                "tumor_vaf_ratio_adj": pd.to_numeric(getattr(row, "tumor_vaf_ratio_adj", np.nan), errors="coerce"),
                "tumor_alt_support": pd.to_numeric(getattr(row, "tumor_alt_support", np.nan), errors="coerce"),
            }
        )
    events.sort(key=lambda e: (0 if e["svtype"] == "DEL" else 1, e["start"]))
    return events


def sv_vaf_label(ev: Dict[str, object]) -> str:
    parts: List[str] = []
    vaf = pd.to_numeric(ev.get("tumor_vaf"), errors="coerce")
    if np.isfinite(vaf):
        parts.append(f"VAF={float(vaf):.3f}")
    vadj = pd.to_numeric(ev.get("tumor_vaf_adj"), errors="coerce")
    if np.isfinite(vadj):
        parts.append(f"CCF={float(vadj):.3f}")
    ratio = pd.to_numeric(ev.get("tumor_vaf_ratio_adj"), errors="coerce")
    if np.isfinite(ratio):
        parts.append(f"ratio={float(ratio):.2f}")
    alt = pd.to_numeric(ev.get("tumor_alt_support"), errors="coerce")
    if not parts and np.isfinite(alt):
        parts.append(f"alt={int(alt)}")
    return " · ".join(parts)


def sv_events_summary_text(events: List[Dict[str, object]]) -> str:
    if not events:
        return ""
    chunks = []
    for e in events:
        chunk = f"{e['svtype']} {e['span_label']}"
        vaf_txt = sv_vaf_label(e)
        if vaf_txt:
            chunk += f" [{vaf_txt}]"
        chunk += f" ({e['sv_id'][:28]})"
        chunks.append(chunk)
    return "; ".join(chunks)


def _draw_hbar(
    ax,
    y: float,
    start: int,
    end: int,
    color: str,
    *,
    label_right: str = "",
    label_above: str = "",
    hatch: str = "",
    alpha: float = 0.85,
    height: float = 0.38,
) -> None:
    s, e = min(start, end), max(start, end)
    width = max(e - s, 1)
    ax.barh(
        y,
        width,
        left=s,
        height=height,
        color=color,
        alpha=alpha,
        edgecolor="0.25",
        linewidth=0.5,
        hatch=hatch,
    )
    if label_above:
        ax.text(
            (s + e) / 2,
            y + height / 2 + 0.06,
            label_above,
            ha="center",
            va="bottom",
            fontsize=7,
            clip_on=False,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", alpha=0.85, edgecolor="none"),
        )
    if label_right:
        ax.text(e, y, f"  {label_right}", va="center", ha="left", fontsize=6.5, clip_on=True)


def draw_segment_track(
    ax,
    segment_label: object,
    window: Tuple[str, int, int],
    y: float,
    *,
    copy_number: Optional[float] = None,
    cn_state: str = "",
    cpe: Optional[float] = None,
) -> Optional[Dict[str, object]]:
    seg = summarize_interval(segment_label)
    if not seg:
        return None
    clipped = _clip_interval(str(seg["chrom"]), int(seg["start"]), int(seg["end"]), window)
    truncated = False
    if clipped:
        cs, ce = clipped
        truncated = cs > int(seg["start"]) or ce < int(seg["end"])
    cn_txt = ""
    if copy_number is not None and np.isfinite(copy_number):
        cn_txt = f" · CN≈{float(copy_number):.1f}"
    if cn_state:
        cn_txt += f" ({cn_state})"
    if cpe is not None and np.isfinite(cpe):
        cn_txt += f" · CPE={float(cpe):.2f}"
    above = f"ASCAT3 segment {seg['span_label']}{cn_txt}"
    if truncated:
        above += " · bar clipped to window"
    if clipped:
        _draw_hbar(
            ax,
            y,
            clipped[0],
            clipped[1],
            "#bdbdbd",
            label_above=above,
            label_right=f"{seg['chrom']}:{seg['start']:,}-{seg['end']:,}",
            hatch="///" if truncated else "",
        )
    return seg


def draw_sv_event_stack(
    ax,
    events: List[Dict[str, object]],
    window: Tuple[str, int, int],
    *,
    y_top: float = 1.8,
    y_step: float = 0.42,
) -> float:
    """Draw one horizontal bar per SV; returns lowest y used."""
    y = y_top
    for ev in events:
        clipped = _clip_interval(str(ev["chrom"]), int(ev["start"]), int(ev["end"]), window)
        if not clipped:
            continue
        color = _SVTYPE_COLORS.get(str(ev["svtype"]), "#7f8c8d")
        vaf_txt = sv_vaf_label(ev)
        label = f"{ev['svtype']} {ev['span_label']}"
        if vaf_txt:
            label += f" · {vaf_txt}"
        label += f" · {ev['sv_id'][:28]}"
        _draw_hbar(ax, y, clipped[0], clipped[1], color, label_right=label)
        y -= y_step
    return y


def schematic_y_layout(n_sv_events: int) -> Tuple[float, float, float, float]:
    """Return (segment_y, sv_y_top, hairpin_y, y_max) for dynamic stacking."""
    sv_y_top = 0.35 + max(n_sv_events, 1) * 0.42
    segment_y = sv_y_top + 1.05
    hairpin_y = -0.15
    y_max = segment_y + 0.75
    return segment_y, sv_y_top, hairpin_y, y_max


def _parse_interval(label: object) -> Optional[Tuple[str, int, int]]:
    if label is None or (isinstance(label, float) and np.isnan(label)):
        return None
    m = _INTERVAL_RE.match(str(label).strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def _ucsc_url(window: str) -> str:
    if not window:
        return ""
    return f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={window}"


def _clip_interval(chrom: str, start: int, end: int, win: Tuple[str, int, int]) -> Optional[Tuple[int, int]]:
    if chrom != win[0]:
        return None
    s, e = min(start, end), max(start, end)
    ws, we = win[1], win[2]
    cs, ce = max(s, ws), min(e, we)
    if ce <= cs:
        return None
    return cs, ce


def plot_case_schematic(
    row: pd.Series,
    sv_rows: pd.DataFrame,
    *,
    out_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    win = _parse_interval(row.get("igv_window"))
    locus = _parse_interval(row.get("igv_locus"))
    if not win:
        return

    wchrom, w0, w1 = win
    span_mb = (w1 - w0) / 1e6
    events = summarize_sv_events(sv_rows)
    seg_y, sv_y_top, hairpin_y, y_max = schematic_y_layout(len(events))

    fig_h = 3.6 + max(0, len(events) - 1) * 0.35
    fig, ax = plt.subplots(figsize=(14, fig_h), dpi=140)
    ax.set_xlim(w0, w1)
    ax.set_ylim(-0.45, y_max)
    ax.set_xlabel(f"{wchrom} (bp)")
    ax.set_yticks([hairpin_y, (sv_y_top + 0.2) / 2 if events else 0.8, seg_y])
    ax.set_yticklabels(["Hairpin locus", f"SV events (n={len(events)})", "ASCAT3 segment"])
    ax.grid(axis="x", alpha=0.25, lw=0.5)

    cn = pd.to_numeric(row.get("copy_number"), errors="coerce")
    cstate = str(row.get("cn_state") or "")
    draw_segment_track(
        ax,
        row.get("overlap_segment"),
        win,
        seg_y,
        copy_number=float(cn) if np.isfinite(cn) else None,
        cn_state=cstate,
        cpe=pd.to_numeric(row.get("CPE"), errors="coerce") if "CPE" in row.index else None,
    )

    if events:
        draw_sv_event_stack(ax, events, win, y_top=sv_y_top)
    else:
        ax.text(w0, sv_y_top, "  (no overlapping SV)", fontsize=7, va="center")

    if locus:
        lh, ls, le = locus
        lname = str(row.get("locus_name", row.get("locus_id", "")))
        mid = (ls + le) / 2
        ax.vlines(mid, hairpin_y - 0.12, hairpin_y + 0.12, colors="#c0392b", linewidth=2.5, zorder=5)
        ax.text(mid, hairpin_y + 0.28, lname, fontsize=7, ha="center", color="#c0392b")

    sv_txt = sv_events_summary_text(events)
    if sv_txt:
        fig.text(0.01, 0.01, f"SVs: {sv_txt}", fontsize=6.5, color="0.3", wrap=True)

    title = (
        f"{row.get('review_queue', '')} | {row.get('participant', '')} | "
        f"{row.get('PAM50_final', '')} | CN={row.get('copy_number', '')} ({cstate}) | "
        f"layer={row.get('dosage_layer', '')}"
    )
    ax.set_title(title, fontsize=9, loc="left")
    note = f"View window {span_mb:.2f} Mb · {wchrom}:{w0}-{w1}"
    fig.text(0.01, 0.04, note, fontsize=7, color="0.35")
    legend_elems = [
        Line2D([0], [0], color=_SVTYPE_COLORS["DEL"], lw=3, label="DEL"),
        Line2D([0], [0], color=_SVTYPE_COLORS["DUP"], lw=3, label="DUP"),
    ]
    ax.legend(handles=legend_elems, loc="upper right", fontsize=7, framealpha=0.9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def _safe_slug(*parts: object) -> str:
    raw = "_".join(str(p) for p in parts if p is not None and str(p) != "nan")
    return re.sub(r"[^\w.-]+", "_", raw)[:120]


def write_case_plot_gallery(
    *,
    queue_dir: Optional[Path] = None,
    hits_path: Optional[Path] = None,
    out_dir: Optional[Path] = None,
    max_per_queue: Optional[int] = None,
    plot_all: bool = False,
    default_large_queue_cap: int = 40,
    queues: Optional[List[str]] = None,
) -> Dict[str, Path]:
    queue_dir = Path(queue_dir or C.MIRNA_LOCUS_SV_DIR / "review_queues")
    hits_path = Path(hits_path or C.MIRNA_LOCUS_SV_DIR / "mirna_locus_sv_hits.tsv.gz")
    out_dir = Path(out_dir or queue_dir / "case_plots")
    out_dir.mkdir(parents=True, exist_ok=True)

    hits = pd.read_csv(hits_path, sep="\t", low_memory=False)
    paths: Dict[str, Path] = {}
    manifest_rows: List[dict] = []

    queue_files = sorted(queue_dir.glob("review_queue_*.tsv.gz"))
    if queues:
        want = {f"review_queue_{q}.tsv.gz" for q in queues}
        queue_files = [p for p in queue_files if p.name in want]

    for qpath in queue_files:
        qname = qpath.name.replace("review_queue_", "").replace(".tsv.gz", "")
        df = pd.read_csv(qpath, sep="\t", low_memory=False)
        if df.empty:
            continue
        if not plot_all and max_per_queue is None and qname != "sv_only_focal_le_1mb":
            cap = default_large_queue_cap
        elif plot_all or max_per_queue is None:
            cap = None
        else:
            cap = max_per_queue

        if cap is not None and len(df) > cap:
            sort_col = "median_sv_span_mb" if "median_sv_span_mb" in df.columns else None
            if sort_col:
                df = df.sort_values(sort_col, na_position="last").head(cap)
            else:
                df = df.head(cap)

        qout = out_dir / qname
        qout.mkdir(parents=True, exist_ok=True)

        for _, row in df.iterrows():
            sub = hits.loc[
                (hits["participant"] == row["participant"]) & (hits["locus_id"].astype(str) == str(row["locus_id"]))
            ]
            slug = _safe_slug(row["participant"], row["locus_id"])
            png = qout / f"{slug}.png"
            events = summarize_sv_events(sub)
            seg = summarize_interval(row.get("overlap_segment"))
            plot_case_schematic(row, sub, out_path=png)
            n_del = sum(1 for e in events if e["svtype"] == "DEL")
            n_dup = sum(1 for e in events if e["svtype"] == "DUP")
            vafs = [e["tumor_vaf"] for e in events if np.isfinite(pd.to_numeric(e.get("tumor_vaf"), errors="coerce"))]
            ratios = [
                e["tumor_vaf_ratio_adj"]
                for e in events
                if np.isfinite(pd.to_numeric(e.get("tumor_vaf_ratio_adj"), errors="coerce"))
            ]
            med_vaf = float(np.median(vafs)) if vafs else np.nan
            med_ratio = float(np.median(ratios)) if ratios else np.nan
            manifest_rows.append(
                {
                    "review_queue": qname,
                    "participant": row["participant"],
                    "locus_id": row["locus_id"],
                    "locus_name": row.get("locus_name", ""),
                    "PAM50_final": row.get("PAM50_final", ""),
                    "dosage_layer": row.get("dosage_layer", ""),
                    "cn_state": row.get("cn_state", ""),
                    "CPE": row.get("CPE", np.nan),
                    "segment_span_label": seg.get("span_label", "") if seg else "",
                    "overlap_segment": row.get("overlap_segment", ""),
                    "sv_events_summary": sv_events_summary_text(events),
                    "median_tumor_vaf": round(med_vaf, 4) if np.isfinite(med_vaf) else np.nan,
                    "median_vaf_ratio_adj": round(med_ratio, 4) if np.isfinite(med_ratio) else np.nan,
                    "n_sv_del": n_del,
                    "n_sv_dup": n_dup,
                    "has_mixed_del_dup": bool(n_del and n_dup),
                    "igv_window": row.get("igv_window", ""),
                    "ucsc_url": _ucsc_url(str(row.get("igv_window", ""))),
                    "plot_png": str(png.relative_to(out_dir)),
                }
            )

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_dir / "case_plot_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    paths["manifest"] = manifest_path

    index_path = out_dir / "index.html"
    index_path.write_text(_render_html_index(manifest, out_dir), encoding="utf-8")
    paths["index"] = index_path

    meta = {
        "n_plots": int(len(manifest)),
        "queues": sorted(manifest["review_queue"].unique().tolist()) if not manifest.empty else [],
        "index_html": str(index_path),
    }
    meta_path = out_dir / "case_plot_gallery.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    paths["meta"] = meta_path

    print(f"[mirna_locus_sv_case_plots] {len(manifest)} schematics -> {out_dir}")
    return paths


def _render_html_index(manifest: pd.DataFrame, out_dir: Path) -> str:
    if manifest.empty:
        return "<html><body><p>No case plots generated.</p></body></html>"

    sections: List[str] = []
    for qname, grp in manifest.groupby("review_queue"):
        cards = []
        for row in grp.itertuples(index=False):
            cards.append(
                f"""
                <div class="card">
                  <img src="{html.escape(row.plot_png)}" alt="{html.escape(str(row.participant))}" loading="lazy"/>
                  <div class="meta">
                    <strong>{html.escape(str(row.participant))}</strong>
                    · {html.escape(str(row.PAM50_final))}
                    · {html.escape(str(row.locus_name or row.locus_id))}<br/>
                    ASCAT3 segment: {html.escape(str(getattr(row, 'segment_span_label', '') or '—'))}
                    · CN {html.escape(str(row.cn_state))}
                    · CPE {getattr(row, 'CPE', '')}
                    · layer {html.escape(str(row.dosage_layer))}<br/>
                    SV: {html.escape(str(getattr(row, 'sv_events_summary', '') or 'none'))}
                    {' · <b>mixed DEL+DUP</b>' if getattr(row, 'has_mixed_del_dup', False) else ''}
                    · med VAF={getattr(row, 'median_tumor_vaf', '')}<br/>
                    <a href="{html.escape(row.ucsc_url)}" target="_blank">UCSC hg38 (no BAM)</a>
                  </div>
                </div>
                """
            )
        sections.append(f"<h2>{html.escape(str(qname))} ({len(grp)})</h2><div class='grid'>{''.join(cards)}</div>")

    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"/>
<title>miRNA SV/CNV review case plots</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 1.5rem; background: #fafafa; }}
h1 {{ font-size: 1.25rem; }}
h2 {{ margin-top: 2rem; font-size: 1.05rem; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(420px, 1fr)); gap: 1rem; }}
.card {{ background: #fff; border: 1px solid #ddd; border-radius: 6px; overflow: hidden; }}
.card img {{ width: 100%; display: block; background: #fff; }}
.meta {{ padding: 0.6rem 0.75rem; font-size: 0.82rem; line-height: 1.35; }}
code {{ font-size: 0.75rem; }}
</style></head>
<body>
<h1>miRNA SV/CNV review — static case schematics (no BAM)</h1>
<p>ASCAT3 segment size is labeled above the gray bar. Each SV is on its own row
(red=DEL, blue=DUP) with genomic span. Folder: <code>{html.escape(str(out_dir))}</code></p>
{''.join(sections)}
</body></html>
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--queue-dir", type=Path, default=None)
    ap.add_argument("--hits", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument(
        "--max-per-queue",
        type=int,
        default=None,
        help="Cap plots per queue (overrides default: all focal + 40 each for other queues)",
    )
    ap.add_argument("--all", action="store_true", help="Plot every row in every queue (can be slow/large)")
    ap.add_argument("--queue", action="append", default=None, help="Restrict to queue stem, e.g. sv_only_focal_le_1mb")
    args = ap.parse_args()

    max_per = args.max_per_queue
    if max_per is None and not args.all:
        max_per = None  # gallery applies per-queue defaults

    write_case_plot_gallery(
        queue_dir=args.queue_dir,
        hits_path=args.hits,
        out_dir=args.out_dir,
        max_per_queue=args.max_per_queue,
        plot_all=args.all,
        queues=args.queue,
    )


if __name__ == "__main__":
    main()
