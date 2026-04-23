#!/usr/bin/env python3
"""
Build per-module sample counts and overlap statistics across TCGA BRCA omics sources.

Run from repo root:
    python analysis/sample_module_coverage.py [--out DIR]

See analysis/README.md for outputs and interpretation.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.config import PATHS  # noqa: E402
from pipeline.sample_ids import count_unique_tcga_sample_types, normalize_tcga_id  # noqa: E402

from analysis.clinical_omics_stratification import run_stratification  # noqa: E402
from analysis.sample_coverage.coverage_module_lists import (  # noqa: E402
    METADATA_MODULES,
    OMICS_MODULES,
    PARTICIPANT_EXTRA,
    PARTICIPANT_MODULES,
    PARTICIPANT_OMICS_ONLY,
    VIAL_MODULES,
)


def _norm_vial(raw: str) -> Optional[str]:
    tid = normalize_tcga_id(str(raw).strip())
    return tid.sample_vial or tid.sample


def _norm_participant(raw: str) -> Optional[str]:
    tid = normalize_tcga_id(str(raw).strip())
    return tid.participant


def tumor_sample_vial_from_gdc_manifest_row(row: pd.Series) -> Optional[str]:
    """
    Tumor sample vial from a GDC-style manifest row (SNV/SV/CNV samples.tsv).

    Prefer Tissue Type list aligned with comma-split Sample ID; else SNV loader rule.
    """
    sid = str(row.get("Sample ID", "")).strip()
    if not sid or sid.lower() == "nan":
        return None
    parts = [p.strip() for p in sid.split(",") if p.strip()]
    if not parts:
        return None
    tt_raw = str(row.get("Tissue Type", "")).strip()
    tt_parts = [p.strip() for p in tt_raw.split(",") if p.strip()]
    chosen: Optional[str] = None
    if tt_parts and len(tt_parts) == len(parts):
        for i, t in enumerate(tt_parts):
            tl = t.lower()
            if "tumor" in tl or "metastatic" in tl:
                chosen = parts[i]
                break
    if chosen is None:
        td = str(row.get("Tumor Descriptor", "")).strip()
        if td.startswith("Not") and len(parts) >= 2:
            chosen = parts[1]
        else:
            chosen = parts[0]
    return _norm_vial(chosen)


def load_gdc_manifest_all_raw_sample_ids(tsv: Path) -> List[str]:
    """
    All TCGA specimen tokens from GDC-style manifests (comma-separated Sample ID per row).
    Used for tumor vs normal (01/10/11) breakdown; not limited to tumor-only pick.
    """
    if not tsv.exists():
        return []
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if "Sample ID" not in df.columns:
        return []
    out: List[str] = []
    for raw in df["Sample ID"].dropna().astype(str):
        out.extend([p.strip() for p in raw.split(",") if p.strip()])
    return out


def load_wide_matrix_tcga_column_headers(path: Path) -> List[str]:
    """TCGA column names from a wide matrix header (same logic as header scan)."""
    hdr = _read_wide_matrix_header(path)
    if hdr is None:
        return []
    return [str(c).strip() for c in hdr.columns if str(c).strip().startswith("TCGA")]


def load_gdc_manifest_vials(tsv: Path) -> Set[str]:
    if not tsv.exists():
        return set()
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    out: Set[str] = set()
    for _, row in df.iterrows():
        v = tumor_sample_vial_from_gdc_manifest_row(row)
        if v:
            out.add(v)
    return out


def load_simple_sample_id_column(tsv: Path, col: str = "Sample ID") -> Set[str]:
    """RPPA / methylation manifests: one tumor vial per row."""
    if not tsv.exists():
        return set()
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if col not in df.columns:
        return set()
    out: Set[str] = set()
    for raw in df[col].astype(str):
        v = _norm_vial(raw)
        if v:
            out.add(v)
    return out


def load_simple_participant_id_column(tsv: Path, col: str = "Sample ID") -> Set[str]:
    """Participant-keyed manifests (e.g. HiCHIP sample columns are TCGA participants)."""
    if not tsv.exists():
        return set()
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if col not in df.columns:
        return set()
    out: Set[str] = set()
    for raw in df[col].astype(str):
        p = _norm_participant(raw)
        if p:
            out.add(p)
    return out


def load_normalized_sample_vials(tsv: Path, col: str = "sample_vial") -> Set[str]:
    if not tsv.exists():
        return set()
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if col not in df.columns:
        return set()
    out: Set[str] = set()
    for raw in df[col].dropna().astype(str):
        v = _norm_vial(raw)
        if v:
            out.add(v)
    return out


def _read_wide_matrix_header(path: Path) -> Optional[pd.DataFrame]:
    """
    Read only the header row. Try tab and comma (and optional index_col).

    Tab-first naive reads can "succeed" on comma-separated files by collapsing the whole
    header into a single column; we pick the parse with the most ``TCGA*`` column names.
    """
    if not path.exists():
        return None
    best_hdr: Optional[pd.DataFrame] = None
    best_score = 0
    for sep in ("\t", ","):
        for kwargs in ({}, {"index_col": 0}):
            try:
                hdr = pd.read_csv(path, sep=sep, nrows=0, low_memory=False, **kwargs)
            except Exception:
                continue
            score = sum(1 for c in hdr.columns if str(c).strip().startswith("TCGA"))
            if score > best_score:
                best_score = score
                best_hdr = hdr
    return best_hdr


def load_tcga_columns_from_wide_matrix(path: Path) -> Set[str]:
    """Wide matrix: any column header starting with TCGA (RNA / ATAC case-level, etc.)."""
    hdr = _read_wide_matrix_header(path)
    if hdr is None:
        return set()
    out: Set[str] = set()
    for c in hdr.columns:
        s = str(c).strip()
        if not s.startswith("TCGA"):
            continue
        v = _norm_vial(s)
        if v:
            out.add(v)
    return out


def vials_to_participants(vials: Iterable[str]) -> Set[str]:
    ps: Set[str] = set()
    for v in vials:
        p = _norm_participant(v)
        if p:
            ps.add(p)
    return ps


def load_immune_advanced_participants(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Normalized file uses sample_id at ...-01 granularity
    col = "sample_id" if "sample_id" in df.columns else df.columns[0]
    out: Set[str] = set()
    for raw in df[col].dropna().astype(str):
        p = _norm_participant(raw)
        if p:
            out.add(p)
    return out


def load_thornsson_participants(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    col = "TCGA Participant Barcode" if "TCGA Participant Barcode" in df.columns else df.columns[0]
    out: Set[str] = set()
    for raw in df[col].dropna().astype(str):
        p = _norm_participant(raw)
        if p:
            out.add(p)
    return out


def load_clinical_unified_participants(path: Path) -> Set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "participant" in df.columns:
        return {str(x).strip() for x in df["participant"].dropna() if str(x).strip().startswith("TCGA")}
    out: Set[str] = set()
    if "sample_id" in df.columns:
        for raw in df["sample_id"].dropna().astype(str):
            p = _norm_participant(raw)
            if p:
                out.add(p)
    return out


def collect_module_raw_tcga_id_lists() -> Dict[str, List[str]]:
    """
    Raw TCGA identifiers per module (all manifest tokens / matrix columns) for 01/10/11 typing.
    Mirrors ``analysis/sample_coverage/sample_type_breakdown.py`` sources.
    """
    ann = PATHS.annotations_dir
    return {
        "SNV": load_gdc_manifest_all_raw_sample_ids(ann / "SNV" / "samples.tsv"),
        "SV": load_gdc_manifest_all_raw_sample_ids(ann / "SV" / "samples.tsv"),
        "CNV": load_gdc_manifest_all_raw_sample_ids(ann / "CNV" / "samples.tsv"),
        "Methylation": load_gdc_manifest_all_raw_sample_ids(PATHS.methylation_sample_manifest),
        "RPPA": load_gdc_manifest_all_raw_sample_ids(ann / "rppa" / "samples.tsv"),
        "miRNA": load_gdc_manifest_all_raw_sample_ids(PATHS.mirna_samples_tsv)
        or load_wide_matrix_tcga_column_headers(PATHS.mirna_expression_tsv),
        "RNA": load_gdc_manifest_all_raw_sample_ids(PATHS.rna_samples_tsv)
        or load_normalized_sample_vials_list(PATHS.rna_expression_sample_metadata)
        or load_wide_matrix_tcga_column_headers(PATHS.rna_expression)
        or load_wide_matrix_tcga_column_headers(PATHS.rna_expression_raw),
        "ATAC": load_gdc_manifest_all_raw_sample_ids(PATHS.atac_samples_tsv)
        or load_normalized_sample_vials_list(PATHS.atac_case_level_sample_metadata)
        or load_wide_matrix_tcga_column_headers(PATHS.atac_case_level_matrix),
        "HLA": load_gdc_manifest_all_raw_sample_ids(PATHS.hla_samples_tsv)
        or load_normalized_sample_vials_list(PATHS.hla_types_normalized),
    }


def load_normalized_sample_vials_list(tsv: Path, col: str = "sample_vial") -> List[str]:
    if not tsv.exists():
        return []
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    if col not in df.columns:
        return []
    return [str(x).strip() for x in df[col].dropna().astype(str)]


def collect_vial_sets() -> Dict[str, Set[str]]:
    ann = PATHS.annotations_dir
    return {
        "SNV": load_gdc_manifest_vials(ann / "SNV" / "samples.tsv"),
        "SV": load_gdc_manifest_vials(ann / "SV" / "samples.tsv"),
        "CNV": load_gdc_manifest_vials(ann / "CNV" / "samples.tsv"),
        "Methylation": load_simple_sample_id_column(PATHS.methylation_sample_manifest),
        "RPPA": load_simple_sample_id_column(ann / "rppa" / "samples.tsv"),
        "miRNA": (
            load_simple_sample_id_column(PATHS.mirna_samples_tsv)
            or load_tcga_columns_from_wide_matrix(PATHS.mirna_expression_tsv)
        ),
        "RNA": (
            load_simple_sample_id_column(PATHS.rna_samples_tsv)
            or load_normalized_sample_vials(PATHS.rna_expression_sample_metadata)
            or load_tcga_columns_from_wide_matrix(PATHS.rna_expression)
            or load_tcga_columns_from_wide_matrix(PATHS.rna_expression_raw)
        ),
        "ATAC": (
            load_simple_sample_id_column(PATHS.atac_samples_tsv)
            or load_normalized_sample_vials(PATHS.atac_case_level_sample_metadata)
            or load_tcga_columns_from_wide_matrix(PATHS.atac_case_level_matrix)
        ),
        "HLA": (
            load_simple_sample_id_column(PATHS.hla_samples_tsv)
            or load_normalized_sample_vials(PATHS.hla_types_normalized)
        ),
    }


def collect_participant_sets(vial_by_mod: Dict[str, Set[str]]) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {}
    for m in VIAL_MODULES:
        out[m] = vials_to_participants(vial_by_mod.get(m, set()))
    # Participant-keyed assay modules
    out["HiCHIP"] = load_simple_participant_id_column(PATHS.hichip_samples_tsv)
    out["Immune_advanced"] = load_immune_advanced_participants(
        PATHS.immune_subtype_annotations_normalized
    )
    out["Immune_thornsson"] = load_thornsson_participants(PATHS.thornsson_immune_table_normalized)
    out["Clinical_unified"] = load_clinical_unified_participants(PATHS.brca_clinical_immune_unified)
    return out


def presence_matrix(ids: List[str], modules: List[str], sets_by_mod: Dict[str, Set[str]]) -> pd.DataFrame:
    data = {m: [sid in sets_by_mod.get(m, set()) for sid in ids] for m in modules}
    return pd.DataFrame(data, index=ids)


def intersection_table(df: pd.DataFrame, modules: List[str]) -> pd.DataFrame:
    """Rows = non-empty module combinations; cols pattern, n."""
    sub = df[list(modules)].astype(bool)
    active = sub.any(axis=1)
    sub = sub.loc[active]
    if sub.empty:
        return pd.DataFrame(columns=["pattern", "n"])
    def pat(row: pd.Series) -> str:
        return "+".join(m for m in modules if row[m])

    g = sub.assign(_p=sub.apply(pat, axis=1)).groupby("_p").size().reset_index(name="n")
    g = g.rename(columns={"_p": "pattern"}).sort_values("n", ascending=False)
    return g


def pairwise_table(sets: Dict[str, Set[str]], modules: List[str]) -> pd.DataFrame:
    rows = []
    for a, b in combinations(modules, 2):
        sa, sb = sets.get(a, set()), sets.get(b, set())
        inter = len(sa & sb)
        uni = len(sa | sb)
        j = inter / uni if uni else 0.0
        rows.append({"module_a": a, "module_b": b, "n_intersection": inter, "n_a": len(sa), "n_b": len(sb), "jaccard": round(j, 6)})
    return pd.DataFrame(rows)


def _symmetric_matrix(pair_df: pd.DataFrame, modules: List[str], value_col: str) -> np.ndarray:
    """Square matrix (len(modules),) with pairwise values on off-diagonal; diagonal = NaN."""
    n = len(modules)
    mat = np.full((n, n), np.nan, dtype=float)
    pos = {m: i for i, m in enumerate(modules)}
    for _, r in pair_df.iterrows():
        i, j = pos[r["module_a"]], pos[r["module_b"]]
        v = float(r[value_col])
        mat[i, j] = mat[j, i] = v
    return mat


def build_sample_type_breakdown_rows(raw_by_mod: Dict[str, List[str]]) -> pd.DataFrame:
    """Per-module counts for TCGA 01 / 10 / 11 / other (vial-level modules + participant-only rows)."""
    rows: List[dict] = []
    for m in VIAL_MODULES:
        ids = raw_by_mod.get(m, [])
        n_t, n01, n10, n11, no = count_unique_tcga_sample_types(ids)
        n_norm = n10 + n11
        rows.append(
            {
                "module": m,
                "n_unique_specimens": n_t,
                "n_primary_tumor_01": n01,
                "n_normal_blood_derived_10": n10,
                "n_normal_solid_tissue_11": n11,
                "n_normal_total_10_plus_11": n_norm,
                "n_other_or_unknown": no,
                "pct_tumor_01": round(100.0 * n01 / max(n_t, 1), 3) if n_t else 0.0,
                "pct_normal_blood_10_of_normals": round(100.0 * n10 / max(n_norm, 1), 3) if n_norm else 0.0,
                "pct_normal_solid_11_of_normals": round(100.0 * n11 / max(n_norm, 1), 3) if n_norm else 0.0,
                "notes": "",
            }
        )

    ann = PATHS.annotations_dir
    participant_only_sources = {
        "HiCHIP": load_simple_participant_id_column(PATHS.hichip_samples_tsv),
        "Immune_advanced": load_immune_advanced_participants(PATHS.immune_subtype_annotations_normalized),
        "Immune_thornsson": load_thornsson_participants(PATHS.thornsson_immune_table_normalized),
        "Clinical_unified": load_clinical_unified_participants(PATHS.brca_clinical_immune_unified),
    }
    for mod, pset in participant_only_sources.items():
        n_p = len(pset)
        rows.append(
            {
                "module": mod,
                "n_unique_specimens": n_p,
                "n_primary_tumor_01": n_p,
                "n_normal_blood_derived_10": 0,
                "n_normal_solid_tissue_11": 0,
                "n_normal_total_10_plus_11": 0,
                "n_other_or_unknown": 0,
                "pct_tumor_01": 100.0 if n_p else 0.0,
                "pct_normal_blood_10_of_normals": 0.0,
                "pct_normal_solid_11_of_normals": 0.0,
                "notes": "participant_only_assumed_primary_tumor",
            }
        )
    return pd.DataFrame(rows).sort_values("module")


def save_sample_type_coverage_figures(
    figures_dir: Path,
    tag: str,
    scope: str,
    breakdown: pd.DataFrame,
    overall: Tuple[int, int, int, int, int],
    *,
    vial_modules: List[str],
) -> None:
    """Stacked bar (per module) + overall union bar for tumor vs normal (10 vs 11)."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"Skipping sample-type figures ({scope}: {exc})")
        return

    figures_dir.mkdir(parents=True, exist_ok=True)
    d = breakdown[breakdown["module"].isin(vial_modules)].copy()
    if d.empty:
        return
    d = d.sort_values("module")
    labels = d["module"].tolist()
    c01 = d["n_primary_tumor_01"].to_numpy(dtype=float)
    c10 = d["n_normal_blood_derived_10"].to_numpy(dtype=float)
    c11 = d["n_normal_solid_tissue_11"].to_numpy(dtype=float)
    co = d["n_other_or_unknown"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(10, max(5.0, 0.38 * len(labels))))
    y = np.arange(len(labels))
    ax.barh(y, c01, label="01 primary tumor", color="#c44e52")
    left = c01
    ax.barh(y, c10, left=left, label="10 blood-derived normal", color="#4c72b0")
    left = left + c10
    ax.barh(y, c11, left=left, label="11 solid-tissue normal", color="#55a868")
    left = left + c11
    ax.barh(y, co, left=left, label="other / unknown", color="#8172b2")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Unique specimens (sample_vial preferred else sample)")
    ax.set_title(f"[{scope}] TCGA sample-type mix per module")
    ax.legend(loc="lower right", fontsize=8)
    plt.tight_layout()
    fig.savefig(figures_dir / f"fig_sample_type_stacked_per_module_{scope}_{tag}.png", dpi=150)
    plt.close(fig)

    n_t, n01, n10, n11, no = overall
    fig2, ax2 = plt.subplots(figsize=(8.5, 2.8))
    ax2.barh([0], [n01], height=0.45, label="01 primary tumor", color="#c44e52")
    x0 = n01
    ax2.barh([0], [n10], height=0.45, left=[x0], label="10 blood-derived normal", color="#4c72b0")
    x0 += n10
    ax2.barh([0], [n11], height=0.45, left=[x0], label="11 solid-tissue normal", color="#55a868")
    x0 += n11
    ax2.barh([0], [no], height=0.45, left=[x0], label="other / unknown", color="#8172b2")
    ax2.set_yticks([0])
    ax2.set_yticklabels([f"Union (omics vial modules, N={n_t} unique specimens)"])
    ax2.set_xlabel("Count")
    ax2.set_title(f"[{scope}] Overall tumor vs normal mix (unique across modules)")
    ax2.legend(loc="lower center", ncol=4, fontsize=7, bbox_to_anchor=(0.5, -0.38))
    plt.tight_layout()
    fig2.savefig(figures_dir / f"fig_sample_type_overall_union_{scope}_{tag}.png", dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"Saved sample-type figures under {figures_dir}")


def write_module_notes(notes_dir: Path) -> None:
    notes = pd.DataFrame(
        [
            {
                "module": "SCREEN_ABC_ChIP",
                "scope": "cell_line_reference",
                "note": "cCRE / SCREEN / ABC / ChIP evidence is reference and cell-line keyed; omitted from TCGA sample overlap tables.",
            },
        ]
    )
    notes_dir.mkdir(parents=True, exist_ok=True)
    notes.to_csv(notes_dir / "module_notes.tsv", sep="\t", index=False)


def _heatmap(
    ax,
    mat: np.ndarray,
    row_labels: List[str],
    col_labels: List[str],
    *,
    vmin: Optional[float],
    vmax: Optional[float],
    cmap: str,
    cbar_label: str,
    annotate: bool = False,
    fmt: str = ".2f",
    annotate_fontsize: Optional[float] = None,
) -> None:
    im = ax.imshow(mat, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=7)
    if annotate:
        ncell = max(mat.shape[0], mat.shape[1])
        fs = annotate_fontsize if annotate_fontsize is not None else max(4.0, min(9.0, 220.0 / max(ncell, 1)))
        vhi = vmax if vmax is not None and np.isfinite(vmax) else float(np.nanmax(mat))
        lo = vmin if vmin is not None and np.isfinite(vmin) else float(np.nanmin(mat))
        span = (vhi - lo) if np.isfinite(vhi) and np.isfinite(lo) and (vhi > lo) else 0.0
        mid = lo + 0.55 * span if span > 0 and np.isfinite(span) else (0.55 * vhi if np.isfinite(vhi) else 0.5)
        pe = None
        try:  # matplotlib is available if we're here; this is for better contrast
            from matplotlib import patheffects as pe  # type: ignore
        except Exception:
            pe = None
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isnan(v):
                    continue
                s = f"{v:{fmt}}"
                t = ax.text(j, i, s, ha="center", va="center", color="white" if v > mid else "black", fontsize=fs)
                if pe is not None:
                    # Outline keeps text readable on bright colormap regions (e.g. high Jaccard ~ yellow).
                    t.set_path_effects([pe.withStroke(linewidth=1.25, foreground="black")])
    ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=cbar_label)


def save_figure_bundle(
    out_dir: Path,
    key: str,
    scope: str,
    *,
    include_vials: bool,
    figures_dir: Path,
    vial_counts: pd.Series,
    participant_counts: pd.Series,
    pair_vial: pd.DataFrame,
    pair_participant: pd.DataFrame,
    inter_vial: pd.DataFrame,
    inter_participant: pd.DataFrame,
    df_vial: pd.DataFrame,
    df_participant: pd.DataFrame,
) -> None:
    """scope is 'omics' or 'metadata' — used in output filenames before the run key."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover
        print(f"Skipping plots ({scope}: {exc})")
        return

    tag = f"{scope}_{key}"
    figures_dir.mkdir(parents=True, exist_ok=True)

    if include_vials and not vial_counts.empty:
        fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(vial_counts))))
        vial_counts.sort_values().plot.barh(ax=ax, color="steelblue")
        ax.set_title(f"[{scope}] Unique tumor sample_vials per module")
        ax.set_xlabel("Count")
        plt.tight_layout()
        fig.savefig(figures_dir / f"fig_module_counts_sample_vial_{tag}.png", dpi=150)
        plt.close(fig)

    if not participant_counts.empty:
        fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(participant_counts))))
        participant_counts.sort_values().plot.barh(ax=ax, color="darkseagreen")
        ax.set_title(f"[{scope}] Unique participants per module")
        ax.set_xlabel("Count")
        plt.tight_layout()
        fig.savefig(figures_dir / f"fig_module_counts_participant_{tag}.png", dpi=150)
        plt.close(fig)

    for label, df, fname in (
        ("sample_vial", df_vial, "fig_n_modules_per_sample_vial"),
        ("participant", df_participant, "fig_n_modules_per_participant"),
    ):
        if not include_vials and label == "sample_vial":
            continue
        if df.empty:
            continue
        n_hit = df.astype(bool).sum(axis=1)
        if n_hit.max() == 0:
            continue
        fig, ax = plt.subplots(figsize=(8, 4.5))
        max_k = int(n_hit.max()) + 1
        bins = np.arange(0.5, max_k + 1.5, 1.0)
        ax.hist(n_hit, bins=bins, color="teal", edgecolor="white")
        ax.set_xticks(range(1, max_k + 1))
        ax.set_xlabel("Number of modules with data")
        ax.set_ylabel(f"Number of {label}s")
        ax.set_title(f"[{scope}] Distribution of module coverage per {label}")
        plt.tight_layout()
        fig.savefig(figures_dir / f"{fname}_{tag}.png", dpi=150)
        plt.close(fig)

    if not pair_participant.empty:
        mod_p = sorted(set(pair_participant["module_a"]) | set(pair_participant["module_b"]))
        mat_j_p = _symmetric_matrix(pair_participant, mod_p, "jaccard")
        fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(mod_p)), max(8, 0.4 * len(mod_p))))
        _heatmap(
            ax,
            mat_j_p,
            mod_p,
            mod_p,
            vmin=0,
            vmax=1,
            cmap="viridis",
            cbar_label="Jaccard",
            annotate=True,
            fmt=".3f",
        )
        ax.set_title(f"[{scope}] Pairwise Jaccard (participant sets)")
        plt.tight_layout()
        fig.savefig(figures_dir / f"fig_jaccard_participant_{tag}.png", dpi=150)
        plt.close(fig)

        mat_n_p = _symmetric_matrix(pair_participant, mod_p, "n_intersection")
        fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(mod_p)), max(8, 0.4 * len(mod_p))))
        vmax = np.nanmax(mat_n_p) if np.isfinite(np.nanmax(mat_n_p)) else 1.0
        _heatmap(ax, mat_n_p, mod_p, mod_p, vmin=0, vmax=vmax, cmap="magma", cbar_label="|A ∩ B|", annotate=True, fmt=".0f")
        ax.set_title(f"[{scope}] Pairwise intersection counts (participants)")
        plt.tight_layout()
        fig.savefig(figures_dir / f"fig_pairwise_intersection_participant_{tag}.png", dpi=150)
        plt.close(fig)

    if include_vials and not pair_vial.empty:
        mod_v = sorted(set(pair_vial["module_a"]) | set(pair_vial["module_b"]))
        mat_j_v = _symmetric_matrix(pair_vial, mod_v, "jaccard")
        fig, ax = plt.subplots(figsize=(max(9, 0.5 * len(mod_v)), max(7, 0.45 * len(mod_v))))
        _heatmap(
            ax,
            mat_j_v,
            mod_v,
            mod_v,
            vmin=0,
            vmax=1,
            cmap="cividis",
            cbar_label="Jaccard",
            annotate=True,
            fmt=".3f",
        )
        ax.set_title(f"[{scope}] Pairwise Jaccard (sample_vial sets)")
        plt.tight_layout()
        fig.savefig(figures_dir / f"fig_jaccard_sample_vial_{tag}.png", dpi=150)
        plt.close(fig)

        mat_n_v = _symmetric_matrix(pair_vial, mod_v, "n_intersection")
        fig, ax = plt.subplots(figsize=(max(9, 0.5 * len(mod_v)), max(7, 0.45 * len(mod_v))))
        vmax = np.nanmax(mat_n_v) if np.isfinite(np.nanmax(mat_n_v)) else 1.0
        _heatmap(ax, mat_n_v, mod_v, mod_v, vmin=0, vmax=vmax, cmap="inferno", cbar_label="|A ∩ B|", annotate=True, fmt=".0f")
        ax.set_title(f"[{scope}] Pairwise intersection counts (sample_vials)")
        plt.tight_layout()
        fig.savefig(figures_dir / f"fig_pairwise_intersection_sample_vial_{tag}.png", dpi=150)
        plt.close(fig)

    for inter_df, title_suffix, fname in (
        (inter_participant, "participants", "fig_top_intersection_patterns_participant"),
        (inter_vial, "sample_vials", "fig_top_intersection_patterns_sample_vial"),
    ):
        if not include_vials and title_suffix == "sample_vials":
            continue
        if inter_df.empty:
            continue
        head = inter_df.head(30).iloc[::-1]
        fig, ax = plt.subplots(figsize=(11, max(6, 0.28 * len(head))))
        ax.barh(head["pattern"], head["n"], color="slategray")
        ax.set_title(f"[{scope}] Top intersection patterns ({title_suffix})")
        ax.set_xlabel("Count")
        plt.tight_layout()
        fig.savefig(figures_dir / f"{fname}_{tag}.png", dpi=150)
        plt.close(fig)

    print(f"Saved {scope} figures under {figures_dir} (matplotlib)")


def main() -> None:
    parser = argparse.ArgumentParser(description="Sample coverage and intersections across modules.")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: analysis/sample_coverage/output/run_<UTC timestamp>)",
    )
    parser.add_argument(
        "--write-current",
        action="store_true",
        help="Also mirror outputs to analysis/sample_coverage/output/current/ for stable downstream joins.",
    )
    args = parser.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    out_dir = args.out or (ROOT / "analysis" / "sample_coverage" / "output" / f"run_{stamp}")
    out_dir.mkdir(parents=True, exist_ok=True)
    key = out_dir.name.replace(" ", "_")

    tables_dir = out_dir / "tables"
    tables_omics = tables_dir / "omics"
    tables_meta = tables_dir / "metadata"
    figures_dir = out_dir / "figures"
    figures_omics = figures_dir / "omics"
    figures_meta = figures_dir / "metadata"
    notes_dir = out_dir / "notes"
    for p in (tables_omics, tables_meta, figures_omics, figures_meta, notes_dir):
        p.mkdir(parents=True, exist_ok=True)

    vial_sets = collect_vial_sets()
    participant_sets = collect_participant_sets(vial_sets)

    all_vials = sorted(set().union(*[vial_sets.get(m, set()) for m in VIAL_MODULES]))
    all_participants = sorted(set().union(*[participant_sets.get(m, set()) for m in PARTICIPANT_MODULES]))

    df_vial = presence_matrix(all_vials, VIAL_MODULES, vial_sets)
    df_vial.index.name = "sample_vial"
    df_vial.to_csv(tables_omics / "sample_vial_presence.tsv", sep="\t")

    df_p = presence_matrix(all_participants, PARTICIPANT_MODULES, participant_sets)
    df_p.index.name = "participant"
    df_p.to_csv(tables_dir / "participant_presence.tsv", sep="\t")
    df_p[OMICS_MODULES].to_csv(tables_omics / "participant_presence_omics.tsv", sep="\t")
    df_p[METADATA_MODULES].to_csv(tables_meta / "participant_presence_metadata.tsv", sep="\t")

    raw_by_mod = collect_module_raw_tcga_id_lists()
    sample_type_df = build_sample_type_breakdown_rows(raw_by_mod)
    sample_type_df.to_csv(tables_omics / "module_sample_type_breakdown.tsv", sep="\t", index=False)

    union_ids: List[str] = []
    for m in VIAL_MODULES:
        union_ids.extend(raw_by_mod.get(m, []))
    n_u, n01_u, n10_u, n11_u, no_u = count_unique_tcga_sample_types(union_ids)
    pd.DataFrame(
        [
            {
                "scope": "union_all_omics_vial_modules",
                "n_unique_specimens": n_u,
                "n_primary_tumor_01": n01_u,
                "n_normal_blood_derived_10": n10_u,
                "n_normal_solid_tissue_11": n11_u,
                "n_normal_total_10_plus_11": n10_u + n11_u,
                "n_other_or_unknown": no_u,
            }
        ]
    ).to_csv(tables_omics / "overall_sample_type_union_omics.tsv", sep="\t", index=False)

    st_lookup = sample_type_df.set_index("module")
    counts_omics = []
    for m in OMICS_MODULES:
        n_v = len(vial_sets.get(m, set())) if m in VIAL_MODULES else ""
        row = {
            "module": m,
            "n_sample_vial": n_v,
            "n_participant": len(participant_sets.get(m, set())),
        }
        if m in st_lookup.index:
            r = st_lookup.loc[m]
            row.update(
                {
                    "n_primary_tumor_01": int(r["n_primary_tumor_01"]),
                    "n_normal_blood_derived_10": int(r["n_normal_blood_derived_10"]),
                    "n_normal_solid_tissue_11": int(r["n_normal_solid_tissue_11"]),
                    "n_normal_total_10_plus_11": int(r["n_normal_total_10_plus_11"]),
                    "n_other_or_unknown": int(r["n_other_or_unknown"]),
                    "sample_type_notes": str(r.get("notes", "")),
                }
            )
        counts_omics.append(row)
    pd.DataFrame(counts_omics).to_csv(tables_omics / "module_counts_omics.tsv", sep="\t", index=False)

    counts_meta = []
    for m in METADATA_MODULES:
        counts_meta.append(
            {
                "module": m,
                "n_sample_vial": "",
                "n_participant": len(participant_sets.get(m, set())),
            }
        )
    pd.DataFrame(counts_meta).to_csv(tables_meta / "module_counts_metadata.tsv", sep="\t", index=False)

    # Legacy combined table (omics + metadata rows)
    pd.concat(
        [pd.DataFrame(counts_omics), pd.DataFrame(counts_meta)],
        ignore_index=True,
    ).to_csv(tables_dir / "module_counts.tsv", sep="\t", index=False)

    inter_v = intersection_table(df_vial, VIAL_MODULES)
    inter_v.to_csv(tables_omics / "intersections_sample_vial.tsv", sep="\t", index=False)

    inter_p_omics = intersection_table(df_p, OMICS_MODULES)
    inter_p_omics.to_csv(tables_omics / "intersections_participant_omics.tsv", sep="\t", index=False)
    inter_p_meta = intersection_table(df_p, METADATA_MODULES)
    inter_p_meta.to_csv(tables_meta / "intersections_participant_metadata.tsv", sep="\t", index=False)

    pair_p_omics = pairwise_table(participant_sets, OMICS_MODULES)
    pair_p_omics.to_csv(tables_omics / "pairwise_overlap_participant_omics.tsv", sep="\t", index=False)
    pair_p_meta = pairwise_table(participant_sets, METADATA_MODULES)
    pair_p_meta.to_csv(tables_meta / "pairwise_overlap_participant_metadata.tsv", sep="\t", index=False)

    pair_v = pairwise_table(vial_sets, VIAL_MODULES)
    pair_v.to_csv(tables_omics / "pairwise_overlap_sample_vial.tsv", sep="\t", index=False)

    write_module_notes(notes_dir)

    save_sample_type_coverage_figures(
        figures_omics,
        key,
        "omics",
        sample_type_df,
        (n_u, n01_u, n10_u, n11_u, no_u),
        vial_modules=list(VIAL_MODULES),
    )

    # Stratify coverage by PAM50 / collapsed stage within the same run folder.
    try:
        strat_tables = tables_omics / "stratification"
        strat_figures = figures_omics / "stratification"
        run_stratification(
            presence_path=tables_omics / "participant_presence_omics.tsv",
            out_tables_dir=strat_tables,
            out_figures_dir=strat_figures,
            clinical_path=None,
        )
    except Exception as exc:
        print(f"Skipping stratification (clinical_omics_stratification): {exc}")

    vial_series = pd.Series({m: len(vial_sets.get(m, set())) for m in VIAL_MODULES})
    part_series_omics = pd.Series({m: len(participant_sets.get(m, set())) for m in OMICS_MODULES})
    part_series_meta = pd.Series({m: len(participant_sets.get(m, set())) for m in METADATA_MODULES})

    save_figure_bundle(
        out_dir,
        key,
        "omics",
        include_vials=True,
        figures_dir=figures_omics,
        vial_counts=vial_series,
        participant_counts=part_series_omics,
        pair_vial=pair_v,
        pair_participant=pair_p_omics,
        inter_vial=inter_v,
        inter_participant=inter_p_omics,
        df_vial=df_vial,
        df_participant=df_p[OMICS_MODULES],
    )
    save_figure_bundle(
        out_dir,
        key,
        "metadata",
        include_vials=False,
        figures_dir=figures_meta,
        vial_counts=pd.Series(dtype=float),
        participant_counts=part_series_meta,
        pair_vial=pd.DataFrame(),
        pair_participant=pair_p_meta,
        inter_vial=pd.DataFrame(),
        inter_participant=inter_p_meta,
        df_vial=pd.DataFrame(),
        df_participant=df_p[METADATA_MODULES],
    )

    print(f"Wrote outputs to {out_dir}")
    print(
        "Key outputs under tables/{omics,metadata}/ and figures/{omics,metadata}/ (plus notes/)."
    )

    if args.write_current:
        import shutil

        current = ROOT / "analysis" / "sample_coverage" / "output" / "current"
        current.parent.mkdir(parents=True, exist_ok=True)
        shutil.rmtree(current, ignore_errors=True)
        shutil.copytree(out_dir, current)
        print(f"Mirrored outputs to canonical folder: {current}")


if __name__ == "__main__":
    main()
