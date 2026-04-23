from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from pipeline.sample_ids import normalize_tcga_id


def _parse_hits_cell(x) -> List[Dict]:
    if x is None:
        return []
    if isinstance(x, list):
        return x
    if isinstance(x, str):
        s = x.strip()
        if s in ("", "[]", "nan", "None"):
            return []
        try:
            v = ast.literal_eval(s)
            return v if isinstance(v, list) else []
        except Exception:
            return []
    return []


def _sv_sample_vial_from_filename(path: Path) -> Optional[str]:
    stem = path.name.split("_strict_sv_set")[0]
    ids = normalize_tcga_id(stem)
    v = ids.sample_vial or ids.sample
    return str(v) if v else None


def _hit_is_breakpoint_like(hit: Dict) -> bool:
    side = str(hit.get("hit_side", "")).lower()
    return side in {"bp1", "bp2", "point"}


def _hit_is_span_like(hit: Dict) -> bool:
    side = str(hit.get("hit_side", "")).lower()
    return side == "span"


def _flag(hit: Dict, key: str) -> bool:
    v = hit.get(key, 0)
    if isinstance(v, bool):
        return v
    try:
        return int(v) == 1
    except Exception:
        return False


@dataclass(frozen=True)
class SvGeneSummaryConfig:
    sv_output_root: Path
    stage_dir: str = "07_final_sv_with_fimo"
    sv_glob: str = "*_strict_sv_set.csv"


def summarize_sv_gene_hits(
    *,
    config: SvGeneSummaryConfig,
    genes_of_interest: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """
    General-purpose SV→gene “mirroring”:
    produces a sample_vial-indexed wide table of gene-centric disruption summaries.

    If genes_of_interest is None, the output includes all genes observed in gene_hits
    across the SV tables (can be wide).
    """
    d = Path(config.sv_output_root) / config.stage_dir
    genes_filter = set(g.upper() for g in genes_of_interest) if genes_of_interest else None

    rows: List[Dict[str, object]] = []
    for csv in sorted(d.glob(config.sv_glob)):
        sid = _sv_sample_vial_from_filename(csv)
        if not sid:
            continue
        try:
            sv = pd.read_csv(csv, low_memory=False)
        except Exception:
            continue
        if "gene_hits" not in sv.columns:
            continue

        per_gene: Dict[str, Dict[str, int]] = {}
        genes_seen: set[str] = set()
        global_stats = {
            "n_hits": 0,
            "n_bp_hits": 0,
            "n_span_hits": 0,
            "n_promoter_hits": 0,
            "n_cds_hits": 0,
            "n_cds_bp_hits": 0,
            "n_exon_hits": 0,
            "n_exon_bp_hits": 0,
            "n_start_stop_hits": 0,
        }

        for cell in sv["gene_hits"]:
            for hit in _parse_hits_cell(cell):
                g = str(hit.get("gene_name", "")).strip()
                if not g:
                    continue
                gU = g.upper()
                genes_seen.add(gU)

                is_bp = _hit_is_breakpoint_like(hit)
                is_span = _hit_is_span_like(hit)

                promoter = _flag(hit, "promoter_flag")
                cds = _flag(hit, "cds_flag")
                exon = _flag(hit, "exon_flag")
                start_stop = _flag(hit, "start_codon_flag") or _flag(hit, "stop_codon_flag")

                global_stats["n_hits"] += 1
                if is_bp:
                    global_stats["n_bp_hits"] += 1
                if is_span:
                    global_stats["n_span_hits"] += 1
                if promoter:
                    global_stats["n_promoter_hits"] += 1
                if cds:
                    global_stats["n_cds_hits"] += 1
                    if is_bp:
                        global_stats["n_cds_bp_hits"] += 1
                if exon:
                    global_stats["n_exon_hits"] += 1
                    if is_bp:
                        global_stats["n_exon_bp_hits"] += 1
                if start_stop:
                    global_stats["n_start_stop_hits"] += 1

                if genes_filter is not None and gU not in genes_filter:
                    continue

                bucket = per_gene.setdefault(
                    gU,
                    {
                        "n_hits": 0,
                        "n_bp_hits": 0,
                        "n_promoter_hits": 0,
                        "n_cds_hits": 0,
                        "n_cds_bp_hits": 0,
                        "n_exon_hits": 0,
                        "n_exon_bp_hits": 0,
                        "n_start_stop_hits": 0,
                    },
                )
                bucket["n_hits"] += 1
                if is_bp:
                    bucket["n_bp_hits"] += 1
                if promoter:
                    bucket["n_promoter_hits"] += 1
                if cds:
                    bucket["n_cds_hits"] += 1
                    if is_bp:
                        bucket["n_cds_bp_hits"] += 1
                if exon:
                    bucket["n_exon_hits"] += 1
                    if is_bp:
                        bucket["n_exon_bp_hits"] += 1
                if start_stop:
                    bucket["n_start_stop_hits"] += 1

        row: Dict[str, object] = {"sample_vial": sid}
        row["SV_any_gene__n_unique_genes_hit"] = int(len(genes_seen))
        for k, v in global_stats.items():
            row[f"SV_any_gene__{k}"] = int(v)

        for gU, stats in per_gene.items():
            row[f"SV_{gU}__any_hit"] = stats["n_hits"] > 0
            row[f"SV_{gU}__any_promoter_hit"] = stats["n_promoter_hits"] > 0
            row[f"SV_{gU}__any_cds_hit"] = stats["n_cds_hits"] > 0
            row[f"SV_{gU}__any_cds_bp_hit"] = stats["n_cds_bp_hits"] > 0
            row[f"SV_{gU}__any_exon_hit"] = stats["n_exon_hits"] > 0
            row[f"SV_{gU}__any_exon_bp_hit"] = stats["n_exon_bp_hits"] > 0
            for kk, vv in stats.items():
                row[f"SV_{gU}__{kk}"] = int(vv)

        rows.append(row)

    out = pd.DataFrame(rows).set_index("sample_vial") if rows else pd.DataFrame(index=pd.Index([], name="sample_vial"))
    out.index = out.index.astype(str)
    out.index.name = "sample_vial"
    return out


def write_sv_gene_summary(
    df: pd.DataFrame,
    *,
    out_path: Path,
    write_csv: bool = True,
) -> Tuple[Path, Optional[Path]]:
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq = out_path.with_suffix(".parquet")
    df.to_parquet(pq)
    csv = None
    if write_csv:
        csv = out_path.with_suffix(".csv")
        df.to_csv(csv)
    return pq, csv

