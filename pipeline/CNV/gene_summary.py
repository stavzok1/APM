"""
Gene-centric CNV summarization.

Takes per-segment CNV tables annotated with `gene_hits` (list[dict]) and derives a
per-(sample,gene) summary table with discrete states (loss/gain/amp/loh_only/etc.).
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class GeneCnvRules:
    """Heuristics for calling gene-level CNV states from segments."""

    min_body_overlap_pct: float = 30.0
    min_promoter_call: bool = True
    # When multiple segments hit a gene, use max score segment:
    # score = abs(cn_total-2) * (overlap_pct/100)


def _safe_parse_list(x: Any) -> List[dict]:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    if isinstance(x, list):
        return x
    s = str(x).strip()
    if not s or s == "[]":
        return []
    try:
        v = ast.literal_eval(s)
        return v if isinstance(v, list) else []
    except Exception:
        return []


def explode_gene_hits(cnv_df: pd.DataFrame) -> pd.DataFrame:
    """
    Explode a CNV segments DataFrame into one row per (segment × gene_hit).
    Expects columns: Chromosome, Start, End, cn_total, loh_flag, cn_state, gene_hits
    """
    df = cnv_df.copy()
    if "gene_hits" not in df.columns:
        return pd.DataFrame()
    df["gene_hits"] = df["gene_hits"].apply(_safe_parse_list)
    rows: List[Dict[str, Any]] = []
    for _, seg in df.iterrows():
        hits = seg["gene_hits"]
        if not hits:
            continue
        for h in hits:
            rows.append(
                {
                    "Chromosome": seg.get("Chromosome"),
                    "Start": int(seg.get("Start")) if pd.notna(seg.get("Start")) else None,
                    "End": int(seg.get("End")) if pd.notna(seg.get("End")) else None,
                    "cn_total": float(seg.get("cn_total")) if pd.notna(seg.get("cn_total")) else np.nan,
                    "cn_state": seg.get("cn_state"),
                    "loh_flag": int(seg.get("loh_flag")) if pd.notna(seg.get("loh_flag")) else 0,
                    "gene_name": h.get("gene_name"),
                    "gene_id": h.get("gene_id"),
                    "strand": h.get("strand"),
                    "signed_dist": h.get("signed_dist"),
                    "overlap_bp": h.get("overlap_bp", 0),
                    "overlap_percent": float(h.get("overlap_percent", 0.0) or 0.0),
                    "promoter_flag": int(h.get("promoter_flag", 0) or 0),
                    "gene_body_flag": int(h.get("gene_body_flag", 0) or 0),
                    "exon_flag": int(h.get("exon_flag", 0) or 0),
                    "intron_only_flag": int(h.get("intron_only_flag", 0) or 0),
                    "region_hit": h.get("region_hit", ""),
                }
            )
    return pd.DataFrame(rows)


def _state_from_cn_total(cn_total: float) -> str:
    if not np.isfinite(cn_total):
        return "NA"
    cn = int(round(cn_total))
    if cn <= 0:
        return "deep_del"
    if cn == 1:
        return "loss"
    if cn == 2:
        return "neutral"
    if cn == 3:
        return "gain"
    if cn >= 4:
        return "amp"
    return "other"


def summarize_gene_cnv(
    cnv_df: pd.DataFrame,
    *,
    rules: Optional[GeneCnvRules] = None,
) -> pd.DataFrame:
    """
    Produce per-gene CNV calls from a per-segment CNV table with gene_hits.
    """
    rules = rules or GeneCnvRules()
    long = explode_gene_hits(cnv_df)
    if long.empty:
        return pd.DataFrame()

    # Only meaningful body overlaps for main calls; promoter calls handled separately.
    long["body_ok"] = (long["gene_body_flag"] == 1) & (long["overlap_percent"] >= rules.min_body_overlap_pct)
    long["prom_ok"] = long["promoter_flag"] == 1

    # score: magnitude × coverage (coverage=overlap_percent)
    long["cn_delta"] = (long["cn_total"] - 2.0).abs()
    long["cnv_score"] = long["cn_delta"] * (long["overlap_percent"] / 100.0)

    out_rows: List[Dict[str, Any]] = []
    for gene, sub in long.groupby("gene_name", dropna=True):
        if not isinstance(gene, str) or not gene:
            continue

        # Candidate segments for body call
        body = sub.loc[sub["body_ok"]].copy()
        prom = sub.loc[sub["prom_ok"]].copy()

        # Determine best body segment
        best = None
        if not body.empty:
            best = body.sort_values("cnv_score", ascending=False).iloc[0]

        # LOH-only: overlap body and cn_total approx 2 but loh_flag set.
        loh_only = False
        if not body.empty:
            loh_only = bool(((body["loh_flag"] == 1) & (body["cn_total"].round() == 2)).any())

        # Promoter-specific calls (use strongest magnitude among promoter overlaps)
        prom_best = None
        if not prom.empty:
            prom = prom.assign(_prom_score=prom["cn_delta"] * 1.0)
            prom_best = prom.sort_values("_prom_score", ascending=False).iloc[0]

        state = "no_call"
        cn_total_call = np.nan
        if best is not None:
            cn_total_call = float(best["cn_total"])
            state = _state_from_cn_total(cn_total_call)
            if loh_only and state == "neutral":
                state = "loh_only"

        promoter_state = ""
        promoter_cn_total = np.nan
        if prom_best is not None:
            promoter_cn_total = float(prom_best["cn_total"])
            promoter_state = _state_from_cn_total(promoter_cn_total)

        out_rows.append(
            {
                "gene_name": gene,
                "gene_id": "" if sub["gene_id"].isna().all() else str(sub["gene_id"].dropna().iloc[0]),
                "cnv_state": state,
                "cn_total_call": cn_total_call,
                "loh_any": int((sub["loh_flag"] == 1).any()),
                "loh_only": int(loh_only),
                "max_overlap_percent": float(sub["overlap_percent"].max()),
                "max_overlap_bp": int(pd.to_numeric(sub["overlap_bp"], errors="coerce").fillna(0).max()),
                "best_cnv_score": float(best["cnv_score"]) if best is not None else 0.0,
                "promoter_state": promoter_state,
                "promoter_cn_total": promoter_cn_total,
                "n_segments_hitting": int(len(sub)),
                "n_body_segments": int(len(body)),
                "n_promoter_segments": int(len(prom)),
            }
        )

    out = pd.DataFrame(out_rows).sort_values(["cnv_state", "best_cnv_score"], ascending=[True, False])
    return out

