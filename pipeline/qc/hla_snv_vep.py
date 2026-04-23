"""
Per-sample HLA-focused SNV summaries from VEP ``gene_hits`` (already filtered to panel genes).

This is **not** copy-number LOH; it flags likely **allelic HLA loss / damage via somatic coding
variants** (LoF / splice / high-impact missense) on ``HLA-A`` / ``HLA-B`` / ``HLA-C`` / ``B2M``.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import numpy as np
import pandas as pd

from pipeline.CNV.loader import extract_sample_id_from_annotations
from pipeline.config import PATHS

HLA_PANEL_GENES: Set[str] = {"HLA-A", "HLA-B", "HLA-C", "B2M"}

# --- Thresholds for “counts toward HLA SNV damage” (combined with CN LOH in QC) ---
# Max population allele frequency (max over gnomADe/g, MAX_AF, AF when present).
# ≤ this value, or **no** population AF fields parsed, counts as eligible (somatic tables often omit AF).
HLA_SNV_MAX_POP_AF_FOR_DAMAGE = 0.01
# Missense: require at least one score when any pathogenicity field is present; otherwise allow HIGH missense alone.
HLA_SNV_CADD_PHRED_MIN = 20.0
HLA_SNV_REVEL_MIN = 0.5
HLA_SNV_SPLICEAI_DS_MIN = 0.5

LOF_CONSEQUENCES: Set[str] = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "transcript_ablation",
}

MISSENSE_TERM = "missense_variant"


def build_tumor_vial_to_snv_per_sample_csv(*, snv_output_dir: Path, snv_manifest_path: Path) -> Dict[str, Path]:
    """
    Map tumor vial id → ``.../per_sample/TCGA-BRCA.<uuid>*_snv_variants.csv`` using the GDC SNV manifest.
    """
    ann = pd.read_csv(snv_manifest_path, sep="\t", low_memory=False)
    per = Path(snv_output_dir) / "per_sample"
    if not per.is_dir():
        return {}
    out: Dict[str, Path] = {}
    uuid_re = re.compile(r"TCGA-BRCA\.([0-9a-f-]{36})\.", flags=re.I)
    for _, row in ann.iterrows():
        fn = str(row.get("File Name", "") or "")
        if "mutect2" not in fn.lower():
            continue
        m = uuid_re.search(fn)
        if not m:
            continue
        uuid = m.group(1)
        try:
            vial = extract_sample_id_from_annotations(ann, fn)
        except Exception:
            continue
        vial = str(vial).strip()
        if not vial:
            continue
        hits = sorted(per.glob(f"TCGA-BRCA.{uuid}*_snv_variants.csv"))
        if hits:
            out[vial] = hits[0]
    return out


def _parse_gene_hits_cell(v: object) -> List[dict]:
    if v is None:
        return []
    if isinstance(v, list):
        return [x for x in v if isinstance(x, dict)]
    if not isinstance(v, str):
        return []
    s = v.strip()
    if not s or s == "[]":
        return []
    try:
        out = ast.literal_eval(s)
        if isinstance(out, list):
            return [x for x in out if isinstance(x, dict)]
    except Exception:
        return []
    return []


def _consequence_terms(hit: Dict[str, Any]) -> Set[str]:
    cons = str(hit.get("Consequence") or "")
    return {t.strip() for t in cons.split("&") if t.strip()}


def _hit_lof_like(hit: Dict[str, Any]) -> bool:
    terms = _consequence_terms(hit)
    if terms & LOF_CONSEQUENCES:
        return True
    lof = str(hit.get("LoF") or "").upper()
    if lof in {"HC", "LC", "OS"}:  # LoFtee-style coarse flags when present
        return True
    return False


def _hit_high_missense(hit: Dict[str, Any]) -> bool:
    imp = str(hit.get("IMPACT") or "").upper()
    terms = _consequence_terms(hit)
    return imp == "HIGH" and MISSENSE_TERM in terms


def _safe_float(x: object) -> Optional[float]:
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return None
        if isinstance(x, str) and (not x.strip() or x.strip() in (".", "-")):
            return None
        return float(x)
    except Exception:
        return None


def _gnomad_max(hit: Dict[str, Any]) -> Optional[float]:
    vals = []
    for k in ("gnomADe_AF", "gnomADg_AF", "gnomAD_AF", "MAX_AF", "AF"):
        v = _safe_float(hit.get(k))
        if v is not None:
            vals.append(v)
    return max(vals) if vals else None


def _pop_af_rare_or_unknown(hit: Dict[str, Any]) -> bool:
    g = _gnomad_max(hit)
    if g is None:
        return True
    return float(g) <= float(HLA_SNV_MAX_POP_AF_FOR_DAMAGE)


def _pathogenic_clin_sig(hit: Dict[str, Any]) -> bool:
    cs = str(hit.get("CLIN_SIG") or "").lower()
    return any(
        x in cs
        for x in (
            "pathogenic",
            "likely_pathogenic",
            "pathogenic/likely_pathogenic",
        )
    )


def _spliceai_max_ds(hit: Dict[str, Any]) -> Optional[float]:
    ag = _safe_float(hit.get("SpliceAI_pred_DS_AG"))
    al = _safe_float(hit.get("SpliceAI_pred_DS_AL"))
    xs = [x for x in (ag, al) if x is not None]
    return max(xs) if xs else None


def _pathogenicity_scores_present(hit: Dict[str, Any]) -> bool:
    return any(
        x is not None
        for x in (
            _safe_float(hit.get("CADD_PHRED")),
            _safe_float(hit.get("REVEL")),
            _spliceai_max_ds(hit),
        )
    )


def _pathogenicity_support(hit: Dict[str, Any]) -> bool:
    if _pathogenic_clin_sig(hit):
        return True
    c = _safe_float(hit.get("CADD_PHRED"))
    if c is not None and c >= HLA_SNV_CADD_PHRED_MIN:
        return True
    r = _safe_float(hit.get("REVEL"))
    if r is not None and r >= HLA_SNV_REVEL_MIN:
        return True
    sm = _spliceai_max_ds(hit)
    if sm is not None and sm >= HLA_SNV_SPLICEAI_DS_MIN:
        return True
    return False


def _hit_lof_damage(hit: Dict[str, Any]) -> bool:
    return bool(_hit_lof_like(hit) and _pop_af_rare_or_unknown(hit))


def _hit_high_missense_damage(hit: Dict[str, Any]) -> bool:
    if not _hit_high_missense(hit) or not _pop_af_rare_or_unknown(hit):
        return False
    if not _pathogenicity_scores_present(hit):
        return True
    return bool(_pathogenicity_support(hit))


def _hit_clin_path_damage(hit: Dict[str, Any]) -> bool:
    return bool(_pathogenic_clin_sig(hit) and _pop_af_rare_or_unknown(hit))


def summarize_hla_snv_from_df(df: pd.DataFrame) -> Dict[str, object]:
    """
    Scan rows for HLA transcript hits and return scalar QC fields.

    Expects columns at least: ``gene_hits``, ``tumor_vaf`` (optional).
    """
    if df is None or df.empty or "gene_hits" not in df.columns:
        return {
            "hla_snv_n_rows": 0,
            "hla_snv_lof_any": False,
            "hla_snv_n_lof_variants": 0,
            "hla_snv_lof_rare_any": False,
            "hla_snv_n_lof_rare_variants": 0,
            "hla_snv_high_missense_any": False,
            "hla_snv_high_missense_damage_any": False,
            "hla_snv_n_high_missense_damage_variants": 0,
            "hla_snv_clin_path_any": False,
            "hla_snv_clin_path_damage_any": False,
            "hla_snv_n_clin_path_damage_variants": 0,
            "hla_snv_max_tumor_vaf": np.nan,
            "hla_snv_min_pop_af_seen": np.nan,
            "hla_snv_genes_hit": "",
        }

    n_lof_variants = 0
    n_lof_rare_variants = 0
    n_high_missense_variants = 0
    n_high_missense_damage_variants = 0
    n_clin_path_variants = 0
    n_clin_path_damage_variants = 0
    genes_hit: Set[str] = set()
    tvafs: List[float] = []
    pop_afs: List[float] = []

    for _, row in df.iterrows():
        hits = _parse_gene_hits_cell(row.get("gene_hits"))
        hla_hits = [
            h
            for h in hits
            if str(h.get("Feature_type") or "") == "Transcript" and str(h.get("SYMBOL") or "") in HLA_PANEL_GENES
        ]
        if not hla_hits:
            continue
        row_lof = any(_hit_lof_like(h) for h in hla_hits)
        row_lof_rare = any(_hit_lof_damage(h) for h in hla_hits)
        row_miss = any(_hit_high_missense(h) for h in hla_hits)
        row_miss_dmg = any(_hit_high_missense_damage(h) for h in hla_hits)
        row_clin = any(_pathogenic_clin_sig(h) for h in hla_hits)
        row_clin_dmg = any(_hit_clin_path_damage(h) for h in hla_hits)
        if row_lof:
            n_lof_variants += 1
        if row_lof_rare:
            n_lof_rare_variants += 1
        if row_miss:
            n_high_missense_variants += 1
        if row_miss_dmg:
            n_high_missense_damage_variants += 1
        if row_clin:
            n_clin_path_variants += 1
        if row_clin_dmg:
            n_clin_path_damage_variants += 1
        for h in hla_hits:
            sym = str(h.get("SYMBOL") or "")
            if sym in HLA_PANEL_GENES:
                genes_hit.add(sym)
            g = _gnomad_max(h)
            if g is not None:
                pop_afs.append(float(g))
        tv = _safe_float(row.get("tumor_vaf"))
        if tv is not None:
            tvafs.append(float(tv))

    return {
        "hla_snv_n_rows": int(len(df)),
        "hla_snv_lof_any": bool(n_lof_variants > 0),
        "hla_snv_n_lof_variants": int(n_lof_variants),
        "hla_snv_lof_rare_any": bool(n_lof_rare_variants > 0),
        "hla_snv_n_lof_rare_variants": int(n_lof_rare_variants),
        "hla_snv_high_missense_any": bool(n_high_missense_variants > 0),
        "hla_snv_high_missense_damage_any": bool(n_high_missense_damage_variants > 0),
        "hla_snv_n_high_missense_damage_variants": int(n_high_missense_damage_variants),
        "hla_snv_clin_path_any": bool(n_clin_path_variants > 0),
        "hla_snv_clin_path_damage_any": bool(n_clin_path_damage_variants > 0),
        "hla_snv_n_clin_path_damage_variants": int(n_clin_path_damage_variants),
        "hla_snv_max_tumor_vaf": float(max(tvafs)) if tvafs else float("nan"),
        "hla_snv_min_pop_af_seen": float(min(pop_afs)) if pop_afs else float("nan"),
        "hla_snv_genes_hit": ",".join(sorted(genes_hit)),
    }


def load_hla_snv_summary_for_vial(
    *,
    tumor_vial: str,
    vial_to_snv_csv: Dict[str, Path],
) -> Dict[str, object]:
    p = vial_to_snv_csv.get(str(tumor_vial).strip())
    if p is None or not p.is_file():
        return {
            "hla_snv_file_present": False,
            "hla_snv_n_rows": 0,
            "hla_snv_lof_any": False,
            "hla_snv_n_lof_variants": 0,
            "hla_snv_lof_rare_any": False,
            "hla_snv_n_lof_rare_variants": 0,
            "hla_snv_high_missense_any": False,
            "hla_snv_high_missense_damage_any": False,
            "hla_snv_n_high_missense_damage_variants": 0,
            "hla_snv_clin_path_any": False,
            "hla_snv_clin_path_damage_any": False,
            "hla_snv_n_clin_path_damage_variants": 0,
            "hla_snv_max_tumor_vaf": np.nan,
            "hla_snv_min_pop_af_seen": np.nan,
            "hla_snv_genes_hit": "",
        }
    df = pd.read_csv(p, sep=",", low_memory=False)
    summ = summarize_hla_snv_from_df(df)
    summ["hla_snv_file_present"] = True
    summ["hla_snv_path"] = str(p)
    return summ


def build_default_vial_to_snv_map(snv_output_dir: Path) -> Dict[str, Path]:
    return build_tumor_vial_to_snv_per_sample_csv(
        snv_output_dir=snv_output_dir,
        snv_manifest_path=PATHS.annotations_dir / "SNV" / "samples.tsv",
    )
