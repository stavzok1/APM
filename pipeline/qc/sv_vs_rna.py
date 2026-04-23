from __future__ import annotations

import ast
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from pipeline.config import PATHS
from pipeline.qc.base import Finding, Level
from pipeline.qc.context import to_tcga_sample_id
from pipeline.sample_ids import normalize_tcga_id
from pipeline.RNA_exp.signatures import read_tpm_wide_subset, GeneSets


def _to_sample_vial(x: str) -> str:
    tid = normalize_tcga_id(str(x))
    return str(tid.sample_vial or "").strip()


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


def compute_sv_disruption_vs_rna(
    *,
    ctx,
    sample_ids: List[str],
    genes: List[str] | None = None,
    out_dir: Path | None = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """
    Tier-3 check: SV disruption (promoter/CDS) vs RNA expression loss for key genes.
    """
    gs = GeneSets()
    key_genes = genes or list(gs.apm_class_i)

    # Region-specific “hit” definitions for SV↔RNA QC.
    # For promoter hits we optionally require a minimum overlap fraction, when the SV tables
    # carry promoter overlap fields (newer SV outputs). Older tables only have promoter_flag.
    promoter_min_overlap_frac = 0.20

    # Expression: genes x samples (sample_vial columns)
    expr = read_tpm_wide_subset(
        PATHS.rna_expression,
        genes=key_genes,
        sample_cols=sample_ids,
        gene_col=PATHS.rna_gene_col,
    )

    # IFN covariate (prefer Thorsson IFN-gamma Response; fallback to RNA IFNG_mean signature)
    ifn_by_sample: Dict[str, float] = {}
    th_path = Path("annotations") / "Thornsson_immune_table.tsv"
    if th_path.exists():
        th = pd.read_csv(th_path, sep="\t", low_memory=False)
        if "TCGA Participant Barcode" in th.columns and "IFN-gamma Response" in th.columns:
            th["tcga_sample_id"] = th["TCGA Participant Barcode"].map(to_tcga_sample_id)
            m = th[["tcga_sample_id", "IFN-gamma Response"]].dropna().drop_duplicates(subset=["tcga_sample_id"])
            mp = dict(zip(m["tcga_sample_id"].astype(str), pd.to_numeric(m["IFN-gamma Response"], errors="coerce")))
            for sid in sample_ids:
                ifn_by_sample[sid] = float(mp.get(to_tcga_sample_id(sid))) if to_tcga_sample_id(sid) in mp else float("nan")

    if not ifn_by_sample:
        # compute IFNG_mean from RNA signatures
        from pipeline.RNA_exp.signatures import compute_gene_set_scores_from_tpm

        sig = compute_gene_set_scores_from_tpm(PATHS.rna_expression, sample_ids=sample_ids, gene_col=PATHS.rna_gene_col)
        if "IFNG_mean" in sig.columns:
            for sid in sample_ids:
                ifn_by_sample[sid] = float(pd.to_numeric(sig.loc[sid, "IFNG_mean"], errors="coerce"))

    # SV tables: prefer ctx.sv_output_root if populated else PATHS.sv_output_root
    sv_root = Path(ctx.sv_output_root)
    sv_dir = sv_root / "07_final_sv_with_fimo"
    if not sv_dir.exists() or not any(sv_dir.glob("*_strict_sv_set.csv")):
        sv_dir = Path(PATHS.sv_output_root) / "07_final_sv_with_fimo"

    per_sample_rows: List[Dict[str, object]] = []
    per_gene_rows: List[Dict[str, object]] = []
    summary: Dict[str, object] = {"sv_dir": str(sv_dir), "sv_dir_exists": sv_dir.exists()}

    if not sv_dir.exists():
        per_sample = pd.DataFrame()
        per_gene = pd.DataFrame()
        return per_sample, per_gene, summary

    # Build per-sample disrupted gene set
    disrupted: Dict[str, set] = {sid: set() for sid in sample_ids}
    disrupted_promoter: Dict[str, set] = {sid: set() for sid in sample_ids}
    disrupted_cds: Dict[str, set] = {sid: set() for sid in sample_ids}
    disrupted_by_svtype: Dict[str, Dict[str, set]] = {sid: {} for sid in sample_ids}
    disrupted_promoter_by_svtype: Dict[str, Dict[str, set]] = {sid: {} for sid in sample_ids}
    disrupted_cds_by_svtype: Dict[str, Dict[str, set]] = {sid: {} for sid in sample_ids}
    disrupted_utr: Dict[str, set] = {sid: set() for sid in sample_ids}
    disrupted_intron_only: Dict[str, set] = {sid: set() for sid in sample_ids}
    promoter_hits_missing_overlap_fields = 0

    for sid in sample_ids:
        f = sv_dir / f"{_to_sample_vial(sid)}_strict_sv_set.csv"
        if not f.exists():
            continue
        for chunk in pd.read_csv(f, chunksize=2000, low_memory=False):
            if "gene_hits" not in chunk.columns:
                continue
            has_svtype = "SVTYPE" in chunk.columns
            # cheap prefilter by substring to avoid parsing most rows
            hitmask = chunk["gene_hits"].astype(str).apply(
                lambda s: any(g in s for g in key_genes) if s and s != "[]" else False
            )
            if not hitmask.any():
                continue
            cols = ["gene_hits"]
            if has_svtype:
                cols.append("SVTYPE")
            sub = chunk.loc[hitmask, cols]
            for _, vr in sub.iterrows():
                gh = vr.get("gene_hits")
                svt = str(vr.get("SVTYPE") or "").strip().upper() if has_svtype else ""
                svt = svt or "NA"
                for h in _parse_gene_hits_cell(gh):
                    g = str(h.get("gene_name") or "").strip()
                    if g not in key_genes:
                        continue
                    disrupted[sid].add(g)
                    disrupted_by_svtype[sid].setdefault(svt, set()).add(g)

                    promoter_flag = int(h.get("promoter_flag") or 0) == 1
                    if promoter_flag:
                        # Prefer overlap fraction threshold if the SV table includes it.
                        pov = h.get("promoter_overlap_frac")
                        if pov is None:
                            promoter_hits_missing_overlap_fields += 1
                            promoter_ok = True
                        else:
                            try:
                                promoter_ok = float(pov) >= float(promoter_min_overlap_frac)
                            except Exception:
                                promoter_ok = True
                        if promoter_ok:
                            disrupted_promoter[sid].add(g)
                            disrupted_promoter_by_svtype[sid].setdefault(svt, set()).add(g)
                    if int(h.get("cds_flag") or 0) == 1:
                        disrupted_cds[sid].add(g)
                        disrupted_cds_by_svtype[sid].setdefault(svt, set()).add(g)
                    if int(h.get("utr_flag") or 0) == 1:
                        disrupted_utr[sid].add(g)
                    if int(h.get("intron_only_flag") or 0) == 1:
                        disrupted_intron_only[sid].add(g)

    # CNV covariate: per gene per sample cnv_log2 from ASCAT3 gene tables
    cnv_log2: Dict[tuple, float] = {}
    for sid in sample_ids:
        p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls_ascat3.csv"
        if not p.exists():
            p = Path(ctx.cnv_gene_tables_dir) / f"{sid}_cnv_gene_calls.csv"
        if not p.exists():
            continue
        df = pd.read_csv(p, usecols=["gene_name", "copy_number"], low_memory=False)
        df["gene_name"] = df["gene_name"].astype(str)
        sub = df[df["gene_name"].isin(key_genes)].copy()
        cn = pd.to_numeric(sub["copy_number"], errors="coerce")
        with np.errstate(divide="ignore", invalid="ignore"):
            log2r = np.log2(cn / 2.0)
        for g, v in zip(sub["gene_name"].tolist(), log2r.tolist()):
            try:
                cnv_log2[(sid, str(g))] = float(v)
            except Exception:
                continue

    # Per-sample SV burden (counts) + IFN covariate
    for sid in sample_ids:
        per_sample_rows.append(
            {
                "sample_id": sid,
                "n_apm_genes_disrupted_any": int(len(disrupted.get(sid, set()))),
                "n_apm_genes_disrupted_promoter": int(len(disrupted_promoter.get(sid, set()))),
                "n_apm_genes_disrupted_cds": int(len(disrupted_cds.get(sid, set()))),
                "ifn_covariate": float(ifn_by_sample.get(sid, float("nan"))),
            }
        )

    per_sample = pd.DataFrame(per_sample_rows)

    # Per-gene tests: stratification + regression adjusted for CNV + IFN
    min_hit = 5
    min_no = 10
    cnv_neutral_abs = 0.2

    # Enumerate SVTYPE strata observed in this cohort (kept small and readable in reports).
    svtypes_seen = set()
    for sid in sample_ids:
        svtypes_seen.update(disrupted_by_svtype.get(sid, {}).keys())
    svtypes_seen = {s for s in svtypes_seen if s and s != "NA"}
    svtypes_order = [s for s in ("DEL", "DUP", "INV", "BND", "INS") if s in svtypes_seen] + sorted(
        [s for s in svtypes_seen if s not in {"DEL", "DUP", "INV", "BND", "INS"}]
    )
    summary["svtypes_seen"] = svtypes_order
    summary["promoter_hit_min_overlap_frac"] = float(promoter_min_overlap_frac)
    summary["n_promoter_hits_missing_overlap_fields"] = int(promoter_hits_missing_overlap_fields)

    for g in key_genes:
        if g not in expr.index:
            continue
        y = pd.to_numeric(expr.loc[g], errors="coerce")
        # Hit strata:
        # - ANY: any SV gene hit
        # - PROMOTER/CDS/UTR/INTRON_ONLY: region-specific flags (not mutually exclusive)
        # - plus SVTYPE stratification for ANY (informational)
        strata = [
            ("ANY", "ANY", disrupted),
            ("REGION", "PROMOTER", disrupted_promoter),
            ("REGION", "CDS", disrupted_cds),
            ("REGION", "UTR", disrupted_utr),
            ("REGION", "INTRON_ONLY", disrupted_intron_only),
        ] + [
            ("SVTYPE", svt, {sid: disrupted_by_svtype.get(sid, {}).get(svt, set()) for sid in sample_ids})
            for svt in svtypes_order
        ]
        for strat_kind, strat_value, dmap in strata:
            hit = pd.Series([g in dmap.get(sid, set()) for sid in y.index], index=y.index).astype(bool)
            n_hit = int(hit.sum())
            n_no = int((~hit).sum())
            rec: Dict[str, object] = {
                "gene": g,
                "strat_kind": strat_kind,
                "strat_value": strat_value,
                "n_sv_hit": n_hit,
                "n_sv_no_hit": n_no,
            }

            # Unadjusted median difference
            if n_hit >= 1 and n_no >= 1:
                rec["expr_median_sv_hit"] = float(y[hit].median()) if y[hit].notna().any() else np.nan
                rec["expr_median_no_sv"] = float(y[~hit].median()) if y[~hit].notna().any() else np.nan
                if rec["expr_median_sv_hit"] == rec["expr_median_sv_hit"] and rec["expr_median_no_sv"] == rec["expr_median_no_sv"]:
                    rec["delta_median_unadjusted"] = float(rec["expr_median_sv_hit"] - rec["expr_median_no_sv"])

            # CNV-neutral stratification
            cnv_s = pd.Series([cnv_log2.get((sid, g), np.nan) for sid in y.index], index=y.index, dtype="float")
            neutral = cnv_s.abs() <= cnv_neutral_abs
            nh = int((hit & neutral).sum())
            nn = int(((~hit) & neutral).sum())
            rec["n_cnv_neutral_sv_hit"] = nh
            rec["n_cnv_neutral_no_sv"] = nn
            if nh >= min_hit and nn >= min_no:
                rec["delta_median_cnv_neutral"] = float(y[hit & neutral].median() - y[(~hit) & neutral].median())

            # Regression: y ~ intercept + sv_hit + cnv_log2 + ifn
            if n_hit >= min_hit and n_no >= min_no:
                ifn_s = pd.Series([ifn_by_sample.get(sid, np.nan) for sid in y.index], index=y.index, dtype="float")
                X = pd.DataFrame(
                    {
                        "sv_hit": hit.astype(int),
                        "cnv_log2": cnv_s,
                        "ifn": ifn_s,
                    },
                    index=y.index,
                )
                m = y.notna() & X["cnv_log2"].notna() & X["ifn"].notna()
                if int(m.sum()) >= (min_hit + min_no):
                    yv = y[m].astype(float).values
                    Xm = X[m].astype(float)
                    A = np.column_stack([np.ones(len(Xm)), Xm["sv_hit"].values, Xm["cnv_log2"].values, Xm["ifn"].values])
                    coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
                    # basic R2
                    yhat = A @ coef
                    ss_res = float(((yv - yhat) ** 2).sum())
                    ss_tot = float(((yv - yv.mean()) ** 2).sum()) if len(yv) else np.nan
                    r2 = 1.0 - ss_res / ss_tot if ss_tot and ss_tot == ss_tot and ss_tot > 0 else np.nan
                    rec["reg_beta_sv_hit"] = float(coef[1])
                    rec["reg_beta_cnv_log2"] = float(coef[2])
                    rec["reg_beta_ifn"] = float(coef[3])
                    rec["reg_r2"] = float(r2) if r2 == r2 else np.nan

            per_gene_rows.append(rec)

    per_gene = pd.DataFrame(per_gene_rows)

    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        per_sample.to_csv(out_dir / "sv_vs_rna__per_sample.csv", index=False)
        per_gene.to_csv(out_dir / "sv_vs_rna__per_gene.csv", index=False)

    return per_sample, per_gene, summary


def assert_sv_disruption_vs_rna(summary: Dict[str, object]) -> List[Finding]:
    findings: List[Finding] = []
    if not summary.get("sv_dir_exists", False):
        findings.append(
            Finding(
                check="sv_vs_rna.inputs_present",
                level=Level.WARN,
                message="SV final tables not found; SV↔RNA check skipped.",
                context={"sv_dir": summary.get("sv_dir")},
            )
        )
        return findings

    # Cohort-level assertion moved to per-gene outputs (the caller should summarize per_gene).
    return findings

