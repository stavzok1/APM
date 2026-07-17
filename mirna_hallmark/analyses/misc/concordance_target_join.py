"""Target-resolved joins for miRNA arms with strong CN→expression concordance."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D


def _participant_series(s: pd.Series) -> pd.Series:
    """One value per participant (median if duplicate index labels)."""
    s = pd.to_numeric(s, errors="coerce")
    if s.index.has_duplicates:
        s = s.groupby(level=0).median()
    return s


def _gene_row(df: pd.DataFrame, gene: str) -> pd.Series:
    """One participant x expression/CN row per gene symbol (median if duplicated)."""
    row = df.loc[gene]
    if isinstance(row, pd.DataFrame):
        row = row.apply(pd.to_numeric, errors="coerce").median(axis=0)
    return row


def _spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    xa = _participant_series(x)
    ya = _participant_series(y)
    df = pd.concat([xa, ya], axis=1, join="inner").dropna()
    if len(df) < 20:
        return np.nan, np.nan, len(df)
    rho, p = spearmanr(df.iloc[:, 0], df.iloc[:, 1])
    return float(rho), float(p), len(df)


def build_concordance_target_join(
    long: pd.DataFrame,
    concordance: pd.DataFrame,
    *,
    mirtar_path: Optional[Path] = None,
    top_n_arms: int = 25,
    fdr_alpha: float = C.FDR_ALPHA,
    min_evidence_score: int = 0,
    min_target_n: int = 20,
) -> pd.DataFrame:
    """
    For top CN↔miRNA-expression concordant arms, join miRTarBase targets and compute
    participant-level Spearman links among miRNA CN, miRNA expr, target CN, target RNA.

    Sign expectations (distinct questions):
    - miRNA CN ↔ miRNA expr: positive concordance is the dosage-propagation read.
    - miRNA CN/expr ↔ target RNA: for functional repression, negative (anti-) correlation
      is the mechanistic prior; positive links often reflect co-CNV, weak regulation in
      bulk, or non-causal coupling — not “expected” repressor support.
    """
    from analysis.cohort_landscapes.cnv.target_gene_cnv_burden import load_mirtar_edges
    from mirna_hallmark import stats as S

    if concordance.empty:
        return pd.DataFrame()

    sig = concordance.loc[
        (concordance["spearman_q"] < fdr_alpha) & (concordance["spearman_rho"] > 0)
    ].sort_values("spearman_rho", ascending=False)
    top_arms = sig.head(top_n_arms)["arm"].astype(str).tolist()
    if not top_arms:
        return pd.DataFrame()

    edges = load_mirtar_edges(
        Path(mirtar_path or C.MIRTARBASE_SUMMARY),
        min_evidence_score=min_evidence_score,
        functional_mti_only=False,
    )
    edges = edges.loc[edges["miRNA_norm"].isin(top_arms)].copy()
    if edges.empty:
        return pd.DataFrame()

    arms = long.loc[long["entity_level"] == "arm"].copy()
    mir_cn = (
        arms.dropna(subset=["copy_number"])
        .groupby(["entity_label", "participant"])["copy_number"]
        .median()
        .unstack(level=0)
    )
    mir_expr = D.load_mirna_arms()
    target_genes = sorted(edges["gene"].astype(str).unique())
    target_cn = D.load_cnv_target_genes(target_genes)
    target_rna = D.load_rna()

    conc_map = concordance.set_index("arm")
    rows: List[dict] = []
    for arm in top_arms:
        arm_edges = edges.loc[edges["miRNA_norm"] == arm]
        if arm not in mir_cn.columns or arm not in mir_expr.index:
            continue
        cn_s = mir_cn[arm]
        expr_s = mir_expr.loc[arm]
        arm_rho = conc_map.loc[arm, "spearman_rho"] if arm in conc_map.index else np.nan
        arm_q = conc_map.loc[arm, "spearman_q"] if arm in conc_map.index else np.nan
        arm_n = conc_map.loc[arm, "n"] if arm in conc_map.index else np.nan

        for gene, erow in arm_edges.groupby("gene"):
            gene = str(gene)
            if gene not in target_cn.index or gene not in target_rna.index:
                continue
            tc = _gene_row(target_cn, gene)
            tr = _gene_row(target_rna, gene)
            ev = float(pd.to_numeric(erow["evidence_score"], errors="coerce").max())
            rho_ce, p_ce, n_ce = _spearman(cn_s, tr)
            rho_ee, p_ee, n_ee = _spearman(expr_s, tr)
            rho_tc, p_tc, n_tc = _spearman(tc, tr)
            rho_ct, p_ct, n_ct = _spearman(cn_s, tc)
            n_use = min(n_ce, n_ee, n_tc)
            if n_use < min_target_n:
                continue
            rows.append(
                {
                    "mirna_arm": arm,
                    "target_gene": gene,
                    "evidence_score": ev,
                    "n_participants": n_use,
                    "arm_spearman_rho": arm_rho,
                    "arm_spearman_q": arm_q,
                    "arm_concordance_n": arm_n,
                    "rho_mirna_cn__target_rna": round(rho_ce, 4) if np.isfinite(rho_ce) else np.nan,
                    "p_mirna_cn__target_rna": p_ce,
                    "rho_mirna_expr__target_rna": round(rho_ee, 4) if np.isfinite(rho_ee) else np.nan,
                    "p_mirna_expr__target_rna": p_ee,
                    "rho_target_cn__target_rna": round(rho_tc, 4) if np.isfinite(rho_tc) else np.nan,
                    "p_target_cn__target_rna": p_tc,
                    "rho_mirna_cn__target_cn": round(rho_ct, 4) if np.isfinite(rho_ct) else np.nan,
                    "p_mirna_cn__target_cn": p_ct,
                    "mean_mirna_cn": round(float(cn_s.mean()), 4),
                    "mean_mirna_expr": round(float(expr_s.mean()), 4),
                    "mean_target_cn": round(float(tc.mean()), 4),
                    "mean_target_rna": round(float(tr.mean()), 4),
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_mirna_cn__target_rna"] = S.bh_fdr(out["p_mirna_cn__target_rna"].values)
    out["q_mirna_expr__target_rna"] = S.bh_fdr(out["p_mirna_expr__target_rna"].values)
    out = out.sort_values(["arm_spearman_rho", "rho_mirna_cn__target_rna"], ascending=[False, False])
    return out
