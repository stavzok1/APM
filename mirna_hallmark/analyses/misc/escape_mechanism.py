"""Subtype-resolved chromatin escape context for decoupled miRNA->gene edges.

For lost / target-specific-escape edges, tests whether the **target gene's promoter is
chromatin-open** (ATAC) in tumor — a plausible non-miRNA driver of target upregulation that
would break partial coupling without the miRNA failing.

Because ATAC n≈74 (pipeline full-panel), results are **PAM50-aggregated** (mean promoter ATAC
per gene × subtype), not per-sample coupling. Promoter methylation (hub-gene matrix, where
available) is joined as a parallel epigenetic layer.

Outputs under ``output/tissue_reference/escape_mechanism/``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D

OUT_DIR = C.TISSUE_REFERENCE_DIR / "decoupling_validation"
STATE_DIR = C.TISSUE_REFERENCE_DIR / "mirna_state_class"
DECOUPLE = OUT_DIR / "decoupling_validation.tsv"
ATAC_HALLMARK = C.REPO_ROOT / "analysis/output/chromatin_per_sample_hallmark/per_sample_atac_promoter.tsv.gz"
ATAC_FULLPANEL = C.REPO_ROOT / "analysis/output/chromatin_per_sample_fullpanel/per_sample_atac_promoter.tsv.gz"
ATAC_PATH = ATAC_HALLMARK if ATAC_HALLMARK.exists() else ATAC_FULLPANEL
HALLMARK_METH_PATH = C.OUTPUT_ROOT / "matrices/hallmark_gene_methylation.tsv.gz"
HUB_METH_PATH = C.OUTPUT_ROOT / "matrices/hub_gene_methylation.tsv.gz"
METH_PATH = HALLMARK_METH_PATH if HALLMARK_METH_PATH.exists() else HUB_METH_PATH
PAM50 = ("LumA", "LumB", "Her2", "Basal")


def _meth_long(path: Path = METH_PATH) -> pd.DataFrame:
    """Normalise the per-sample methylation long table to (participant, gene, beta)."""
    if not path.exists():
        return pd.DataFrame()
    m = pd.read_csv(path, sep="\t")
    gene_col = "gene_name" if "gene_name" in m.columns else ("gene" if "gene" in m.columns else None)
    beta_col = next((c for c in ("meth_beta_mean", "promoter_beta", "beta_mean", "beta")
                     if c in m.columns), None)
    if gene_col is None or beta_col is None or "participant" not in m.columns:
        return pd.DataFrame()
    out = m[["participant", gene_col, beta_col]].rename(columns={gene_col: "gene", beta_col: "beta"})
    out["beta"] = pd.to_numeric(out["beta"], errors="coerce")
    return out.dropna(subset=["beta"])


def load_meth_per_sample(genes: Sequence[str], path: Path = METH_PATH) -> pd.DataFrame:
    """gene x participant promoter/body methylation beta (per-sample, Hallmark genes)."""
    long = _meth_long(path)
    if long.empty:
        return pd.DataFrame()
    long = long[long["gene"].isin(set(genes))]
    return long.pivot_table(index="gene", columns="participant", values="beta", aggfunc="median")


def load_atac_per_sample(genes: Sequence[str], *, path: Path = ATAC_PATH,
                         signal_col: str = "mean_atac_signal_promoter") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-sample ATAC promoter signal with **subtype-representative fallback**.

    Returns (value_matrix, imputed_mask) as gene x participant over ALL PAM50-typed tumor
    participants: a participant's value is its DIRECT ATAC where measured (n~74), otherwise the
    mean of same-PAM50 participants that DO have ATAC. `imputed_mask` is True where filled by the
    subtype representative (so downstream can flag that those carry only between-subtype variation).
    """
    if not path.exists():
        return pd.DataFrame(), pd.DataFrame()
    atac = pd.read_csv(path, sep="\t")
    atac = atac[atac["gene"].isin(set(genes))].copy()
    atac["sig"] = pd.to_numeric(atac[signal_col], errors="coerce")
    clin = D.load_clinical_strata().dropna(subset=["PAM50_final"])
    pam = clin.set_index("participant")["PAM50_final"].to_dict()
    atac["PAM50"] = atac["participant"].map(pam)

    direct = atac.pivot_table(index="gene", columns="participant", values="sig", aggfunc="median")
    sub_mean = (atac[atac["PAM50"].isin(PAM50)]
                .groupby(["gene", "PAM50"])["sig"].mean().unstack())  # gene x PAM50

    # all PAM50-typed participants (the universe we want a per-sample value for)
    all_parts = [p for p, s in pam.items() if s in PAM50]
    genes_idx = sub_mean.index
    val = pd.DataFrame(index=genes_idx, columns=all_parts, dtype=float)
    imp = pd.DataFrame(False, index=genes_idx, columns=all_parts)
    direct_parts = set(direct.columns)
    for p in all_parts:
        s = pam[p]
        if p in direct_parts:
            col = direct[p].reindex(genes_idx)
            # backfill any gene the direct sample is missing with the subtype mean
            miss = col.isna()
            if miss.any() and s in sub_mean.columns:
                col[miss] = sub_mean.loc[miss[miss].index, s]
                imp.loc[miss[miss].index, p] = True
            val[p] = col
        else:
            if s in sub_mean.columns:
                val[p] = sub_mean[s].reindex(genes_idx)
                imp[p] = True
    return val, imp


def _load_atac_by_pam50(path: Path = ATAC_PATH) -> pd.DataFrame:
    """(gene, PAM50) mean promoter ATAC signal."""
    if not path.exists():
        return pd.DataFrame()
    atac = pd.read_csv(path, sep="\t")
    clin = D.load_clinical_strata().dropna(subset=["PAM50_final"])
    pam = clin.set_index("participant")["PAM50_final"]
    atac["PAM50"] = atac["participant"].map(pam)
    atac = atac[atac["PAM50"].isin(PAM50)]
    sig = pd.to_numeric(atac["mean_atac_signal_promoter"], errors="coerce")
    atac["atac_promoter"] = sig
    return (atac.groupby(["gene", "PAM50"], as_index=False)
            .agg(n_atac=("participant", "nunique"),
                 mean_atac_promoter=("atac_promoter", "mean"),
                 median_atac_promoter=("atac_promoter", "median")))


def _load_meth_by_pam50(path: Path = METH_PATH) -> pd.DataFrame:
    """(gene, PAM50) mean methylation beta from the per-sample long matrix."""
    long = _meth_long(path)
    if long.empty:
        return pd.DataFrame()
    clin = D.load_clinical_strata().dropna(subset=["PAM50_final"])
    pam = clin.set_index("participant")["PAM50_final"]
    long["PAM50"] = long["participant"].map(pam)
    long = long[long["PAM50"].isin(PAM50)]
    return (long.groupby(["gene", "PAM50"], as_index=False)
            .agg(n_meth=("participant", "nunique"),
                 mean_promoter_beta=("beta", "mean")))


def run(*, decouple_path: Path = DECOUPLE, out_dir: Path = OUT_DIR,
        classes: Sequence[str] = ("lost", "nat_decoupled"),
        mechanisms: Sequence[str] = ("target_specific_escape", "arm_wide_failure", "intermediate")) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    if not decouple_path.exists():
        raise FileNotFoundError(f"Run decoupling_validation first: {decouple_path}")
    dc = pd.read_csv(decouple_path, sep="\t")
    dc = dc[dc["joint_edge_class"].isin(classes)]
    if "decoupling_mechanism" in dc.columns and mechanisms:
        dc = dc[dc["decoupling_mechanism"].isin(mechanisms)]
    print(f"[escape] {len(dc)} decoupled edges for chromatin context")

    atac = _load_atac_by_pam50()
    meth = _load_meth_by_pam50()
    if atac.empty:
        print("[escape] WARN: ATAC matrix missing — skipping")
        return pd.DataFrame()

    # gene-level ATAC reference (cohort median across subtypes) for open/closed call
    gene_med = atac.groupby("gene")["mean_atac_promoter"].median()

    rows: List[dict] = []
    for _, e in dc.iterrows():
        g = e["gene"]
        base = {"miRNA": e["miRNA"], "gene": g,
                "decoupling_mechanism": e.get("decoupling_mechanism"),
                "nat_anchored_confident": e.get("nat_anchored_confident"),
                "regulator_share_rank_pct": e.get("regulator_share_rank_pct"),
                "n_competitors_coupled_tumor": e.get("n_competitors_coupled_tumor")}
        g_med = float(gene_med.get(g, np.nan))
        for sub in PAM50:
            row = dict(base)
            row["PAM50"] = sub
            a = atac.loc[(atac["gene"] == g) & (atac["PAM50"] == sub)]
            if not a.empty:
                val = float(a["mean_atac_promoter"].iloc[0])
                row["mean_atac_promoter"] = val
                row["n_atac"] = int(a["n_atac"].iloc[0])
                row["atac_open_vs_gene_median"] = bool(val > g_med) if np.isfinite(g_med) else np.nan
            else:
                row["mean_atac_promoter"] = row["n_atac"] = np.nan
                row["atac_open_vs_gene_median"] = np.nan
            if not meth.empty:
                m = meth.loc[(meth["gene"] == g) & (meth["PAM50"] == sub)]
                if not m.empty:
                    row["mean_promoter_beta"] = float(m["mean_promoter_beta"].iloc[0])
                    row["n_meth"] = int(m["n_meth"].iloc[0])
                else:
                    row["mean_promoter_beta"] = row["n_meth"] = np.nan
            rows.append(row)

    long = pd.DataFrame(rows)
    long.to_csv(out_dir / "escape_atac_by_pam50.tsv", sep="\t", index=False)

    # edge-level summary: is target ATAC-open in ANY subtype (promoter accessibility escape)?
    if not long.empty:
        agg_spec = {"any_subtype_atac_open": ("atac_open_vs_gene_median", lambda s: bool(s.any())),
                    "max_mean_atac": ("mean_atac_promoter", "max"),
                    "n_subtypes_atac": ("mean_atac_promoter", lambda s: int(s.notna().sum()))}
        if "mean_promoter_beta" in long.columns:
            agg_spec["mean_promoter_beta"] = ("mean_promoter_beta", "mean")
        edge_sum = long.groupby(["miRNA", "gene"], as_index=False).agg(**agg_spec)
        edge_sum.to_csv(out_dir / "escape_atac_edge_summary.tsv", sep="\t", index=False)
        n_open = int(edge_sum["any_subtype_atac_open"].sum())
        print(f"[escape] target promoter ATAC-open (any subtype): {n_open}/{len(edge_sum)} edges "
              f"({100*n_open/max(len(edge_sum),1):.0f}%)")
        conf = edge_sum.merge(dc[dc.get("nat_anchored_confident", False) == True][["miRNA", "gene"]],
                              on=["miRNA", "gene"], how="inner")
        if not conf.empty:
            print(f"[escape] NAT-confident subset ATAC-open: "
                  f"{int(conf['any_subtype_atac_open'].sum())}/{len(conf)}")

    manifest = {"module": "mirna_hallmark.analyses.misc.escape_mechanism",
                "generated_utc": datetime.now(timezone.utc).isoformat(),
                "n_edges": int(len(dc)), "atac_path": str(ATAC_PATH),
                "meth_path": str(METH_PATH) if METH_PATH.exists() else None,
                "aggregation": "PAM50 mean (ATAC n~74)"}
    (out_dir / "escape_method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return long


def main() -> None:
    ap = argparse.ArgumentParser(description="ATAC/meth escape context for decoupled edges")
    ap.add_argument("--decouple-path", type=Path, default=DECOUPLE)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    run(decouple_path=args.decouple_path, out_dir=args.out_dir)


if __name__ == "__main__":
    main()
