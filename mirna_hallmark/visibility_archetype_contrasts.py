"""Visibility-archetype contrasts: actual cold-Basal / hot-Luminal (APM D30/D76).

Uses the parent APM visibility-archetype assignments
(``analysis/output/visibility_archetypes/archetype_assignments.tsv``) with the
same group logic as ``analysis/mechanisms/synthesis_second_pass._load_archetype_panel``:

- **cold_Basal**     = PAM50 Basal + cold_silenced archetype
- **hot_Basal**      = PAM50 Basal + hot_IFN_MHCI
- **other_Basal**    = remaining Basal
- **hot_Luminal**    = LumA/LumB + hot_IFN_MHCI
- **other_Luminal**  = remaining LumA/LumB

Normal / Her2 / intermediate archetypes fall in ``other`` and are not primary
contrast arms (same rationale as setting Normal-like aside in PAM50 contrasts).

Unlike ``subtype_contrasts`` (bulk luminal vs non-luminal), this module tests
the **translational archetypes** from APM synthesis.

Outputs (``output/visibility_archetype_contrasts/``):
- ``hallmark_embodiment_by_archetype.tsv``  -- all 50 Hallmarks × group
- ``hallmark_contrast_coldBasal_vs_hotLuminal.tsv`` -- primary contrast deltas
- ``gene_route_focus.tsv``                    -- per-gene IRF1-style panel (immune +
                                              tumor-suppressor/apoptosis nodes)
- ``gene_route_top_mirnas.tsv``               -- top high-evidence driver miRNAs
                                              per gene × group (realized load)
- ``differential_gene_pressure_archetype.tsv`` -- per-gene pressure by archetype group
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.hallmark_interaction import (
    compute_gene_pressure,
    hallmark_expression_matrix,
    hallmark_pressure_matrix,
)

# Parent APM visibility archetype table (not duplicated here).
ARCHETYPE_ASSIGNMENTS = C.REPO_ROOT / "analysis/output/visibility_archetypes/archetype_assignments.tsv"

# Primary contrast groups (ordered for tables).
ARCHETYPE_GROUPS = (
    "cold_Basal",
    "hot_Basal",
    "other_Basal",
    "hot_Luminal",
    "other_Luminal",
)

# Per-gene route panel: immune cascade (APM cold-Basal) + tumor-suppressor/apoptosis
# nodes flagged in bulk basal differential-pressure analysis.
GENE_ROUTE_PANEL: Dict[str, Sequence[str]] = {
    "IRF1": ("hsa-miR-23a-3p", "hsa-miR-106b-5p"),  # APM D76
    "PTEN": ("hsa-miR-21-5p", "hsa-miR-221-3p", "hsa-miR-222-3p"),
    "FOXO1": ("hsa-miR-27a-3p", "hsa-miR-182-5p", "hsa-miR-183-5p"),
    "CDKN1A": ("hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106b-5p", "hsa-miR-93-5p"),
    "BCL2L11": ("hsa-miR-17-5p", "hsa-miR-20a-5p", "hsa-miR-106a-5p", "hsa-miR-106b-5p"),
    "TGFBR2": ("hsa-miR-21-5p", "hsa-miR-93-5p", "hsa-miR-106b-5p"),
    "TP53": ("hsa-miR-125b-5p", "hsa-miR-504-5p", "hsa-miR-34a-5p"),
    "STAT1": ("hsa-miR-145-5p", "hsa-miR-221-3p"),
    "NLRC5": ("hsa-miR-145-5p", "hsa-miR-34a-5p"),
    "CIITA": ("hsa-miR-145-5p", "hsa-miR-155-5p"),
}


def load_archetype_panel() -> pd.DataFrame:
    """Participant table with ``group`` column (APM visibility-archetype logic)."""
    if not ARCHETYPE_ASSIGNMENTS.exists():
        raise FileNotFoundError(
            f"Missing visibility archetype assignments: {ARCHETYPE_ASSIGNMENTS}. "
            "Run: .venv/bin/python3 -m analysis.synthesis.visibility_archetypes"
        )
    assignments = pd.read_csv(ARCHETYPE_ASSIGNMENTS, sep="\t")
    clinical = D.load_clinical_strata()
    df = assignments.merge(clinical, on="participant", how="left", suffixes=("", "_clin"))
    if "PAM50_final" not in df.columns and "PAM50_final_clin" in df.columns:
        df["PAM50_final"] = df["PAM50_final_clin"]
    cold = {a for a in df["archetype"].dropna().unique() if str(a).startswith("cold")}
    hot = {a for a in df["archetype"].dropna().unique() if str(a).startswith("hot")}
    df["group"] = np.select(
        [
            (df["PAM50_final"].eq("Basal")) & df["archetype"].isin(cold),
            (df["PAM50_final"].eq("Basal")) & df["archetype"].isin(hot),
            df["PAM50_final"].eq("Basal"),
            (df["PAM50_final"].isin(["LumA", "LumB"])) & df["archetype"].isin(hot),
            df["PAM50_final"].isin(["LumA", "LumB"]),
        ],
        list(ARCHETYPE_GROUPS),
        default="other",
    )
    return df.drop_duplicates("participant")


def _group_members(panel: pd.DataFrame) -> Dict[str, List[str]]:
    return {g: panel.loc[panel["group"] == g, "participant"].tolist() for g in ARCHETYPE_GROUPS}


def hallmark_embodiment_all(
    hp_gated: pd.DataFrame,
    he: pd.DataFrame,
    edges: pd.DataFrame,
    members: Dict[str, List[str]],
    *,
    mirna_expr: pd.DataFrame,
) -> pd.DataFrame:
    """All Hallmarks × archetype group: gated pressure, target expr, reg-miRNA expr."""
    hev = edges.loc[edges["high_evidence"], ["miRNA", "hallmark_set"]].drop_duplicates()
    rows = []
    for hset in hp_gated.index:
        mirs = hev.loc[hev["hallmark_set"] == hset, "miRNA"]
        mirs = [m for m in mirs if m in mirna_expr.index]
        for label, samples in members.items():
            cols = [c for c in samples if c in hp_gated.columns]
            if len(cols) < 10:
                continue
            ecols = [c for c in samples if c in he.columns]
            mcols = [c for c in samples if c in mirna_expr.columns]
            pvals = hp_gated.loc[hset, cols].to_numpy(dtype=float) if hset in hp_gated.index else np.array([])
            evals = he.loc[hset, ecols].to_numpy(dtype=float) if hset in he.index and ecols else np.array([])
            mvals = (
                mirna_expr.loc[mirs, mcols].to_numpy(dtype=float).ravel()
                if mirs and mcols else np.array([])
            )
            rows.append({
                "hallmark_set": hset,
                "group": label,
                "n": len(cols),
                "mean_gated_pressure": round(float(np.nanmean(pvals)), 4) if pvals.size else np.nan,
                "mean_target_log2tpm": round(float(np.nanmean(evals)), 4) if evals.size else np.nan,
                "mean_reg_mirna_log2rpm": round(float(np.nanmean(mvals)), 4) if mvals.size else np.nan,
            })
    return pd.DataFrame(rows)


def hallmark_primary_contrast(emb: pd.DataFrame) -> pd.DataFrame:
    """Pivot key contrasts: cold_Basal vs hot_Luminal (+ vs other_Basal, hot_Basal)."""
    piv_p = emb.pivot_table(index="hallmark_set", columns="group", values="mean_gated_pressure")
    piv_e = emb.pivot_table(index="hallmark_set", columns="group", values="mean_target_log2tpm")
    piv_m = emb.pivot_table(index="hallmark_set", columns="group", values="mean_reg_mirna_log2rpm")
    out = pd.DataFrame(index=piv_p.index)
    for col in ARCHETYPE_GROUPS:
        if col in piv_p.columns:
            out[f"pressure_{col}"] = piv_p[col]
        if col in piv_e.columns:
            out[f"target_expr_{col}"] = piv_e[col]
        if col in piv_m.columns:
            out[f"reg_mirna_{col}"] = piv_m[col]
    if "cold_Basal" in piv_p.columns and "hot_Luminal" in piv_p.columns:
        out["delta_pressure_coldBasal_minus_hotLuminal"] = piv_p["cold_Basal"] - piv_p["hot_Luminal"]
        out["delta_target_coldBasal_minus_hotLuminal"] = piv_e["cold_Basal"] - piv_e["hot_Luminal"]
    if "cold_Basal" in piv_p.columns and "other_Basal" in piv_p.columns:
        out["delta_pressure_coldBasal_minus_otherBasal"] = piv_p["cold_Basal"] - piv_p["other_Basal"]
    if "cold_Basal" in piv_p.columns and "hot_Basal" in piv_p.columns:
        out["delta_pressure_coldBasal_minus_hotBasal"] = piv_p["cold_Basal"] - piv_p["hot_Basal"]
    if "hot_Luminal" in piv_p.columns and "other_Luminal" in piv_p.columns:
        out["delta_pressure_hotLuminal_minus_otherLuminal"] = piv_p["hot_Luminal"] - piv_p["other_Luminal"]
    return out.reset_index().sort_values("delta_pressure_coldBasal_minus_hotLuminal")


def gene_route_focus(
    gp_gated: pd.DataFrame,
    edges: pd.DataFrame,
    members: Dict[str, List[str]],
    *,
    mirna_expr: pd.DataFrame,
    rna: pd.DataFrame,
) -> pd.DataFrame:
    """IRF1-style per-gene panel: gated pressure + route miRNA expr + target RNA by group."""
    rows = []
    for gene, route_mirs in GENE_ROUTE_PANEL.items():
        for label, samples in members.items():
            cols = [c for c in samples if c in gp_gated.columns]
            mcols = [c for c in samples if c in mirna_expr.columns]
            rcols = [c for c in samples if c in rna.columns]
            row = {
                "gene": gene,
                "group": label,
                "n": len(cols),
                "mean_gated_pressure": (
                    round(float(np.nanmean(gp_gated.loc[gene, cols].to_numpy(dtype=float))), 4)
                    if gene in gp_gated.index and cols else np.nan
                ),
                "mean_target_log2tpm": (
                    round(float(np.nanmean(rna.loc[gene, rcols].to_numpy(dtype=float))), 4)
                    if gene in rna.index and rcols else np.nan
                ),
            }
            for m in route_mirs:
                present = bool(((edges["miRNA"] == m) & (edges["gene"] == gene)).any())
                he = bool(
                    ((edges["miRNA"] == m) & (edges["gene"] == gene) & (edges["high_evidence"])).any()
                )
                row[f"{m}_log2rpm"] = (
                    round(float(mirna_expr.loc[m, mcols].mean()), 4)
                    if m in mirna_expr.index and mcols else np.nan
                )
                row[f"{m}_edge_present"] = present
                row[f"{m}_high_evidence"] = he
            rows.append(row)
    return pd.DataFrame(rows)


def gene_route_top_mirnas(
    edges: pd.DataFrame,
    members: Dict[str, List[str]],
    *,
    mirna_expr: pd.DataFrame,
    genes: Sequence[str],
    top_n: int = 5,
) -> pd.DataFrame:
    """Top high-evidence driver miRNAs per gene × group by realized structural load."""
    he = edges.loc[edges["high_evidence"]].copy()
    he["w"] = np.log1p(pd.to_numeric(he["evidence_score"], errors="coerce").fillna(0.0))
    z = S.zscore_rows(mirna_expr).fillna(0.0)
    rows = []
    for gene in genes:
        sub = he.loc[he["gene"] == gene]
        if sub.empty:
            continue
        load = sub.groupby("miRNA")["w"].sum()
        for label, samples in members.items():
            s = [c for c in samples if c in z.columns]
            if not s:
                continue
            gmean = z[s].mean(axis=1)
            for mir, w in load.sort_values(ascending=False).head(top_n).items():
                rows.append({
                    "gene": gene,
                    "group": label,
                    "miRNA": mir,
                    "structural_weight": round(float(w), 4),
                    "mean_z_mirna": round(float(gmean.get(mir, np.nan)), 4),
                    "realized_load": round(float(w * gmean.get(mir, 0.0)), 4),
                })
    return pd.DataFrame(rows)


def differential_gene_pressure_archetype(gp_gated: pd.DataFrame, members: Dict[str, List[str]]) -> pd.DataFrame:
    cols = {}
    for label, samples in members.items():
        s = [c for c in samples if c in gp_gated.columns]
        cols[label] = gp_gated[s].mean(axis=1) if s else pd.Series(np.nan, index=gp_gated.index)
    df = pd.DataFrame(cols)
    if "cold_Basal" in df and "hot_Luminal" in df:
        df["delta_coldBasal_minus_hotLuminal"] = df["cold_Basal"] - df["hot_Luminal"]
    if "cold_Basal" in df and "other_Basal" in df:
        df["delta_coldBasal_minus_otherBasal"] = df["cold_Basal"] - df["other_Basal"]
    return df.rename_axis("gene").reset_index().sort_values(
        "delta_coldBasal_minus_hotLuminal", ascending=False, na_position="last"
    )


def run(*, out_dir: Path = C.VISIBILITY_ARCHETYPE_DIR, include_tnrc6: Optional[bool] = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    panel = load_archetype_panel()
    members = _group_members(panel)
    print("[visibility_archetype] group sizes: " + ", ".join(f"{k}={len(v)}" for k, v in members.items()))

    hs = HallmarkSets.load()
    edges = D.load_hallmark_edges()
    rna = D.load_rna()
    mirna_expr = D.load_mirna_arms()

    gp = compute_gene_pressure(hs.universe)
    from dataclasses import replace
    params = C.AGO_GATE if include_tnrc6 is None else replace(C.AGO_GATE, include_tnrc6=include_tnrc6)
    gate = compute_ago_gate(rna, params=params)["ago_gate"]
    shared = gp.columns.intersection(gate.index)
    gp_gated = gp[shared].mul(gate.reindex(shared), axis=1)
    hp_gated = hallmark_pressure_matrix(gp_gated, hs)
    he = hallmark_expression_matrix(rna, hs)

    emb = hallmark_embodiment_all(hp_gated, he, edges, members, mirna_expr=mirna_expr)
    emb.to_csv(out_dir / "hallmark_embodiment_by_archetype.tsv", sep="\t", index=False)

    contrast = hallmark_primary_contrast(emb)
    contrast.to_csv(out_dir / "hallmark_contrast_coldBasal_vs_hotLuminal.tsv", sep="\t", index=False)

    route = gene_route_focus(gp_gated, edges, members, mirna_expr=mirna_expr, rna=rna)
    route.to_csv(out_dir / "gene_route_focus.tsv", sep="\t", index=False)

    top_m = gene_route_top_mirnas(
        edges, members, mirna_expr=mirna_expr, genes=list(GENE_ROUTE_PANEL.keys()), top_n=5
    )
    top_m.to_csv(out_dir / "gene_route_top_mirnas.tsv", sep="\t", index=False)

    dgp = differential_gene_pressure_archetype(gp_gated, members)
    dgp.to_csv(out_dir / "differential_gene_pressure_archetype.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.visibility_archetype_contrasts",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "archetype_source": str(ARCHETYPE_ASSIGNMENTS),
        "group_logic": "synthesis_second_pass._load_archetype_panel (cold_silenced/hot_IFN_MHCI)",
        "group_sizes": {k: len(v) for k, v in members.items()},
        "gene_route_panel": {g: list(m) for g, m in GENE_ROUTE_PANEL.items()},
        "primary_contrast": "cold_Basal vs hot_Luminal (all 50 Hallmarks + gene panel)",
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[visibility_archetype] wrote outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include-tnrc6", action="store_true", default=None)
    ap.add_argument("--out-dir", type=Path, default=C.VISIBILITY_ARCHETYPE_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, include_tnrc6=args.include_tnrc6)


if __name__ == "__main__":
    main()
