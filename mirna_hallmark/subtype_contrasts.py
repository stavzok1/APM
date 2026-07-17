"""Subtype-contrast deep-dive: who is pressured, by which miRNAs, and the immune axis.

Adds the finer-grained, subtype-resolved views requested on top of the cohort
`hallmark_interaction` results, with **Normal-like set aside** (small n +
normal-tissue admixture) and HER2+Basal grouped as *non-luminal* vs LumA+LumB as
*luminal* (plus one-vs-rest singletons):

1. **Grouped anti-correlation** -- pressure↔target-expression Spearman within
   luminal / non-luminal (and per subtype), with Δrho = nonluminal - luminal.
2. **Differential gene pressure** -- per target gene, mean AGO-gated pressure in
   non-luminal vs luminal (+ per subtype); which genes are *differently pressured*.
3. **Per-Hallmark miRNA impact** -- decomposition of each Hallmark's pressure
   into per-miRNA structural targeting load (Σ log1p(evidence) over member-gene
   targets) and the *realized* pressure each miRNA exerts within each group
   (load × group-mean miRNA z). Reveals which miRNA is most impactful per
   Hallmark and whether the dominant regulator differs by subtype.
4. **Immune-axis embodiment (cold-Basal / hot-Luminal)** -- for immune Hallmarks,
   gated pressure / target expression / regulatory-miRNA expression by Basal vs
   Luminal, plus an IRF1 focus matching APM D76 (miR-23a-3p, miR-106b-5p route).

Outputs (``output/subtype_contrasts/``):
- ``hallmark_anticorr_by_group.tsv`` / ``hallmark_anticorr_luminal_vs_nonluminal.tsv``
- ``differential_gene_pressure.tsv``
- ``hallmark_mirna_impact.tsv``
- ``immune_axis_embodiment.tsv`` / ``irf1_route_focus.tsv``
- ``method_manifest.json``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

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


def _subtype_membership(clinical: pd.DataFrame) -> Dict[str, List[str]]:
    """Map group/subtype label -> participant list (Normal-like excluded)."""
    pam = clinical.set_index("participant")["PAM50_final"]
    members: Dict[str, List[str]] = {}
    for grp, subs in C.SUBTYPE_GROUPS.items():
        members[grp] = pam.index[pam.isin(subs)].tolist()
    for sub in C.SUBTYPE_SINGLETONS:
        members[sub] = pam.index[pam == sub].tolist()
    return members


def _anticorr_group(hp, he, clinical, samples, label):
    clin = clinical.set_index("participant")
    conf = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    rows = []
    cols = hp.columns.intersection(samples).intersection(he.columns)
    if len(cols) < 20:
        return pd.DataFrame()
    for hset in hp.index.intersection(he.index):
        df = pd.concat([hp.loc[hset, cols].rename("p"), he.loc[hset, cols].rename("e")], axis=1).dropna()
        if len(df) < 20:
            continue
        cov = D.augment_tcga_batch(clin.reindex(df.index)[conf]) if conf else None
        st = S.correlation_pair(df["e"], df["p"], cov)
        rows.append({
            "group": label,
            "hallmark_set": hset,
            "n": st["n"],
            "spearman_rho": round(float(st["spearman_rho"]), 4),
            "spearman_p": st["spearman_p"],
            "pearson_r": round(float(st["pearson_r"]), 4),
            "pearson_p": st["pearson_p"],
            "partial_rho": (
                round(float(st["partial_rho"]), 4) if np.isfinite(st["partial_rho"]) else np.nan
            ),
            "partial_pearson_r": (
                round(float(st["partial_pearson_r"]), 4)
                if np.isfinite(st["partial_pearson_r"])
                else np.nan
            ),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["spearman_q"] = S.bh_fdr(out["spearman_p"].values)
        out["pearson_q"] = S.bh_fdr(out["pearson_p"].values)
    return out


def grouped_anticorrelation(hp_gated, he, clinical, members) -> pd.DataFrame:
    parts = []
    for label, samples in members.items():
        parts.append(_anticorr_group(hp_gated, he, clinical, samples, label))
    out = pd.concat([p for p in parts if not p.empty], ignore_index=True)
    return out


def differential_gene_pressure(gp_gated, hs, members) -> pd.DataFrame:
    """Per gene: mean gated pressure per group + nonluminal-luminal delta.

    Raw pressure magnitude scales with a gene's miRNA in-degree, so we also report
    a **per-gene z-scored** pressure (standardized across samples). The
    ``zdelta_nonluminal_minus_luminal`` column is in SD units and is the
    cross-gene-comparable effect size; raw means are kept for transparency.
    """
    g2sets = {g: ";".join(s) for g, s in hs.gene_to_sets.items()}
    gpz = S.zscore_rows(gp_gated)  # per-gene z across samples
    raw, z = {}, {}
    for label, samples in members.items():
        s = [c for c in samples if c in gp_gated.columns]
        raw[label] = gp_gated[s].mean(axis=1) if s else pd.Series(np.nan, index=gp_gated.index)
        z[label] = gpz[s].mean(axis=1) if s else pd.Series(np.nan, index=gpz.index)
    df = pd.DataFrame({f"pressure_{k}": v for k, v in raw.items()})
    for k, v in z.items():
        df[f"zpressure_{k}"] = v.round(4)
    df["delta_nonluminal_minus_luminal"] = (df["pressure_nonluminal"] - df["pressure_luminal"]).round(4)
    df["zdelta_nonluminal_minus_luminal"] = (df["zpressure_nonluminal"] - df["zpressure_luminal"]).round(4)
    df.insert(0, "hallmark_sets", df.index.map(g2sets))
    df = df.rename_axis("gene").reset_index()
    return df.sort_values("zdelta_nonluminal_minus_luminal")


def hallmark_mirna_impact(edges, hs, members, *, high_evidence_only: bool = True) -> pd.DataFrame:
    """Per (hallmark, miRNA): structural targeting load + realized pressure by group.

    structural_weight = Σ over member-gene targets of log1p(evidence).
    realized_pressure(group) = structural_weight × mean_z(miRNA within group),
      i.e. how much pressure the miRNA actually exerts on the program in that group.
    """
    e = edges.loc[edges["high_evidence"]] if high_evidence_only else edges
    e = e.copy()
    e["w"] = np.log1p(pd.to_numeric(e["evidence_score"], errors="coerce").fillna(0.0))
    agg = e.groupby(["hallmark_set", "miRNA"]).agg(
        structural_weight=("w", "sum"), n_targets_in_hallmark=("gene", "nunique")
    ).reset_index()

    # group-mean miRNA z (cohort z first, then mean within each group's samples)
    z = S.zscore_rows(D.load_mirna_arms()).fillna(0.0)
    agg = agg.loc[agg["miRNA"].isin(z.index)].copy()
    for label, samples in members.items():
        s = [c for c in samples if c in z.columns]
        gmean = z[s].mean(axis=1) if s else pd.Series(0.0, index=z.index)
        agg[f"mean_z_{label}"] = agg["miRNA"].map(gmean)
        agg[f"realized_pressure_{label}"] = (agg["structural_weight"] * agg[f"mean_z_{label}"]).round(4)
    agg["structural_weight"] = agg["structural_weight"].round(4)
    if "realized_pressure_luminal" in agg and "realized_pressure_nonluminal" in agg:
        agg["realized_delta_nonlum_minus_lum"] = (
            agg["realized_pressure_nonluminal"] - agg["realized_pressure_luminal"]
        ).round(4)
    return agg.sort_values(["hallmark_set", "structural_weight"], ascending=[True, False])


def immune_axis(hp_gated, he, edges, clinical, members, *, mirna_expr) -> pd.DataFrame:
    """Immune-Hallmark gated pressure / target expr / regulatory-miRNA expr by group."""
    # regulatory-miRNA mean expression per immune hallmark per group
    hev = edges.loc[edges["high_evidence"], ["miRNA", "hallmark_set"]].drop_duplicates()
    rows = []
    for hm in C.IMMUNE_HALLMARKS:
        if hm not in hp_gated.index:
            continue
        mirs = hev.loc[hev["hallmark_set"] == hm, "miRNA"]
        mirs = [m for m in mirs if m in mirna_expr.index]
        for label, samples in members.items():
            cols = [c for c in samples if c in hp_gated.columns]
            if len(cols) < 20:
                continue
            ecols = [c for c in samples if c in he.columns]
            mcols = [c for c in samples if c in mirna_expr.columns]
            rows.append({
                "hallmark_set": hm, "group": label, "n": len(cols),
                "mean_gated_pressure": round(float(hp_gated.loc[hm, cols].mean()), 4),
                "mean_target_log2tpm": round(float(he.loc[hm, ecols].mean()), 4) if hm in he.index else np.nan,
                "mean_reg_mirna_log2rpm": round(float(mirna_expr.loc[mirs, mcols].mean().mean()), 4) if mirs else np.nan,
            })
    return pd.DataFrame(rows)


def irf1_focus(gp_gated, edges, clinical, members, *, mirna_expr) -> pd.DataFrame:
    """IRF1-specific pressure by group + the APM D76 miRNA route presence."""
    rows = []
    have = "IRF1" in gp_gated.index
    route_present = {
        m: bool(((edges["miRNA"] == m) & (edges["gene"] == "IRF1")).any()) for m in C.IRF1_ROUTE_MIRNAS
    }
    route_he = {
        m: bool(((edges["miRNA"] == m) & (edges["gene"] == "IRF1") & (edges["high_evidence"])).any())
        for m in C.IRF1_ROUTE_MIRNAS
    }
    for label, samples in members.items():
        cols = [c for c in samples if c in gp_gated.columns]
        mcols = [c for c in samples if c in mirna_expr.columns]
        row = {
            "group": label, "n": len(cols),
            "irf1_mean_gated_pressure": (round(float(gp_gated.loc["IRF1", cols].mean()), 4)
                                         if have and cols else np.nan),
        }
        for m in C.IRF1_ROUTE_MIRNAS:
            row[f"{m}_mean_log2rpm"] = (round(float(mirna_expr.loc[m, mcols].mean()), 4)
                                        if m in mirna_expr.index and mcols else np.nan)
        rows.append(row)
    df = pd.DataFrame(rows)
    df.attrs["route_present"] = route_present
    df.attrs["route_high_evidence"] = route_he
    return df


def run(*, out_dir: Path = C.CONTRASTS_DIR, include_tnrc6: Optional[bool] = None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    edges = D.load_hallmark_edges()
    members = _subtype_membership(clinical)
    print("[contrasts] group sizes: " + ", ".join(f"{k}={len(v)}" for k, v in members.items()))

    print("[contrasts] computing gated gene pressure ...")
    gp = compute_gene_pressure(hs.universe)
    rna = D.load_rna()
    from dataclasses import replace
    params = C.AGO_GATE if include_tnrc6 is None else replace(C.AGO_GATE, include_tnrc6=include_tnrc6)
    gate = compute_ago_gate(rna, params=params)["ago_gate"]
    shared = gp.columns.intersection(gate.index)
    gp_gated = gp[shared].mul(gate.reindex(shared), axis=1)

    hp_gated = hallmark_pressure_matrix(gp_gated, hs)
    he = hallmark_expression_matrix(rna, hs)
    mirna_expr = D.load_mirna_arms()

    # 1) grouped anti-correlation + luminal-vs-nonluminal delta
    ac = grouped_anticorrelation(hp_gated, he, clinical, members)
    ac.to_csv(out_dir / "hallmark_anticorr_by_group.tsv", sep="\t", index=False)
    piv = ac.pivot_table(index="hallmark_set", columns="group", values="spearman_rho")
    if {"luminal", "nonluminal"}.issubset(piv.columns):
        piv["delta_nonluminal_minus_luminal"] = piv["nonluminal"] - piv["luminal"]
        piv.sort_values("delta_nonluminal_minus_luminal").to_csv(
            out_dir / "hallmark_anticorr_luminal_vs_nonluminal.tsv", sep="\t"
        )

    # 2) differential gene pressure
    dgp = differential_gene_pressure(gp_gated, hs, members)
    dgp.to_csv(out_dir / "differential_gene_pressure.tsv", sep="\t", index=False)

    # 3) per-Hallmark miRNA impact
    imp = hallmark_mirna_impact(edges, hs, members)
    imp.to_csv(out_dir / "hallmark_mirna_impact.tsv", sep="\t", index=False)

    # 4) immune axis + IRF1 route
    iax = immune_axis(hp_gated, he, edges, clinical, members, mirna_expr=mirna_expr)
    iax.to_csv(out_dir / "immune_axis_embodiment.tsv", sep="\t", index=False)
    irf1 = irf1_focus(gp_gated, edges, clinical, members, mirna_expr=mirna_expr)
    irf1.to_csv(out_dir / "irf1_route_focus.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.subtype_contrasts",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "groups": {k: list(v) for k, v in C.SUBTYPE_GROUPS.items()},
        "singletons": list(C.SUBTYPE_SINGLETONS),
        "excluded": list(C.SUBTYPE_EXCLUDE_FROM_CONTRAST),
        "group_sizes": {k: len(v) for k, v in members.items()},
        "miRNA_impact": "structural_weight=Σ log1p(evidence) over member-gene targets; "
                        "realized_pressure(group)=structural_weight × group-mean miRNA z",
        "irf1_route": {"mirnas": list(C.IRF1_ROUTE_MIRNAS),
                       "present_any_evidence": irf1.attrs.get("route_present"),
                       "present_high_evidence": irf1.attrs.get("route_high_evidence")},
        "ago_gate_include_tnrc6": params.include_tnrc6,
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    print(f"[contrasts] wrote outputs under {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--include-tnrc6", action="store_true", default=None)
    ap.add_argument("--out-dir", type=Path, default=C.CONTRASTS_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, include_tnrc6=args.include_tnrc6)


if __name__ == "__main__":
    main()
