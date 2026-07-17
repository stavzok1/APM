"""S0 vs S1 claim audit: downstream diffs + Basal bootstrap on coupling borderlines.

Run:
    .venv/bin/python3 -m mirna_hallmark.spine_claim_audit
    .venv/bin/python3 -m mirna_hallmark.spine_claim_audit --n-bootstrap 500
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.dual_spine_comparison import (
    S1_BASE_SCORER,
    _build_s0_edges,
    _build_s1_edges,
    _edge_weight_table,
    _hub_rank_rho,
)
from mirna_hallmark.encori_edges import enc_depth_lookup
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_engine import compute_gene_pressure
from mirna_hallmark.robustness_checks import (
    HUB_ROUTES,
    KEY_HALLMARKS,
    _hallmark_cn_matrix,
    _hallmark_expr_matrix,
    _proliferation_proxies,
    _zscore_within,
    gene_pressure_route_partial_corr,
    hallmark_coupling_by_subtype,
    hub_route_partial_corr,
)

OUT_DIR = C.OUTPUT_ROOT / "dual_spine_comparison" / "claim_audit"
S1_ALPHA = 0.5
S1_RUN_ID = f"S1__alpha={S1_ALPHA}__base={S1_BASE_SCORER}"
Q_THRESH = 0.10
SURV_P = 0.05


def _attach_high_evidence(edges: pd.DataFrame) -> pd.DataFrame:
    flags = D.load_hallmark_edges()[["miRNA", "gene", "high_evidence"]].drop_duplicates()
    out = edges.merge(flags, on=["miRNA", "gene"], how="left")
    out["high_evidence"] = out["high_evidence"].fillna(False).astype(bool)
    return out


def _gated_pressure(edges: pd.DataFrame, mirna: pd.DataFrame, gate: pd.Series, hs: HallmarkSets):
    gp = compute_gene_pressure(
        edges, mirna, genes=list(hs.universe),
        expr_mode="softmax_z_logrpm", target_norm="evidence_mass", aggregate="sum",
    )
    shared = gp.columns.intersection(gate.index)
    gp_g = gp[shared].mul(gate.reindex(shared), axis=1)
    hp = hallmark_pressure_matrix(gp_g, hs)
    return gp_g, hp


def _pressure_diff(hp0: pd.DataFrame, hp1: pd.DataFrame) -> pd.DataFrame:
    rows = []
    shared_cols = hp0.columns.intersection(hp1.columns)
    for hset in hp0.index.intersection(hp1.index):
        a = pd.to_numeric(hp0.loc[hset, shared_cols], errors="coerce")
        b = pd.to_numeric(hp1.loc[hset, shared_cols], errors="coerce")
        mask = a.notna() & b.notna()
        if mask.sum() < 20:
            rho = float("nan")
        else:
            rho, _ = spearmanr(a[mask], b[mask])
        rows.append({
            "hallmark_set": hset,
            "key_hallmark": hset in KEY_HALLMARKS,
            "spearman_pressure_s0_s1": float(rho) if np.isfinite(rho) else float("nan"),
            "mean_abs_delta": float((b - a).abs().mean()),
        })
    return pd.DataFrame(rows).sort_values("spearman_pressure_s0_s1")


def _hub_arm_rank_diff(s0_edges: pd.DataFrame, s1_edges: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in HUB_ROUTES:
        e0 = _edge_weight_table(s0_edges.loc[s0_edges["gene"] == gene], "evidence_mass")
        e1 = _edge_weight_table(s1_edges.loc[s1_edges["gene"] == gene], "evidence_mass")
        if e0.empty or e1.empty:
            continue
        m = e0.merge(e1, on=["miRNA", "gene"], suffixes=("_s0", "_s1"))
        if len(m) < 2:
            continue
        r0 = m["edge_w_s0"].rank(ascending=False)
        r1 = m["edge_w_s1"].rank(ascending=False)
        rho, _ = spearmanr(r0, r1)
        top_s0 = set(m.nlargest(5, "edge_w_s0")["miRNA"])
        top_s1 = set(m.nlargest(5, "edge_w_s1")["miRNA"])
        rows.append({
            "gene": gene,
            "n_arms": len(m),
            "rank_spearman": float(rho) if np.isfinite(rho) else float("nan"),
            "top5_jaccard": len(top_s0 & top_s1) / len(top_s0 | top_s1) if top_s0 | top_s1 else float("nan"),
            "top5_s0_only": ",".join(sorted(top_s0 - top_s1)),
            "top5_s1_only": ",".join(sorted(top_s1 - top_s0)),
        })
    return pd.DataFrame(rows)


def _survival_flip_table(
    s0_df: pd.DataFrame,
    s1_df: pd.DataFrame,
    *,
    scope: str,
    rho_col: str,
    p_col: str,
    id_cols: List[str],
) -> pd.DataFrame:
    m = s0_df.merge(s1_df, on=id_cols + ["scope", "n"], suffixes=("_s0", "_s1"))
    m = m.loc[m["scope"] == scope]
    m["surv_s0"] = (m[f"{rho_col}_s0"] < 0) & (m[f"{p_col}_s0"] < SURV_P)
    m["surv_s1"] = (m[f"{rho_col}_s1"] < 0) & (m[f"{p_col}_s1"] < SURV_P)
    m["flip"] = m["surv_s0"] != m["surv_s1"]
    return m


def _basal_coupling_row(
    hset: str,
    cols: pd.Index,
    hp: pd.DataFrame,
    he: pd.DataFrame,
    hc: pd.DataFrame,
    clin: pd.DataFrame,
    pro: pd.Series,
) -> Tuple[float, float]:
    from analysis.utils.common.loaders import partial_spearman

    p = pd.to_numeric(hp.loc[hset].reindex(cols), errors="coerce")
    e = pd.to_numeric(he.loc[hset].reindex(cols), errors="coerce")
    cn = pd.to_numeric(hc.loc[hset].reindex(cols), errors="coerce") if hset in hc.index else None
    cpe = pd.to_numeric(clin["CPE"].reindex(cols), errors="coerce")
    hrd = pd.to_numeric(clin["thornsson_hrd_score"].reindex(cols), errors="coerce")
    pro_s = pd.to_numeric(pro.reindex(cols), errors="coerce")
    ez, pz = _zscore_within(e), _zscore_within(p)
    cov = pd.concat([cpe, hrd, pro_s], axis=1)
    if cn is not None:
        cov = pd.concat([cov, cn], axis=1)
    rho, pval, _ = partial_spearman(ez, pz, cov)
    return float(rho), float(pval)


def bootstrap_basal_coupling(
    hp_s0: pd.DataFrame,
    hp_s1: pd.DataFrame,
    rna: pd.DataFrame,
    clinical: pd.DataFrame,
    proxies: Dict[str, pd.Series],
    hs: HallmarkSets,
    cnv: pd.DataFrame,
    *,
    n_boot: int,
    seed: int,
) -> pd.DataFrame:
    clin = clinical.set_index("participant")
    basal = sorted(clin.index[clin["PAM50_final"].eq("Basal")])
    he = _hallmark_expr_matrix(rna, hs)
    hc = _hallmark_cn_matrix(cnv, hs)
    pro = proxies["e2f_g2m"]
    hallmarks = sorted(hp_s0.index.intersection(he.index))
    rng = np.random.default_rng(seed)

    records = []
    for b in range(n_boot):
        cols = pd.Index(rng.choice(basal, size=len(basal), replace=True))
        rows_s0, rows_s1 = [], []
        for hset in hallmarks:
            rho0, p0 = _basal_coupling_row(hset, cols, hp_s0, he, hc, clin, pro)
            rho1, p1 = _basal_coupling_row(hset, cols, hp_s1, he, hc, clin, pro)
            rows_s0.append({"hallmark_set": hset, "rho": rho0, "p": p0})
            rows_s1.append({"hallmark_set": hset, "rho": rho1, "p": p1})
        for label, rows in [("S0", rows_s0), ("S1", rows_s1)]:
            df = pd.DataFrame(rows)
            df["q"] = S.bh_fdr(df["p"].values)
            df["neg_sig"] = (df["rho"] < 0) & (df["q"] < Q_THRESH)
            records.append({
                "bootstrap": b,
                "spec": label,
                "n_neg_sig": int(df["neg_sig"].sum()),
                "n_key_neg_sig": int(df.loc[df.hallmark_set.isin(KEY_HALLMARKS), "neg_sig"].sum()),
                "heme_neg_sig": bool(df.loc[df.hallmark_set == "HALLMARK_HEME_METABOLISM", "neg_sig"].iloc[0]),
            })
        records.append({
            "bootstrap": b,
            "spec": "diff",
            "s0_minus_s1_neg_sig": int(
                pd.DataFrame(rows_s0).assign(
                    q=S.bh_fdr(pd.DataFrame(rows_s0)["p"].values),
                    neg_sig=lambda d: (d["rho"] < 0) & (d["q"] < Q_THRESH),
                )["neg_sig"].sum()
                - pd.DataFrame(rows_s1).assign(
                    q=S.bh_fdr(pd.DataFrame(rows_s1)["p"].values),
                    neg_sig=lambda d: (d["rho"] < 0) & (d["q"] < Q_THRESH),
                )["neg_sig"].sum()
            ),
        })

    boot = pd.DataFrame([r for r in records if r["spec"] != "diff"])
    diff = pd.DataFrame([r for r in records if r["spec"] == "diff"])

    summary_rows = []
    for spec in ("S0", "S1"):
        g = boot.loc[boot.spec == spec]
        summary_rows.append({
            "metric": f"{spec}_basal_neg_sig",
            "mean": float(g["n_neg_sig"].mean()),
            "p5": float(g["n_neg_sig"].quantile(0.05)),
            "p50": float(g["n_neg_sig"].median()),
            "p95": float(g["n_neg_sig"].quantile(0.95)),
            "p_s0_gt_s1": float((diff["s0_minus_s1_neg_sig"] > 0).mean()) if spec == "S0" else float(
                (diff["s0_minus_s1_neg_sig"] < 0).mean()
            ),
        })
        summary_rows.append({
            "metric": f"{spec}_heme_neg_sig_rate",
            "mean": float(g["heme_neg_sig"].mean()),
            "p5": float("nan"),
            "p50": float(g["heme_neg_sig"].median()),
            "p95": float("nan"),
            "p_s0_gt_s1": float("nan"),
        })
    summary = pd.DataFrame(summary_rows)

    # pairwise heme flip rate
    s0h = boot.loc[boot.spec == "S0"].set_index("bootstrap")["heme_neg_sig"]
    s1h = boot.loc[boot.spec == "S1"].set_index("bootstrap")["heme_neg_sig"]
    heme_flip = float((s0h != s1h).mean())

    summary = pd.concat([
        summary,
        pd.DataFrame([{
            "metric": "heme_neg_sig_flip_rate",
            "mean": heme_flip,
            "p5": float("nan"),
            "p50": heme_flip,
            "p95": float("nan"),
            "p_s0_gt_s1": float((s0h & ~s1h).mean()),
        }]),
    ], ignore_index=True)

    return boot, summary


def _registry_checklist(
    coupling_s0: pd.DataFrame,
    coupling_s1: pd.DataFrame,
    gene_s0: pd.DataFrame,
    gene_s1: pd.DataFrame,
    hub_s0: pd.DataFrame,
    hub_s1: pd.DataFrame,
) -> pd.DataFrame:
    def basal_neg_sig(coupling: pd.DataFrame) -> int:
        g = coupling.loc[coupling.subtype == "Basal"]
        return int(((g.rho_prolif_cn_wsd_adj < 0) & (g.q_prolif_cn_wsd_adj < Q_THRESH)).sum())

    def basal_key(coupling: pd.DataFrame) -> int:
        g = coupling.loc[(coupling.subtype == "Basal") & coupling.key_hallmark]
        return int(((g.rho_prolif_cn_wsd_adj < 0) & (g.q_prolif_cn_wsd_adj < Q_THRESH)).sum())

    def gene_surv_basal_cn(df: pd.DataFrame, target: str) -> bool:
        sub = df.loc[(df.target == target) & (df.scope == "Basal")]
        if sub.empty:
            return False
        r = sub.iloc[0]
        return bool(r.get("survives_e2f_g2m_CN", r.get("survives_e2f_g2m", False)))

    def hub_surv_basal_cn(df: pd.DataFrame, target: str, mirna: str) -> bool:
        sub = df.loc[(df.target == target) & (df.miRNA == mirna) & (df.scope == "Basal")]
        if sub.empty:
            return False
        r = sub.iloc[0]
        return bool(r.get("survives_e2f_g2m_CN", r.get("survives_e2f_g2m", False)))

    claims = [
        ("MH-30", "Basal neg-sig Hallmarks", f"{basal_neg_sig(coupling_s0)}/50", f"{basal_neg_sig(coupling_s1)}/50"),
        ("MH-17", "Basal key 8/8 prolif-adj", f"{basal_key(coupling_s0)}/8", f"{basal_key(coupling_s1)}/8"),
        ("MH-22", "CDKN1A gene pressure Basal CN", gene_surv_basal_cn(gene_s0, "CDKN1A"), gene_surv_basal_cn(gene_s1, "CDKN1A")),
        ("MH-22", "TGFBR2 gene pressure Basal CN", gene_surv_basal_cn(gene_s0, "TGFBR2"), gene_surv_basal_cn(gene_s1, "TGFBR2")),
        ("MH-22", "VIM gene pressure Basal CN", gene_surv_basal_cn(gene_s0, "VIM"), gene_surv_basal_cn(gene_s1, "VIM")),
        ("MH-22", "PTEN gene pressure Basal CN", gene_surv_basal_cn(gene_s0, "PTEN"), gene_surv_basal_cn(gene_s1, "PTEN")),
        ("MH-22", "IRF1 gene pressure Basal CN", gene_surv_basal_cn(gene_s0, "IRF1"), gene_surv_basal_cn(gene_s1, "IRF1")),
        ("MH-19", "PTEN miR-106b route Basal CN", hub_surv_basal_cn(hub_s0, "PTEN", "hsa-miR-106b-5p"), hub_surv_basal_cn(hub_s1, "PTEN", "hsa-miR-106b-5p")),
        ("MH-19", "CDKN1A miR-17 route Basal CN", hub_surv_basal_cn(hub_s0, "CDKN1A", "hsa-miR-17-5p"), hub_surv_basal_cn(hub_s1, "CDKN1A", "hsa-miR-17-5p")),
    ]
    rows = []
    for cid, label, v0, v1 in claims:
        rows.append({
            "claim_id": cid,
            "claim": label,
            "S0": v0,
            "S1": v1,
            "unchanged": str(v0) == str(v1),
        })
    return pd.DataFrame(rows)


def run(*, n_bootstrap: int = 500, seed: int = 42) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[claim_audit] bootstrap n={n_bootstrap}")

    hs = HallmarkSets.load()
    clinical = D.load_clinical_strata()
    rna = D.load_rna().groupby(level=0).mean()
    mirna = D.load_mirna_arms()
    proxies = _proliferation_proxies(rna, hs)
    hub_genes = sorted(set(HUB_ROUTES) | {g for gs in hs.sets.values() for g in gs})
    cnv = D.load_cnv_target_genes(hub_genes)
    enc_lookup = enc_depth_lookup(genes=hs.universe)
    gate = compute_ago_gate(
        D.load_rna(),
        params=AgoGateParams(
            include_tnrc6=False,
            gate_min=C.AGO_GATE.gate_min,
            gate_k=C.AGO_GATE.gate_k,
            gate_midpoint=C.AGO_GATE.gate_midpoint,
        ),
    )["ago_gate"]

    s0_edges = _attach_high_evidence(_build_s0_edges(hs.universe, mirna))
    s1_edges = _attach_high_evidence(
        _build_s1_edges(hs.universe, mirna, s0_edges, enc_lookup, alpha=S1_ALPHA, base_scorer=S1_BASE_SCORER)
    )

    gp_s0, hp_s0 = _gated_pressure(s0_edges, mirna, gate, hs)
    gp_s1, hp_s1 = _gated_pressure(s1_edges, mirna, gate, hs)

    # --- downstream diffs ---
    pressure_diff = _pressure_diff(hp_s0, hp_s1)
    pressure_diff.to_csv(OUT_DIR / "hallmark_pressure_correlation.tsv", sep="\t", index=False)

    hub_arms = _hub_arm_rank_diff(s0_edges, s1_edges)
    hub_arms.to_csv(OUT_DIR / "hub_top_arm_rank_stability.tsv", sep="\t", index=False)

    coupling_s0 = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp_s0, cnv)
    coupling_s1 = hallmark_coupling_by_subtype(rna, clinical, proxies, hs, hp_s1, cnv)
    coupling_s0["spec"] = "S0"
    coupling_s1["spec"] = S1_RUN_ID
    pd.concat([coupling_s0, coupling_s1]).to_csv(OUT_DIR / "coupling_recheck.tsv", sep="\t", index=False)

    gene_s0 = gene_pressure_route_partial_corr(rna, clinical, proxies, gp_s0, cnv)
    gene_s1 = gene_pressure_route_partial_corr(rna, clinical, proxies, gp_s1, cnv)
    gene_flip = _survival_flip_table(
        gene_s0, gene_s1, scope="Basal", rho_col="rho_e2f_g2m_CN", p_col="p_e2f_g2m_CN",
        id_cols=["target"],
    )
    gene_flip.to_csv(OUT_DIR / "gene_pressure_survival_flip_basal.tsv", sep="\t", index=False)

    hub_s0 = hub_route_partial_corr(rna, mirna, clinical, proxies, s0_edges, cnv)
    hub_s1 = hub_route_partial_corr(rna, mirna, clinical, proxies, s1_edges, cnv)
    hub_flip = _survival_flip_table(
        hub_s0, hub_s1, scope="Basal", rho_col="rho_e2f_g2m_CN", p_col="p_e2f_g2m_CN",
        id_cols=["target", "miRNA"],
    )
    hub_flip.to_csv(OUT_DIR / "hub_route_survival_flip_basal.tsv", sep="\t", index=False)

    registry = _registry_checklist(coupling_s0, coupling_s1, gene_s0, gene_s1, hub_s0, hub_s1)
    registry.to_csv(OUT_DIR / "registry_claim_checklist.tsv", sep="\t", index=False)

    # --- bootstrap ---
    boot, boot_summary = bootstrap_basal_coupling(
        hp_s0, hp_s1, rna, clinical, proxies, hs, cnv,
        n_boot=n_bootstrap, seed=seed,
    )
    boot.to_csv(OUT_DIR / "bootstrap_basal_coupling.tsv", sep="\t", index=False)
    boot_summary.to_csv(OUT_DIR / "bootstrap_summary.tsv", sep="\t", index=False)

    manifest = {
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "s1_alpha": S1_ALPHA,
        "n_bootstrap": n_bootstrap,
        "pressure_median_spearman": float(pressure_diff["spearman_pressure_s0_s1"].median()),
        "pressure_min_spearman_key8": float(
            pressure_diff.loc[pressure_diff.key_hallmark, "spearman_pressure_s0_s1"].min()
        ),
        "gene_basal_survival_flips": int(gene_flip["flip"].sum()),
        "hub_basal_survival_flips": int(hub_flip["flip"].sum()),
        "registry_unchanged": int(registry["unchanged"].sum()),
        "registry_total": int(len(registry)),
        "bootstrap": boot_summary.set_index("metric")["mean"].to_dict(),
    }
    (OUT_DIR / "audit_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    # README
    lines = [
        "# S0 vs S1 claim audit",
        "",
        f"Generated: {manifest['finished_utc']}",
        "",
        "## Pressure stability",
        f"- Median Spearman(sample pressure) across Hallmarks: **{manifest['pressure_median_spearman']:.4f}**",
        f"- Key-8 minimum: **{manifest['pressure_min_spearman_key8']:.4f}**",
        "",
        "## Registry claims",
        f"- Unchanged: **{manifest['registry_unchanged']}/{manifest['registry_total']}**",
        "",
        "## Gene / hub survival flips (Basal, e2f_g2m+CN)",
        f"- Gene aggregate pressure: **{manifest['gene_basal_survival_flips']}** flips",
        f"- Hub miRNA routes: **{manifest['hub_basal_survival_flips']}** flips",
        "",
        "## Bootstrap (Basal neg-sig, primary metric)",
    ]
    for _, r in boot_summary.iterrows():
        lines.append(f"- `{r['metric']}`: mean={r['mean']:.3f}, median={r['p50']:.3f}")
    (OUT_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[claim_audit] wrote {OUT_DIR}")


def main() -> None:
    ap = argparse.ArgumentParser(description="S0 vs S1 downstream claim audit + bootstrap")
    ap.add_argument("--n-bootstrap", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    run(n_bootstrap=args.n_bootstrap, seed=args.seed)


if __name__ == "__main__":
    main()
