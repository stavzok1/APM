"""Roadmap Phase 3 (outcome) — does the framework's PRESSURE predict survival, and at which RESOLUTION?

Companion to `outcome_survival` (raw arm expression). Here we screen the framework's pressure construct
against TCGA outcome across resolutions — expression (arm/miRNA/seed-family), AGO gate capacity,
gene-target pressure, hub-gene pressure, program pressure, sample-aggregate pressure, and role-aggregated
pressure (TF / high-evidence-target). Core measurement (`hallmark_interaction`); each feature z-scored,
Cox vs DFI (recurrence, primary) + OS, crude (age+TNM) and composition-adjusted (epi/imm/str/prolif),
BH-FDR within each resolution family. Headline: does pressure beat expression, and where?

Pressure is TCGA-only (needs matched mRNA+miRNA+AGO) → no cross-cohort pressure validation until
METABRIC-full ([[pending-controlled-data-requests]]); only expression resolutions replicate in GSE22216 (MH-60).

Run: ``.venv/bin/python3 -m mirna_hallmark.analyses.outcome.outcome_pressure_survival``
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome

OUT_DIR = C.OUTPUT_ROOT / "outcome_pressure_survival"
HI = C.OUTPUT_ROOT / "hallmark_interaction"
HYB = C.OUTPUT_ROOT / "hybrid_pressure" / "m2"
AGO = C.OUTPUT_ROOT / "ago_gate" / "per_sample_ago_capacity.tsv.gz"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
TF_TSV = C.REPO_ROOT / "annotations" / "humantfs_lambert2018_tf.tsv"
KEY_PROGRAMS = ["HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"]


def _load(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def _screen(mat: pd.DataFrame, fam: str, cl: dict) -> pd.DataFrame:
    """mat: features x samples. Cox each feature vs DFI (crude+adj) + OS (adj)."""
    base = cl["base"]; adj = pd.concat([base, cl["comp"]], axis=1)
    rows = []
    for feat in mat.index:
        x = mat.loc[feat].astype(float)
        if x.std(skipna=True) == 0 or x.notna().sum() < 50:
            continue
        dcb, dpb = _cox_one(cl["dfi_t"], cl["dfi_e"], x, base)
        dca, dpa = _cox_one(cl["dfi_t"], cl["dfi_e"], x, adj)
        oca, opa = _cox_one(cl["os_t"], cl["os_e"], x, adj)
        rows.append({"resolution": fam, "feature": str(feat),
                     "dfi_hr_adj": np.exp(dca), "dfi_p_adj": dpa, "dfi_p_base": dpb,
                     "os_hr_adj": np.exp(oca), "os_p_adj": opa})
    r = pd.DataFrame(rows)
    for col in ("dfi_p_adj", "os_p_adj"):
        ok = r[col].notna()
        r.loc[ok, col.replace("_p_", "_q_")] = bh_fdr(r.loc[ok, col].to_numpy())
    return r


def _feature_matrices() -> dict:
    mats = {}
    mirna = D.load_mirna_arms(); mirna = mirna[~mirna.index.duplicated(keep="first")]
    mats["expr_arm"] = mirna
    base = mirna.groupby(lambda a: re.sub(r"-(3p|5p)$", "", a)).mean()
    mats["expr_miRNA"] = base
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    e["arm"] = e["key"].str.split(r"\|\|").str[0]
    fam_map = e.dropna(subset=["family"]).drop_duplicates("arm").set_index("arm")["family"]
    mats["expr_seedfamily"] = mirna.join(fam_map.rename("fam")).groupby("fam").mean() if "family" in e else None
    mats["ago_gate"] = _load(AGO).T                                   # 11 features x samples
    gene = _load(HYB / "gene_pressure_per_sample.tsv.gz")
    mats["pressure_gene"] = gene
    mats["pressure_hubgene"] = _load(HYB / "hub_gene_pressure_per_sample.tsv.gz")
    prog = _load(HI / "hallmark_pressure_per_sample.tsv.gz")
    mats["pressure_program"] = prog
    agg = pd.DataFrame({"total_pressure": prog.mean(axis=0),
                        "composite_HYP_EMT_E2F_G2M": prog.reindex(KEY_PROGRAMS).dropna(how="all").mean(axis=0)}).T
    mats["pressure_aggregate"] = agg
    # role-aggregated: mean gene pressure over TF / high-evidence-target genes, per sample
    tfset = set(pd.read_csv(TF_TSV, sep="\t").iloc[:, 0].astype(str))
    he_genes = set(e[e["q_by"] < 0.05]["key"].str.split(r"\|\|").str[1])   # high-evidence (BY-sig) target genes
    role = pd.DataFrame({"TF_pressure": gene.loc[gene.index.intersection(tfset)].mean(axis=0),
                         "HEtarget_pressure": gene.loc[gene.index.intersection(he_genes)].mean(axis=0)}).T
    mats["pressure_role"] = role
    return {k: v for k, v in mats.items() if v is not None and len(v)}


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cl = load_tcga_outcome()
    print(f"[pressure-surv] DFI events {int(cl['dfi_e'].sum())} | OS events {int(cl['os_e'].sum())}")
    mats = _feature_matrices()
    res = pd.concat([_screen(m, fam, cl) for fam, m in mats.items()], ignore_index=True)
    res.to_csv(out_dir / "pressure_resolution_survival.tsv", sep="\t", index=False)

    fam_rows = []
    for fam, g in res.groupby("resolution"):
        top = g.reindex(g["dfi_p_adj"].sort_values().index).head(1)
        fam_rows.append({"resolution": fam, "n_features": int(len(g)),
                         "n_DFI_FDR": int((g["dfi_q_adj"] < C.FDR_ALPHA).sum()),
                         "n_OS_FDR": int((g["os_q_adj"] < C.FDR_ALPHA).sum()),
                         "n_DFI_nominal": int((g["dfi_p_adj"] < 0.05).sum()),
                         "best_DFI": f"{top['feature'].iloc[0]} HR={top['dfi_hr_adj'].iloc[0]:.2f} q={top['dfi_q_adj'].iloc[0]:.3f}" if len(top) else None})
    fam = pd.DataFrame(fam_rows).sort_values("n_DFI_FDR", ascending=False)
    summary = {
        "module": "mirna_hallmark.analyses.outcome.outcome_pressure_survival",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint_primary": "TCGA DFI (recurrence), composition-adjusted; OS secondary",
        "by_resolution": fam.to_dict("records"),
        "top_pressure_hits_DFI": res[res["resolution"].str.startswith("pressure")]
            .reindex(res[res["resolution"].str.startswith("pressure")]["dfi_q_adj"].sort_values().index).head(12)
            .apply(lambda r: f"{r['resolution']}:{r['feature']} HR={r['dfi_hr_adj']:.2f} q={r['dfi_q_adj']:.3f}", axis=1).tolist(),
        "expression_vs_pressure": {
            "expr_DFI_FDR": int(((res["resolution"].str.startswith("expr")) & (res["dfi_q_adj"] < C.FDR_ALPHA)).sum()),
            "pressure_DFI_FDR": int(((res["resolution"].str.startswith("pressure")) & (res["dfi_q_adj"] < C.FDR_ALPHA)).sum()),
        },
        "caveats": ["TCGA-only (pressure needs matched mRNA+miRNA+AGO); DFI sparse (~100 events)",
                    "core measurement only — layer/model robustness deferred to hits",
                    "gene resolution = many features; FDR within each resolution family"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print("[pressure-surv] by resolution (n_DFI_FDR / n_OS_FDR / best DFI):")
    for r in summary["by_resolution"]:
        print(f"   {r['resolution']:22s} feats={r['n_features']:5d} DFI_FDR={r['n_DFI_FDR']:3d} OS_FDR={r['n_OS_FDR']:3d}  {r['best_DFI']}")
    print(f"[pressure-surv] expr DFI-FDR {summary['expression_vs_pressure']['expr_DFI_FDR']} vs pressure {summary['expression_vs_pressure']['pressure_DFI_FDR']}")
    print(f"[pressure-surv] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
