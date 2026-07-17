"""Program-level miRNA-pressure aggregation stratified by gene ROLE (TF / high-evidence target).

Motivation (the marker-vs-target lesson, MH-56): single activity markers and naive program-average
pressure are misleading — SLC2A1's program-average pressure was +0.24 but its gene-level pressure was
−2.47 (a masked brake-release). So aggregate the framework's *own* per-gene pressure
(`gene_pressure_by_state.delta_tumor_gtex`; <0 = released/brake-release, >0 = gained) **per Hallmark
program**, and stratify by gene role using existing annotations:
  - **TF** vs non-TF (Lambert 2018 human-TF census, `annotations/humantfs_lambert2018_tf.tsv`)
  - **high-evidence miRNA-target** vs other (the curated edge table)

Questions answered: (1) which programs are net brake-released vs net pressure-gained, aggregated
properly? (2) Is brake-release CONCENTRATED on the program's regulators (TFs) / well-supported
targets, or diffuse? — i.e. are miRNAs releasing the *drivers* of a program.

Run: ``.venv/bin/python3 -m mirna_hallmark.program_pressure_by_role``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr

GP = C.TISSUE_REFERENCE_DIR / "cross_state_landscape" / "gene_pressure_by_state.tsv"
EDGES = C.OUTPUT_ROOT / "edges" / "mirna_hallmark_edges.tsv.gz"
TF = C.REPO_ROOT / "annotations" / "humantfs_lambert2018_tf.tsv"
OUT_DIR = C.OUTPUT_ROOT / "program_pressure_by_role"


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    gp = pd.read_csv(GP, sep="\t")[["gene", "delta_tumor_gtex"]]
    edges = pd.read_csv(EDGES, sep="\t")
    tfset = set(pd.read_csv(TF, sep="\t").iloc[:, 0].astype(str))
    he_genes = set(edges[edges["high_evidence"]]["gene"])
    gene_hallmarks = edges[["gene", "hallmark_set"]].drop_duplicates()

    g = gp.merge(gene_hallmarks, on="gene", how="inner")
    g["is_tf"] = g["gene"].isin(tfset)
    g["is_he"] = g["gene"].isin(he_genes)
    g["released"] = g["delta_tumor_gtex"] < 0      # brake-release

    # --- per-program aggregation, stratified ---
    rows = []
    for hm, sub in g.groupby("hallmark_set"):
        if len(sub) < 8:
            continue
        tf = sub[sub["is_tf"]]; nontf = sub[~sub["is_tf"]]
        rows.append({
            "hallmark": hm, "n_genes": len(sub), "n_tf": int(len(tf)),
            "median_delta_all": round(float(sub["delta_tumor_gtex"].median()), 3),
            "frac_released_all": round(float(sub["released"].mean()), 3),
            "median_delta_TF": round(float(tf["delta_tumor_gtex"].median()), 3) if len(tf) >= 3 else None,
            "median_delta_nonTF": round(float(nontf["delta_tumor_gtex"].median()), 3),
            "median_delta_HEtarget": round(float(sub[sub["is_he"]]["delta_tumor_gtex"].median()), 3) if sub["is_he"].sum() >= 3 else None,
        })
    prog = pd.DataFrame(rows).sort_values("median_delta_all")
    prog.to_csv(out_dir / "program_pressure_by_role.tsv", sep="\t", index=False)

    # --- global: is brake-release concentrated on TFs? (gene-level, dedup; dropna) ---
    gu = g.drop_duplicates("gene").dropna(subset=["delta_tumor_gtex"])
    tf_d = gu[gu["is_tf"]]["delta_tumor_gtex"]; nontf_d = gu[~gu["is_tf"]]["delta_tumor_gtex"]
    # NB negative delta = released; "less" tests TFs MORE released (more negative) than non-TF
    p_tf_rel = mannwhitneyu(tf_d, nontf_d, alternative="less")[1] if len(tf_d) and len(nontf_d) else np.nan
    p_tf_gain = mannwhitneyu(tf_d, nontf_d, alternative="greater")[1] if len(tf_d) and len(nontf_d) else np.nan
    summary = {
        "module": "mirna_hallmark.program_pressure_by_role",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_gene_hallmark_pairs": int(len(g)), "n_unique_genes": int(len(gu)), "n_programs": int(len(prog)),
        "TF_vs_nonTF": {"n_tf": int(len(tf_d)), "median_TF": round(float(tf_d.median()), 3),
                        "median_nonTF": round(float(nontf_d.median()), 3),
                        "TF_more_released_p": float(f"{p_tf_rel:.3g}") if np.isfinite(p_tf_rel) else None,
                        "TF_more_GAINED_p": float(f"{p_tf_gain:.3g}") if np.isfinite(p_tf_gain) else None},
        "most_brake_released_programs": prog.nsmallest(6, "median_delta_all")[
            ["hallmark", "median_delta_all", "median_delta_TF", "frac_released_all"]].to_dict("records"),
        "most_pressure_gained_programs": prog.nlargest(6, "median_delta_all")[
            ["hallmark", "median_delta_all", "median_delta_TF"]].to_dict("records"),
        "caveats": ["pressure-by-state is bulk TCGA (composition-confounded); delta is healthy(GTEx)->tumor",
                    "TF = Lambert2018 curated; HE = curated high-evidence edge target",
                    "aggregation corrects the marker/program-average pitfall (MH-56)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[prog_role] {len(g)} gene-hallmark pairs, {len(prog)} programs")
    print(f"[prog_role] TF vs nonTF pressure: median {summary['TF_vs_nonTF']['median_TF']} vs {summary['TF_vs_nonTF']['median_nonTF']} "
          f"(TF-more-released p={summary['TF_vs_nonTF']['TF_more_released_p']}; TF-more-GAINED p={summary['TF_vs_nonTF']['TF_more_GAINED_p']})")
    print("[prog_role] most brake-RELEASED programs (median_delta_all, TF):")
    for r in summary["most_brake_released_programs"]:
        print(f"    {r['hallmark']:38s} {r['median_delta_all']:+.2f}  TF={r['median_delta_TF']}  frac_rel={r['frac_released_all']}")
    print("[prog_role] most pressure-GAINED programs:")
    for r in summary["most_pressure_gained_programs"]:
        print(f"    {r['hallmark']:38s} {r['median_delta_all']:+.2f}  TF={r['median_delta_TF']}")
    print(f"[prog_role] wrote {out_dir}")
    return prog


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
