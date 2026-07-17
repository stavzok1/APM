"""Cell-line <-> tumour-subtype fidelity across the framework's measurement layers.

Instead of asking whether the framework's miRNA->target edges anti-correlate across cell lines (the
cross-line repression test, MH-67 -- which the NanoString miRNA + co-regulation axis can't deliver),
this asks the orthogonal, more answerable question the user raised: **does each breast cell line look
like its tumour subtype, in the framework's own quantities, and does that hold across modalities?**

For each TCGA PAM50 subtype we build a centroid in two comparable Hallmark-gene spaces:
  (1) **mRNA** target expression (the readout layer; TCGA RNA vs DepMap RNA-seq -- reliable on both sides)
  (2) **miRNA PRESSURE** (the regulation layer; `compute_gene_pressure` projects each platform's miRNA
       onto the SAME gene space, so TCGA-MIMAT and CCLE-NanoString become directly comparable)
Then each curated breast line (TAD-aligned, reliable labels) is Spearman-correlated to all 4 centroids;
its **best-match** subtype is compared to its **nominal** subtype. Per-modality accuracy + confusion.

Reading: mRNA is the positive control (cell lines are known to track tumour subtype transcriptionally);
if **pressure** also places lines with their subtype, the framework's miRNA-regulation layer carries
subtype identity in cell lines too -- and any modality where the match BREAKS localises the weak link
(expected: pressure < mRNA, because the cell-line miRNA is NanoString, cf. MH-67).

Run: ``.venv/bin/python3 -m mirna_hallmark.cell_line_subtype_fidelity``
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import pressure_build as PB
from mirna_hallmark.analyses.cnv_locus.ccle_mirna_cn_concordance import _load_model_lineage
from mirna_hallmark.analyses.cnv_locus.ccle_breast_target_anticorr import _mrna_breast_matrix, _curated_labels

OUT_DIR = C.OUTPUT_ROOT / "cell_line_subtype_fidelity"
SUBTYPES = ["LumA", "LumB", "HER2", "Basal"]
PAM50_NORM = {"Her2": "HER2", "HER2": "HER2", "LumA": "LumA", "LumB": "LumB", "Basal": "Basal"}


def _ccle_arm_name_matrix(model: pd.DataFrame) -> pd.DataFrame:
    """CCLE NanoString miRNA as arm-NAME x ModelID, log2(50+nSolver) -- engine-ready for pressure."""
    gct = pd.read_csv(C.CCLE_MIRNA_GCT, sep="\t", skiprows=2, index_col=0)
    ccle2model = model.set_index("CCLEName")["ModelID"].to_dict()
    cols = [c for c in gct.columns if c != "Description" and c in ccle2model]
    desc = gct["Description"].astype(str)
    mat = gct[cols].apply(pd.to_numeric, errors="coerce")
    mat.index = desc.values
    mat = np.log2(C.CCLE_MIRNA_EXPR_PSEUDOCOUNT + mat)
    mat = mat.groupby(mat.index).mean()
    mat.columns = [ccle2model[c] for c in cols]
    return mat


def _centroids(mat: pd.DataFrame, labels: pd.Series) -> pd.DataFrame:
    """gene x subtype centroid = per-subtype mean of per-gene z-scored values (pattern, not level)."""
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
    out = {}
    for st in SUBTYPES:
        cols = [c for c in z.columns if labels.get(c) == st]
        if cols:
            out[st] = z[cols].mean(axis=1)
    return pd.DataFrame(out)


def _match_lines(line_mat: pd.DataFrame, centroids: pd.DataFrame, nominal: pd.Series) -> pd.DataFrame:
    """Spearman each cell line to each subtype centroid over shared genes; best-match vs nominal."""
    shared = sorted(set(line_mat.index) & set(centroids.index))
    rows = []
    for line in line_mat.columns:
        x = line_mat[line].reindex(shared)
        ok = x.notna()
        if ok.sum() < 30:
            continue
        cors = {}
        for st in centroids.columns:
            y = centroids[st].reindex(shared)
            m = ok & y.notna()
            if m.sum() >= 30:
                cors[st] = float(spearmanr(x[m], y[m]).statistic)
        if not cors:
            continue
        best = max(cors, key=cors.get)
        nom = nominal.get(line)
        own = cors.get(nom, np.nan)
        others = [v for k, v in cors.items() if k != nom]
        rows.append({
            "line": line, "nominal": nom, "best_match": best, "agree": best == nom,
            "corr_own": round(own, 3) if pd.notna(own) else np.nan,
            "margin_own_minus_bestother": round(own - max(others), 3) if pd.notna(own) and others else np.nan,
            **{f"corr_{st}": round(cors.get(st, np.nan), 3) for st in SUBTYPES},
        })
    return pd.DataFrame(rows)


def _tcga_self_fidelity(mat: pd.DataFrame, labels: pd.Series, seed: int = 0) -> dict:
    """Control: do TUMOURS themselves separate by subtype in this modality? 50/50 split, build
    centroids on train, best-match accuracy on held-out test. Disambiguates 'modality not subtype-
    informative' from 'cell-line measurement broken' (if TCGA separates but cell lines don't)."""
    rng = np.random.default_rng(seed)
    cols = [c for c in mat.columns if c in labels.index]
    rng.shuffle(cols)
    train, test = cols[: len(cols) // 2], cols[len(cols) // 2:]
    cent = _centroids(mat[train], labels)
    z = mat[test].sub(mat[test].mean(axis=1), axis=0).div(mat[test].std(axis=1).replace(0, np.nan), axis=0)
    m = _match_lines(z, cent, labels)
    if m.empty:
        return {"n": 0}
    return {"n": int(len(m)), "accuracy": round(float(m["agree"].mean()), 3),
            "mean_corr_own": round(float(m["corr_own"].mean()), 3)}


def run(*, out_dir: Path = OUT_DIR, gene_cap: int | None = None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    model = _load_model_lineage()
    breast_ids = set(model.loc[model["OncotreeLineage"] == "Breast", "ModelID"])

    # ---- nominal labels for curated breast lines (reliable ground truth) ----------------------
    curated = _curated_labels(model, breast_ids).map(lambda v: PAM50_NORM.get(v, "Unknown"))
    curated = curated[curated.isin(SUBTYPES)]

    # ---- TCGA subtype labels ------------------------------------------------------------------
    clin = D.load_clinical_strata().set_index("participant")
    tcga_lab = clin["PAM50_final"].map(lambda v: PAM50_NORM.get(v, "Unknown"))
    tcga_lab = tcga_lab[tcga_lab.isin(SUBTYPES)]

    # ---- gene universe = framework Hallmark target genes with miRTarBase edges ----------------
    edges = pd.read_csv(C.OUTPUT_ROOT / "edges" / "mirna_hallmark_edges.tsv.gz", sep="\t",
                        usecols=["gene"])
    genes = sorted(edges["gene"].dropna().unique())
    if gene_cap:
        genes = genes[:gene_cap]

    # ================= MODALITY 1: mRNA target expression =====================================
    rna = D.load_rna()
    rna = rna[~rna.index.duplicated(keep="first")]         # collapse duplicate gene symbols
    mrna_ccle = _mrna_breast_matrix(breast_ids)
    mrna_ccle = mrna_ccle[~mrna_ccle.index.duplicated(keep="first")]
    g_mrna = sorted(set(genes) & set(rna.index) & set(mrna_ccle.index))
    tcga_mrna_mat = rna.loc[g_mrna, [c for c in rna.columns if c in tcga_lab.index]]
    tcga_cent_mrna = _centroids(tcga_mrna_mat, tcga_lab)
    ccle_lines_mrna = mrna_ccle.loc[g_mrna, [c for c in mrna_ccle.columns if c in curated.index]]
    # z-score cell lines on the same gene axis (pattern), then match
    z_ccle_mrna = ccle_lines_mrna.sub(ccle_lines_mrna.mean(axis=1), axis=0).div(
        ccle_lines_mrna.std(axis=1).replace(0, np.nan), axis=0)
    match_mrna = _match_lines(z_ccle_mrna, tcga_cent_mrna, curated)
    match_mrna["modality"] = "mRNA"

    # ================= MODALITY 2: miRNA pressure =============================================
    print("[fidelity] computing TCGA gene pressure ...")
    tcga_press = PB.compute_gene_pressure(g_mrna)
    print("[fidelity] computing cell-line gene pressure ...")
    ccle_arms = _ccle_arm_name_matrix(model)
    ccle_press = PB.compute_gene_pressure(g_mrna, mirna=ccle_arms)
    g_press = sorted(set(tcga_press.index) & set(ccle_press.index))
    tcga_press_mat = tcga_press.loc[g_press, [c for c in tcga_press.columns if c in tcga_lab.index]]
    tcga_cent_press = _centroids(tcga_press_mat, tcga_lab)
    ccle_lines_press = ccle_press.loc[g_press, [c for c in ccle_press.columns if c in curated.index]]
    z_ccle_press = ccle_lines_press.sub(ccle_lines_press.mean(axis=1), axis=0).div(
        ccle_lines_press.std(axis=1).replace(0, np.nan), axis=0)
    match_press = _match_lines(z_ccle_press, tcga_cent_press, curated)
    match_press["modality"] = "pressure"

    allm = pd.concat([match_mrna, match_press], ignore_index=True)
    name = {m: str(model.set_index("ModelID")["CCLEName"].get(m, m)).split("_")[0] for m in allm["line"]}
    allm.insert(1, "line_name", allm["line"].map(name))
    # coarse (Luminal = LumA+LumB) -- the LumA/LumB split is hard even for tumours
    coarse = {"LumA": "Luminal", "LumB": "Luminal", "HER2": "HER2", "Basal": "Basal"}
    allm["agree_coarse"] = allm["nominal"].map(coarse) == allm["best_match"].map(coarse)
    allm.to_csv(out_dir / "line_subtype_match.tsv", sep="\t", index=False)

    # ---- per-modality summary + confusion ----------------------------------------------------
    summary = {"module": "mirna_hallmark.cell_line_subtype_fidelity",
               "generated_utc": datetime.now(timezone.utc).isoformat(),
               "n_genes_mrna": len(g_mrna), "n_genes_pressure": len(g_press),
               "n_curated_lines": int(curated.isin(SUBTYPES).sum()),
               "tcga_subtype_n": tcga_lab.value_counts().to_dict(),
               "tcga_self_fidelity_control": {
                   "mRNA": _tcga_self_fidelity(tcga_mrna_mat, tcga_lab),
                   "pressure": _tcga_self_fidelity(tcga_press_mat, tcga_lab),
                   "note": "tumour-only 50/50 split: does the modality separate subtypes among TUMOURS? "
                           "If pressure separates here but fails for cell lines -> the cell-line miRNA "
                           "(NanoString) is the weak link, not pressure being subtype-uninformative.",
               },
               "by_modality": {}}
    for mod, sub in allm.groupby("modality"):
        valid = sub.dropna(subset=["best_match"])
        conf = pd.crosstab(valid["nominal"], valid["best_match"]).reindex(index=SUBTYPES, columns=SUBTYPES).fillna(0).astype(int)
        summary["by_modality"][mod] = {
            "n_lines": int(len(valid)),
            "accuracy_best_match_eq_nominal": round(float(valid["agree"].mean()), 3),
            "accuracy_coarse_luminal_her2_basal": round(float(valid["agree_coarse"].mean()), 3),
            "mean_corr_own_subtype": round(float(valid["corr_own"].mean()), 3),
            "mean_margin": round(float(valid["margin_own_minus_bestother"].mean()), 3),
            "confusion_nominal_x_bestmatch": conf.to_dict(),
            "misassigned": valid.loc[~valid["agree"], ["line_name", "nominal", "best_match", "corr_own"]].to_dict("records"),
        }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    ctrl = summary["tcga_self_fidelity_control"]
    print(f"[fidelity] TCGA-self control (tumour 50/50): mRNA acc={ctrl['mRNA'].get('accuracy')} "
          f"(own ρ {ctrl['mRNA'].get('mean_corr_own')}) | pressure acc={ctrl['pressure'].get('accuracy')} "
          f"(own ρ {ctrl['pressure'].get('mean_corr_own')})")
    for mod, d in summary["by_modality"].items():
        print(f"[fidelity] {mod}: accuracy={d['accuracy_best_match_eq_nominal']:.0%} "
              f"(coarse {d['accuracy_coarse_luminal_her2_basal']:.0%}; {d['n_lines']} lines), "
              f"own-subtype ρ={d['mean_corr_own_subtype']}, margin={d['mean_margin']}")
        if d["misassigned"]:
            print(f"[fidelity]   misassigned: " +
                  "; ".join(f"{m['line_name']} {m['nominal']}->{m['best_match']}" for m in d["misassigned"]))
    print(f"[fidelity] wrote {out_dir}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--gene-cap", type=int, default=None, help="debug: cap gene universe for speed")
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, gene_cap=args.gene_cap)


if __name__ == "__main__":
    main()
