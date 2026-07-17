"""Independent-cohort validation of miRNA→target coupling in the Buffa 2011 cohort.

The subproject's anti-correlation findings (miRNA pressure ↔ lower target mRNA) and
the TargetScan-orphan coupling story (MH-23/26) live in TCGA-BRCA. CPTAC adds an
orthogonal *protein* layer (`cptac_validation`). This module adds the missing leg:
**replication in an independent breast cohort of new patients** — Buffa et al. 2011
(GEO GSE22216 miRNA + GSE22219 mRNA, 210 paired primary breast tumours, Illumina),
which has **no overlap with TCGA patients**.

Mirrors the CPTAC "RNA part" (per-edge miRNA-arm ↔ target-mRNA anti-correlation,
raw + proliferation-adjusted partial Spearman) for two edge sets:
  1. **High-evidence miRTarBase edges** — compared edge-by-edge to the TCGA partial-ρ
     (`edge_partial_corr_panels`): sign concordance + replication rate of TCGA
     neg-sig edges.
  2. **TargetScan-orphan edges** (TS-predicted, not HE; MH-23) — does the HE > orphan
     coupling gap reproduce independently?

Honest scope / caveats:
  - Illumina miRNA uses **legacy bare arm names** (`hsa-miR-21`, no -5p/-3p); edge arms
    are mapped by bare stem, so the measured arm is assumed to be the edge's arm — true
    for guide-arm edges, approximate for passenger-arm edges (flagged).
  - Buffa has no purity/HRD; the only adjustment is a **proliferation proxy** (E2F/G2M
    mean-z from the cohort's own mRNA). Raw ρ reported alongside.
  - This validates *direction/replication of coupling*, an independent-patient leg the
    same-patient CPTAC protein layer cannot provide.

Run:
  .venv/bin/python3 -m mirna_hallmark.eval.buffa_validation
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr
from analysis.utils.common.loaders import partial_spearman

GEO_DIR = Path(__file__).resolve().parent.parent / "data" / "external_cache" / "geo"
BUFFA_MIRNA = GEO_DIR / "GSE22216_mirna_arm_matrix_paired.tsv.gz"
BUFFA_MRNA = GEO_DIR / "GSE22219_mrna_gene_matrix.tsv.gz"
TCGA_CORR = (C.TISSUE_REFERENCE_DIR / "edge_partial_corr_panels"
             / "edge_partial_corr_panel_corr_table.tsv")
CPTAC_ORPHAN = (C.OUTPUT_ROOT / "cptac_validation" / "orphan_discovery"
                / "orphan_candidates.tsv")
OUT_DIR = C.OUTPUT_ROOT / "buffa_validation"
MIN_N = 30


# --------------------------------------------------------------------------- #
# Load cohort
# --------------------------------------------------------------------------- #

BUFFA_MIRNA_SERIES = GEO_DIR / "GSE22216_miRNA_series_matrix.txt.gz"


def load_buffa() -> Tuple[pd.DataFrame, pd.DataFrame]:
    mi = pd.read_csv(BUFFA_MIRNA, sep="\t", index_col=0)
    rna = pd.read_csv(BUFFA_MRNA, sep="\t", index_col=0)
    shared = mi.columns.intersection(rna.columns)
    print(f"[buffa] miRNA {mi.shape}  mRNA {rna.shape}  paired samples {len(shared)}")
    return mi[shared], rna[shared]


def load_buffa_er() -> pd.Series:
    """ER status (0/1) per miRNA-GSM patient key, from the GSE22216 characteristics."""
    import gzip
    gsms = None
    er = None
    with gzip.open(BUFFA_MIRNA_SERIES, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                gsms = [x.strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and "er status" in line:
                er = [x.strip('"').replace("er status:", "").strip() for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                break
    s = pd.Series(pd.to_numeric(pd.Series(er, index=gsms), errors="coerce"), name="er")
    return s


def _stem(a: str) -> str:
    return re.sub(r"-(5p|3p)$", "", str(a))


def prolif_proxy(rna: pd.DataFrame, hs: HallmarkSets) -> pd.Series:
    """Per-sample mean z of E2F_TARGETS ∪ G2M_CHECKPOINT genes (the MH-17 confound)."""
    genes = set(hs.sets.get("HALLMARK_E2F_TARGETS", [])) | set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))
    present = [g for g in genes if g in rna.index]
    sub = rna.loc[present].apply(pd.to_numeric, errors="coerce")
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1) + 1e-9, axis=0)
    return z.mean(axis=0).rename("prolif")


# --------------------------------------------------------------------------- #
# Edge coupling in Buffa
# --------------------------------------------------------------------------- #

def edge_coupling(
    edges: pd.DataFrame,
    mirna: pd.DataFrame,
    rna: pd.DataFrame,
    cov: pd.DataFrame,
) -> pd.DataFrame:
    """Per-edge raw + proliferation-adjusted partial Spearman (arm vs target) in Buffa."""
    arm_by_stem: Dict[str, str] = {}
    for a in mirna.index:
        arm_by_stem.setdefault(_stem(a), a)
    rna_genes = set(rna.index)

    rows: List[dict] = []
    for arm, gene in zip(edges["miRNA"], edges["gene"]):
        if gene not in rna_genes:
            continue
        b_arm = arm_by_stem.get(_stem(arm))
        if b_arm is None:
            continue
        x = pd.to_numeric(mirna.loc[b_arm], errors="coerce")
        y = pd.to_numeric(rna.loc[gene] if not isinstance(rna.loc[gene], pd.DataFrame)
                          else rna.loc[gene].iloc[0], errors="coerce")
        from scipy.stats import spearmanr
        df = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
        if len(df) < MIN_N:
            continue
        raw_rho, raw_p = spearmanr(df["x"], df["y"])
        prho, pp, pn = partial_spearman(y, x, cov)
        rows.append({
            "miRNA": arm, "gene": gene, "buffa_arm": b_arm,
            "n": len(df), "raw_rho": round(float(raw_rho), 4), "raw_p": float(raw_p),
            "partial_rho": round(float(prho), 4) if np.isfinite(prho) else np.nan,
            "partial_p": pp, "partial_n": pn,
            "passenger_arm_caveat": str(arm).endswith("-3p"),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["raw_q"] = bh_fdr(out["raw_p"].fillna(1.0).values)
        valid = out["partial_p"].notna()
        out["partial_q"] = np.nan
        if valid.any():
            out.loc[valid, "partial_q"] = bh_fdr(out.loc[valid, "partial_p"].values)
    return out


def _neg_sig_frac(tbl: pd.DataFrame, rho_col: str, q_col: str) -> float:
    if tbl.empty:
        return float("nan")
    return float(((tbl[rho_col] < 0) & (tbl[q_col] < 0.05)).mean())


# --------------------------------------------------------------------------- #
# Triple-cohort validation of genuine orphans (absent from miRTarBase)
# --------------------------------------------------------------------------- #

def validate_orphans_in_buffa(
    mirna: pd.DataFrame, rna: pd.DataFrame, cov: pd.DataFrame, out_dir: Path,
) -> Dict[str, object]:
    """Test the CPTAC genuine-orphan edges (absent from miRTarBase, coupled in TCGA+CPTAC)
    in the independent Buffa mRNA cohort → triple-cohort (TCGA·CPTAC·Buffa) validation.

    'Genuine orphan' = `mirtar_any==False` (no miRTarBase study at all): TargetScan/ENCORI-
    predicted edges with NO curation, yet protein-coupled in CPTAC and mRNA-coupled in TCGA.
    Buffa adds a third, independent-patient mRNA readout. Annotated with miRGeneDB arm class
    (is the arm a high-confidence guide?).
    """
    if not CPTAC_ORPHAN.exists():
        print(f"[buffa] CPTAC orphan table absent ({CPTAC_ORPHAN}) — skipping triple-cohort")
        return {}
    cand = pd.read_csv(CPTAC_ORPHAN, sep="\t")
    genuine = cand[~cand["mirtar_any"].astype(bool)].copy()
    print(f"[buffa] CPTAC orphan candidates: {len(cand)}; genuine (absent from miRTarBase): {len(genuine)}")

    edges = genuine[["miRNA", "gene"]].drop_duplicates()
    buffa = edge_coupling(edges, mirna, rna, cov)[
        ["miRNA", "gene", "raw_rho", "raw_q", "partial_rho", "partial_q", "n"]
    ].rename(columns={"raw_rho": "buffa_raw_rho", "raw_q": "buffa_raw_q",
                      "partial_rho": "buffa_partial_rho", "partial_q": "buffa_partial_q",
                      "n": "buffa_n"})
    m = genuine.merge(buffa, on=["miRNA", "gene"], how="left")

    # miRGeneDB arm class (the high-confidence "~65% of miRTarBase arms" reference)
    try:
        from mirna_hallmark.analyses.edge_panels.edge_prior_refinement import mirgenedb_arm_status
        arm = mirgenedb_arm_status(do_fetch=False)
        armmap = dict(zip(arm["miRNA"], arm["arm_class"]))
        m["arm_class"] = m["miRNA"].map(armmap).fillna("unknown")
    except Exception:  # noqa: BLE001
        m["arm_class"] = "unknown"

    m["buffa_neg"] = m["buffa_partial_rho"] < 0
    m["buffa_neg_sig"] = (m["buffa_partial_rho"] < 0) & (m["buffa_partial_q"] < 0.05)
    # triple cohort: CPTAC-sig (protein_q) AND TCGA-replicated AND Buffa same-sign (neg)
    m["cptac_sig"] = m["protein_q"] < 0.05
    m["triple_validated"] = m["cptac_sig"] & m["tcga_replicated"].astype(bool) & m["buffa_neg"]
    m["triple_validated_sig"] = m["cptac_sig"] & m["tcga_replicated"].astype(bool) & m["buffa_neg_sig"]

    m = m.sort_values(["triple_validated", "buffa_partial_rho"], ascending=[False, True])
    m.to_csv(out_dir / "orphan_triple_cohort_validation.tsv", sep="\t", index=False)

    core = m[m["cptac_sig"] & m["tcga_replicated"].astype(bool) & m["buffa_partial_rho"].notna()]
    n_core = len(core)
    n_buffa_neg = int((core["buffa_partial_rho"] < 0).sum())
    n_triple = int(m["triple_validated"].sum())
    n_triple_sig = int(m["triple_validated_sig"].sum())
    res = {
        "n_genuine_orphans": int(len(genuine)),
        "n_testable_in_buffa": int(m["buffa_partial_rho"].notna().sum()),
        "n_cptac_sig_and_tcga_repl_testable": n_core,
        "n_buffa_same_sign": n_buffa_neg,
        "buffa_replication_rate_of_cptac_tcga_orphans": round(n_buffa_neg / n_core, 3) if n_core else float("nan"),
        "n_triple_validated": n_triple,
        "n_triple_validated_sig_in_buffa": n_triple_sig,
        "arm_class_of_triple": m.loc[m["triple_validated"], "arm_class"].value_counts().to_dict(),
    }
    print(f"[buffa] genuine orphans testable in Buffa: {res['n_testable_in_buffa']}")
    print(f"[buffa] of CPTAC-sig+TCGA-repl orphans testable in Buffa ({n_core}): "
          f"{n_buffa_neg} same-sign ({res['buffa_replication_rate_of_cptac_tcga_orphans']:.0%})")
    print(f"[buffa] TRIPLE-cohort validated (TCGA·CPTAC·Buffa, same sign): {n_triple}  "
          f"(neg-sig in Buffa: {n_triple_sig})")
    print("[buffa] top triple-validated orphans:")
    for _, r in m[m["triple_validated"]].head(12).iterrows():
        print(f"    {r['miRNA']:<18}→{r['gene']:<8} "
              f"buffa ρ={r['buffa_partial_rho']:+.2f} cptac ρ={r['protein_rho']:+.2f} "
              f"tcga ρ={r['tcga_protein_rho']:+.2f}  arm={r['arm_class']}")
    return res


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def run(*, out_dir: Path = OUT_DIR) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    mirna, rna = load_buffa()
    cov = prolif_proxy(rna, hs).to_frame()
    er = load_buffa_er().reindex(mirna.columns)
    cov["er"] = er  # ER status = subtype proxy (Buffa lacks PAM50/purity/HRD)
    print(f"[buffa] covariates: prolif (E2F/G2M) + ER status (n_ER={int(er.notna().sum())})")

    # ---- HE edges ---- #
    edges_all = D.load_hallmark_edges()
    he = D.high_evidence_edges(edges_all)[["miRNA", "gene"]].drop_duplicates()
    print(f"[buffa] HE edges: {len(he)} ; computing Buffa coupling …")
    he_buffa = edge_coupling(he, mirna, rna, cov)
    print(f"[buffa] HE edges with Buffa coupling: {len(he_buffa)}")

    # compare to TCGA partial-ρ
    tcga = pd.read_csv(TCGA_CORR, sep="\t")[["miRNA", "gene", "rho_adj", "q_adj"]].dropna(subset=["rho_adj"])
    cmp = he_buffa.merge(tcga, on=["miRNA", "gene"], how="inner", suffixes=("", "_tcga"))
    cmp = cmp.rename(columns={"rho_adj": "tcga_rho_adj", "q_adj": "tcga_q_adj"})
    cmp.to_csv(out_dir / "he_edge_coupling_buffa_vs_tcga.tsv", sep="\t", index=False)
    he_buffa.to_csv(out_dir / "he_edge_coupling_buffa.tsv", sep="\t", index=False)

    # replication metrics on the edges with both
    from scipy.stats import spearmanr as _sp
    both = cmp.dropna(subset=["partial_rho", "tcga_rho_adj"])
    rho_concord = float(_sp(both["partial_rho"], both["tcga_rho_adj"])[0]) if len(both) > 10 else float("nan")
    tcga_negsig = both[(both["tcga_rho_adj"] < 0) & (both["tcga_q_adj"] < 0.05)]
    repl_rate = float((tcga_negsig["partial_rho"] < 0).mean()) if len(tcga_negsig) else float("nan")
    repl_rate_sig = float(((tcga_negsig["partial_rho"] < 0) & (tcga_negsig["partial_q"] < 0.05)).mean()) \
        if len(tcga_negsig) else float("nan")
    sign_concord = float((np.sign(both["partial_rho"]) == np.sign(both["tcga_rho_adj"])).mean()) if len(both) else float("nan")

    # ---- orphan TargetScan edges ---- #
    from mirna_hallmark.eval.targetscan_orphan_coupling import build_orphan_edge_table, MIN_TS_WEIGHT
    he_pairs: Set[Tuple[str, str]] = set(zip(he["miRNA"], he["gene"]))
    orphan = build_orphan_edge_table(sorted(hs.universe), he_pairs, min_ts=MIN_TS_WEIGHT)
    orphan = orphan[["miRNA", "gene"]].drop_duplicates()
    print(f"[buffa] orphan TS edges: {len(orphan)} ; computing Buffa coupling …")
    orphan_buffa = edge_coupling(orphan, mirna, rna, cov)
    orphan_buffa.to_csv(out_dir / "orphan_edge_coupling_buffa.tsv", sep="\t", index=False)
    print(f"[buffa] orphan edges with Buffa coupling: {len(orphan_buffa)}")

    he_negsig = _neg_sig_frac(he_buffa, "partial_rho", "partial_q")
    orphan_negsig = _neg_sig_frac(orphan_buffa, "partial_rho", "partial_q")

    # ---- triple-cohort validation of genuine orphans (TCGA·CPTAC·Buffa) ---- #
    triple = validate_orphans_in_buffa(mirna, rna, cov, out_dir)

    summary = {
        "module": "mirna_hallmark.eval.buffa_validation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort": "Buffa 2011 (GSE22216 miRNA + GSE22219 mRNA), 210 paired breast tumours, "
                  "independent of TCGA",
        "confounder": "proliferation proxy (E2F/G2M mean-z from cohort mRNA); no purity/HRD in Buffa",
        "arm_mapping": "legacy bare-name → edge arm by bare stem (guide-arm-accurate)",
        "he_edges": {
            "n_with_buffa_coupling": int(len(he_buffa)),
            "n_compared_to_tcga": int(len(both)),
            "buffa_vs_tcga_rho_concordance": round(rho_concord, 3),
            "sign_concordance": round(sign_concord, 3),
            "tcga_negsig_replication_rate_same_sign": round(repl_rate, 3),
            "tcga_negsig_replication_rate_sig": round(repl_rate_sig, 3),
            "buffa_neg_sig_frac": round(he_negsig, 3),
        },
        "orphan_edges": {
            "n_with_buffa_coupling": int(len(orphan_buffa)),
            "buffa_neg_sig_frac": round(orphan_negsig, 3),
            "he_vs_orphan_neg_sig": f"HE {he_negsig:.3f} vs orphan {orphan_negsig:.3f} "
                                    "(replicates MH-23 HE>orphan gap if HE higher)",
        },
        "genuine_orphan_triple_cohort": triple,
        "outputs": [
            "he_edge_coupling_buffa.tsv", "he_edge_coupling_buffa_vs_tcga.tsv",
            "orphan_edge_coupling_buffa.tsv", "orphan_triple_cohort_validation.tsv",
        ],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")

    # replication scatter: Buffa partial-ρ vs TCGA partial-ρ
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.5, 7))
        sig = (both["tcga_q_adj"] < 0.05)
        ax.scatter(both.loc[~sig, "tcga_rho_adj"], both.loc[~sig, "partial_rho"],
                   s=6, c="#cccccc", alpha=0.4, label="TCGA q≥0.05")
        ax.scatter(both.loc[sig, "tcga_rho_adj"], both.loc[sig, "partial_rho"],
                   s=10, c="#4E79A7", alpha=0.55, label="TCGA q<0.05")
        ax.axhline(0, color="#888", lw=0.7, ls="--"); ax.axvline(0, color="#888", lw=0.7, ls="--")
        lim = 0.6
        ax.plot([-lim, lim], [-lim, lim], color="#E15759", lw=0.8, ls=":", label="y=x")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xlabel("TCGA partial ρ (CPE+HRD+prolif)", fontsize=10)
        ax.set_ylabel("Buffa partial ρ (ER+prolif)", fontsize=10)
        ax.set_title(
            f"Independent-cohort replication of edge coupling\n"
            f"n={len(both)} edges · concordance ρ={rho_concord:+.2f} · "
            f"TCGA neg-sig same-sign in Buffa {repl_rate:.0%}",
            fontsize=10)
        ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.15, ls="--")
        fig.tight_layout()
        (out_dir / "figures").mkdir(exist_ok=True)
        fig.savefig(out_dir / "figures" / "buffa_vs_tcga_replication.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[buffa] figure -> {out_dir/'figures'/'buffa_vs_tcga_replication.png'}")
    except Exception as e:  # noqa: BLE001
        print(f"[buffa] figure skipped: {e}")

    print("\n=== Buffa independent-cohort validation ===")
    print(f"HE edges: {len(he_buffa)} coupled; {len(both)} also in TCGA")
    print(f"  Buffa-vs-TCGA partial-ρ concordance (Spearman): {rho_concord:+.3f}")
    print(f"  sign concordance: {sign_concord:.3f}")
    print(f"  TCGA neg-sig edges replicating (same sign in Buffa): {repl_rate:.3f}")
    print(f"  TCGA neg-sig edges replicating (neg-sig in Buffa):   {repl_rate_sig:.3f}")
    print(f"HE neg-sig frac {he_negsig:.3f}  vs  orphan neg-sig frac {orphan_negsig:.3f}")
    return {"he_buffa": he_buffa, "cmp": cmp, "orphan_buffa": orphan_buffa, "summary": summary}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--min-purity", type=float, default=None, help="runner compat (no-op)")
    args = ap.parse_args()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
