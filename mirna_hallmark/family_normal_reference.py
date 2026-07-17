"""miR-301/130/454 family: tumor-adjacent-normal (NAT) reference contrasts.

Answers the design question from FORMULAS.md / chat: should the abundance anchor
(cohort median) be replaced or augmented with a *healthy-tissue* reference, and is
the "miR-130a-3p is the abundant-but-uncoupled family hub" picture **tumor-acquired**
or already present in normal breast?

It deliberately keeps full TCGA barcodes (the default `mirna_hallmark` loader
truncates to 12-char participant and *averages tumor + NAT together*, which would
silently dilute the contrast). Sample types: ``01`` primary tumor, ``11`` solid
tissue normal (NAT), ``06`` metastatic (ignored).

Three outputs:
1. ``family_arm_tumor_vs_nat.tsv`` — per family arm: median log2(RPM+1) in tumor vs
   NAT, log2 fold change, Mann-Whitney U + BH-FDR, and abundance rank within the
   family in each compartment (is 130a the hub in normal too?).
2. ``family_normal_anchored_absratio.tsv`` — a normal-anchored abundance ratio
   (tumor median / NAT median, linear RPM) as an alternative to the cohort-median
   ``absratio`` denominator.
3. ``family_nat_coupling.tsv`` — DESCRIPTIVE raw Spearman of each family arm vs the
   focus genes within matched NAT samples (n~100, no purity/CN covariates — normal
   tissue), with the tumor-cohort CN-adjusted rho alongside and a sign-flip flag.

NAT coupling is under-powered and confounder-naive by construction; it is a
qualitative "same sign in normal?" check, not a claim.

Run:
  .venv/bin/python3 -m mirna_hallmark.family_normal_reference
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from mirna_hallmark import config as C
from mirna_hallmark.analyses.mir301.mir301_family_depth import FAMILY_ARMS, OUT_DIR as FAM_DEPTH_DIR, _focus_genes
from mirna_hallmark.stats import bh_fdr


OUT_DIR = C.TISSUE_REFERENCE_DIR / "family_normal_reference"
MIN_NAT_COUPLING_N = 30

# GTEx v10 true-healthy breast (small RNA-seq introduced in V10). Downloaded under
# data/GTEx/ (see scripts note in catalog).
#   miRNA: TCGA = RPM (GDC/MIMAT), GTEx = "TPM" (RNAcentral/URS) -> both per-million
#          but DIFFERENT pipelines + annotations; compare RANK / SHARE / coupling SIGN.
#   mRNA : TCGA = STAR TPM, GTEx = v10 TPM -> both TPM, broadly comparable (batch aside).
GTEX_DIR = C.MIRNA_EXPRESSION.parent.parent / "GTEx"
GTEX_MIRNA_TPM = GTEX_DIR / "miRNA_TPM_matrix_v10.txt.gz"
GTEX_SMALLRNA_ANNOT = GTEX_DIR / "smallRNA_annotation_v10.txt"
GTEX_SAMPLE_ATTR = GTEX_DIR / "GTEx_v10_SampleAttributesDS.txt"
GTEX_BREAST_TPM = GTEX_DIR / "gene_tpm_v10_breast.gct.gz"
GTEX_TISSUE = "Breast - Mammary Tissue"


def _gtex_donor(sampid: str) -> str:
    return "-".join(sampid.split("-")[:2])


def _sample_type(barcode: str) -> str:
    parts = barcode.split("-")
    return parts[3][:2] if len(parts) >= 4 else "?"


def _participant(barcode: str) -> str:
    return barcode[:12]


def _load_full_mirna() -> pd.DataFrame:
    """arm x FULL-barcode log2(RPM+1) (no participant collapse)."""
    from analysis.expression.mirna_target_integration import load_mimat_to_arm

    mimat_to_arm = load_mimat_to_arm(C.MIRNA_MATURE_LOCI)
    df = pd.read_csv(C.MIRNA_EXPRESSION, sep="\t", index_col=0)
    df.index = df.index.astype(str).map(lambda i: mimat_to_arm.get(i, i))
    df = df.groupby(level=0).mean()
    return df


def _load_full_rna(genes: Sequence[str]) -> pd.DataFrame:
    """gene x FULL-barcode log2(TPM+1), restricted to ``genes`` (symbol-indexed file)."""
    df = pd.read_csv(C.RNA_EXPRESSION, sep="\t", index_col=0)
    if "Ensembl_ID" in df.columns:
        df = df.drop(columns=["Ensembl_ID"])
    df = df.loc[df.index.isin(set(genes))]
    df = df.apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    return np.log2(df.clip(lower=0) + 1)


def _split_types(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    types = pd.Series({c: _sample_type(c) for c in df.columns})
    tumor = df.loc[:, types.eq("01").values]
    nat = df.loc[:, types.eq("11").values]
    # Cross-state matrices split by sample type, not PAM50, so the global Normal-like drop
    # (in load_clinical_strata) does not reach them — filter Normal-like primary tumours here.
    if getattr(C, "EXCLUDE_NORMAL_LIKE", False):
        from mirna_hallmark import data_loaders as _D
        norm = _D.normal_like_participants()
        if norm:
            tumor = tumor.loc[:, [c for c in tumor.columns if _participant(c) not in norm]]
    return {"tumor": tumor, "nat": nat}


def _family_abundance_contrast(mirna_full: pd.DataFrame) -> pd.DataFrame:
    split = _split_types(mirna_full)
    tumor, nat = split["tumor"], split["nat"]
    rows: List[dict] = []
    for arm in FAMILY_ARMS:
        if arm not in mirna_full.index:
            continue
        t = pd.to_numeric(tumor.loc[arm], errors="coerce").dropna()
        n = pd.to_numeric(nat.loc[arm], errors="coerce").dropna()
        if len(t) < 5 or len(n) < 5:
            continue
        try:
            u, p = mannwhitneyu(t, n, alternative="two-sided")
        except ValueError:
            u, p = np.nan, np.nan
        rows.append(
            {
                "arm": arm,
                "n_tumor": int(len(t)),
                "n_nat": int(len(n)),
                "tumor_median_log2rpm": float(t.median()),
                "nat_median_log2rpm": float(n.median()),
                "log2fc_tumor_minus_nat": float(t.median() - n.median()),
                "mwu_p": float(p),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["mwu_q"] = bh_fdr(out["mwu_p"].fillna(1.0).values)
    # abundance rank within family (1 = most abundant) in each compartment
    out["tumor_family_rank"] = out["tumor_median_log2rpm"].rank(ascending=False, method="min")
    out["nat_family_rank"] = out["nat_median_log2rpm"].rank(ascending=False, method="min")
    return out.sort_values("tumor_median_log2rpm", ascending=False).reset_index(drop=True)


def _normal_anchored_absratio(contrast: pd.DataFrame) -> pd.DataFrame:
    """tumor/NAT abundance ratio in linear RPM (alternative to cohort-median absratio)."""
    if contrast.empty:
        return contrast
    out = contrast[["arm", "tumor_median_log2rpm", "nat_median_log2rpm", "log2fc_tumor_minus_nat"]].copy()
    out["tumor_rpm"] = np.expm1(np.log(2.0) * out["tumor_median_log2rpm"])
    out["nat_rpm"] = np.expm1(np.log(2.0) * out["nat_median_log2rpm"])
    out["normal_anchored_absratio"] = out["tumor_rpm"] / out["nat_rpm"].replace(0, np.nan)
    return out


def _tumor_cohort_rho() -> Dict[tuple, float]:
    """(arm, gene) -> tumor cohort CN-adjusted single-arm rho from the depth module."""
    path = FAM_DEPTH_DIR / "family_all_partials.tsv"
    if not path.is_file():
        return {}
    p = pd.read_csv(path, sep="\t")
    sa = p.loc[p["predictor_type"].eq("single_arm") & p["scope"].eq("cohort")]
    col = "rho_CPE_HRD_CN" if "rho_CPE_HRD_CN" in sa.columns else "rho_CPE_HRD"
    return {(r["predictor"], r["gene"]): float(r[col]) for _, r in sa.iterrows() if pd.notna(r[col])}


def _nat_coupling(mirna_full: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    """DESCRIPTIVE raw Spearman family arm vs focus gene in matched NAT samples."""
    mir_nat = _split_types(mirna_full)["nat"]
    rna_nat = _split_types(_load_full_rna(genes))["nat"]
    # match on 12-char participant
    mir_nat = mir_nat.rename(columns={c: _participant(c) for c in mir_nat.columns})
    rna_nat = rna_nat.rename(columns={c: _participant(c) for c in rna_nat.columns})
    mir_nat = mir_nat.loc[:, ~mir_nat.columns.duplicated()]
    rna_nat = rna_nat.loc[:, ~rna_nat.columns.duplicated()]
    shared = sorted(set(mir_nat.columns) & set(rna_nat.columns))
    tumor_rho = _tumor_cohort_rho()
    rows: List[dict] = []
    for arm in FAMILY_ARMS:
        if arm not in mir_nat.index:
            continue
        x = pd.to_numeric(mir_nat.loc[arm, shared], errors="coerce")
        for gene in genes:
            if gene not in rna_nat.index:
                continue
            y = pd.to_numeric(rna_nat.loc[gene, shared], errors="coerce")
            ok = x.notna() & y.notna()
            if int(ok.sum()) < MIN_NAT_COUPLING_N:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            t_rho = tumor_rho.get((arm, gene), np.nan)
            rows.append(
                {
                    "arm": arm,
                    "gene": gene,
                    "n_nat": int(ok.sum()),
                    "rho_nat_raw": float(rho),
                    "p_nat_raw": float(p),
                    "rho_tumor_cohort_cn": t_rho,
                    "sign_flip_vs_tumor": bool(
                        pd.notna(t_rho) and np.sign(rho) != np.sign(t_rho)
                    ),
                }
            )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_nat_raw"] = bh_fdr(out["p_nat_raw"].fillna(1.0).values)
    return out.sort_values(["arm", "rho_nat_raw"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# GTEx v10 true-healthy breast reference
# --------------------------------------------------------------------------- #
def _gtex_available() -> bool:
    return all(p.is_file() for p in (GTEX_MIRNA_TPM, GTEX_SMALLRNA_ANNOT,
                                     GTEX_SAMPLE_ATTR, GTEX_BREAST_TPM))


def _gtex_urs_to_arm() -> Dict[str, str]:
    annot = pd.read_csv(GTEX_SMALLRNA_ANNOT, sep="\t", usecols=["id", "mirbase_name"],
                        low_memory=False)
    annot = annot.dropna(subset=["mirbase_name"])
    annot = annot.loc[annot["mirbase_name"].astype(str).str.startswith(("hsa-miR", "hsa-let"))]
    return dict(zip(annot["id"].astype(str), annot["mirbase_name"].astype(str)))


def _gtex_breast_samples(smafrze: str) -> List[str]:
    sa = pd.read_csv(GTEX_SAMPLE_ATTR, sep="\t", usecols=["SAMPID", "SMTSD", "SMAFRZE"],
                     low_memory=False)
    sub = sa.loc[sa["SMTSD"].eq(GTEX_TISSUE) & sa["SMAFRZE"].eq(smafrze)]
    return sub["SAMPID"].astype(str).tolist()


def _gtex_family_mirna() -> pd.DataFrame:
    """family arm x GTEx-breast-donor TPM (true healthy)."""
    urs_to_arm = _gtex_urs_to_arm()
    fam_urs = {u for u, a in urs_to_arm.items() if a in set(FAMILY_ARMS)}
    breast_smlrna = set(_gtex_breast_samples("SMLRNA"))
    rows = []
    for chunk in pd.read_csv(GTEX_MIRNA_TPM, sep="\t", index_col=0, chunksize=400, low_memory=False):
        keep = chunk.loc[chunk.index.astype(str).isin(fam_urs)]
        if not keep.empty:
            rows.append(keep)
    if not rows:
        return pd.DataFrame()
    mat = pd.concat(rows)
    mat.index = mat.index.astype(str).map(urs_to_arm)
    breast_cols = [c for c in mat.columns if c in breast_smlrna]
    mat = mat[breast_cols].groupby(level=0).mean()
    mat.columns = [_gtex_donor(c) for c in mat.columns]
    return mat.loc[:, ~mat.columns.duplicated()]


def _gtex_breast_rna(genes: Sequence[str]) -> pd.DataFrame:
    """gene x GTEx-breast-donor log2(TPM+1) for focus genes."""
    gset = set(genes)
    rows = []
    for chunk in pd.read_csv(GTEX_BREAST_TPM, sep="\t", skiprows=2, chunksize=4000, low_memory=False):
        sub = chunk.loc[chunk["Description"].isin(gset)]
        if not sub.empty:
            rows.append(sub)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows)
    df = df.drop(columns=["Name"]).set_index("Description")
    df = df.apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    df.columns = [_gtex_donor(c) for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return np.log2(df.clip(lower=0) + 1)


def _gtex_family_abundance(gtex_mir: pd.DataFrame) -> pd.DataFrame:
    if gtex_mir.empty:
        return gtex_mir
    med = gtex_mir.median(axis=1)
    out = pd.DataFrame({"arm": med.index, "gtex_median_tpm": med.values})
    out["gtex_family_rank"] = out["gtex_median_tpm"].rank(ascending=False, method="min")
    out["gtex_family_share"] = out["gtex_median_tpm"] / out["gtex_median_tpm"].sum()
    out["n_gtex_breast"] = int(gtex_mir.shape[1])
    return out.sort_values("gtex_median_tpm", ascending=False).reset_index(drop=True)


def _gtex_coupling(gtex_mir: pd.DataFrame, genes: Sequence[str]) -> pd.DataFrame:
    rna = _gtex_breast_rna(genes)
    if gtex_mir.empty or rna.empty:
        return pd.DataFrame()
    shared = sorted(set(gtex_mir.columns) & set(rna.columns))
    tumor_rho = _tumor_cohort_rho()
    rows: List[dict] = []
    for arm in FAMILY_ARMS:
        if arm not in gtex_mir.index:
            continue
        x = pd.to_numeric(gtex_mir.loc[arm, shared], errors="coerce")
        for gene in genes:
            if gene not in rna.index:
                continue
            y = pd.to_numeric(rna.loc[gene, shared], errors="coerce")
            ok = x.notna() & y.notna()
            if int(ok.sum()) < MIN_NAT_COUPLING_N:
                continue
            rho, p = spearmanr(x[ok], y[ok])
            t_rho = tumor_rho.get((arm, gene), np.nan)
            rows.append({
                "arm": arm, "gene": gene, "n_gtex": int(ok.sum()),
                "rho_gtex_raw": float(rho), "p_gtex_raw": float(p),
                "rho_tumor_cohort_cn": t_rho,
                "sign_flip_vs_tumor": bool(pd.notna(t_rho) and np.sign(rho) != np.sign(t_rho)),
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q_gtex_raw"] = bh_fdr(out["p_gtex_raw"].fillna(1.0).values)
    return out.sort_values(["arm", "rho_gtex_raw"]).reset_index(drop=True)


def run(*, out_dir: Path = OUT_DIR, top_orphans: int = 25, include_gtex: bool = True) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = _focus_genes(top_orphans=top_orphans)
    mirna_full = _load_full_mirna()

    contrast = _family_abundance_contrast(mirna_full)
    absratio = _normal_anchored_absratio(contrast)
    coupling = _nat_coupling(mirna_full, genes)

    contrast.to_csv(out_dir / "family_arm_tumor_vs_nat.tsv", sep="\t", index=False)
    absratio.to_csv(out_dir / "family_normal_anchored_absratio.tsv", sep="\t", index=False)
    coupling.to_csv(out_dir / "family_nat_coupling.tsv", sep="\t", index=False)

    gtex_abund = pd.DataFrame()
    gtex_coupling = pd.DataFrame()
    if include_gtex and _gtex_available():
        gtex_mir = _gtex_family_mirna()
        gtex_abund = _gtex_family_abundance(gtex_mir)
        gtex_coupling = _gtex_coupling(gtex_mir, genes)
        gtex_abund.to_csv(out_dir / "gtex_family_abundance.tsv", sep="\t", index=False)
        gtex_coupling.to_csv(out_dir / "gtex_family_coupling.tsv", sep="\t", index=False)

    # cross-reference: family abundance rank across the three references
    xref = contrast[["arm", "tumor_family_rank", "nat_family_rank",
                     "tumor_median_log2rpm", "nat_median_log2rpm",
                     "log2fc_tumor_minus_nat"]].copy()
    if not gtex_abund.empty:
        xref = xref.merge(
            gtex_abund[["arm", "gtex_family_rank", "gtex_median_tpm", "gtex_family_share"]],
            on="arm", how="left",
        )
    xref.to_csv(out_dir / "family_abundance_rank_xref.tsv", sep="\t", index=False)

    manifest = {
        "module": "mirna_hallmark.family_normal_reference",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "family_arms": list(FAMILY_ARMS),
        "sample_types": {"tumor": "01", "nat": "11", "ignored": "06 metastatic"},
        "n_focus_genes": len(genes),
        "nat_coupling": "raw Spearman, no covariates, descriptive only",
        "gtex": {
            "available": bool(include_gtex and _gtex_available()),
            "tissue": GTEX_TISSUE,
            "version": "v10 small RNA-seq miRNA TPM; mRNA TPM v10",
            "note": "true healthy (not adjacent). miRNA cross-pipeline/annotation (RPM/MIMAT vs TPM/URS) -> rank/share/sign; mRNA TPM-vs-TPM broadly comparable",
        },
        "caveats": [
            "NAT is adjacent-normal (field effect), not true healthy tissue",
            "NAT n~100; coupling is confounder-naive and under-powered",
            "cohort default loader averages tumor+NAT per participant; this module does not",
            "miRNA GTEx vs TCGA: different pipeline + annotation (RPM/MIMAT vs TPM/URS); use rank/share/sign, not absolute levels. mRNA is TPM in both.",
        ],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print("[family_normal_reference] Tumor vs NAT abundance:\n", contrast.to_string(index=False))
    if not coupling.empty:
        flips = coupling.loc[coupling["sign_flip_vs_tumor"]]
        print(f"\n[family_normal_reference] NAT coupling rows={len(coupling)}; "
              f"sign-flip vs tumor={len(flips)}/{len(coupling)}; matched NAT n="
              f"{int(coupling['n_nat'].max())}")
    if not gtex_abund.empty:
        print("\n[family_normal_reference] GTEx healthy-breast family abundance:\n",
              gtex_abund.to_string(index=False))
        print("\n[family_normal_reference] Family abundance rank xref (tumor/NAT/GTEx):\n",
              xref.to_string(index=False))
    if not gtex_coupling.empty:
        gf = gtex_coupling.loc[gtex_coupling["sign_flip_vs_tumor"]]
        for arm, g in gtex_coupling.groupby("arm"):
            negsig = g[(g.rho_gtex_raw < 0) & (g.q_gtex_raw < 0.05)]
            print(f"[gtex-coupling] {arm:18s} n={int(g.n_gtex.max())} median_rho={g.rho_gtex_raw.median():.3f} "
                  f"neg_sig={len(negsig)}/{len(g)} sign_flip_vs_tumor={int(g.sign_flip_vs_tumor.sum())}")
    return {"contrast": contrast, "absratio": absratio, "coupling": coupling,
            "gtex_abundance": gtex_abund, "gtex_coupling": gtex_coupling, "xref": xref}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--top-orphans", type=int, default=25)
    ap.add_argument("--skip-gtex", action="store_true")
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, top_orphans=args.top_orphans, include_gtex=not args.skip_gtex)


if __name__ == "__main__":
    main()
