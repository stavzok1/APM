#!/usr/bin/env python3
"""Within-cohort / per-PAM50 OLS explainability for hub-route target genes.

Nested linear models on log2(TPM+1) quantify how much variance route miRNA
(arm / family / HE pressure) explains **beyond** CN, CPE, proliferation,
promoter methylation, and (optional) SV hits.

Mirrors APM: ``cn_adjusted_rna_states``, ``per_gene_explainability``,
``effector_composition_score._incremental_r2``.

Scopes:
* ``cohort`` — all primary tumors; PAM50 dummies in base model.
* ``pam50_<subtype>`` — fit **within** LumA / LumB / Her2 / Basal (no PAM50 dummies).

Also fits an **RNA-low** subset (bottom quartile of expression within scope).

Outputs under ``mirna_hallmark/output/gene_expression_explainability/``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.hallmark_sets import load_hallmark_sets
from mirna_hallmark.robustness_checks import BASAL_MIRNA_FAMILIES, HUB_ROUTES

from mirna_hallmark.mirna_arm_resolve import resolve_arm
from mirna_hallmark.pressure_build import compute_gene_pressure, load_mirtar_edges  # noqa: E402
from analysis.utils.sample_filters import add_purity_filter_argument, apply_purity_filter  # noqa: E402
from mirna_hallmark.analyses.misc.hub_gene_methylation import load_hub_gene_methylation_matrix  # noqa: E402
from pipeline.config import PATHS  # noqa: E402

OUTPUT_DIR = C.OUTPUT_ROOT / "gene_expression_explainability"
METH_DIR = PATHS.methylation_output_dir / "per_sample"
SV_BURDEN = C.REPO_ROOT / "analysis/output/sv_per_sample/per_sample_sv_gene_burden.tsv.gz"
PAM50_FIT_SUBTYPES = ("LumA", "LumB", "Her2", "Basal")

# IRF1 route variants: curated APM arms, TargetScan-nominated alts, and their union.
IRF1_CURATED_ARMS = ["hsa-miR-23a-3p", "hsa-miR-106b-5p"]
IRF1_TARGETSCAN_ARMS = ["hsa-miR-130b-3p", "hsa-miR-301a-3p", "hsa-miR-301b-3p"]
IRF1_COMBINED_ARMS = list(dict.fromkeys([*IRF1_CURATED_ARMS, *IRF1_TARGETSCAN_ARMS]))

GENE_ROUTE_VARIANTS: Dict[str, Dict[str, List[str]]] = {
    "IRF1": {
        "curated": IRF1_CURATED_ARMS,
        "targetscan_alt": IRF1_TARGETSCAN_ARMS,
        "curated_plus_targetscan": IRF1_COMBINED_ARMS,
    },
}


def _route_sets_for_gene(gene: str) -> Dict[str, List[str]]:
    if gene in GENE_ROUTE_VARIANTS:
        return GENE_ROUTE_VARIANTS[gene]
    if gene in HUB_ROUTES:
        return {"hub": list(HUB_ROUTES[gene])}
    return {}


def _min_samples(scope_label: str, *, rna_low: bool = False) -> int:
    if scope_label == "cohort":
        return 25
    if rna_low:
        return 15
    return 30


def _align_series(s: pd.Series, index: pd.Index) -> pd.Series:
    if s.empty:
        return pd.Series(np.nan, index=index)
    if s.index.duplicated().any():
        s = s.groupby(level=0).mean()
    return s.reindex(index)


def _matrix_row_as_series(mat: pd.DataFrame, key: str, index: pd.Index) -> pd.Series:
    if key not in mat.index:
        return pd.Series(np.nan, index=index)
    row = mat.loc[key]
    if isinstance(row, pd.DataFrame):
        row = row.mean(axis=0)
    return _align_series(row, index)


def _pam50_dummies(pam: pd.Series) -> pd.DataFrame:
    s = pam.astype(str).replace("nan", "NA")
    return pd.get_dummies(s, prefix="pam50", drop_first=True, dtype=float)


def _family_groups_for_mirnas(mirnas: Sequence[str]) -> Dict[str, List[str]]:
    mir_set = set(mirnas)
    groups: Dict[str, List[str]] = {}
    for fam, arms in BASAL_MIRNA_FAMILIES.items():
        hit = [m for m in arms if m in mir_set]
        if hit:
            groups[fam] = hit
    assigned = {m for xs in groups.values() for m in xs}
    other = [m for m in mirnas if m not in assigned]
    if other:
        groups["other_hub_arms"] = other
    return groups


def _proliferation_metagene(rna: pd.DataFrame) -> pd.Series:
    hs = load_hallmark_sets()
    e2f = hs.get("HALLMARK_E2F_TARGETS", ())
    g2m = hs.get("HALLMARK_G2M_CHECKPOINT", ())
    genes = [g for g in (*e2f, *g2m) if g in rna.index]
    if len(genes) < 10:
        return pd.Series(dtype=float)
    return S.zscore_rows(rna.loc[genes]).mean(axis=0)


def _load_sv_hits(genes: Sequence[str]) -> Optional[pd.DataFrame]:
    if not SV_BURDEN.exists():
        print(f"[explain] SV burden missing ({SV_BURDEN}); SV layer skipped.")
        return None
    keep = ["participant", "gene", "n_with_regulatory_hit", "n_coding_disruptive", "n_bp_promoter"]
    sv = pd.read_csv(SV_BURDEN, sep="\t", usecols=lambda c: c in keep, low_memory=False)
    for col in ["n_with_regulatory_hit", "n_coding_disruptive", "n_bp_promoter"]:
        if col not in sv.columns:
            sv[col] = 0
    sv = sv.loc[sv["gene"].isin(set(genes))].copy()
    sv["sv_hit"] = (
        (pd.to_numeric(sv["n_with_regulatory_hit"], errors="coerce").fillna(0) > 0)
        | (pd.to_numeric(sv["n_coding_disruptive"], errors="coerce").fillna(0) > 0)
        | (pd.to_numeric(sv["n_bp_promoter"], errors="coerce").fillna(0) > 0)
    ).astype(float)
    sv = sv.rename(columns={"gene": "gene_name"})
    return sv.groupby(["participant", "gene_name"], as_index=False)["sv_hit"].max()


def _build_design(
    df: pd.DataFrame,
    numeric_cols: Sequence[str],
    *,
    include_pam50: bool,
) -> pd.DataFrame:
    parts: List[pd.DataFrame] = [pd.DataFrame({"intercept": 1.0}, index=df.index)]
    for col in numeric_cols:
        if col not in df.columns:
            continue
        parts.append(pd.DataFrame({col: pd.to_numeric(df[col], errors="coerce")}))
    if include_pam50 and "PAM50_final" in df.columns:
        parts.append(_pam50_dummies(df["PAM50_final"]))
    X = pd.concat(parts, axis=1)
    X = X.replace([np.inf, -np.inf], np.nan)
    med = X.median(numeric_only=True)
    X = X.fillna(med)
    keep = [c for c in X.columns if c == "intercept" or X[c].std(ddof=0) > 1e-12]
    return X[keep]


def _ols_fit(
    y: pd.Series,
    X: pd.DataFrame,
    *,
    min_n: int,
) -> Tuple[float, int, Dict[str, float], pd.Series]:
    dat = pd.concat([y.rename("_y"), X], axis=1).dropna()
    pred = pd.Series(np.nan, index=y.index, dtype=float)
    if len(dat) < min_n or len(dat) < X.shape[1] + 5:
        return np.nan, len(dat), {}, pred
    yv = dat["_y"].to_numpy(dtype=float)
    Xm = dat[X.columns].to_numpy(dtype=float)
    try:
        beta, *_ = np.linalg.lstsq(Xm, yv, rcond=1e-8)
    except np.linalg.LinAlgError:
        return np.nan, len(dat), {}, pred
    sst = float(((yv - yv.mean()) ** 2).sum())
    if sst <= 0:
        return np.nan, len(dat), {}, pred
    yhat = Xm @ beta
    pred.loc[dat.index] = yhat
    r2 = float(1.0 - ((yv - yhat) ** 2).sum() / sst)
    coefs = {name: float(b) for name, b in zip(X.columns, beta)}
    return r2, len(dat), coefs, pred


def _prepare_gene_frame(
    gene: str,
    *,
    participants: Sequence[str],
    rna: pd.DataFrame,
    cnv: pd.DataFrame,
    clinical: pd.DataFrame,
    meth: pd.DataFrame,
    mirna: pd.DataFrame,
    route_pressure: pd.Series,
    route_arms: Sequence[str],
    proliferation: Optional[pd.Series],
    sv_hits: Optional[pd.DataFrame],
) -> pd.DataFrame:
    parts = pd.Index(participants)
    out = pd.DataFrame(index=parts)
    out["rna_log2tpm"] = _matrix_row_as_series(rna, gene, parts)
    out["copy_number"] = _matrix_row_as_series(cnv, gene, parts)
    clin = clinical.drop_duplicates("participant").set_index("participant").reindex(parts)
    out["CPE"] = pd.to_numeric(clin["CPE"], errors="coerce")
    out["PAM50_final"] = clin["PAM50_final"]
    out["promoter_beta_mean"] = _matrix_row_as_series(meth, gene, parts)
    out["route_he_pressure"] = _align_series(route_pressure, parts)
    if proliferation is not None and not proliferation.empty:
        out["proliferation_z"] = _align_series(proliferation, parts)
    if sv_hits is not None and not sv_hits.empty:
        sv_g = sv_hits.loc[sv_hits["gene_name"].eq(gene)].set_index("participant")["sv_hit"]
        out["sv_hit"] = _align_series(sv_g, parts)
    else:
        out["sv_hit"] = np.nan
    for arm in route_arms:
        resolved, _ = resolve_arm(arm, mirna.index)
        if resolved in mirna.index:
            out[f"arm:{arm}"] = _align_series(mirna.loc[resolved], parts)
    for fam, fam_arms in _family_groups_for_mirnas(route_arms).items():
        present = []
        for a in fam_arms:
            resolved, _ = resolve_arm(a, mirna.index)
            if resolved in mirna.index:
                present.append(resolved)
        if present:
            short = fam.replace("miR-17~92/106b~25/106a~363", "miR-17~92")
            fam_df = pd.DataFrame({a: _align_series(mirna.loc[a], parts) for a in present})
            out[f"family:{short}"] = fam_df.mean(axis=1)
    return out


def _rna_low_mask(y: pd.Series) -> pd.Series:
    v = pd.to_numeric(y, errors="coerce")
    ok = v.dropna()
    if len(ok) < 8:
        return pd.Series(False, index=y.index)
    q1 = ok.quantile(0.25)
    return v.notna() & v.lt(q1)


def explain_one_gene(
    gene: str,
    *,
    frame: pd.DataFrame,
    route_set: str,
    include_proliferation: bool,
    include_sv: bool,
    include_pam50: bool,
    scope_label: str,
    fit_set: str = "all",
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return (incremental_r2, coefficients, residuals)."""
    if fit_set == "rna_low":
        sub = frame.loc[_rna_low_mask(frame["rna_log2tpm"])].copy()
    else:
        sub = frame
    min_n = _min_samples(scope_label, rna_low=(fit_set == "rna_low"))

    base_num = ["copy_number", "CPE"]
    if include_proliferation and "proliferation_z" in sub.columns:
        base_num.append("proliferation_z")
    meth_num = [*base_num, "promoter_beta_mean"]
    sv_num = [*meth_num, "sv_hit"] if include_sv and sub["sv_hit"].notna().any() else meth_num
    family_cols = [c for c in sub.columns if c.startswith("family:")]
    arm_cols = [c for c in sub.columns if c.startswith("arm:")]

    specs: List[Tuple[str, List[str]]] = [
        ("base_cn_cpe", base_num),
        ("base_plus_meth", meth_num),
    ]
    if sv_num != meth_num:
        specs.append(("base_meth_plus_sv", sv_num))
    core = sv_num
    specs.extend(
        [
            ("core_plus_route_pressure", [*core, "route_he_pressure"]),
            ("core_plus_families", [*core, *family_cols]),
            ("core_plus_route_arms", [*core, *arm_cols]),
            ("full_route_block", [*core, "route_he_pressure", *family_cols, *arm_cols]),
        ]
    )

    y = pd.to_numeric(sub["rna_log2tpm"], errors="coerce")
    r2_rows: List[dict] = []
    coef_rows: List[dict] = []
    resid_rows: List[dict] = []
    r2_by_model: Dict[str, float] = {}
    preds: Dict[str, pd.Series] = {}
    core_model_name = "base_meth_plus_sv" if include_sv and sub["sv_hit"].notna().any() else "base_plus_meth"

    for model_name, num_cols in specs:
        if model_name.endswith("_arms") and not arm_cols:
            continue
        if model_name.endswith("_families") and not family_cols:
            continue
        if "route_pressure" in model_name and sub["route_he_pressure"].notna().sum() < min_n:
            continue
        if "sv" in model_name and "sv_hit" not in num_cols:
            continue
        X = _build_design(sub, num_cols, include_pam50=include_pam50)
        r2, n, coefs, pred = _ols_fit(y, X, min_n=min_n)
        r2_by_model[model_name] = r2
        preds[model_name] = pred
        base_key = "base_cn_cpe"
        meth_key = "base_plus_meth"
        delta_base = r2 - r2_by_model.get(base_key, np.nan) if model_name != base_key else 0.0
        delta_meth = (
            r2 - r2_by_model.get("base_plus_meth", np.nan)
            if model_name not in {"base_cn_cpe", "base_plus_meth"}
            else 0.0
        )
        r_core = r2_by_model.get(core_model_name, np.nan)
        delta_core = (
            r2 - r_core
            if model_name.startswith("core_plus") or model_name == "full_route_block"
            else 0.0
        )
        r2_rows.append(
            {
                "gene": gene,
                "route_set": route_set,
                "scope": scope_label,
                "fit_set": fit_set,
                "model": model_name,
                "r2": round(r2, 4) if pd.notna(r2) else np.nan,
                "delta_r2_vs_base": round(delta_base, 4) if pd.notna(delta_base) else np.nan,
                "delta_r2_vs_meth": round(delta_meth, 4) if pd.notna(delta_meth) else np.nan,
                "delta_r2_vs_core": round(delta_core, 4) if pd.notna(delta_core) else np.nan,
                "n": int(n),
                "n_meth_assayed": int(sub["promoter_beta_mean"].notna().sum()),
                "n_sv_assayed": int(sub["sv_hit"].notna().sum()) if include_sv else 0,
                "n_route_arms": len(arm_cols),
                "include_proliferation": include_proliferation,
                "include_sv": include_sv,
            }
        )
        for term, beta in coefs.items():
            if term in {"intercept", "copy_number", "CPE", "promoter_beta_mean", "proliferation_z", "sv_hit"}:
                continue
            if not term.startswith(("arm:", "family:", "route_he_pressure", "pam50_")):
                continue
            coef_rows.append(
                {
                    "gene": gene,
                    "route_set": route_set,
                    "scope": scope_label,
                    "fit_set": fit_set,
                    "model": model_name,
                    "term": term,
                    "beta": round(beta, 5),
                    "n": int(n),
                }
            )

    meth_pred = preds.get("base_plus_meth", pd.Series(dtype=float))
    arm_pred = preds.get("core_plus_route_arms", pd.Series(dtype=float))
    for pid in sub.index:
        yv = y.get(pid, np.nan)
        if pd.isna(yv):
            continue
        resid_rows.append(
            {
                "participant": pid,
                "gene": gene,
                "route_set": route_set,
                "scope": scope_label,
                "fit_set": fit_set,
                "PAM50_final": sub.loc[pid, "PAM50_final"] if "PAM50_final" in sub.columns else np.nan,
                "rna_log2tpm": float(yv),
                "yhat_base_meth": float(meth_pred.get(pid, np.nan)),
                "yhat_core_arms": float(arm_pred.get(pid, np.nan)),
                "residual_base_meth": float(yv - meth_pred.get(pid, np.nan))
                if pd.notna(meth_pred.get(pid, np.nan))
                else np.nan,
                "residual_core_arms": float(yv - arm_pred.get(pid, np.nan))
                if pd.notna(arm_pred.get(pid, np.nan))
                else np.nan,
            }
        )

    return pd.DataFrame(r2_rows), pd.DataFrame(coef_rows), pd.DataFrame(resid_rows)


def _scope_out_dir(out_dir: Path, scope_label: str) -> Path:
    if scope_label == "cohort":
        return out_dir / "cohort"
    return out_dir / f"pam50_{scope_label}"


def run_explainability_scope(
    *,
    genes: Sequence[str],
    out_dir: Path,
    clinical: pd.DataFrame,
    rna: pd.DataFrame,
    mirna: pd.DataFrame,
    cnv: pd.DataFrame,
    meth: pd.DataFrame,
    edges: pd.DataFrame,
    sv_hits: Optional[pd.DataFrame],
    proliferation: Optional[pd.Series],
    scope_label: str,
    include_proliferation: bool,
    include_sv: bool,
    fit_rna_low: bool,
) -> None:
    scoped_out = _scope_out_dir(out_dir, scope_label)
    scoped_out.mkdir(parents=True, exist_ok=True)
    include_pam50 = scope_label == "cohort"

    participants = D.common_participants(clinical["participant"], rna.columns, mirna.columns)
    if scope_label != "cohort":
        clin_idx = clinical.set_index("participant")
        participants = [
            p
            for p in participants
            if p in clin_idx.index and clin_idx.loc[p, "PAM50_final"] == scope_label
        ]

    all_r2: List[pd.DataFrame] = []
    all_coef: List[pd.DataFrame] = []
    all_resid: List[pd.DataFrame] = []

    for gene in genes:
        route_sets = _route_sets_for_gene(gene)
        if not route_sets:
            continue
        for route_set, route_arms in route_sets.items():
            present_arms = []
            for m in route_arms:
                resolved, _ = resolve_arm(m, mirna.index)
                if resolved in mirna.index and resolved not in present_arms:
                    present_arms.append(resolved)
            if not present_arms:
                print(f"[explain] skip {gene}/{route_set}: no arms in expression matrix")
                continue
            route_edges = edges.loc[
                edges["gene"].eq(gene) & edges["miRNA"].isin(present_arms)
            ].drop_duplicates(["miRNA", "gene"])
            if route_edges.empty:
                pressure = pd.Series(0.0, index=participants)
            else:
                pressure = compute_gene_pressure([gene], edges=route_edges, mirna=mirna).loc[gene]

            frame = _prepare_gene_frame(
                gene,
                participants=participants,
                rna=rna,
                cnv=cnv,
                clinical=clinical,
                meth=meth,
                mirna=mirna,
                route_pressure=pressure,
                route_arms=present_arms,
                proliferation=proliferation,
                sv_hits=sv_hits,
            )
            for fit_set in (["all", "rna_low"] if fit_rna_low else ["all"]):
                r2_df, coef_df, resid_df = explain_one_gene(
                    gene,
                    frame=frame,
                    route_set=route_set,
                    include_proliferation=include_proliferation,
                    include_sv=include_sv,
                    include_pam50=include_pam50,
                    scope_label=scope_label,
                    fit_set=fit_set,
                )
                if r2_df.empty:
                    continue
                all_r2.append(r2_df)
                all_coef.append(coef_df)
                all_resid.append(resid_df)
                arms = r2_df.loc[r2_df["model"].eq("core_plus_route_arms")]
                if not arms.empty:
                    row = arms.iloc[0]
                    print(
                        f"[explain] {gene}/{route_set} ({scope_label}, {fit_set}): "
                        f"ΔR²(arms|core)={row['delta_r2_vs_core']:.3f} "
                        f"R²={row['r2']:.3f} n={row['n']}"
                    )

    r2_out = pd.concat(all_r2, ignore_index=True) if all_r2 else pd.DataFrame()
    coef_out = pd.concat(all_coef, ignore_index=True) if all_coef else pd.DataFrame()
    resid_out = pd.concat(all_resid, ignore_index=True) if all_resid else pd.DataFrame()
    r2_out.to_csv(scoped_out / "incremental_r2_hub_genes.tsv", sep="\t", index=False)
    coef_out.to_csv(scoped_out / "mirna_term_coefficients.tsv", sep="\t", index=False)
    if not resid_out.empty:
        resid_out.to_csv(scoped_out / "per_sample_residuals.tsv.gz", sep="\t", index=False, compression="gzip")

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "scope": scope_label,
        "include_pam50_in_base": include_pam50,
        "include_proliferation": include_proliferation,
        "include_sv": include_sv,
        "fit_rna_low": fit_rna_low,
        "genes": list(genes),
        "n_participants": len(participants),
    }
    (scoped_out / "method_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[explain] wrote {scoped_out}")


def build_subtype_summary(out_dir: Path) -> pd.DataFrame:
    """Pivot ΔR²(arms|core) across PAM50 subtypes for each gene/route_set."""
    rows: List[dict] = []
    for sub in PAM50_FIT_SUBTYPES:
        path = _scope_out_dir(out_dir, sub) / "incremental_r2_hub_genes.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        arms = df.loc[
            df["model"].eq("core_plus_route_arms") & df["fit_set"].eq("all")
        ]
        for _, r in arms.iterrows():
            rows.append(
                {
                    "gene": r["gene"],
                    "route_set": r["route_set"],
                    "subtype": sub,
                    "delta_r2_arms_vs_core": r["delta_r2_vs_core"],
                    "r2_arms_model": r["r2"],
                    "n": r["n"],
                }
            )
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary.to_csv(out_dir / "incremental_r2_by_subtype_summary.tsv", sep="\t", index=False)
    return summary


def run_explainability(
    *,
    genes: Sequence[str],
    out_dir: Path,
    scopes: Sequence[str],
    min_purity: Optional[float] = None,
    include_proliferation: bool = True,
    include_sv: bool = True,
    fit_rna_low: bool = True,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = sorted({g for g in genes if g in HUB_ROUTES or g in GENE_ROUTE_VARIANTS})
    if not genes:
        raise SystemExit("No hub-route genes to fit.")

    clinical = D.load_clinical_strata().drop_duplicates("participant")
    if min_purity is not None:
        clinical = apply_purity_filter(clinical, min_purity=min_purity)
    rna = D.load_rna()
    mirna = D.load_mirna_arms()
    cnv = D.load_cnv_target_genes(genes)
    meth = load_hub_gene_methylation_matrix(genes, per_sample_dir=METH_DIR)
    if meth.empty:
        print("[explain] hub methylation cache empty; meth layer inactive.")
    edges = load_mirtar_edges(genes, resolve_arms=True)
    sv_hits = _load_sv_hits(genes) if include_sv else None
    prolif = _proliferation_metagene(rna) if include_proliferation else None

    for scope in scopes:
        run_explainability_scope(
            genes=genes,
            out_dir=out_dir,
            clinical=clinical,
            rna=rna,
            mirna=mirna,
            cnv=cnv,
            meth=meth,
            edges=edges,
            sv_hits=sv_hits,
            proliferation=prolif,
            scope_label=scope,
            include_proliferation=include_proliferation,
            include_sv=include_sv,
            fit_rna_low=fit_rna_low,
        )

    if any(s in PAM50_FIT_SUBTYPES for s in scopes):
        build_subtype_summary(out_dir)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", nargs="*", default=None, help="Hub genes (default: all)")
    ap.add_argument("--out-dir", type=Path, default=OUTPUT_DIR)
    ap.add_argument(
        "--scope",
        choices=["cohort", "all-pam50", "all"],
        default="all",
        help="cohort | all four PAM50 subtypes | both (default)",
    )
    ap.add_argument("--pam50", default=None, help="Single PAM50 subtype only (overrides --scope)")
    ap.add_argument("--no-proliferation", action="store_true")
    ap.add_argument("--no-sv", action="store_true")
    ap.add_argument("--no-rna-low", action="store_true")
    add_purity_filter_argument(ap)
    args = ap.parse_args()

    genes = args.genes or sorted(set(HUB_ROUTES) | set(GENE_ROUTE_VARIANTS))
    if args.pam50:
        scopes = [args.pam50]
    elif args.scope == "cohort":
        scopes = ["cohort"]
    elif args.scope == "all-pam50":
        scopes = list(PAM50_FIT_SUBTYPES)
    else:
        scopes = ["cohort", *PAM50_FIT_SUBTYPES]

    run_explainability(
        genes=genes,
        out_dir=args.out_dir,
        scopes=scopes,
        min_purity=args.min_purity,
        include_proliferation=not args.no_proliferation,
        include_sv=not args.no_sv,
        fit_rna_low=not args.no_rna_low,
    )


if __name__ == "__main__":
    main()
