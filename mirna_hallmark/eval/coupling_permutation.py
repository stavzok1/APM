"""Permutation-based coupling inference at EVERY resolution.

Every coupling claim the framework makes -- edge, gene, target-set (miRNA module),
program (Hallmark), and universe -- is re-tested here against an empirical
**Freedman-Lane permutation null** instead of (only) the asymptotic partial-Spearman
p, and corrected with BOTH senses of "family":

- the *statistical* testing family: BH within each resolution, made dependence-robust
  with Benjamini-Yekutieli (``q_by``);
- the *biological* seed family: paralogue tests collapsed via Simes (``q_seedfamily``)
  and a per-seed-family min-statistic FWER (``p_seedfamily_fwer``), applied where the
  test unit is a miRNA (edge / target-set). Justified empirically in
  ``seed_family_justification``.

Predictor by resolution follows MODELING_FRAMEWORK §6.3 (abundance at edge/module,
aggregate pressure at gene/program). Confounder adjustment = CPE + HRD + TCGA batch
(the sample-level stack shared across all rows; the per-gene target-CN rung stays in
the parametric ladder). Outputs under ``output/coupling_permutation/``.

Run: ``.venv/bin/python3 -m mirna_hallmark.coupling_permutation [--n-perm N] [--smoke] [--gated]``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.coupling_inference import permutation_pvalues, seed_family_map
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark import pressure_build as PB
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.hallmark_interaction import (
    hallmark_pressure_matrix,
    hallmark_expression_matrix,
)

OUT_DIR = C.OUTPUT_ROOT / "coupling_permutation"


def _covariates(clinical: pd.DataFrame, samples: pd.Index) -> Optional[pd.DataFrame]:
    clin = clinical.set_index("participant")
    cols = [c for c in C.CONFOUNDER_NUMERIC if c in clin.columns]
    if not cols:
        return None
    cov = D.augment_tcga_batch(clin.reindex(samples)[cols])
    return cov.reindex(samples)


def _summ(res: pd.DataFrame, *, alpha: float, fam: bool) -> Dict[str, object]:
    neg = res["rho"] < 0
    out = {
        "n_tests": int(res["rho"].notna().sum()),
        "n_neg": int(neg.sum()),
        "n_perm_sig": int((neg & (res["p_perm"] < alpha)).sum()),
        "n_bh_sig": int((neg & (res["q_bh"] < alpha)).sum()),
        "n_by_sig": int((neg & (res["q_by"] < alpha)).sum()),
    }
    if fam and "q_seedfamily" in res.columns:
        sig_fam = neg & (res["q_seedfamily"] < alpha)
        out["n_edges_in_sig_seedfamily"] = int(sig_fam.sum())
        # the honest, de-duplicated discovery count: distinct seed families, not edges
        out["n_seedfamilies_total"] = int(res["family"].nunique())
        out["n_seedfamilies_sig"] = int(res.loc[sig_fam, "family"].nunique())
        out["n_seedfamily_fwer_sig"] = int((neg & (res["p_seedfamily_fwer"] < alpha)).sum())
    return out


def _run_resolution(
    name: str,
    predictor: pd.DataFrame,
    target: pd.DataFrame,
    cov: Optional[pd.DataFrame],
    *,
    seed_families=None,
    n_perm: int,
    alpha: float,
    out_dir: Path,
    return_null: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, object], Optional[np.ndarray]]:
    print(f"[perm] {name}: {predictor.shape[0]} tests x {predictor.shape[1]} samples, {n_perm} perms ...")
    r = permutation_pvalues(
        predictor, target, covariates=cov, families=seed_families,
        n_perm=n_perm, tail="neg", seed=0, return_null=return_null,
    )
    null_mat = None
    if return_null:
        r, null_mat = r
    if seed_families is not None:
        r = r.rename(columns={"q_family": "q_seedfamily", "p_family_fwer": "p_seedfamily_fwer"})
    r.to_csv(out_dir / f"coupling_{name}.tsv.gz", sep="\t", compression="gzip")
    return r, _summ(r, alpha=alpha, fam=seed_families is not None), null_mat


def run(
    *,
    out_dir: Path = OUT_DIR,
    n_perm: int = 2000,
    gated: bool = False,
    smoke: bool = False,
) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    alpha = C.FDR_ALPHA
    hs = HallmarkSets.load()
    genes = list(hs.universe)
    clinical = D.load_clinical_strata()  # primary tumours, PAM50 Normal-like already excluded
    mirna = D.load_mirna_arms()
    rna = D.load_rna()
    # dedup matrix indices so positional .loc gathers stay 1:1 with the edge list
    mirna = mirna[~mirna.index.duplicated(keep="first")]
    rna = rna[~rna.index.duplicated(keep="first")]
    # restrict to the framework cohort (tumour-only AND Normal-like dropped): the raw
    # miRNA/RNA matrices still carry the 36 Normal-like participants that clinical drops.
    cohort = set(clinical["participant"])
    mirna = mirna.loc[:, mirna.columns.isin(cohort)]
    rna = rna.loc[:, rna.columns.isin(cohort)]
    print(f"[perm] cohort: {mirna.shape[1]} miRNA x {rna.shape[1]} RNA participants "
          f"(primary tumour, Normal-like excluded)")

    # high-evidence spine edges, resolved onto the expression matrix, deduped to (arm,gene)
    edges = D.high_evidence_edges()
    edges = edges[["miRNA", "gene", "evidence_score"]].drop_duplicates()
    edges = PB.resolve_pressure_edges(edges, mirna)
    edges = edges[edges["gene"].isin(set(rna.index)) & edges["miRNA"].isin(set(mirna.index))]
    edges = edges.drop_duplicates(subset=["miRNA", "gene"])
    if smoke:  # cap work for a quick wiring check
        edges = edges.groupby("gene", group_keys=False).head(3).head(600)
        genes = sorted(set(edges["gene"]))

    # arm-name -> seed family (for the biological-family corrections)
    from mirna_hallmark.seed_family_justification import _arm_name_index
    arm_name = _arm_name_index(mirna)
    fam_of_name = seed_family_map(arm_name.unique())
    key_to_family = {k: fam_of_name[arm_name[k]] for k in mirna.index}

    summaries: Dict[str, object] = {}

    # ---------------- aggregate pressure (gene-level predictor) --------------
    gp = PB.compute_gene_pressure(genes, edges=edges if not smoke else None, mirna=mirna)
    gate = None
    if gated:
        gd = compute_ago_gate(rna)["ago_gate"]
        shared = gp.columns.intersection(gd.index)
        gp = gp[shared].mul(gd.reindex(shared), axis=1)

    # ---------------- EDGE: arm abundance vs gene expression -----------------
    eid = edges["miRNA"].astype(str) + "||" + edges["gene"].astype(str)
    Xe = mirna.loc[edges["miRNA"]].copy(); Xe.index = eid
    Ye = rna.loc[edges["gene"]].copy(); Ye.index = eid
    fam_e = [key_to_family.get(m, str(m)) for m in edges["miRNA"]]
    cov = _covariates(clinical, Xe.columns.intersection(Ye.columns))
    _, summaries["edge"], _ = _run_resolution(
        "edge", Xe, Ye, cov, seed_families=fam_e, n_perm=n_perm, alpha=alpha, out_dir=out_dir)

    # ---------------- GENE: aggregate pressure vs own expression -------------
    g_common = [g for g in gp.index if g in rna.index]
    Xg = gp.loc[g_common]; Yg = rna.loc[g_common]; Yg.index = Xg.index
    cov = _covariates(clinical, Xg.columns.intersection(Yg.columns))
    _, summaries["gene"], _ = _run_resolution(
        "gene", Xg, Yg, cov, seed_families=None, n_perm=n_perm, alpha=alpha, out_dir=out_dir)

    # ---------------- TARGET-SET / miRNA: arm abundance vs module composite ---
    arm_targets = edges.groupby("miRNA")["gene"].apply(lambda s: sorted(set(s)))
    rows_x, rows_y, mods = [], [], []
    z_rna = S.zscore_rows(rna)
    for arm, tgts in arm_targets.items():
        tg = [g for g in tgts if g in z_rna.index]
        if len(tg) < 3 or arm not in mirna.index:
            continue
        rows_x.append(mirna.loc[arm]); rows_y.append(z_rna.loc[tg].mean(axis=0)); mods.append(arm)
    Xm = pd.DataFrame(rows_x, index=mods); Ym = pd.DataFrame(rows_y, index=mods)
    fam_m = [key_to_family.get(m, str(m)) for m in mods]
    cov = _covariates(clinical, Xm.columns.intersection(Ym.columns))
    _, summaries["target_set"], _ = _run_resolution(
        "target_set", Xm, Ym, cov, seed_families=fam_m, n_perm=n_perm, alpha=alpha, out_dir=out_dir)

    # ---------------- PROGRAM: hallmark pressure vs hallmark expression -------
    hp = hallmark_pressure_matrix(gp, hs)
    he = hallmark_expression_matrix(rna, hs)
    common_h = hp.index.intersection(he.index)
    hp, he = hp.loc[common_h], he.loc[common_h]
    cov = _covariates(clinical, hp.columns.intersection(he.columns))
    res_prog, summaries["program"], null_prog = _run_resolution(
        "program", hp, he, cov, seed_families=None, n_perm=n_perm, alpha=alpha,
        out_dir=out_dir, return_null=True)

    # ---------------- UNIVERSE: count of coupled programs vs count-null -------
    obs_count = int(((res_prog["rho"] < 0) & (res_prog["p_perm"] < alpha)).sum())
    # per-perm count of programs whose permuted stat is as negative as each program's
    # own observed stat (a within-resolution count statistic on the joint null)
    obs_rho = res_prog["rho"].to_numpy()
    thr = np.nanquantile(null_prog, alpha, axis=0)  # per-program alpha-quantile of its null
    null_counts = np.nansum(null_prog <= thr[None, :], axis=1)  # ~alpha*n_programs under null
    universe = {
        "n_programs": int(res_prog["rho"].notna().sum()),
        "observed_coupled_programs_perm_sig": obs_count,
        "null_count_mean": round(float(np.mean(null_counts)), 3),
        "null_count_p95": float(np.nanpercentile(null_counts, 95)),
        "universe_p_perm": float((1.0 + np.sum(null_counts >= obs_count)) / (len(null_counts) + 1.0)),
    }
    summaries["universe"] = universe

    manifest = {
        "module": "mirna_hallmark.coupling_permutation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_perm": n_perm, "gated": gated, "smoke": smoke, "alpha": alpha,
        "confounders": list(C.CONFOUNDER_NUMERIC) + ["tcga_batch"],
        "predictor_by_resolution": {
            "edge": "arm log2(RPM)", "gene": "aggregate pressure", "target_set": "arm log2(RPM)",
            "program": "hallmark pressure", "universe": "count of coupled programs",
        },
        "family_senses": {
            "statistical_testing_family": "BH per resolution + Benjamini-Yekutieli (q_by)",
            "biological_seed_family": "Simes q_seedfamily + min-stat p_seedfamily_fwer (edge/target_set)",
        },
        "summaries": summaries,
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(summaries, indent=2))
    print(f"[perm] wrote outputs under {out_dir}")
    return summaries


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--gated", action="store_true", help="apply AGO/RISC gate to pressure predictors")
    ap.add_argument("--smoke", action="store_true", help="tiny subset + few perms for a wiring check")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, n_perm=(200 if args.smoke else args.n_perm),
        gated=args.gated, smoke=args.smoke)


if __name__ == "__main__":
    main()
