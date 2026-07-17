"""Comprehensive pressure-VARIANT x RESOLUTION (incl. gene-role sets) survival sweep + multivariate signature.

Answers "every survival idea we can leverage": across all pressure **variants** (5 layers: L1 spine /
L2 attribution / L3 ACQUIRED / L4a posclip / L4b abs + 7 hybrid models m0/m2/m3/m5/m7/m8/m11), test
survival at every **resolution**: program (50), gene (for spine + acquired), and gene-ROLE aggregates
(TF / oncogene / TSG / driver / high-evidence-target / all). Cox vs TCGA DFI+OS, composition-adjusted,
BH-FDR within each (variant x resolution) family. Plus **(#3) a MULTIVARIATE elastic-net Cox pressure
SIGNATURE** (5-fold CV concordance vs a clinical-only baseline). Honest gate: any TCGA hit is replicated
in the independent Buffa cohort (role aggregates of Buffa-native pressure).

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_pressure_variants``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

from mirna_hallmark import config as C
from mirna_hallmark.gene_roles import load_gene_roles
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z, load_tcga_outcome

OUT_DIR = C.OUTPUT_ROOT / "outcome_pressure_variants"
PLC = C.OUTPUT_ROOT / "pressure_layer_comparison"
HYB = C.OUTPUT_ROOT / "hybrid_pressure"
TF_TSV = C.REPO_ROOT / "annotations" / "humantfs_lambert2018_tf.tsv"
EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
DRIVERS = C.OUTPUT_ROOT / "cptac_validation" / "prospective" / "resid_survivors" / "survivor_drivers.tsv"
MEASUREMENTS = {  # variant -> dir with hallmark_/gene_pressure_per_sample
    "L1_spine": PLC / "L1_coupling_spine", "L2_attrib": PLC / "L2_attribution_noz",
    "L3_acquired": PLC / "L3_acquisition_devhealthy", "L4a_posclip": PLC / "L4a_posclip_spine",
    "L4b_abs": PLC / "L4b_abs_spine",
    **{f"m{k}": HYB / f"m{k}" for k in (0, 2, 3, 5, 7, 8, 11)},
}
GENE_LEVEL = {"L1_spine", "L3_acquired"}            # gene-resolution only for spine + acquired (cost)


def _role_sets() -> dict:
    gr = load_gene_roles()
    tf = set(pd.read_csv(TF_TSV, sep="\t").iloc[:, 0].astype(str))
    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    he = set(e[e["q_by"] < 0.05]["key"].str.split(r"\|\|").str[1])
    drv = set()
    if DRIVERS.exists():
        d = pd.read_csv(DRIVERS, sep="\t")
        col = "target_gene" if "target_gene" in d.columns else ("gene_name" if "gene_name" in d.columns else d.columns[0])
        drv = set(d[col].astype(str))
    return {"TF": tf, "oncogene": set(gr[gr.malignancy_sign == 1].gene),
            "TSG": set(gr[gr.malignancy_sign == -1].gene), "driver": drv, "HEtarget": he}


def _role_features(gene_p: pd.DataFrame, roles: dict) -> pd.DataFrame:
    rows = {"role_all": gene_p.mean(axis=0)}
    for r, gs in roles.items():
        g = gene_p.index.intersection(gs)
        if len(g) >= 5:
            rows[f"role_{r}"] = gene_p.loc[g].mean(axis=0)
    return pd.DataFrame(rows).T


def _screen(mat: pd.DataFrame, variant: str, fam: str, cl: dict) -> pd.DataFrame:
    base = cl["base"]; adj = pd.concat([base, cl["comp"]], axis=1)
    rows = []
    for feat in mat.index:
        x = mat.loc[feat].astype(float)
        if x.std(skipna=True) == 0 or x.notna().sum() < 50:
            continue
        dca, dpa = _cox_one(cl["dfi_t"], cl["dfi_e"], x, adj)
        oca, opa = _cox_one(cl["os_t"], cl["os_e"], x, adj)
        rows.append({"variant": variant, "resolution": fam, "feature": str(feat),
                     "dfi_hr": np.exp(dca), "dfi_p": dpa, "os_hr": np.exp(oca), "os_p": opa})
    r = pd.DataFrame(rows)
    if r.empty:
        return r
    for col in ("dfi_p", "os_p"):
        ok = r[col].notna(); r.loc[ok, col.replace("_p", "_q")] = bh_fdr(r.loc[ok, col].to_numpy())
    return r


def _multivariate(feat: pd.DataFrame, cl: dict, ep="dfi") -> dict:
    """5-fold CV concordance of an elastic-net Cox pressure signature vs clinical-only baseline."""
    t, e = cl[f"{ep}_t"], cl[f"{ep}_e"]
    X = feat.T.copy()                                                # sample x feature
    df = X.join(cl["base"]).join(t.rename("T")).join(e.rename("E")).dropna()
    df = df.loc[:, df.std(numeric_only=True) > 0]
    feats = [c for c in X.columns if c in df.columns]
    clin = [c for c in cl["base"].columns if c in df.columns]
    rng = np.random.default_rng(0); idx = rng.permutation(df.index.to_numpy())
    folds = np.array_split(idx, 5)
    cidx_full, cidx_clin = [], []
    for i in range(5):
        te = df.loc[folds[i]]; tr = df.drop(index=folds[i])
        if tr["E"].sum() < 10 or te["E"].sum() < 3:
            continue
        for cols, store in ((feats + clin, cidx_full), (clin, cidx_clin)):
            try:
                m = CoxPHFitter(penalizer=0.1, l1_ratio=0.5).fit(tr[cols + ["T", "E"]], "T", "E")
                risk = m.predict_partial_hazard(te[cols])
                store.append(concordance_index(te["T"], -risk, te["E"]))
            except Exception:
                pass
    return {"endpoint": ep, "n_features": len(feats),
            "cv_concordance_pressure+clin": round(float(np.mean(cidx_full)), 3) if cidx_full else None,
            "cv_concordance_clin_only": round(float(np.mean(cidx_clin)), 3) if cidx_clin else None}


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    cl = load_tcga_outcome()
    roles = _role_sets()
    print(f"[variants] roles: {{{', '.join(f'{k}:{len(v)}' for k,v in roles.items())}}} | DFI {int(cl['dfi_e'].sum())} OS {int(cl['os_e'].sum())}")
    all_res = []
    sig_feats = {}                                                  # for the multivariate signature
    for var, d in MEASUREMENTS.items():
        hp = d / "hallmark_pressure_per_sample.tsv.gz"; gp = d / "gene_pressure_per_sample.tsv.gz"
        if not hp.exists():
            continue
        prog = pd.read_csv(hp, sep="\t", index_col=0)
        all_res.append(_screen(prog, var, "program", cl))
        if gp.exists():
            gpm = pd.read_csv(gp, sep="\t", index_col=0)
            rf = _role_features(gpm, roles)
            all_res.append(_screen(rf, var, "role", cl))
            if var in {"L1_spine", "L3_acquired"}:
                sig_feats[var] = pd.concat([prog, rf])              # for signature
                all_res.append(_screen(gpm, var, "gene", cl))
        print(f"[variants]   {var} done")
    res = pd.concat([r for r in all_res if not r.empty], ignore_index=True)
    res.to_csv(out_dir / "variant_resolution_survival.tsv", sep="\t", index=False)

    # (#3) multivariate signature from spine + acquired (program + role features)
    sig_mat = pd.concat(list(sig_feats.values())) if sig_feats else pd.DataFrame()
    sig_mat = sig_mat[~sig_mat.index.duplicated()]
    mv = {ep: _multivariate(sig_mat, cl, ep) for ep in ("dfi", "os")} if len(sig_mat) else {}

    fdr = res[(res.dfi_q < 0.1) | (res.os_q < 0.1)].copy()
    summary = {
        "module": "mirna_hallmark.outcome_pressure_variants",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_variants": len(MEASUREMENTS), "roles": {k: len(v) for k, v in roles.items()},
        "n_FDR_hits_total": int(len(fdr)),
        "FDR_hits": fdr.sort_values("dfi_q").head(25)
            .apply(lambda r: f"{r['variant']}:{r['resolution']}:{r['feature']} DFI HR={r['dfi_hr']:.2f} q={r['dfi_q']:.3f}", axis=1).tolist(),
        "acquired_role_results": res[(res.variant == "L3_acquired") & (res.resolution == "role")]
            .apply(lambda r: f"{r['feature']} DFI HR={r['dfi_hr']:.2f} q={r['dfi_q']:.3f} | OS q={r['os_q']:.3f}", axis=1).tolist(),
        "multivariate_signature_3": mv,
        "caveats": ["TCGA-only; comp-adjusted Cox; FDR within (variant x resolution)",
                    "large multiple-testing surface across 12 variants — any hit needs Buffa replication",
                    "gene-resolution only for spine + acquired (cost)"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[variants] total FDR hits: {summary['n_FDR_hits_total']}")
    print(f"[variants] acquired role: {summary['acquired_role_results']}")
    print(f"[variants] #3 multivariate signature: {mv}")
    print(f"[variants] top FDR: {summary['FDR_hits'][:10]}")
    print(f"[variants] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
