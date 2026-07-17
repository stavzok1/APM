"""Phase 1-refined (roadmap): does the CORE framework's headline coupling survive composition adjustment
with the REAL BayesPrism cell-type fractions (MH-70) — the gold-standard upgrade of MH-57 (metagenes)?

Same design as `core_coupling_composition_retest` (MH-57a), but the composition covariates are the
deconvolution fractions (Cancer/Normal Epithelial, CAFs, T/B/Myeloid/Plasmablasts, Endothelial, PVL)
instead of epi/immune/stroma metagene proxies. Baseline CPE+HRD+batch; test adds the 8 fractions (one
dropped for collinearity) + proliferation. Partial Spearman per headline BY-neg edge; report how many
stay negative — vs MH-57's 93% (metagene). apm-rigor-protocol B.

Run: ``.venv/bin/python3 -m mirna_hallmark.core_coupling_deconv_retest``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.utils.common.loaders import partial_spearman
from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark.analyses.cross_state.cross_state_coupling import _prolif_metagene
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr

EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
FRACTIONS = C.OUTPUT_ROOT / "brca_deconvolution" / "tcga_bayesprism_fractions.tsv"
OUT_DIR = C.OUTPUT_ROOT / "core_coupling_deconv_retest"


def run(*, out_dir: Path = OUT_DIR, max_edges: int = 1500, fractions: Path = FRACTIONS) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    rna = D.load_rna(); mirna = D.load_mirna_arms(); clin = D.load_clinical_strata()
    rna = rna[~rna.index.duplicated(keep="first")]; mirna = mirna[~mirna.index.duplicated(keep="first")]
    clin = clin.set_index("participant") if "participant" in clin.columns else clin

    frac = pd.read_csv(fractions, sep="\t", index_col=0)
    frac = frac.loc[[s for s in frac.index if not str(s).endswith("-NAT")]]      # tumours only
    frac = frac.drop(columns=[frac.mean().idxmin()])                            # drop smallest -> avoid sum=1 collinearity
    prolif = _prolif_metagene(rna, HallmarkSets.load())
    print(f"[deconv-retest] BayesPrism fractions: {frac.shape[1]} cell-type covariates, {frac.shape[0]} tumours")

    base = clin.reindex(frac.index)[[c for c in ("CPE", "thornsson_hrd_score") if c in clin.columns]].apply(pd.to_numeric, errors="coerce")
    try:
        from mirna_hallmark.tcga_batch import tcga_batch_dummies
        bd = tcga_batch_dummies(list(frac.index), kind=C.TCGA_BATCH_KIND).reindex(frac.index)
        base = base.join(bd)
        print(f"[deconv-retest] baseline = CPE+HRD + {bd.shape[1]} batch dummies")
    except Exception as ex:
        print(f"[deconv-retest] WARN no batch ({ex})")
    cov_base = base
    cov_deconv = base.join(frac).join(prolif.rename("prolif"))

    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    neg = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()
    neg["arm"] = neg["key"].str.split(r"\|\|").str[0]; neg["gene"] = neg["key"].str.split(r"\|\|").str[1]
    neg = neg[neg["arm"].isin(mirna.index) & neg["gene"].isin(rna.index)].head(max_edges)
    print(f"[deconv-retest] testing {len(neg)} headline BY-neg edges")

    rows = []
    for r in neg.itertuples(index=False):
        y, x = rna.loc[r.gene], mirna.loc[r.arm]
        rb, pb, _ = partial_spearman(y, x, cov_base)
        rc, pc, _ = partial_spearman(y, x, cov_deconv)
        rows.append({"arm": r.arm, "gene": r.gene, "rho_framework": r.rho, "rho_base": rb, "rho_deconv": rc, "p_deconv": pc})
    res = pd.DataFrame(rows)
    ok = res["p_deconv"].notna(); res["q_deconv"] = np.nan
    res.loc[ok, "q_deconv"] = bh_fdr(res.loc[ok, "p_deconv"].to_numpy())
    res.to_csv(out_dir / "headline_edges_deconv_retest.tsv", sep="\t", index=False)

    bneg = res["rho_base"] < 0
    summary = {
        "module": "mirna_hallmark.core_coupling_deconv_retest",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "composition_covariates": "BayesPrism cell-type fractions (8) + proliferation",
        "n_edges_tested": int(res["rho_base"].notna().sum()),
        "median_rho_base": round(float(res["rho_base"].median()), 4),
        "median_rho_deconv": round(float(res["rho_deconv"].median()), 4),
        "frac_kept_negative_after_deconv": round(float((res.loc[bneg, "rho_deconv"] < 0).mean()), 3),
        "n_base_neg": int(bneg.sum()),
        "n_still_FDR_neg_after_deconv": int(((res["rho_deconv"] < 0) & (res["q_deconv"] < C.FDR_ALPHA)).sum()),
        "median_attenuation": round(float(res.loc[bneg, "rho_base"].median() - res.loc[bneg, "rho_deconv"].median()), 4),
        "most_attenuated": [{"edge": f"{x.arm}->{x.gene}", "base": round(x.rho_base, 3), "deconv": round(x.rho_deconv, 3)}
                            for x in res.assign(d=res.rho_deconv - res.rho_base).nlargest(10, "d").itertuples()],
        "vs_MH57_metagene": "MH-57a (metagene): 93% kept negative, median rho -0.148 -> -0.124",
        "caveats": ["BayesPrism deconvolution fractions (MH-70) as composition covariates — gold-standard vs metagene",
                    "baseline CPE+HRD+batch; partial Spearman; tumour samples only"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[deconv-retest] median rho base {summary['median_rho_base']} -> deconv {summary['median_rho_deconv']} "
          f"(attenuation {summary['median_attenuation']})")
    print(f"[deconv-retest] kept negative after deconvolution: {summary['frac_kept_negative_after_deconv']:.0%} of {summary['n_base_neg']} "
          f"(MH-57 metagene was 93%); still FDR-neg {summary['n_still_FDR_neg_after_deconv']}")
    print(f"[deconv-retest] most attenuated: {[m['edge'] for m in summary['most_attenuated'][:6]]}")
    print(f"[deconv-retest] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--max-edges", type=int, default=1500)
    ap.add_argument("--fractions", type=Path, default=FRACTIONS)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, max_edges=args.max_edges, fractions=args.fractions)


if __name__ == "__main__":
    main()
