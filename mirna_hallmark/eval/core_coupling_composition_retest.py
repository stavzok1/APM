"""Phase 1 (roadmap): does the CORE framework's headline coupling survive composition adjustment?

The framework's tumour couplings (MH-1..47) adjust for CPE+HRD+batch but NOT cell composition — the
exact gap that flipped miR-29c (MH-54). This re-tests the headline **BY-significant negative edges**
(`coupling_permutation/coupling_edge`) on the TCGA spine: baseline partial Spearman (CPE+HRD) vs
+composition (epi/immune/stroma metagenes + proliferation), and asks which survive. This is the
primary-tumour version of the DCIS composition lesson, with the metagene route (a fuller scRNA-
deconvolution version is roadmap Phase 0/1-refined). `apm-rigor-protocol` B/C.

Run: ``.venv/bin/python3 -m mirna_hallmark.core_coupling_composition_retest --max-edges 250``
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
from mirna_hallmark.analyses.cross_state.cross_state_coupling import (
    EPI_MARKERS, IMMUNE_MARKERS, STROMA_MARKERS, _metagene, _prolif_metagene,
)
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.stats import bh_fdr

EDGES = C.OUTPUT_ROOT / "coupling_permutation" / "coupling_edge.tsv.gz"
OUT_DIR = C.OUTPUT_ROOT / "core_coupling_composition_retest"


def run(*, out_dir: Path = OUT_DIR, max_edges: int = 250) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    rna = D.load_rna(); mirna = D.load_mirna_arms(); clin = D.load_clinical_strata()
    rna = rna[~rna.index.duplicated(keep="first")]            # dup gene symbols -> .loc would return a frame
    mirna = mirna[~mirna.index.duplicated(keep="first")]
    clin = clin.set_index("participant") if "participant" in clin.columns else clin
    print(f"[retest] rna {rna.shape}, mirna {mirna.shape}, clinical {clin.shape}")

    # composition + proliferation covariates from TCGA RNA (same markers as the DCIS work)
    cov = pd.DataFrame({"epi": _metagene(rna, EPI_MARKERS), "immune": _metagene(rna, IMMUNE_MARKERS),
                        "stroma": _metagene(rna, STROMA_MARKERS), "prolif": _prolif_metagene(rna, HallmarkSets.load())})
    base = clin.reindex(cov.index)[[c for c in ("CPE", "thornsson_hrd_score") if c in clin.columns]].apply(pd.to_numeric, errors="coerce")
    # maximal-depth baseline = framework's CPE+HRD+sequencing-batch (plate_both); test adds composition
    try:
        from mirna_hallmark.tcga_batch import tcga_batch_dummies
        bd = tcga_batch_dummies(list(cov.index), kind=C.TCGA_BATCH_KIND).reindex(cov.index)
        base = base.join(bd)
        print(f"[retest] baseline = CPE+HRD + {bd.shape[1]} batch dummies ({C.TCGA_BATCH_KIND})")
    except Exception as ex:
        print(f"[retest] WARN no batch ({ex}); baseline CPE+HRD only")
    cov_base = base
    cov_comp = base.join(cov)

    e = pd.read_csv(EDGES, sep="\t").rename(columns={"Unnamed: 0": "key"})
    neg = e[(e["rho"] < 0) & (e["q_by"] < 0.05)].copy()
    neg["arm"] = neg["key"].str.split(r"\|\|").str[0]; neg["gene"] = neg["key"].str.split(r"\|\|").str[1]
    neg = neg[neg["arm"].isin(mirna.index) & neg["gene"].isin(rna.index)].head(max_edges)
    print(f"[retest] testing {len(neg)} headline BY-neg edges (cap {max_edges})")

    rows = []
    for r in neg.itertuples(index=False):
        y, x = rna.loc[r.gene], mirna.loc[r.arm]
        rb, pb, _ = partial_spearman(y, x, cov_base)
        rc, pc, _ = partial_spearman(y, x, cov_comp)
        rows.append({"arm": r.arm, "gene": r.gene, "rho_framework": r.rho,
                     "rho_base": rb, "rho_comp": rc, "p_comp": pc})
    res = pd.DataFrame(rows)
    ok = res["p_comp"].notna()
    res["q_comp"] = np.nan
    res.loc[ok, "q_comp"] = bh_fdr(res.loc[ok, "p_comp"].to_numpy())
    res.to_csv(out_dir / "headline_edges_composition_retest.tsv", sep="\t", index=False)

    base_neg = res["rho_base"] < 0
    summary = {
        "module": "mirna_hallmark.core_coupling_composition_retest",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_edges_tested": int(res["rho_base"].notna().sum()),
        "covariates_base": list(cov_base.columns), "covariates_added": list(cov.columns),
        "median_rho_base": round(float(res["rho_base"].median()), 4),
        "median_rho_comp": round(float(res["rho_comp"].median()), 4),
        "frac_kept_negative_after_composition": round(float((res.loc[base_neg, "rho_comp"] < 0).mean()), 3),
        "n_base_neg": int(base_neg.sum()),
        "n_still_FDR_neg_after_comp": int(((res["rho_comp"] < 0) & (res["q_comp"] < C.FDR_ALPHA)).sum()),
        "median_attenuation": round(float(res.loc[base_neg, "rho_base"].median() - res.loc[base_neg, "rho_comp"].median()), 4),
        "most_attenuated": [{"edge": f"{x.arm}->{x.gene}", "base": round(x.rho_base, 3), "comp": round(x.rho_comp, 3)}
                            for x in res.assign(d=res.rho_comp - res.rho_base).nlargest(8, "d").itertuples()],
        "caveats": ["metagene composition route (scRNA-deconvolution is the refined Phase-0/1 version)",
                    "capped batch of headline edges; baseline CPE+HRD (framework also uses batch)",
                    "bulk TCGA; partial Spearman; the primary-tumour analogue of MH-49/54"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[retest] median rho base {summary['median_rho_base']} -> comp {summary['median_rho_comp']} "
          f"(attenuation {summary['median_attenuation']})")
    print(f"[retest] kept negative after composition: {summary['frac_kept_negative_after_composition']:.0%} of {summary['n_base_neg']}; "
          f"still FDR-neg {summary['n_still_FDR_neg_after_comp']}")
    print(f"[retest] most attenuated: {[m['edge'] for m in summary['most_attenuated'][:5]]}")
    print(f"[retest] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--max-edges", type=int, default=250)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, max_edges=args.max_edges)


if __name__ == "__main__":
    main()
