"""Does the promiscuity denominator matter on the ATTRIBUTION axis?

The held-out tuning showed the denominator (``target_norm``) is coupling-neutral, so
it cannot be justified by "sharpest coupling". Its *design intent* is promiscuity
control on the attribution axis — but intent is not evidence. This sweep MEASURES the
denominator's effect on attribution so §3.3 can be justified (or dropped) by data:

For ``target_norm in {none, degree, evidence_mass (spine), combined_mass}`` it recomputes
the per-(gene,arm) shares in the **attribution mode** (``softmax_logrpm``, no-z, where
share/specificity actually live) and reports, relative to the spine and to bare ``none``:

- **dominant-arm concordance** — fraction of genes whose dominant regulator (argmax
  share) is unchanged. Low concordance ⇒ the denominator re-assigns credit.
- **specificity** — median ``spec`` and its rank-correlation across denominators.
- **hub concentration** — genes-dominated-per-arm: Gini, top-10-arm share, #arms to reach
  50% of all dominances. Promiscuity control should *lower* concentration vs ``none``
  (a promiscuous high-evidence arm dominating many genes is exactly what it penalizes).

Run: ``.venv/bin/python3 -m mirna_hallmark.denominator_attribution_sweep [--smoke]``.
Outputs under ``output/denominator_attribution_sweep/``.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import pressure_build as PB
from mirna_hallmark.hallmark_sets import HallmarkSets

OUT_DIR = C.OUTPUT_ROOT / "denominator_attribution_sweep"
DENOMS = ["none", "degree", "evidence_mass", "combined_mass"]  # evidence_mass = spine
SPINE_DENOM = "evidence_mass"
ATTR_MODE = "softmax_logrpm"  # no-z attribution mode (abs_share == signed_share)


def _gini(x: np.ndarray) -> float:
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    if n == 0 or x.sum() == 0:
        return np.nan
    return float((2.0 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


def _dominant_and_spec(contrib: pd.DataFrame) -> Dict[str, object]:
    """From a contributions frame, derive per-gene dominant arm + per-edge specificity."""
    c = contrib.copy()
    # dominant regulator per gene (max gene-level magnitude share)
    dom = c.loc[c.groupby("gene")["global_abs_share"].idxmax(), ["gene", "miRNA"]]
    dom = dom.set_index("gene")["miRNA"]
    # specificity: arm's mass on a gene / arm's total mass across its genes
    arm_mass = c.groupby("miRNA")["mean_abs_contribution"].transform("sum")
    c["spec"] = c["mean_abs_contribution"] / arm_mass.replace(0, np.nan)
    spec = c.set_index(["gene", "miRNA"])["spec"]
    # hub concentration: genes-dominated per arm
    dom_counts = dom.value_counts()
    return {"dominant": dom, "spec": spec, "dom_counts": dom_counts}


def run(*, out_dir: Path = OUT_DIR, smoke: bool = False) -> Dict[str, object]:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    genes = list(hs.universe)
    clinical = D.load_clinical_strata()
    mirna = D.load_mirna_arms()
    mirna = mirna[~mirna.index.duplicated(keep="first")]
    cohort = sorted(set(clinical["participant"]) & set(mirna.columns))
    mirna = mirna[cohort]
    if smoke:
        genes = genes[:400]
    print(f"[attr] {len(genes)} genes, {len(cohort)} participants, mode={ATTR_MODE}")

    per_denom: Dict[str, Dict[str, object]] = {}
    rows: List[dict] = []
    for tn in DENOMS:
        print(f"[attr] contributions under target_norm={tn} ...")
        contrib = PB.compute_gene_pressure_contributions(
            genes, expr_mode=ATTR_MODE, target_norm=tn, mirna=mirna,
        )
        if contrib.empty:
            continue
        per_denom[tn] = _dominant_and_spec(contrib)
        dc = per_denom[tn]["dom_counts"].to_numpy()
        rows.append({
            "target_norm": tn,
            "is_spine": tn == SPINE_DENOM,
            "n_genes_attributed": int(per_denom[tn]["dominant"].notna().sum()),
            "n_distinct_dominant_arms": int(per_denom[tn]["dom_counts"].size),
            "hub_gini": round(_gini(dc), 4),
            "top10_arm_dominance_share": round(float(np.sort(dc)[::-1][:10].sum() / dc.sum()), 4),
            "n_arms_for_50pct_dominance": int(
                (np.cumsum(np.sort(dc)[::-1]) >= 0.5 * dc.sum()).argmax() + 1
            ),
            "median_spec": round(float(per_denom[tn]["spec"].median()), 4),
        })

    summary_tbl = pd.DataFrame(rows)

    # concordance + spec-rank correlation vs the spine and vs `none`
    spine_dom = per_denom[SPINE_DENOM]["dominant"]
    spine_spec = per_denom[SPINE_DENOM]["spec"]
    for tn in per_denom:
        d = per_denom[tn]["dominant"]
        common = spine_dom.index.intersection(d.index)
        conc_spine = float((spine_dom.reindex(common) == d.reindex(common)).mean())
        sp = per_denom[tn]["spec"]
        ci = spine_spec.index.intersection(sp.index)
        spec_corr = float(
            pd.concat([spine_spec.reindex(ci), sp.reindex(ci)], axis=1)
            .dropna().corr(method="spearman").iloc[0, 1]
        ) if len(ci) > 10 else np.nan
        summary_tbl.loc[summary_tbl["target_norm"] == tn, "dominant_concordance_vs_spine"] = round(conc_spine, 4)
        summary_tbl.loc[summary_tbl["target_norm"] == tn, "spec_rankcorr_vs_spine"] = round(spec_corr, 4)

    summary_tbl.to_csv(out_dir / "attribution_by_denominator.tsv", sep="\t", index=False)

    # ---- verdict: does the denominator materially change attribution? --------
    none_row = summary_tbl.loc[summary_tbl["target_norm"] == "none"].iloc[0]
    spine_row = summary_tbl.loc[summary_tbl["target_norm"] == SPINE_DENOM].iloc[0]
    dom_conc_none_vs_spine = float(none_row["dominant_concordance_vs_spine"])
    gini_drop = float(none_row["hub_gini"]) - float(spine_row["hub_gini"])  # >0 ⇒ spine de-concentrates
    materially_changes = bool((dom_conc_none_vs_spine < 0.9) or (abs(gini_drop) > 0.03))

    summary = {
        "module": "mirna_hallmark.denominator_attribution_sweep",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "attribution_mode": ATTR_MODE, "n_genes": len(genes), "n_cohort": len(cohort), "smoke": smoke,
        "spine_denominator": SPINE_DENOM,
        "none_vs_spine_dominant_concordance": round(dom_conc_none_vs_spine, 4),
        "hub_gini_none": round(float(none_row["hub_gini"]), 4),
        "hub_gini_spine": round(float(spine_row["hub_gini"]), 4),
        "gini_change_none_to_spine": round(gini_drop, 4),
        "verdict_denominator_materially_changes_attribution": materially_changes,
        "interpretation": (
            "If True, the denominator re-assigns dominant-regulator credit and/or de-concentrates "
            "the hub structure → justified on the attribution axis (not coupling). If False, it is "
            "near-vacuous on attribution too → prefer the simplest (none) and say so."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(summary_tbl.to_string(index=False))
    print(json.dumps(summary, indent=2))
    print(f"[attr] verdict: denominator "
          f"{'MATERIALLY changes attribution' if materially_changes else 'is near-vacuous on attribution too'}"
          f"; wrote {out_dir}")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir, smoke=args.smoke)


if __name__ == "__main__":
    main()
