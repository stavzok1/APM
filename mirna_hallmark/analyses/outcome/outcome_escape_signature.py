"""Roadmap Phase 3 outcome — regulatory ESCAPE + functional TS-miRNA-panel signature (push on MH-67).

MH-67 showed functional repression (anti-corr) out-prognoses abundance. Two follow-ups:

  ESCAPE (per miRNA, per patient): escape = -residual(anti_corr ~ abundance) — repression *shortfall*
    given the miRNA level (miRNA abundant but targets NOT repressed = escaped control). Tests whether
    *dysfunction* beats either level alone, per miRNA and panel-aggregated.

  FUNCTIONAL PANEL signature: mean anti-corr over the TS-miRNA panel = per-patient "TS-miRNA functional
    activity", vs the abundance-based panel and vs an oncomiR panel. Does functional > abundance panel?

TCGA (DFI+OS) + Buffa (DRFS), composition/ER-adjusted, head-to-head.

Run: ``.venv/bin/python3 -m mirna_hallmark.outcome_escape_signature``
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr
from mirna_hallmark.analyses.outcome.outcome_survival import _cox_one, _z
from mirna_hallmark.analyses.outcome.outcome_famous_mirnas import FAMOUS, _edges_by_arm
from mirna_hallmark.analyses.outcome.outcome_famous_compare import _tcga, _buffa

OUT_DIR = C.OUTPUT_ROOT / "outcome_escape_signature"


def _per_mirna(coh, eba) -> dict:
    out = {}
    for arm, role in FAMOUS.items():
        a = coh["arm_of"](arm)
        if a is None or a not in coh["mi"].index:
            continue
        targets, _ = eba.get(arm, ([], None))
        tg = [g for g in targets if g in coh["rna"].index]
        if len(tg) < 5:
            continue
        sub = coh["rna"].loc[tg]
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
        ab = _z(coh["mi"].loc[a]); anti = -z.mean(axis=0)
        df = pd.concat([ab.rename("ab"), anti.rename("anti")], axis=1).dropna()
        if len(df) < 30:
            continue
        b = np.polyfit(df["ab"], df["anti"], 1)
        escape = -(anti - (b[0] * ab + b[1]))                          # high = abundant but under-repressing
        out[arm] = {"role": role, "abundance": ab, "anti_corr": anti, "escape": escape}
    return out


def _cox(coh, x, fam, feat, rows):
    for ep, (t, e) in coh["endpoints"].items():
        c, p = _cox_one(t, e, x.astype(float), coh["cov"])
        rows.append({"cohort": coh["name"], "family": fam, "feature": feat, "endpoint": ep, "hr": np.exp(c), "p": p})


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    eba = _edges_by_arm()
    rows = []
    for loader in (_tcga, _buffa):
        coh = loader(); pm = _per_mirna(coh, eba)
        # per-miRNA escape (+ abundance/anti for reference)
        for arm, d in pm.items():
            for metric in ("abundance", "anti_corr", "escape"):
                _cox(coh, d[metric], f"permirna_{metric}", arm, rows)
        # panels
        for pname, sel in {"TS": [m for m in pm if pm[m]["role"].startswith("TS")],
                           "onco": [m for m in pm if pm[m]["role"].startswith("onco")]}.items():
            if len(sel) < 3:
                continue
            for metric in ("abundance", "anti_corr", "escape"):
                score = pd.concat([pm[m][metric] for m in sel], axis=1).mean(axis=1)
                _cox(coh, score, "panel", f"{pname}_{metric}", rows)
        print(f"[escape] {coh['name']}: {len(pm)} miRNAs")
    res = pd.DataFrame(rows)
    for key, idx in res.groupby(["cohort", "family", "endpoint"]).groups.items():
        ok = res.loc[idx, "p"].notna(); res.loc[res.loc[idx].index[ok], "q"] = bh_fdr(res.loc[idx, "p"].dropna().to_numpy())
    res.to_csv(out_dir / "escape_signature.tsv", sep="\t", index=False)

    pe = res[res.family.str.startswith("permirna")]
    win = pe.pivot_table(index=["cohort", "feature", "endpoint"], columns="family", values="p")
    wcols = [c for c in ("permirna_abundance", "permirna_anti_corr", "permirna_escape") if c in win.columns]
    summary = {
        "module": "mirna_hallmark.outcome_escape_signature",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "permirna_FDR_by_metric": {m.split("_", 1)[1]: int(((pe.family == m) & (pe.q < 0.1)).sum()) for m in pe.family.unique()},
        "permirna_median_neglog10p": {m.split("_", 1)[1]: round(float(np.median(-np.log10(pe[pe.family == m]["p"].dropna()))), 3) for m in pe.family.unique()},
        "escape_wins_head2head(min p)": (win[wcols].idxmin(axis=1).value_counts().to_dict() if len(wcols) == 3 else {}),
        "panel_results": res[res.family == "panel"].sort_values("p").apply(lambda r: f"{r['cohort']}/{r['feature']}/{r['endpoint']} HR={r['hr']:.2f} p={r['p']:.3f} q={r.get('q',np.nan):.3f}", axis=1).tolist(),
        "escape_FDR_hits": pe[(pe.family == "permirna_escape") & (pe.q < 0.1)].sort_values("p").apply(lambda r: f"{r['cohort']}/{r['feature']}/{r['endpoint']} HR={r['hr']:.2f} q={r['q']:.3f}", axis=1).tolist(),
        "caveats": ["escape = -resid(anti~abundance); panel = mean over TS/onco miRNAs; TCGA+Buffa"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"[escape] per-miRNA FDR by metric: {summary['permirna_FDR_by_metric']}")
    print(f"[escape] per-miRNA median -log10p: {summary['permirna_median_neglog10p']}")
    print(f"[escape] escape head-to-head wins: {summary['escape_wins_head2head(min p)']}")
    print(f"[escape] escape FDR hits: {summary['escape_FDR_hits']}")
    print(f"[escape] PANEL: {summary['panel_results']}")
    print(f"[escape] wrote {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
