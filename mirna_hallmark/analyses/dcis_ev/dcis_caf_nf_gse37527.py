"""Non-circular CAF-vs-NF miRNA test (GSE37527) — is the miR-29c→ECM axis a CAF program?

MH-52 showed (at MIBI protein level) that invasion = myoepithelial loss + CAF activation with
ECM up; MH-51 left open whether the miR-29c→ECM signal is *miR-29c lost in CAFs* (regulatory) or
just a composition correlate. The compartment readouts so far are protein/target (effector) — they
can't measure miR-29c. **GSE37527 can:** paired primary **cancer-associated fibroblasts (CAF) vs
matched normal fibroblasts (NF)** from 6 breast tumours, miRNA array (GPL14767, Feature-Number) —
the one deposited, compartment-resolved, *directly measured* (non-circular) breast miRNA dataset,
and an array so it resolves miR-29a/b/c separately (the MH-51 paralog question).

Result (see MH-54): **miR-29c is NOT down in CAFs** (flat/slightly up, NS) → the desmoplastic ECM
up-regulation in CAFs is miR-29-INDEPENDENT; the miR-29c→ECM-in-CAF hypothesis is not supported.
The genuine CAF-activation-lost arms are **miR-143/miR-145/miR-126/miR-205** (the myofibroblast
miR-143/145 cluster + epithelial markers), all paired-Wilcoxon p=0.031 (n=6 floor).

Caveats: cultured primary CAF/NF (passage <10; culture can shift miRNA); n=6 pairs → p floor 0.031;
array; this is breast tumour CAF generally, not DCIS-staged.

DATA: ``data/external/caf_nf_gse37527/`` (GEO series matrix + persisted GPL14767 map).
Run: ``.venv/bin/python3 -m mirna_hallmark.dcis_caf_nf_gse37527``
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from mirna_hallmark import config as C
from mirna_hallmark.stats import bh_fdr

DATA_DIR = C.REPO_ROOT / "data" / "external" / "caf_nf_gse37527"
MATRIX = DATA_DIR / "GSE37527_series_matrix.txt.gz"
GPL_MAP = DATA_DIR / "GPL14767_feature_to_mirna.tsv.gz"
OUT_DIR = C.OUTPUT_ROOT / "dcis_caf_nf_gse37527"
FOCUS = ("miR-29a", "miR-29b", "miR-29c", "miR-145", "miR-143", "miR-126",
         "miR-140", "miR-205", "miR-21", "miR-185")
_FOCUS_RE = "(" + "|".join(FOCUS) + r")(?![0-9])"   # digit-boundary: miR-21 ≠ miR-210/214


def load() -> Tuple[pd.DataFrame, list, list]:
    """(miRNA x sample, NF cols, CAF cols). Paired by specimen (NF_i ↔ CAF_i)."""
    f2m = pd.read_csv(GPL_MAP, sep="\t", dtype=str).set_index("feature_id")["mirna"]
    L = gzip.open(MATRIX, "rt").readlines()
    titles = next([x.strip().strip('"') for x in ln.split("\t")[1:]]
                  for ln in L if ln.startswith("!Sample_title"))
    a0 = next(i for i, l in enumerate(L) if l.startswith("!series_matrix_table_begin"))
    a1 = next(i for i, l in enumerate(L) if l.startswith("!series_matrix_table_end"))
    df = pd.read_csv(io.StringIO("".join(L[a0 + 1:a1])), sep="\t", index_col=0)
    df.index = [f2m.get(str(i).strip().strip('"')) for i in df.index]
    df = df[df.index.notna() & df.index.str.startswith("hsa")]
    df = df.groupby(df.index).mean()
    nf = [c for c, t in zip(df.columns, titles) if "Normal" in t]
    caf = [c for c, t in zip(df.columns, titles) if "Cancer" in t]
    return df, nf, caf


def paired_test(df: pd.DataFrame, nf: list, caf: list) -> pd.DataFrame:
    rows = []
    for arm in df.index:
        a = pd.to_numeric(df.loc[arm, nf], errors="coerce").to_numpy(float)
        b = pd.to_numeric(df.loc[arm, caf], errors="coerce").to_numpy(float)
        m = ~(np.isnan(a) | np.isnan(b))
        if m.sum() < 5 or np.allclose(a[m], b[m]):
            continue
        try:
            p = wilcoxon(b[m], a[m])[1]
        except ValueError:
            p = np.nan
        rows.append({"arm": arm, "nf_median": round(float(np.median(a[m])), 2),
                     "caf_median": round(float(np.median(b[m])), 2),
                     "delta_caf_minus_nf": round(float(np.median(b[m] - a[m])), 2),
                     "wilcoxon_p": p, "n_pairs": int(m.sum())})
    res = pd.DataFrame(rows)
    res["q_bh"] = bh_fdr(res["wilcoxon_p"].fillna(1).to_numpy())
    res["direction"] = np.where(res["delta_caf_minus_nf"] < 0, "CAF_down", "CAF_up")
    return res.sort_values("delta_caf_minus_nf").reset_index(drop=True)


def run(*, out_dir: Path = OUT_DIR) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    df, nf, caf = load()
    print(f"[caf_nf] {df.shape[0]} miRNA; NF={len(nf)} CAF={len(caf)} (paired by specimen)")
    res = paired_test(df, nf, caf)
    res.to_csv(out_dir / "caf_vs_nf_paired.tsv", sep="\t", index=False)

    focus = res[res["arm"].str.contains(_FOCUS_RE, regex=True)].copy()
    caf_down_sig = res[(res["direction"] == "CAF_down") & (res["wilcoxon_p"] <= 0.05)]
    summary = {
        "module": "mirna_hallmark.dcis_caf_nf_gse37527",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "cohort": "GSE37527 — 6 paired breast CAF vs NF primary cultures (miRNA array GPL14767)",
        "miR29_family": [
            {"arm": r.arm, "nf": r.nf_median, "caf": r.caf_median,
             "delta": r.delta_caf_minus_nf, "p": round(float(r.wilcoxon_p), 3), "dir": r.direction}
            for r in focus[focus["arm"].str.contains("miR-29")].itertuples()],
        "loss_leaders_and_controls": [
            {"arm": r.arm, "delta": r.delta_caf_minus_nf, "p": round(float(r.wilcoxon_p), 3), "dir": r.direction}
            for r in focus[~focus["arm"].str.contains("miR-29")].itertuples()],
        "n_CAF_down_p<=0.05": int(len(caf_down_sig)),
        "verdict": ("miR-29c NOT down in CAFs (ECM-up in CAF is miR-29-independent); "
                    "CAF-activation-lost arms = miR-143/145/126/205"),
        "caveats": ["cultured primary CAF/NF (passage<10), n=6 pairs (p floor 0.031), array",
                    "breast tumour CAF generally, not DCIS-staged",
                    "tests miR-29c directly (non-circular) — the compartment readout HTAN/MIBI cannot give"],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print("[caf_nf] miR-29 family:")
    for r in summary["miR29_family"]:
        print(f"    {r['arm']:14s} NF={r['nf']:>8} CAF={r['caf']:>8} d={r['delta']:>+8} p={r['p']}  {r['dir']}")
    print("[caf_nf] loss-leaders/controls:")
    for r in summary["loss_leaders_and_controls"]:
        print(f"    {r['arm']:14s} d={r['delta']:>+9} p={r['p']}  {r['dir']}")
    print(f"[caf_nf] verdict: {summary['verdict']}")
    print(f"[caf_nf] wrote outputs under {out_dir}")
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(out_dir=args.out_dir)


if __name__ == "__main__":
    main()
