"""Phase 3 — program-wise learned-M run (Design §6 Phase 3 §4): fit the model WITHIN each Hallmark program
(not a flat universe run), on the 5 priority programs EMT / P53 / G2M / IFN-γ / HYPOXIA.

Per HE gene (family-collapsed, ledger prior): the OOF coupling gate (Bars 1–2 — beats raw abundance /
matches curated fixed-M) + the persisted learned M. Per program: aggregate gate + the identifiability
object (hierarchical posterior |β/sd|>2, Decision F). One assemble per gene (no redundant loads).

Writes: output/learned/programs/{program}_{genes,summary}.tsv + identified_edges.tsv.
CLI: `python -m mirna_hallmark.learned.programs`
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import data_loaders as D
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.learned import hierarchical as H
from mirna_hallmark.learned import mvp as M

PRIORITY = ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_P53_PATHWAY",
            "HALLMARK_G2M_CHECKPOINT", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_HYPOXIA"]
OUT = Path("mirna_hallmark/output/learned/programs")


def _program_genes() -> dict:
    hs = HallmarkSets.load()
    sets = hs.sets if hasattr(hs, "sets") else {}
    he = set(D.high_evidence_edges()["gene"].dropna().astype(str))
    return {p.replace("HALLMARK_", ""): sorted(set(sets.get(p, set())) & he) for p in PRIORITY}


def run(*, alpha: float = 0.005, deconv: bool = False, batch: bool = True) -> pd.DataFrame:
    """Program-wise gate + identifiability. deconv=True adds CIBERSORTx non-malignant composition to C, so
    the identified edges are CELL-INTRINSIC (the stroma-driven miR-29→collagen etc. drop out) — the
    composition pass (#0, 2026-07-05: those edges retain only 6–44% under deconv)."""
    tag = "_deconv" if deconv else ""
    OUT.mkdir(parents=True, exist_ok=True)
    progs = _program_genes()
    prog_summ, ident_all = [], []
    for pname, genes in progs.items():
        rows = []
        for i, g in enumerate(genes):
            if i % 30 == 0:
                print(f"[{pname}{tag}] {i}/{len(genes)} gate...", flush=True)
            try:
                rows.append(M.oof_gate(g, alpha=alpha, w_prior_source="ledger", family=True,
                                       deconv=deconv, batch=batch))
            except Exception as e:
                rows.append({"gene": g, "error": repr(e)[:50]})
        gate = pd.DataFrame(rows)
        gate.to_csv(OUT / f"{pname}{tag}_genes.tsv", sep="\t", index=False)
        ok = gate[gate.get("rho_model").notna()] if "rho_model" in gate else gate.iloc[:0]
        # identifiability: hierarchical posterior sd on this program (one shared fit)
        try:
            ident = H.uncertainty(genes, family=True, deconv=deconv)
            ident["program"] = pname
            ident_all.append(ident[ident["z"].abs() > 2])
            n_ident = int((ident["z"].abs() > 2).sum())
        except Exception as e:
            print(f"[{pname}] uncertainty failed: {repr(e)[:60]}"); n_ident = -1
        prog_summ.append({
            "program": pname, "n_genes": len(genes), "n_fit": len(ok),
            "mean_rho_model": round(float(ok["rho_model"].mean()), 3) if len(ok) else np.nan,
            "beats_abundance": f"{int(ok['vs_abund'].sum())}/{len(ok)}" if len(ok) else "0/0",
            "matches_curated": f"{int(ok['vs_curated'].sum())}/{len(ok)}" if len(ok) else "0/0",
            "mean_stability": round(float(ok["stability"].mean()), 2) if len(ok) else np.nan,
            "n_identified_edges": n_ident,
        })
        print(f"  → {pname}: {prog_summ[-1]}", flush=True)
    summ = pd.DataFrame(prog_summ)
    summ.to_csv(OUT / f"program_summary{tag}.tsv", sep="\t", index=False)
    if ident_all:
        idf = pd.concat(ident_all, ignore_index=True).sort_values("beta")
        idf.to_csv(OUT / f"identified_edges{tag}.tsv", sep="\t", index=False)
    print(f"\n=== PROGRAM-WISE M SUMMARY{tag} (5 priority Hallmarks) ===")
    print(summ.to_string(index=False))
    print(f"\nwrote {OUT}/ (per-program genes + summary + identified_edges, tag='{tag}')")
    return summ


if __name__ == "__main__":
    import sys
    run(deconv="--deconv" in sys.argv)
