"""Measure the card layer: dimensions, registry coverage, per-column fill.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.card_scan

Writes `output/learned/card_coverage.json` and prints a summary. This is the
producer behind `docs/derived/DOSSIER.html`'s column index — the dossier must
never carry a hand-typed coverage number (the producer-less shape that rotted
the literature sets, MH-196/MH-218).

What it checks, beyond fill:
  * every card column is registered in `card_registry.tsv`  (unregistered -> gap)
  * every registry entry points at a column that exists     (registry-only -> stale)
Both were 0 on 2026-08-18. `gen_cards.py --check` is the stricter gate; this one
exists to produce NUMBERS, not to pass/fail.

Fill is `notna().sum() / n_rows`. Read it beside `domain`, which states which rows
a column is DEFINED on: a low fill with a domain entry is the domain doing its job,
a low fill WITHOUT one is a real hole (CARD_RUNG_DOCTRINE §2).
"""
import json
import pathlib
import statistics as st

import pandas as pd

BASE = pathlib.Path(__file__).resolve().parents[2] / "output" / "learned"

CARDS = {
    "arm":         (BASE / "arm_card.tsv",             ["arm"]),
    "edge":        (BASE / "realization/edge_card.tsv", ["gene", "arm"]),
    "gene":        (BASE / "realization/gene_card.tsv", ["gene"]),
    "gene_family": (BASE / "gene_family_card.tsv",      ["gene", "family"]),
    "seed_family": (BASE / "seed_family_card.tsv",      ["seed_family"]),
}


def scan():
    reg = pd.read_csv(BASE / "card_registry.tsv", sep="\t", dtype=str).fillna("")
    regidx = {(r.card, r.column): r for r in reg.itertuples()}
    out = {"cards": {}, "columns": []}

    for name, (path, key) in CARDS.items():
        df = pd.read_csv(path, sep="\t", low_memory=False)
        n = len(df)
        out["cards"][name] = {"path": str(path.relative_to(BASE.parents[1])),
                              "key": key, "rows": n, "cols": df.shape[1],
                              "bytes": path.stat().st_size}
        for c in df.columns:
            nn = int(df[c].notna().sum())
            r = regidx.get((name, c))
            out["columns"].append({
                "card": name, "col": c, "n": nn,
                "frac": round(nn / n, 4) if n else 0.0,
                "rung": (r.rung if r is not None else ""),
                "agg_of": (r.agg_of if r is not None else ""),
                "domain": (r.domain if r is not None else ""),
                "registered": r is not None,
            })

    have = {(c["card"], c["col"]) for c in out["columns"]}
    out["registry_only"] = [{"card": r.card, "col": r.column}
                            for r in reg.itertuples() if (r.card, r.column) not in have]
    out["unregistered"] = [{"card": c["card"], "col": c["col"]}
                           for c in out["columns"] if not c["registered"]]
    return out


def main():
    out = scan()
    (BASE / "card_coverage.json").write_text(json.dumps(out))

    print("cards:", {k: (v["rows"], v["cols"]) for k, v in out["cards"].items()})
    print(f"columns: {len(out['columns'])} | unregistered: {len(out['unregistered'])} "
          f"| registry-only: {len(out['registry_only'])}")
    for nm in CARDS:
        cs = [c for c in out["columns"] if c["card"] == nm]
        wd = [c["frac"] for c in cs if c["domain"]]
        nd = [c["frac"] for c in cs if not c["domain"]]
        print(f"  {nm:12s} with-domain n={len(wd):3d} med_fill="
              f"{st.median(wd) if wd else float('nan'):.3f}  |  no-domain n={len(nd):3d} med_fill="
              f"{st.median(nd) if nd else float('nan'):.3f}")


if __name__ == "__main__":
    main()
