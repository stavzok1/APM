"""Build a compact per-gene DepMap CRISPR dependency summary for the continuous (all-gene) role weight.

Source: DepMap 26Q1 `CRISPRGeneEffect.csv` (Chronos gene-effect; rows = cell-line ModelID, columns =
"SYMBOL (ENTREZ)"). Chronos gene-effect is signed: strongly **negative** = essential (knockout depletes
the line) = the tumour *depends* on the gene (oncogene-like); ~0 = non-essential; **positive** = knockout
*increases* fitness (tumour-suppressor-like).

We collapse the matrix to per-gene means over (i) the 96-line BREAST panel (`OncotreeLineage=="Breast"`)
and (ii) the full pan-cancer panel, and emit `data/CCLE/depmap_gene_effect_summary.tsv`. The downstream
role weight orients it so that **r(g) = gene_effect(g)** reads as "repressing g is pro-tumor" (>0, TSG-like)
vs "anti-tumor / hits a dependency" (<0, oncogene-like) — the continuous, all-gene analogue of the binary
`malignancy_sign`. CRISPR essentiality is one-sided (strong on dependencies, weak on TSG direction), which
the consumer documents.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd

from mirna_hallmark import config as C

CCLE = C.REPO_ROOT / "data/CCLE"
SRC = CCLE / "CRISPRGeneEffect.csv"
MODEL = CCLE / "Model.csv"
OUT = CCLE / "depmap_gene_effect_summary.tsv"


def main() -> None:
    model = pd.read_csv(MODEL)
    breast = set(model.loc[model["OncotreeLineage"] == "Breast", "ModelID"])
    print(f"[depmap] {len(breast)} breast ModelIDs in metadata")

    df = pd.read_csv(SRC, index_col=0)
    df.index = df.index.astype(str)
    print(f"[depmap] gene-effect matrix: {df.shape[0]} lines x {df.shape[1]} genes")

    breast_rows = df.index.intersection(breast)
    print(f"[depmap] {len(breast_rows)} breast lines with CRISPR data")

    pan_mean = df.mean(axis=0)
    breast_mean = df.loc[breast_rows].mean(axis=0)

    def _sym(col: str) -> str:
        return col.split(" (")[0].strip()

    out = pd.DataFrame({
        "gene": [_sym(c) for c in df.columns],
        "dep_effect_breast": breast_mean.values,
        "dep_effect_pan": pan_mean.values,
    })
    out["n_breast_lines"] = len(breast_rows)
    # collapse any duplicate symbols (rare) by mean
    out = out.groupby("gene", as_index=False).mean(numeric_only=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"[depmap] wrote {len(out)} genes -> {OUT}")
    print(out["dep_effect_breast"].describe().round(3).to_dict())


if __name__ == "__main__":
    main()
