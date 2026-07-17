"""Build a compact curated human transcription-factor census for the within-set hierarchy weight.

Source: Lambert et al. 2018, "The Human Transcription Factors" (Cell) — humantfs.ccbr.utoronto.ca
`DatabaseExtract_v_1.01.csv` (saved as `annotations/humantfs_lambert2018_raw.csv`). We keep the
1,639 genes flagged `Is TF? == Yes` and emit `annotations/humantfs_lambert2018_tf.tsv` with the DBD
family + TF assessment.

Why this complements the OmniPath-induced topology: the architecture `w_arch` (reverse-PageRank on the
induced signed subgraph) only *weakly* tracks curated TF identity — ~48% of curated TFs fall below the
median `w_arch` because transcriptional regulation is under-represented in the signed-PPI network. A
curated identity flag (`is_tf`) is therefore an independent, coverage-robust hierarchy annotation that
can be manually weighted on top of the topology.
"""
from __future__ import annotations
import pandas as pd

from mirna_hallmark import config as C

RAW = C.REPO_ROOT / "annotations/humantfs_lambert2018_raw.csv"
OUT = C.REPO_ROOT / "annotations/humantfs_lambert2018_tf.tsv"


def main() -> None:
    df = pd.read_csv(RAW)
    tf = df[df["Is TF?"] == "Yes"].copy()
    out = pd.DataFrame({
        "gene": tf["HGNC symbol"].astype(str),
        "is_tf": True,
        "dbd": tf["DBD"].astype(str),
        "tf_assessment": tf["TF assessment"].astype(str),
    }).dropna(subset=["gene"]).drop_duplicates("gene")
    out = out[out["gene"] != "nan"]
    out.to_csv(OUT, sep="\t", index=False)
    print(f"[tf_census] wrote {len(out)} curated TFs -> {OUT}")


if __name__ == "__main__":
    main()
