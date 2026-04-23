"""Filter ENCODE TF (and histone) ChIP JSON metadata for the APM project.

Inputs
------
data/CHIP/ENCODE/encode_tf_chip.json
data/CHIP/ENCODE/encode_histone_chip.json

Outputs (data/CHIP/ENCODE/filtered/)
------------------------------------
panel_tf_chip__all_biosamples.csv     -- panel TFs across every biosample (decision aid)
panel_tf_chip__breast_only.csv        -- panel TFs in breast/mammary biosamples
panel_tf_coverage_matrix.csv          -- TF x biosample (count of experiments) for breast biosamples
panel_tf_missing_in_breast.csv        -- panel TFs with zero breast experiments (consider non-breast surrogates)
breast_histone_chip.csv               -- histone marks in breast/mammary biosamples (full set, decision aid)

Panel definition is split into "core required" and "context optional" so it is
easy to prune later. Sources:
- Prior-session priority table (IRF1, NF-YA, RFX5, RELA, IRF2, STAT2+IRF9,
  EZH2, AP-1, BRD4, BATF3/SPI1/CEBPB).
- pipeline/config.py PRIMARY_GENES + EXTENDED_PRIMARY_GENES TFs.
- research_plan/01_gene_panel_extended.md and 02_hypothesis_catalog.md
  regulators referenced by H1, H21, H22, H41, H43.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ENCODE_DIR = ROOT / "data" / "CHIP" / "ENCODE"
OUT_DIR = ENCODE_DIR / "filtered"
OUT_DIR.mkdir(exist_ok=True)

# ---- Panel definition -----------------------------------------------------
# Aliases the user might encounter for the same factor (ENCODE label first).
TF_PANEL_CORE: dict[str, list[str]] = {
    # SXY enhanceosome (HLA-I, B2M, TAPBP, CIITA targets)
    "NFYA":   ["NFYA", "NF-YA"],
    "NFYB":   ["NFYB", "NF-YB"],
    "NFYC":   ["NFYC", "NF-YC"],
    "RFX5":   ["RFX5"],
    "RFXAP":  ["RFXAP"],
    "RFXANK": ["RFXANK"],
    "CIITA":  ["CIITA"],
    "NLRC5":  ["NLRC5"],
    # IFN-γ / JAK-STAT
    "STAT1":  ["STAT1"],
    "STAT2":  ["STAT2"],
    "STAT3":  ["STAT3"],
    "IRF1":   ["IRF1"],
    "IRF2":   ["IRF2"],
    "IRF3":   ["IRF3"],
    "IRF7":   ["IRF7"],
    "IRF9":   ["IRF9"],
    # NF-κB / Enhancer A (TNF synergy)
    "RELA":   ["RELA", "p65"],
    "NFKB1":  ["NFKB1", "p50"],
    "RELB":   ["RELB"],
    # Coactivators / chromatin
    "EP300":  ["EP300", "p300"],
    "CREBBP": ["CREBBP", "CBP"],
    "BRD4":   ["BRD4"],
    # PRC2 silencing of MHC-I (Burr 2019)
    "EZH2":   ["EZH2"],
    "SUZ12":  ["SUZ12"],
}

TF_PANEL_CONTEXT: dict[str, list[str]] = {
    # AP-1 (TGF-β interface, EMT/APM)
    "JUN":    ["JUN", "c-JUN"],
    "JUNB":   ["JUNB"],
    "JUND":   ["JUND"],
    "FOS":    ["FOS", "c-FOS"],
    "FOSL1":  ["FOSL1"],
    "FOSL2":  ["FOSL2"],
    # TGF-β / EMT
    "SMAD3":  ["SMAD3"],
    "SMAD4":  ["SMAD4"],
    # Immune-compartment (skip for bulk-tumor; useful if deconvolving)
    "BATF3":  ["BATF3"],
    "SPI1":   ["SPI1", "PU.1"],
    "CEBPB":  ["CEBPB"],
    "IRF4":   ["IRF4"],
    "IRF5":   ["IRF5"],
    "IRF8":   ["IRF8"],
    # Counter-regulators sometimes relevant at APM enhancers
    "MYC":    ["MYC", "c-MYC"],
    "BCL6":   ["BCL6"],
    "FOXM1":  ["FOXM1"],
    # CTCF (insulation context for SXY/Enhancer A)
    "CTCF":   ["CTCF"],
}


def expand_label_set(panel: dict[str, list[str]]) -> dict[str, str]:
    """Map every ENCODE-side label variant -> canonical panel name."""
    out = {}
    for canonical, labels in panel.items():
        for lab in labels:
            out[lab] = canonical
    return out


# ---- Biosample whitelist --------------------------------------------------
# Pipeline canonical breast biosamples (biosample_names.py):
#   MCF7, T47D, MCF10A, MDA-MB-231, MDA-MB-468, HCC1954, ZR751, HMEC.
# ENCODE-side spelling differs; match liberally then normalize.
BREAST_PATTERN = re.compile(
    r"\b(MCF[\s-]?7|MCF[\s-]?10A|T[\s-]?47D|MDA[\s-]?MB[\s-]?\d+|HCC\d+|"
    r"ZR[\s-]?75|HMEC|breast|mammary)\b",
    re.IGNORECASE,
)


def is_breast_row(row: pd.Series) -> bool:
    for col in ("biosample", "organ_slims", "biosample_summary"):
        v = row.get(col) or ""
        if BREAST_PATTERN.search(str(v)):
            return True
    return False


# ---- IO -------------------------------------------------------------------
def load_chip_json(path: Path) -> pd.DataFrame:
    with open(path) as f:
        data = json.load(f)
    rows = []
    for exp in data.get("@graph", []):
        target = exp.get("target", {})
        target_label = target.get("label") if isinstance(target, dict) else None
        bs = exp.get("biosample_ontology", {})
        biosample = bs.get("term_name") if isinstance(bs, dict) else None
        classification = bs.get("classification") if isinstance(bs, dict) else None
        organ_slims = bs.get("organ_slims") if isinstance(bs, dict) else None
        if isinstance(organ_slims, list):
            organ_slims = ";".join(organ_slims)
        rows.append({
            "accession": exp.get("accession"),
            "target": target_label,
            "biosample": biosample,
            "classification": classification,
            "organ_slims": organ_slims,
            "biosample_summary": exp.get("biosample_summary"),
            "assay_title": exp.get("assay_title"),
            "status": exp.get("status"),
        })
    return pd.DataFrame(rows)


# ---- Main -----------------------------------------------------------------
def filter_tf(tf_df: pd.DataFrame) -> None:
    full_panel = {**TF_PANEL_CORE, **TF_PANEL_CONTEXT}
    label_map = expand_label_set(full_panel)
    tier_map = {
        canonical: ("core" if canonical in TF_PANEL_CORE else "context")
        for canonical in full_panel
    }

    # Keep only released experiments
    tf = tf_df[tf_df["status"] == "released"].copy()
    tf["panel_name"] = tf["target"].map(label_map)
    panel = tf.dropna(subset=["panel_name"]).copy()
    panel["tier"] = panel["panel_name"].map(tier_map)

    # All biosamples
    panel_all = panel.sort_values(["tier", "panel_name", "biosample"])
    panel_all.to_csv(OUT_DIR / "panel_tf_chip__all_biosamples.csv", index=False)

    # Breast/mammary subset
    is_breast = panel.apply(is_breast_row, axis=1)
    panel_breast = panel[is_breast].sort_values(["tier", "panel_name", "biosample"])
    panel_breast.to_csv(OUT_DIR / "panel_tf_chip__breast_only.csv", index=False)

    # Coverage matrix (panel TF x breast biosample, count of experiments)
    cov = (
        panel_breast.assign(n=1)
        .pivot_table(index=["tier", "panel_name"], columns="biosample",
                     values="n", aggfunc="sum", fill_value=0)
        .sort_index()
    )
    cov.to_csv(OUT_DIR / "panel_tf_coverage_matrix__breast.csv")

    # Missing-in-breast (which TFs would need a non-breast surrogate)
    present_breast = set(panel_breast["panel_name"])
    missing_rows = []
    for canonical in full_panel:
        if canonical in present_breast:
            continue
        non_breast_n = (panel["panel_name"] == canonical).sum()
        non_breast_bs = sorted(
            panel.loc[panel["panel_name"] == canonical, "biosample"]
            .dropna().unique().tolist()
        )
        missing_rows.append({
            "panel_name": canonical,
            "tier": tier_map[canonical],
            "n_experiments_any_biosample": int(non_breast_n),
            "non_breast_biosamples": ";".join(non_breast_bs),
        })
    pd.DataFrame(missing_rows).sort_values(
        ["tier", "n_experiments_any_biosample"], ascending=[True, False]
    ).to_csv(OUT_DIR / "panel_tf_missing_in_breast.csv", index=False)

    print("\n--- TF FILTER SUMMARY ---")
    print(f"input experiments (released): {len(tf)}")
    print(f"panel hits (all biosamples):  {len(panel_all)}  "
          f"covering {panel_all['panel_name'].nunique()} TFs")
    print(f"panel hits (breast/mammary):  {len(panel_breast)}  "
          f"covering {panel_breast['panel_name'].nunique()} TFs")
    print("breast biosamples touched:",
          sorted(panel_breast["biosample"].dropna().unique().tolist()))
    print("\nbreast coverage matrix:")
    print(cov.to_string())
    print(f"\npanel TFs with NO breast experiment: "
          f"{len(full_panel) - panel_breast['panel_name'].nunique()}")


def filter_histone(hist_df: pd.DataFrame) -> None:
    is_breast = hist_df.apply(is_breast_row, axis=1)
    breast = hist_df[is_breast].copy()
    breast.to_csv(OUT_DIR / "breast_histone_chip.csv", index=False)
    cov = (
        breast.assign(n=1)
        .pivot_table(index="target", columns="biosample",
                     values="n", aggfunc="sum", fill_value=0)
    )
    cov.to_csv(OUT_DIR / "breast_histone_coverage_matrix.csv")
    print("\n--- HISTONE FILTER SUMMARY ---")
    print(f"breast histone experiments: {len(breast)}")
    print(f"distinct marks: {breast['target'].nunique()}")
    print(f"distinct breast biosamples: {breast['biosample'].nunique()}")
    print(cov.to_string())


def main():
    tf_df = load_chip_json(ENCODE_DIR / "encode_tf_chip.json")
    hist_df = load_chip_json(ENCODE_DIR / "encode_histone_chip.json")
    filter_tf(tf_df)
    filter_histone(hist_df)
    print(f"\nOutputs in: {OUT_DIR}")


if __name__ == "__main__":
    main()
