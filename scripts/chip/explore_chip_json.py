"""Quick exploration of ENCODE TF and histone ChIP json files."""
import json
import pandas as pd
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TF_JSON = ROOT / "data" / "CHIP" / "ENCODE" / "encode_tf_chip.json"
HIST_JSON = ROOT / "data" / "CHIP" / "ENCODE" / "encode_histone_chip.json"


def load(path: Path) -> pd.DataFrame:
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
        disease = exp.get("biosample_summary")
        rows.append({
            "accession": exp.get("accession"),
            "target": target_label,
            "biosample": biosample,
            "classification": classification,
            "organ_slims": organ_slims,
            "biosample_summary": disease,
            "assay_title": exp.get("assay_title"),
            "status": exp.get("status"),
        })
    return pd.DataFrame(rows)


def main():
    tf = load(TF_JSON)
    hist = load(HIST_JSON)
    print("=== TF ===")
    print(f"experiments={len(tf)} targets={tf['target'].nunique()} biosamples={tf['biosample'].nunique()}")
    print("status counts:", tf["status"].value_counts().to_dict())
    print("\nbreast biosample summaries:")
    breast = tf[
        tf["biosample"].fillna("").str.contains("breast|MCF|T47|HCC|HMEC|MDA|SUM|ZR-75|BT-|MB-", case=False, regex=True)
        | tf["organ_slims"].fillna("").str.contains("breast|mammary", case=False, regex=True)
        | tf["biosample_summary"].fillna("").str.contains("breast|mammary", case=False, regex=True)
    ]
    print(breast["biosample"].value_counts().head(40).to_string())

    print("\n=== HIST ===")
    print(f"experiments={len(hist)} targets={hist['target'].nunique()} biosamples={hist['biosample'].nunique()}")
    print("targets:", sorted(hist["target"].dropna().unique().tolist())[:30])
    print("breast histone biosamples:")
    breast_h = hist[
        hist["biosample"].fillna("").str.contains("breast|MCF|T47|HCC|HMEC|MDA|SUM|ZR-75|BT-|MB-", case=False, regex=True)
        | hist["organ_slims"].fillna("").str.contains("breast|mammary", case=False, regex=True)
        | hist["biosample_summary"].fillna("").str.contains("breast|mammary", case=False, regex=True)
    ]
    print(breast_h.groupby("biosample")["target"].apply(lambda s: ";".join(sorted(set(s)))).to_string())

    # Check the prior-session TF panel
    panel = ["IRF1", "NFYA", "NF-YA", "RFX5", "RELA", "p65", "IRF2", "STAT1", "STAT2",
             "IRF9", "EZH2", "JUN", "FOS", "BRD4", "BATF3", "SPI1", "PU.1", "CEBPB",
             "CIITA", "NLRC5", "JUND", "FOSL1", "FOSL2"]
    print("\nTF panel coverage in TF json (any biosample):")
    for tfn in panel:
        n_total = (tf["target"] == tfn).sum()
        n_breast = ((breast["target"] == tfn)).sum()
        if n_total or n_breast:
            biosamples_breast = sorted(breast.loc[breast["target"] == tfn, "biosample"].dropna().unique().tolist())
            print(f"  {tfn:8s} total={n_total:4d} breast={n_breast:3d} biosamples={biosamples_breast}")
        else:
            print(f"  {tfn:8s} MISSING in TF json")

if __name__ == "__main__":
    main()
