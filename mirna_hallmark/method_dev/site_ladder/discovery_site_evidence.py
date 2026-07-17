"""Site-level chimeric evidence for the DISCOVERY candidates (MH-155 wire).

Extends the site ladder (`site_filter_ladder` → `site_genomic_l5`) from the curated HE universe to the discovery
orphan candidates, so each candidate edge carries **site-level** occupancy — does its *predicted 3'UTR site*
physically coincide with an observed chimeric duplex (Manakov eCLIP / TarBase CLASH) of the SAME miRNA in the
SAME gene? This is the sharper form of `chimeric_evidence.py`'s edge-level weight (which only asks whether a
duplex exists *somewhere* on the gene).

WHY IT MATTERS (measured 2026-07-17): site-level Manakov overlap tracks coupling at **MWU p=1.9e−20** vs the
edge-level p=4.6e−8 — localising to the site sharpens the convergent-evidence signal by orders of magnitude. It
does not rescue per-edge FDR (the honest null is unbeatable per-edge), but it is a claim *against the null by
concept*: a site-free arm has no predicted site to coincide with a duplex, and the coincidence rate climbs
monotonically up the site-confidence ladder (7mer-A1 Manakov 9% → 8mer+conserved+3'-supp 27%).

Run: `.venv/bin/python3 -m mirna_hallmark.method_dev.site_ladder.discovery_site_evidence`
Out: `output/learned/discoveries_sitelevel.tsv` (discoveries + site_manakov / site_clip_any / best_type / n_sites).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.method_dev.site_ladder import site_filter_ladder as SFL
from mirna_hallmark.method_dev.site_ladder import site_genomic_l5 as L5

_DISC = C.OUTPUT_ROOT / "learned" / "discoveries.tsv"
_OUT = C.OUTPUT_ROOT / "learned" / "discoveries_sitelevel.tsv"


def build(discoveries: str | Path = _DISC) -> pd.DataFrame:
    d = pd.read_csv(discoveries, sep="\t")
    edges = d[["arm", "gene"]].drop_duplicates().rename(columns={"arm": "miRNA"})
    print(f"[disc-site] scanning {len(edges):,} discovery candidate edges through the site ladder", flush=True)
    # extend the ladder to the orphan edges — SEPARATE output names so the HE artifacts are preserved
    SFL.build(edges=edges, out_name="utr_site_ladder_discovery.tsv.gz")
    L5.build(in_name="utr_site_ladder_discovery.tsv.gz", out_name="utr_site_ladder_genomic_discovery.tsv.gz")
    g = pd.read_csv(SFL.OUT / "utr_site_ladder_genomic_discovery.tsv.gz", sep="\t")
    # collapse sites → per-edge site-level occupancy (a candidate has support if ANY of its sites overlaps)
    site = (g.groupby(["miRNA", "gene"])
              .agg(site_manakov=("site_manakov", "max"), site_tarbase=("site_tarbase", "max"),
                   site_clip_any=("site_clip_any", "max"), best_type=("type", "first"), n_sites=("type", "size"))
              .reset_index().rename(columns={"miRNA": "arm"}))
    m = d.merge(site, on=["arm", "gene"], how="left")
    for c in ("site_manakov", "site_tarbase", "site_clip_any"):
        m[c] = m[c].fillna(False).astype(bool)
    m["has_pred_site"] = m["n_sites"].notna()
    m.to_csv(_OUT, sep="\t", index=False)
    print(f"[disc-site] {int(m['has_pred_site'].sum()):,}/{len(m):,} candidates have a predicted 7mer+ site | "
          f"site overlaps a Manakov duplex {int(m['site_manakov'].sum()):,} | any chimeric "
          f"{int(m['site_clip_any'].sum()):,} → wrote {_OUT}")
    return m


if __name__ == "__main__":
    build()
