"""THE DISCOVERY GOLD SET — the A1∩site-Manakov convergent-evidence table, with a producer.

    .venv/bin/python3 -m mirna_hallmark.eval.discovery_gold_set [--report]

⭐ WHY THIS EXISTS. The discovery lane's actual deliverable is **not** a per-edge FDR result — that is
genuinely empty (0 survivors at `q_family_emp<0.05`, MH-155, and MH-197 confirmed widening the null cannot
change it). The deliverable is **convergent evidence**: edges that clear three *independent* hurdles at once.
Until now that table existed only as a hand-merge of `discovery_fall_diagnosis.tsv` and
`discoveries_sitelevel.tsv` — reconstructable, but with no producer, no provenance and no stated caveats.
**That is exactly how the five literature ground-truth sets rotted (MH-196): a headline set nobody could
regenerate.** So it gets a module, a stamp, and its caveats attached to the artifact itself.

THE THREE HURDLES (all required)
--------------------------------
  1. **Bulk coupling beyond curation** — `partial_coupling` < −0.15, partial on C + the curated HE aggregate,
     i.e. signal curation missed. (`discovery.scan_all`)
  2. **A1 = composition-robust and BH-passing** — survives the deconvolution C block and clears a *standard*
     BH; it fails only the honest heavy-tailed empirical FDR. (`fall_reason == "A1_BY_multiplicity"`)
  3. ⭐ **SITE-LEVEL physical duplex** — the arm's *predicted 3'UTR site* coincides with a **Manakov chimeric
     eCLIP duplex**, not merely "a duplex exists somewhere on this gene". (`site_manakov`)

Hurdle 3 is what makes this a claim against the null *by concept*: site-coincident edges couple stronger at
**MWU p=2.3e−27** (re-derived 2026-08-02 over the full site-level table; the handoff's 1.9e−20 was a
narrower scope), and the coincidence rate climbs monotonically up the site-confidence ladder — a gradient a
site-free arm cannot produce.

⚠⚠ TWO CAVEATS THAT MUST TRAVEL WITH THE NUMBER
-----------------------------------------------
**(1) IT IS NOT 157 INDEPENDENT FINDINGS.** By seed family: **96/157 (61%) are ONE family**
(miR-17-5p/20-5p/93-5p/106-5p/519-3p), and **139/157 (89%) belong to the miR-17~92 polycistron and its
paralogue clusters** (miR-106b~25, miR-106a~363) — co-transcribed, so their doses are not independent draws.
The honest unit is closer to **one oncogenic cluster's realized target set** than to 157 discoveries. Any
count quoted per-edge overstates the evidence; quote families (11) alongside it.

**(2) THE HANDOFF'S FLAGSHIP QUARTET IS WRONG ON ONE MEMBER.** `HANDOFF_DISCOVERY_NULL` §6 reads
*"miR-18a→{STAM2, KIF3B, MAP3K1, NEDD4}, all with Manakov duplexes."* **KIF3B has `site_manakov = False`** —
edge-level CLIP support only, no site-coincident duplex — so it belongs to the 296-edge EDGE-level set and
**drops out of the 157-edge SITE-level one**. STAM2, MAP3K1 and NEDD4 hold.

WHAT IS GENUINELY STRONG HERE
-----------------------------
**94 / 157 sit on genes with NO curated regulator at all** (`no_he_gene`) — novel gene territory, not a
re-derivation of curation. And only **3 / 157** are same-seed paralogues of a regulator already curated for
that gene (`same_seed_he`), so the set is not curation leaking back in through family membership.

⬜ **STATUS: the registry FINDING is deliberately NOT written.** `HANDOFF_DISCOVERY_NULL` §7.3 requires the
framing (convergent evidence vs per-edge significance; cluster-level vs edge-level headline) to be decided
**with the user** first. This module builds and characterises the table; it does not make the claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "discovery_gold_set.tsv"
PROV = OUT / "discovery_gold_set_provenance.json"
FALL = OUT / "discovery_fall_diagnosis.tsv"
SITE = OUT / "discoveries_sitelevel.tsv"

A1 = "A1_BY_multiplicity"
KEEP = ["gene", "arm", "seed_family", "partial_coupling", "partial_deconv", "retention", "null_z",
        "p_sitefree", "q_bh", "scanmir_rep", "arm_abundance", "chimeric_wt", "chimeric_src",
        "ev_classes", "ev_npmid", "same_seed_he", "no_he_gene", "fam_size"]
SITE_COLS = ["site_manakov", "site_clip_any", "best_type", "n_sites", "has_pred_site"]


def build() -> pd.DataFrame:
    fall = pd.read_csv(FALL, sep="\t", low_memory=False)
    site = pd.read_csv(SITE, sep="\t", low_memory=False)
    a1 = fall[fall["fall_reason"] == A1].copy()
    cols = [c for c in KEEP if c in a1.columns]
    m = a1[cols].merge(site[["gene", "arm"] + [c for c in SITE_COLS if c in site.columns]],
                       on=["gene", "arm"], how="left")
    gold = m[m["site_manakov"].fillna(False).astype(bool)].copy()
    return gold.sort_values("partial_coupling").reset_index(drop=True)


def _concentration(gold: pd.DataFrame) -> dict:
    """The independence caveat, computed rather than asserted."""
    vc = gold["seed_family"].value_counts()
    # the miR-17~92 polycistron + its two paralogue clusters (miR-106b~25, miR-106a~363)
    cluster = [f for f in vc.index if any(k in str(f) for k in
                                          ("miR-17-5p/20-5p", "miR-19-3p", "miR-25-3p", "miR-18-5p"))]
    n_cluster = int(vc[cluster].sum()) if cluster else 0
    return {"n_edges": int(len(gold)), "n_genes": int(gold.gene.nunique()),
            "n_arms": int(gold.arm.nunique()), "n_families": int(gold.seed_family.nunique()),
            "top_family": str(vc.index[0]), "top_family_n": int(vc.iloc[0]),
            "top_family_frac": round(float(vc.iloc[0] / len(gold)), 3),
            "mir17_92_cluster_n": n_cluster,
            "mir17_92_cluster_frac": round(float(n_cluster / len(gold)), 3),
            "no_he_gene_n": int(gold.no_he_gene.fillna(False).astype(bool).sum()),
            "same_seed_he_n": int(gold.same_seed_he.fillna(False).astype(bool).sum())}


def report(gold: pd.DataFrame) -> None:
    c = _concentration(gold)
    print(f"\n{'='*92}\nTHE DISCOVERY GOLD SET — A1 ∩ site-coincident Manakov duplex\n{'='*92}")
    print(f"  {c['n_edges']} edges · {c['n_genes']} genes · {c['n_arms']} arms · "
          f"⭐ {c['n_families']} SEED FAMILIES (the honest unit)\n")
    print("  --- ⚠ CONCENTRATION: this is NOT N independent findings ---")
    for k, v in gold["seed_family"].value_counts().items():
        print(f"    {v:4d}  {k}")
    print(f"\n    top family {c['top_family_n']}/{c['n_edges']} = {c['top_family_frac']:.0%}"
          f"  ·  miR-17~92 polycistron + paralogue clusters {c['mir17_92_cluster_n']}/{c['n_edges']} = "
          f"{c['mir17_92_cluster_frac']:.0%}")
    print(f"    ⇒ quote FAMILIES alongside edges; a per-edge count overstates the evidence.\n")
    print("  --- ✅ WHAT IS STRONG ---")
    print(f"    {c['no_he_gene_n']}/{c['n_edges']} on genes with NO curated regulator at all (novel territory)")
    print(f"    only {c['same_seed_he_n']}/{c['n_edges']} are same-seed paralogues of an already-curated "
          f"regulator of that gene")
    print(f"\n  --- top edges by coupling ---")
    show = [x for x in ("gene", "arm", "partial_coupling", "null_z", "retention", "best_type", "n_sites")
            if x in gold.columns]
    print(gold.head(10)[show].to_string(index=False))
    # the concept-level claim, re-derived
    site = pd.read_csv(SITE, sep="\t", low_memory=False).dropna(subset=["partial_coupling", "site_manakov"])
    a = site[site.site_manakov.astype(bool)].partial_coupling
    b = site[~site.site_manakov.astype(bool)].partial_coupling
    p = stats.mannwhitneyu(a, b, alternative="less").pvalue
    print(f"\n  --- the concept-level claim (re-derived) ---")
    print(f"    site-coincident {a.mean():+.4f} (n={len(a):,}) vs not {b.mean():+.4f} (n={len(b):,})  "
          f"MWU p={p:.3g}")
    print(f"    ⚠ HANDOFF §6 records 1.9e−20 for this — a narrower scope; the full-table value is above.")
    print(f"\n  --- ⚠ handoff flagship check ---")
    for g in ("STAM2", "KIF3B", "MAP3K1", "NEDD4"):
        hit = gold[(gold.gene == g) & gold.arm.str.contains("18a", na=False)]
        verdict = ("IN the site-level gold set" if len(hit)
                   else "⛔ NOT in it (edge-level CLIP only, no site-coincident duplex)")
        print(f"    miR-18a→{g:7s} {verdict}")
    print(f"\n  ⬜ THE REGISTRY FINDING IS DELIBERATELY UNWRITTEN — HANDOFF §7.3 requires the framing "
          f"decision with the user first.")


def run() -> pd.DataFrame:
    gold = build()
    gold.to_csv(DEST, sep="\t", index=False)
    prov = {"rule": f"fall_reason=={A1} AND site_manakov==True",
            **_concentration(gold),
            "sources": {p.name: {"bytes": p.stat().st_size,
                                 "mtime": pd.Timestamp(p.stat().st_mtime, unit="s").isoformat()}
                        for p in (FALL, SITE)},
            "sha256": hashlib.sha256(gold.to_csv(sep="\t", index=False).encode()).hexdigest()}
    PROV.write_text(json.dumps(prov, indent=2))
    print(f"-> {DEST} ({len(gold)} edges)\n-> {PROV} (sha256 {prov['sha256'][:16]}…)")
    return gold


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    g = run()
    if a.report or True:
        report(g)
