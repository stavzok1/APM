"""THE LITERATURE GROUND TRUTH — a VERSIONED, AUDITABLE, REPRODUCIBLE canonical-regulator set.

    .venv/bin/python3 -m mirna_hallmark.eval.lit_ground_truth [--min-pmid 3] [--report]

⭐ WHY THIS EXISTS. The attribution axis — *can the model NAME the miRNA family that represses this gene?* —
is the one place the program tests identity rather than coupling, and its verdict (MH-124 §4b, MH-126c) rests
on a hand-curated literature set. **That set cannot be reproduced.** Audited 2026-08-02:

    literature_headtohead.tsv        n= 16   producer: ABSENT
    between_family_literature.tsv    n= 19   producer: ABSENT
    variance_attribution.tsv         n= 21   producer: ABSENT       <- MH-124 §4b's n
    cn_family_literature.tsv         n= 23   producer: ABSENT
    literature_power_test.tsv        n=173   producer: ABSENT
    MH-126c TEST 1                   n= 32   NOT ON DISK AT ALL     <- the REGISTERED headline

All five were written 2026-07-13 by `/tmp` scratchpad scripts that no longer exist. Five mutually
inconsistent lists, no definition, no version, no way to regenerate any of them. **That is not a power
problem, it is an auditability problem**, and it is why the board's caveat — *"β does not attribute" is
MEASURED; "the model cannot attribute" is UNDERPOWERED* — could never be settled.

THE RULE (one line, mechanical, re-derivable from the repo)
-----------------------------------------------------------
    canonical family of gene G  =  argmax over seed families of
                                   #DISTINCT PubMed IDs carrying LOW-THROUGHPUT FUNCTIONAL evidence
                                   for any member arm -> G

`LT_FUNC = {reporter, western, proteomics}` — luciferase reporter, immunoblot, mass-spec — with miRTarBase
"Weak" and TarBase "Positive" rows excluded. **Deliberately NOT `qpcr_rna`** (it lumps high-throughput
RNA-seq/microarray in with qPCR) and **NOT `ago_clip`/`chimeric`** (binding is not repression, and those are
the very orthogonal channels the discovery lane uses as INDEPENDENT evidence — spending them here would
couple the two arcs).

Source of record: `learned/evidence/ledger.py` — PMID-deduped globally across miRTarBase × TarBase-v9 by
(edge × PMID × assay_class), so a paper reported by both databases counts once. Versioned by its two input
files; this module stamps their size+mtime and a content hash of its own output into
`lit_ground_truth_provenance.json` so any number quoted from the set can be traced to the exact pull.

⚠⚠ THREE CONSTRAINTS, EACH LOAD-BEARING
---------------------------------------
1. **CIRCULARITY — β and `identity` PASS; `prior_pi`/`pip_discovery` are BANNED.** This ledger *is* the
   source of `w`, so scoring a w-informed quantity against it is circular. **Verified in code, not inherited**
   (`readouts.py:190-202`): `collapse_by_family` returns `wf`, `_prep` never reads it, and the dense chain
   runs `_gibbs_posterior(Xz, yr, np.ones(p))` — **π≡1** — while `_bagged_nnls_meansd(Xz, yr)` takes no w at
   all. So `beta` and `identity` are computed with **zero** reference to w (matching MH-126c's empirical gate:
   max|Δβ| = 0.0 under shuffled *and* constant w). The only w-dependent path is `_evidence_pi -> pip_disc`.
   ⇒ `attribution_rank.py` refuses to score `prior_pi` / `pip_discovery`, in code.
2. **THE UNIT IS THE FAMILY, NOT THE GENE.** Abundance is family-constant, so genes sharing a canonical
   family are not independent draws; the cluster-bootstrap unit is `lit_family` (MH-126c established this at
   13 families — the reason its n=32 was really an n=13).
3. **⚠ STUDY BIAS IS NOT REMOVED BY SCALING n.** PMID depth tracks abundance (MH-126c: spearman(abundance, w)
   within gene = +0.171 — well-studied miRNAs *are* abundant). A bigger set does not neutralise that; it makes
   **abundance the baseline to beat**, and abundance currently BEATS β. This set is expected to make the
   negative verdict SHARPER, not softer. It is built because it is the version that can fail.

ADMISSION CRITERIA (all four, per gene)
---------------------------------------
  (a) the gene is in the learned design (`gene_family_card.tsv`)
  (b) the top family by LT-func PMID depth is PRESENT among that gene's design families
  (c) the design offers >= `min_competitors` families — with one family there is nothing to attribute
  (d) the top is UNAMBIGUOUS: `margin` = n_pmid(top) - n_pmid(runner-up) >= 1, so ties are excluded rather
      than broken arbitrarily

TIERS (emitted together; the caller picks, and the trade-off stays visible)
    tier2 = >=2 PMIDs   tier3 = >=3   tier4 = >=4
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "lit_ground_truth.tsv"
PROV = OUT / "lit_ground_truth_provenance.json"
FAMILY_CARD = OUT / "gene_family_card.tsv"

#: low-throughput FUNCTIONAL assays — repression demonstrated, not binding observed.
LT_FUNC = ("reporter", "western", "proteomics")
MIN_COMPETITORS = 2


def _clean_arms(s: pd.Series) -> pd.Series:
    """Same source hygiene `ledger._clean_edge_arms` applies to the HE edge list — unicode dashes and
    comma-joined multi-arm cells are curation artifacts present in the RAW ledger too, and an uncleaned
    `hsa-miR‑15a` silently fails to map to a family (it becomes a phantom singleton)."""
    return (s.astype(str).str.strip().str.replace(r"[‐-―−]", "-", regex=True))


def build(*, min_pmid: int = 2, min_competitors: int = MIN_COMPETITORS) -> pd.DataFrame:
    """The pull. One row per admitted gene; `tier` records the deepest PMID threshold it satisfies."""
    from mirna_hallmark.learned import families as FAM
    from mirna_hallmark.learned.evidence import ledger as LG

    lg = LG.build_ledger()
    f = lg[lg["assay_class"].isin(LT_FUNC) & (~lg["weak"])].copy()
    f["arm"] = _clean_arms(f["arm"])
    f = f[f["arm"].str.startswith("hsa-")]
    fam = FAM.family_of(pd.Index(f["arm"].unique()))
    f["family"] = f["arm"].map(fam)
    unmapped = f["family"].isna().mean()
    f = f.dropna(subset=["family"])

    # family-level depth = DISTINCT PMIDs across the family's member arms (a paper on miR-17 and a paper on
    # miR-20a are two independent lines of evidence for the SAME family; the same paper on both is one).
    dep = f.groupby(["gene", "family"])["pmid"].nunique().reset_index(name="n_pmid")
    cls = (f.groupby(["gene", "family"])["assay_class"]
             .agg(lambda s: ",".join(sorted(set(s)))).reset_index(name="assay_classes"))
    pm = (f.groupby(["gene", "family"])["pmid"]
            .agg(lambda s: ",".join(str(int(x)) for x in sorted(set(s))[:12])).reset_index(name="pmids"))
    dep = dep.merge(cls, on=["gene", "family"]).merge(pm, on=["gene", "family"])

    card = pd.read_csv(FAMILY_CARD, sep="\t", usecols=["gene", "family"])
    design = card.groupby("gene")["family"].apply(set).to_dict()
    n_design = card.groupby("gene")["family"].nunique()

    rows, drop = [], {"gene_not_in_design": 0, "top_not_in_design": 0,
                      "too_few_competitors": 0, "ambiguous_tie": 0, "below_min_pmid": 0}
    for g, sub in dep.groupby("gene"):
        if g not in design:
            drop["gene_not_in_design"] += 1
            continue
        sub = sub.sort_values("n_pmid", ascending=False)
        top = sub.iloc[0]
        runner = int(sub["n_pmid"].iloc[1]) if len(sub) > 1 else 0
        if top["family"] not in design[g]:
            drop["top_not_in_design"] += 1
            continue
        if int(n_design[g]) < min_competitors:
            drop["too_few_competitors"] += 1
            continue
        if int(top["n_pmid"]) - runner < 1:
            drop["ambiguous_tie"] += 1
            continue
        if int(top["n_pmid"]) < min_pmid:
            drop["below_min_pmid"] += 1
            continue
        n = int(top["n_pmid"])
        rows.append({"gene": g, "lit_family": top["family"], "n_pmid": n,
                     "runner_up_pmid": runner, "margin": n - runner,
                     "n_lit_families": int(len(sub)), "n_design_families": int(n_design[g]),
                     "assay_classes": top["assay_classes"], "pmids": top["pmids"],
                     "tier": 4 if n >= 4 else (3 if n >= 3 else 2)})
    R = pd.DataFrame(rows).sort_values(["n_pmid", "gene"], ascending=[False, True]).reset_index(drop=True)
    R.attrs["dropped"] = drop
    R.attrs["unmapped_arm_frac"] = float(unmapped)
    return R


def _stamp(R: pd.DataFrame) -> dict:
    """Provenance: the two source files by size+mtime, the ledger cache, and a content hash of THIS table.
    Enough to answer 'which pull produced the number in that sentence?' without re-running anything."""
    from mirna_hallmark.learned.evidence import ledger as LG
    src = {}
    for name, p in (("mirtarbase", LG.MIRTAR), ("tarbase", LG.TARBASE), ("ledger_cache", LG.CACHE)):
        p = Path(p)
        src[name] = ({"path": str(p), "bytes": p.stat().st_size,
                      "mtime": pd.Timestamp(p.stat().st_mtime, unit="s").isoformat()}
                     if p.exists() else {"path": str(p), "missing": True})
    body = R.to_csv(sep="\t", index=False).encode()
    return {"rule": "argmax family by #distinct low-throughput-functional PMIDs; "
                    f"LT_FUNC={list(LT_FUNC)}; weak excluded",
            "n_genes": int(len(R)), "n_families": int(R["lit_family"].nunique()),
            "tiers": {f"tier{t}": int((R.tier >= t).sum()) for t in (2, 3, 4)},
            "dropped": R.attrs.get("dropped", {}),
            "unmapped_arm_frac": R.attrs.get("unmapped_arm_frac"),
            "sources": src, "sha256": hashlib.sha256(body).hexdigest()}


def load(min_tier: int = 2) -> pd.DataFrame:
    """Read the persisted set. `min_tier` selects PMID depth (2/3/4)."""
    R = pd.read_csv(DEST, sep="\t")
    return R[R["tier"] >= min_tier].reset_index(drop=True)


def report(R: pd.DataFrame) -> None:
    d = R.attrs.get("dropped", {})
    print(f"\n=== LITERATURE GROUND TRUTH — {len(R):,} genes / {R['lit_family'].nunique()} distinct "
          f"canonical families ===\n")
    print(f"  rule: argmax family by DISTINCT low-throughput functional PMIDs {LT_FUNC}, weak excluded")
    print(f"  admitted: in-design AND canonical-family-in-design AND >={MIN_COMPETITORS} competitors "
          f"AND unambiguous top\n")
    print("  --- TIERS (deeper = purer, smaller) ---")
    for t in (2, 3, 4):
        s = R[R.tier >= t]
        print(f"    tier{t} (>={t} PMIDs): {len(s):4d} genes  {s['lit_family'].nunique():3d} families  "
              f"median competitors {s.n_design_families.median():.0f}  median margin {s.margin.median():.0f}")
    print("\n  --- DROPPED (why a gene with literature did not make it) ---")
    for k, v in d.items():
        print(f"    {k:24s} {v:5d}")
    print(f"\n  --- FOR COMPARISON, the sets this replaces (all producer-less, 2026-07-13) ---")
    print(f"    MH-124 §4b n=21 · MH-126c TEST 1 n=32 (not on disk)  ->  this: "
          f"{int((R.tier>=3).sum())} at tier3, {len(R)} at tier2")
    print(f"\n  --- DEPTH ---")
    print(f"    n_pmid: median {R.n_pmid.median():.0f}  max {R.n_pmid.max():.0f}  "
          f"| assay mix: {R.assay_classes.str.contains('reporter').mean():.0%} reporter, "
          f"{R.assay_classes.str.contains('western').mean():.0%} western, "
          f"{R.assay_classes.str.contains('proteomics').mean():.0%} proteomics")
    top = R.head(8)[["gene", "lit_family", "n_pmid", "runner_up_pmid", "n_design_families"]]
    print(f"\n  --- deepest-evidence genes ---\n{top.to_string(index=False)}")
    print("\n  ⚠ STUDY BIAS IS NOT REMOVED BY SIZE — PMID depth tracks abundance (+0.171 within gene).")
    print("    This set makes ABUNDANCE the baseline to beat; it is expected to sharpen the negative verdict.")


def run(*, min_pmid: int = 2) -> pd.DataFrame:
    R = build(min_pmid=min_pmid)
    prov = _stamp(R)
    DEST.parent.mkdir(parents=True, exist_ok=True)
    R.to_csv(DEST, sep="\t", index=False)
    PROV.write_text(json.dumps(prov, indent=2))
    print(f"-> {DEST}  ({len(R):,} genes)\n-> {PROV}  (sha256 {prov['sha256'][:16]}…)")
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-pmid", type=int, default=2)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    R = run(min_pmid=a.min_pmid)
    report(R)
