"""Attribution-IDENTITY — the abundance-removed "who OWNS this gene" readout, sequence-grounded.

Distinct from three neighbours:
  * attribution-MAGNITUDE (who exerts the most force; abundance-IN) — the loud generalist wins.
  * IDENTIFIABILITY (can the DATA separate this edge from collinear ones) — a data-geometry property.
  * the shipped `structural_identity.structural_share` (evidence-mass fraction × potential × confidence) — same
    spirit but the evidence-mass fraction SATURATES to ~0 for hub arms with huge curated targetomes; this readout
    swaps in **seq-spec** (per-edge TargetScan affinity concentration), which keeps resolution there, and gates by
    **AGO loading** (a passenger arm cannot OWN a gene however specific its sequence).

    identity_raw(m,g) = seq_spec(m,g) · loading(m)                 # sequence dedication × guide-arm
    identity_share(m,g) = softmax_{m'∈covered(g)}( identity_raw / T )   # within-gene ownership distribution

Uncovered arms (no seq-spec) get share 0 (no identity claim). Non-circular (sequence + AGO measurement, never TCGA
X/Y, never coupling). Validated: recovers textbook specialists abundance-removed (miR-19/25→PTEN, miR-15/16/34→
CCND1, miR-200→ZEB1) while the loud generalists (miR-21, miR-17, miR-182) rank low. See
LEARNED_MODEL_DISCOVERY_SYNTHESIS §6a; [[ago-loading-arm-axis]].

CLI: .venv/bin/python3 -m mirna_hallmark.learned.attribution_identity PTEN CCND1 ZEB1
"""
from __future__ import annotations

import io
from contextlib import redirect_stdout

import numpy as np
import pandas as pd

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM
from mirna_hallmark.learned import ago_loading as AG
from mirna_hallmark.learned import seq_specificity as SQ
from mirna_hallmark.learned import chimeric_evidence as CE


def attribution_identity(gene: str, *, temperature: float = 1.0, w_prior_source: str = "ledger",
                         chim_sources=None, arm_level: bool = False, spec_source: str = "context++") -> pd.DataFrame:
    """TWO-LEVEL ownership of `gene`, matching the model's estimand (M is per seed-family; per-arm split = nomination).

    LEVEL 1 — FAMILY ownership (`ownership`): softmax over seed-families of the z-scored **affinity percentile**
    (is g among the family's strongest targets). This is the identity unit: co-seed arms target identically (shared
    seed → shared sites), so ownership is a SEED/family property — an arm-level ownership score would be spurious
    resolution. Percentile (not concentration-share) is scale-free and captures a hub-family-with-a-famous-target
    (miR-200 is a hub, but ZEB1 is its #1). Loading is NOT in the ownership score (it's an arm-level signal).

    LEVEL 2 — within-family ARM nomination (`nom_arm`): which specific mature carries the family = the member with
    the highest expression × AGO loading. Sequence/affinity cannot resolve this (co-seed = same seed); it needs the
    non-seed signals. (The empirical alternative: the member actually seen in chimeric reads — chimeric carries the
    full arm sequence, so it resolves the mature; a TODO hook.)"""
    _, X, _, w = LD.assemble_gene(gene, w_prior_source=w_prior_source)    # ARM-level
    arms = X.columns
    fam_map = FAM.family_of(arms)                                        # arm → TargetScan seed family
    Xf, wf, _ = FAM.collapse_by_family(X, w, fam_map)
    fam_cols = list(Xf.columns)
    fam_to_arms: dict = {}
    for a in arms:
        f = fam_map.get(a) if hasattr(fam_map, "get") else fam_map[a]
        fam_to_arms.setdefault(f, []).append(a)
    # LEVEL 1: family ownership (seed-pooled affinity percentile; loading excluded).
    # spec_source='kd' → validated per-arm GENOME-WIDE scanMiR K_D (kd_fair_bench: recovers specialists 0.89 vs
    # context++ 0.79), aggregated to family = strongest member (the "Family = strongest member" convention).
    if spec_source == "kd":
        kd_arm = SQ.affinity_percentile_kd(gene, list(arms))
        ded = pd.Series({f: (np.nanmax(v) if len(v := [kd_arm.get(a, np.nan) for a in fam_to_arms.get(f, [])])
                                              and not np.all(np.isnan(v)) else np.nan) for f in fam_cols}).reindex(fam_cols)
    else:
        ded = SQ.affinity_percentile(gene, pd.Index(fam_cols))
    dv = ded.to_numpy(dtype=float); cov = ~np.isnan(dv)
    own = np.zeros(len(fam_cols))
    if cov.sum() > 1 and dv[cov].std() > 1e-12:
        z = np.where(cov, (dv - dv[cov].mean()) / dv[cov].std() / max(temperature, 1e-6), -np.inf)
        z = z - np.nanmax(z[cov]); e = np.where(cov, np.exp(z), 0.0)
        own = e / e.sum() if e.sum() > 0 else e
    # LEVEL 2: within-family arm DISTRIBUTION (a family's action can be split across several expressed members,
    # not winner-take-all). Weights = POOLED chimeric reads (empirical, arm-resolved) where present, else
    # expression×loading. Reported as "arm:share" list + the chimeric source flag.
    arm_load = AG.loading(arms).reindex(arms); arm_expr = X.median()
    arm_dist, nom_src, fam_share = {}, {}, {}
    for f in fam_cols:
        mem = fam_to_arms.get(f, [])
        if not mem:
            arm_dist[f] = "-"; fam_share[f] = pd.Series(dtype=float); continue
        chim = CE.evidence(gene, mem, sources=chim_sources)             # pooled chimeric weight per member
        if chim.sum() > 0:
            wv = chim; nom_src[f] = "chimera"
        else:
            wv = pd.Series({m: float(arm_expr.get(m, 0.0)) * float(arm_load.get(m, 1.0)) for m in mem})
            nom_src[f] = "expr×load"
        wv = wv[wv > 0]
        share = (wv / wv.sum()).sort_values(ascending=False) if wv.sum() > 0 else wv
        fam_share[f] = share                                            # L2 within-family distribution (raw)
        arm_dist[f] = ", ".join(f"{m.replace('hsa-','')}:{s:.2f}" for m, s in share.head(3).items()) or "-"
    own_by_fam = dict(zip(fam_cols, own))
    if arm_level:
        # HIERARCHICAL COMBINE: arm_ownership = ownership(family) × share(arm|family) → joint over ALL arms
        rows = []
        for f in fam_cols:
            sh = fam_share[f]
            if sh.empty:                                               # no resolvable member → family keeps its own weight
                rows.append((f.split("/")[0], f.split("/")[0], round(own_by_fam[f], 4),
                             round(own_by_fam[f], 3), 1.0, nom_src.get(f, "unresolved")))
                continue
            for m, s in sh.items():
                rows.append((m.replace("hsa-", ""), f.split("/")[0], round(own_by_fam[f] * s, 4),
                             round(own_by_fam[f], 3), round(s, 3), nom_src.get(f, "-")))
        arm = pd.DataFrame(rows, columns=["arm", "family", "arm_ownership", "fam_ownership", "within_fam", "src"])
        return arm.sort_values("arm_ownership", ascending=False).set_index("arm")
    out = pd.DataFrame({
        "ownership": np.round(own, 3),                                   # L1 family ownership (the identity)
        "aff_pctile": ded.round(3).to_numpy(),                          # family seed-level dedication
        "fam_abund": Xf.median().round(2).to_numpy(),                    # MAGNITUDE contrast
        "arm_dist": [arm_dist.get(f, "-") for f in fam_cols],          # L2 DISTRIBUTION over members (arm:share)
        "nom_src": [nom_src.get(f, "-") for f in fam_cols],
    }, index=[f.split("/")[0] for f in fam_cols])
    return out.sort_values("ownership", ascending=False)


def _report(gene: str):
    d = attribution_identity(gene)
    print(f"\n=== {gene}: L1 FAMILY ownership + L2 arm nomination (abundance-removed) ===")
    print(d.head(6).to_string())
    top_id = d.index[0]
    top_mag = d.sort_values("fam_abund", ascending=False).index[0]
    print(f"  owner family: {top_id} → arms [{d.loc[top_id,'arm_dist']}] (own {d.loc[top_id,'ownership']}) "
          f"| loudest family: {top_mag} (own {d.loc[top_mag,'ownership']})")
    from scipy.stats import spearmanr
    cov = d["aff_pctile"].notna()
    if cov.sum() > 4:
        r = spearmanr(d.loc[cov, "ownership"], d.loc[cov, "fam_abund"]).correlation
        print(f"  Spearman(ownership, fam_abund) = {r:+.3f}  (want ~0 = abundance-removed)")


if __name__ == "__main__":
    import sys
    for g in (sys.argv[1:] or ["PTEN", "CCND1", "ZEB1"]):
        _report(g)
