"""⭐ THE SEED-FAMILY CARD — the FIFTH rung, one row per seed family (MH-212).

    .venv/bin/python3 -m mirna_hallmark.learned.seed_family_card

**WHY THIS EXISTS — the same gap as the arm card, one rung up (user-identified, 2026-08-03).**
`card_rungs.CARDS["family"]` is keyed **`['gene', 'family']`** — it is a *gene × family* card (5,117
rows over 592 families), which can say "how family F behaves at gene G" but **cannot say anything about
family F itself**. So per-FAMILY facts — how concentrated the family's dose is, whether its members are
seed-homogeneous, how redundant its members are as tests — had nowhere to live, exactly as per-ARM facts
had nowhere to live before MH-209. `family_context.tsv` (573 rows) is the only one-row-per-family
artifact in the subproject and it was never registered as a card.

⛔ **NAME COLLISION AVOIDED DELIBERATELY:** the artifact is `seed_family_card.tsv`, NOT `family_card.tsv`
— that name is already taken by the (gene, family) card. Two different rungs must never share a filename;
`--check`'s cross-card report exists because that class of collision is this subproject's most common
false bug (axiom 6).

## What a FAMILY row is for

A seed family is ONE targeting unit by sequence, but its members are wildly unequal by dose. Every
column here is a property of the unit; the *member's* position inside it is the arm card's `famrole_*`.
The pair is the point: `bc_fam_*` on the arm card broadcasts these DOWN, and `famrole_*` says which
member the family actually is.

⚠ **`fam_*` counts and dose statistics are aggregates over the family's MEMBER ARMS restricted to the
arm card's universe** — a family whose members are unexpressed in TCGA will have a small `fam_n_members`
here and that is a coverage fact, not biology. `fam_n_members_context` is the independent count from
`genomic_context.family_context()`; compare them before reading either.
⚠ **Dose pooling is the TRUE RPM pool** `log2(1 + Σ(2^x − 1))`, never a mean of logs.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "seed_family_card.tsv"
ARM_CARD = OUT / "arm_card.tsv"


def _hhi(v) -> float:
    """Concentration, FLOOR-CORRECTED. Raw HHI is bounded below by 1/k so it falls as k grows whatever
    the data does — that cost MH-201 a wrong claim (−0.667 raw vs −0.075 corrected)."""
    v = np.abs(np.asarray(v, float))
    v = v[np.isfinite(v)]
    k = len(v)
    if k < 2 or v.sum() <= 0:
        return np.nan
    h = float(((v / v.sum()) ** 2).sum())
    return float((h - 1 / k) / (1 - 1 / k))


def build() -> pd.DataFrame:
    if not ARM_CARD.exists():
        raise SystemExit("build the arm card first: python -m mirna_hallmark.learned.arm_card")
    a = pd.read_csv(ARM_CARD, sep="\t", low_memory=False)
    a = a[a["seed_family"].notna()].copy()
    lin = (np.power(2.0, pd.to_numeric(a["abund_med"], errors="coerce")) - 1.0).clip(lower=0)
    a["_lin"] = lin
    g = a.groupby("seed_family")

    out = pd.DataFrame({
        # ── STRUCTURE
        "fam_n_members": g.size(),
        "fam_n_expressed": g["cov_expr"].sum(),
        # ── DOSE: level, concentration, and WHO carries it
        "fam_dose_total_log2": np.log2(1.0 + g["_lin"].sum()),          # TRUE RPM pool
        "fam_dose_med_log2": g["abund_med"].median(),
        "fam_dose_max_log2": g["abund_med"].max(),
        "fam_dose_hhi": g["_lin"].apply(_hhi),                          # floor-corrected
        # ⚠ `idxmax` is NaN for a family whose members all lack abundance — indexing with it raises.
        "fam_dominant_arm": a.loc[g["_lin"].idxmax().dropna(), ["seed_family", "arm"]]
                             .set_index("seed_family")["arm"],
        "fam_dominant_share": g["_lin"].apply(lambda s: float(s.max() / s.sum()) if s.sum() > 0 else np.nan),
        # ── VARIANCE: dynamic range is what predicts, MEAN level is not (axiom 8)
        "fam_var_med": g["abund_sd"].median(),
        "fam_var_max": g["abund_sd"].max(),
        "fam_var_hhi": g["abund_sd"].apply(_hhi),
    })

    # ── SEQUENCE HETEROGENEITY — is the collapse even valid for this family?
    if "iso_seed_shift_frac" in a.columns:
        # ⚠ gate on n_samples: the seed-shift median is n-DEPENDENT (0.056 raw vs 0.037 at n>=300)
        w = a[pd.to_numeric(a.get("iso_n_samples"), errors="coerce") >= 300]
        if len(w):
            out["fam_seed_shift_max_gated"] = w.groupby("seed_family")["iso_seed_shift_frac"].max()
            out["fam_n_shift_gt20_gated"] = (w[w.iso_seed_shift_frac > 0.20]
                                             .groupby("seed_family").size())
            out["fam_n_shift_gt20_gated"] = out["fam_n_shift_gt20_gated"].fillna(0)

    # ── TARGETOME + EVIDENCE, aggregated over members
    for src, dst, how in (("seq_n_genes_strong", "fam_targetome_seq_med", "median"),
                          ("sdsz_N_eff", "fam_targetome_seed_med", "median"),
                          ("fame_npmid", "fam_fame_npmid_sum", "sum"),
                          ("cur_he_degree_expr", "fam_cur_degree_max", "max"),
                          ("real_n_c10", "fam_real_n_c10_sum", "sum"),
                          ("real_n_scored", "fam_real_n_scored_sum", "sum"),
                          ("arb_n_edges", "fam_n_model_edges", "sum")):
        if src in a.columns:
            out[dst] = getattr(g[src], how)()

    # ── AGO / guide structure
    if "ago_dom" in a.columns:
        out["fam_n_guide"] = a[a.get("ago_guide_class") == "guide"].groupby("seed_family").size()
        out["fam_n_passenger"] = a[a.get("ago_guide_class") == "passenger"].groupby("seed_family").size()
        for c in ("fam_n_guide", "fam_n_passenger"):
            out[c] = out[c].fillna(0)

    # ── LOCUS homogeneity (the independent per-family artifact)
    p = OUT / "family_context.tsv"
    if p.exists():
        fc = pd.read_csv(p, sep="\t").set_index("family")
        keep = [c for c in ("n_members", "context_homogeneous", "dose_coding_host_frac",
                            "dose_cotranscribed_frac", "n_distinct_coding_hosts", "composition")
                if c in fc.columns]
        fc = fc[keep].rename(columns={c: f"fam_ctx_{c}" for c in keep})
        fc = fc.rename(columns={"fam_ctx_n_members": "fam_n_members_context"})
        out = out.join(fc, how="left")

    # ── REDUNDANCY: how many INDEPENDENT tests the family's members really are
    p = C.REPO_ROOT / "mirna_hallmark/output/coupling_inference/seed_family_justification/per_family_redundancy.tsv"
    if p.exists():
        r = pd.read_csv(p, sep="\t").set_index("seed_family")
        keep = [c for c in ("m_eff", "overcount_frac", "median_within_rho") if c in r.columns]
        out = out.join(r[keep].rename(columns={c: f"fam_redund_{c}" for c in keep}), how="left")

    # ⚠ DEGENERACY FLAG — with ONE member every share/rank/HHI is 1 or NaN BY CONSTRUCTION. Split on
    # this before pooling anything (gene_axes.mask_degenerate's rule, at the family rung).
    out["fam_single_member"] = out["fam_n_members"] <= 1

    out = out.rename_axis("seed_family").reset_index()
    out.to_csv(DEST, sep="\t", index=False)
    print(f"[seed_family_card] {len(out):,} families × {out.shape[1]} columns -> {DEST}")
    print(f"  single-member (degenerate): {int(out.fam_single_member.sum()):,} "
          f"({100*out.fam_single_member.mean():.1f}%)")
    multi = out[~out.fam_single_member]
    print(f"  multi-member: {len(multi):,} | median members {multi.fam_n_members.median():.0f} | "
          f"median dominant share {multi.fam_dominant_share.median():.3f}")
    return out


if __name__ == "__main__":
    build()
