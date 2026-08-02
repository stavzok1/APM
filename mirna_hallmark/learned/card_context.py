"""DESIGN & COMPETENCE CONTEXT for the gene and edge cards (user-asked, 2026-08-01).

    .venv/bin/python3 -m mirna_hallmark.learned.card_context              # build context tables only
    .venv/bin/python3 -m mirna_hallmark.learned.card_context --annotate   # + join onto the cards

WHY. The cards answer "what does this edge/gene do". They could not answer the PRIOR question —
**"can this gene's regulation be measured at all, and is its real-vs-fake advantage real here?"** —
because the competence lens (`gene_atlas`), the decoy bench (`decoy_bench`) and the CPTAC protein
retention (`dossier` gold) all lived in separate artifacts. An edge sitting in a gene whose OOF ceiling
is <=0 is **uninterpretable**, and nothing on the card said so. This block puts that on the card.

WHAT IT ADDS — all `ctx_` prefixed, all A PRIORI or already-measured elsewhere, nothing re-estimated here.

GENE level (join key `gene`)
    ctx_ceiling            5-fold OOF R^2 of an UNRESTRICTED OLS on the gene's real design = the upper
                           bound for ANY linear weighting of its regulators (`gene_atlas`).
    ctx_measurable         ⚠ gated at **> 0.02, NOT > 0** — MH-144: `ceiling <= 0` is a threshold count
                           evaluated exactly where the mass piles up (42.8% of genes sit within +/-0.01),
                           so a 0.06%-of-scale shift reclassifies ~8% of the universe. 0.02 reproduces.
    ctx_apriori_class      A_COMPETENT / B_PARTIAL / C_WEAK / D_UNDEFINED. ⚠ MH-147/MH-169: the evidence
                           axis is UNDERPOWERED given width — treat this as `n_fam >= 3` plus noise.
    ctx_gap_core/_deconv   the gene's real-vs-fake decoy gap (`decoy_bench`), core and composition blocks.
    ctx_gap_retention      gap_deconv/gap_core, **GATED** at |gap_core| >= 0.005 (axiom 5) else NaN.
    ctx_n_abund            real regulator arms above the canonical floor log2(11) (RPM>=10).
    ctx_n_fam_abund        seed families carrying at least one abundant arm.
    ctx_collapse           n_arm/n_fam — >1 means curation bundles same-seed mates onto this gene.
    ctx_d_fam              fake minus real design width (MH-167; real but gap-INERT, kept for audit).

EDGE level (join key `gene`,`arm`) — the gene block is broadcast, plus:
    ctx_arm_dose           this arm's median log2(RPM+1).
    ctx_arm_abundant       ctx_arm_dose >= log2(11).
    ctx_protein_*          per-edge CPTAC protein coupling on the INDEPENDENT prospective cohort, both
                           composition-ADJUSTED and the historical RAW, plus their gated retention and
                           class. Gold/cross-candidacy edges only (`discovery_dossier_gold`); NaN elsewhere.

⚠⚠ READ THIS BEFORE USING `ctx_gap_*` PER GENE. The decoy gap is a **SET-LEVEL** quantity. A single gene's
gap is one noisy draw — MH-123/124 measured the per-edge null as 3-4x too narrow, and MH-137/139's headline
is a mean over ~1,350 genes (−0.0129). **Use these columns to STRATIFY and to EXCLUDE, never to call an
individual gene "validated".** Same for `ctx_protein_*`: MH-171 measured median retention 0.414 on gold
edges with **241/495 `composition_explained`**, so a raw protein coupling is not evidence on its own.

⚠ ADDITIVE ONLY. `--annotate` adds `ctx_` columns and rewrites the cards; every pre-existing column is
verified BIT-IDENTICAL first and the write is aborted otherwise (axiom 2 — the cards have 3 consumers:
`canonical_card.py`, `discovery.py`, `realization.py`).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
GENE_CTX = OUT / "gene_design_context.tsv"
EDGE_CTX = OUT / "edge_design_context.tsv"
# ⭐ THERE IS EXACTLY ONE EDGE CARD (consolidated 2026-08-01, user-directed).
# `edge_card.tsv` is a strict SUPERSET of `edge_card_base.tsv` — measured: canonical has
# **0** columns the progression card lacks, all 108 are shared, and the progression card adds 14
# within-patient realization columns (MH-158). Annotating both duplicated ~108 columns across two files
# and gave two answers to "where is the edge card".
# `edge_card_base.tsv` is NOT a deliverable — it is the BUILD INTERMEDIATE that `realization.py` reads
# (3 call sites) to construct the edge card. It is deliberately left UNANNOTATED so there is one place
# to look. ⚠ Its filename still says "canonical", which is now misleading; renaming it (and the edge
# card to `edge_card.tsv`) is a clean follow-up.
# ⭐ RENAMED 2026-08-01 (MH-181): the artifacts are now `edge_card.tsv` / `gene_card.tsv` /
# `edge_card_base.tsv` — previously `progression_edge_card` / `progression_gene_card` / `canonical_card`,
# names that mislabelled the complete card as progression-specific and the intermediate as "canonical".
# ⚠ The old names REMAIN in `ANALYSIS_RUN_LEDGER` / `DISCOVERY_REGISTRY` rows ON PURPOSE: those are
# HISTORICAL records of what ran under what name, and rewriting them would falsify the record.
CARDS_EDGE = [OUT / "realization/edge_card.tsv"]
CARDS_GENE = [OUT / "realization/gene_card.tsv"]
INTERMEDIATE = OUT / "edge_card_base.tsv"

_FLOOR = np.log2(11)          # RPM>=10, the canonical arm-expression floor
GAP_GATE = 0.005              # axiom 5: below this the retention ratio is a coin-flip with decimals
RHO_GATE = 0.05               # same, for the protein retention denominator
CEILING_ROBUST = 0.02         # MH-144: NOT 0 — the mass piles up at 0
# every prefix this module writes onto a card; used to strip-before-rejoin AND to select
ANNOT_PREFIXES = ("ctx_", "cptac_", "tcga_", "comp_", "cal_")
# ⛔ EXPLICIT NAMES, NOT A PREFIX, for the arm-rung block. Adding `"arm_"` to ANNOT_PREFIXES on
# 2026-08-01 silently STRIPPED 8 pre-existing card columns (`arm_lfc_*`, `arm_med_rpm`, `arm_pct_floor`,
# `arm_iqr`, `arm_id_status`, …) — the card went 104 -> 96 base columns. The additivity check could NOT
# catch it: it compares against `before`, which had already been stripped. A prefix is a blunt instrument
# on a 140-column table; name what you own.
ARM_RUNG_COLS = ("n_arm_in_cell", "beta_arm", "sd_arm", "z_arm", "arm_dbeta", "arm_sep_z",
                 "oof_rho_arm", "oof_rho_fam", "oof_drho", "arm_resolvable", "coupling_fam")

# ⭐ POSTERIOR-WIDTH CALIBRATION CONSTANT (MH-185, 2026-08-01). `calibration.posterior_calibration`
# compares the reported posterior SD against the TRUE sampling SD from INDEPENDENT half-cohorts.
# MEASURED over 4 seeds (200 genes, 3 splits): Gibbs Gaussian ratio **0.815 ± 0.084** ⇒ the shipped
# widths are ~1.23x too NARROW. ⚠⚠ THE CONSTANT IS ITSELF IMPRECISE — the seed SD is 0.084 (range
# 0.74–0.93), so the inflation is 1.23 with a plausible range ~1.07–1.35. It is emitted WITH that
# range rather than as a bare multiplier, per axiom 5.
# ⛔ Do NOT swap in Student-t instead: ν=7 reaches only 0.849 ± 0.051 — better, still not calibrated.
# ⛔ And do NOT model a subtype random effect for this: subtype-STRATIFIED half-splits were measured
# and did NOT reduce the true sampling SD (Δ −0.026 ± 0.051, p=0.63) ⇒ the heterogeneity is finer
# grained than PAM50, so the honest fix is REPORTING, which is what these columns are.
WIDTH_RATIO, WIDTH_RATIO_SD = 0.815, 0.084


def _read(p: Path, **kw):
    return pd.read_csv(p, sep="\t", **kw) if p.exists() else None


def gene_context() -> pd.DataFrame:
    atlas = _read(OUT / "gene_competence_map.tsv")
    bench = _read(OUT / "decoy_bench.tsv")
    prof = _read(OUT / "decoy_design_profile.tsv")
    if atlas is None:
        raise FileNotFoundError("gene_competence_map.tsv missing — run learned.analyses.gene_atlas")

    g = atlas[["gene", "n_fam", "n_arms", "w_max", "ceiling_R2_oof_fam_core", "apriori_class"]].copy()
    g = g.rename(columns={"ceiling_R2_oof_fam_core": "ctx_ceiling", "apriori_class": "ctx_apriori_class"})
    g["ctx_measurable"] = g.ctx_ceiling > CEILING_ROBUST          # ⚠ 0.02, not 0 (MH-144)

    if bench is not None:
        b = bench[["gene", "gap_core", "gap_deconv", "d_dose"]].copy()
        b = b.rename(columns={"gap_core": "ctx_gap_core", "gap_deconv": "ctx_gap_deconv",
                              "d_dose": "ctx_gap_d_dose"})
        # ⭐ GATED RATIO (axiom 5) — a retention over a vanishing gap is not a number
        b["ctx_gap_retention"] = np.where(b.ctx_gap_core.abs() >= GAP_GATE,
                                          b.ctx_gap_deconv / b.ctx_gap_core, np.nan)
        if {"n_fam_fake", "collin_real", "collin_fake"} <= set(bench.columns):
            b["ctx_d_fam"] = bench.n_fam_fake - bench.n_fam
            b["ctx_d_collin"] = bench.collin_fake - bench.collin_real
        g = g.merge(b, on="gene", how="left")

    if prof is not None:
        p = prof[["gene", "n_abund", "n_fam_abund", "frac_abund", "dose_max", "dose_med",
                  "collapse", "n_fam_multi"]].copy()
        p.columns = ["gene"] + [f"ctx_{c}" for c in p.columns[1:]]
        g = g.merge(p, on="gene", how="left")

    g.to_csv(GENE_CTX, sep="\t", index=False)
    print(f"[card_context] gene context: {len(g):,} genes, {g.shape[1]-1} cols -> {GENE_CTX.name}")
    return g


def edge_context() -> pd.DataFrame:
    """Per-(gene, arm): arm abundance + the gold CPTAC protein retention. Gene block is joined at annotate."""
    from mirna_hallmark.learned import data as LD

    X = LD._load()["X"]
    dose = pd.Series(np.nanmedian(X.to_numpy(float), axis=1), index=X.index, name="ctx_arm_dose")

    gold = _read(OUT / "discovery_dossier_gold.tsv")
    if gold is not None:
        e = gold[["gene", "arm", "protein_coupling", "protein_coupling_RAW", "prot_n"]].copy()
        e = e.rename(columns={"protein_coupling": "ctx_protein_rho",
                              "protein_coupling_RAW": "ctx_protein_rho_RAW", "prot_n": "ctx_protein_n"})
        # ⭐ GATED (axiom 5) + the FLAG-don't-delete class (axiom 2a)
        e["ctx_protein_retention"] = np.where(e.ctx_protein_rho_RAW.abs() >= RHO_GATE,
                                              e.ctx_protein_rho / e.ctx_protein_rho_RAW, np.nan)
        r = e.ctx_protein_retention
        cls = pd.Series(pd.NA, index=e.index, dtype="object")   # object dtype: NA and str coexist
        cls[r >= 0.7] = "cell_intrinsic"
        cls[(r >= 0.4) & (r < 0.7)] = "partial"
        cls[r < 0.4] = "composition_explained"
        e["ctx_protein_class"] = cls
        e = e.drop_duplicates(subset=["gene", "arm"])
    else:
        e = pd.DataFrame(columns=["gene", "arm"])

    # ⚠ the ARM block is kept SEPARATE and keyed on `arm` alone. Folding it into the gold-edge frame
    # (an outer merge keyed on gene+arm) would have restricted `ctx_arm_dose` to the 747 GOLD edges —
    # it must reach all ~5,650 card edges, most of which have no gold row.
    arm = dose.rename_axis("arm").reset_index()
    arm["ctx_arm_abundant"] = arm.ctx_arm_dose >= _FLOOR

    e.to_csv(EDGE_CTX, sep="\t", index=False)
    print(f"[card_context] edge context: {len(e):,} gold edges "
          f"({e.ctx_protein_rho.notna().sum() if 'ctx_protein_rho' in e else 0} with protein rho) "
          f"+ {len(arm):,} arms with dose -> {EDGE_CTX.name}")
    return e, arm


def _annotate(card_path: Path, blocks: list, keys: list) -> None:
    """Join `blocks` onto a card. ADDITIVE ONLY — aborts if any pre-existing column would change."""
    if not card_path.exists():
        print(f"  ⚠ absent, skipped: {card_path}")
        return
    card = pd.read_csv(card_path, sep="\t")
    # ⚠ strip EVERY annotated prefix, not just ctx_ — otherwise a re-run appends `_x`/`_y` duplicates
    # of the CPTAC block instead of replacing it, and the additivity control below would pass anyway.
    _owned = set(ARM_RUNG_COLS)
    before = card.drop(columns=[c for c in card.columns
                                if c.startswith(ANNOT_PREFIXES) or c in _owned], errors="ignore")
    out = before.copy()
    for blk, key in zip(blocks, keys):
        cols = [c for c in blk.columns if c.startswith(ANNOT_PREFIXES) or c in _owned or c in key]
        out = out.merge(blk[cols].drop_duplicates(subset=key), on=key, how="left")
    # ⭐ THE CONTROL: every pre-existing column must be bit-identical, in the same order
    chk = out[before.columns]
    if not chk.equals(before):
        diff = [c for c in before.columns if not chk[c].equals(before[c])]
        raise SystemExit(f"⛔ ABORT — annotation is NOT additive for {card_path.name}: {diff[:8]}")
    if len(out) != len(before):
        raise SystemExit(f"⛔ ABORT — row count changed for {card_path.name}: {len(before)} -> {len(out)}")
    out.to_csv(card_path, sep="\t", index=False)
    added = [c for c in out.columns if c.startswith(ANNOT_PREFIXES) or c in _owned]
    print(f"  ✅ {card_path.name}: {len(before.columns)} -> {out.shape[1]} cols "
          f"(+{len(added)} annotated), {len(out):,} rows, pre-existing columns bit-identical")


def _calibrate_widths(card_path: Path) -> None:
    """Emit HONEST posterior widths next to the raw ones — the model is overconfident by a MEASURED factor.

    `beta_sd` is ~1.23x too narrow (MH-185), so every `z = beta/sd` is ~1.23x too large and the
    `identified` (|z|>2) count is inflated. Rather than silently rescale the shipped column, this adds a
    parallel set so a reader can see both and the correction is auditable:
        cal_beta_sd        beta_sd / WIDTH_RATIO                 the honest width
        cal_z              beta / cal_beta_sd                    the honest z
        cal_identified     |cal_z| > 2                           the honest identifiability call
        cal_beta_sd_lo/hi  the width under the constant's OWN uncertainty (ratio +/- 1 seed SD)
        cal_identified_lo/hi   identifiability at those two ends — the honest RANGE, not a point
    ⚠ This does NOT make the model calibrated; it makes its overconfidence explicit and correctable.
    """
    if not card_path.exists():
        return
    d = pd.read_csv(card_path, sep="\t", low_memory=False)
    if "beta_sd" not in d.columns or "beta" not in d.columns:
        return
    raw_ident = d["identified"].mean() if "identified" in d else np.nan
    for tag, r in (("", WIDTH_RATIO), ("_lo", WIDTH_RATIO + WIDTH_RATIO_SD),
                   ("_hi", WIDTH_RATIO - WIDTH_RATIO_SD)):     # lower ratio => WIDER sd => fewer identified
        sd = d["beta_sd"] / r
        z = d["beta"] / sd.replace(0, np.nan)
        if tag == "":
            d["cal_beta_sd"], d["cal_z"] = sd, z
        else:
            d[f"cal_beta_sd{tag}"] = sd
        d[f"cal_identified{tag}"] = (z.abs() > 2)
    d.to_csv(card_path, sep="\t", index=False)
    print(f"     ↳ calibrated widths: identified {raw_ident:.1%} (raw) -> "
          f"{d.cal_identified.mean():.1%} at ratio {WIDTH_RATIO:.3f} "
          f"[{d.cal_identified_hi.mean():.1%}-{d.cal_identified_lo.mean():.1%}] over its +/-1SD range")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--annotate", action="store_true", help="join the context onto the cards in place")
    a = ap.parse_args()
    g = gene_context()
    e, arm = edge_context()
    if not a.annotate:
        print("\n(context tables written; pass --annotate to join them onto the cards)")
        return
    # ⛔ THE GOLD PROTEIN BLOCK IS DELIBERATELY *NOT* JOINED ONTO THE CARDS (measured 2026-08-01).
    # `discovery_dossier_gold` scores ORPHAN edges; the cards are the HE universe. They are DISJOINT
    # BY CONSTRUCTION — measured: |card ∩ gold| = **0** of 5,648 vs 747, even though 466 genes and 115
    # arms are shared. This is MH-142's structural fact ("0 of dossier's 6,744 candidate pairs are HE
    # edges"), not a join bug: joining anyway shipped 5 columns that were 0.0% populated, which reads as
    # "no protein data for this edge" when the truth is "this edge is not in that lane's universe."
    # The per-edge protein retention stays in `edge_design_context.tsv`, where the ORPHAN lane can use it.
    # ⭐ CPTAC block — the same coupling question in an INDEPENDENT cohort, plus the protein and
    # protein-vs-mRNA layers TCGA cannot supply (`learned/cptac_card.py`). Joined on (gene, arm); unlike
    # the gold-protein block above this one IS in the HE universe, so it actually populates (79% / 70%).
    cptac = _read(OUT / "edge_cptac_card.tsv")
    blocks, keys = [g, arm], [["gene"], ["arm"]]
    # ⭐ ARM RUNG (MH-186): an arm-level β + an IDENTIFIABILITY verdict, on the 20% of edges where a
    # family contributes >1 arm. The family β stays the CANONICAL default — these columns say WHERE the
    # equal-β broadcast is an assumption you can check, and whether the arms are separable at all.
    # ⭐ FAMILY-rung coupling (MH-187) — the card's `coupling_*` is single-ARM while `beta` is FAMILY.
    # `coupling_fam` puts the family rung beside it. ⛔ `identity_arm` from the same module is
    # DELIBERATELY NOT JOINED: it FAILED its pre-registered coherence check against `arm_resolvable`
    # (concentration was HIGHER in non-resolvable cells, MWU p=1) and the failure is unexplained
    # (abundance span does not account for it, ρ=+0.018 p=0.74). It stays in `rung_parity.tsv` only.
    rp = _read(OUT / "rung_parity.tsv")
    if rp is not None and "coupling_fam" in rp.columns:
        blocks.append(rp[["gene", "arm", "coupling_fam"]]); keys.append(["gene", "arm"])
    armr = _read(OUT / "arm_rung.tsv")
    if armr is not None:
        armr = armr.drop(columns=[c for c in ("seed_family",) if c in armr.columns])
        blocks.append(armr); keys.append(["gene", "arm"])
    if cptac is not None:
        cptac = cptac.rename(columns={c: c for c in cptac.columns})
        blocks.append(cptac); keys.append(["gene", "arm"])
    print("\n[card_context] annotating EDGE cards (gene block + arm block"
          f"{' + CPTAC block' if cptac is not None else ''}):")
    for p in CARDS_EDGE:
        _annotate(p, blocks, keys)
        _calibrate_widths(p)
    # ⭐ GENE-level CPTAC AGGREGATE (`cptac_card.build_gene`): the gene's whole miRNA budget —
    # SUM_a beta_a*X_a with beta the DENSE GIBBS mean fit on TCGA mRNA at n~1000 — scored against each
    # CPTAC layer, alongside an UNWEIGHTED sum-abundance reference and a TCGA reference built the same way.
    gcp = _read(OUT / "gene_cptac_card.tsv")
    gblocks, gkeys = [g], [["gene"]]
    if gcp is not None:
        gblocks.append(gcp); gkeys.append(["gene"])
    # ⭐ COMPARTMENT block (`learned/compartment_card.py`): WHICH compartment drives a gene's apparent
    # coupling — the CAUSE, where `retention`/`composition_class` only carried the CONSEQUENCE.
    cmp_ = _read(OUT / "compartment_summary.tsv")
    if cmp_ is not None:
        gblocks.append(cmp_); gkeys.append(["gene"])
    print(f"[card_context] annotating GENE cards"
          f"{' (design block + CPTAC aggregate block)' if gcp is not None else ''}:")
    for p in CARDS_GENE:
        _annotate(p, gblocks, gkeys)


if __name__ == "__main__":
    main()
