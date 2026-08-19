"""Cross-card LADDERS — closes the three rung gaps the MH-215 cross-card audit found.

    .venv/bin/python3 -m mirna_hallmark.learned.card_ladders            # report only
    .venv/bin/python3 -m mirna_hallmark.learned.card_ladders --annotate # + join onto the cards

**WHY.** MH-215 put a realized-coupling threshold LADDER on the arm card ("how many genes does this arm
target at coupling < −0.10?"). The audit then found the mirror question has no home:

  * ⬜ **GENE side** — the gene card has 7 realization columns but **no ladder**: nothing anywhere answers
    *"how many of THIS GENE's regulators realize at depth X"*. Added here as `greal_`.
  * ⬜ **EDGE side** — PAM50 subtype heterogeneity is natively per-(gene, arm) yet lived only in the side
    artifact `subtype_edge_heterogeneity.tsv`, never on a card. Added here as `esub_` — NOT `sub_`,
    which the arm card already owns for its per-arm aggregate of the same source (domains match by
    prefix globally, so sharing it would misfile the edge block).
  * ✅ **AGO fabrication — FIXED, and fixed WITHOUT changing model behaviour (MH-214).** The edge card's
    `ago_loading` carries `ago_loading.loading()`'s `.fillna(1.0)`: a measured guide (1.0) is
    indistinguishable from an arm AGO never measured (also 1.0), for **~53% of the card's arms**.
    ⚠ The fill is NOT simply wrong — `loading()` has two consumers with opposite needs.
    `analyses/attribution_identity` MULTIPLIES by it, where "no discount for an unmeasured gate" is the
    defensible default; `readouts` PERSISTS it as a label, where the same default destroys information.
    ⇒ changing the return value would silently alter a model computation, so instead the fabrication is
    made **visible**: `ago_loading.measured()` (new) reports who was really measured, and `esub_`'s
    sibling column `ago_loading_measured` is annotated onto the edge card here.

⚠ **`ago_loading_measured` deliberately carries no module prefix** — it must sit adjacent to the column
it qualifies. That needs no special handling: `_annotate` derives what it owns from the blocks, so an
un-prefixed column is stripped and re-added exactly like any other.

⚠ **ADDITIVE ONLY, and idempotent.** `_annotate` strips exactly the columns the blocks supply before
re-joining (so a re-run replaces rather than appending `_x`/`_y`), then asserts **every pre-existing
column is bit-identical and the row count is unchanged** — aborting otherwise. Same contract as
`card_context._annotate`, which established the pattern; the owned-set derivation is stricter.

⚠ **RE-RUN THIS AFTER ANY `gene_card` / `edge_card` REBUILD** — those builders do not know about these
columns and will drop them.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
EDGE_CARD = OUT / "realization/edge_card.tsv"
GENE_CARD = OUT / "realization/gene_card.tsv"
SUBTYPE = OUT / "subtype_edge_heterogeneity.tsv"
# ⛔ NO MODULE-LEVEL PREFIX LIST. `_annotate` derives what it owns from the blocks it is given —
# see its docstring. A global list here is what destroyed the arm card's `chim_` block (MH-218).

# the SAME cuts the arm card uses — a ladder is only comparable across rungs if the rungs agree on it
COUPLING_CUTS = (-0.05, -0.10, -0.15, -0.20, -0.30)

#: ⛔ THE RETIRED-HEURISTIC BLOCK ON THE GENE CARD. All five arrive from
#: `analyses/misc/mirna_comovement.py` -> `tissue_reference/mirna_comovement/gene_corepression.tsv`, i.e.
#: the **§6b-RETIRED pressure heuristic**, NOT the learned model. `heur_n_regulators` in particular is a
#: different count from the model's `n_arms` (they differ on 306 of 1,409 genes, in both directions).
#: Prefixing them makes the provenance visible at a glance and groups them as ONE block.
HEUR_RENAME = {"n_regulators": "heur_n_regulators",
               "gene_repression_class": "heur_repression_class",
               "gene_net_repressed_tumor": "heur_net_repressed_tumor",
               "rho_gene_pressure_tumor": "heur_rho_pressure_tumor",
               "delta_tumor_nat": "heur_delta_tumor_nat"}


def heur_col(d, name: str):
    """Tolerant accessor for the retired-heuristic block — accepts either the old or the `heur_` name.
    Lets a consumer work against a card written before or after the 2026-08-19 rename."""
    new = HEUR_RENAME.get(name, name)
    for c in (new, name):
        if c in getattr(d, "columns", []):
            return d[c]
    return None


def gene_realization_ladder() -> pd.DataFrame:
    """Per GENE: how many of its regulators realize at each coupling depth (the arm ladder's mirror).

    ⚠ Denominator is `greal_n_scored` — the gene's regulators that COULD be scored, not its full design.
    ⚠ Rates are gated at n≥3 regulators. A "fraction of regulators realized" on a 1-regulator gene is
    0 or 1 and nothing else, which is axiom 5's coin-flip-with-decimals.
    ⚠ `coupling_tum` is bulk-observational: absence ≠ not-real (MH-166, ~13% sensitive).
    """
    d = pd.read_csv(EDGE_CARD, sep="\t", low_memory=False)
    if "coupling_tum" not in d.columns:
        return pd.DataFrame()
    d = d[d["coupling_tum"].notna()]
    g = d.groupby("gene")
    out = pd.DataFrame({"greal_n_scored": g.size(),
                        "greal_med_coupling": g["coupling_tum"].median(),
                        "greal_best_coupling": g["coupling_tum"].min()})
    for cut in COUPLING_CUTS:
        tag = f"c{abs(int(cut * 100)):02d}"
        out[f"greal_n_{tag}"] = (d[d.coupling_tum < cut].groupby("gene").size()
                                 .reindex(out.index).fillna(0))
    if "coupling_p_tum" in d.columns:
        out["greal_n_sig"] = (d[d.coupling_p_tum < 0.05].groupby("gene").size()
                              .reindex(out.index).fillna(0))
    ok = out["greal_n_scored"] >= 3
    out["greal_frac_c10"] = (out["greal_n_c10"] / out["greal_n_scored"]).where(ok)
    return out.reset_index()


def gene_arm_resolution() -> pd.DataFrame:
    """⭐ IS THIS GENE MODELABLE AT ARM RESOLUTION? — the gene-level roll-up of arm identifiability, which
    did not exist (column review, unit 1, user-asked 2026-08-19).

    `arm_rung.py` answers *can these same-seed arms be told apart* per (gene, seed_family) CELL, and
    `aid_*` rolls it up per ARM. Nobody rolled it up per GENE, so *"can I trust an arm-level answer for
    THIS gene"* had no column — even though it decides whether an arm-resolved claim about the gene is
    admissible at all.

    ⚠ THE DENOMINATOR IS THE POINT. Only genes with at least one MULTI-ARM cell can be asked: elsewhere the
    arm rung IS the family rung by construction and the question is undefined, not negative. Measured:
    **273 of 1,549 genes (17.6%)** have any multi-arm cell — so `armres_class` is `not_applicable` for 82%,
    and that is a design fact, not a failure.

    Among the 273: `armres_frac_resolvable` **median 0.00**, ALL cells resolvable in **68** genes and NONE in
    **168** ⇒ the distribution is strongly BIMODAL, so do not read its mean. `armres_med_drho` median
    −0.0012 with a best of −0.206 — the gain is a TAIL, concentrated in the miR-29 family.
    """
    p = OUT / "realization/edge_card.tsv"
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t", low_memory=False)
    need = {"gene", "seed_family", "n_arm_in_cell"}
    if not need <= set(d.columns):
        return pd.DataFrame()
    num = lambda s: pd.to_numeric(s.map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0})
                                  if s.dtype == object else s, errors="coerce")
    multi = d[num(d["n_arm_in_cell"]) > 1]
    if multi.empty:
        return pd.DataFrame()
    cells = multi.drop_duplicates(["gene", "seed_family"])
    out = pd.DataFrame({
        "armres_n_multi_cells": cells.groupby("gene").size(),
        "armres_n_arms_split": multi.groupby("gene")["arm"].nunique(),
        "armres_frac_resolvable": multi.groupby("gene")["arm_resolvable"].apply(lambda s: num(s).mean()),
        "armres_med_sep_z": multi.groupby("gene")["arm_sep_z"].apply(lambda s: num(s).median()),
        "armres_best_drho": multi.groupby("gene")["oof_drho"].apply(lambda s: num(s).min()),
        "armres_med_drho": multi.groupby("gene")["oof_drho"].apply(lambda s: num(s).median()),
    })
    fr = out["armres_frac_resolvable"]
    out["armres_class"] = np.where(fr >= 0.999, "arm_modelable",
                        np.where(fr <= 0.001, "family_only", "partial"))
    return out.reset_index()


# ⭐⭐ THE NaN-BOOLEAN REPAIR TABLE (MH-256, 2026-08-19). Each entry is (flag, inputs, mode).
# ⛔ THE CLASS: `NaN < x`, `NaN > x`, `NaN >= x` and `series.eq(v)` on NaN all return **False**, never NaN,
# so a flag built as a bare comparison of possibly-missing numbers is defined on EVERY row and each
# unmeasurable row silently joins the NEGATIVE class. Thirteen instances were found in one session; the
# eight below were found in one pass by `analyses/ops/nan_bool_audit.py` rather than one block at a time.
# ⚠ Fixed at SOURCE as well — this table exists because several sources (`card.py`) only re-run under the
# expensive canonical rebuild, while every input is already on the delivered card, so the repair is a pure
# derivation. `mode="any"` masks where ALL inputs are missing (an `|` of comparisons); `mode="group_max"`
# masks on the gene-level max, because that flag is a GENE property and a NaN row inside a measured gene is
# legitimate — testing the row's own value overstated it 10× (1,570 rows vs the true 154).
_NAN_FLAG_REPAIRS = {
    "edge": [("gene_dominated", ["share_TUM"], "group_max"),
             ("healthy_uninformative", ["healthy_potential"], "all"),
             ("arm_spiker", ["arm_pct_floor", "arm_iqr"], "all"),
             ("ctx_measurable", ["ctx_ceiling"], "all")],
    "gene": [("ctx_measurable", ["ctx_ceiling"], "all")],
    "family": [("ctx_measurable", ["ctx_ceiling"], "all")],
    "arm": [("famrole_is_dominant", ["famrole_abund_share"], "all"),
            ("dose_confounded", ["dose_comp_retention", "dose_prolif_retention"], "any"),
            ("hly_from_seedmate", ["hly_baseline_src"], "all"),
            ("arm_spiker", ["arm_pct_floor", "arm_iqr"], "all")],
}


def _repair_nan_flags(d: pd.DataFrame, card: str) -> int:
    """Mask the flags in `_NAN_FLAG_REPAIRS` wherever their inputs are missing. Returns rows unmasked."""
    moved = 0
    for flag, ins, mode in _NAN_FLAG_REPAIRS.get(card, []):
        if flag not in d.columns or not all(i in d.columns for i in ins):
            continue
        ok = None
        for i in ins:
            col = d[i]
            s = col.notna() if col.dtype == object else pd.to_numeric(col, errors="coerce").notna()
            ok = s if ok is None else (ok | s if mode == "any" else ok & s)
        if mode == "group_max" and "gene" in d.columns:
            v = pd.to_numeric(d[ins[0]], errors="coerce")
            ok = v.groupby(d["gene"]).transform("max").notna()
        was = d[flag].notna().sum()
        d[flag] = d[flag].where(ok)
        if d[flag].notna().sum() != was:
            moved += was - int(d[flag].notna().sum())
            print(f"  ✅ {card}.{flag}: {was} -> {int(d[flag].notna().sum())} defined "
                  f"({was - int(d[flag].notna().sum())} unmeasurable rows unmasked from a silent False)")
    return moved


# ⭐ RATIO GATES APPLIED IN PLACE (MH-257) — same reasoning as `_NAN_FLAG_REPAIRS`: fixed at source, but
# the source only re-runs under the expensive canonical rebuild while the denominator is already on the
# delivered card, so the repair is a pure derivation. ⚠ Only columns whose DENOMINATOR is on the card can
# be repaired here — `retention_rho` and `own_specific_frac` cannot, and stay ungated until the rebuild.
_RATIO_GATE_REPAIRS = {"gene": [("realized_retention", "realized_rho_raw")]}


def _repair_ratio_gates(d: pd.DataFrame, card: str) -> None:
    from mirna_hallmark.config import RHO_GATE
    for col, den in _RATIO_GATE_REPAIRS.get(card, []):
        if col not in d.columns or den not in d.columns:
            continue
        r = pd.to_numeric(d[den], errors="coerce")
        x = pd.to_numeric(d[col], errors="coerce")
        was, wmax = int(x.notna().sum()), float(x.abs().max())
        d[col] = x.where(r.abs() >= RHO_GATE)
        now = pd.to_numeric(d[col], errors="coerce")
        if int(now.notna().sum()) != was:
            print(f"  ✅ {card}.{col}: gated at |{den}| >= {RHO_GATE} — {was} -> {int(now.notna().sum())} "
                  f"rows, max|x| {wmax:.1f} -> {float(now.abs().max()):.2f}")


def normalise_flag_cards() -> None:
    """Apply `_NAN_FLAG_REPAIRS` to the EDGE and ARM cards, which have no name-normaliser of their own.

    ⚠ Kept separate from `normalise_gene_card` / `normalise_family_card` because those also rename and
    prune; this one ONLY masks flags, so it is safe to run on any card at any time.
    """
    for card, rel in (("edge", "realization/edge_card.tsv"), ("arm", "arm_card.tsv")):
        path = OUT / rel
        if not path.exists():
            continue
        d = pd.read_csv(path, sep="\t", low_memory=False)
        moved = _repair_nan_flags(d, card)
        # ⛔ unit 21: drop the pruned constant from the delivered card too
        if card == "edge" and "echim_any" in d.columns:
            d = d.drop(columns=["echim_any"]); moved += 1
            print("  ✅ edge.echim_any DROPPED (constant True, 0 bits beyond echim_n_sources.notna())")
        if moved:
            d.to_csv(path, sep="\t", index=False)


def normalise_gene_card() -> None:
    """⭐ THE GENE CARD's PUBLIC NAMES — applied in place, no refit (column review unit 3, 2026-08-19).

    The edge card got `realization._normalise_edge_names` as the ONE place its delivered vocabulary is
    decided. The gene card had no such funnel, so its two schema debts were wrongly parked behind the
    Gibbs refit queue — both are **pure derivations** over columns already on disk.

      `n_discovered` -> `n_pip_disc_gt50`   ⛔ the old name reads as a RESULT. It is `pip_discovery > 0.5`,
          a posterior-inclusion count, and it is >0 for **725 genes (46.8%)** — while per-edge and
          per-family discovery are EMPTY under the honest empirical FDR. Reading it as "genes with
          discoveries" inverts the axis. Renamed at source in `readouts.py` too, so the two agree.
      `n_dense_included`  DROPPED           ⛔ bit-identical to `n_fam` on every row (verified). It once
          produced duplicate scan results that looked like two independent axes agreeing.
    """
    if not GENE_CARD.exists():
        return
    g = pd.read_csv(GENE_CARD, sep="\t", low_memory=False)
    before = g.shape[1]
    if "n_discovered" in g.columns and "n_pip_disc_gt50" not in g.columns:
        g = g.rename(columns={"n_discovered": "n_pip_disc_gt50"})
    if "n_dense_included" in g.columns:
        g = g.drop(columns=["n_dense_included"])
    # ⛔ `n_regulators` -> `heur_n_regulators` (user-directed 2026-08-19). It is NOT the model's count:
    # it is `edges.groupby('gene')['miRNA'].nunique()` from `mirna_comovement`, i.e. the §6b-RETIRED
    # heuristic pressure lane's edge table. `n_arms` is the fit's own width, and the two differ on 306 of
    # 1,409 genes IN BOTH DIRECTIONS. Two counts reading as one name is axiom 6's collision class.
    # ⭐ AND ITS FOUR SIBLINGS TOO (user-directed 2026-08-19) — they arrive from the SAME retired lane, so
    # the whole block is renamed together and reads as one provenance unit rather than five loose columns.
    # ⚠ RENAMED ON THE CARD ONLY. `gene_corepression.tsv` and its direct readers
    # (`pressure_attribution_validation`, `decoupling_validation`) keep the original names — the rename is
    # about what the CARD delivers, not about rewriting a retired lane's own outputs.
    for _old, _new in HEUR_RENAME.items():
        if _old in g.columns and _new not in g.columns:
            g = g.rename(columns={_old: _new})
    # ⛔⛔ REPAIR THE CPTAC "beats abundance" FLAGS IN PLACE (column review unit 17, 2026-08-19).
    # `cptac_card` built them with a bare `rp < ab`, and `NaN < NaN` is **False**, not NaN — so a gene whose
    # protein coupling was never computed read as "beta does NOT beat abundance". Fixed at source, but the
    # source is an expensive CPTAC re-scoring, and the repair is a PURE DERIVATION over columns already on
    # disk: mask wherever either input is missing. ✅ The comparison itself was correct on measurable rows
    # (verified 100% match), so only the mask changes.
    _repair_nan_flags(g, "gene")
    _repair_ratio_gates(g, "gene")
    for _pref in ("cptac_prosp", "cptac_t105"):
        _f, _a, _g = f"{_pref}_agg_beats_abund_prot", f"{_pref}_abund_rho_prot", f"{_pref}_agg_rho_prot"
        if _f in g.columns and _a in g.columns and _g in g.columns:
            _ok = pd.to_numeric(g[_a], errors="coerce").notna() & pd.to_numeric(g[_g], errors="coerce").notna()
            _was = g[_f].notna().sum()
            g[_f] = g[_f].where(_ok)
            if g[_f].notna().sum() != _was:
                print(f"  ✅ {_f}: {_was} -> {g[_f].notna().sum()} defined "
                      f"({_was - g[_f].notna().sum()} unmeasurable rows unmasked from a silent False)")
    if g.shape[1] != before or "n_pip_disc_gt50" in g.columns:
        g.to_csv(GENE_CARD, sep="\t", index=False)
        print(f"  ✅ gene_card normalised: {before} -> {g.shape[1]} cols "
              f"(n_discovered->n_pip_disc_gt50, n_dense_included dropped)")


def normalise_family_card() -> None:
    """⭐ THE GENE_FAMILY CARD's PUBLIC NAMES — and a rename I reported as done that reached ONE of two cards.

    ⛔ **THE DEFECT (found in the bare-column review, 2026-08-19).** Column-review unit C pruned `p_fam` and
    renamed `n` -> `n_samples`. That landed via `realization._normalise_edge_names`, which funnels the EDGE
    card only — so the gene_family card still carried BOTH `n` and `p_fam` while I had recorded the change
    as applied. This is doctrine §4.5 biting exactly as written: **`gene_family_card.tsv` has no
    `_finish_card` call site**, so it inherits none of the edge card's schema discipline. It now has a
    funnel of its own, which is the actual fix.

      `n` -> `n_samples`   the sample count the cell was fit on — CONSTANT 1,040 across the card. Bare `n`
          is the single most overloaded name in the codebase; `n_samples` says which n.
      `p_fam`  DROPPED     ⛔ it is NOT a p-value. It is the DESIGN DIMENSION (predictor count), and on this
          card it is bit-identical to `n_fam`. A name that reads as a p-value on a card full of real
          p-values is axiom 6's collision class, and here it was also redundant.
    """
    p = OUT / "gene_family_card.tsv"
    if not p.exists():
        return
    d = pd.read_csv(p, sep="\t", low_memory=False)
    before = d.shape[1]
    _repair_nan_flags(d, "family")
    if "n" in d.columns and "n_samples" not in d.columns:
        d = d.rename(columns={"n": "n_samples"})
    if "p_fam" in d.columns:
        d = d.drop(columns=["p_fam"])
    if d.shape[1] != before or "n_samples" in d.columns:
        d.to_csv(p, sep="\t", index=False)
        print(f"  ✅ gene_family_card normalised: {before} -> {d.shape[1]} cols (n->n_samples, p_fam dropped)")


def discovery_queue_rollup(level: str) -> pd.DataFrame:
    """⭐ MAKE THE 157-EDGE DISCOVERY QUEUE VISIBLE FROM THE CARDS (user-asked 2026-08-19).

    ⛔ **WHY IT CANNOT BE AN EDGE FLAG.** The gold set and the curated edge card are **DISJOINT — 0 of the
    157 pairs appear on the edge card** (verified). They are discovery CANDIDATES on pairs curation never
    admitted, so there is no edge row to flag; inventing one would put un-curated pairs into the design.

    ⇒ surfaced as a COUNT at the two rungs that do overlap, so a reader of any gene or arm can see that the
    discovery lane touches it:
        `disc_n_gold_edges`   how many of the 157 name this gene / this arm
        `disc_gold_families`  which seed families they belong to (arm rung: usually one)

    ⚠⚠ **COVERAGE IS ASYMMETRIC AND THE GENE SIDE IS BADLY INCOMPLETE.** All **18 of 18** gold arms are on
    the arm card, but only **39 of the 90** gold genes are on the gene card — the rest have **no curated
    regulator at all**, which is exactly what makes them discovery candidates (MH-198: 94/157 sit on genes
    with no curated regulator). ⇒ **a 0 on the gene card means "no gold edge OR this gene's gold edges are
    invisible here"; it does NOT mean the queue is empty for it.** Read the queue from
    `discovery_gold_set.tsv`, and treat these columns as a pointer, never as the denominator.

    ⚠ The set is **11 seed families and 61% is ONE family** — always quote the families beside the edges.
    ⚠ It is a ranked VALIDATION QUEUE, not 157 findings: per-edge FDR is empty by construction.
    """
    p = OUT / "discovery_gold_set.tsv"
    if not p.exists():
        return pd.DataFrame()
    gs = pd.read_csv(p, sep="\t")
    key = "gene" if level == "gene" else "arm"
    if key not in gs.columns:
        return pd.DataFrame()
    out = pd.DataFrame({"disc_n_gold_edges": gs.groupby(key).size()})
    if "seed_family" in gs.columns:
        out["disc_gold_families"] = gs.groupby(key)["seed_family"].apply(
            lambda s: ";".join(sorted(set(s))))
    out = out.reset_index()
    # ⛔ FIXED 2026-08-19 (column review unit 8): the count was emitted as NaN for rows the gold set does
    # not name, which is WRONG — **every** row on both cards WAS checked against the set; the join simply
    # found nothing. A NaN there is indistinguishable from "not scanned", which is the exact ambiguity the
    # `cov_` flags exist to remove, so this block was manufacturing it. Emit a MEASURED ZERO instead.
    # ⚠ `disc_gold_families` stays NaN — there is genuinely no family to name, and "" would read as a value.
    card = GENE_CARD if key == "gene" else OUT / "arm_card.tsv"
    if card.exists():
        universe = pd.read_csv(card, sep="\t", usecols=[key], low_memory=False)[key].dropna().unique()
        out = (pd.DataFrame({key: universe}).merge(out, on=key, how="left"))
        out["disc_n_gold_edges"] = out["disc_n_gold_edges"].fillna(0).astype(int)
    return out


def edge_hly_leg_concordance() -> pd.DataFrame:
    """⭐ DO THE TWO HLY→TUM LEGS AGREE? — derived, no refit (column review unit 5, 2026-08-19).

    ⛔ **THE PROBLEM.** `arm_lfc_HLY_TUM_QN` and `arm_lfc_HLY_TUM_raw` are nominally the same comparison,
    and they are not: spearman **+0.544**, median |diff| **1.49**, and they **DISAGREE ON SIGN FOR 23.9%
    of edges**. Which leg a reader picks flips the direction for ~1 edge in 4 — silently, because both
    columns look equally authoritative.

    ⚠ Compounding it, **only `arm_lfc_NAT_TUM` is same-platform (TCGA→TCGA); every HLY leg crosses
    GTEx→TCGA**, the boundary the `fst_` work measured as producing 5–8× more near-total flips.

    ⇒ `hly_leg_concordant` marks the edges where the two legs at least AGREE ON SIGN, so a cross-state
    claim can be gated on it instead of resting on an arbitrary choice of leg. **It is not a correctness
    flag** — agreement does not make either leg right — it is a *disagreement* flag, which is the honest
    thing the data supports. NaN where either leg is missing.
    """
    if not EDGE_CARD.exists():
        return pd.DataFrame()
    d = pd.read_csv(EDGE_CARD, sep="\t", low_memory=False)
    need = {"gene", "arm", "arm_lfc_HLY_TUM_QN", "arm_lfc_HLY_TUM_raw"}
    if not need <= set(d.columns):
        return pd.DataFrame()
    q = pd.to_numeric(d["arm_lfc_HLY_TUM_QN"], errors="coerce")
    r = pd.to_numeric(d["arm_lfc_HLY_TUM_raw"], errors="coerce")
    ok = q.notna() & r.notna()
    conc = (np.sign(q) == np.sign(r)).where(ok)
    return pd.DataFrame({"gene": d["gene"], "arm": d["arm"], "hly_leg_concordant": conc})


def gene_concentration_adj() -> pd.DataFrame:
    """⭐ THE FLOOR-CORRECTED CONCENTRATION — derived, no refit (column review unit 3, 2026-08-19).

    ⛔ **`concentration` = `beta.max() / beta.sum()` is a TOP SHARE, and a top share is BOUNDED BELOW BY
    1/k** where k = the gene's family count. So it falls as designs widen *whatever the biology does* —
    this is the moving-support trap (axiom 8), and here it is at its most extreme:

        spearman(concentration, n_fam)                  = −0.9399   almost entirely mechanical
        spearman(floor-corrected, n_fam)                = −0.2028   the honest number

    The floor is not theoretical — the data sits ON it: at `n_fam == 1` concentration is **exactly 1.0000
    for all 730 genes (47%)**, and at `n_fam == 2` the minimum is **0.5005** against a floor of 0.5000.

    ⇒ `concentration_adj = (c − 1/k) / (1 − 1/k)` ∈ [0,1] — "how concentrated is this gene RELATIVE to the
    most diffuse its design allows". **Undefined (NaN) at k = 1**, which is correct and important: with one
    family there is nothing to concentrate, so the question does not exist for 47% of genes. That NaN is
    the honest answer, not a gap — `gene_axes.mask_degenerate` exists for exactly this.

    ⚠ Ships BESIDE the raw column, never replacing it (axiom 2a: flag, don't delete). Report the raw value
    only with `n_fam` next to it; report `concentration_adj` when comparing genes of different widths.
    """
    if not GENE_CARD.exists():
        return pd.DataFrame()
    g = pd.read_csv(GENE_CARD, sep="\t", low_memory=False)
    if not {"gene", "concentration", "n_fam"} <= set(g.columns):
        return pd.DataFrame()
    c = pd.to_numeric(g["concentration"], errors="coerce")
    k = pd.to_numeric(g["n_fam"], errors="coerce")
    adj = ((c - 1.0 / k) / (1.0 - 1.0 / k)).where(k > 1)
    return pd.DataFrame({"gene": g["gene"], "concentration_adj": adj.clip(0, 1).round(4)})


def gene_top_identity_gated() -> pd.DataFrame:
    """⭐ THE GENE CARD's `top_identity`, GATED — the last ungated survivor (column review unit 11).

    ⛔ `readouts` computes `top_identity = float(d.identity.max())` per gene, ungated. `identity` is a
    SIGNED share (≈10% negative), so the max is unbounded: on the gene card it still reaches **+740.007**,
    and 80 genes exceed 1. The `identity_reliable` gate shipped to the edge and gene_family cards in the
    unit-1/2 review but **never reached this DERIVED column** — a gate on the inputs does not gate a
    statistic someone already computed from them.

    ⇒ recompute the max over `identity_reliable` edges only. Pure derivation from the edge card, so no
    refit. Emits `top_identity_gated` (bounded) and `top_identity_n_reliable` (its denominator — a max over
    2 families is not the same statement as a max over 12, and axiom 5 says print the denominator).
    ⚠ Ships BESIDE the raw column rather than replacing it: `readouts` still owns `top_identity`, and
    silently changing a value another module wrote is how provenance rots.
    """
    if not EDGE_CARD.exists():
        return pd.DataFrame()
    d = pd.read_csv(EDGE_CARD, sep="\t", low_memory=False)
    if not {"gene", "identity", "identity_reliable"} <= set(d.columns):
        return pd.DataFrame()
    rel = d[d["identity_reliable"].map({True: True, "True": True}).fillna(False)]
    if rel.empty:
        return pd.DataFrame()
    i = pd.to_numeric(rel["identity"], errors="coerce")
    g = rel.assign(_i=i).groupby("gene")["_i"]
    return pd.DataFrame({"top_identity_gated": g.max().round(4),
                         "top_identity_n_reliable": g.size()}).reset_index()


def _identity_gate(d: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """⭐ THE IDENTITY GATE, DERIVED — no refit required (user-asked 2026-08-19).

    `readouts.add_reliability` gained `identity_coherence` / `identity_abs` / `identity_reliable` in the
    unit-1 review, but **that code has never reached a card**: both the edge and gene_family cards are
    written upstream of it and neither has been regenerated since, so `identity` ships UNGATED on both
    (verified 2026-08-19 — the columns are absent from each). Reporting the gate as shipped was wrong.

    ⭐ **AND NO REFIT IS NEEDED TO FIX IT.** The gate is a PURE FUNCTION of the `identity` column already
    on the card plus its per-gene sum — it re-estimates nothing. A Gibbs refit is required only when the
    ESTIMATES change; here they do not. So this derives the three columns in place, and `annotate()` will
    keep deriving them after every rebuild.

        identity_coherence = 1 / sum|identity| per gene   in (0,1]; 1 when nothing cancels, -> 0 as
                                                          suppressor cancellation grows (identity sums
                                                          to exactly 1 per gene, so |sum| is 1)
        identity_abs       = |identity| / sum|identity|   the BOUNDED companion, always in [0,1]
        identity_reliable  = coherence >= 0.5 AND |identity| <= 1

    ⚠ The grouping key is the GENE on both cards, because identity is normalised per gene — NOT per
    (gene, family) and NOT per arm. Grouping on the card's own key would give every row coherence 1.0.
    """
    if "identity" not in d.columns or "gene" not in d.columns:
        return pd.DataFrame()
    i = pd.to_numeric(d["identity"], errors="coerce")
    tot = i.abs().groupby(d["gene"]).transform("sum")
    with np.errstate(divide="ignore", invalid="ignore"):
        coh = (1.0 / tot.replace(0, np.nan)).clip(0, 1)
        iab = (i.abs() / tot.replace(0, np.nan)).clip(0, 1)
    out = d[keys].copy()
    out["identity_coherence"] = coh.round(4)
    out["identity_abs"] = iab.round(4)
    # ⛔⛔ MASKED 2026-08-19 (column review unit 19): a bare `&` of two comparisons cannot say "unknown"
    # — `NaN >= 0.5` and `NaN <= 1.0` are both False, so an edge whose `identity` was NEVER COMPUTED read
    # as "unreliable". Measured: **218 of the 377 False rows (57.8%) had NaN identity** ⇒ anyone filtering
    # `identity_reliable == False` to study unreliable attribution got a set that was 58% "not computed",
    # and the unreliable RATE was overstated 2.3× (6.7% -> 2.9% on the measurable set). Fifth instance of
    # this class in one session. **A three-state quantity needs three states: True / False / NaN.**
    out["identity_reliable"] = ((coh >= 0.5) & (i.abs() <= 1.0)).where(i.notna() & coh.notna())
    return out


def edge_identity_gate() -> pd.DataFrame:
    """The identity gate for the EDGE card. See `_identity_gate`."""
    if not EDGE_CARD.exists():
        return pd.DataFrame()
    return _identity_gate(pd.read_csv(EDGE_CARD, sep="\t", low_memory=False), ["gene", "arm"])


def family_identity_gate() -> pd.DataFrame:
    """The identity gate for the GENE_FAMILY card. See `_identity_gate`."""
    p = OUT / "gene_family_card.tsv"
    if not p.exists():
        return pd.DataFrame()
    return _identity_gate(pd.read_csv(p, sep="\t", low_memory=False), ["gene", "family"])


def _paired_dshare(d: pd.DataFrame, lin: dict, a: str, b: str) -> pd.Series:
    """⛔⛔ THE SHARE DELTAS HAD A MOVING DENOMINATOR — measured, and it was manufacturing the finding.

    `fst_share_{HLY,NAT,TUM}` come from THREE DIFFERENT arm-card abundance columns with very different
    fill: `arm_med_rpm` **19.8%**, `hly_nat_median` **74.0%**, `hly_gtex_median` **17.0%**. Each state's
    share therefore normalises over a DIFFERENT member set, so a naive `share_TUM - share_NAT` differences
    two fractions with **different denominators**.

    Measured on the 320 multi-member families (2026-08-19): the TUM denominator holds a mean of **0.66**
    members against NAT's **1.88** (2.8×), only **21.2%** of families use the same member count, and the
    raw delta came out **mean +0.156 / median +0.032, 67.6% positive, Wilcoxon p=1.2e-12** — i.e. it read
    as *"arms systematically gain family share in tumour"* when the sparser state simply has fewer members
    to divide by. **A statistic whose denominator moves with coverage is measuring the coverage** (axiom 5).

    ⇒ both shares are RE-NORMALISED over the members measured in BOTH states before differencing. The
    delta is NaN where fewer than 2 such members exist, because a one-member common set forces both shares
    to 1.0 and the difference to a structural 0 (the degeneracy trap, axiom 8) — a fake zero is worse than
    a blank. `fst_n_common_{A}_{B}` ships alongside so the denominator is visible.
    """
    if a not in lin or b not in lin:
        return pd.Series(np.nan, index=d.index)
    la, lb = lin[a], lin[b]
    both = la.notna() & lb.notna()
    fam = d["seed_family"]
    n_common = both.groupby(fam).transform("sum")
    ok = both & (n_common >= 2)
    out = pd.Series(np.nan, index=d.index)
    for src, sign in ((lb, 1.0), (la, -1.0)):
        masked = src.where(ok)
        tot = masked.groupby(fam).transform("sum")
        share = masked / tot.replace(0, np.nan)
        out = out.add(sign * share, fill_value=0.0) if sign > 0 else out.sub(-sign * share)
    d[f"fst_n_common_{a}_{b}"] = n_common.where(both)
    return out.where(ok).round(4)


def family_state_shift() -> pd.DataFrame:
    """⭐ THE ARM'S SHARE OF ITS OWN FAMILY, ACROSS STATES — gene-free (user-asked 2026-08-19).

    ⛔ THE GAP THIS FILLS. Every cross-state share/rank in the system is anchored to a GENE
    (`share_HLY/NAT/TUM`, `rank_*`, `d_rank_*`) or to the whole cohort (`grank_*`, `dGlobal_*`). The arm's
    standing INSIDE ITS OWN SEED FAMILY exists only as a single-state snapshot (`famrole_abund_share`,
    `famrole_is_dominant`) — verified 2026-08-19: no HLY/NAT/TUM variant of any `famrole_` column exists.
    So *"does the miR-29 family's dominant member change between healthy and tumour"* had no answer.

    ⭐ WHY IT IS WORTH HAVING, beyond symmetry. The gene-level `regulatory_handoff` is **mechanically
    monotone in design width** (1.6% at n_fam==1 vs 38.6% at n_fam>=2), so every use of it needs
    conditioning. A family-internal switch has NO such dependence — the denominator is the family itself,
    which is fixed regardless of how many genes it regulates. It is the confound-free version of the same
    question, and it speaks directly to MH-215's structural fact that one member carries ~62% of a typical
    multi-member family's dose.

    Emits per ARM (arm rung, gene-free): `fst_share_{HLY,NAT,TUM}` — the arm's fraction of its family's
    total linear dose in each state; `fst_d_share_HLY_TUM` / `_NAT_TUM`; `fst_rank_{HLY,TUM}` within the
    family; `fst_is_dominant_{HLY,TUM}` (>0.5 of family dose) and `fst_dominance_switch`.

    ⚠ SINGLETON FAMILIES ARE DEGENERATE — share ≡ 1.0 in every state by construction, so the switch flag is
    structurally False. **66.9% of arms (1,639 of 2,450)** are singletons; pooling them is the degeneracy
    trap (`gene_axes.mask_degenerate`). Mask on `fst_n_members > 1`.

    ⛔⛔ **USE THE NAT LEG. THE HEALTHY LEG IS CROSS-PLATFORM AND CARRIES AN ARTIFACT — measured, and the
    obvious guard does NOT fix it.** `_HLY` comes from GTEx while `_TUM` comes from TCGA, so a share
    computed across them compares two quantification pipelines, on top of MH-210's multi-mapping collapse.
    I predicted that requiring EVERY family member to have a measured GTEx level would clean it up. **It did
    not — it made it worse**, which is the useful result:

        guard                      n     dominance switch   near-total flips (|Δshare|>0.90)
        this arm measured        177          46.3%                 7.1%
        ALL members measured      68          52.9%                11.6%      <- stricter, WORSE
        ── same platform ──
        NAT -> TUM               210          26.7%                 1.4%      <- 5-8x fewer flips

    A near-total flip (share 0.0002 -> 1.0000) is not biology; it is one pipeline seeing an arm the other
    cannot. It is **5–8× rarer** the moment both legs come from TCGA. ⇒ the missing-sibling story was wrong;
    the defect is the platform boundary itself, and no membership guard can repair it.

    ⛔⛔ **AND THE RAW 26.7% SWITCH RATE IS MOSTLY MINOR ARMS — GATE IT (user-caught, measured).** The
    unfiltered NAT→TUM dominance switch is 56/210 = 26.7%, but it is an axiom-5 threshold count sitting
    exactly where the mass piles up:
      * switching arms have median **log2 RPM 0.80 vs 3.30** for stable ones, and **73% sit BELOW the
        expression floor** (log2 11 ≈ 3.46) vs 51%;
      * they sit a median **0.133 from the 0.5 boundary**, with **29% within 0.05** — coin flips;
      * the rate is monotone in abundance, arm tertile **38.0% → 26.8% → 14.7%**, and monotone in the
        FAMILY's own abundance, **47.1% → 20.0% → 12.9%**.
    ✅ **GATED — both arms above the floor AND the switch margin > 0.10 — the rate is 10.7% (8 of 75).**
    Quote **~11%**, never 26.7%. The eight survivors are legible and include a clean RECIPROCAL pair
    (miR-196b-5p 0.75→0.32 against miR-196a-5p 0.25→0.68, both well above floor), plus
    miR-30e-3p 0.29→1.00 (RPM 12.4), miR-200a-3p 0.28→1.00, miR-96-5p, miR-19a-3p, miR-190b, miR-193a-3p.

    ✅ **WHAT IS DEFENSIBLE:** the same-platform malignant step, gated as above. Within-family dominance
    moves for ~11% of well-measured multi-member arms — and unlike the gene-level `regulatory_handoff` it
    cannot be a design-width artifact, because the denominator is the family.
    ⚠ Still DESCRIPTIVE: no decoy, no null. A registry CANDIDATE.
    """
    try:
        from mirna_hallmark.learned import data as LD, families as FAM
        from mirna_hallmark import healthy_anchor as HA          # noqa: F401  (presence check)
    except Exception as ex:
        print(f"  ⚠ family_state_shift unavailable ({ex})")
        return pd.DataFrame()
    a_path = OUT / "arm_card.tsv"
    if not a_path.exists():
        return pd.DataFrame()
    A = pd.read_csv(a_path, sep="\t", low_memory=False)
    num = lambda s: pd.to_numeric(s, errors="coerce")
    lev = {"TUM": "arm_med_rpm", "NAT": "hly_nat_median", "HLY": "hly_gtex_median"}
    have = {k: v for k, v in lev.items() if v in A.columns}
    if "TUM" not in have or "seed_family" not in A.columns:
        return pd.DataFrame()
    d = A[["arm", "seed_family"] + list(have.values())].copy()
    d = d[d.seed_family.notna()]
    lin_by_state: dict[str, pd.Series] = {}
    for st, col in have.items():
        lin = np.power(2.0, num(d[col])) - 1.0        # log2 RPM -> linear, matching the dose convention
        lin = lin.clip(lower=0)
        tot = lin.groupby(d.seed_family).transform("sum")
        d[f"fst_share_{st}"] = (lin / tot.replace(0, np.nan)).round(4)
        d[f"fst_rank_{st}"] = lin.groupby(d.seed_family).rank(ascending=False, method="min")
        lin_by_state[st] = lin
    d["fst_n_members"] = d.groupby("seed_family")["arm"].transform("size")
    # ⭐ THE BLOCK'S OWN COVERAGE FLAG (column review unit 8). 485 of 2,450 arms have a family share; the
    # rest were NEVER COMPUTED because the abundance input is absent — that is genuinely UNSCANNED, not a
    # zero share, and without a flag the two are indistinguishable. ⚠ Only 16 of the arm card's 42 blocks
    # carry such a flag; this is one of the two added rather than inherited.
    d["cov_fst"] = d["fst_share_TUM"].notna()
    if "HLY" in have:
        d["fst_hly_measured"] = num(d[have["HLY"]]).notna()
        # the stricter guard, emitted because it is INFORMATIVE even though it does not repair the leg:
        # requiring every sibling to be measured RAISES the flip rate (7.1% -> 11.6%), which is what
        # localised the defect to the platform boundary rather than to missing members.
        d["fst_hly_family_measured"] = d.groupby(d["seed_family"])["fst_hly_measured"].transform("all")
        d["fst_d_share_HLY_TUM"] = _paired_dshare(d, lin_by_state, "HLY", "TUM")
        # ⛔ FIXED 2026-08-19: `NaN > 0.5` is **False**, not NaN — so an arm with NO measured share was
        # reading as "measured, and not dominant" on 100% of rows. That is the `echim_any` failure again
        # (a flag that can never say "unknown"), and it was in a block I had just written. Mask explicitly.
        d["fst_is_dominant_HLY"] = (d["fst_share_HLY"] > 0.5).where(d["fst_share_HLY"].notna())
        d["fst_is_dominant_TUM"] = (d["fst_share_TUM"] > 0.5).where(d["fst_share_TUM"].notna())
        # ⚠ a switch is only MEANINGFUL where the family has >1 member AND the healthy leg is measured
        ok = (d["fst_n_members"] > 1) & d["fst_hly_measured"]
        d["fst_dominance_switch"] = np.where(
            ok, d["fst_is_dominant_HLY"].ne(d["fst_is_dominant_TUM"]), np.nan)
    if "NAT" in have:
        d["fst_d_share_NAT_TUM"] = _paired_dshare(d, lin_by_state, "NAT", "TUM")
    return d.drop(columns=list(have.values()) + ["seed_family"])


def edge_subtype() -> pd.DataFrame:
    """Per (gene, arm): which PAM50 subtype the edge couples in.

    ⚠ MH-165: subtype coupling heterogeneity is REAL but **DISTRIBUTIONAL** — 16% of edges vs a
    marginal-preserving null's 5%, net of composition-slope modification. It is NOT a validated per-edge
    subtype assignment, and the "coupling surfaces in a subtype" discovery hope tested AT CHANCE
    (337 vs 287±35, z=+1.4). Read `esub_best_sub` as descriptive, never as a claim about that edge.
    """
    if not SUBTYPE.exists():
        return pd.DataFrame()
    d = pd.read_csv(SUBTYPE, sep="\t")
    keep = [c for c in ("gene", "arm", "rho_pool", "rho_LumA", "rho_LumB", "rho_Basal", "rho_Her2",
                        "best_sub", "rho_best", "q_adj") if c in d.columns]
    # ⛔ prefix is `esub_`, NOT `sub_`: the ARM card already owns `sub_` for its per-arm AGGREGATE
    # of this same source, and `card_rungs.domain_of` matches by prefix GLOBALLY across cards —
    # sharing it would file the edge block under the arm block's domain string.
    d = d[keep].rename(columns={c: f"esub_{c}" for c in keep if c not in ("gene", "arm")})
    return d.drop_duplicates(subset=["gene", "arm"])


def edge_ago_measured() -> pd.DataFrame:
    """Per (gene, arm): was `ago_loading` MEASURED, or is it the `.fillna(1.0)` default? (MH-214)"""
    if not EDGE_CARD.exists():
        return pd.DataFrame()
    from mirna_hallmark.learned import ago_loading as AG
    d = pd.read_csv(EDGE_CARD, sep="\t", usecols=["gene", "arm"], low_memory=False)
    arms = sorted(d["arm"].dropna().astype(str).unique())
    m = AG.measured(arms)
    d["ago_loading_measured"] = d["arm"].astype(str).map(m)
    return d.drop_duplicates(subset=["gene", "arm"])


def edge_admissibility() -> pd.DataFrame:
    """⭐⭐ THE MOST IMPORTANT MISSING COLUMN: does this edge even have a SEED SITE? (MH-216)

    The edge card carried β, calibrated coupling, realization, subtype and CN — and **no way to tell
    whether the arm has a seed match in the target's 3′UTR at all**. That is the project's single most
    load-bearing conditioning variable:

      ⭐ **MH-136, the decoy arc's first real positive control and its strongest result:** using the
      genome-wide 6mer map, **187 genes whose curated edges are ALL SEEDLESS give gap +0.0006 — exactly
      zero** — against **−0.0325** for 810 all-seeded genes. *No site → no repression → no coupling → no
      gap.* A priori, sequence-only, non-circular, correct on BOTH sides. It also closed MH-132's 21%
      positive-coupling mystery (`frac_seedless` predicts positive coupling, +0.042, p=2.2e-04).

    Without `adm_has_site` on the card, a reader can look at a strong coupling on a **seedless** edge and
    have no way to know. Source `admissibility_edge_tags.tsv` (MH-161): TargetScan context++ ∪ scanMiR
    duplex ∪ Poisson-gated 6mer for the site leg, TCGA-tumour ∨ GTEx-breast for the expression leg.

    ⚠ **`adm_admissible` is a DIAGNOSTIC FLAG, NOT A FILTER (MH-161's own ruling).** The pre-registered
    prediction that gating on it sharpens the decoy gap was **REFUTED** — n.s. at every width, pooled
    p=0.39 — and the gate shrinks `n_fam` 2→1, so filtering on it silently changes the design width the
    gap scales with. Census: has_site 89.9%, expressed 63.9%, admissible 60.4%.
    """
    p = OUT / "admissibility_edge_tags.tsv"
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    keep = [c for c in ("has_site", "expr_tcga", "expr_gtex", "expressed", "admissible")
            if c in d.columns]
    return (d[["gene", "arm"] + keep]
            .rename(columns={c: f"adm_{c}" for c in keep})
            .drop_duplicates(subset=["gene", "arm"]))


def gene_lit_ground_truth() -> pd.DataFrame:
    """⭐ WHO THE LITERATURE SAYS OWNS THIS GENE — the missing counterpart to the model's owner columns.

    The gene card carries FOUR "who owns this gene" columns — `top_family_identity`,
    `top_family_magnitude`, `realization_owner`, `static_owner_family` — and **nothing to compare them
    against**. `lit_ground_truth.tsv` (MH-196) is that comparator: canonical family = argmax of DISTINCT
    low-throughput functional PMIDs (reporter/western/proteomics), admitted only where the gene is in the
    design, the family is in that gene's design, ≥2 families compete, and the top is unambiguous.
    **329 genes / 92 families, sha256-stamped** — it replaced five mutually inconsistent, producer-less
    hand-curated lists.

    ⛔⛔ **`lit_agrees_*` IS DESCRIPTIVE AND IS NOT THE ATTRIBUTION RESULT. DO NOT READ IT AS ONE.**
    MH-196's verdict is a **cluster-bootstrapped RANK test with a family-FAME null**, and it is:
    *the model RANKS the canonical family above chance (β 0.436, p=0.021; identity 0.416, p=0.004) but
    **ARGMAX IS AT CHANCE AT EVERY k***. An agreement RATE is exactly that argmax statistic, so it will
    look unimpressive by construction and says nothing the rank test has not already said properly.
    ⚠ And the head-to-head β-vs-abundance is **unresolvable in principle** — only ~330 canonical families
    exist in the entire curated literature; separating a 0.033 effect needs ~1,241 clusters. These
    columns exist to make a per-gene case INSPECTABLE, never to be aggregated into a score.
    """
    p = OUT / "lit_ground_truth.tsv"
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    keep = [c for c in ("lit_family", "n_pmid", "margin", "n_lit_families", "tier") if c in d.columns]
    out = d[["gene"] + keep].rename(columns={c: (c if c.startswith("lit_") else f"lit_{c}")
                                             for c in keep})
    g = pd.read_csv(GENE_CARD, sep="\t", low_memory=False,
                    usecols=lambda c: c in ("gene", "top_family_identity", "top_family_magnitude"))
    out = out.merge(g, on="gene", how="left")
    for src, dst in (("top_family_identity", "lit_agrees_identity"),
                     ("top_family_magnitude", "lit_agrees_magnitude")):
        if src in out.columns:
            out[dst] = (out[src] == out["lit_family"]).where(out[src].notna())
    return out.drop(columns=[c for c in ("top_family_identity", "top_family_magnitude")
                             if c in out.columns]).drop_duplicates(subset=["gene"])


def arm_admissibility_rollup() -> pd.DataFrame:
    """Per-arm rollup of the site/expression legs — how much of this arm's design is even site-backed."""
    p = OUT / "admissibility_edge_tags.tsv"
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    g = d.groupby("arm")
    out = pd.DataFrame({"adm_n_edges": g.size()})
    for c, nm in (("has_site", "adm_n_with_site"), ("admissible", "adm_n_admissible")):
        if c in d.columns:
            out[nm] = d[d[c].astype(bool)].groupby("arm").size().reindex(out.index).fillna(0)
    # ⚠ gated: a fraction over 1-2 edges is a coin-flip with decimals (axiom 5)
    ok = out["adm_n_edges"] >= 3
    if "adm_n_with_site" in out.columns:
        out["adm_frac_with_site"] = (out["adm_n_with_site"] / out["adm_n_edges"]).where(ok)
    return out.reset_index()


def edge_chimeric() -> pd.DataFrame:
    """⭐ CHIMERIC (AGO-ligation) evidence AT ITS NATIVE RUNG — per (gene, arm). (MH-218)

    This was a real gap: chimeric ligation is **the one evidence type that RESOLVES THE MATURE ARM**, it is
    natively per (arm, gene), and only the ARM ROLLUP was carded. So an edge-level question — "does THIS
    pair have a physical duplex?" — had no card answer, which is exactly why MH-217's seedless-strong
    candidates could be counted but not found.

    ⛔ Prefix is `echim_`, NOT `chim_` — the ARM card owns `chim_` for its own per-arm rollup, and
    `PREFIXES` is stripped GLOBALLY by `_annotate`. Adding `chim_` here silently DESTROYED the arm card's
    native 5-column block on the first run (293 cols -> 288). Same convention as `esub_`.

    ⚠⚠ **PER SOURCE, NEVER POOLED.** Four sources with INCOMPATIBLE scales (Manakov weight = reads;
    TarBase = record count) and ~52% pair overlap between Manakov and TarBase qCLASH — a naive sum
    double-counts and lets Manakov dominate. `chim_n_sources` is the honest cross-source column.
    ⚠⚠ **Manakov is ~90% of all rows and is 100% HEK293T — a CELL LINE, not breast. No breast chimeric
    data exists in ANY source** (TarBase CLASH is umbilical vein / stomach / intestine / kidney). This is
    cross-tissue corroboration of a physical duplex, never breast-context validation.
    """
    p = OUT / "chimeric_evidence.parquet"
    if not p.exists():
        return pd.DataFrame()
    d = pd.read_parquet(p)
    d = d.rename(columns={"mature": "arm"})
    d["_grp"] = np.where(d["source"].astype(str).str.startswith("manakov"), "manakov", "tarbase")
    w = (d.groupby(["gene", "arm", "_grp"])["weight"].sum().unstack("_grp")
         .rename(columns={"manakov": "echim_manakov_w", "tarbase": "echim_tarbase_w"}))
    n = d.groupby(["gene", "arm"])["source"].nunique().rename("echim_n_sources")
    out = w.join(n).reset_index()
    # ⛔ PRUNED 2026-08-19 (column review unit 21): `out["echim_any"] = True` — a literal constant on
    # a frame that only ever holds edges WITH chimeric evidence, so it carried **0 bits**: its
    # non-null mask is IDENTICAL to `echim_n_sources`'s (verified on all 1,449 rows) and its only
    # value is True. It was the original exemplar of "a flag that can never say unknown" (MH-256);
    # the honest fix is deletion, not masking — `echim_n_sources.notna()` already answers "any?".
    # ⚠ Do NOT re-add it as a True/False over ALL edges without deciding what a False means:
    # "no chimeric evidence" and "not scanned for chimeric evidence" are different claims.
    return out.drop_duplicates(subset=["gene", "arm"])


def edge_affinity_pct() -> pd.DataFrame:
    """⭐ "IS THIS GENE AMONG THE ARM'S STRONGEST TARGETS?" — persisted at last. (MH-218)

    `kd.genome_affinity_pct()` is the scale-free per-arm specificity percentile over the arm's GENOME-WIDE
    scanMiR K_D targetome, validated (kd_fair_bench) to recover canonical specialists better than context++
    — and it was **computed on the fly and never persisted anywhere**, so no card could answer it.

    ⚠ It is a WITHIN-ARM rank, so it says nothing about absolute affinity: a promiscuous arm's 99th
    percentile and a specialist's 99th percentile are not comparable in K_D. `kd_repression` (the raw
    scanMiR value) ships beside it for exactly that reason.
    ⚠ Coverage is the K_D-scanned arm set, so a NaN means UNSCANNED, not "weak target".
    """
    try:
        from mirna_hallmark.learned import kd as KD
        pct = KD.genome_affinity_pct().rename("kd_affinity_pct").reset_index()
        rep = KD.genome_affinity().rename(columns={"repression": "kd_repression"})
        out = pct.merge(rep, on=["arm", "gene"], how="left")
        return out[["gene", "arm", "kd_affinity_pct", "kd_repression"]].drop_duplicates(
            subset=["gene", "arm"])
    except Exception as ex:
        print(f"  ⚠ kd affinity unavailable ({ex})")
        return pd.DataFrame()


def seedless_chimeric_candidates(strong: float = -0.10) -> pd.DataFrame:
    """⭐ GIVE MH-217's SEEDLESS+STRONG+CHIMERIC EDGES A HOME (MH-218).

    MH-217 reported that ~8% of seedless-but-strongly-coupled edges carry chimeric evidence — a 2.8×
    enrichment over seedless-weak — and identified them as *the only genuinely interesting subset* of the
    seedless class. **But only the COUNT was ever produced; the edges themselves lived in a `/tmp`
    scratchpad and were never written down.** That is precisely the producer-less shape that rotted
    MH-196's five literature sets and MH-201's modifier scan, called out one row earlier and then repeated.

    These edges are the candidate **non-canonical (seed-free) binding** class: no seed match, real
    anti-correlation, and a physical AGO duplex. ⚠ They are a HYPOTHESIS QUEUE, not findings — n is tiny,
    the duplex evidence is HEK293T, and MH-217's own conclusion is that the seedless class as a whole is a
    depleted tail (18.3% vs 39.0% strong-coupling rate) needing no special mechanism.
    """
    if not EDGE_CARD.exists():
        return pd.DataFrame()
    e = pd.read_csv(EDGE_CARD, sep="\t", low_memory=False)
    need = {"adm_has_site", "coupling_tum"}
    if not need <= set(e.columns):
        print("  ⚠ run --annotate first (needs adm_ + chim_ on the edge card)")
        return pd.DataFrame()
    s = e[e["adm_has_site"].notna() & e["coupling_tum"].notna()]
    s = s[~s["adm_has_site"].astype(bool) & (s["coupling_tum"] < strong)]
    if "echim_any" in s.columns:
        s = s[s["echim_any"].fillna(False).astype(bool)]
    keep = [c for c in ("gene", "arm", "seed_family", "coupling_tum", "coupling_p_tum", "beta",
                        "adm_expressed", "echim_manakov_w", "echim_tarbase_w", "echim_n_sources",
                        "kd_affinity_pct", "kd_repression", "retention_rho") if c in s.columns]
    return s[keep].sort_values("coupling_tum")


def _annotate(card_path: Path, blocks: list) -> None:
    """Join EVERY block for this card in ONE call. ADDITIVE ONLY — aborts if a pre-existing column changes.

    ⛔⛔ **`blocks` is a LIST, and that is not cosmetic.** This function strips everything it is about to
    re-add, so calling it twice on the same card makes the second call DELETE the first call's
    columns — measured when exactly that happened here: edge_card went 165→173 (+esub_) and then
    165→166, silently dropping all 8 `esub_` columns. Idempotency and multi-block joins are in tension;
    the resolution is one call per card, never one call per block.
    """
    blocks = [(b, k) for b, k in blocks if b is not None and len(b)]
    if not card_path.exists() or not blocks:
        print(f"  ⚠ absent/empty, skipped: {card_path.name}")
        return
    card = pd.read_csv(card_path, sep="\t", low_memory=False)
    # ⭐⭐ WHAT THIS CALL OWNS IS DERIVED FROM THE BLOCKS, NOT FROM A GLOBAL PREFIX LIST (MH-220).
    # The old form stripped a module-level `PREFIXES` tuple, which is GLOBAL while columns are
    # CARD-SPECIFIC — so adding `chim_` for the edge card silently DELETED the arm card's own native
    # 5-column `chim_` rollup (293 -> 288), because the arm call re-added only `adm_`. Deriving the
    # owned set from the blocks makes that failure UNREPRESENTABLE: this function can now only ever
    # strip columns it is about to re-add, on whatever card it is pointed at.
    owned = {c for b, k in blocks for c in b.columns if c not in k}
    before = card.drop(columns=[c for c in card.columns if c in owned], errors="ignore")
    out = before.copy()
    for block, key in blocks:
        cols = [c for c in block.columns if c in owned or c in key]
        out = out.merge(block[cols].drop_duplicates(subset=key), on=key, how="left")
    # ⭐ THE CONTROL: every pre-existing column bit-identical, same order, same row count
    chk = out[before.columns]
    if not chk.equals(before):
        diff = [c for c in before.columns if not chk[c].equals(before[c])]
        raise SystemExit(f"⛔ ABORT — not additive for {card_path.name}: {diff[:8]}")
    if len(out) != len(before):
        raise SystemExit(f"⛔ ABORT — row count changed for {card_path.name}: {len(before)} -> {len(out)}")
    out.to_csv(card_path, sep="\t", index=False)
    added = sorted(owned & set(out.columns))
    print(f"  ✅ {card_path.name}: {len(before.columns)} -> {out.shape[1]} cols (+{len(added)}), "
          f"{len(out):,} rows, pre-existing columns bit-identical")


def annotate() -> None:
    """Build the ladder/subtype blocks and JOIN them onto the cards.

    ⭐ Callable, not CLI-only (MH-227): this is a SECOND annotation pass, separate from
    `card_context.annotate()`, and only the latter was wired into the card builders by MH-222. Skipping it
    silently costs **37 columns** (20 edge `adm_*`/`echim_*`/`esub_*`, 17 gene `greal_*`/`lit_*`) — measured
    when the tumour-only rebuild dropped them and nothing errored. `realization._finish_card` now calls both.
    """
    _run(annotate_cards=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotate", action="store_true", help="join the blocks onto the cards")
    a = ap.parse_args()
    _run(annotate_cards=a.annotate)


def _run(*, annotate_cards: bool) -> None:
    gl = gene_realization_ladder()
    st = edge_subtype()
    if len(gl):
        n = gl["greal_n_scored"]
        print(f"[gene ladder] {len(gl):,} genes | median regulators scored {n.median():.0f} | "
              f"totals " + " ".join(f"c{abs(int(c*100)):02d}={int(gl[f'greal_n_c{abs(int(c*100)):02d}'].sum()):,}"
                                    for c in COUPLING_CUTS))
        print(f"[gene ladder] genes with >=1 regulator at c10: {int((gl.greal_n_c10 > 0).sum()):,} "
              f"| median gated frac_c10 {gl.greal_frac_c10.median():.3f}")
    if len(st):
        print(f"[edge subtype] {len(st):,} (gene,arm) rows")
    if not annotate_cards:
        print("\n(report only — pass --annotate to join onto the cards)")
        return
    print("\n[annotate]")
    normalise_gene_card()          # rename/prune BEFORE annotating, so the registry sees final names
    normalise_family_card()        # doctrine §4.5: this card has no _finish_card funnel of its own
    normalise_flag_cards()         # MH-256: mask flags that cannot say 'unknown' (edge + arm)
    _annotate(GENE_CARD, [(gl, ["gene"]), (gene_lit_ground_truth(), ["gene"]),
                          (gene_arm_resolution(), ["gene"]), (gene_concentration_adj(), ["gene"]),
                          (discovery_queue_rollup("gene"), ["gene"]),
                          (gene_top_identity_gated(), ["gene"])])
    _annotate(EDGE_CARD, [(st, ["gene", "arm"]), (edge_ago_measured(), ["gene", "arm"]),
                          (edge_admissibility(), ["gene", "arm"]),
                          (edge_chimeric(), ["gene", "arm"]), (edge_affinity_pct(), ["gene", "arm"]),
                          (edge_identity_gate(), ["gene", "arm"]),
                          (edge_hly_leg_concordance(), ["gene", "arm"])])
    # ⭐ the gene_family card gets its identity gate here too — it is the ONLY annotator that touches it,
    # and the gate needs no refit (it is derived from `identity`, see `_identity_gate`).
    _annotate(OUT / "gene_family_card.tsv", [(family_identity_gate(), ["gene", "family"])])
    _annotate(OUT / "arm_card.tsv", [(arm_admissibility_rollup(), ["arm"]),
                                     (family_state_shift(), ["arm"]),
                                     (discovery_queue_rollup("arm"), ["arm"])])
    # ⭐ MH-218: give MH-217's candidates a FILE, not just a count in a chat log.
    cand = seedless_chimeric_candidates()
    if len(cand):
        dest = OUT / "seedless_chimeric_candidates.tsv"
        cand.to_csv(dest, sep="\t", index=False)
        print(f"\n  ⭐ {len(cand)} seedless+strong+chimeric candidate edge(s) -> {dest.name}")
        print(cand.head(8).to_string(index=False))
    print("\n⚠ re-run `card_rungs --check`; and re-run THIS after any gene_card/edge_card rebuild.")


if __name__ == "__main__":
    main()
