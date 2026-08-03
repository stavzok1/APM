"""Cross-card LADDERS — closes the three rung gaps the MH-212 cross-card audit found.

    .venv/bin/python3 -m mirna_hallmark.learned.card_ladders            # report only
    .venv/bin/python3 -m mirna_hallmark.learned.card_ladders --annotate # + join onto the cards

**WHY.** MH-212 put a realized-coupling threshold LADDER on the arm card ("how many genes does this arm
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

⚠ **`ago_loading_measured` is the ONE annotated column without a module prefix** — deliberately, because
it must sit adjacent to the column it qualifies. It is added to `PREFIXES` explicitly so the idempotent
re-run still strips it.

⚠ **ADDITIVE ONLY, and idempotent.** `_annotate` strips this module's own prefixes before re-joining (so
a re-run replaces rather than appending `_x`/`_y`), then asserts **every pre-existing column is
bit-identical and the row count is unchanged** — aborting otherwise. Same contract as
`card_context._annotate`, which established the pattern.

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
PREFIXES = ("greal_", "esub_", "ago_loading_measured")

# the SAME cuts the arm card uses — a ladder is only comparable across rungs if the rungs agree on it
COUPLING_CUTS = (-0.05, -0.10, -0.15, -0.20, -0.30)


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


def _annotate(card_path: Path, blocks: list) -> None:
    """Join EVERY block for this card in ONE call. ADDITIVE ONLY — aborts if a pre-existing column changes.

    ⛔⛔ **`blocks` is a LIST, and that is not cosmetic.** This function strips ALL of `PREFIXES` before
    re-joining, so calling it twice on the same card makes the second call DELETE the first call's
    columns — measured when exactly that happened here: edge_card went 165→173 (+esub_) and then
    165→166, silently dropping all 8 `esub_` columns. Idempotency and multi-block joins are in tension;
    the resolution is one call per card, never one call per block.
    """
    blocks = [(b, k) for b, k in blocks if b is not None and len(b)]
    if not card_path.exists() or not blocks:
        print(f"  ⚠ absent/empty, skipped: {card_path.name}")
        return
    card = pd.read_csv(card_path, sep="\t", low_memory=False)
    # ⚠ strip this module's prefixes first, so a RE-RUN replaces instead of appending `_x`/`_y`
    before = card.drop(columns=[c for c in card.columns if c.startswith(PREFIXES)], errors="ignore")
    out = before.copy()
    for block, key in blocks:
        cols = [c for c in block.columns if c.startswith(PREFIXES) or c in key]
        out = out.merge(block[cols].drop_duplicates(subset=key), on=key, how="left")
    # ⭐ THE CONTROL: every pre-existing column bit-identical, same order, same row count
    chk = out[before.columns]
    if not chk.equals(before):
        diff = [c for c in before.columns if not chk[c].equals(before[c])]
        raise SystemExit(f"⛔ ABORT — not additive for {card_path.name}: {diff[:8]}")
    if len(out) != len(before):
        raise SystemExit(f"⛔ ABORT — row count changed for {card_path.name}: {len(before)} -> {len(out)}")
    out.to_csv(card_path, sep="\t", index=False)
    added = [c for c in out.columns if c.startswith(PREFIXES)]
    print(f"  ✅ {card_path.name}: {len(before.columns)} -> {out.shape[1]} cols (+{len(added)}), "
          f"{len(out):,} rows, pre-existing columns bit-identical")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotate", action="store_true", help="join the blocks onto the cards")
    a = ap.parse_args()

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
    if not a.annotate:
        print("\n(report only — pass --annotate to join onto the cards)")
        return
    print("\n[annotate]")
    _annotate(GENE_CARD, [(gl, ["gene"])])
    _annotate(EDGE_CARD, [(st, ["gene", "arm"]), (edge_ago_measured(), ["gene", "arm"])])
    print("\n⚠ re-run `card_rungs --check`; and re-run THIS after any gene_card/edge_card rebuild.")


if __name__ == "__main__":
    main()
