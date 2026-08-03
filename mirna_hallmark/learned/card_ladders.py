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
PREFIXES = ("greal_", "esub_", "ago_loading_measured", "adm_", "lit_", "echim_", "kd_")

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
    out["echim_any"] = True
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
    _annotate(GENE_CARD, [(gl, ["gene"]), (gene_lit_ground_truth(), ["gene"])])
    _annotate(EDGE_CARD, [(st, ["gene", "arm"]), (edge_ago_measured(), ["gene", "arm"]),
                          (edge_admissibility(), ["gene", "arm"]),
                          (edge_chimeric(), ["gene", "arm"]), (edge_affinity_pct(), ["gene", "arm"])])
    _annotate(OUT / "arm_card.tsv", [(arm_admissibility_rollup(), ["arm"])])
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
