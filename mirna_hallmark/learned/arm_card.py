"""⭐ THE ARM CARD — the missing fourth rung, and the one home for every per-ARM property (MH-209).

    .venv/bin/python3 -m mirna_hallmark.learned.arm_card [--rebuild]

**WHY THIS EXISTS.** `card_rungs.py` registered exactly three cards — `edge` (gene,arm), `family`
(gene,family), `gene` (gene) — and **no ARM card**. So every genuinely per-arm property in this subproject
had nowhere to live and ended up scattered across modules, each re-deriving its own: `genomewide_promiscuity`
(targetome breadth), `kd` (scanMiR affinity), `genomic_context` (host locus), `ago_loading` (RISC loading),
`rarity_bench` (seed rarity), `arm_rung` (arm-level β). MH-208 is what that costs: a promiscuity annotation
sat mis-specified for five weeks with one consumer, because there was no rung where it would have been read
next to the abundance and fame columns that explain it.

**⛔⛔ THE THREE TARGETOME UNIVERSES ARE DIFFERENT AND MUST NEVER BE BLENDED.** This is why every column
carries a provenance prefix rather than a friendly name:

| prefix | source | universe | arms covered |
|---|---|---|---|
| `seq_`  | scanMiR RBNS K_D, genome-wide scan | **16,761 breast-expressed genes** | **746** |
| `site_` | scanMiR `kd_repression` site-type counts | ⚠ **HALLMARK-SCOPED, 1,432 genes only** | 771 |
| `ts_`   | TargetScan context++ default predictions | per-SITE, human rows only | ⚠ **321** |
| `cur_` / `fame_` | miRTarBase curation / ledger PMIDs | curated literature | 962 / 2,643 |

A `site_frac_8mer` is therefore *"the 8mer share of this arm's sites among Hallmark genes"*, NOT of its
targetome. Reading it as the latter is the MH-208 error one level down.

**⛔ `cur_*` AND `fame_*` ARE FAME AXES, NOT BIOLOGY (MH-208).** Curated degree correlates with an arm's
distinct-PMID count at **ρ=+0.736** and with abundance at **+0.556**, but with the sequence targetome at only
**+0.124** — their top-10 lists share nothing. They are carried here so a reader can CONTROL for fame, never
so it can be used as targetome breadth.

**⚠ MISSING ≠ ZERO.** An arm absent from a scan is UNSCANNED. Every block ships a `cov_*` boolean; nothing
is imputed. TargetScan covers 321 arms — do not read a NaN `ts_frac_8mer` as "no 8mer sites".

⚠ **SITE-TYPE CODING IS VERIFIED, NOT ASSUMED**: TargetScan `Site Type` 3/2/1 → 8mer / 7mer-m8 / 7mer-A1,
confirmed by the monotone context++ ordering (−0.253 / −0.169 / −0.127 mean), not by the documentation.

**COVERAGE — read it on the denominator that MATTERS** (measured at build):

| arm set | n | `seq_` | `site_` | `ts_` | `cur_` |
|---|---|---|---|---|---|
| all (union of blocks) | 3,241 | 23.0% | 23.8% | 9.9% | 29.7% |
| expressed in TCGA | 2,236 | 32.5% | 31.0% | 13.4% | 37.7% |
| **model arms (≥1 HE edge)** | **803** | **65.8%** | **96.0%** | **33.7%** | **99.0%** |

✅ **DEFINITIONAL CONTROL — the card re-derives MH-208 TO THE DIGIT**, on numbers measured before it
existed: curated↔fame **+0.736**, curated↔abundance **+0.556**, curated↔sequence **+0.124**,
sequence↔fame **+0.197**, sequence↔abundance **−0.004**. An aggregate artifact that cannot reproduce a
known answer is not trustworthy; this one does.

⚠⚠ **TWO THINGS THE DATA CONTRADICTED IN THIS CARD'S OWN DESIGN — do not re-assume them:**
  1. **"SITES vs GENES" IS NOT THE INDEPENDENT AXIS IT LOOKS LIKE.** `site_sites_per_gene` spans only
     **1.00–2.88 (median 1.34)** and correlates with `site_n_genes_canonical` at **ρ=+0.927** — so
     "many sites in few genes" vs "few sites in many genes" barely exists here; the two are ~redundant.
     ⚠ It is also a ratio carrying its own denominator (axiom 5), which would *depress* the correlation,
     so the true association is stronger still. Use `site_n_canonical`; treat per-gene density as thin.
  2. ⛔ **THE 8mer SHARE IS TOOL-DEPENDENT — NEVER COMPARE IT ACROSS SOURCES.** The same concept measured
     two ways agrees at only **ρ=+0.620**, and the medians differ **2.2×** (`site_frac_8mer` 0.150 vs
     `ts_frac_8mer` 0.329) — TargetScan's default predictions are a filtered, 8mer-enriched set. Pick ONE
     source per analysis and say which.
  ✅ Sanity that DOES hold: site-type counts are monotone as biology requires —
     8mer 56,967 < 7mer 305,467 < 6mer 1,282,442 < non-canonical 5,493,916.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "arm_card.tsv"
TS_FILE = C.REPO_ROOT / "data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt"
KD_SITES = C.REPO_ROOT / "data/external_cache/scanmir/kd_repression.tsv.gz"
TS_SITE_TYPE = {3: "8mer", 2: "7mer_m8", 1: "7mer_A1"}   # verified via context++ ordering, see docstring

_MEM: dict = {}


# --------------------------------------------------------------------------- #
# blocks
# --------------------------------------------------------------------------- #
def _abundance() -> pd.DataFrame:
    """TCGA arm expression: LEVEL and DYNAMIC RANGE. Dispersion matters, level alone does not (axiom 8)."""
    from mirna_hallmark.learned import data as LD
    X = LD._load()["X"].apply(pd.to_numeric, errors="coerce")
    return pd.DataFrame({
        "abund_med": X.median(axis=1), "abund_sd": X.std(axis=1),
        "abund_iqr": X.quantile(0.75, axis=1) - X.quantile(0.25, axis=1),
        "abund_frac_detected": (X > 0).mean(axis=1),
    }).rename_axis("arm")


def _sequence_targetome() -> pd.DataFrame:
    """GENOME-WIDE breadth (breast-expressed). The only curation-free promiscuity — see MH-208."""
    from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
    if not GP.SEQ_CACHE.exists():
        GP.build_sequence_promiscuity()
    d = pd.read_csv(GP.SEQ_CACHE, sep="\t").set_index("arm")
    return d.rename(columns={"seq_degree_strong": "seq_n_genes_strong",
                             "seq_degree_any": "seq_n_genes_any",
                             "promisc_seq_strong": "seq_promisc",
                             "promisc_seq_any": "seq_promisc_any"})


def _site_composition() -> pd.DataFrame:
    """⚠ HALLMARK-SCOPED site-type composition (scanMiR `kd_repression`, 1,432 genes only).

    Counts SITES, not genes — the distinction the targetome columns cannot make: an arm with many sites in
    few genes is a different object from one with single sites in many genes (`site_sites_per_gene`).
    """
    d = pd.read_csv(KD_SITES, sep="\t")
    n = d.groupby("arm")[["n_8mer", "n_7mer", "n_6mer", "n_noncanon"]].sum()
    n.columns = ["site_n_8mer", "site_n_7mer", "site_n_6mer", "site_n_noncanon"]
    canon = n["site_n_8mer"] + n["site_n_7mer"]
    tot = canon + n["site_n_6mer"] + n["site_n_noncanon"]
    n["site_n_canonical"] = canon
    n["site_n_total"] = tot
    n["site_frac_8mer"] = (n["site_n_8mer"] / canon.replace(0, np.nan))       # ⚠ gated denominator (axiom 5)
    n["site_frac_canonical"] = canon / tot.replace(0, np.nan)
    g = d[(d.n_8mer + d.n_7mer) > 0].groupby("arm")["gene"].nunique()
    n["site_n_genes_canonical"] = g.reindex(n.index).fillna(0).astype(int)
    n["site_sites_per_gene"] = canon / n["site_n_genes_canonical"].replace(0, np.nan)
    r = d.groupby("arm")["repression"]
    n["site_repression_med"] = r.median()
    n["site_repression_min"] = r.min()                                        # the arm's strongest edge
    return n.rename_axis("arm")


def _targetscan() -> pd.DataFrame:
    """⚠ TargetScan per-SITE context++ (321 human arms). Site type verified by context++ ordering."""
    d = pd.read_csv(TS_FILE, sep="\t", low_memory=False,
                    usecols=["Gene Tax ID", "Gene Symbol", "miRNA", "Site Type",
                             "context++ score", "Predicted relative KD"])
    d = d[(d["Gene Tax ID"] == 9606) & d["miRNA"].str.startswith("hsa-", na=False)]
    d = d[d["Site Type"].isin(TS_SITE_TYPE)]
    d["ctx"] = pd.to_numeric(d["context++ score"], errors="coerce")
    d["kd"] = pd.to_numeric(d["Predicted relative KD"], errors="coerce")
    g = d.groupby("miRNA")
    out = pd.DataFrame({
        "ts_n_sites": g.size(),
        "ts_n_genes": g["Gene Symbol"].nunique(),
        "ts_ctx_mean": g["ctx"].mean(),
        "ts_ctx_best": g["ctx"].min(),          # most negative = strongest predicted repression
        "ts_kd_med": g["kd"].median(),
    })
    for code, name in TS_SITE_TYPE.items():
        out[f"ts_n_{name}"] = d[d["Site Type"] == code].groupby("miRNA").size().reindex(out.index).fillna(0).astype(int)
    out["ts_frac_8mer"] = out["ts_n_8mer"] / out["ts_n_sites"].replace(0, np.nan)
    out["ts_sites_per_gene"] = out["ts_n_sites"] / out["ts_n_genes"].replace(0, np.nan)
    return out.rename_axis("arm")


def _curation_and_fame() -> pd.DataFrame:
    """⛔ FAME AXES (MH-208). Carried to CONTROL for study bias — never as targetome breadth."""
    from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
    from mirna_hallmark.learned.evidence import ledger as LG
    cur = pd.read_csv(GP.CACHE, sep="\t").set_index("arm")[["he_degree", "he_degree_expr"]]
    cur.columns = ["cur_he_degree", "cur_he_degree_expr"]
    L = LG.build_ledger()
    fame = pd.DataFrame({"fame_npmid": L.groupby("arm")["pmid"].nunique(),
                         "fame_n_genes_curated": L.groupby("arm")["gene"].nunique()})
    return cur.join(fame, how="outer").rename_axis("arm")


def _locus_and_loading() -> pd.DataFrame:
    """Host-locus class (strand-aware) + AGO/RISC loading — both genuinely per-arm."""
    blocks = []
    try:
        from mirna_hallmark.learned import genomic_context as GC
        hm = GC.locus_host_map()
        rows = [{"arm": a, "host_class": v.get("mir_class") if isinstance(v, dict) else v,
                 "clustered": (v.get("clustered") if isinstance(v, dict) else np.nan),
                 "mirgenedb": (v.get("mirgenedb") if isinstance(v, dict) else np.nan)}
                for a, v in hm.items()]
        blocks.append(pd.DataFrame(rows).set_index("arm"))
    except Exception as e:
        print(f"  ⚠ genomic_context unavailable ({e}) — host_class omitted")
    try:
        ago = pd.read_parquet(OUT / "ago_loading.parquet")
        blocks.append(ago.rename(columns={"ago_reads": "ago_reads"}).rename_axis("arm"))
    except Exception as e:
        print(f"  ⚠ ago_loading unavailable ({e})")
    if not blocks:
        return pd.DataFrame().rename_axis("arm")
    out = blocks[0]
    for b in blocks[1:]:
        out = out.join(b, how="outer")
    return out.rename_axis("arm")


def _model_footprint() -> pd.DataFrame:
    """How much of the LEARNED MODEL this arm actually touches — the design, not the literature."""
    from mirna_hallmark.learned.evidence import ledger as LG
    he = LG.pooled_he_edges()
    g = he.groupby("miRNA")
    return pd.DataFrame({"model_n_genes": g["gene"].nunique(),
                         "model_n_edges": g.size()}).rename_axis("arm")


# --------------------------------------------------------------------------- #
def build() -> pd.DataFrame:
    print("[arm_card] building blocks…")
    blocks = {"abundance": _abundance(), "sequence targetome": _sequence_targetome(),
              "site composition": _site_composition(), "targetscan": _targetscan(),
              "curation/fame": _curation_and_fame(), "locus/loading": _locus_and_loading(),
              "model footprint": _model_footprint()}
    card = None
    for name, b in blocks.items():
        if b is None or not len(b):
            print(f"  ⚠ {name}: EMPTY, skipped")
            continue
        b = b[~b.index.duplicated(keep="first")]
        print(f"  {name:20s}: {len(b):5,d} arms × {b.shape[1]:2d} cols")
        card = b if card is None else card.join(b, how="outer")

    # seed family + size (the family rung's key, so the two cards join cleanly)
    try:
        from mirna_hallmark.coupling_inference import seed_family_map
        fam = pd.Series(seed_family_map(list(card.index)), name="seed_family")
        card = card.join(fam)
        # ⚠ NOT `family_size` — the EDGE card already has that name at the FAMILY rung, counting arms in
        # the MODEL DESIGN. This counts arms in THIS CARD's universe (every arm any block knows about),
        # a different denominator. Same name, different unit is exactly axiom 6's collision, and the
        # `--check` cross-card report flagged it on the first run. Renamed rather than footnoted.
        card["family_n_arms_card"] = card.groupby("seed_family")["seed_family"].transform("size")
    except Exception as e:
        print(f"  ⚠ seed family map unavailable ({e})")

    # ⚠ COVERAGE FLAGS — a missing block means UNSCANNED, never zero.
    card["cov_seq_scan"] = card.get("seq_n_genes_strong", pd.Series(index=card.index)).notna()
    card["cov_site_types"] = card.get("site_n_total", pd.Series(index=card.index)).notna()
    card["cov_targetscan"] = card.get("ts_n_sites", pd.Series(index=card.index)).notna()
    card["cov_curated"] = card.get("cur_he_degree", pd.Series(index=card.index)).notna()

    card = card.rename_axis("arm").reset_index()
    DEST.parent.mkdir(parents=True, exist_ok=True)
    card.to_csv(DEST, sep="\t", index=False)
    print(f"\n[arm_card] {len(card):,} arms × {card.shape[1]} columns -> {DEST}")
    for c in ("cov_seq_scan", "cov_site_types", "cov_targetscan", "cov_curated"):
        print(f"    {c:18s}: {int(card[c].sum()):5,d} arms ({100*card[c].mean():.1f}%)")
    return card


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rebuild", action="store_true", help="rebuild the sequence-promiscuity cache too")
    a = ap.parse_args()
    if a.rebuild:
        from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
        GP.build_sequence_promiscuity()
    build()


if __name__ == "__main__":
    main()
