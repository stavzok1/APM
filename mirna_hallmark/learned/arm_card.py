"""⭐ THE ARM CARD — the fourth rung, and the one home for every per-ARM property (MH-209, v2 MH-211).

    .venv/bin/python3 -m mirna_hallmark.learned.arm_card [--rebuild]

**WHY THIS EXISTS.** `card_rungs.py` registered exactly three cards — `edge` (gene,arm), `family`
(gene,family), `gene` (gene) — and **no ARM card**. So every genuinely per-arm property had nowhere to
live and each module re-derived its own: `genomewide_promiscuity` (targetome breadth), `kd` (affinity),
`genomic_context` (host locus), `ago_loading` (RISC loading), `seed_rarity`, `arm_rung` (arm-level β).
MH-208 is what that costs: a promiscuity annotation sat mis-specified for five weeks with one consumer,
because there was no rung where it would have been read next to the abundance and fame columns that
explain it.

## ⛔ v1 DEFECT, FIXED HERE (MH-211) — READ BEFORE TRUSTING ANY v1 NUMBER

v1's locus block called `genomic_context.locus_host_map()`, which is keyed by **locus_id (`MI*`
precursor accession), NOT by arm**, and whose dicts carry no `clustered`/`mirgenedb` keys. Measured on
the shipped v1 card: **3,241 rows = 2,730 real arms + 506 `MI0*` rows + 5 unicode-hyphen rows**;
`host_class` non-null = 506 **entirely on the bogus rows**; `clustered` and `mirgenedb` **0 non-null**.
⇒ v1's published coverage percentages were computed on an inflated denominator and three columns were
dead. Fixed by reading the real per-arm artifact `output/learned/genomic_context.tsv` (714 arms × 20
cols) and by gating the key space (`_arm_key` + the `build()` assertions).
⚠ The **15 `MIMAT*` rows are REAL arms** (15/15 carry `abund_med`), merely unnamed in miRBase — they
stay. **The true denominator is 2,730.**

## ⛔⛔ THE FOUR TARGETOME UNIVERSES ARE DIFFERENT OBJECTS — NEVER BLEND THEM

This is why every column carries a provenance prefix rather than a friendly name.

| prefix | source | universe | arms |
|---|---|---|---|
| `seq_`  | scanMiR RBNS K_D, genome-wide scan | 16,761 breast-expressed genes | 746 |
| `site_` | scanMiR `kd_repression` site-type counts | ⚠ **HALLMARK-SCOPED, 1,432 genes only** | 771 |
| `ts_`   | TargetScan context++ default predictions | per-SITE, human rows only | ⚠ 321 |
| `sdsz_` | ⭐ **NEW (v2)** `seed_rarity` MANE 3′UTR 7/8-mer scan | all MANE 3′UTRs, seed-keyed | 2,656 |

A `site_frac_8mer` is *"the 8mer share of this arm's sites among Hallmark genes"*, NOT of its targetome.

## PREFIX TAXONOMY — three orthogonal facets, all greppable

| token | meaning | the test that proves it |
|---|---|---|
| *(none)* | arm-native measured | — |
| `bc_` | **BROADCAST** from a coarser unit (seed / family / hairpin) — inherited, not measured here | constant within its broadcast key |
| `hx_` | **HEURISTIC / retired-estimator** derivation — may not sit unlabelled beside a learned β | provenance grep of the producing call |
| `cur_` / `fame_` | ⛔ **FAME AXES** — curation depth, carried to CONTROL for study bias, never as biology | MH-208: ρ=+0.736 with #PMIDs, +0.124 with the sequence targetome |
| `cov_<block>` | boolean; **MISSING = UNSCANNED, never zero. Nothing is imputed.** | block columns all-NaN outside it |

`bc_` is **not** redundant with `card_rungs.AGG_OF`: AGG_OF records *"summarises a LOWER unit"*, `bc_`
records *"inherited from a HIGHER unit"*. Opposite directions.

## PROVENANCE AUDIT — what was deliberately EXCLUDED (do not re-litigate)

- ⛔ **`mirna_state_class.tsv` contributes NOTHING.** Its `mass_X(m) = Σ_g w(m,g)·f_softmax_logrpm`
  (`mirna_state_class.py:27`) is the **§6b-RETIRED pressure heuristic (MH-115)**, so `share_*`/`rank_*`
  and everything downstream (`dHN/dNT/dHT`, `primary_class`, `acquired_*`) are heuristic-derived. And
  `subtype_tau` is `_tau()` computed **over `subshare_*`, which are shares of `mass_*`** ⇒ also
  heuristic. The whole table fails the audit.
- ⛔ **`comovement`'s `dHT` / `primary_class` / `trajectory_cluster`** — `dHT` is the exact quantity
  **MH-210 proved was being fabricated by the GTEx multi-mapping collapse**. Only the pure
  co-expression columns are taken.
- ⛔ **`arm_survival.gse_*`** — its producer matches arms **by base name with no arm suffix**; that is
  the identical bare-stem failure that makes Buffa's 501 probes unusable at arm rung.
- ⛔ **`gtex_qn_baseline.floor0` is never emitted as a value** (1,137 arms) — that 0-vs-NaN collapse is
  precisely what MH-210 removed from `budget_shift`; re-importing it reinstates the bug one rung up.
- ⛔ **`ago_loading.loading()` is NOT used** — it ends `.fillna(1.0)` ("uncovered → 1.0"), fabricating a
  guide call for every unmeasured arm. `ago_dom` preserves NaN.
- ⛔ **`arm_variability.tsv` skipped** — verified byte-identical to `abund_med/sd/iqr` (max|d| = 0,
  n=2,236). A tertile is derived in-card instead of importing a second copy.
- ⛔ **`arm_id_status` not lifted** from the edge card — it varies within arm (263/736), so it is
  arm-in-family rung despite sitting in `_ARM_DESC`.

## COVERAGE — on the denominator that MATTERS (axiom 5)

v2 = **2,450 arms**, of which **2,236 = the ENTIRE TCGA matrix** (no expressed arm is lost).

| arm set | n | `site_` | `cur_` | `sdsz_` | `arb_` | `gctx_` | `cnv_` | `seq_` | `foc_` | `ago_` | `ts_` | `bc_meth_` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | 2,450 | 29% | 36% | 97% | 30% | 29% | 35% | 30% | 30% | 28% | 13% | 14% |
| expressed | 2,236 | 31% | 38% | 98% | 32% | 32% | 37% | 33% | 33% | 28% | 14% | 15% |
| **model (≥1 HE edge)** | **741** | **97%** | **99%** | **99%** | **98%** | **96%** | **73%** | **70%** | **68%** | **47%** | **38%** | **32%** |

## VERIFICATION — what was pre-registered, and what the data said back

✅ **no-impute** (0 values outside any `cov_`) · ✅ **key-join cross-check** — `foc_partial_rho` vs
`cnvc_partial_rho` **ρ=+0.982** (n=733): two CN measures of the same link agree, so the CNV joins are
right · ✅ **broadcast constancy** — `bc_meth_*` / `bc_fam_*` take exactly one value per hairpin /
family · ⚠ **anti-blend**: the four universes correlate but are not interchangeable (seed-scan vs
scanMiR +0.825, vs TargetScan +0.645; scanMiR-gw vs Hallmark-sites +0.799).

⚠ **THE DEFINITIONAL CONTROL WAS REFUTED IN THE LETTER — reported, not buried.** Four of MH-208's five
correlations reproduce to the digit (abundance +0.556, sequence +0.124, seq↔fame +0.197, seq↔abundance
−0.004), but **curated↔fame moved +0.736 → +0.752**. Diagnosed rather than waved away: v1's universe
held **76 bare-stem curation rows** (`hsa-miR-1183`, `hsa-let-7a`, …) that are not arms and that v2's
key gate drops. On v1's raw keys the value reproduces **exactly (+0.7361, n=955)**; on the arm-only
universe it is **+0.7518 (n=879)**; and **the dropped rows alone carry only +0.370**. ⇒ the stems were
*diluting* the estimate, so **MH-208's conclusion is STRENGTHENED, not weakened.**

⛔ **A NAME PATTERN IS NOT THE ARM UNIVERSE** (learned twice while building this). A first cut of
`_arm_key` required a `-(3p|5p)$` suffix and silently lost **582 of 2,236 real matures** — the
canonical index is 1,654 suffixed + 582 unsuffixed. And the authority must be consulted **before** any
rewriting: `hsa-miR-496.1` is a real canonical name that `.N`-stripping destroys.
"""
from __future__ import annotations

import argparse
import re
import unicodedata

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

OUT = C.REPO_ROOT / "mirna_hallmark/output/learned"
DEST = OUT / "arm_card.tsv"
TS_FILE = C.REPO_ROOT / "data/miRNA/Predicted_Targets_Context_Scores.default_predictions.txt"
KD_SITES = C.REPO_ROOT / "data/external_cache/scanmir/kd_repression.tsv.gz"
TS_SITE_TYPE = {3: "8mer", 2: "7mer_m8", 1: "7mer_A1"}   # verified via context++ ordering, see below

_MEM: dict = {}
_DROPS: dict = {}     # source -> {reason: [keys]} — printed at build, so drops are never silent


# --------------------------------------------------------------------------- #
# THE KEY NORMALIZER — load-bearing; six sources half-join without it
# --------------------------------------------------------------------------- #
def _mimat_map() -> dict:
    if "mimat" not in _MEM:
        try:
            from analysis.expression.mirna_target_integration import load_mimat_to_arm
            _MEM["mimat"] = load_mimat_to_arm(C.MIRNA_MATURE_LOCI)
        except Exception as e:                                   # never fatal
            print(f"  ⚠ MIMAT→arm bridge unavailable ({e}) — MIMAT-keyed sources will under-join")
            _MEM["mimat"] = {}
    return _MEM["mimat"]


_ARM_RE = re.compile(r"^hsa-(?:miR|let|mir)-.+-(?:3p|5p)$", re.I)


def _universe() -> set:
    """⭐ THE AUTHORITY for what counts as an arm: the TCGA mature-miRNA matrix index (2,236).

    ⚠ MEASURED, and it is why a name-pattern rule alone is wrong: the canonical universe is
    **1,654 arm-suffixed + 582 NOT suffixed** — 582 real mature molecules (unnamed `MIMAT*`
    accessions and matures miRBase never gave a -3p/-5p) that a `-(3p|5p)$` regex silently rejects.
    A first cut of this function did exactly that and lost 25% of the card.
    """
    if "universe" not in _MEM:
        try:
            from mirna_hallmark.learned import data as LD
            _MEM["universe"] = {str(a) for a in LD._load()["X"].index}
        except Exception as e:
            print(f"  ⚠ canonical arm universe unavailable ({e}) — falling back to the name pattern")
            _MEM["universe"] = set()
    return _MEM["universe"]


def _arm_key(x, *, src: str = "?"):
    """Normalize any source's miRNA label onto the card's arm key, or return None with a logged reason.

    Order: normalize → **canonical-universe membership wins** → MIMAT bridge → arm-suffix pattern.

    THREE REJECTION classes, all deliberate and all logged:
      * **non-human** (`ebv-`, `gga-`, `sja-` appear in miRTarBase) — not this cohort's biology.
      * ⛔ **bare stem with no -5p/-3p AND not in the universe** (`hsa-miR-4723`, `MiR-210`) — the two
        arms of a hairpin are different molecules with different seeds. Broadcasting a stem-level value
        to both is the exact failure that makes Buffa's 501 probes unusable at arm rung and that
        disqualifies `arm_survival.gse_*`. ⚠ The universe check must come FIRST, or genuinely
        unsuffixed matures get caught by this rule.
      * **non-ASCII residue** — a key that still holds a non-ASCII character after normalisation is a
        source encoding artifact (miRTarBase ships `MicroRNA‑3907` with U+2011), not an arm.
    """
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    s = unicodedata.normalize("NFKC", str(x)).strip()
    s = re.sub(r"[‐-―−]", "-", s)          # unicode dashes → ASCII hyphen
    if not s or s.lower() == "nan":
        return None
    # ⭐⭐ THE AUTHORITY IS CONSULTED BEFORE ANY REWRITING. A canonical name may itself contain the
    # patterns the rewrites target: the TCGA universe really does hold `hsa-miR-496.1`, and stripping
    # `.N` first turned it into a non-existent `hsa-miR-496` that was then dropped as a bare stem —
    # from the abundance block, i.e. from the universe's own source. Check first, rewrite second.
    if s in _universe():
        return s
    # ⭐ TargetScan multi-seed variants (`hsa-miR-101-3p.2`) — the `.N` suffix names an alternative seed
    # of the SAME arm. Memory `universe-redefinition-pending-refresh` measured that normalising it
    # recovers ~18 HE arms and ~300 edges; without it TargetScan silently loses 20 arms here.
    s = re.sub(r"\.\d+$", "", s)
    # `learned/methylation._hairpin_loci` keys arms WITHOUT a species prefix and lower-cased
    # (`mir-200b-5p`) — un-prefixed, all 353 of its keys were being rejected as non-human.
    if re.match(r"^(mir|let)-", s, flags=re.I):
        s = "hsa-" + s
    if s in _universe():
        return s
    if s.upper().startswith("MIMAT"):
        m = _mimat_map().get(s) or _mimat_map().get(s.upper())
        return m if m else _drop(src, "unmapped_MIMAT", s)
    s = re.sub(r"^hsa-", "hsa-", s, flags=re.I)           # case-fix the species prefix
    s = re.sub(r"^(hsa-)(microrna|mir|miR)-", r"\1miR-", s, flags=re.I)
    if any(ord(ch) > 127 for ch in s):
        return _drop(src, "non_ascii_residue", s)
    if not s.lower().startswith("hsa-"):
        return _drop(src, "non_human", s)
    if s in _universe():
        return s
    if not _ARM_RE.match(s):
        return _drop(src, "bare_stem_no_arm", s)
    return s


def _drop(src: str, reason: str, val: str):
    _DROPS.setdefault(src, {}).setdefault(reason, set()).add(val)
    return None


def _key_index(df: pd.DataFrame, col: str, *, src: str) -> pd.DataFrame:
    """Re-key a frame onto normalized arm keys, dropping unmappable rows (logged) and duplicates."""
    k = df[col].map(lambda v: _arm_key(v, src=src))
    out = df.loc[k.notna()].copy()
    out.index = pd.Index(k[k.notna()].to_numpy(), name="arm")
    return out.drop(columns=[col], errors="ignore")[~out.index.duplicated(keep="first")]


def _norm_index(b: pd.DataFrame, *, src: str) -> pd.DataFrame:
    """Rewrite a block's index THROUGH `_arm_key` and drop what does not survive.

    ⚠ This is the single place the key space is enforced. v1's gate tested whether an index entry
    *could* be normalised but then kept the **original** string — so a unicode-hyphen name that mapped
    onto a valid arm was retained in its raw form and tripped the ASCII assertion. Normalise here, once,
    for every block; collapse any duplicates normalisation creates (first non-null wins).
    """
    k = [_arm_key(a, src=src) for a in b.index]
    keep = [i for i, v in enumerate(k) if v is not None]
    out = b.iloc[keep].copy()
    out.index = pd.Index([k[i] for i in keep], name="arm")
    return out.groupby(level=0).first() if out.index.has_duplicates else out


def _read(path, **kw) -> pd.DataFrame:
    """Module-level memoized read — axiom 3a.1: every static file is read ONCE, never per arm."""
    key = str(path)
    if key not in _MEM:
        _MEM[key] = pd.read_parquet(path) if str(path).endswith(".parquet") else pd.read_csv(path, **kw)
    return _MEM[key]


# --------------------------------------------------------------------------- #
# BLOCKS — each returns an arm-indexed frame
# --------------------------------------------------------------------------- #
def _abundance() -> pd.DataFrame:
    """TCGA arm expression: LEVEL and DYNAMIC RANGE. Dispersion predicts where level does not (axiom 8)."""
    from mirna_hallmark.learned import data as LD
    X = LD._load()["X"].apply(pd.to_numeric, errors="coerce")
    d = pd.DataFrame({
        "abund_med": X.median(axis=1), "abund_sd": X.std(axis=1),
        "abund_iqr": X.quantile(0.75, axis=1) - X.quantile(0.25, axis=1),
        "abund_frac_detected": (X > 0).mean(axis=1),
    })
    d["abund_sd_tertile"] = pd.qcut(d["abund_sd"], 3, labels=["low", "mid", "high"])
    return d.rename_axis("arm")


def _sequence_targetome() -> pd.DataFrame:
    """GENOME-WIDE breadth (breast-expressed). The curation-free promiscuity — MH-208."""
    from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
    if not GP.SEQ_CACHE.exists():
        GP.build_sequence_promiscuity()
    d = _read(GP.SEQ_CACHE, sep="\t").set_index("arm")
    return d.rename(columns={"seq_degree_strong": "seq_n_genes_strong",
                             "seq_degree_any": "seq_n_genes_any",
                             "promisc_seq_strong": "seq_promisc",
                             "promisc_seq_any": "seq_promisc_any"})


def _site_composition() -> pd.DataFrame:
    """⚠ HALLMARK-SCOPED site-type composition (scanMiR `kd_repression`, 1,432 genes only).

    Counts SITES, not genes. ⚠ v1 measured that `site_sites_per_gene` spans only 1.00–2.88 (median
    1.34) and correlates with `site_n_genes_canonical` at ρ=+0.927 — "many sites in few genes" barely
    exists here, so prefer `site_n_canonical` and treat per-gene density as thin.
    """
    d = _read(KD_SITES, sep="\t")
    n = d.groupby("arm")[["n_8mer", "n_7mer", "n_6mer", "n_noncanon"]].sum()
    n.columns = ["site_n_8mer", "site_n_7mer", "site_n_6mer", "site_n_noncanon"]
    canon = n["site_n_8mer"] + n["site_n_7mer"]
    tot = canon + n["site_n_6mer"] + n["site_n_noncanon"]
    n["site_n_canonical"], n["site_n_total"] = canon, tot
    n["site_frac_8mer"] = n["site_n_8mer"] / canon.replace(0, np.nan)        # gated denominator (axiom 5)
    n["site_frac_canonical"] = canon / tot.replace(0, np.nan)
    g = d[(d.n_8mer + d.n_7mer) > 0].groupby("arm")["gene"].nunique()
    n["site_n_genes_canonical"] = g.reindex(n.index).fillna(0).astype(int)
    n["site_sites_per_gene"] = canon / n["site_n_genes_canonical"].replace(0, np.nan)
    r = d.groupby("arm")["repression"]
    n["site_repression_med"], n["site_repression_min"] = r.median(), r.min()
    return n.rename_axis("arm")


def _targetscan() -> pd.DataFrame:
    """⚠ TargetScan per-SITE context++ (321 human arms).

    ⚠ Site-type coding is VERIFIED, not assumed: `Site Type` 3/2/1 → 8mer / 7mer-m8 / 7mer-A1, confirmed
    by the monotone context++ ordering (−0.253 / −0.169 / −0.127 mean), not by the documentation.
    ⛔ The 8mer SHARE is TOOL-DEPENDENT — `site_frac_8mer` vs `ts_frac_8mer` agree at only ρ=+0.620 with
    medians 2.2× apart. Never compare an 8mer fraction across sources; pick one and say which.
    """
    d = _read(TS_FILE, sep="\t", low_memory=False,
              usecols=["Gene Tax ID", "Gene Symbol", "miRNA", "Site Type",
                       "context++ score", "Predicted relative KD"])
    d = d[(d["Gene Tax ID"] == 9606) & d["miRNA"].str.startswith("hsa-", na=False)]
    d = d[d["Site Type"].isin(TS_SITE_TYPE)].copy()
    d["ctx"] = pd.to_numeric(d["context++ score"], errors="coerce")
    d["kd"] = pd.to_numeric(d["Predicted relative KD"], errors="coerce")
    g = d.groupby("miRNA")
    out = pd.DataFrame({"ts_n_sites": g.size(), "ts_n_genes": g["Gene Symbol"].nunique(),
                        "ts_ctx_mean": g["ctx"].mean(), "ts_ctx_best": g["ctx"].min(),
                        "ts_kd_med": g["kd"].median()})
    for code, name in TS_SITE_TYPE.items():
        out[f"ts_n_{name}"] = (d[d["Site Type"] == code].groupby("miRNA").size()
                               .reindex(out.index).fillna(0).astype(int))
    out["ts_frac_8mer"] = out["ts_n_8mer"] / out["ts_n_sites"].replace(0, np.nan)
    out["ts_sites_per_gene"] = out["ts_n_sites"] / out["ts_n_genes"].replace(0, np.nan)
    return out.rename_axis("arm")


def _seed_scan() -> pd.DataFrame:
    """⭐ THE FOURTH UNIVERSE — MANE 3′UTR 7/8-mer seed scan (`learned/seed_rarity.py`), 2,656 arms.

    Not `seq_`, not `site_`, not `ts_`. Effective size `N_eff = N_8mer + 0.5·(N_7mer_plus − N_8mer)`.
    ⚠ `seed_rarity.rarity()` min-max scales over WHATEVER LIST IS PASSED, so its value is not comparable
    across calls (axiom 5). We persist the ABSOLUTE sizes and scale once, here, over the full universe.
    """
    p = OUT / "seed_targetome_size.parquet"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p)
    d = _key_index(d, "arm", src="seed_targetome")
    out = pd.DataFrame(index=d.index)
    out["sdsz_N_7mer_plus"] = pd.to_numeric(d["N_7mer_plus"], errors="coerce")
    out["sdsz_N_8mer"] = pd.to_numeric(d["N_8mer"], errors="coerce")
    out["sdsz_N_eff"] = out["sdsz_N_8mer"] + 0.5 * (out["sdsz_N_7mer_plus"] - out["sdsz_N_8mer"])
    r = 1.0 / np.log1p(out["sdsz_N_eff"].clip(lower=0))
    out["sdsz_rarity"] = (r - r.min()) / (r.max() - r.min()) if r.notna().any() else np.nan
    return out


def _genomic_context() -> pd.DataFrame:
    """⭐ THE v1 FIX — the real per-arm locus table (714 arms), not the locus-keyed `locus_host_map()`."""
    p = OUT / "genomic_context.tsv"
    if not p.exists():
        try:
            from mirna_hallmark.learned import genomic_context as GC
            GC.classify_he_arms()
        except Exception as e:
            print(f"  ⚠ genomic_context unavailable ({e})")
            return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="genomic_context")
    keep = ["mirgenedb", "coord_source", "n_loci", "mir_class", "host", "host_type", "strand_rel",
            "in_exon", "host_region", "context_mixed", "clustered", "promoter_class",
            "promoter_gene", "promoter_gene_type", "promoter_coord_class", "promoter_host_tss_kb"]
    d = d[[c for c in keep if c in d.columns]]
    return d.rename(columns={c: f"gctx_{c}" for c in d.columns})


def _expression_tier() -> pd.DataFrame:
    """Detectability tier (RPM≥10 floor). ⛔ `median_log2rpm` is byte-identical to `abund_med` — dropped."""
    p = C.REPO_ROOT / "mirna_hallmark/output/matrices/arm_expression_tiers.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="expression_tiers")
    return d[["tier", "expressed", "frac_expressed", "p90_log2rpm", "max_log2rpm"]].rename(
        columns={"tier": "tier_class", "expressed": "tier_expressed",
                 "frac_expressed": "tier_frac_expressed", "p90_log2rpm": "tier_p90",
                 "max_log2rpm": "tier_max"})


def _within_patient_shift() -> pd.DataFrame:
    """Within-PATIENT paired NAT→tumour dose shift (2,236 arms — supersedes the 736-arm edge-card lift)."""
    p = OUT / "realization/dose_shift_arm.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="dose_shift_arm")
    return d.rename(columns={c: f"wshift_{c}" for c in d.columns})


def _isomir() -> pd.DataFrame:
    """5′-isomiR seed shift — the arm's true seed vs its family's.

    ⚠⚠ AXIOM 5: the median is n-DEPENDENT (all 687 arms 0.056 → n≥300 0.037). `iso_n_samples` MUST be
    carried and gated on; the biggest raw shifters are measured in 3–17 samples.
    """
    p = OUT / "within_family/seed_shift.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t").copy()
    d["_k"] = d["arm"].where(d["arm"].notna(), d.get("mimat"))
    d = _key_index(d, "_k", src="seed_shift")
    keep = [c for c in ("seed_shift_frac", "seed_shift_sd", "n_samples", "top_shift_nt", "total_rpm")
            if c in d.columns]
    return d[keep].rename(columns={c: f"iso_{c}" for c in keep})


def _isomir_composition() -> pd.DataFrame:
    """⭐ ISOMIR OCCURRENCE per arm — where an arm's reads actually GO once 5′-isomiRs are resolved.

    `seed_shift.tsv` (the `iso_` block) says how much an arm's 5′ end MOVES. This says what that movement
    COSTS: from `seed_composition.tsv.gz` (participant × arm × destination family, with `frac`, `n_reads`,
    `is_canonical`), per arm and averaged over participants:
      * `isoc_canonical_frac` — share of reads still assigned to the arm's OWN seed family.
      * ⭐ `isoc_orphan_frac` — share landing in the **`orphan`** bucket: reads whose shifted seed matches
        NO family. **These are DROPPED from every family**, so `X_fam` does not conserve mass — MH-153
        measured the loss at **4.1% of total RPM (654.3M → 627.4M) with nothing flagging it**, and the
        board's open item is *"audit RECEIVERS and the ORPHAN MASS"*. This column is that flag, per arm.
      * `isoc_n_dest_families` — how many DIFFERENT families the arm donates to (1 = clean).
      * `isoc_donated_frac` — share going anywhere other than its own family or the orphan bucket.

    ⚠ MH-153: high-shift **donors** in multi-arm families are a NON-event (pooling buffers them: miR-17~92
    dose corr 0.95 under `isomir=True`, PTEN's identity argmax unchanged). **The real exposure is the
    mirror image** — small families with no pool to dilute a donation (`miR-543` gains **+351%**), and the
    orphan mass. Read `isoc_orphan_frac` and the RECEIVING side, not the donor's shift alone.
    """
    p = OUT / "seed_composition.tsv.gz"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t")
    need = {"arm", "target_family", "frac"}
    if not need <= set(d.columns):
        return pd.DataFrame().rename_axis("arm")
    d = d.assign(_k=d["arm"].map(lambda v: _arm_key(v, src="seed_composition")))
    d = d[d._k.notna()].copy()
    d["frac"] = pd.to_numeric(d["frac"], errors="coerce")
    is_orph = d["target_family"].astype(str).str.lower().eq("orphan")
    is_canon = d["is_canonical"].astype(bool) if "is_canonical" in d.columns else ~is_orph
    # mean over participants of each per-participant fraction
    out = pd.DataFrame({
        "isoc_canonical_frac": d[is_canon].groupby("_k")["frac"].mean().round(4),
        "isoc_orphan_frac": d[is_orph].groupby("_k")["frac"].mean().round(4),
        "isoc_donated_frac": d[~is_canon & ~is_orph].groupby("_k")["frac"].mean().round(4),
        "isoc_n_dest_families": d[~is_canon & ~is_orph].groupby("_k")["target_family"].nunique(),
    })
    if "n_reads" in d.columns:
        out["isoc_med_reads"] = d.groupby("_k")["n_reads"].median()
    for c in ("isoc_orphan_frac", "isoc_donated_frac", "isoc_n_dest_families"):
        out[c] = out[c].fillna(0)                       # measured-and-zero, not unscanned: the arm IS in
    return out.rename_axis("arm")                       # the table, it simply donates nothing


def _locus_cnv() -> pd.DataFrame:
    """miRNA-locus copy number, cohort stratum (858 arms). ⚠ Multi-locus arms are DILUTED here — the
    paralog-decontaminated view is `foc_` (focal locus)."""
    p = C.REPO_ROOT / "mirna_hallmark/output/mirna_locus_cnv/mirna_cnv_by_stratum.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t")
    d = d[(d.get("entity_level") == "arm") & (d.get("stratification_layer") == "cohort")].copy()
    if d.empty:
        return pd.DataFrame().rename_axis("arm")
    d["_k"] = d["entity_label"].where(d["entity_label"].notna(), d["entity_id"])
    d = _key_index(d, "_k", src="locus_cnv")
    keep = [c for c in ("mean_copy_number", "median_copy_number", "sd_copy_number", "pct_gain",
                        "pct_loss", "pct_amp", "pct_deep_del", "pct_loh", "pct_altered")
            if c in d.columns]
    return d[keep].rename(columns={c: f"cnv_{c.replace('copy_number', 'cn')}" for c in keep})


def _cnv_concordance() -> pd.DataFrame:
    """Does the arm's locus CN actually propagate to its mature expression (756 arms)."""
    p = C.REPO_ROOT / "mirna_hallmark/output/mirna_locus_cnv/mirna_cnv_expr_concordance.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="cnv_concordance")
    keep = [c for c in ("n", "spearman_rho", "partial_rho", "partial_p", "partial_q") if c in d.columns]
    return d[keep].rename(columns={c: f"cnvc_{c}" for c in keep})


def _focal_locus() -> pd.DataFrame:
    """⭐ Paralog-decontaminated CN→expression concordance at the arm's OWN focal locus (733 arms).

    ⛔ This is DESCRIPTIVE ONLY. Do NOT call it "instrument strength", "admissibility" or "F>10" — both
    named CN instruments are RETRACTED (MH-124r/126, MH-133) and the hard F>10 gate was itself replaced
    by a soft null-corrected weight. There is no persisted per-arm F anywhere in the repo.
    """
    p = OUT / "instrument/focal_locus_concordance.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="focal_locus")
    keep = [c for c in ("focal_locus", "n_loci", "multilocus", "n", "partial_rho", "partial_p",
                        "spearman_rho") if c in d.columns]
    return d[keep].rename(columns={c: f"foc_{c}" for c in keep})


def _healthy_anchor() -> pd.DataFrame:
    """GTEx-healthy anchor — MEASURED LEG ONLY.

    ⛔ `gtex_qn_baseline` mixes four legs (floor0 1,137 · mited_anchor 634 · gtex 417 · gtex_family 49).
    `hly_gtex_median`/`hly_gtex_pct` are emitted **only where `source=="gtex"`**; the transferred value
    goes to `hx_hly_baseline_imputed` with `hly_baseline_src`; and **`floor0` is NEVER emitted as a
    value** — that 0-vs-NaN collapse is exactly what MH-210 removed from `budget_shift`.
    """
    p = C.REPO_ROOT / "mirna_hallmark/output/tissue_reference/healthy_anchor/gtex_qn_baseline.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="healthy_anchor")
    src = d.get("source")
    out = pd.DataFrame(index=d.index)
    meas = src.eq("gtex") if src is not None else pd.Series(False, index=d.index)
    out["hly_gtex_median"] = pd.to_numeric(d.get("gtex_median"), errors="coerce").where(meas)
    out["hly_gtex_pct"] = pd.to_numeric(d.get("gtex_pct"), errors="coerce").where(meas)
    out["hly_nat_median"] = pd.to_numeric(d.get("nat_median"), errors="coerce")
    out["hly_nat_detect_n"] = pd.to_numeric(d.get("nat_detect_n"), errors="coerce")
    out["hly_tumor_prev"] = pd.to_numeric(d.get("tumor_prev"), errors="coerce")
    out["hly_baseline_src"] = src
    # ⛔ labelled hx_: this is the miTED rank-transfer / seed-mate value, NOT a GTEx measurement.
    imp = pd.to_numeric(d.get("healthy_baseline_tcga"), errors="coerce")
    out["hx_hly_baseline_imputed"] = imp.where(src.ne("gtex") & src.ne("floor0")) if src is not None else imp
    # ⛔ NOT `bc_` — the broadcast-constancy test (MH-217) caught this on its first run: it
    # varies WITHIN seed_family in 38 of 1,959 groups, and correctly so. It records whether
    # THIS ARM fell back to the seed-mate baseline, which differs between members of one
    # family (one measured in GTEx, the other imputed). Arm-native provenance, not inherited.
    # ⛔ MASKED 2026-08-19 (nan_bool_audit, MH-256): a bare comparison reads False on NaN. `.eq()` on NaN is False too, so the 213 arms with NO recorded
    # baseline source read as "not from a seed-mate" rather than "source unknown".
    out["hly_from_seedmate"] = (src.eq("gtex_family").where(src.notna())
                                if src is not None else np.nan)
    return out


def _survival() -> pd.DataFrame:
    """Per-arm Cox outcome, TCGA legs only (145 arms).

    ⛔ `gse_*` EXCLUDED: its producer matches arms by BASE NAME with no arm suffix — the same bare-stem
    conflation that disqualifies Buffa's 501 probes at arm rung.
    """
    p = C.REPO_ROOT / "mirna_hallmark/output/outcome_survival/arm_survival.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="survival")
    keep = [c for c in d.columns if c.startswith("tcga_")]
    return d[keep].rename(columns={c: f"surv_{c[5:]}" for c in keep})


def _comovement() -> pd.DataFrame:
    """Co-expression module membership (371 arms) — the PURE co-expression columns only.

    ⛔ `dHT`, `primary_class`, `secondary_class`, `trajectory_cluster` EXCLUDED: they come from the
    §6b-retired state-class heuristic, and `dHT` is the exact quantity MH-210 showed the GTEx collapse
    was fabricating.
    """
    p = C.REPO_ROOT / "mirna_hallmark/output/tissue_reference/mirna_comovement/comovement_modules.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "miRNA", src="comovement")
    keep = [c for c in ("comove_module", "top_coexpressed") if c in d.columns]
    return d[keep].rename(columns={"comove_module": "comv_module", "top_coexpressed": "comv_top_partner"})


def _curation_and_fame() -> pd.DataFrame:
    """⛔ FAME AXES (MH-208). Carried to CONTROL for study bias — never as targetome breadth."""
    from mirna_hallmark.analyses.misc import genomewide_promiscuity as GP
    from mirna_hallmark.learned.evidence import ledger as LG
    cur = _read(GP.CACHE, sep="\t").set_index("arm")[["he_degree", "he_degree_expr"]]
    cur.columns = ["cur_he_degree", "cur_he_degree_expr"]
    L = LG.build_ledger()
    fame = pd.DataFrame({"fame_npmid": L.groupby("arm")["pmid"].nunique(),
                         "fame_n_genes_curated": L.groupby("arm")["gene"].nunique()})
    return cur.join(fame, how="outer").rename_axis("arm")


def _ledger_assay_profile() -> pd.DataFrame:
    """⛔ FAME: the arm's curation QUALITY profile, not just its volume.

    Splits LOW-THROUGHPUT FUNCTIONAL (reporter / western / proteomics — what MH-196 defines canonical
    regulators by) from HIGH-THROUGHPUT (ago_clip / qpcr_rna / degradome, 91.6% of ledger rows).
    """
    from mirna_hallmark.learned.evidence import ledger as LG
    L = LG.build_ledger()
    lt = set(getattr(LG, "_LT_FUNC", ("reporter", "western", "proteomics")))
    L = L.assign(_lt=L["assay_class"].isin(lt))
    g = L.groupby("arm")
    out = pd.DataFrame({
        "fame_led_n_pmid": g["pmid"].nunique(),
        "fame_led_n_ltfunc_pmid": L[L._lt].groupby("arm")["pmid"].nunique(),
        "fame_led_n_classes": g["assay_class"].nunique(),
    })
    out["fame_led_n_ltfunc_pmid"] = out["fame_led_n_ltfunc_pmid"].fillna(0)
    out["fame_led_frac_ltfunc"] = out["fame_led_n_ltfunc_pmid"] / out["fame_led_n_pmid"].replace(0, np.nan)
    if "weak" in L.columns:
        out["fame_led_n_weak"] = L[L["weak"].astype(bool)].groupby("arm").size()
    # ⛔ PRUNED 2026-08-19: `fame_led_n_pmid` was BIT-IDENTICAL to `fame_npmid` from `_curation_and_fame`
    # — the same `groupby("arm")["pmid"].nunique()` on the same ledger, computed twice. It is still needed
    # as the denominator of `fame_led_frac_ltfunc` above, so it is computed and then dropped.
    out = out.drop(columns=["fame_led_n_pmid"])
    return out.rename_axis("arm")


def _curation_profile() -> pd.DataFrame:
    """⛔ FAME: per-assay-class miRTarBase study counts (widest arm coverage of any evidence source)."""
    p = C.REPO_ROOT / "mirna_hallmark/output/edges/mirtarbase/mirtar_mirna_summary.csv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p)
    kcol = "miRNA" if "miRNA" in d.columns else d.columns[0]
    d = _key_index(d, kcol, src="mirtar_summary")
    keep = [c for c in d.columns if c.startswith("n_") and "json" not in c.lower()]
    # ⛔ PRUNED 2026-08-19: this block took EVERY `n_*` column, and five of them are ALL-ZERO and
    # bit-identical — miRTarBase carries no `perturbation` assay class for these arms, nor any
    # `binding__nonfunctional_mti_weak` row. They were 5 dead columns on the largest card (the `fame_`
    # block is 43 of the arm card's 297) for an axis the doctrine calls a confound-to-control, never
    # biology. Dropping CONSTANT columns is measured, not curated, so a class that gains rows later
    # reappears automatically.
    d = d[keep]
    # ⛔ THREE MEASURED DROP RULES (column review unit 4, 2026-08-19; user chose the conservative cut).
    # All three are MEASURED, never a hardcoded name list — so a class that gains signal later reappears
    # on its own. The block is 31 columns with ZERO consumers in any module (verified), and the user
    # elected to prune only what is provably dead rather than the whole cross-tabulation.
    #   (1) CONSTANT            — carries no information at all.
    #   (2) NEAR-DEGENERATE     — >95% zero or <=4 distinct values. Eight qualify, e.g.
    #       `binding__nonfunctional_mti` (99.9% zero, 2 values) and `protein__nonfunctional_mti_weak`
    #       (99.8%, 2 values). A count that is zero for 99% of arms cannot separate anything.
    # ⛔⛔ (3) THE CORRELATION RULE IS REMOVED — it was WRONG, and checking the ingestion is what showed it.
    #     I dropped two columns at r>=0.9999 (`unique_targets` ~ `total_studies`, `binding_studies` ~
    #     `binding__functional_mti_weak_studies`) calling them redundant. **Both were DATA COINCIDENCES,
    #     not definitional identities**, and `pipeline/genes/mirtarbase.py` settles it:
    #         n_unique_targets  = group["gene"].nunique()          <- DISTINCT GENES
    #         n_total_studies   = group["study_id"].nunique()      <- DISTINCT STUDIES
    #     They differ on **293 of 3,039 arms**; they track because most (miRNA, gene) pairs carry exactly
    #     ONE study. Likewise `binding_studies` vs its functional-weak cell differ on 12 arms.
    #     ⇒ `unique_targets` was NOT a misnomer — it is the only true TARGET count in this block, and I
    #     removed it on a rule that cannot tell "identical by definition" from "identical in this pull".
    #     **A correlation rule prunes columns that will diverge the moment curation changes.** Only EXACT
    #     equality is safe, and the duplicate audit already covers that.
    dropped: dict[str, str] = {}
    for c in list(d.columns):
        s = pd.to_numeric(d[c], errors="coerce").dropna()
        if not len(s) or s.nunique() <= 1:
            dropped[c] = "constant"
        elif (s == 0).mean() > 0.95 or s.nunique() <= 4:
            dropped[c] = f"degenerate ({(s == 0).mean():.0%} zero, {s.nunique()} distinct)"
    if dropped:
        print(f"  fame_assay: dropping {len(dropped)} of {len(d.columns)} column(s)")
        for c, why in sorted(dropped.items()):
            print(f"      {c[2:]:<46} {why}")
        d = d.drop(columns=list(dropped))
    return d.rename(columns={c: f"fame_assay_{c[2:]}" for c in d.columns})


def _model_footprint() -> pd.DataFrame:
    """How much of the LEARNED MODEL DESIGN this arm touches — the design, not the literature."""
    from mirna_hallmark.learned.evidence import ledger as LG
    he = LG.pooled_he_edges()
    g = he.groupby("miRNA")
    # ⛔ PRUNED 2026-08-19: `model_n_edges` was bit-identical to `model_n_genes` — same structural reason
    # as `arb_n_genes`: one row per (arm, gene), so `size()` == `gene.nunique()` by construction.
    return pd.DataFrame({"model_n_genes": g["gene"].nunique()}).rename_axis("arm")


def _beta_rollup() -> pd.DataFrame:
    """⭐ Per-arm rollup of the LEARNED β (733 arms) — did not exist anywhere before.

    ⛔ `identity` spans −2,791…+2,792 on this table — a ratio evaluated where its denominator vanishes
    (axiom 5). `arb_max_identity` is computed ONLY over `beta_frac_reliable` rows and ships its own
    gated denominator `arb_n_identity_reliable`.
    """
    p = OUT / "readouts_arm_edges.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t")
    kcol = "arm" if "arm" in d.columns else "miRNA"
    d = d.assign(_k=d[kcol].map(lambda v: _arm_key(v, src="readouts_arm_edges")))
    d = d[d._k.notna()]
    g = d.groupby("_k")
    # ⛔ PRUNED 2026-08-19: `arb_n_genes` was BIT-IDENTICAL to `arb_n_edges`. Not a bug in the arithmetic —
    # `readouts_arm_edges.tsv` holds exactly one row per (arm, gene) (verified: 0 duplicate pairs), so
    # `size()` and `gene.nunique()` are the same number BY CONSTRUCTION. The column promised a second
    # quantity and delivered the first.
    out = pd.DataFrame({
        "arb_n_edges": g.size(),
        "arb_mean_abs_beta": g["beta"].apply(lambda s: s.abs().mean()) if "beta" in d.columns else np.nan,
        "arb_max_abs_beta": g["beta"].apply(lambda s: s.abs().max()) if "beta" in d.columns else np.nan,
    })
    if "identified" in d.columns:
        out["arb_n_identified"] = g["identified"].sum()
        out["arb_frac_identified"] = g["identified"].mean()
    if "composition_class" in d.columns:
        out["arb_n_comp_explained"] = d[d.composition_class.eq("composition_explained")].groupby("_k").size()
        out["arb_n_comp_explained"] = out["arb_n_comp_explained"].fillna(0)
    if "identity" in d.columns:
        # ⛔⛔ FIXED 2026-08-19 — THIS GATE WAS INERT AND THE DOCSTRING ABOVE RELIED ON IT.
        # It used `beta_frac_reliable`, which watches **β** sign-cancellation. MH-124 fixed the
        # `_rtnorm_pos` sampler bug that produced negative βs, so β is now strictly positive
        # (0 negatives / 5,802) ⇒ `net_pressure ≡ 1.0` ⇒ the gate admitted **100.0%** of rows and
        # excluded **0 of the 124** rows with |identity| > 1. `arb_max_identity` therefore reached
        # **+740.0** while its docstring claimed it was protected. The live cancellation is in
        # **identity** (9.9% negative — legitimate suppressor contributions), so gate on THAT.
        # `readouts.add_reliability` now emits `identity_reliable`; recomputed here so the card is
        # correct against artifacts written before that change.
        if "identity_reliable" in d.columns:
            keep = d["identity_reliable"].astype(bool)
        else:
            ia = d.groupby("gene")["identity"].transform(lambda s: s.abs().sum())
            keep = ((1.0 / ia.replace(0, np.nan)) >= 0.5) & (d["identity"].abs() <= 1.0)
            keep = keep.fillna(False)
        rel = d[keep]
        out["arb_max_identity"] = rel.groupby("_k")["identity"].max()
        out["arb_n_identity_reliable"] = rel.groupby("_k").size().reindex(out.index).fillna(0)
    # ⭐ THE β LADDER — the second reading of "how many genes does this arm target with such-and-such
    # coupling". `real_*` counts REALIZED (observed bulk) coupling; this counts the MODEL's coupling.
    # They are different estimands and must not be read as one: MH-166 — coupling is REALIZATION, β is
    # the de-confounded joint estimate. An arm can have high β on a gene it does not visibly realize.
    # ⚠ Cuts are set from the measured distribution, not by taste: β is half-normal-slabbed so it is
    # strictly ≥0 (min 0.0003, max 1.52, 0 negatives / 5,646), and it is HEAVILY right-skewed —
    # median 0.005, p75 0.016, p90 0.044, p95 0.091, p99 0.719. A cut at 0.10 therefore selects roughly
    # the top 5% of edges; anything below ~0.01 is the bulk and separates nothing.
    # ⛔⛔ AND THE UNGATED β LADDER IS MISLEADING ON ITS OWN — caught by reading its own output. The arms
    # that TOP it (miR-302a/b/d-3p, miR-603, miR-663a: 9-12 edges all above β=0.10) have
    # `arb_n_identified == 0` and no realized coupling at all. Mechanism, and it is structural, not a
    # bug: the half-normal slab has a STRICTLY POSITIVE MEAN, so an un-informed family **cannot be
    # zeroed** — measured genome-wide at `beta == 0` in 0/5,117 (MH-138). A large β with no
    # identification is the prior showing through, not a strong edge. ⇒ the `_ident` ladder is the one
    # to read; the raw one is kept only so the gap between them is visible.
    if "beta" in d.columns:
        ident = d[d["identified"].astype(bool)] if "identified" in d.columns else None
        for cut in (0.01, 0.05, 0.10, 0.25):
            tag = f"b{int(cut * 100):03d}"
            out[f"arb_n_{tag}"] = (d[d["beta"] > cut].groupby("_k").size()
                                   .reindex(out.index).fillna(0))
            if ident is not None:
                out[f"arb_n_{tag}_ident"] = (ident[ident["beta"] > cut].groupby("_k").size()
                                             .reindex(out.index).fillna(0))
    return out.rename_axis("arm")


def _arm_identifiability() -> pd.DataFrame:
    """Arm-rung identifiability rollup. DOMAIN = multi-arm (gene,family) cells only (20.3% of edges).

    ⚠ Folds the IDENTIFIABILITY VERDICT, not `beta_arm` — MH-191 showed the card's `beta` is already
    fit per arm, so `arm_rung`'s β partly re-derives it. What is unique here is `arm_resolvable`/OOF.
    """
    p = OUT / "arm_rung.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t")
    d = d.assign(_k=d["arm"].map(lambda v: _arm_key(v, src="arm_rung")))
    d = d[d._k.notna()]
    g = d.groupby("_k")
    out = pd.DataFrame({"aid_n_cells": g.size()})
    for c, name in (("arm_resolvable", "aid_frac_resolvable"), ("arm_sep_z", "aid_med_arm_sep_z"),
                    ("oof_drho", "aid_med_oof_drho")):
        if c in d.columns:
            out[name] = g[c].mean() if c == "arm_resolvable" else g[c].median()
    return out.rename_axis("arm")


def _chimeric() -> pd.DataFrame:
    """AGO-ligation (chimeric) evidence PER SOURCE — never pooled.

    ⚠⚠ Four sources with INCOMPATIBLE scales (Manakov weight = reads, max ~14k; TarBase = record count,
    max ~90) and ~52% pair overlap between Manakov and TarBase qCLASH ⇒ a naive sum double-counts and
    lets Manakov dominate. `chim_n_sources` is the honest cross-source column (corroboration).
    ⚠⚠ **Manakov is ~90% of all rows and is 100% HEK293T — a cell line, NOT breast.** No breast chimeric
    data exists in any source. Do not read `chim_` as breast targeting evidence.
    """
    p = OUT / "chimeric_evidence.parquet"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p)
    d = d.assign(_k=d["mature"].map(lambda v: _arm_key(v, src="chimeric")))
    d = d[d._k.notna()].copy()
    d["_grp"] = np.where(d["source"].astype(str).str.startswith("manakov"), "manakov", "tarbase")
    out = pd.DataFrame(index=pd.Index(sorted(d._k.unique()), name="arm"))
    for grp in ("manakov", "tarbase"):
        s = d[d._grp == grp]
        out[f"chim_{grp}_n_genes"] = s.groupby("_k")["gene"].nunique()
        out[f"chim_{grp}_weight"] = s.groupby("_k")["weight"].sum()
    out["chim_n_sources"] = d.groupby("_k")["source"].nunique()
    return out


def _ago() -> pd.DataFrame:
    """RISC loading. ⛔ `ago_loading.loading()` is NOT used — it ends `.fillna(1.0)`, fabricating a guide
    call for every unmeasured arm (verified `ago_loading.py:154`). `ago_dom` preserves NaN, and
    `cov_ago_dom` says who was actually measured. A 1.0 here means MEASURED dominance, not "unknown".
    """
    blocks = []
    try:
        from mirna_hallmark.learned import ago_loading as AL
        eng = AL._manakov_engagement()
        dom = AL._dominance(eng)
        blocks.append(pd.DataFrame({"ago_dom": dom, "ago_dom_src": "manakov"}))
    except Exception as e:
        print(f"  ⚠ AGO dominance unavailable ({e})")
    try:
        ago = _read(OUT / "ago_loading.parquet")
        blocks.append(ago.rename_axis("arm"))
    except Exception as e:
        print(f"  ⚠ ago_loading.parquet unavailable ({e})")
    if not blocks:
        return pd.DataFrame().rename_axis("arm")
    out = blocks[0]
    for b in blocks[1:]:
        out = out.join(b, how="outer")
    # ⭐ THE CATEGORICAL GUIDE/PASSENGER CALL — the operational rule (`dominance < 0.2 ⇒ passenger`) that
    # `ago_loading._report` prints but never materialises anywhere in the subproject.
    # ⚠ It is DEFINED ONLY where dominance was MEASURED (678 arms); everywhere else it is NaN, not
    # "guide". That is the whole point of not using `loading()`, which fabricates 1.0 for the unmeasured.
    # ⛔ The miRBase -5p/-3p label is NOT the designation source — `ago_loading`'s docstring rejects it,
    # and miR-181a-5p is the guide on abundance (10.4 vs 7.8 RPM) regardless of its suffix.
    if "ago_dom" in out.columns:
        v = pd.to_numeric(out["ago_dom"], errors="coerce")
        out["ago_guide_class"] = pd.Series(
            np.where(v.isna(), None, np.where(v < 0.2, "passenger",
                                       np.where(v > 0.8, "guide", "co_loaded"))), index=out.index)
    return out.rename_axis("arm")          # index normalised centrally by `_norm_index` at join time


def _methylation() -> pd.DataFrame:
    """⚠ BROADCAST from the pri-miRNA HAIRPIN — 5p and 3p share one promoter, so they share Δβ.

    ⛔ Does NOT call `promoter_probes(all_arms)`: that ends `.drop_duplicates("probeID")` ACROSS arms,
    silently dropping the second arm's copy of a shared probe. The arm→probe map is built directly.
    ⚠ Only the GROUP tumour-vs-normal contrast is available; the per-sample lane is BLOCKED on the
    UUID→TCGA-barcode map.
    """
    try:
        from mirna_hallmark.learned import methylation as ME
        hp, pr, mat = ME._hairpin_loci(), ME._probe_ref(), ME._beta_cache()
        if mat is None or hp.empty or pr.empty:
            return pd.DataFrame().rename_axis("arm")
        tcols = [c for c in mat.columns if str(c).startswith("T::")]
        ncols = [c for c in mat.columns if str(c).startswith("N::")]
        rows = {}
        for r in hp.drop_duplicates("pre_gene_id").itertuples():
            lo = r.tss - 2000 if r.strand == "+" else r.tss - 500
            hi = r.tss + 500 if r.strand == "+" else r.tss + 2000
            hit = pr[(pr["chrom"] == r.chrom) & (pr["center"] >= lo) & (pr["center"] <= hi)]
            if hit.empty:
                continue
            sub = mat.reindex(hit["probeID"].tolist())
            t, n = float(sub[tcols].mean(axis=0).mean()), float(sub[ncols].mean(axis=0).mean())
            # ⚠ NO PARTIAL ROWS. A locus whose probes carry no β data would otherwise ship
            # `n_probes` with a NaN Δβ — a value outside its own coverage flag, i.e. exactly the
            # "missing = unscanned" contract this card exists to keep. Emit all-or-nothing.
            if not (np.isfinite(t) and np.isfinite(n)):
                continue
            rows[r.arm_key] = {"bc_meth_n_probes": len(hit),
                               "bc_meth_n_cgi": int(hit["in_CGI"].sum()) if "in_CGI" in hit else np.nan,
                               "bc_meth_tumour_beta": round(t, 4), "bc_meth_normal_beta": round(n, 4),
                               "bc_meth_dbeta": round(t - n, 4)}
        if not rows:
            return pd.DataFrame().rename_axis("arm")
        out = pd.DataFrame.from_dict(rows, orient="index").rename_axis("arm_key")
        out["bc_meth_direction"] = np.where(out["bc_meth_dbeta"] >= 0.15, "hyper",
                                     np.where(out["bc_meth_dbeta"] <= -0.15, "hypo", "none"))
        return out
    except Exception as e:
        print(f"  ⚠ methylation unavailable ({e})")
        return pd.DataFrame().rename_axis("arm")


def _family_context(card_index, fam_series) -> pd.DataFrame:
    """⚠ BROADCAST from the seed FAMILY (`genomic_context.family_context()`)."""
    p = OUT / "family_context.tsv"
    if not p.exists() or fam_series is None:
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t")
    kcol = "family" if "family" in d.columns else d.columns[0]
    d = d.set_index(kcol)
    keep = [c for c in ("n_members", "context_homogeneous", "dose_coding_host_frac",
                        "dose_cotranscribed_frac", "n_distinct_coding_hosts") if c in d.columns]
    if not keep:
        return pd.DataFrame().rename_axis("arm")
    m = d[keep].rename(columns={c: f"bc_fam_{c}" for c in keep})
    return m.reindex(fam_series.to_numpy()).set_index(pd.Index(card_index, name="arm"))


# --------------------------------------------------------------------------- #
# v3 — ATTRIBUTION, REALIZATION and WITHIN-FAMILY ROLE (MH-215)
# --------------------------------------------------------------------------- #
def _attribution() -> pd.DataFrame:
    """⭐ THE ARM SHAPLEY / ATTRIBUTION ROLL-UP — the axis v2 was missing entirely.

    v2's only Shapley-derived column was `arb_max_identity`. The dedicated per-arm decomposition
    (`task3c_shapley_arm.tsv`, joint/within/nested/budget share) lives at (gene, arm) grain, so it needs
    a roll-up — which is why it fell outside v2's Tier 1, not because it was unwanted.

    ⚠ Every column here is an AGGREGATE OVER THE ARM'S EDGES (declared `agg_of="edge"`): "how this arm
    is credited across the genes it regulates", never a single-gene attribution.
    ⚠ `attr_*_share` are gated on the source's own `reliable` flag where it exists — a share is a ratio,
    and MH-119/MH-146 are what happens when you average one over rows whose denominator vanishes.
    """
    out = []
    p = OUT / "task3c_shapley_arm.tsv"
    if p.exists():
        d = _read(p, sep="\t")
        d = d.assign(_k=d["arm"].map(lambda v: _arm_key(v, src="shapley")))
        d = d[d._k.notna()]
        rel = d[d["reliable"].astype(bool)] if "reliable" in d.columns else d
        g, gr = d.groupby("_k"), rel.groupby("_k")
        b = pd.DataFrame({"attr_n_edges": g.size(), "attr_n_reliable": gr.size()})
        for c, nm in (("joint_share", "joint"), ("within_share", "within"),
                      ("nested_share", "nested"), ("budget_share", "budget")):
            if c in d.columns:
                b[f"attr_med_{nm}_share"] = gr[c].median()
                b[f"attr_max_{nm}_share"] = gr[c].max()
        if "phi_joint_raw" in d.columns:
            b["attr_med_phi_joint"] = gr["phi_joint_raw"].median()
        out.append(b)
    p = OUT / "rung_parity.tsv"
    if p.exists():
        d = _read(p, sep="\t")
        d = d.assign(_k=d["arm"].map(lambda v: _arm_key(v, src="rung_parity")))
        d = d[d._k.notna()].groupby("_k")
        out.append(pd.DataFrame({"attr_med_credit_share": d["arm_credit_share"].median(),
                                 "attr_max_credit_share": d["arm_credit_share"].max()}))
    p = OUT / "arm_attribution_from_family.tsv"
    if p.exists():
        d = _read(p, sep="\t")
        d = d.assign(_k=d["miRNA"].map(lambda v: _arm_key(v, src="arm_attribution")))
        d = d[d._k.notna()].groupby("_k")
        out.append(pd.DataFrame({"attr_med_dose_share": d["dose_share"].median()}))
    if not out:
        return pd.DataFrame().rename_axis("arm")
    o = out[0]
    for b in out[1:]:
        o = o.join(b, how="outer")
    return o.rename_axis("arm")


def _subtype() -> pd.DataFrame:
    """Which PAM50 subtype this arm's edges couple in (roll-up over its edges).

    ⚠ MH-165: subtype coupling heterogeneity is REAL but **distributional** — a shift in the spread of
    per-edge coupling, NOT a validated per-edge subtype assignment. `sub_best_mode` is descriptive.
    """
    p = OUT / "subtype_edge_heterogeneity.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t")
    d = d.assign(_k=d["arm"].map(lambda v: _arm_key(v, src="subtype")))
    d = d[d._k.notna()]
    g = d.groupby("_k")
    out = pd.DataFrame({"sub_n_edges": g.size(), "sub_med_rho_pool": g["rho_pool"].median(),
                        "sub_med_rho_best": g["rho_best"].median()})
    if "best_sub" in d.columns:
        out["sub_best_mode"] = g["best_sub"].agg(lambda s: s.mode().iat[0] if len(s.mode()) else np.nan)
        out["sub_best_frac"] = g["best_sub"].agg(lambda s: s.value_counts(normalize=True).max()
                                                 if s.notna().any() else np.nan)
    if "q_adj" in d.columns:
        out["sub_n_q05"] = d[d.q_adj < 0.05].groupby("_k").size()
        out["sub_n_q05"] = out["sub_n_q05"].fillna(0)
    return out.rename_axis("arm")


def _field_effect() -> pd.DataFrame:
    """How much of the arm's TUMOUR level is its own patient's NAT baseline (field effect) vs acquired."""
    p = OUT / "realization/mirna_nat_retention.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="nat_retention")
    keep = [c for c in ("n", "r_own", "r_perm", "retention", "r_adj") if c in d.columns]
    return d[keep].rename(columns={c: f"field_{c}" for c in keep})


def _compartment() -> pd.DataFrame:
    """Which cell compartment the arm's DOSE tracks, and how tightly (`arm_lock`).

    ⚠⚠ PROVENANCE: `arm_compartment_loading.tsv` has **NO PRODUCER in the repo** (MH-111 ad-hoc) — the
    same shape that rotted the MH-196 literature sets. Carried because the axis is real and MH-201 showed
    composition axes are the ones that predict SIGN, but flagged: re-derive before it carries a claim.
    """
    p = OUT / "arm_compartment_loading.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _key_index(_read(p, sep="\t"), "arm", src="compartment")
    keep = [c for c in ("arm_lock", "arm_top_compartment", "dose_rpm") if c in d.columns]
    # ⛔ prefix is `cload_`, NOT `comp_`: the GENE and FAMILY cards already carry a `comp_tcga_mrna_*`
    # block (the CPTAC/TCGA compartment analysis), and `card_rungs.domain_of` matches by prefix
    # GLOBALLY across cards — a bare `comp_` here mislabelled 10 of their columns with THIS
    # block's domain string. Third instance of that root cause (cf. `esub_`, `echim_`).
    return d[keep].rename(columns={c: f"cload_{c.replace('arm_', '')}" for c in keep})


# ⭐ THE REALIZATION LADDER. axiom 5: a threshold count is only honest as a SWEEP — one cut is a
# coin-flip with decimals, so every cut is shipped and the reader sees where the mass sits.
_COUPLING_CUTS = (-0.05, -0.10, -0.15, -0.20, -0.30)


def _realization_ladder() -> pd.DataFrame:
    """⭐ HOW MUCH OF THIS ARM'S TARGETOME IS ACTUALLY REALIZED — a threshold LADDER, never one cut.

    Counts, per arm, how many of its scored target genes reach each realized-coupling depth
    (`coupling_tum`, the per-state site-free-null-CALIBRATED coupling — MH-166, not the retired −0.1 cut),
    plus how many clear calibrated significance.

    ⚠⚠ **THE RATE COLUMNS ARE THE DANGEROUS ONES AND ARE GATED.** `real_frac_*` divides by the number of
    SCORED edges (curated design), and `real_rate_vs_seq_*` divides by the arm's genome-wide targetome —
    two completely different denominators, and the second spans 85–7,666. Both are emitted ONLY where the
    denominator is non-trivial, and the raw numerator + denominator always ship beside them so the ratio
    can be re-derived. ⛔ Never quote a rate without `real_n_scored` / `seq_n_genes_strong`.
    ⚠ `coupling_tum` is a BULK observational correlation: presence = dominant realized regulator,
    ABSENCE ≠ not-real (MH-166: coupling is REALIZATION, not existence, and is ~13% sensitive).
    """
    p = OUT / "realization/edge_card.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t", low_memory=False)
    if "coupling_tum" not in d.columns:
        return pd.DataFrame().rename_axis("arm")
    d = d.assign(_k=d["arm"].map(lambda v: _arm_key(v, src="realization_ladder")))
    d = d[d._k.notna() & d["coupling_tum"].notna()]
    g = d.groupby("_k")
    out = pd.DataFrame({"real_n_scored": g.size(),
                        "real_med_coupling": g["coupling_tum"].median(),
                        "real_best_coupling": g["coupling_tum"].min()})
    for cut in _COUPLING_CUTS:
        tag = f"c{abs(int(cut * 100)):02d}"
        out[f"real_n_{tag}"] = d[d.coupling_tum < cut].groupby("_k").size().reindex(out.index).fillna(0)
    if "coupling_p_tum" in d.columns:
        out["real_n_sig"] = d[d.coupling_p_tum < 0.05].groupby("_k").size().reindex(out.index).fillna(0)
    # gated rates — denominator must be non-trivial, and the numerator ships beside it
    ok = out["real_n_scored"] >= 5
    out["real_frac_c10"] = (out["real_n_c10"] / out["real_n_scored"]).where(ok)
    if "real_n_sig" in out.columns:
        out["real_frac_sig"] = (out["real_n_sig"] / out["real_n_scored"]).where(ok)
    return out.rename_axis("arm")


def _family_role(card: pd.DataFrame) -> pd.DataFrame:
    """⭐ THE ARM'S ROLE INSIDE ITS SEED FAMILY — family structure dissected, then annotated per arm.

    A seed family is one targeting unit whose members are interchangeable by SEQUENCE but wildly
    different by DOSE: measured on this card, members of a multi-arm family sit a median ~2.4 log2 apart.
    So "which member IS the family" is an arm-level fact the family rung cannot express, and it is exactly
    what `bc_fam_*` (a broadcast, constant within family) cannot say.

    Computed here from the card's own `abund_*` so it needs no new source:
      * `famrole_abund_share` — this arm's share of its family's TOTAL LINEAR RPM (a true pool,
        `2^x − 1`, never a mean of logs — the pooling defect `decoy_bench.pool_family` documents).
      * `famrole_abund_rank` / `famrole_var_rank` — 1 = the family's loudest / most variable member.
      * `famrole_is_dominant` — this arm carries ≥50% of the family's dose.
      * `famrole_class` — sole / dominant / co / minor.
    ⚠ `famrole_*` are DEGENERATE for single-arm families (share ≡ 1, rank ≡ 1 by construction) —
    `famrole_n_members` is shipped so that stratum can be split off, exactly as `gene_axes.mask_degenerate`
    demands for the one-unit case.
    """
    if "seed_family" not in card.columns or "abund_med" not in card.columns:
        return pd.DataFrame().rename_axis("arm")
    d = card[["seed_family", "abund_med", "abund_sd"]].copy()
    d = d[d["seed_family"].notna()]
    lin = np.power(2.0, pd.to_numeric(d["abund_med"], errors="coerce")) - 1.0
    d["_lin"] = lin.clip(lower=0)
    g = d.groupby("seed_family")
    tot = g["_lin"].transform("sum")
    out = pd.DataFrame(index=d.index)
    out["famrole_n_members"] = g["_lin"].transform("size")
    out["famrole_abund_share"] = (d["_lin"] / tot.replace(0, np.nan)).round(4)
    out["famrole_abund_rank"] = g["_lin"].rank(ascending=False, method="min")
    out["famrole_var_rank"] = g["abund_sd"].rank(ascending=False, method="min")
    # ⛔ MASKED 2026-08-19 (nan_bool_audit, MH-256): a bare comparison reads False on NaN. 219 arms with no abundance share read "not dominant".
    out["famrole_is_dominant"] = (out["famrole_abund_share"] >= 0.5).where(
        out["famrole_abund_share"].notna())
    out["famrole_class"] = np.where(out["famrole_n_members"] == 1, "sole",
                             np.where(out["famrole_abund_share"] >= 0.5, "dominant",
                               np.where(out["famrole_abund_share"] >= 0.2, "co", "minor")))
    return out.rename_axis("arm")


def _cptac_cov() -> pd.Series:
    try:
        from mirna_hallmark.learned import cptac_data as CD
        idx = [_arm_key(a, src="cptac") for a in CD.arms().index]
        return pd.Series(True, index=pd.Index(sorted({a for a in idx if a}), name="arm"), name="cov_cptac")
    except Exception as e:
        print(f"  ⚠ CPTAC coverage unavailable ({e})")
        return pd.Series(dtype=bool, name="cov_cptac")


# --------------------------------------------------------------------------- #
# EDGE-CARD LIFT
# --------------------------------------------------------------------------- #
_LIFT = ("arm_med_rpm", "arm_pct_floor", "arm_iqr", "spiker", "detection", "healthy_leg",
         "healthy_potential", "healthy_uninformative", "surrogate_instrument", "surrogate_corr",
         "dose_comp_retention", "dose_prolif_retention", "dose_confounded",
         "arm_lfc_NAT_TUM", "arm_lfc_HLY_TUM_QN", "arm_lfc_HLY_TUM_raw", "arm_lfc_HLY_NAT_raw",
         "grank_HLY", "grank_NAT", "grank_TUM",
         "dGlobal_HLY_NAT", "dGlobal_NAT_TUM", "dGlobal_HLY_TUM",
         "ctx_arm_dose", "ctx_arm_abundant")


def _lift_edge_card(verify: bool = True) -> pd.DataFrame:
    """Lift the ARM-rung columns already sitting on the edge card (dedupe by arm, 736 arms).

    ⛔ `arm_id_status` is NOT lifted — measured to vary within arm in 263/736 arms, so it is
    arm-in-family rung despite sitting in `card_rungs._ARM_DESC`.
    ⛔ `ago_loading` is NOT lifted — it carries the fabricated 1.0 (see `_ago`).
    ⛔ `mean_own_shift`/`own_specific_frac`/`mean_dGlobalRank` are NOT lifted — `_within_patient_shift`
    supersedes them at 2,236 arms vs 736 (verified identical on the overlap).
    """
    p = OUT / "realization/edge_card.tsv"
    if not p.exists():
        return pd.DataFrame().rename_axis("arm")
    d = _read(p, sep="\t", low_memory=False)
    cols = [c for c in _LIFT if c in d.columns]
    if verify:                        # a test that passes everything is not a test
        bad = {c: int((d.groupby("arm")[c].nunique(dropna=True) > 1).sum()) for c in cols}
        viol = {c: n for c, n in bad.items() if n}
        if viol:
            print(f"  ⚠ NOT arm-invariant, dropped from the lift: {viol}")
            cols = [c for c in cols if not bad[c]]
        if "arm_id_status" in d.columns:
            n = int((d.groupby("arm")["arm_id_status"].nunique(dropna=True) > 1).sum())
            print(f"  [lift] control — arm_id_status varies in {n} arms (expected ~263; a test that "
                  f"flags nothing is not a test)")
    out = d.drop_duplicates("arm")[["arm"] + cols]
    return _key_index(out, "arm", src="edge_card")


# --------------------------------------------------------------------------- #
def build() -> pd.DataFrame:
    print("[arm_card] building blocks…")
    blocks = {
        "abundance": _abundance(), "expression tier": _expression_tier(),
        "sequence targetome": _sequence_targetome(), "site composition": _site_composition(),
        "targetscan": _targetscan(), "seed scan (4th universe)": _seed_scan(),
        "genomic context": _genomic_context(), "within-patient shift": _within_patient_shift(),
        "isomiR seed shift": _isomir(), "isomiR composition": _isomir_composition(), "locus CNV": _locus_cnv(),
        "CN concordance": _cnv_concordance(), "focal locus": _focal_locus(),
        "healthy anchor": _healthy_anchor(), "survival": _survival(), "comovement": _comovement(),
        "curation/fame": _curation_and_fame(), "ledger assay profile": _ledger_assay_profile(),
        "mirtar curation profile": _curation_profile(), "model footprint": _model_footprint(),
        "beta rollup": _beta_rollup(), "arm identifiability": _arm_identifiability(),
        "chimeric": _chimeric(), "AGO": _ago(), "methylation (bc)": _methylation(),
        "edge-card lift": _lift_edge_card(),
        # ── v3 (MH-215): attribution, realization and within-family role
        "attribution/shapley": _attribution(), "subtype": _subtype(),
        "field effect": _field_effect(), "compartment": _compartment(),
        "realization ladder": _realization_ladder(),
    }
    card = None
    for name, b in blocks.items():
        if b is None or not len(b):
            print(f"  ⚠ {name:26s}: EMPTY, skipped")
            continue
        b = _norm_index(b, src=name)
        print(f"  {name:26s}: {len(b):5,d} arms × {b.shape[1]:2d} cols")
        card = b if card is None else card.join(b, how="outer")

    # ⚠ HARD UNIVERSE GATE — v1 shipped 506 locus-ids and 5 unicode rows as if they were arms.
    # Every block index is already normalised by `_norm_index`; this is the belt-and-braces assertion.
    print(f"\n[arm_card] {len(card):,} arms after normalised join")
    assert not any(str(a).startswith("MI0") for a in card.index), "MI0* locus ids leaked into the card"
    assert not any(any(ord(ch) > 127 for ch in str(a)) for a in card.index), "non-ASCII key leaked"

    # seed family + the family broadcast
    try:
        from mirna_hallmark.coupling_inference import seed_family_map
        fam = pd.Series(seed_family_map(list(card.index)), name="seed_family").reindex(card.index)
        card["seed_family"] = fam
        # ⛔ PRUNED 2026-08-19: `family_n_arms_card` was bit-identical to `famrole_n_members`, which is
        # computed below and belongs to a coherent block (the arm's role inside its family). One name
        # kept, and it is the one that sits with its siblings. Nothing read this column
        # (grep: written here, read nowhere).
        # ⚠ Neither is `family_size` — the EDGE card owns that name at FAMILY rung over the MODEL DESIGN.
        fc = _family_context(card.index, fam)
        if len(fc):
            card = card.join(fc)
        # ⭐ v3: the arm's ROLE inside its family — must run AFTER seed_family exists, and it reads the
        # card's own abundance rather than a new source.
        fr = _family_role(card)
        if len(fr):
            card = card.join(fr)
            print(f"  {'family role (v3)':26s}: {len(fr):5,d} arms × {fr.shape[1]:2d} cols")
    except Exception as e:
        print(f"  ⚠ seed family map unavailable ({e})")

    cc = _cptac_cov()
    card["cov_cptac"] = card.index.isin(cc.index) if len(cc) else False

    # ⚠ COVERAGE FLAGS — a missing block means UNSCANNED, never zero. Nothing is imputed.
    for flag, probe in (("cov_seq_scan", "seq_n_genes_strong"), ("cov_site_types", "site_n_total"),
                        ("cov_targetscan", "ts_n_sites"), ("cov_seed_scan", "sdsz_N_eff"),
                        ("cov_curated", "cur_he_degree"), ("cov_gctx", "gctx_mir_class"),
                        ("cov_cnv", "cnv_mean_cn"), ("cov_focal", "foc_partial_rho"),
                        ("cov_isomir", "iso_seed_shift_frac"), ("cov_isocomp", "isoc_canonical_frac"), ("cov_beta", "arb_n_edges"),
                        ("cov_ago_dom", "ago_dom"), ("cov_meth", "bc_meth_dbeta"),
                        ("cov_expr", "abund_med"),
                        # v3
                        ("cov_attr", "attr_n_edges"), ("cov_subtype", "sub_n_edges"),
                        ("cov_field", "field_r_own"), ("cov_comp", "cload_lock"),
                        ("cov_real", "real_n_scored"), ("cov_famrole", "famrole_abund_share")):
        card[flag] = card[probe].notna() if probe in card.columns else False

    card = card.rename_axis("arm").reset_index()
    DEST.parent.mkdir(parents=True, exist_ok=True)
    card.to_csv(DEST, sep="\t", index=False)

    if _DROPS:
        print("\n[arm_card] key-normalisation drops (visible, never silent):")
        for src, reasons in sorted(_DROPS.items()):
            for why, vals in sorted(reasons.items()):
                ex = ", ".join(sorted(vals)[:3])
                print(f"    {src:22s} {why:18s} {len(vals):4d}   e.g. {ex}")

    print(f"\n[arm_card] {len(card):,} arms × {card.shape[1]} columns -> {DEST}")
    _coverage_report(card)
    return card


def _coverage_report(card: pd.DataFrame) -> None:
    """⚠ AXIOM 5 — report coverage on the denominator that MATTERS, not just the union."""
    sets = [("ALL arms", card),
            ("EXPRESSED (TCGA)", card[card.get("cov_expr", False)]),
            # ⚠ 2026-08-19: was keyed on `model_n_edges`, PRUNED as bit-identical to `model_n_genes`.
            # The `else card.iloc[:0]` fallback would have silently reported this denominator as EMPTY
            # rather than failing — a pruned column becoming a silent zero is exactly the ripple
            # axiom 2 exists for. Repointed to the surviving twin.
            ("MODEL arms (>=1 HE edge)", card[card.get("model_n_genes").notna()]
             if "model_n_genes" in card.columns else card.iloc[:0])]
    flags = [c for c in card.columns if c.startswith("cov_")]
    print(f"\n{'denominator':<26s}{'n':>7s}  " + "  ".join(f"{f[4:]:>10s}" for f in flags))
    for lab, sub in sets:
        if not len(sub):
            continue
        print(f"{lab:<26s}{len(sub):>7,d}  " +
              "  ".join(f"{100*sub[f].mean():9.1f}%" for f in flags))


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
