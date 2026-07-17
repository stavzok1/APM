"""CPTAC proteome validation of the mirna_hallmark pressure spine.

miRNA action represses **protein** (mRNA destabilization *and* translational block).
The subproject's findings rest on miRNA pressure anti-correlating with target-gene
**mRNA**; CPTAC adds the orthogonal **protein** layer for the *same TCGA-BRCA patients*
(CPTAC TCGA-105 iTRAQ, Mertins 2016), so we can test the pressure claim at the readout
miRNAs actually control.

Cohort: **CPTAC TCGA-105** only (``analysis.mechanisms.cptac_common`` ``tcga105`` config).
These 105 proteome columns are the *same patients* as TCGA, joined deterministically via
``PATHS.cptac_tcga_participant_map`` (no RNA-correlation matching). The 122-column
prospective cohort is a *different* patient set and is handled separately (Task B / its own
miRNA-seq), not here. See ``pipeline/md/CPTAC_PROTEOME_DATA.md``.

Pressure variants (all four reported co-equal):
  edge set ∈ {all-edge spine, high-evidence only} × gate ∈ {ungated, AGO-gated}.

Three layers (every test: marginal Spearman/Pearson + partial adjusting CPE + HRD,
per the subproject conduct guardrails; raw reported alongside adjusted):

  L1  rna_protein_gap   gap = cptac_rna_z − cptac_protein_z; pressure should be **positively**
                        associated with the gap (protein repressed beyond mRNA). Directly
                        re-runs the question of parent D63 (``cptac_mirna_protein_repression``,
                        pooled-null on legacy pressure) with the gated/high-evidence pressure.
  L2  protein_anticorr  pressure vs cptac_protein_z; expected **negative** (the core mRNA
                        anti-correlation finding, re-tested at protein).
  L3  per-Hallmark / per-stratum   L1+L2 aggregated to Hallmark programs (mean over member
                        genes) and split by PAM50 / TNBC / immune / stage strata, to test
                        whether the most heavily-targeted programs couple most strongly.

Outputs (``output/cptac_validation/``):
  - ``gene_level_associations.tsv.gz``   variant × layer × gene (cohort), raw+partial, FDR.
  - ``pooled_associations.tsv``          variant × layer × stratum, within-gene-centered pooled.
  - ``hallmark_associations.tsv``        variant × layer × hallmark × stratum, ranked.
  - ``method_manifest.json``

Run::

    .venv/bin/python3 -m mirna_hallmark.eval.cptac_validation
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D
from mirna_hallmark import stats as S
from mirna_hallmark.ago_gate import compute_ago_gate
from mirna_hallmark.config import AgoGateParams
from mirna_hallmark.hallmark_interaction import hallmark_pressure_matrix
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.pressure_build import compute_gene_pressure, method_blurb

from analysis.mechanisms.cptac_common import (
    cohort_z,
    get_cohort_config,
    load_cct,
    load_mapped,
    sid_participant_maps,
)
from analysis.expression.mirna_target_integration import (
    load_mimat_to_arm,
    load_mirna_expression,
)
from pipeline.config import PATHS

OUT_DIR = C.OUTPUT_ROOT / "cptac_validation"
MIN_N = 15

# CPTAC-2 prospective breast miRNA (built by scripts/cptac/ingest_cptac2_mirna_arms.py).
# Lives alongside the LinkedOmics matrices under data/CPTAC/.
PROSPECTIVE_MIMAT_LOG2RPM = PATHS.cptac_proteome_cct.parent / "cptac2_brca_mirna_isoform_log2rpm.tsv"

# (edge set, gate) -> variant label. "all" = evidence-weighted spine; "highev" = high-evidence only.
PRESSURE_VARIANTS: Tuple[Tuple[str, bool, str], ...] = (
    ("all", False, "all_edge|ungated"),
    ("all", True, "all_edge|gated"),
    ("highev", False, "highev|ungated"),
    ("highev", True, "highev|gated"),
)


# --------------------------------------------------------------------------- #
# Pressure construction
# --------------------------------------------------------------------------- #
_PRESSURE_CACHE_DIR = C.OUTPUT_ROOT / "cptac_validation" / "_pressure_cache"


def _df_sig(df: Optional[pd.DataFrame]) -> str:
    """Content signature of a frame (shape + row-hash sum), NaN/dtype sensitive."""
    if df is None:
        return "none"
    h = int(pd.util.hash_pandas_object(df, index=True).sum() & 0xFFFFFFFFFFFF)
    return f"{df.shape[0]}x{df.shape[1]}:{h}"


def _cached_gene_pressure(universe: Sequence[str], edges: pd.DataFrame,
                          mirna: Optional[pd.DataFrame], *, tag: str) -> pd.DataFrame:
    """``compute_gene_pressure`` with a content-hashed disk cache.

    The all-edge variant builds pressure over ~914k edges (the dominant CPTAC cost).
    The result is a pure function of (edges, miRNA matrix, gene universe, pressure
    mode, abundance floor), so the cache key hashes all of them -- any input change
    misses and recomputes, so a stale hit is impossible. Set env CPTAC_FORCE_PRESSURE=1
    to bypass. Stored as pickle (exact dtype round-trip)."""
    m = mirna if mirna is not None else D.load_mirna_arms()
    gene_sig = hashlib.sha1("\n".join(map(str, sorted(universe))).encode()).hexdigest()[:12]
    key_src = "|".join([tag, _df_sig(edges), _df_sig(m), gene_sig,
                        method_blurb(), str(C.PRESSURE_ABUNDANCE_FLOOR)])
    key = hashlib.sha1(key_src.encode()).hexdigest()[:16]
    fp = _PRESSURE_CACHE_DIR / f"gp_{tag}_{key}.pkl"
    if not os.environ.get("CPTAC_FORCE_PRESSURE") and fp.exists():
        print(f"[cptac_val] gene pressure: {tag} (cache hit {fp.name})")
        return pd.read_pickle(fp)
    gp = compute_gene_pressure(universe, edges=edges, mirna=m)
    try:
        _PRESSURE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
        gp.to_pickle(fp)
    except Exception as e:  # noqa: BLE001 -- cache write must never break the run
        print(f"[cptac_val]   (pressure cache write skipped: {e})")
    return gp


def build_pressure_matrices(
    hs: HallmarkSets,
    *,
    include_tnrc6: Optional[bool] = None,
    mirna: Optional[pd.DataFrame] = None,
    rna_for_gate: Optional[pd.DataFrame] = None,
) -> Dict[str, pd.DataFrame]:
    """Return the four gene×sample pressure matrices keyed by variant label.

    ``mirna`` / ``rna_for_gate`` default to the subproject TCGA matrices. Pass the
    CPTAC-2 arm matrix + CPTAC RNA to build pressure **within** the prospective
    cohort (sample space = those columns).
    """
    universe = list(hs.universe)
    edges_all = D.load_hallmark_edges()
    edges_he = D.high_evidence_edges(edges_all)

    print(f"[cptac_val] gene pressure: all-edge ({len(edges_all):,} edges) ...")
    gp_all = _cached_gene_pressure(universe, edges_all, mirna, tag="all")
    print(f"[cptac_val] gene pressure: high-evidence ({len(edges_he):,} edges) ...")
    gp_he = _cached_gene_pressure(universe, edges_he, mirna, tag="he")

    print("[cptac_val] AGO/RISC gate ...")
    rna = rna_for_gate if rna_for_gate is not None else D.load_rna()
    from dataclasses import replace
    params = C.AGO_GATE if include_tnrc6 is None else replace(C.AGO_GATE, include_tnrc6=include_tnrc6)
    gate = compute_ago_gate(rna, params=params)["ago_gate"]

    out: Dict[str, pd.DataFrame] = {}
    for edge_set, gated, label in PRESSURE_VARIANTS:
        gp = gp_all if edge_set == "all" else gp_he
        if gated:
            shared = gp.columns.intersection(gate.index)
            gp = gp[shared].mul(gate.reindex(shared), axis=1)
        out[label] = gp
    return out


# --------------------------------------------------------------------------- #
# CPTAC protein / RNA z-scores (participant-keyed)
# --------------------------------------------------------------------------- #
def load_cptac_layers(cohort: str = "tcga105") -> Dict[str, pd.DataFrame]:
    """Sample-keyed cohort z-scores of CPTAC protein & RNA + the rna−protein gap.

    ``tcga105`` columns are remapped to TCGA participants (same patients as the
    pressure spine). ``prospective`` keeps the LinkedOmics ``X<case>`` columns —
    pressure for that cohort is built on the *same* columns, so no remap.
    """
    if cohort == "tcga105":
        cfg = get_cohort_config(cohort)
        mapped = load_mapped(cohort)
        sids = mapped["cptac_sample_id"].tolist()
        sid_to_part, _ = sid_participant_maps(mapped)
        prot = load_cct(cfg.proteome_cct)
        rna = load_cct(cfg.rna_cct)
        common = sorted(set(prot.index) & set(rna.index))
        prot_z = cohort_z(prot.loc[common], sids).rename(columns=sid_to_part)
        rna_z = cohort_z(rna.loc[common], sids).rename(columns=sid_to_part)
    else:  # prospective: self-contained, X<case> columns
        cfg = get_cohort_config("pancan122")
        prot = load_cct(cfg.proteome_cct)
        rna = load_cct(cfg.rna_cct)
        sids = sorted(set(prot.columns) & set(rna.columns))
        common = sorted(set(prot.index) & set(rna.index))
        prot_z = cohort_z(prot.loc[common], sids)
        rna_z = cohort_z(rna.loc[common], sids)
    gap = rna_z - prot_z  # positive = protein repressed relative to mRNA (assumes slope 1)
    resid = _protein_on_rna_residual(prot_z, rna_z)  # protein beyond mRNA (per-gene OLS)
    return {"protein_z": prot_z, "rna_z": rna_z, "gap": gap, "protein_resid": resid}


def _protein_on_rna_residual(prot_z: pd.DataFrame, rna_z: pd.DataFrame) -> pd.DataFrame:
    """Per-gene residual of protein_z regressed on its own rna_z (slope-free gap).

    For each gene g: resid = protein_z − (a + b·rna_z) with (a,b) the per-gene OLS fit
    across samples. Negative residual ~ high pressure = protein repressed **beyond**
    what mRNA explains — the rigorous translational-repression readout (vs the
    slope-1 ``gap``). Mean-centered, so low resid = protein below its mRNA-predicted level.
    """
    common = prot_z.index.intersection(rna_z.index)
    cols = prot_z.columns
    out = pd.DataFrame(index=common, columns=cols, dtype=float)
    Y = prot_z.loc[common]
    X = rna_z.loc[common, cols]
    for gene in common:
        y = Y.loc[gene].astype(float)
        x = X.loc[gene].astype(float)
        m = y.notna() & x.notna()
        if m.sum() < MIN_N or x[m].std(ddof=0) == 0:
            continue
        xv, yv = x[m], y[m]
        b = np.cov(xv, yv, ddof=0)[0, 1] / np.var(xv)
        a = yv.mean() - b * xv.mean()
        out.loc[gene, m[m].index] = (yv - (a + b * xv)).values
    return out


def load_prospective_mirna_arms() -> pd.DataFrame:
    """CPTAC-2 breast miRNA arm × X<case> log2(RPM+1), in the subproject arm space."""
    if not PROSPECTIVE_MIMAT_LOG2RPM.exists():
        raise FileNotFoundError(
            f"Missing {PROSPECTIVE_MIMAT_LOG2RPM}. Build it first:\n"
            "  .venv/bin/python3 scripts/cptac/fetch_cptac2_breast_mirna.py\n"
            "  .venv/bin/python3 scripts/cptac/ingest_cptac2_mirna_arms.py"
        )
    mimat_to_arm = load_mimat_to_arm(C.MIRNA_MATURE_LOCI)
    return load_mirna_expression(PROSPECTIVE_MIMAT_LOG2RPM, mimat_to_arm=mimat_to_arm)


def load_prospective_clinical() -> pd.DataFrame:
    """Prospective-cohort strata + confounders from HS_CPTAC_BRCA_2018 CLI .tsi.

    Mapped to the subproject's stratum column names where possible. HRD is not
    available; the confounder set is purity (ESTIMATE.TumorPurity) + chromosomal
    instability (CIN) — documented in the manifest as the TCGA CPE+HRD substitute.
    """
    cli = pd.read_csv(PATHS.cptac_cli_tsi, sep="\t", index_col=0)
    cli = cli[~cli.index.isin(["IDX"])].copy()
    out = pd.DataFrame(index=cli.index)
    out["participant"] = cli.index.astype(str)
    out["PAM50_final"] = cli.get("PAM50")
    tnbc = cli.get("TNBC.Updated.Clinical.Status")
    out["tnbc_subtype_4"] = tnbc.map(lambda v: "TNBC" if str(v).strip() in {"Yes", "1", "TNBC"} else "non-TNBC") \
        if tnbc is not None else np.nan
    out["pathologic_stage_collapsed"] = cli.get("Stage.reduced", cli.get("Stage"))
    out["purity"] = pd.to_numeric(cli.get("ESTIMATE.TumorPurity"), errors="coerce")
    out["cin"] = pd.to_numeric(cli.get("Chromosome.INstability.index.CIN."), errors="coerce")
    return out.reset_index(drop=True)


# Layer -> (CPTAC response matrix key, expected sign of association with pressure).
LAYERS: Tuple[Tuple[str, str, str], ...] = (
    ("rna_protein_gap", "gap", "positive"),         # pressure -> bigger gap (protein < mRNA), slope-1
    ("protein_resid", "protein_resid", "negative"), # pressure -> protein below its mRNA-predicted level (slope-free)
    ("protein_anticorr", "protein_z", "negative"),  # pressure -> lower protein (mRNA+protein combined)
)


# --------------------------------------------------------------------------- #
# Gene-level (cohort) associations
# --------------------------------------------------------------------------- #
def gene_level_associations(
    pressures: Dict[str, pd.DataFrame],
    layers: Dict[str, pd.DataFrame],
    cov: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[dict] = []
    for variant, gp in pressures.items():
        for layer, resp_key, _sign in LAYERS:
            resp = layers[resp_key]
            genes = sorted(set(gp.index) & set(resp.index))
            samples = gp.columns.intersection(resp.columns)
            for gene in genes:
                x = gp.loc[gene, samples]
                y = resp.loc[gene, samples]
                st = S.correlation_pair(y, x, cov)
                if st["n"] < MIN_N:
                    continue
                rows.append({
                    "variant": variant,
                    "layer": layer,
                    "gene_name": gene,
                    "n": st["n"],
                    "spearman_rho": _r(st["spearman_rho"]),
                    "spearman_p": st["spearman_p"],
                    "partial_rho": _r(st["partial_rho"]),
                    "partial_p": st["partial_p"],
                    "pearson_r": _r(st["pearson_r"]),
                    "pearson_p": st["pearson_p"],
                })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # FDR within each (variant, layer).
    out["spearman_q"] = np.nan
    out["partial_q"] = np.nan
    for _, idx in out.groupby(["variant", "layer"]).groups.items():
        sl = out.loc[idx]
        out.loc[idx, "spearman_q"] = S.bh_fdr(sl["spearman_p"].fillna(1.0).values)
        if sl["partial_p"].notna().any():
            out.loc[idx, "partial_q"] = S.bh_fdr(sl["partial_p"].fillna(1.0).values)
    return out.sort_values(["variant", "layer", "partial_rho"])


# --------------------------------------------------------------------------- #
# Pooled within-gene-centered test (cohort + per stratum) -- comparable to D63 pooled.
# --------------------------------------------------------------------------- #
def _pooled_one(
    gp: pd.DataFrame, resp: pd.DataFrame, samples: Sequence[str]
) -> Tuple[float, float, int, int]:
    """Within-gene-center pressure & response across the given samples, then Spearman."""
    genes = sorted(set(gp.index) & set(resp.index))
    recs = []
    for gene in genes:
        x = gp.loc[gene, samples]
        y = resp.loc[gene, samples]
        df = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
        if len(df) < 5:
            continue
        df = df - df.mean()  # within-gene center removes per-gene baselines
        df["gene"] = gene
        recs.append(df)
    if not recs:
        return (np.nan, np.nan, 0, 0)
    pooled = pd.concat(recs)
    if len(pooled) < MIN_N:
        return (np.nan, np.nan, int(len(pooled)), int(pooled["gene"].nunique()))
    rho, p = stats.spearmanr(pooled["x"], pooled["y"])
    return (float(rho), float(p), int(len(pooled)), int(pooled["gene"].nunique()))


def pooled_associations(
    pressures: Dict[str, pd.DataFrame],
    layers: Dict[str, pd.DataFrame],
    clinical: pd.DataFrame,
) -> pd.DataFrame:
    strata = _stratum_groups(clinical)
    rows: List[dict] = []
    for variant, gp in pressures.items():
        for layer, resp_key, sign in LAYERS:
            resp = layers[resp_key]
            base = list(gp.columns.intersection(resp.columns))
            for group, value, members in strata:
                samples = [s for s in base if s in members]
                if len(samples) < 20:
                    continue
                rho, p, n_gp, n_genes = _pooled_one(gp, resp, samples)
                rows.append({
                    "variant": variant, "layer": layer, "expected_sign": sign,
                    "stratum_group": group, "stratum_value": value,
                    "n_samples": len(samples), "n_gene_participant": n_gp,
                    "n_genes": n_genes, "pooled_spearman": _r(rho), "pooled_p": p,
                })
    return pd.DataFrame(rows).sort_values(["variant", "layer", "stratum_group", "stratum_value"])


# --------------------------------------------------------------------------- #
# Hallmark-level (cohort + per stratum) associations
# --------------------------------------------------------------------------- #
def _member_mean_matrix(value_gene_x_part: pd.DataFrame, hs: HallmarkSets) -> pd.DataFrame:
    """Hallmark × participant = mean of a gene×participant matrix over present members."""
    return hallmark_pressure_matrix(value_gene_x_part, hs)


def hallmark_associations(
    pressures: Dict[str, pd.DataFrame],
    layers: Dict[str, pd.DataFrame],
    clinical: pd.DataFrame,
    hs: HallmarkSets,
    cov: pd.DataFrame,
    enr: pd.DataFrame,
) -> pd.DataFrame:
    strata = _stratum_groups(clinical)
    fold = enr.set_index("hallmark_set")["fold_enrichment"].to_dict()
    resp_hm = {key: _member_mean_matrix(layers[key], hs) for _, key, _ in LAYERS}

    rows: List[dict] = []
    for variant, gp in pressures.items():
        hp = _member_mean_matrix(gp, hs)
        for layer, resp_key, sign in LAYERS:
            rh = resp_hm[resp_key]
            base = hp.columns.intersection(rh.columns)
            hsets = hp.index.intersection(rh.index)
            for group, value, members in strata:
                samples = base.intersection(pd.Index(list(members)))
                if len(samples) < 20:
                    continue
                for hset in hsets:
                    st = S.correlation_pair(rh.loc[hset, samples], hp.loc[hset, samples], cov)
                    if st["n"] < MIN_N:
                        continue
                    rows.append({
                        "variant": variant, "layer": layer, "expected_sign": sign,
                        "stratum_group": group, "stratum_value": value,
                        "hallmark_set": hset, "n": st["n"],
                        "fold_highev_enrichment": fold.get(hset, np.nan),
                        "spearman_rho": _r(st["spearman_rho"]), "spearman_p": st["spearman_p"],
                        "partial_rho": _r(st["partial_rho"]), "partial_p": st["partial_p"],
                    })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["spearman_q"] = np.nan
    for _, idx in out.groupby(["variant", "layer", "stratum_group", "stratum_value"]).groups.items():
        sl = out.loc[idx]
        out.loc[idx, "spearman_q"] = S.bh_fdr(sl["spearman_p"].fillna(1.0).values)
    return out.sort_values(["variant", "layer", "stratum_group", "stratum_value", "partial_rho"])


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _r(v) -> float:
    return round(float(v), 4) if v is not None and np.isfinite(v) else np.nan


def _stratum_groups(clinical: pd.DataFrame) -> List[Tuple[str, str, set]]:
    """List of (group, value, member-participant-set), with a cohort-wide 'all' row first."""
    clin = clinical.drop_duplicates("participant")
    groups: List[Tuple[str, str, set]] = [("cohort", "all", set(clin["participant"]))]
    for col, layer_name in C.STRATUM_SPECS:
        if col not in clin.columns:
            continue
        for value, sub in clin.dropna(subset=[col]).groupby(col):
            groups.append((layer_name, str(value), set(sub["participant"])))
    return groups


def _covariates(clinical: pd.DataFrame, cohort: str, *, composition: bool = True) -> pd.DataFrame:
    """Numeric confounder frame, participant-indexed. tcga105=CPE+HRD; prospective=purity+CIN — **plus the
    CELL-COMPOSITION block** (8 Wu-major lineages + malignant proliferation) via `confounders.build_C`.

    ⚠⚠ **COMPOSITION IS NOT OPTIONAL — omitting it was a REAL CONFOUND (MH-107, 2026-07-12).**
    Until now this block was `purity + CIN` only, while **the TCGA learned model that this module VALIDATES
    already conditions on the 8 Wu-major lineages** (`data._DECONV_COLS`) — i.e. the validation ran with a
    WEAKER confounder block than the model it was validating. Measured consequence on the prospective cohort:
      * gene-level FDR-negative hits **66 → 3** (median |ρ| 0.311 → 0.163, −48%; **9 genes SIGN-FLIP**);
      * the flagship **EMT** hallmark **−0.123 → −0.009** (collapses) — EMT proteins correlate with the **CAF
        fraction at +0.509** (vs +0.033 for all other proteins) and **ZEB1 protein at +0.768**, while miR-200 (an
        EPITHELIAL miRNA) correlates with CAF at −0.35 ⇒ the "miR-200 represses EMT protein" signal in bulk is
        **compartment arithmetic**, not cell-autonomous repression;
      * **but this is a REFRAME, not a blanket retraction** — tumour-cell-INTRINSIC programs SURVIVE
        (ESTROGEN_RESPONSE_EARLY −0.214→−0.141, PI3K_AKT_MTOR −0.068→−0.087) while stroma-loaded ones collapse.
    The underlying **edge is still real**: miR-200→ZEB1 at TCGA mRNA, composition-adjusted, n=1041, **ρ=−0.209,
    p=8.7e-12** — so composition is a CONFOUND (it inflates), not a mediator whose removal kills the effect.

    `composition=False` restores the historical (confounded) block for provenance/reproduction only.
    """
    clin = clinical.drop_duplicates("participant").set_index("participant")
    wanted = C.CONFOUNDER_NUMERIC if cohort == "tcga105" else ("purity", "cin")
    cols = [c for c in wanted if c in clin.columns]
    out = clin[cols].apply(pd.to_numeric, errors="coerce")
    if not composition:
        return out
    try:
        from mirna_hallmark.learned import confounders as CF
        key = "tcga" if cohort == "tcga105" else "cptac"
        comp = CF.build_C(key, list(out.index), purity=False)      # lineages + prolif; purity already above
        comp = comp.reindex(out.index)
        out = out.join(comp[[c for c in comp.columns if c not in out.columns]])
    except Exception as exc:                                        # never silently drop the block
        raise RuntimeError(f"composition block unavailable for {cohort!r} — refusing to run a CONFOUNDED "
                           f"validation (pass composition=False to reproduce the historical numbers): {exc}")
    return out


def _target_enrichment(hs: HallmarkSets) -> pd.DataFrame:
    from mirna_hallmark.hallmark_interaction import target_enrichment
    return target_enrichment(D.load_hallmark_edges(), hs)


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #
def run(*, cohort: str = "tcga105", out_dir: Path = OUT_DIR, include_tnrc6: Optional[bool] = None,
        batch_kind: str = "none") -> None:
    if cohort not in ("tcga105", "prospective"):
        raise ValueError(f"cohort must be 'tcga105' or 'prospective', got {cohort!r}")
    # tcga105 stays at the top-level dir; prospective writes to a subdir.
    if out_dir == OUT_DIR and cohort != "tcga105":
        out_dir = OUT_DIR / cohort
    if batch_kind != "none":
        out_dir = out_dir / f"batch_{batch_kind}"
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()

    if cohort == "tcga105":
        clinical = D.load_clinical_strata()
        pressures = build_pressure_matrices(hs, include_tnrc6=include_tnrc6)
        conf_label = "CPE + thornsson_hrd_score"
    else:  # prospective: pressure built within the CPTAC-2 cohort (X<case> columns)
        clinical = load_prospective_clinical()
        cptac_mirna = load_prospective_mirna_arms()
        cptac_rna_gate = load_cct(get_cohort_config("pancan122").rna_cct)
        print(f"[cptac_val] prospective miRNA arms: {cptac_mirna.shape[0]} arms × {cptac_mirna.shape[1]} samples")
        pressures = build_pressure_matrices(
            hs, include_tnrc6=include_tnrc6, mirna=cptac_mirna, rna_for_gate=cptac_rna_gate
        )
        conf_label = "ESTIMATE.TumorPurity + CIN (HRD unavailable in CPTAC CLI)"

    cov = _covariates(clinical, cohort)
    from mirna_hallmark import cptac_batch as B
    cov, bcols = B.augment_cov(cov, cohort, batch_kind)
    if bcols:
        print(f"[cptac_val] +batch ({batch_kind}): {len(bcols)} dummies -> {out_dir}")
    # The manifest must name the covariates that were ACTUALLY used, not a literal written before the
    # composition join — otherwise it under-reports the block and a reader mistakes an adjusted run for a
    # confounded one (this bit us: the first MH-107 re-run wrote "purity + CIN" while fitting 11 covariates).
    conf_label = f"{conf_label} | fitted: {', '.join(map(str, cov.columns))}"
    print(f"[cptac_val] covariates ({cov.shape[1]}): {list(cov.columns)}")
    layers = load_cptac_layers(cohort)
    enr = _target_enrichment(hs)

    print("[cptac_val] gene-level associations ...")
    gene_tbl = gene_level_associations(pressures, layers, cov)
    print("[cptac_val] pooled within-gene-centered ...")
    pooled_tbl = pooled_associations(pressures, layers, clinical)
    print("[cptac_val] Hallmark-level (cohort + strata) ...")
    hm_tbl = hallmark_associations(pressures, layers, clinical, hs, cov, enr)

    gene_tbl.to_csv(out_dir / "gene_level_associations.tsv.gz", sep="\t", index=False, compression="gzip")
    pooled_tbl.to_csv(out_dir / "pooled_associations.tsv", sep="\t", index=False)
    hm_tbl.to_csv(out_dir / "hallmark_associations.tsv", sep="\t", index=False)

    n_part = int(len(set(layers["protein_z"].columns) & set(next(iter(pressures.values())).columns)))
    cohort_desc = (
        "CPTAC TCGA-105 iTRAQ — same patients as TCGA, pressure from TCGA miRNA"
        if cohort == "tcga105" else
        "CPTAC-2 prospective breast — independent patients, pressure from CPTAC-2 miRNA-seq (self-contained)"
    )
    manifest = {
        "built_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "module": "mirna_hallmark.eval.cptac_validation",
        "doc": "pipeline/md/CPTAC_PROTEOME_DATA.md",
        "cohort": cohort,
        "cohort_description": cohort_desc,
        "n_participants_with_pressure_and_protein": n_part,
        "pressure_method": method_blurb(),
        "pressure_variants": [lbl for *_x, lbl in PRESSURE_VARIANTS],
        "layers": {name: f"{key} ~ pressure (expect {sign})" for name, key, sign in LAYERS},
        "confounders": f"{conf_label} (partial Spearman/Pearson; raw reported)",
        "batch_kind": batch_kind,
        "compare_to": "analysis.mechanisms.cptac_mirna_protein_repression (D63, legacy pressure pooled-null)",
        "outputs": [
            "gene_level_associations.tsv.gz",
            "pooled_associations.tsv",
            "hallmark_associations.tsv",
        ],
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    # Headline: cohort pooled per variant/layer (compare gap-row to D63's ρ≈-0.013 null).
    print(f"[cptac_val] wrote {out_dir} (n={n_part})")
    cohort_pool = pooled_tbl[(pooled_tbl["stratum_group"] == "cohort")]
    for _, r in cohort_pool.iterrows():
        print(f"  {r['variant']:<18} {r['layer']:<16} pooled ρ={r['pooled_spearman']:+.3f} "
              f"(p={r['pooled_p']:.1e}, n_gp={r['n_gene_participant']})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", default="tcga105", choices=["tcga105", "prospective"],
                    help="tcga105 = same TCGA patients; prospective = independent CPTAC-2 cohort")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--include-tnrc6", action="store_true", default=None,
                    help="Include TNRC6A/B/C in the AGO gate (default: follow config.AGO_GATE)")
    ap.add_argument("--batch-kind", default="none", choices=["none", "site", "plex", "auto"],
                    help="Add MS-plex/site batch dummies to the cov (writes to batch_<kind>/ subdir)")
    ap.add_argument("--min-purity", type=float, default=None, help="Accepted for runner compatibility (no-op)")
    args = ap.parse_args()
    run(cohort=args.cohort, out_dir=args.out_dir, include_tnrc6=args.include_tnrc6, batch_kind=args.batch_kind)


if __name__ == "__main__":
    main()
