"""Assemble aligned (Y, X, C, w_prior) for a gene from the spine loaders.

Reuses (no rewrite): data_loaders.load_mirna_arms (arm×participant log2 RPM+1),
load_rna (gene×participant log2 TPM+1), load_clinical_strata (participant-keyed; already
drops PAM50 Normal-like via config.EXCLUDE_NORMAL_LIKE), high_evidence_edges / load_hallmark_edges.

Confounders here are the current spine minimal set (CPE + HRD). The design's transcription-control
block and the malignant-compartment (CIBERSORTx) proliferation / deconvolution composition
(Decision B) are Phase-0/Phase-3 upgrades wired later — flagged in `confounder_columns()`.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd

from pathlib import Path

from mirna_hallmark import config as C
from mirna_hallmark import data_loaders as D

_CACHE: dict = {}

# The user's breast single-cell-atlas CIBERSORTx run (participant-keyed cell fractions). Design §Decision B:
# use the NON-malignant compartments as the composition block; NOT Cancer Epithelial (conditioning on the
# malignant fraction over-controls the signal, since the target is expressed there).
_DECONV_PATH = Path("mirna_hallmark/output/brca_deconvolution/tcga_cibersortx_fractions.tsv")
_DECONV_COLS = ["CAFs", "T-cells", "Myeloid", "B-cells", "Endothelial", "PVL",
                "Normal Epithelial", "Plasmablasts"]


def _latent(k: int = 10):
    """Global transcriptional-state proxy: top-k PCs of the most-variable genes (Design §B transcription
    control). Captures shared programs (proliferation, ER, immune) that co-drive a miRNA and its target.

    ⚠ EMPIRICALLY OVER-CORRECTS (tested 2026-07-04): adding k=10 PCs to C **eliminates real, causally-
    confirmed edges** — miR-21→PTEN (−0.17→+0.03), miR-17~92→BCL2 (−0.29→−0.08), miR-206→ESR1 (−0.32→−0.03)
    — because the top PCs ARE the programs the miRNAs regulate (proliferation/ER/EMT), so residualising on
    them removes the signal itself (Design §B "over-residualising removes real signal"). ⇒ **Do NOT use
    global PCs as the transcription proxy.** Use an ORTHOGONAL proxy: intronic/pre-mRNA reads (direct
    transcription rate) or the gene's specific TF-regulon activity. Kept here only as the documented
    negative control."""
    key = f"latent{k}"
    if key not in _CACHE:
        from sklearn.decomposition import PCA
        M = _load()["Y"].T                                      # participant × gene
        top = M.var().nlargest(2000).index                     # most-variable genes
        Ms = M[top]
        Ms = ((Ms - Ms.mean()) / Ms.std(ddof=0)).fillna(0.0)
        pcs = PCA(n_components=k, random_state=0).fit_transform(Ms.to_numpy())
        _CACHE[key] = pd.DataFrame(pcs, index=M.index, columns=[f"PC{i+1}" for i in range(k)])
    return _CACHE[key]


_COLLECTRI = Path("data/external_cache/collectri/collectri_human.tsv")
_TS_CONTEXT = Path("data/external_cache/targetscan/Predicted_Targets_Context_Scores.default_predictions.txt")
_TS_CACHE = Path("mirna_hallmark/output/learned/targetscan_context.parquet")


def _targetscan_context() -> Optional[pd.DataFrame]:
    """Per-(arm, gene) cumulative TargetScan context++ magnitude (Design §E: continuous, sequence-based
    prior for EVERY predicted edge incl. orphans). weighted context++ is negative (stronger = more
    negative); we sum over sites and flip sign → positive `ts_mag` (larger = stronger predicted repression)."""
    if "ts" not in _CACHE:
        if _TS_CACHE.exists():
            _CACHE["ts"] = pd.read_parquet(_TS_CACHE)
        elif _TS_CONTEXT.exists():
            arms = set(_load()["X"].index)
            keep = []
            for ch in pd.read_csv(_TS_CONTEXT, sep="\t",
                                  usecols=["Gene Symbol", "miRNA", "weighted context++ score"],
                                  chunksize=500_000):
                ch = ch[ch["miRNA"].isin(arms)]                # exact arm-name match (hsa-miR-..-5p/3p)
                if len(ch):
                    keep.append(ch)
            ts = pd.concat(keep, ignore_index=True).rename(
                columns={"Gene Symbol": "gene", "miRNA": "arm", "weighted context++ score": "wcs"})
            ts = ts.groupby(["arm", "gene"])["wcs"].sum().reset_index()
            ts["ts_mag"] = (-ts["wcs"]).clip(lower=0.0)         # positive magnitude, larger = stronger
            _TS_CACHE.parent.mkdir(parents=True, exist_ok=True)
            _CACHE["ts"] = ts[["arm", "gene", "ts_mag"]]
            _CACHE["ts"].to_parquet(_TS_CACHE)
        else:
            _CACHE["ts"] = None
    return _CACHE["ts"]


def _tf_regulon():
    """target gene → DataFrame[tf, strength] of its TF regulators, with **edge STRENGTH** (CollecTRI
    `n_references`; complexes expanded, strength inherited). A gene's TFs are not equal — TP53→PTEN (79
    refs) is far better-supported than a 2-ref edge — so the transcription proxy weights each TF by its
    regulatory strength, not just its presence (Design §B). TFs act at transcription, orthogonal to the
    miRNA's post-transcriptional action."""
    if "tf_reg" not in _CACHE:
        if _COLLECTRI.exists():
            c = pd.read_csv(_COLLECTRI, sep="\t", usecols=["source_genesymbol", "target_genesymbol",
                            "references", "is_stimulation", "is_inhibition"])
            c["strength"] = np.log1p(c["references"].fillna("").apply(
                lambda s: len(str(s).split(";")) if s else 0))
            c["mor"] = c["is_stimulation"].astype(int) - c["is_inhibition"].astype(int)  # +1/-1/0(dual)
            c = c.assign(tf=c["source_genesymbol"].str.split("_")).explode("tf")
            _CACHE["tf_reg"] = {g: sub.groupby("tf")[["strength", "mor"]].first()
                                for g, sub in c.groupby("target_genesymbol")}
        else:
            _CACHE["tf_reg"] = None
    return _CACHE["tf_reg"]


def _tf_activity(gene: str, parts) -> Optional[pd.DataFrame]:
    """decoupleR-style SIGNED TF ACTIVITY (the 'proper' proxy): one covariate =
    Σ_TF mor(TF→g)·strength(TF→g)·z(TF_expr) — the net predicted transcriptional drive on g, using
    CollecTRI's activation/repression sign AND reference-count strength. Contrast with `_tf_proxy` (raw
    strength-weighted expression PCs). Dual-sign (mor=0) TFs are dropped (no orientation)."""
    reg = _tf_regulon()
    if reg is None or gene not in reg:
        return None
    Yall = _load()["Y"]
    r = reg[gene]
    r = r[(r["mor"] != 0) & r.index.isin(Yall.index) & (r.index != gene)]
    if len(r) < 2:
        return None
    E = Yall.loc[r.index, parts].T                              # participant × TF
    E = ((E - E.mean()) / E.std(ddof=0)).fillna(0.0)
    w = (r["mor"] * r["strength"]).to_numpy()
    act = pd.Series(E.to_numpy() @ w, index=parts, name="TFact")
    return act.to_frame()


def _tf_proxy(gene: str, parts, k: int = 3) -> Optional[pd.DataFrame]:
    """STRENGTH-WEIGHTED top-k PCs of the gene's TF-regulon expression. Each TF's standardised expression
    is scaled by its CollecTRI edge strength (log1p n_references) BEFORE PCA, so the components are driven
    by the gene's *strongest* regulators (its master TFs), not merely its highest-variance ones — the
    surgical transcription proxy (Design §B). Signs are learned by the downstream residualisation."""
    reg = _tf_regulon()
    if reg is None or gene not in reg:
        return None
    Yall = _load()["Y"]
    s = reg[gene]["strength"]
    tfs = [t for t in s.index if t in Yall.index and t != gene]
    if len(tfs) < 2:
        return None
    from sklearn.decomposition import PCA
    E = Yall.loc[tfs, parts].T                                  # participant × TF
    E = ((E - E.mean()) / E.std(ddof=0)).fillna(0.0)
    E = E.mul(s.reindex(tfs).to_numpy(), axis=1)               # scale each TF by its edge STRENGTH
    kk = min(k, E.shape[1])
    pcs = PCA(n_components=kk, random_state=0).fit_transform(E.to_numpy())
    return pd.DataFrame(pcs, index=parts, columns=[f"TF{i+1}" for i in range(kk)])


def _deconv():
    if "deconv" not in _CACHE:
        if _DECONV_PATH.exists():
            f = pd.read_csv(_DECONV_PATH, sep="\t").drop_duplicates("Mixture").set_index("Mixture")
            _CACHE["deconv"] = f[[c for c in _DECONV_COLS if c in f.columns]].apply(
                pd.to_numeric, errors="coerce")
        else:
            _CACHE["deconv"] = None
    return _CACHE["deconv"]


def confounder_columns(*, drop_hrd: bool = True) -> list[str]:
    """Design §4.1 minimal-core scalar confounders. Default = **CPE only** (Design: HRD/batch off by
    default). target-CN, malignant proliferation, TF-regulon transcription proxy, and program-conditional
    deconvolution fractions are added by `assemble_gene` (below), not here."""
    cols = list(C.CONFOUNDER_NUMERIC)                          # ('CPE', 'thornsson_hrd_score')
    return [c for c in cols if c != "thornsson_hrd_score"] if drop_hrd else cols


def _arm_name_map(Xall) -> dict:
    """Edge arm name -> its ACTUAL name in the miRNA matrix `X`. Three transforms, all safe:
    (1) `.N` duplicate-suffix — `X` stores `hsa-miR-101-3p.1` while edges use plain `hsa-miR-101-3p`, so the old
        `m in Xall.index` test SILENTLY DROPPED them (incl. canonical miR-101/124/126-3p, +~300 edges);
    (2) case-insensitive;
    (3) old-style SUFFIXLESS name -> the GUIDE (most-abundant) arm of the hairpin — the pre-arm-convention `hsa-miR-181a`
        means the guide (miR-181a-5p, median 10.4 RPM), NOT the -3p passenger (7.8). Resolving by abundance is
        well-defined (not the ambiguous blind -3p/-5p guess an earlier version wrongly used). NB the two arms share
        the LOCUS, so for CN/localization the choice is moot anyway.
    See memory `universe-redefinition-pending-refresh` / DATA_SOURCES arm→locus canonicalization.

    ⚠ The cache is keyed on a FINGERPRINT of `Xall.index`, not a bare literal: this map is now built for
    MULTIPLE cohorts (TCGA `X`, CPTAC prospective arms, GTEx/NAT), whose arm indices differ. A single global
    key silently returned the FIRST cohort's map to every later caller (verified: CPTAC-first poisoned TCGA,
    6272 → 5963 entries). Same index ⇒ same key ⇒ same object, so this is a strict no-op for existing callers."""
    import hashlib as _hl
    idx = list(Xall.index)
    # fingerprint = index + n_columns: the index alone would collide `X` with a column-SUBSET `X[parts]`, whose
    # `max(axis=1)` (the guide-arm tiebreak below) can differ. Distinct shapes ⇒ distinct maps.
    _key = ("arm_name_map_" + _hl.md5("\x00".join(map(str, idx)).encode()).hexdigest()[:16]
            + f"_{getattr(Xall, 'shape', (0, 0))[1]}")
    if _key not in _CACHE:
        import re as _re
        m: dict = {}
        for a in idx:
            m.setdefault(a, a)                                     # exact
            base = _re.sub(r"\.[0-9]+$", "", a)                    # strip trailing .N duplicate suffix
            m.setdefault(base, a)
            m.setdefault(base.lower(), a)                          # case-insensitive
        abund = Xall.max(axis=1)                                   # guide = most-abundant arm of the hairpin
        hp: dict = {}
        for a in idx:
            mm = _re.match(r"(.+?)-([35]p)(\.[0-9]+)?$", a)
            if mm:
                hp.setdefault(mm.group(1), []).append(a)
        for b, arms in hp.items():
            g = max(arms, key=lambda x: abund.get(x, -1.0))
            m.setdefault(b, g)                                     # suffixless old-style name -> guide arm
            m.setdefault(b.lower(), g)
        _CACHE[_key] = m
    return _CACHE[_key]


def _target_cn(genes) -> pd.DataFrame:
    """target-gene ASCAT3 copy number × participant (Design §B C-block: the gene's OWN genomic state)."""
    key = "tcn"
    if key not in _CACHE:
        _CACHE[key] = {}
    need = [g for g in genes if g not in _CACHE[key]]
    if need:
        if "tcn_full" not in _CACHE:                           # read the cached gene×participant CN matrix ONCE
            try:                                               # (was a full re-read per gene via load_cnv_target_genes)
                cp = Path(C.OUTPUT_ROOT) / "matrices" / "cnv_target_genes.tsv.gz"
                _CACHE["tcn_full"] = pd.read_csv(cp, sep="\t", index_col=0) if cp.exists() else None
            except Exception:
                _CACHE["tcn_full"] = None
        full = _CACHE["tcn_full"]
        if full is not None and all(g in full.index for g in need):
            for g in need:
                _CACHE[key][g] = full.loc[g]
        else:                                                  # some genes absent from the cache → source parse (rare), refresh
            cn = D.load_cnv_target_genes(sorted(set(need)))    # gene × participant
            _CACHE.pop("tcn_full", None)                       # invalidate: next miss re-reads the refreshed cache file
            for g in need:
                _CACHE[key][g] = cn.loc[g] if g in getattr(cn, "index", []) else None
    return _CACHE[key]


def _malignant_prolif(parts) -> Optional[pd.Series]:
    """One malignant-compartment proliferation (Design §B): bulk E2F+G2M metagene RESIDUALISED on the
    Cancer-Epithelial CIBERSORTx fraction (the design's immediately-available approximation — removes the
    composition-driven part, leaving per-malignant-cell proliferation)."""
    if "mprolif" not in _CACHE:
        from mirna_hallmark.hallmark_sets import HallmarkSets
        Y = _load()["Y"]
        hs = HallmarkSets.load()
        pg = sorted((set(hs.sets.get("HALLMARK_E2F_TARGETS", [])) |
                     set(hs.sets.get("HALLMARK_G2M_CHECKPOINT", []))) & set(Y.index))
        sub = Y.loc[pg]
        z = sub.sub(sub.mean(1), axis=0).div(sub.std(1) + 1e-9, axis=0)
        prolif = z.mean(0)                                     # bulk proliferation score / participant
        try:
            f = pd.read_csv(_DECONV_PATH, sep="\t").drop_duplicates("Mixture").set_index("Mixture")
            cep = pd.to_numeric(f.get("Cancer Epithelial"), errors="coerce").dropna()
            common = prolif.index.intersection(cep.index)
            lr = LinearRegression().fit(cep.loc[common].to_numpy().reshape(-1, 1), prolif.loc[common])
            resid = prolif.copy()
            resid.loc[common] = prolif.loc[common] - lr.predict(cep.loc[common].to_numpy().reshape(-1, 1))
            _CACHE["mprolif"] = resid                          # composition-adjusted (malignant) proliferation
        except Exception:
            _CACHE["mprolif"] = prolif
    return _CACHE["mprolif"].reindex(parts)


def _mycaf(parts) -> Optional[pd.Series]:
    """Fibroblast-subtype (myCAF vs iCAF) axis the lumped 'CAFs' fraction misses (user 2026-07-05): a
    myofibroblast/contractile marker metagene RESIDUALISED on the total CAF fraction → the TGF-β-driven
    myCAF state *beyond* total fibroblast content. Structural markers (not TGF-β-signaling genes) ⇒
    non-circular for the TGF-β-target robustness check."""
    if "mycaf" not in _CACHE:
        Y = _load()["Y"]
        mk = [g for g in ["ACTA2", "TAGLN", "POSTN", "FAP", "COL11A1", "THBS2", "MYH11", "CNN1"] if g in Y.index]
        sub = Y.loc[mk]; z = sub.sub(sub.mean(1), axis=0).div(sub.std(1) + 1e-9, axis=0); score = z.mean(0)
        try:
            f = pd.read_csv(_DECONV_PATH, sep="\t").drop_duplicates("Mixture").set_index("Mixture")
            caf = pd.to_numeric(f.get("CAFs"), errors="coerce").dropna()
            common = score.index.intersection(caf.index)
            lr = LinearRegression().fit(caf.loc[common].to_numpy().reshape(-1, 1), score.loc[common])
            resid = score.copy()
            resid.loc[common] = score.loc[common] - lr.predict(caf.loc[common].to_numpy().reshape(-1, 1))
            _CACHE["mycaf"] = resid
        except Exception:
            _CACHE["mycaf"] = score
    return _CACHE["mycaf"].reindex(parts)


def _mal_emt(parts) -> Optional[pd.Series]:
    """Within-malignant EMT/mesenchymal STATE (the axis `mal_prolif` doesn't capture, user 2026-07-05):
    structural EMT markers RESIDUALISED on the non-malignant (stroma+immune) CIBERSORTx fractions → the
    cancer-cell mesenchymal state beyond stromal mesenchyme. Structural markers (VIM/FN1/CDH2…), NOT
    TGF-β-signaling genes ⇒ non-circular for the TGF-β-target robustness check."""
    if "malemt" not in _CACHE:
        Y = _load()["Y"]
        mk = [g for g in ["VIM", "FN1", "CDH2", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "SPARC"] if g in Y.index]
        sub = Y.loc[mk]; z = sub.sub(sub.mean(1), axis=0).div(sub.std(1) + 1e-9, axis=0); score = z.mean(0)
        try:
            dv = _deconv()
            common = score.index.intersection(dv.index)
            lr = LinearRegression().fit(dv.loc[common].to_numpy(float), score.loc[common])
            resid = score.copy()
            resid.loc[common] = score.loc[common] - lr.predict(dv.loc[common].to_numpy(float))
            _CACHE["malemt"] = resid
        except Exception:
            _CACHE["malemt"] = score
    return _CACHE["malemt"].reindex(parts)


def _load() -> dict:
    if "X" not in _CACHE:                                       # robust to other loaders sharing _CACHE
        X = D.load_mirna_arms()                                 # arm × participant
        Y = D.load_rna()                                        # gene × participant
        _CACHE["X"] = X[~X.index.duplicated(keep="first")]
        _CACHE["Y"] = Y[~Y.index.duplicated(keep="first")]      # dedupe repeated gene symbols
        _CACHE["C"] = D.load_clinical_strata().drop_duplicates("participant").set_index("participant")
    return _CACHE



def _by_gene(df: pd.DataFrame, key: str) -> dict:
    """Memoized {gene -> sub-frame} index for the per-gene lookups in `assemble_gene`.

    Those lookups used to be `df["gene"] == gene` — an OBJECT-dtype (string) comparison scanning the WHOLE
    frame on every call. With `edge_weights` at ~10^5-10^6 rows and `assemble_gene` called once per gene, a
    1.5k-gene sweep re-scanned them 1.5k times: cProfile put `comp_method_OBJECT_ARRAY` at 63% of the
    per-gene cost (after the ledger parquet memo). Grouping ONCE makes each lookup a dict hit.
    Sub-frames are returned as-is and callers only read / re-index them (never mutate).
    """
    k = f"_bygene::{key}"
    if k not in _CACHE:
        _CACHE[k] = {g: sub for g, sub in df.groupby("gene", sort=False)}
    return _CACHE[k]


def assemble_gene(
    gene: str,
    *,
    he_only: bool = True,
    edges: Optional[pd.DataFrame] = None,
    w_prior_source: str = "evidence_score",
    deconv: bool = False,
    n_latent: int = 0,
    n_tf: int = 0,
    tf_activity: bool = False,
    orphans: bool = False,
    orphan_source: str = "targetscan",   # orphan candidacy: 'targetscan' (context++ sites) or 'kd' (scanMiR K_D)
    target_cn: bool = True,       # Design §B core (verified: neutral); the gene's own CN
    prolif: bool = True,          # Design §B core (verified: STRENGTHENS coupling); malignant proliferation
    mycaf: bool = False,          # finer stroma: myCAF-vs-iCAF axis beyond the lumped CAF fraction
    mal_emt: bool = False,        # within-malignant EMT state (beyond stromal mesenchyme)
    batch: bool = False,          # naive plate dummies over-parametrise → off; needs ComBat-style pre-adj
) -> Tuple[pd.Series, pd.DataFrame, pd.DataFrame, pd.Series]:
    """Return (Y, X, C, w_prior) aligned on common tumour participants.

    Y : Series[participant]            target gene log2(TPM+1)
    X : DataFrame[participant × arm]    regulator arm log2(RPM+1)
    C : DataFrame[participant × conf]   numeric confounders
    w_prior : Series[arm]               adaptive-penalty prior weight (larger = stronger edge)

    w_prior_source: 'evidence_score' (hand-weighted, un-deduped; the baseline) or 'ledger'
    (the PMID-deduped, method-centric fused weight from learned.evidence.ledger, Design §E).
    """
    d = _load()
    Xall, Yall, clin = d["X"], d["Y"], d["C"]
    if gene not in Yall.index:
        raise KeyError(f"{gene} not in RNA matrix")

    if edges is not None:
        ed = edges
    elif he_only:                                             # MIGRATED 2026-07-06: HE → POOLED-HE
        from mirna_hallmark.learned.evidence import ledger as LG   # miRTarBase-HE ∪ TarBase-v9 low-throughput
        ed = LG.pooled_he_edges()                                  # functional (~5,940, +~640 over miRTarBase-HE)
    else:
        if "raw_edges" not in _CACHE:                        # cached read-only view (avoids load_hallmark_edges' per-call full copy)
            _CACHE["raw_edges"] = D.load_hallmark_edges()
        ed = _CACHE["raw_edges"]
    if edges is None:                                          # shared frame -> use the memoized gene index
        _ek = "pooled_he" if he_only else "raw_edges"
        ed = _by_gene(ed, _ek).get(gene, ed.iloc[:0])
    else:
        ed = ed.loc[ed["gene"] == gene]
    ed = ed.loc[:, ["miRNA", "evidence_score"]].drop_duplicates("miRNA")
    _nm = _arm_name_map(Xall)                                           # `.N`/case/suffixless-guide edge->X resolution
    regs = [m for m in ed["miRNA"] if (_nm.get(m) or _nm.get(str(m).lower()))]
    if not regs:
        raise ValueError(f"no regulators of {gene} present in the arm matrix (he_only={he_only})")

    fallback = np.log1p(ed.set_index("miRNA")["evidence_score"].reindex(regs).astype(float))
    if w_prior_source in ("ledger", "ledger_mrna"):
        from mirna_hallmark.learned.evidence import ledger as LG
        from mirna_hallmark.learned.evidence.methods import CLASS_WEIGHT_MRNA
        lw = LG.edge_weights(weights=CLASS_WEIGHT_MRNA if w_prior_source == "ledger_mrna" else None)
        _lwk = f"lw_{w_prior_source}"
        _sub = _by_gene(lw, _lwk).get(gene, lw.iloc[:0])
        lw_g = _sub.set_index("arm")["ledger_weight"]
        w_prior = lw_g.reindex(regs).fillna(fallback)          # ledger weight; fall back to log-evidence
    elif w_prior_source == "scanmir":                          # biochemical affinity as the adaptive prior
        from mirna_hallmark.learned import kd as KD
        aff = KD.affinity()
        ag = aff.loc[aff["gene"] == gene].set_index("arm")["repression"]
        w_prior = (-ag.reindex(regs)).clip(lower=0.0)          # −repression = affinity magnitude (≥0)
        w_prior = w_prior.where(w_prior > 0).fillna(fallback)  # no site → fall back to log-evidence
    elif w_prior_source == "fused":                            # three-roles: ledger EXISTENCE + scanMiR MAGNITUDE
        from mirna_hallmark.learned.evidence import ledger as LG
        from mirna_hallmark.learned import kd as KD
        lw = LG.edge_weights()
        e = lw.loc[lw["gene"] == gene].set_index("arm")["ledger_weight"].reindex(regs).fillna(fallback)
        aff = KD.affinity()
        a = (-aff.loc[aff["gene"] == gene].set_index("arm")["repression"]).clip(lower=0.0).reindex(regs).fillna(0.0)
        ze = (e - e.mean()) / (e.std() or 1.0)                 # standardize each prior among THIS gene's regulators
        za = (a - a.mean()) / (a.std() or 1.0)
        # NOTE (2026-07-04): on the curated HE-only set this is NON-ADDITIVE — every HE regulator is already
        # existence-supported, so ledger + scanMiR are redundant magnitude proxies and any admixture breaks
        # the weak-coupling CDKN1A abundance gate (4/5) without reaching scanMiR-alone stability (0.90).
        # Symmetric blend kept un-tuned; the three-roles fusion's real slot is existence-gated ORPHAN magnitude.
        w_prior = pd.Series(np.exp(0.5 * (ze + za)), index=regs)   # equal-weight standardized blend (>0)
    else:
        w_prior = fallback
    w_prior = w_prior.clip(lower=1e-3)                          # adaptive-penalty weight (>0)

    if orphans:                                                 # add ORPHAN candidate edges (Design §E discovery)
        if orphan_source == "kd":                               # scanMiR K_D candidacy: biophysical, catches non-canonical seeds
            from mirna_hallmark.learned import kd as _KD
            aff = _KD.affinity()
            ka = aff[(aff["gene"] == gene) & aff["arm"].isin(Xall.index) & ~aff["arm"].isin(regs)
                     & (aff["repression"] < 0)].drop_duplicates("arm")
            ka = ka.nsmallest(80, "repression")                # TOP-80 strongest predicted duplexes (bounds count vs TS)
            if len(ka):
                ow = np.log1p((-ka.set_index("arm")["repression"]).clip(lower=0)).clip(lower=1e-3)
                regs = regs + list(ow.index)                   # orphans: K_D magnitude prior
                w_prior = pd.concat([w_prior, ow])
        else:                                                   # TargetScan context++ candidacy (default)
            ts = _targetscan_context()
            if ts is not None:
                _tsg = _by_gene(ts, "ts_context").get(gene, ts.iloc[:0])
                to = _tsg[_tsg["arm"].isin(Xall.index) & ~_tsg["arm"].isin(regs)]
                if len(to):
                    ow = np.log1p(to.set_index("arm")["ts_mag"]).clip(lower=1e-3)
                    regs = regs + list(ow.index)               # orphans: no curated evidence, context++ prior
                    w_prior = pd.concat([w_prior, ow])

    conf = [c for c in confounder_columns() if c in clin.columns]
    parts = clin.index.intersection(Yall.columns).intersection(Xall.columns)
    Y = Yall.loc[gene, parts].astype(float)
    # undetected arm in a sample → NaN in the RPM matrix → treat as 0 abundance (log2(RPM+1)=0),
    # not a reason to drop the participant (requiring all-arms-present decimates the cohort).
    regs_x = [_nm.get(m) or _nm.get(str(m).lower()) or m for m in regs]  # edge names -> actual X-index names
    X = Xall.loc[regs_x, parts].T.astype(float).fillna(0.0)     # participant × arm (X rows...
    X.columns = list(regs)                                      # ...renamed to canonical edge arm names)
    Cc = clin.loc[parts, conf].apply(pd.to_numeric, errors="coerce")
    if deconv:                                                  # add non-malignant composition (Design §B)
        dv = _deconv()
        if dv is not None:
            dvj = dv.reindex(Cc.index)
            dvj = dvj.fillna(dvj.median())                      # preserve n; impute missing fractions
            Cc = Cc.join(dvj)
    if n_latent > 0:                                            # global-PC transcription proxy (OVER-CORRECTS)
        Cc = Cc.join(_latent(n_latent).reindex(Cc.index))
    if n_tf > 0:                                                # TF-regulon proxy: strength-weighted expr PCs
        tfp = _tf_proxy(gene, Cc.index, k=n_tf)
        if tfp is not None:
            Cc = Cc.join(tfp)
    if tf_activity:                                             # TF-regulon proxy: signed activity (the 'proper' version)
        ta = _tf_activity(gene, Cc.index)
        if ta is not None:
            Cc = Cc.join(ta)
    if target_cn:                                              # Design §B: the gene's OWN copy number
        tcn = _target_cn([gene]).get(gene)
        if tcn is not None:
            Cc = Cc.join(pd.to_numeric(tcn.reindex(Cc.index), errors="coerce").rename("target_cn"))
    if prolif:                                                 # Design §B: malignant-compartment proliferation
        mp = _malignant_prolif(Cc.index)
        if mp is not None:
            Cc = Cc.join(mp.rename("mal_prolif"))
    if mycaf:                                                   # finer fibroblast subtype (myCAF vs iCAF)
        mc = _mycaf(Cc.index)
        if mc is not None:
            Cc = Cc.join(mc.rename("mycaf"))
    if mal_emt:                                                 # within-malignant EMT state
        me = _mal_emt(Cc.index)
        if me is not None:
            Cc = Cc.join(me.rename("mal_emt"))
    if batch:                                                  # TCGA analyte plate (miRNA+RNA), user-requested
        from mirna_hallmark import tcga_batch as TB
        Cc = TB.augment_cov(Cc, kind="plate_both")
    Cc = Cc.apply(lambda s: s.fillna(s.median()) if s.notna().any() else s.fillna(0.0))

    keep = Y.notna() & Cc[conf].notna().all(axis=1)
    return Y[keep], X.loc[keep], Cc.loc[keep], w_prior
