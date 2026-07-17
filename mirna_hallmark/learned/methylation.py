"""miRNA-locus promoter methylation gate — the mechanism confirmation for `lost_specialist`.

Answers "is this miRNA's promoter epigenetically silenced (hyper-methylated) in tumour?" **directly**, via
CpG probes that overlap the miRNA promoter — NOT via the cCRE aggregation the pipeline used as a shortcut.

    promoter window   = strand-aware [TSS − UP, TSS + DOWN] around the pri-miRNA hairpin TSS
                        (hairpin span = min/max over its mature arms in mirna_mature_loci, grouped by pre_gene_id).
    promoter probes    = probe_annotations rows whose (chrom,center) fall in that window (49,869 promoter/CGI-
                        curated probes; the regulatory-relevant subset — miRNA silencing CpG islands live here).
    Δβ (hyper-meth)    = mean β(TCGA-BRCA tumour) − mean β(normal-adjacent), group means over the promoter probes
                        (raw_samples/Tumor 802 vs raw_samples/Normal 99 — no barcode map needed for the contrast).
    gate               = Δβ ≥ DELTA_BETA_HYPER (default 0.15) AND tumour mean β ≥ BETA_METHYLATED (0.3) — a
                        real hyper-methylated promoter, consistent with epigenetic silencing of the arm.

This is the LOSS-mechanism leg of the silenced-specialist detector (structural + baseline-active + tumour-
silenced + **hyper-methylated promoter**). Reuses parent methylation data (data/Methylation/raw_samples).

CLI: `python -m mirna_hallmark.learned.methylation miR-200b-3p miR-34c-5p [--full]`  (--full = all 901 samples)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from mirna_hallmark import config as C

UP, DOWN = 2000, 500                 # promoter window upstream / downstream of the hairpin TSS (bp)
DELTA_BETA_HYPER = 0.15              # tumour−normal β gain to call hyper-methylation
BETA_METHYLATED = 0.20              # tumour promoter must actually be methylated (mean β diluted by shore probes)
_MDIR = Path("data/Methylation/raw_samples")
_BETA_CACHE = Path("mirna_hallmark/output/learned/mirna_promoter_betas.parquet")   # probe × sample β (built once)
_CACHE: dict = {}


def _norm(name: str) -> str:
    """model arm name → canonical key: drop hsa-, lower, strip the trailing copy suffix (miR-124-3p.1/.2/.3
    are identical mature arms from distinct precursors → one key, all their promoters collected)."""
    return re.sub(r"\.\d+$", "", str(name).lower().replace("hsa-", "")).strip()


def _hairpin_loci() -> pd.DataFrame:
    """Per-precursor (hairpin) genomic span, from mature-arm loci grouped by pre_gene_id. Columns:
    chrom,start,end,strand,tss + one row per mature arm key (norm mirbase id) → precursor, so an arm maps in."""
    if "hairpin" in _CACHE:
        return _CACHE["hairpin"]
    m = pd.read_csv(C.MIRNA_MATURE_LOCI)
    m = m.dropna(subset=["pre_gene_id", "chrom", "start", "end"])
    grp = m.groupby("pre_gene_id").agg(chrom=("chrom", "first"), start=("start", "min"),
                                       end=("end", "max"), strand=("strand", "first")).reset_index()
    grp["tss"] = np.where(grp["strand"] == "+", grp["start"], grp["end"])
    # arm key → precursor row (mirbase_mature_id is 'hsa-miR-200b-3p')
    key = m[["pre_gene_id", "mirbase_mature_id"]].copy()
    key["arm_key"] = key["mirbase_mature_id"].map(_norm)
    out = key.merge(grp, on="pre_gene_id", how="left").dropna(subset=["chrom"])
    _CACHE["hairpin"] = out
    return out


def _probe_ref() -> pd.DataFrame:
    """FULL sesame HM450 manifest (485,577 probes, hg38) — NOT the gene-centric 49,869-probe curated subset,
    which has no coverage below chr1:38 Mb and so misses miRNA loci entirely (that gene-promoter blind spot
    is exactly why the pipeline fell back to cCRE). Columns: probeID, chrom, center, in_CGI, CGI_context."""
    if "probes" in _CACHE:
        return _CACHE["probes"]
    from pipeline.config import PATHS
    df = pd.read_csv(PATHS.methylation_probe_reference, sep="\t",
                     usecols=lambda c: c in ("probeID", "CpG_chrm", "CpG_beg", "CpG_end", "CGI", "CGIposition"))
    df = df.rename(columns={"CpG_chrm": "chrom", "CpG_beg": "start", "CpG_end": "end"})
    df["chrom"] = df["chrom"].astype(str)
    df["chrom"] = df["chrom"].where(df["chrom"].str.startswith("chr"), "chr" + df["chrom"])
    df["center"] = ((df["start"] + df["end"]) // 2)
    df = df.dropna(subset=["center"])
    df["center"] = df["center"].astype("int64")
    df["in_CGI"] = df.get("CGI", pd.Series(index=df.index)).notna() & (df.get("CGI", "") != "")
    df["CGI_context"] = df.get("CGIposition", pd.Series(index=df.index)).fillna("OpenSea").replace("", "OpenSea")
    _CACHE["probes"] = df[["probeID", "chrom", "center", "in_CGI", "CGI_context"]]
    return _CACHE["probes"]


def promoter_probes(arms, *, up: int = UP, down: int = DOWN) -> pd.DataFrame:
    """CpG probes overlapping the promoter window of the given miRNA arm(s). DIRECT probe→promoter overlap
    (probe center within [TSS−up, TSS+down] on the hairpin strand); no cCRE intermediary. Columns:
    probeID, arm_key, chrom, center, dist_to_tss, in_CGI, CGI_context."""
    hp = _hairpin_loci()
    keys = {_norm(a) for a in ([arms] if isinstance(arms, str) else arms)}
    loci = hp[hp["arm_key"].isin(keys)].drop_duplicates("pre_gene_id")
    if loci.empty:
        return pd.DataFrame()
    pr = _probe_ref()
    rows = []
    for r in loci.itertuples():
        lo = r.tss - up if r.strand == "+" else r.tss - down
        hi = r.tss + down if r.strand == "+" else r.tss + up
        hit = pr[(pr["chrom"] == r.chrom) & (pr["center"] >= lo) & (pr["center"] <= hi)].copy()
        if hit.empty:
            continue
        hit["arm_key"] = r.arm_key
        hit["dist_to_tss"] = (hit["center"] - r.tss) * (1 if r.strand == "+" else -1)
        rows.append(hit[["probeID", "arm_key", "chrom", "center", "dist_to_tss", "in_CGI", "CGI_context"]])
    return pd.concat(rows, ignore_index=True).drop_duplicates("probeID") if rows else pd.DataFrame()


def all_promoter_probes() -> pd.DataFrame:
    """Union of promoter probes across the WHOLE miRNA arm universe (one map, cached) — the fixed probe set
    to extract once so per-gene calls never re-read raw files."""
    if "allpp" in _CACHE:
        return _CACHE["allpp"]
    hp = _hairpin_loci()
    pp = promoter_probes(hp["arm_key"].dropna().unique().tolist())
    _CACHE["allpp"] = pp
    return pp


def _read_group(probes: list[str], group: str, *, limit: int | None = None) -> pd.DataFrame:
    """probe × sample β for a raw-sample group (Tumor|Normal), subset to `probes`. Reads the 2-col level3
    beta txt files, subsetting to the probe union so each file is touched once."""
    files = sorted((_MDIR / group).glob("*.level3betas.txt"))   # exclude stray gdc_manifest.txt etc.
    if limit:
        files = files[:limit]
    pidx = pd.Index(probes)
    cols = {}
    for f in files:
        s = pd.read_csv(f, sep="\t", header=None, names=["probe", "beta"], usecols=[0, 1],
                        dtype={"probe": str}, na_values=["NA"]).set_index("probe")["beta"]
        cols[f.stem] = pd.to_numeric(s, errors="coerce").reindex(pidx)
    return pd.DataFrame(cols, index=pidx)


def build_beta_cache(*, force: bool = False) -> Path:
    """One-time hoisted read: extract the promoter-probe union β across ALL 802 tumour + 99 normal raw files,
    write probe × sample parquet (col MultiIndex-free; a companion _group row marks Tumor/Normal). ~7–10 min;
    run in background. After this, every `locus_methylation` call is a cache slice (no file I/O)."""
    if _BETA_CACHE.exists() and not force:
        print(f"[methylation] beta cache exists: {_BETA_CACHE}")
        return _BETA_CACHE
    probes = all_promoter_probes()["probeID"].tolist()
    print(f"[methylation] building beta cache for {len(probes)} promoter probes × (802 tum + 99 nor) ...")
    tum = _read_group(probes, "Tumor"); nor = _read_group(probes, "Normal")
    tum.columns = [f"T::{c}" for c in tum.columns]; nor.columns = [f"N::{c}" for c in nor.columns]
    mat = pd.concat([tum, nor], axis=1)
    _BETA_CACHE.parent.mkdir(parents=True, exist_ok=True)
    mat.to_parquet(_BETA_CACHE)
    print(f"[methylation] cached {mat.shape[0]} probes × {mat.shape[1]} samples → {_BETA_CACHE}")
    return _BETA_CACHE


def _beta_cache() -> pd.DataFrame | None:
    if "betamat" in _CACHE:
        return _CACHE["betamat"]
    if not _BETA_CACHE.exists():
        return None
    _CACHE["betamat"] = pd.read_parquet(_BETA_CACHE)
    return _CACHE["betamat"]


def locus_methylation(arms, *, limit: int | None = None) -> dict:
    """Tumour-vs-normal promoter methylation for a miRNA arm/family. Returns dict with the promoter probes,
    per-group mean β, Δβ, and the hyper-methylation gate. Uses the cached probe×sample matrix when present
    (instant); else reads raw files directly (`limit` caps files per group, for profiling)."""
    pp = promoter_probes(arms)
    if pp.empty:
        return {"n_probes": 0, "hypermethylated": False, "reason": "no promoter probes for arm(s)"}
    probes = pp["probeID"].tolist()
    mat = _beta_cache()
    if mat is not None and limit is None:
        sub = mat.reindex(probes)
        tcols = [c for c in sub.columns if c.startswith("T::")]; ncols = [c for c in sub.columns if c.startswith("N::")]
        tum, nor = sub[tcols], sub[ncols]
    else:
        tum = _read_group(probes, "Tumor", limit=limit); nor = _read_group(probes, "Normal", limit=limit)
    # per-sample mean over promoter probes, then group mean (probe-robust to missingness)
    tum_s, nor_s = tum.mean(axis=0, skipna=True), nor.mean(axis=0, skipna=True)
    tmean, nmean = float(tum_s.mean()), float(nor_s.mean())
    dbeta = tmean - nmean
    # SYMMETRIC gate. hyper = tumour GAINED promoter methylation → arm SILENCED → repression LOST (target ↑).
    #                 hypo  = tumour LOST a normal-present promoter methylation → arm DE-repressed → repression GAINED (target ↓).
    hyper = bool(dbeta >= DELTA_BETA_HYPER and tmean >= BETA_METHYLATED)
    hypo = bool(dbeta <= -DELTA_BETA_HYPER and nmean >= BETA_METHYLATED)
    direction = "hyper" if hyper else ("hypo" if hypo else "none")
    return {"n_probes": len(probes), "n_cgi": int(pp["in_CGI"].sum()),
            "tumour_beta": round(tmean, 3), "normal_beta": round(nmean, 3), "delta_beta": round(dbeta, 3),
            "n_tumour": int(tum.shape[1]), "n_normal": int(nor.shape[1]),
            "hypermethylated": hyper, "hypomethylated": hypo, "direction": direction, "probes": pp}


def run(arms, *, limit: int | None = None) -> None:
    for a in ([arms] if isinstance(arms, str) else arms):
        r = locus_methylation(a, limit=limit)
        if r["n_probes"] == 0:
            print(f"\n{a}: {r['reason']}"); continue
        flag = {"hyper": "⚑ HYPER-methylated → repression LOST",
                "hypo": "⚑ HYPO-methylated → repression GAINED", "none": "—"}[r["direction"]]
        print(f"\n=== {a} — promoter methylation ({r['n_probes']} probes, {r['n_cgi']} in CGI; "
              f"n_tum={r['n_tumour']} n_nor={r['n_normal']}) ===")
        print(f"  tumour β={r['tumour_beta']}  normal β={r['normal_beta']}  Δβ={r['delta_beta']:+.3f}   {flag}")


if __name__ == "__main__":
    if "--build" in sys.argv:
        build_beta_cache(force="--force" in sys.argv)
    else:
        lim = None if ("--full" in sys.argv or _BETA_CACHE.exists()) else 40
        a = [x for x in sys.argv[1:] if not x.startswith("-")] or ["miR-200b-3p", "miR-34c-5p", "miR-124-3p", "miR-21-5p"]
        src = "cache" if _BETA_CACHE.exists() else (f"raw (limit={lim}/group)" if lim else "raw FULL")
        print(f"[methylation] source: {src}")
        run(a, limit=lim)
