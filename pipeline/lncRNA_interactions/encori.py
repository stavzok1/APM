from __future__ import annotations

import hashlib
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
import pyarrow.dataset as ds

from pipeline.config import PATHS


ENCORI_BASE = "https://rnasysu.com/encori/api"


@dataclass(frozen=True)
class EncoriQuery:
    """
    Canonical query parameters for ENCORI endpoints.

    Notes:
    - ENCORI returns tab-delimited files (no JSON) for most endpoints.
    - Some queries can be huge if you use target=all and RBP/miRNA=all.
      Prefer restricting by either RBP list or by target gene list.
    """

    assembly: str = "hg38"
    gene_type: str = "lncRNA"  # ENCORI "geneType"
    cell_type: str = "all"

    # Filters
    clip_exp_num: int = 1
    degra_exp_num: int = 0
    pancancer_num: int = 0

    # miRNATarget-specific
    program_num: int = 1
    program: str = "TargetScan,miRanda,PITA"


def _cache_key(url: str) -> str:
    return hashlib.sha256(url.encode("utf-8")).hexdigest()[:16]


def _download_tsv(url: str, out_tsv: Path, timeout_s: int = 120) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "APM-pipeline/encori (academic use)"})
    with urllib.request.urlopen(req, timeout=timeout_s) as r:
        data = r.read()
    out_tsv.write_bytes(data)


def _encori_url(endpoint: str, params: dict[str, str | int]) -> str:
    q = urllib.parse.urlencode(params, doseq=True)
    return f"{ENCORI_BASE}/{endpoint}/?{q}"


def load_modelled_lncrna_gene_names(
    lncrna_gene_table: Path,
) -> list[str]:
    """
    Get lncRNA gene_name list from the pipeline-produced lncRNA multifeature table.

    Expected input is `PATHS.lncrnas_all_features` (e.g. `data/lncRNAs_genes_all_features.csv`).
    """
    df = pd.read_csv(lncrna_gene_table, low_memory=False)
    if "gene_name" not in df.columns:
        raise ValueError(f"Missing gene_name in {lncrna_gene_table}")
    names = (
        df["gene_name"]
        .dropna()
        .astype(str)
        .map(lambda s: s.strip())
        .loc[lambda s: s.ne("")]
        .drop_duplicates()
        .tolist()
    )
    return names


def build_encori_lncrna_target_list(
    *,
    panel_lncrnas: Iterable[str],
    lncrna_proximity_pairs_csv: Path,
    primary_genes: Iterable[str],
    n_extra_close_lncrnas: int = 20,
) -> list[str]:
    """
    Deterministic lncRNA target list for ENCORI miRNA queries:
    - always include the lncRNAs explicitly in the panel
    - add N additional lncRNAs that are closest (by min_distance_bp) to any of the frozen 66 primary genes
    """
    panel = [str(x).strip() for x in panel_lncrnas if str(x).strip()]

    pairs = pd.read_csv(lncrna_proximity_pairs_csv, low_memory=False)
    needed = {"gene_name", "lncRNA_name", "min_distance_bp"}
    if not needed.issubset(set(pairs.columns)):
        raise ValueError(f"{lncrna_proximity_pairs_csv} missing columns {sorted(needed - set(pairs.columns))}")

    prim = {str(x).strip() for x in primary_genes if str(x).strip()}
    pairs = pairs[pairs["gene_name"].astype(str).isin(prim)].copy()
    if pairs.empty:
        return panel

    pairs["min_distance_bp"] = pd.to_numeric(pairs["min_distance_bp"], errors="coerce")
    best = (
        pairs.dropna(subset=["min_distance_bp"])
        .sort_values(["min_distance_bp", "lncRNA_name", "gene_name"], ascending=[True, True, True])
        .groupby("lncRNA_name", sort=False)["min_distance_bp"]
        .min()
        .sort_values()
    )

    extra = [ln for ln in best.index.astype(str).tolist() if ln not in set(panel)][: int(n_extra_close_lncrnas)]
    return panel + extra


def fetch_encori_rbp_targets(
    *,
    rbps: Iterable[str],
    target: str = "all",
    query: EncoriQuery = EncoriQuery(),
    cache_dir: Path,
) -> pd.DataFrame:
    """
    Fetch RBP-Target rows for a set of RBPs.

    This implements the *practical* strategy for "all targets, all cell types":
    - iterate an RBP list and request target=all, cellType=all
    - keep geneType=lncRNA so output stays focused and manageable
    """
    frames: list[pd.DataFrame] = []
    for rbp in rbps:
        params = {
            "assembly": query.assembly,
            "geneType": query.gene_type,
            "RBP": rbp,
            "clipExpNum": int(query.clip_exp_num),
            "pancancerNum": int(query.pancancer_num),
            "target": target,
            "cellType": query.cell_type,
        }
        url = _encori_url("RBPTarget", params)
        out_tsv = cache_dir / "encori" / "RBPTarget" / f"{rbp}.{_cache_key(url)}.tsv"
        if not out_tsv.exists():
            _download_tsv(url, out_tsv)
        # ENCORI responses often include leading "#please cite" comment lines, and may include
        # one-line error strings for invalid RBPs; parse robustly.
        df = pd.read_csv(
            out_tsv,
            sep="\t",
            low_memory=False,
            comment="#",
            on_bad_lines="skip",
        )
        # ENCORI sometimes returns a one-line error string that slips through parsing; keep only real rows.
        if "geneID" in df.columns and "geneName" in df.columns:
            # Some ENCORI releases omit geneID for lncRNA targets; geneName is still populated.
            df = df[df["geneID"].notna() | df["geneName"].notna()].copy()
        elif "geneID" in df.columns:
            df = df[df["geneID"].notna()].copy()
        elif "geneName" in df.columns:
            df = df[df["geneName"].notna()].copy()
        df["__encori_url"] = url
        if not df.empty:
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def fetch_encori_rbp_targets_for_lncrnas(
    *,
    rbps: Iterable[str],
    lncrna_targets: Iterable[str],
    query: EncoriQuery = EncoriQuery(),
    cache_dir: Path,
) -> pd.DataFrame:
    """
    Fetch ENCORI RBPTarget rows for (RBP × lncRNA target) pairs.

    This is the practical alternative to `RBP=all&target=all`, which is typically too large to download.
    """
    frames: list[pd.DataFrame] = []
    for rbp in rbps:
        for tgt in lncrna_targets:
            t = str(tgt).strip()
            if not t:
                continue
            df = fetch_encori_rbp_targets(rbps=[rbp], target=t, query=query, cache_dir=cache_dir)
            if df.empty:
                continue
            df = df.copy()
            df["__encori_lncrna_target_query"] = t
            df["__encori_rbp_query"] = str(rbp)
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def fetch_encori_rbp_targets_for_lncrnas_with_fallback(
    *,
    rbps: Iterable[str],
    lncrna_targets: Iterable[str],
    query: EncoriQuery = EncoriQuery(),
    cache_dir: Path,
    gencode_parquet: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Like `fetch_encori_rbp_targets_for_lncrnas`, but for each lncRNA target it tries:
    - target=<symbol or ENSG as provided>
    - if empty and target is a symbol: try GENCODE gene_id (versioned), then unversioned ENSG

    Returns:
    - combined rows
    - diagnostics per (rbp, query_symbol)
    """
    gencode_parquet = gencode_parquet or Path(PATHS.gencode_gtf_full_pq)
    frames: list[pd.DataFrame] = []
    rows_diag: list[dict[str, object]] = []

    for rbp in rbps:
        r = str(rbp).strip()
        if not r:
            continue
        for tgt in lncrna_targets:
            s = str(tgt).strip()
            if not s:
                continue

            df_sym = fetch_encori_rbp_targets(rbps=[r], target=s, query=query, cache_dir=cache_dir)
            n_sym = int(len(df_sym))

            gid_ver: Optional[str] = None
            gid_unver: Optional[str] = None
            df_gid = pd.DataFrame()
            n_gid = 0
            used = "symbol"

            if n_sym == 0 and (not s.startswith("ENSG")):
                gid_ver = _gencode_gene_id_for_symbol(s, gencode_parquet)
                if gid_ver:
                    df_gid = fetch_encori_rbp_targets(rbps=[r], target=gid_ver, query=query, cache_dir=cache_dir)
                    n_gid = int(len(df_gid))
                    used = "gene_id" if n_gid else used
                    if n_gid == 0 and "." in gid_ver:
                        gid_unver = gid_ver.split(".", 1)[0]
                        df_gid2 = fetch_encori_rbp_targets(rbps=[r], target=gid_unver, query=query, cache_dir=cache_dir)
                        n_gid2 = int(len(df_gid2))
                        if n_gid2:
                            df_gid = df_gid2
                            n_gid = n_gid2
                            used = "gene_id"

            rows_diag.append(
                {
                    "rbp_query": r,
                    "query_symbol": s,
                    "gencode_gene_id": gid_ver,
                    "gencode_gene_id_unversioned": gid_unver,
                    "n_rows_symbol_query": n_sym,
                    "n_rows_gene_id_query": n_gid,
                    "n_rows_total": n_sym + n_gid,
                    "used": used if (n_sym + n_gid) else "none",
                }
            )

            if n_sym:
                out = df_sym.copy()
                out["__encori_lncrna_target_query"] = s
                out["__encori_lncrna_target_query_type"] = "symbol_or_input"
                out["__encori_rbp_query"] = r
                frames.append(out)
            if n_gid:
                out = df_gid.copy()
                out["__encori_lncrna_target_query"] = s
                out["__encori_lncrna_target_query_type"] = "gene_id"
                out["__encori_rbp_query"] = r
                frames.append(out)

    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    diag = pd.DataFrame(rows_diag)
    return df, diag


def fetch_encori_miRNA_targets(
    *,
    targets: Iterable[str],
    miRNA: str = "all",
    query: EncoriQuery = EncoriQuery(),
    cache_dir: Path,
    clip_exp_num: Optional[int] = None,
) -> pd.DataFrame:
    """
    Fetch miRNA-Target rows restricted to a target gene list (recommended for lncRNAs).
    """
    frames: list[pd.DataFrame] = []
    clip_exp = int(query.clip_exp_num if clip_exp_num is None else clip_exp_num)
    for t in targets:
        params = {
            "assembly": query.assembly,
            "geneType": query.gene_type,
            "miRNA": miRNA,
            "clipExpNum": clip_exp,
            "degraExpNum": int(query.degra_exp_num),
            "pancancerNum": int(query.pancancer_num),
            "programNum": int(query.program_num),
            "program": query.program,
            "target": t,
            "cellType": query.cell_type,
        }
        url = _encori_url("miRNATarget", params)
        out_tsv = cache_dir / "encori" / "miRNATarget" / f"{t}.{_cache_key(url)}.tsv"
        if not out_tsv.exists():
            _download_tsv(url, out_tsv)
        df = pd.read_csv(
            out_tsv,
            sep="\t",
            low_memory=False,
            comment="#",
            on_bad_lines="skip",
        )
        # miRNATarget rows should always have a target gene symbol; keep rows even if geneID is missing.
        if "geneName" in df.columns:
            df = df[df["geneName"].notna()].copy()
        df["__encori_url"] = url
        if not df.empty:
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _gencode_gene_id_for_symbol(symbol: str, gencode_parquet: Path) -> Optional[str]:
    """
    Resolve a human gene symbol to a GENCODE gene_id (versioned ENSG…), using gene rows only.
    """
    sym = str(symbol).strip()
    if not sym:
        return None

    d = ds.dataset(str(gencode_parquet), format="parquet")
    cols = set(d.schema.names)
    if "gene_name" not in cols or "gene_id" not in cols or "feature" not in cols:
        return None

    table = d.to_table(
        columns=["gene_name", "gene_id", "gene_type", "feature"],
        filter=(ds.field("feature") == "gene") & (ds.field("gene_name") == sym),
    )
    df = table.to_pandas()
    if df.empty:
        return None

    # Prefer lncRNA rows when present (important for ambiguous symbols / pseudo-genes).
    if "gene_type" in df.columns:
        lnc = df[df["gene_type"].astype(str).str.contains("lncRNA", na=False)]
        if not lnc.empty:
            df = lnc

    gid = str(df.iloc[0]["gene_id"])
    return gid


def fetch_encori_miRNA_targets_with_symbol_fallback(
    *,
    targets: Iterable[str],
    miRNA: str = "all",
    query: EncoriQuery = EncoriQuery(),
    cache_dir: Path,
    clip_exp_num: Optional[int] = None,
    gencode_parquet: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fetch miRNA-target rows for a list of targets, trying **gene symbol first**, then falling back to
    **GENCODE gene_id** if the symbol query yields no rows.

    Returns:
    - combined dataframe (concatenated)
    - per-target diagnostics (symbol, gene_id tried, n_rows_symbol, n_rows_gene_id)
    """
    gencode_parquet = gencode_parquet or Path(PATHS.gencode_gtf_full_pq)

    rows_diag: list[dict[str, object]] = []
    frames: list[pd.DataFrame] = []

    for sym in targets:
        s = str(sym).strip()
        if not s:
            continue

        df_sym = fetch_encori_miRNA_targets(
            targets=[s],
            miRNA=miRNA,
            query=query,
            cache_dir=cache_dir,
            clip_exp_num=clip_exp_num,
        )
        n_sym = int(len(df_sym))

        df_gid = pd.DataFrame()
        n_gid = 0
        gid: Optional[str] = None
        gid_ver: Optional[str] = None
        gid_unver: Optional[str] = None
        if n_sym == 0 and (not s.startswith("ENSG")):
            gid_ver = _gencode_gene_id_for_symbol(s, gencode_parquet)
            gid = gid_ver
            if gid_ver:
                df_gid = fetch_encori_miRNA_targets(
                    targets=[gid_ver],
                    miRNA=miRNA,
                    query=query,
                    cache_dir=cache_dir,
                    clip_exp_num=clip_exp_num,
                )
                n_gid = int(len(df_gid))
                if n_gid == 0 and "." in gid_ver:
                    gid_unver = gid_ver.split(".", 1)[0]
                    df_gid2 = fetch_encori_miRNA_targets(
                        targets=[gid_unver],
                        miRNA=miRNA,
                        query=query,
                        cache_dir=cache_dir,
                        clip_exp_num=clip_exp_num,
                    )
                    n_gid2 = int(len(df_gid2))
                    if n_gid2:
                        df_gid = df_gid2
                        n_gid = n_gid2
                        gid = gid_unver

        rows_diag.append(
            {
                "query_symbol": s,
                "gencode_gene_id": gid_ver,
                "gencode_gene_id_unversioned": gid_unver,
                "n_rows_symbol_query": n_sym,
                "n_rows_gene_id_query": n_gid,
                "n_rows_total": n_sym + n_gid,
            }
        )

        if n_sym:
            df_sym = df_sym.copy()
            df_sym["__encori_target_query"] = s
            df_sym["__encori_target_query_type"] = "symbol"
            frames.append(df_sym)
        if n_gid:
            df_gid = df_gid.copy()
            df_gid["__encori_target_query"] = s
            df_gid["__encori_target_query_type"] = "gene_id"
            frames.append(df_gid)

    diag = pd.DataFrame(rows_diag)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return out, diag

