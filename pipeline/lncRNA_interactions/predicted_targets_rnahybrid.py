from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from pipeline.config import PATHS, PRIMARY_GENES, TIER1_LNCRNA_GENES
from pipeline.lncRNA_interactions.encori import build_encori_lncrna_target_list


def _reference_fasta_uses_chr_prefix(*, genome_fa: Path) -> bool:
    """
    Infer whether a reference FASTA uses UCSC-style ``chr1`` contig names.

    Ensembl/NCBI primary assemblies often use ``1``/``MT`` without a ``chr`` prefix, while many
    GENCODE-derived BEDs use ``chr*``. ``bedtools getfasta`` requires exact string matches.
    """
    with genome_fa.open("rb") as f:
        chunk = f.read(4096)
    text = chunk.decode("utf-8", errors="replace")
    for line in text.splitlines():
        if line.startswith(">"):
            return line.startswith(">chr")
    raise ValueError(f"Could not read a FASTA header from {genome_fa}")


def _normalize_bed_chrom_for_reference(chrom: str, *, ref_has_chr_prefix: bool) -> str:
    c = str(chrom)
    if ref_has_chr_prefix:
        return c if c.startswith("chr") else f"chr{c}"

    # Reference without ``chr`` prefix (common for Ensembl FASTAs).
    if c == "chrM":
        return "MT"
    if c.startswith("chr"):
        return c[len("chr") :]
    return c


@dataclass(frozen=True)
class RnaHybridPaths:
    out_dir: Path = Path(PATHS.working_dir) / "lncRNA_interactions" / "predicted_targets" / "rnahybrid"
    exon_bed: Path = out_dir / "selected_lncrnas.exons.bed"
    exons_fa: Path = out_dir / "selected_lncrnas.exons.fa"
    spliced_fa: Path = out_dir / "selected_lncrnas.spliced_exons.fa"
    spliced_fa_trunc: Path = out_dir / "selected_lncrnas.spliced_exons.3prime_trunc.fa"
    raw_txt: Path = out_dir / "rnahybrid.raw.txt"
    hits_parquet: Path = out_dir / "rnahybrid.hits.parquet"
    mirna_subset_fa: Path = out_dir / "mirnas.subset.fa"


def _write_lncrna_exon_bed(
    *,
    lncrnas: Iterable[str],
    gencode_parquet: Path,
    out_bed: Path,
    genome_fasta: Path,
) -> pd.DataFrame:
    """
    Build a BED6 file of exon intervals for selected lncRNAs.

    Name field is `gene_name|gene_id|exon_start-exon_end|strand` so downstream FASTA headers preserve context.
    """
    wanted = {str(x).strip() for x in lncrnas if str(x).strip()}
    if not wanted:
        raise ValueError("No lncRNAs provided.")

    d = ds.dataset(str(gencode_parquet), format="parquet")
    cols = set(d.schema.names)
    chrom_col = "chrom" if "chrom" in cols else "seqname" if "seqname" in cols else None
    if chrom_col is None:
        raise ValueError(f"Cannot find chrom/seqname column in {gencode_parquet}")

    # Exon rows; filter to desired gene_name. (We load all exon rows then filter in pandas to avoid large IN filters.)
    cols = set(d.schema.names)
    base_cols = [chrom_col, "start", "end", "strand", "gene_id", "feature"]
    if "gene_name" in cols:
        base_cols.insert(4, "gene_name")
    elif "gene" in cols:
        base_cols.insert(4, "gene")
    else:
        raise ValueError(f"Missing gene symbol column (expected gene_name/gene) in {gencode_parquet}")
    if "transcript_id" in cols:
        base_cols.append("transcript_id")

    table = d.to_table(columns=base_cols, filter=(ds.field("feature") == "exon"))
    ex = table.to_pandas()
    if "gene_name" in ex.columns:
        sym = ex["gene_name"].astype(str)
    else:
        sym = ex["gene"].astype(str)
    ex = ex[sym.isin(wanted)].copy()
    ex = ex.rename(columns={chrom_col: "chrom"})
    if "gene_name" not in ex.columns and "gene" in ex.columns:
        ex["gene_name"] = ex["gene"].astype(str)
    ex["start"] = pd.to_numeric(ex["start"], errors="coerce")
    ex["end"] = pd.to_numeric(ex["end"], errors="coerce")
    ex = ex.dropna(subset=["start", "end"])
    ex["start"] = ex["start"].astype(int)
    ex["end"] = ex["end"].astype(int)
    ex["strand"] = ex["strand"].astype(str).replace({"nan": "+"})

    ref_has_chr = _reference_fasta_uses_chr_prefix(genome_fa=genome_fasta)
    ex["chrom"] = ex["chrom"].map(lambda x: _normalize_bed_chrom_for_reference(str(x), ref_has_chr_prefix=ref_has_chr))

    # Sort exons in transcription order for later splicing.
    # Avoid ``groupby(...).apply`` (pandas version differences around grouping columns).
    if "transcript_id" in ex.columns:
        span = (
            ex.groupby(["gene_name", "transcript_id"], sort=False)["end"]
            .max()
            - ex.groupby(["gene_name", "transcript_id"], sort=False)["start"].min()
        )
        span = span.rename("tx_span").reset_index()
        pick = span.sort_values(["gene_name", "tx_span"], ascending=[True, False]).drop_duplicates("gene_name")
        ex = ex.merge(pick[["gene_name", "transcript_id"]], on=["gene_name", "transcript_id"], how="inner")

    strand_minus = ex["strand"].astype(str).eq("-")
    tx_key = np.where(strand_minus, -ex["start"].to_numpy(), ex["start"].to_numpy())
    ex = ex.assign(_tx_key=tx_key)
    sort_cols = ["gene_name"] + (["transcript_id"] if "transcript_id" in ex.columns else []) + ["_tx_key", "start"]
    ex = ex.sort_values(sort_cols, kind="mergesort").drop(columns=["_tx_key"]).reset_index(drop=True)
    ex["name"] = ex.apply(
        lambda r: f"{r['gene_name']}|{r.get('gene_id','')}|{int(r['start'])}-{int(r['end'])}|{r['strand']}",
        axis=1,
    )

    bed = pd.DataFrame(
        {
            "chrom": ex["chrom"].astype(str),
            "chromStart": ex["start"].astype(int),
            "chromEnd": ex["end"].astype(int),
            "name": ex["name"].astype(str),
            "score": 0,
            "strand": ex["strand"].astype(str),
        }
    )
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    bed.to_csv(out_bed, sep="\t", header=False, index=False)
    return ex


def _bedtools_getfasta(*, genome_fa: Path, bed: Path, out_fa: Path) -> None:
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["bedtools", "getfasta", "-fi", str(genome_fa), "-bed", str(bed), "-s", "-name"]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    out_fa.write_text(p.stdout, encoding="utf-8")
    stderr_path = out_fa.with_suffix(out_fa.suffix + ".bedtools.stderr.txt")
    stderr_path.write_text(p.stderr or "", encoding="utf-8")
    if p.returncode != 0:
        raise RuntimeError(f"bedtools getfasta failed (exit {p.returncode}). stderr: {stderr_path}")
    if out_fa.stat().st_size == 0:
        raise RuntimeError(
            "bedtools getfasta produced an empty FASTA. "
            f"This is commonly a chr-prefix mismatch between the BED chrom column and `{genome_fa}`. "
            f"See stderr log: {stderr_path}"
        )


def _splice_exon_fasta(exons_fa: Path, out_spliced_fa: Path) -> None:
    """
    Convert exon-level FASTA into one spliced sequence per lncRNA gene_name.
    """
    out_spliced_fa.parent.mkdir(parents=True, exist_ok=True)

    cur_name: Optional[str] = None
    cur_seq: list[str] = []
    per_gene: dict[str, list[str]] = {}

    with exons_fa.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_name is not None:
                    gene = cur_name.split("|", 1)[0]
                    per_gene.setdefault(gene, []).append("".join(cur_seq))
                cur_name = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line.upper())
        if cur_name is not None:
            gene = cur_name.split("|", 1)[0]
            per_gene.setdefault(gene, []).append("".join(cur_seq))

    with out_spliced_fa.open("w", encoding="utf-8") as w:
        for gene, exon_seqs in per_gene.items():
            seq = "".join(exon_seqs)
            if not seq:
                continue
            w.write(f">{gene}\n")
            for i in range(0, len(seq), 80):
                w.write(seq[i : i + 80] + "\n")


def _read_fasta_dict(path: Path) -> dict[str, str]:
    """
    Minimal FASTA reader: returns {header_id: sequence} where header_id is the first token after '>'.
    """
    out: dict[str, str] = {}
    cur: Optional[str] = None
    buf: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
        if cur is not None:
            out[cur] = "".join(buf).upper()
    return out


def _truncate_fasta_sequences_in_place(*, in_fa: Path, out_fa: Path, max_bp: int, from_end: bool = True) -> dict[str, int]:
    """
    Write a new FASTA with sequences truncated to ``max_bp`` nucleotides.

    Returns per-sequence lengths after truncation.
    """
    if max_bp <= 0:
        raise ValueError("max_bp must be positive.")

    seqs = _read_fasta_dict(in_fa)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    lens: dict[str, int] = {}
    with out_fa.open("w", encoding="utf-8") as w:
        for sid, s in seqs.items():
            if from_end:
                t = s[-max_bp:] if len(s) > max_bp else s
            else:
                t = s[:max_bp] if len(s) > max_bp else s
            lens[sid] = len(t)
            w.write(f">{sid}\n")
            for i in range(0, len(t), 80):
                w.write(t[i : i + 80] + "\n")
    return lens


def _max_fasta_seq_len(fa: Path) -> int:
    lens = _read_fasta_dict(fa)
    return max((len(s) for s in lens.values()), default=0)


def write_mirna_subset_fasta_from_encori(
    *,
    encori_mirna_parquet: Path,
    mature_mirna_fasta: Path,
    out_fa: Path,
) -> pd.DataFrame:
    """
    Build a small miRNA FASTA containing only miRNAs observed in ENCORI miRNA-target output.
    """
    m = pd.read_parquet(encori_mirna_parquet, columns=["miRNAname"])
    names = sorted({str(x).strip() for x in m["miRNAname"].dropna().astype(str).tolist() if str(x).strip()})
    mature = _read_fasta_dict(mature_mirna_fasta)

    rows = []
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w", encoding="utf-8") as w:
        for n in names:
            seq = mature.get(n)
            rows.append({"miRNAname": n, "found_in_mature_fa": bool(seq), "seq_len": len(seq) if seq else 0})
            if not seq:
                continue
            w.write(f">{n}\n")
            for i in range(0, len(seq), 80):
                w.write(seq[i : i + 80] + "\n")
    return pd.DataFrame(rows)


_RNAHYBRID_COMPACT = re.compile(
    r"^(?P<target>[^:]+):(?P<pos>\d+):(?P<mirna>[^:]+):(?P<seedlen>\d+):(?P<mfe>-?\d+(?:\.\d+)?):"
)


def _parse_rnahybrid_raw(raw_txt: Path) -> pd.DataFrame:
    """
    Parser for RNAhybrid ``-c`` compact output lines (colon-delimited), plus passthrough for diagnostics.
    Always keeps `raw_line` for auditability.
    """
    rows = []
    with raw_txt.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            m = _RNAHYBRID_COMPACT.match(s)
            if m:
                rows.append(
                    {
                        "target": m.group("target"),
                        "target_hit_start": int(m.group("pos")),
                        "miRNA": m.group("mirna"),
                        "seed_len": int(m.group("seedlen")),
                        "mfe": float(m.group("mfe")),
                        "raw_line": s,
                    }
                )
            else:
                rows.append(
                    {
                        "target": None,
                        "target_hit_start": None,
                        "miRNA": None,
                        "seed_len": None,
                        "mfe": None,
                        "raw_line": s,
                    }
                )
    return pd.DataFrame(rows)


def run_rnahybrid_predictions(
    *,
    mirna_fasta: Path,
    genome_fasta: Path | None = None,
    gencode_parquet: Path | None = None,
    n_extra_close_lncrnas: int = 20,
    paths: RnaHybridPaths = RnaHybridPaths(),
    rnahybrid_extra_args: Optional[list[str]] = None,
    encori_mirna_parquet: Path | None = None,
    mirna_subset_from_encori: bool = False,
    max_target_bp: int | None = 8000,
) -> None:
    """
    Build spliced exon FASTA for selected lncRNAs and run RNAhybrid miRNA target predictions.

    Requirements:
    - `bedtools` in PATH
    - `RNAhybrid` in PATH
    - A genome FASTA (default: `PATHS.sv_reference_fasta`)
    - A miRNA FASTA of **mature sequences** (e.g. `>hsa-miR-21-5p`)
    """
    genome_fasta = genome_fasta or Path(PATHS.sv_reference_fasta)
    gencode_parquet = gencode_parquet or Path(PATHS.gencode_gtf_full_pq)

    selected_lncrnas = build_encori_lncrna_target_list(
        panel_lncrnas=TIER1_LNCRNA_GENES,
        lncrna_proximity_pairs_csv=Path(PATHS.working_dir) / "lncRNA_matching" / "genes_lncRNAs_1000000bp_distances.csv",
        primary_genes=PRIMARY_GENES,
        n_extra_close_lncrnas=int(n_extra_close_lncrnas),
    )

    paths.out_dir.mkdir(parents=True, exist_ok=True)
    _write_lncrna_exon_bed(
        lncrnas=selected_lncrnas,
        gencode_parquet=gencode_parquet,
        out_bed=paths.exon_bed,
        genome_fasta=genome_fasta,
    )
    _bedtools_getfasta(genome_fa=genome_fasta, bed=paths.exon_bed, out_fa=paths.exons_fa)
    _splice_exon_fasta(paths.exons_fa, paths.spliced_fa)

    q_fa = mirna_fasta
    if mirna_subset_from_encori:
        encori_mirna_parquet = encori_mirna_parquet or (Path(PATHS.working_dir) / "lncRNA_interactions" / "encori_mirna_targets.parquet")
        diag = write_mirna_subset_fasta_from_encori(
            encori_mirna_parquet=encori_mirna_parquet,
            mature_mirna_fasta=mirna_fasta,
            out_fa=paths.mirna_subset_fa,
        )
        diag.to_csv(paths.out_dir / "mirna_subset_from_encori.csv", index=False)
        q_fa = paths.mirna_subset_fa

    target_fa = paths.spliced_fa
    trunc_note = None
    if max_target_bp is not None:
        _truncate_fasta_sequences_in_place(
            in_fa=paths.spliced_fa,
            out_fa=paths.spliced_fa_trunc,
            max_bp=int(max_target_bp),
            from_end=True,
        )
        target_fa = paths.spliced_fa_trunc
        trunc_note = f"3prime_last_{int(max_target_bp)}bp"

    max_query_len = _max_fasta_seq_len(q_fa)
    max_target_len = _max_fasta_seq_len(target_fa)
    if max_query_len <= 0 or max_target_len <= 0:
        raise RuntimeError("RNAhybrid input FASTA appears empty after building targets/queries.")

    cmd = ["RNAhybrid", "-q", str(q_fa), "-t", str(target_fa)]
    # Sensible defaults: one best hit per target, compact output, human 3'UTR energy model, energy cutoff.
    # Users can override by passing explicit `rnahybrid_extra_args`.
    if rnahybrid_extra_args is None:
        cmd.extend(
            [
                "-b",
                "1",
                "-c",
                "-e",
                "-20",
                "-s",
                "3utr_human",
                "-n",
                str(max(32, max_query_len)),
                "-m",
                str(max(2000, max_target_len)),
            ]
        )
    else:
        cmd.extend(list(rnahybrid_extra_args))
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    paths.raw_txt.write_text(p.stdout, encoding="utf-8")

    hits = _parse_rnahybrid_raw(paths.raw_txt)
    hits["stderr"] = p.stderr
    hits["exit_code"] = int(p.returncode)
    hits["cmd"] = " ".join(cmd)
    hits["max_target_bp_truncation"] = trunc_note
    hits.to_parquet(paths.hits_parquet, index=False)

