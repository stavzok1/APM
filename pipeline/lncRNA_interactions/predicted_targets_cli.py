from __future__ import annotations

import argparse
from pathlib import Path

from pipeline.lncRNA_interactions.predicted_targets_rnahybrid import RnaHybridPaths, run_rnahybrid_predictions


def main() -> None:
    ap = argparse.ArgumentParser(description="Predict miRNA→lncRNA binding using RNAhybrid (spliced exons).")
    ap.add_argument("--mirna-fasta", type=Path, required=True, help="FASTA of mature miRNA sequences (arm-level names).")
    ap.add_argument(
        "--mirna-subset-from-encori",
        action="store_true",
        help="Restrict miRNAs to those observed in data/lncRNA_interactions/encori_mirna_targets.parquet (recommended).",
    )
    ap.add_argument("--encori-mirna-parquet", type=Path, default=None, help="Override ENCORI miRNA parquet path.")
    ap.add_argument("--genome-fasta", type=Path, default=None, help="Genome FASTA (default: PATHS.sv_reference_fasta).")
    ap.add_argument(
        "--max-target-bp",
        type=int,
        default=8000,
        help="Truncate each spliced lncRNA target to the last N bp (3' end) before RNAhybrid. Use 0 to disable.",
    )
    ap.add_argument("--n-extra-close-lncrnas", type=int, default=20)
    ap.add_argument("--out-dir", type=Path, default=None, help="Override output dir (default under data/lncRNA_interactions/).")
    ap.add_argument(
        "--rnahybrid-arg",
        action="append",
        default=[],
        help="Extra RNAhybrid args (repeatable), e.g. --rnahybrid-arg=-b --rnahybrid-arg=1 --rnahybrid-arg=-c",
    )
    args = ap.parse_args()

    paths = RnaHybridPaths()
    if args.out_dir is not None:
        paths = RnaHybridPaths(out_dir=args.out_dir)

    run_rnahybrid_predictions(
        mirna_fasta=args.mirna_fasta,
        genome_fasta=args.genome_fasta,
        n_extra_close_lncrnas=int(args.n_extra_close_lncrnas),
        paths=paths,
        rnahybrid_extra_args=list(args.rnahybrid_arg) if args.rnahybrid_arg else None,
        encori_mirna_parquet=args.encori_mirna_parquet,
        mirna_subset_from_encori=bool(args.mirna_subset_from_encori),
        max_target_bp=None if int(args.max_target_bp) == 0 else int(args.max_target_bp),
    )
    print(f"Wrote: {paths.hits_parquet}")


if __name__ == "__main__":
    main()

