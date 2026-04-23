from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

from pipeline.config import PATHS
from pipeline.sample_ids import (
    add_tcga_id_columns_from_index_inplace,
    tcga_best_join_key,
    tcga_sample_type_two_digit,
)

from .io import CovariatesRunPaths, default_run_id, resolve_output_paths, write_covariates_outputs
from .providers import (
    ClinicalImmuneProvider,
    CoverageProvider,
    ExternalCovariateSpec,
    ExternalTableProvider,
    ProviderResult,
    RnaSignatureProvider,
    RppaProvider,
    MutationMafProvider,
    DdrScoresProvider,
    HiChipParticipantProvider,
    SnvNativeProvider,
    SvDisruptionProvider,
)


def _index_union(dfs: Iterable[pd.DataFrame], *, name: str) -> pd.Index:
    idx: List[str] = []
    for df in dfs:
        if df is None or df.empty:
            continue
        idx.extend([str(x) for x in df.index.astype(str).tolist()])
    out = pd.Index(sorted(set(idx)), name=name, dtype="string")
    return out


def _lift_sample_to_sample_vial(
    df_sample: pd.DataFrame,
    *,
    target_vials: pd.Index,
) -> pd.DataFrame:
    """
    Lift a sample-keyed table (…-01) onto a sample_vial index (…-01A)
    by joining on the derived 'sample' column from normalized IDs.
    """
    if df_sample.empty:
        return pd.DataFrame(index=target_vials)

    tmp = df_sample.copy()
    if "sample" not in tmp.columns:
        add_tcga_id_columns_from_index_inplace(tmp, overwrite=False)

    if "sample" not in tmp.columns:
        # Cannot lift; return empty aligned frame.
        return pd.DataFrame(index=target_vials)

    # Build a mapping sample_vial -> sample (string)
    vials = pd.DataFrame(index=target_vials)
    add_tcga_id_columns_from_index_inplace(vials, overwrite=True)
    sample_for_vial = vials["sample"].astype(str)

    # tmp is keyed by sample; make sure that's the index.
    tmp["sample"] = tmp["sample"].astype("string")
    tmp.index = tmp["sample"].where(tmp["sample"].notna(), tmp.index.astype(str)).astype(str)
    tmp.index.name = "sample"
    if tmp.index.has_duplicates:
        # Sample-level sources can have multiple rows per sample; keep the first to allow lifting.
        tmp = tmp[~tmp.index.duplicated(keep="first")].copy()

    lifted = tmp.reindex(sample_for_vial.values)
    lifted.index = target_vials
    lifted.index.name = "sample_vial"
    return lifted


def _lift_participant_to_sample_vial(
    df_participant: pd.DataFrame,
    *,
    target_vials: pd.Index,
) -> pd.DataFrame:
    if df_participant.empty:
        return pd.DataFrame(index=target_vials)

    tmp = df_participant.copy()
    if "participant" not in tmp.columns:
        add_tcga_id_columns_from_index_inplace(tmp, overwrite=False)
    if "participant" not in tmp.columns:
        return pd.DataFrame(index=target_vials)

    tmp["participant"] = tmp["participant"].astype("string")
    tmp.index = tmp["participant"].where(tmp["participant"].notna(), tmp.index.astype(str)).astype(str)
    tmp.index.name = "participant"
    if tmp.index.has_duplicates:
        tmp = tmp[~tmp.index.duplicated(keep="first")].copy()

    vials = pd.DataFrame(index=target_vials)
    add_tcga_id_columns_from_index_inplace(vials, overwrite=True)
    part_for_vial = vials["participant"].astype(str)

    lifted = tmp.reindex(part_for_vial.values)
    lifted.index = target_vials
    lifted.index.name = "sample_vial"
    return lifted


def build_cohort_covariates(
    *,
    run_id: Optional[str] = None,
    out_dir: Optional[Path] = None,
    coverage_run_id: Optional[str] = None,
    rppa_combined_path: Optional[Path] = None,
    rna_tpm_path: Optional[Path] = None,
    clinical_unified_path: Optional[Path] = None,
    external_tables: Optional[List[ExternalCovariateSpec]] = None,
    mutation_maf: Optional[Path] = None,
    ddr_scores_path: Optional[Path] = None,
    include_hichip_participant: bool = True,
    write_csv: bool = True,
) -> Tuple[pd.DataFrame, Dict[str, Any], CovariatesRunPaths]:
    """
    Build a unified covariates table keyed by TCGA ``sample_vial``.

    This function is intentionally separate from ``pipeline.main``.
    """
    rid = run_id or default_run_id("run")

    # Output under the pipeline working dir (data/), unless explicit out_dir given.
    repo_data_dir = Path(PATHS.working_dir)
    paths = resolve_output_paths(repo_data_dir=repo_data_dir, run_id=rid)
    if out_dir is not None:
        out_dir = Path(out_dir)
        paths = CovariatesRunPaths(
            out_dir=out_dir,
            table_parquet=out_dir / paths.table_parquet.name,
            table_csv=out_dir / paths.table_csv.name,
            metadata_json=out_dir / paths.metadata_json.name,
        )

    providers = []
    providers.append(CoverageProvider(analysis_root=Path("."), run_id=coverage_run_id))
    providers.append(
        ClinicalImmuneProvider(
            unified_tsv=Path(clinical_unified_path or PATHS.brca_clinical_immune_unified)
        )
    )
    providers.append(
        RppaProvider(
            combined_path=Path(rppa_combined_path or (PATHS.rppa_processed_dir / "rppa_analysis_combined.parquet"))
        )
    )
    providers.append(
        RnaSignatureProvider(
            tpm_path=Path(rna_tpm_path or PATHS.rna_expression),
            gene_col=PATHS.rna_gene_col,
        )
    )

    # Pipeline-native SNV per-sample summaries (best-effort).
    providers.append(SnvNativeProvider(snv_output_dir=PATHS.snv_output_dir))

    # Canonical external covariates: DDR/HRD/TP53 score table (BRCA-only filtered).
    ddr_path = Path(ddr_scores_path) if ddr_scores_path else (Path(PATHS.working_dir) / "covariates" / "DDRscores.txt")
    if ddr_path.exists():
        providers.append(DdrScoresProvider(path=ddr_path, cancer_acronym="BRCA"))

    # SV → gene mirroring (gene-centric disruption summaries for stratifiers + BRCA1/2).
    providers.append(
        SvDisruptionProvider(
            sv_output_root=PATHS.sv_output_root,
            stage_dir="07_final_sv_with_fimo",
            genes_of_interest=["TP53", "PTEN", "PIK3CA", "BRCA1", "BRCA2"],
        )
    )

    if include_hichip_participant:
        providers.append(HiChipParticipantProvider(path=PATHS.hichip_tcga_processed_csv))

    if mutation_maf is not None:
        providers.append(MutationMafProvider(maf_path=Path(mutation_maf)))

    for spec in list(external_tables or []):
        providers.append(ExternalTableProvider(spec))

    results: List[ProviderResult] = [p.build() for p in providers]
    diag: Dict[str, Any] = {"providers": {r.name: r.diagnostics for r in results}}

    # Determine canonical index from strongest sample_vial-keyed sources.
    sample_vial_frames = [r.df for r in results if r.index_key == "sample_vial"]
    base_index = _index_union(sample_vial_frames, name="sample_vial")
    if len(base_index) == 0:
        # fall back: union over any provider index (raw), then normalize to best_join_key
        raw = _index_union([r.df for r in results], name="raw_id")
        base_index = pd.Index([x for x in raw if x], name="sample_vial", dtype="string")

    base = pd.DataFrame(index=base_index)
    add_tcga_id_columns_from_index_inplace(base, overwrite=True)
    base["tcga_sample_type"] = base.index.map(lambda x: tcga_sample_type_two_digit(str(x)))
    base["tcga_best_join_key"] = base.index.map(lambda x: tcga_best_join_key(str(x)) or "")

    # Vial plurality / replication context.
    if "sample" in base.columns:
        sample_counts = base["sample"].astype(str).value_counts(dropna=False)
        base["n_vials_per_sample"] = base["sample"].astype(str).map(sample_counts).astype("Int64")
    if "participant" in base.columns:
        part_counts = base["participant"].astype(str).value_counts(dropna=False)
        base["n_vials_per_participant"] = base["participant"].astype(str).map(part_counts).astype("Int64")

    # Merge in providers.
    out = base
    for r in results:
        if r.df.empty:
            continue
        df_r = r.df
        if df_r.index.has_duplicates:
            # Avoid hard failure on reindex; keep first occurrence and record in diagnostics.
            n_dupe = int(df_r.index.duplicated(keep="first").sum())
            diag.setdefault("provider_deduped_index", {})[r.name] = {"n_duplicates_dropped": n_dupe}
            df_r = df_r[~df_r.index.duplicated(keep="first")].copy()
        if r.index_key == "sample_vial":
            aligned = df_r.reindex(out.index)
        elif r.index_key == "sample":
            aligned = _lift_sample_to_sample_vial(df_r, target_vials=out.index)
        elif r.index_key == "participant":
            aligned = _lift_participant_to_sample_vial(df_r, target_vials=out.index)
        else:
            # For now: skip unknown granularity (participant/aliquot/raw) unless user provides
            # an explicit mapping strategy. Keep the slot to avoid accidental mis-joins.
            aligned = pd.DataFrame(index=out.index)
        # Prefix columns to prevent collisions.
        prefix = r.name.replace("external:", "external__").replace("/", "_")
        aligned = aligned.copy()
        aligned.columns = [f"{prefix}__{c}" for c in aligned.columns]
        out = out.join(aligned, how="left")

    # Write.
    write_covariates_outputs(
        out,
        paths=paths,
        metadata=diag,
        write_csv=write_csv,
    )
    return out, diag, paths


def main() -> None:
    ap = argparse.ArgumentParser(description="Build cohort-level covariates table (post-processing).")
    ap.add_argument("--run-id", type=str, default="", help="Run id (default: timestamped).")
    ap.add_argument("--out", type=str, default="", help="Output directory (default: data/covariates/<run_id>).")
    ap.add_argument("--coverage-run-id", type=str, default="", help="analysis/sample_coverage run id to use.")
    ap.add_argument("--rppa-combined", type=str, default="", help="Path to rppa_analysis_combined.parquet/csv.")
    ap.add_argument("--rna-tpm", type=str, default="", help="Path to processed wide TPM table.")
    ap.add_argument("--clinical-unified", type=str, default="", help="Path to BRCA clinical+immune unified TSV.")
    ap.add_argument("--external", action="append", default=[], help="External covariate spec: name=path[:sep][:id_col]")
    ap.add_argument("--maf", type=str, default="", help="MAF-like somatic mutation table for TP53/PIK3CA covariates.")
    ap.add_argument("--ddr-scores", type=str, default="", help="Path to DDRscores.txt (defaults to data/covariates/DDRscores.txt).")
    ap.add_argument("--no-hichip", action="store_true", help="Skip HiChIP participant covariates.")
    ap.add_argument("--no-csv", action="store_true", help="Skip writing CSV (parquet only).")
    args = ap.parse_args()

    ext_specs: List[ExternalCovariateSpec] = []
    for raw in args.external:
        # Format: name=path[:sep][:id_col]
        name, rest = raw.split("=", 1) if "=" in raw else (f"ext{len(ext_specs)+1}", raw)
        parts = rest.split(":")
        path = Path(parts[0])
        sep = parts[1] if len(parts) >= 2 and parts[1] else "\t"
        id_col = parts[2] if len(parts) >= 3 and parts[2] else "sample_id"
        ext_specs.append(ExternalCovariateSpec(name=name, path=path, sep=sep, id_col=id_col))

    df, diag, paths = build_cohort_covariates(
        run_id=args.run_id or None,
        out_dir=Path(args.out) if args.out else None,
        coverage_run_id=args.coverage_run_id or None,
        rppa_combined_path=Path(args.rppa_combined) if args.rppa_combined else None,
        rna_tpm_path=Path(args.rna_tpm) if args.rna_tpm else None,
        clinical_unified_path=Path(args.clinical_unified) if args.clinical_unified else None,
        external_tables=ext_specs,
        mutation_maf=Path(args.maf) if args.maf else None,
        ddr_scores_path=Path(args.ddr_scores) if args.ddr_scores else None,
        include_hichip_participant=not args.no_hichip,
        write_csv=not args.no_csv,
    )
    print(f"Wrote: {paths.table_parquet}")
    if not args.no_csv:
        print(f"Wrote: {paths.table_csv}")
    print(f"Metadata: {paths.metadata_json}")
    print(f"Rows: {len(df)}  Cols: {len(df.columns)}")


if __name__ == "__main__":
    main()

