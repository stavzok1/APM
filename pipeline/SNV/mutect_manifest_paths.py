"""
Map GDC Mutect2 manifest ``File Name`` entries to on-disk VEP outputs.

Layout (typical for this repo):
  ``data/SNV/vcfs_snv`` — raw Mutect2 names as in the GDC manifest
  ``data/SNV/vep_vcfs`` — VEP + 1Mb-window subset, e.g.
      ``…raw_somatic_mutation.vcf`` → ``…raw_somatic_mutation.APM_1Mb.vep.vcf``
      (``…raw_somatic_mutation.APM_1Mb.vcf`` is the pre-VEP 1Mb subset and
       **does not carry VEP CSQ**; it must NOT be chosen for SNV processing.)
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional


def vep_apm1mb_filename_candidates(manifest_file_name: str) -> List[str]:
    """
    Possible on-disk VEP / APM-window filenames for one GDC Mutect2 manifest name.

    The resolver will try these in order and return the first hit. VEP-annotated
    names (``…APM_1Mb.vep.vcf[.gz]``) come first so downstream SNV processing
    actually sees ``CSQ`` and can populate gene_hits / regulatory_hits / motif_hits
    / has_* flags. Plain ``…APM_1Mb.vcf`` is the **pre-VEP** 1Mb subset and is
    intentionally LAST (only used when no VEP output exists on disk) — with a
    graceful-empty column fallback in the SNV loader.
    """
    fn = str(manifest_file_name).strip()
    if not fn or fn.lower() == "nan":
        return []
    out: List[str] = []
    if fn.endswith(".vcf.gz"):
        base = fn[: -len(".vcf.gz")]
        out.append(base + ".APM_1Mb.vep.vcf.gz")
        out.append(base + ".APM_1Mb.vep.vcf")
        out.append(base + ".APM_1Mb.vcf.gz")
        out.append(base + ".APM_1Mb.vcf")
    elif fn.endswith(".vcf"):
        base = fn[: -len(".vcf")]
        out.append(base + ".APM_1Mb.vep.vcf")
        out.append(base + ".APM_1Mb.vcf")
    out.append(fn)
    return out


def first_existing_child(rel_name: str, dirs: Iterable[Optional[Path]]) -> Optional[Path]:
    """Return ``dir / rel_name`` for the first directory where that file exists."""
    for d in dirs:
        if d is None:
            continue
        root = Path(d)
        if not root.is_dir():
            continue
        p = root / rel_name
        if p.is_file():
            return p
    return None


def resolve_vep_mutect_vcf(
    manifest_file_name: str,
    *,
    vep_dirs: Iterable[Optional[Path]],
) -> Optional[Path]:
    """
    Locate the VEP-processed Mutect2 VCF for a manifest ``File Name``.

    Tries each ``vep_apm1mb_filename_candidates`` under each directory in ``vep_dirs``.
    """
    for cand in vep_apm1mb_filename_candidates(manifest_file_name):
        hit = first_existing_child(cand, vep_dirs)
        if hit is not None:
            return hit
    return None
