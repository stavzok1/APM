"""miRNA x Hallmark regulation subproject.

A self-contained analysis subproject (separate from ``analysis/``) studying how
miRNAs -- weighted by miRTarBase experimental evidence and gated by AGO/RISC
availability -- regulate the 50 MSigDB Hallmark gene sets, with CNV + RNA
expression context and stratum characterization (PAM50 / TNBC / immune / stage).

See ``mirna_hallmark/AGENTS.md`` and ``mirna_hallmark/docs/`` for onboarding and
its own catalog / discovery-registry / run-ledger.
"""

from __future__ import annotations

__all__ = ["config"]
