"""
Cohort-level covariates builder.

This package is intentionally **not** invoked by ``pipeline.main``.
It is designed to run after cohort processing to produce a unified, sample-keyed
covariates table for downstream modeling.
"""

from .build_covariates import build_cohort_covariates

__all__ = ["build_cohort_covariates"]

