#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import List

ROOT = Path(__file__).resolve().parents[2]


def _run(cmd: List[str]) -> None:
    print("\n$ " + " ".join(cmd))
    subprocess.check_call(cmd, cwd=str(ROOT))


def main() -> None:
    ap = argparse.ArgumentParser(description="Run the QC suite (build RPPA outputs + QC reports).")
    ap.add_argument("--scratch-json", type=str, required=True)
    ap.add_argument("--skip-rppa", action="store_true", help="Skip building cohort RPPA outputs.")
    ap.add_argument("--skip-cross-modal", action="store_true", help="Skip analysis/qc/run_cross_modal_qc.py.")
    ap.add_argument("--skip-full", action="store_true", help="Skip analysis/qc/run_full_qc.py.")
    ap.add_argument("--venv-python", type=str, default=".venv/bin/python3", help="Python executable to use.")
    args = ap.parse_args()

    py = args.venv_python
    scratch = args.scratch_json

    if not args.skip_rppa:
        _run([py, "analysis/qc/build_rppa_outputs_for_cohort.py", "--scratch-json", scratch])

    if not args.skip_cross_modal:
        _run([py, "analysis/qc/run_cross_modal_qc.py", "--scratch-json", scratch])

    if not args.skip_full:
        _run([py, "analysis/qc/run_full_qc.py", "--scratch-json", scratch])

    print("\n[OK] QC suite complete. See: analysis/qc/output/")


if __name__ == "__main__":
    main()

