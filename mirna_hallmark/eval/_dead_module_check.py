"""DEAD-MODULE GUARD — which modules in this subproject CANNOT run, and would present as "stale output"?

    .venv/bin/python3 -m mirna_hallmark.eval._dead_module_check            # sweep everything
    .venv/bin/python3 -m mirna_hallmark.eval._dead_module_check --probe mirna_hallmark.learned.card

WHY THIS EXISTS (2026-08-01). **Twice now, an artifact that looked merely STALE was produced by a module
that could not run at all**, and in both cases the module executed its FULL compute and died only at the
final write — so the failure was indistinguishable from neglect, and the lane sat "pending" for weeks:

  * **MH-38** — `eval/buffa_validation`: `GEO_DIR` was a `__file__`-relative hop, correct at the subproject
    top level, dead after the module moved into `eval/`. *"That is why this lane was never re-run: it
    could not run."*
  * **MH-168** — `learned/analyses/gene_atlas.py`: `parents[2]` was the repo root at `learned/gene_atlas.py`;
    commit `e5d5d84` relocated it one level deeper, so every output resolved to
    `mirna_hallmark/mirna_hallmark/output/learned/` — **dead for 15 days**, blocking the competence-map
    refresh that STATE_OF_PLAY listed as a top-4 next action.

⇒ The shared fingerprint is cheap to detect: **a module-level `Path` whose PARENT DIRECTORY does not
exist.** Nothing can be written there and nothing can be read from there. This guard finds those, plus
modules that cannot even be imported, plus the move-bug's literal signature (a duplicated path segment
such as `mirna_hallmark/mirna_hallmark`).

WHAT IT DOES NOT CLAIM. A clean sweep does **not** mean every module runs — paths built INSIDE functions,
missing input FILES, and runtime errors are out of scope. It rules out exactly the failure mode that has
actually bitten this subproject twice. Treat a hit as "verify this module", not as an automatic bug.

⚠ Each module is imported in its OWN subprocess: import-time side effects and crashes cannot contaminate
  the sweep or each other, and a hung import is reported rather than hanging the run.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUB = ROOT / "mirna_hallmark"
SKIP_DIRS = {"__pycache__", "archive", ".ipynb_checkpoints"}
TIMEOUT = 120


def _modules() -> list:
    out = []
    for p in sorted(SUB.rglob("*.py")):
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        rel = p.relative_to(ROOT).with_suffix("")
        parts = list(rel.parts)
        if parts[-1] == "__init__":
            parts = parts[:-1]
            if not parts:
                continue
        out.append(".".join(parts))
    return out


def probe(name: str) -> dict:
    """Import ONE module and report every module-level Path and its directory status."""
    import importlib

    res = {"module": name, "import": "ok", "paths": []}
    try:
        m = importlib.import_module(name)
    except BaseException as e:                       # noqa: BLE001 - any failure is a finding
        res["import"] = f"{type(e).__name__}: {e}"
        return res
    for attr in dir(m):
        if attr.startswith("__"):
            continue
        try:
            v = getattr(m, attr)
        except Exception:
            continue
        if not isinstance(v, Path):
            continue
        s = str(v)
        res["paths"].append({
            "attr": attr, "path": s,
            "exists": v.exists(),
            "parent_exists": v.parent.exists(),
            "dup_segment": "mirna_hallmark/mirna_hallmark" in s,
        })
    return res


def _run_one(name: str) -> dict:
    # ⚠ the target module name travels by ENV, not argv: several modules in this subproject parse
    # `sys.argv` AT IMPORT TIME, so a "--probe <mod>" flag leaked into them and they raised
    # (`KeyError: '--probe not in RNA matrix'` for the two `_niter_*` tests) — a harness artifact
    # masquerading as a dead module. An empty argv is the only neutral thing to hand them.
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               _DMC_PROBE=name)
    cmd = [sys.executable, "-m", "mirna_hallmark.eval._dead_module_check"]
    try:
        r = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {"module": name, "import": f"TIMEOUT >{TIMEOUT}s (import-time side effects?)", "paths": []}
    for line in reversed(r.stdout.splitlines()):
        if line.startswith("{"):
            try:
                return json.loads(line)
            except Exception:
                pass
    return {"module": name, "import": f"probe failed rc={r.returncode}: "
                                      f"{(r.stderr or r.stdout).strip().splitlines()[-1:] or ['?']}",
            "paths": []}


def main() -> None:
    mods = _modules()
    print(f"[dead_module_check] probing {len(mods)} modules (own subprocess each, {TIMEOUT}s cap)\n")
    with ThreadPoolExecutor(max_workers=8) as ex:
        results = list(ex.map(_run_one, mods))

    bad_import = [r for r in results if r["import"] != "ok"]
    dead_dir, dup = [], []
    for r in results:
        for p in r["paths"]:
            if p["dup_segment"]:
                dup.append((r["module"], p))
            elif not p["parent_exists"]:
                dead_dir.append((r["module"], p))

    print("=" * 100)
    print(f"⛔ DEAD OUTPUT/INPUT DIRECTORY — the MH-38 / MH-168 fingerprint  ({len(dead_dir)})")
    print("   a module-level Path whose PARENT does not exist: nothing can be written or read there")
    print("=" * 100)
    for mod, p in sorted(dead_dir):
        print(f"  {mod}\n      {p['attr']} = {p['path']}")
    if not dead_dir:
        print("  (none)")

    print(f"\n{'=' * 100}\n⛔ DUPLICATED PATH SEGMENT 'mirna_hallmark/mirna_hallmark'  ({len(dup)})")
    print("   the literal signature of a __file__ hop that outlived its module's location")
    print("=" * 100)
    for mod, p in sorted(dup):
        print(f"  {mod}\n      {p['attr']} = {p['path']}")
    if not dup:
        print("  (none)")

    print(f"\n{'=' * 100}\n⚠ CANNOT IMPORT  ({len(bad_import)})\n{'=' * 100}")
    for r in sorted(bad_import, key=lambda r: r["module"]):
        print(f"  {r['module']:62s} {r['import'][:110]}")
    if not bad_import:
        print("  (none)")

    n_paths = sum(len(r["paths"]) for r in results)
    print(f"\n{'=' * 100}")
    print(f"  {len(results)} modules · {n_paths} module-level Paths inspected · "
          f"{len(dead_dir)} dead dirs · {len(dup)} dup-segment · {len(bad_import)} import failures")
    print("  ⚠ A clean sweep does NOT mean every module runs — paths built inside functions, missing")
    print("    input FILES, and runtime errors are out of scope. This rules out ONE failure mode.")
    print("=" * 100)

    out = ROOT / "mirna_hallmark/output/learned/dead_module_check.tsv"
    import pandas as pd
    rows = [{"module": r["module"], "import": r["import"], "attr": p["attr"], "path": p["path"],
             "exists": p["exists"], "parent_exists": p["parent_exists"], "dup_segment": p["dup_segment"]}
            for r in results for p in (r["paths"] or [{"attr": "", "path": "", "exists": None,
                                                       "parent_exists": None, "dup_segment": None}])]
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"-> {out}")


if __name__ == "__main__":
    _p = os.environ.get("_DMC_PROBE")
    if _p:
        sys.argv = [sys.argv[0]]        # neutral argv for modules that parse it at import time
        print(json.dumps(probe(_p)))
    else:
        argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args()
        main()
