"""⭐ FIND EVERY CODE REFERENCE TO A CARD COLUMN THAT NO LONGER EXISTS — the prune-ripple check.

Third of the standing checks, after `nan_bool_audit` (booleans that cannot say "unknown") and
`ratio_blowup_audit` (ratios with an ungated denominator). This one exists because **three separate silent
failures came from MY OWN prunes in one session** (MH-266):

  · `seedless_chimeric_candidates()` read `if "echim_any" in s.columns: s = s[...]`. The prune did not
    raise — the **filter simply stopped applying**, and the function began returning un-filtered rows.
  · `gene_landscape/axes_build.py` built the `chim_frac` axis from `echim_any` and listed `n_regulators`
    and `n_dense_included` in its `CARD` axis set. All guarded ⇒ **the axes silently stopped being built**,
    which in an FDR scan reads as *tested and null*.
  · `adm_expressed` was pruned in a normaliser, printed `✅ DROPPED`, and was still delivered — twice.

⛔ **THE MECHANISM, and it is the whole reason a grep is not enough.** Defensive access is idiomatic and
usually right — `if c in d.columns`, `.get(c)`, `[c for c in COLS if c in d.columns]`. But it converts a
**removed column** from a loud `KeyError` into a **silently skipped branch**. The safer the code looks, the
quieter the failure. ⇒ a prune is not done when the tests pass; it is done when nothing still names it.

**HOW IT KNOWS WHAT WAS REMOVED — and why the obvious two ways both fail.**
  ⛔ *Prefix heuristic* (first attempt): flag any literal carrying a live block prefix that is not a current
     column. **5,393 hits, essentially all false** — `is_` and `w_` are live card prefixes, so every
     internal `is_directed` / `w_arch` matched. Same failure as the first versions of `nan_bool_audit` and
     `ratio_blowup_audit`: a blunt rule buries the real hits, and a report nobody trusts is worse than none.
  ⛔ *Git history of the registry*: `output/` is **gitignored**, so the registry has **0 commits** — there
     is no history to diff against.
  ✅ **A TRACKED SNAPSHOT.** `docs/derived/card_columns.txt` holds one line per delivered column and IS in
     git. The audit diffs it against the live registry: names in the snapshot but absent now are **exactly
     the retired set**, with no heuristic and no hand-maintained list. Pruning a column then shows up as a
     line removed in a reviewable diff, and `--update` re-baselines once the ripple is cleared.

    .venv/bin/python3 -m mirna_hallmark.analyses.ops.column_ref_audit
    .venv/bin/python3 -m mirna_hallmark.analyses.ops.column_ref_audit --all   # include adjudicated
"""
from __future__ import annotations

import ast
import pathlib
import re
import sys

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
REGISTRY = ROOT / "output" / "learned" / "card_registry.tsv"

#: columns deliberately removed or renamed — the reference is EXPECTED and has been dealt with.
#: Each entry records where it is still legitimately named, so the row can be re-checked, not just muted.
RETIRED = {
    "echim_any":            "pruned unit 21; refs fixed to echim_n_sources.notna() (MH-266)",
    "adm_expressed":        "pruned unit 25 -> adm_expr_tcga; producer + post-annotation both fixed",
    "n_regulators":         "card copy renamed heur_n_regulators; the RETIRED heuristic lane still owns "
                            "the name in gene_corepression.tsv and its direct readers — expected",
    "gene_nreg":            "renamed heur_gene_nreg unit 24; rename machinery still names it",
    "n_discovered":         "renamed n_pip_disc_gt50; rename machinery still names it",
    "n_dense_included":     "pruned (bit-identical to n_fam); rename machinery still names it",
    "p_fam":                "pruned unit C; drop machinery still names it",
    "beta_frac_reliable":   "pruned unit 7; readouts still BUILDS it upstream of the card - expected",
    "net_pressure":         "pruned unit 7; readouts still builds it upstream - expected",
    "fame_led_n_pmid":      "pruned; arm_card computes then drops it - expected",
    "arb_n_genes":          "pruned (bit-identical to arb_n_edges)",
    "model_n_edges":        "pruned (bit-identical to model_n_genes)",
    "family_n_arms_card":   "pruned (bit-identical to famrole_n_members)",
    "fame_assay_perturbation_studies": "pruned - miRTarBase carries no perturbation assay",
    # ── MH-269 explicit renames. Each old name is still legitimately referenced in ONE of two ways:
    # a deliberate both-names transition shim, or a DIFFERENT frame that owns the name in its own scope.
    "family_size":       "renamed n_family_in_design (MH-269); only the rename table names it",
    "family_role":       "renamed arm_role_in_family; canonical_card accepts both through the transition",
    "family_dose_share": "renamed arm_share_of_family_dose; rename table + transition shim only",
    "arm_resolvable":    "renamed cell_arms_resolvable; card_context/arm_card/rung_parity/_alias accept both",
    "arm_dbeta":         "renamed cell_beta_spread (MH-279); rename table + both-names registrations only",
    "arm_sep_z":         "renamed cell_arm_sep_z (MH-279); rename table + both-names registrations only",
    "arm_pct_floor":     "renamed arm_pct_above_floor; card.py + _alias accept both",
    "field_retention":   "renamed field_excess_over_perm (it is r_own - r_perm, not a ratio)",
    "iso_total_rpm":     "renamed iso_rpm_summed (a SUM across samples)",
    # ⚠ `role` is NOT a stale card reference in most of its hits: `gene_roles.py` and every analysis reading
    # ITS frame own the name `role` in their own scope, correctly. Only the CARD copy became
    # `gene_cancer_role`. Renaming a source to fix a card is how a second collision starts.
    "role":              "CARD copy renamed gene_cancer_role (MH-269); gene_roles.py and its readers keep "
                         "`role` in their own frame — correct in that scope, not a stale card reference",
}
SKIP_DIRS = {"__pycache__", "archive", "docs"}
DEFENSIVE = re.compile(r"(\bin\s+\w+\.columns\b|\.get\(|\bif\b.*\bin\b|\bhasattr\()")


SNAPSHOT = ROOT / "docs" / "derived" / "card_columns.txt"


def _live_columns() -> set[str]:
    if not REGISTRY.exists():
        raise SystemExit(f"registry not found: {REGISTRY} — run `card_rungs --check` first")
    return set(pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("").column)


def _snapshot() -> set[str]:
    return set(SNAPSHOT.read_text().split()) if SNAPSHOT.exists() else set()


def _write_snapshot(cols: set[str]) -> None:
    SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
    SNAPSHOT.write_text("\n".join(sorted(cols)) + "\n")


def audit() -> list[dict]:
    live = _live_columns()
    retired = _snapshot() - live          # in the tracked snapshot, gone from the card = RETIRED
    retired |= set(RETIRED)               # plus the hand-adjudicated ones, so the first run is not empty
    out: list[dict] = []
    for path in sorted(ROOT.rglob("*.py")):
        if SKIP_DIRS & set(path.parts) or "card_glossary" in path.name:
            continue
        try:
            src = path.read_text(errors="ignore")
            tree = ast.parse(src)
        except (SyntaxError, OSError):
            continue
        lines = src.split("\n")
        for node in ast.walk(tree):
            if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
                continue
            v = node.value
            if v not in retired:          # ⭐ the ONLY rule: it was a column, and it is not one now
                continue
            ln = getattr(node, "lineno", 0)
            text = lines[ln - 1] if 0 < ln <= len(lines) else ""
            if text.lstrip().startswith("#"):
                continue
            out.append({"file": str(path.relative_to(ROOT.parent)), "line": ln, "name": v,
                        "guarded": bool(DEFENSIVE.search(text)),
                        "retired": RETIRED.get(v, ""), "src": text.strip()[:88]})
    return out


def selftest() -> int:
    """⭐ Simulate a prune and confirm the snapshot diff catches it. Run this before trusting a clean report.

    ⚠ Validating by un-adjudicating a name in `RETIRED` does NOT work and looks like a tool failure: those
    keys are the DETECTION set as well as the reason map, so removing one removes it from detection too.
    The honest test is to remove a column from the LIVE set while the snapshot still lists it — which is
    exactly what a prune does.
    """
    real = _live_columns
    probe = {"adm_has_site", "coupling_tum"}
    globals()["_live_columns"] = lambda: real() - probe
    try:
        hits = [h for h in audit() if h["name"] in probe and not h["retired"]]
    finally:
        globals()["_live_columns"] = real
    guarded = sum(h["guarded"] for h in hits)
    print(f"simulated prune of {sorted(probe)}:")
    print(f"  call sites that would break: {len(hits)}  ({guarded} behind a defensive guard)")
    clean = len([h for h in audit() if not h["retired"]])
    print(f"  open hits with the columns restored: {clean}")
    ok = len(hits) > 5 and guarded > 0 and clean == 0
    print("  ✅ VALIDATED — the snapshot diff catches a prune, guarded and direct, and is quiet otherwise."
          if ok else "  ⛔ SELFTEST FAILED — do not trust a clean report until this passes.")
    return 0 if ok else 1


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    if "--update" in sys.argv:
        cols = _live_columns()
        _write_snapshot(cols)
        print(f"✅ snapshot re-baselined: {len(cols)} columns -> {SNAPSHOT.relative_to(ROOT.parent)}")
        print("   commit it — the diff IS the record of what was pruned.")
        return 0
    hits = audit()
    show_all = "--all" in sys.argv
    adjudicated = [h for h in hits if h["retired"]]
    open_ = [h for h in hits if not h["retired"]]
    print(f"live column-name references not present on any card: {len(hits)}")
    print(f"  adjudicated (in RETIRED, with a reason): {len(adjudicated)}")
    print(f"  ⚠ OPEN — not on a card and not adjudicated: {len(open_)}   "
          f"({sum(h['guarded'] for h in open_)} behind a DEFENSIVE guard, i.e. failing silently)")
    rows = hits if show_all else open_
    if not rows:
        print("\n✅ nothing OPEN — every stale column reference has a recorded reason.")
        return 0
    print("\n  ⚠ these may be INTERNAL frame columns that never reach a card — read before acting:")
    for h in sorted(rows, key=lambda x: (not x["guarded"], x["file"], x["line"])):
        tag = "⛔ GUARDED" if h["guarded"] else "  direct  "
        print(f"  {tag} {h['file']}:{h['line']}  {h['name']}")
        print(f"             {h['src']}")
    print("\n  ⛔ GUARDED is the dangerous one: the reference cannot raise, so a removed column becomes a\n"
          "     skipped branch — a missing filter, or an axis that silently stops being built.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
