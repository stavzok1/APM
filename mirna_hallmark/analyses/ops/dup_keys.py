"""⛔ REPO-WIDE DUPLICATE-KEY SCAN — is the silent-overwrite hazard structural, or was it my lapse?

    PYTHONPATH=. .venv/bin/python3 -m mirna_hallmark.analyses.ops.dup_keys [--all]

WHY. A Python dict literal accepts a repeated key **silently** — resolved at construction, LAST one wins,
no error and no warning. In `card_glossary` that killed four freshly-written entries AND, once a guard
existed, immediately exposed an older live instance (a circular stub had been overwriting the real
`identity_` description for every `identity_*` column). Two instances in one file is a pattern, not a
lapse — but "two in one file" is not evidence about the REPO. This scans every module so the claim can be
made or withdrawn on measurement instead of on impression.

⚠ SET AND DICT-COMPREHENSION keys are NOT checked — only literal `{k: v}` dicts and `set` literals, where
the duplicate is unambiguous and always a mistake. A key built at runtime can legitimately repeat.

⛔⛔ **BLIND SPOT CLOSED 2026-08-19 — the within-literal check was NOT ENOUGH, and it was missing live
shadowing at the time it was written.** A key can also be duplicated ACROSS two literals that are later
merged into one mapping (`D.update({...})`, `{**A, **B}`), and the merge resolves it exactly as silently as
a literal does. `card_glossary.COLUMNS` is built from a literal at line 531 plus a `COLUMNS.update({...})`
at line 878 — **three keys were shadowed across that boundary while this module reported the file CLEAN**
(`n_arms`, `n_cell_intrinsic`, `n_composition_explained`, two of which disagreed about whether the unit was
edges or families). ⇒ `scan_merges()` now tracks dict literals by the NAME they are merged into and reports
keys appearing in more than one. **The lesson is the tool's own: a guard that checks one shape reports
"clean" for every other shape, and that reads as coverage.**

⚠ A duplicate is not automatically a BUG — `{"a": 1, "a": 1}` is harmless. What makes it dangerous is a
duplicate with DIFFERENT values, which is silently resolved in favour of the later one. Both are reported,
the differing ones first, because only those can change behaviour.
"""
from __future__ import annotations

import ast
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]


def scan_merges(path: pathlib.Path) -> list[dict]:
    """Keys duplicated ACROSS literals that are merged into the same name — the blind spot above.

    Tracks two shapes: `NAME = {...}` / `NAME: T = {...}` and `NAME.update({...})`. Anything merged at
    runtime from a non-literal is out of scope, deliberately — a key built at runtime may legitimately
    repeat, and guessing would put false positives in front of a reader who then stops trusting the report.
    """
    try:
        tree = ast.parse(path.read_text(errors="ignore"))
    except SyntaxError:
        return []
    byname: dict[str, list[tuple[ast.Dict, str]]] = {}
    # ⛔ MODULE LEVEL ONLY. The first version walked the whole tree and scored **85 hits, almost all false**:
    # a local `row = {...}` in one function and another `row = {...}` in the next are unrelated mappings that
    # merely share a variable name. Only a module-level name can actually be merged across literals, which is
    # the hazard this looks for. (Same lesson as `nan_bool_audit`'s first version — a blunt first heuristic
    # buries the real hits under noise, and a report nobody trusts is worse than no report.)
    for node in ast.iter_child_nodes(tree):
        if isinstance(node, ast.Expr):
            node = node.value
        name = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) \
                and isinstance(node.value, ast.Dict):
            name, dnode = node.targets[0].id, node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) \
                and isinstance(node.value, ast.Dict):
            name, dnode = node.target.id, node.value
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == "update" and isinstance(node.func.value, ast.Name) \
                and len(node.args) == 1 and isinstance(node.args[0], ast.Dict):
            name, dnode = node.func.value.id, node.args[0]
        if name is None:
            continue
        byname.setdefault(name, []).append((dnode, name))
    out: list[dict] = []
    for name, lits in byname.items():
        if len(lits) < 2:
            continue
        seen: dict[str, tuple[int, str]] = {}
        for dnode, _ in lits:
            for k, v in zip(dnode.keys, dnode.values):
                if k is None:
                    continue
                try:
                    tag = repr(ast.literal_eval(k))
                except Exception:
                    continue
                val = ast.dump(v)
                if tag in seen and seen[tag][0] != getattr(k, "lineno", -1):
                    out.append({"file": path, "key": tag, "name": name,
                                "first": seen[tag][0], "second": getattr(k, "lineno", -1),
                                "differs": seen[tag][1] != val, "cross": True})
                else:
                    seen[tag] = (getattr(k, "lineno", -1), val)
    return out


def scan_file(path: pathlib.Path) -> list[dict]:
    try:
        tree = ast.parse(path.read_text(errors="ignore"))
    except SyntaxError:
        return []
    out: list[dict] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Dict):
            continue
        seen: dict[str, tuple[int, str]] = {}
        for k, v in zip(node.keys, node.values):
            if k is None:                      # {**spread}
                continue
            try:
                tag = repr(ast.literal_eval(k))
            except Exception:
                continue                       # runtime-built key: may legitimately repeat
            try:
                val = ast.unparse(v)[:80]
            except Exception:
                val = "?"
            if tag in seen:
                first_line, first_val = seen[tag]
                out.append({"file": str(path.relative_to(ROOT.parent)), "key": tag,
                            "first": first_line, "second": k.lineno,
                            "same_value": first_val == val})
            else:
                seen[tag] = (k.lineno, val)
    return out


def main() -> None:
    show_all = "--all" in sys.argv[1:]
    files = [p for p in ROOT.rglob("*.py") if "__pycache__" not in p.parts]
    hits: list[dict] = []
    merges: list[dict] = []
    for p in files:
        hits.extend(scan_file(p))
        merges.extend(scan_merges(p))
    differing = [h for h in hits if not h["same_value"]]
    identical = [h for h in hits if h["same_value"]]

    print(f"scanned {len(files)} modules under {ROOT.name}/")
    print(f"  duplicate dict keys: {len(hits)}  "
          f"({len(differing)} with DIFFERENT values — these silently change behaviour; "
          f"{len(identical)} harmless repeats)")
    m_diff = [h for h in merges if h["differs"]]
    print(f"  keys shadowed ACROSS merged literals: {len(merges)}  ({len(m_diff)} with DIFFERENT values)")
    if differing:
        print("\n⛔ DIFFERING VALUES — the later one wins and the earlier is DEAD CODE:")
        for h in sorted(differing, key=lambda x: (x["file"], x["first"])):
            print(f"   {h['file']}:{h['first']} -> overwritten at :{h['second']}   key={h['key']}")
    if m_diff:
        print("\n⛔ SHADOWED ACROSS A MERGE — same mapping, two literals, the later wins silently:")
        for h in sorted(m_diff, key=lambda x: (x["file"], x["first"])):
            print(f"   {h['file']}:{h['first']} -> shadowed at :{h['second']}   "
                  f"{h['name']}[{h['key']}]")
    else:
        print("\n✅ no duplicate key with differing values anywhere in the subproject.")
    if identical and show_all:
        print(f"\n(harmless exact repeats, {len(identical)}:)")
        for h in identical[:40]:
            print(f"   {h['file']}:{h['first']} == :{h['second']}   key={h['key']}")
    if not show_all and identical:
        print(f"\n(+{len(identical)} harmless exact repeats — pass --all to list)")
    raise SystemExit(1 if differing else 0)


if __name__ == "__main__":
    main()
