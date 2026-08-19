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

⚠ A duplicate is not automatically a BUG — `{"a": 1, "a": 1}` is harmless. What makes it dangerous is a
duplicate with DIFFERENT values, which is silently resolved in favour of the later one. Both are reported,
the differing ones first, because only those can change behaviour.
"""
from __future__ import annotations

import ast
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]


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
    for p in files:
        hits.extend(scan_file(p))
    differing = [h for h in hits if not h["same_value"]]
    identical = [h for h in hits if h["same_value"]]

    print(f"scanned {len(files)} modules under {ROOT.name}/")
    print(f"  duplicate dict keys: {len(hits)}  "
          f"({len(differing)} with DIFFERENT values — these silently change behaviour; "
          f"{len(identical)} harmless repeats)")
    if differing:
        print("\n⛔ DIFFERING VALUES — the later one wins and the earlier is DEAD CODE:")
        for h in sorted(differing, key=lambda x: (x["file"], x["first"])):
            print(f"   {h['file']}:{h['first']} -> overwritten at :{h['second']}   key={h['key']}")
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
