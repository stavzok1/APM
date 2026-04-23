#!/usr/bin/env python3
"""Compare HOCOMOCO H14CORE MEME vs JASPAR-format bundle against SV_TARGET_TF_SYMBOLS."""
from __future__ import annotations

import re
import sys
from pathlib import Path

MEME_DEFAULT = Path("/home/stavz/HOCOMOCO/H14CORE_meme_format.meme")
JASPAR_DEFAULT = Path("/home/stavz/HOCOMOCO/H14CORE_jaspar_format.txt")

TARGETS = [
    "STAT1",
    "STAT2",
    "STAT3",
    "IRF1",
    "IRF2",
    "IRF3",
    "NFKB1",
    "NFKB2",
    "RELA",
    "p65",
    "RELB",
    "REL",
    "FOS",
    "FOSL1",
    "FOSL2",
    "JUN",
    "JUNB",
    "JUND",
    "ATF3",
    "CTCF",
    "CTCFL",
    "NLRC5",
    "MYC",
]


def motif_tf_token(name: str) -> str:
    return name.split(".", 1)[0] if name else ""


def matches_target(name: str, patterns: list[tuple[str, re.Pattern]]) -> set[str]:
    found: set[str] = set()
    for sym, pat in patterns:
        if pat.search(name):
            found.add(sym)
    return found


def main() -> None:
    meme_path = Path(sys.argv[1]) if len(sys.argv) > 1 else MEME_DEFAULT
    jaspar_path = Path(sys.argv[2]) if len(sys.argv) > 2 else JASPAR_DEFAULT

    patterns = [(sym, re.compile(sym, re.IGNORECASE)) for sym in TARGETS]

    meme_names: list[str] = []
    if meme_path.is_file():
        for line in meme_path.read_text(errors="replace").splitlines():
            if line.startswith("MOTIF "):
                parts = line.split()
                if len(parts) >= 2:
                    meme_names.append(parts[1])

    jaspar_names: list[str] = []
    if jaspar_path.is_file():
        for line in jaspar_path.read_text(errors="replace").splitlines():
            if line.startswith(">"):
                tok = line[1:].strip().split()
                if tok:
                    jaspar_names.append(tok[0])

    meme_tokens = {motif_tf_token(n) for n in meme_names}
    jaspar_tokens = {motif_tf_token(n) for n in jaspar_names}

    meme_hit: set[str] = set()
    jaspar_hit: set[str] = set()
    for n in meme_names:
        meme_hit |= matches_target(n, patterns)
    for n in jaspar_names:
        jaspar_hit |= matches_target(n, patterns)

    missing_both = [s for s in TARGETS if s not in meme_hit and s not in jaspar_hit]
    only_meme = sorted(meme_hit - jaspar_hit)
    only_jaspar = sorted(jaspar_hit - meme_hit)
    both = sorted(meme_hit & jaspar_hit)

    print("=== File stats ===")
    print(f"MEME: {meme_path} — {len(meme_names)} MOTIF blocks, {len(meme_tokens)} unique name-prefix tokens")
    print(f"JASPAR-format: {jaspar_path} — {len(jaspar_names)} entries, {len(jaspar_tokens)} unique name-prefix tokens")

    print("\n=== SV_TARGET_TF_SYMBOLS (substring on full motif name) ===")
    print("In BOTH:", ", ".join(both) or "(none)")
    print("MEME only:", ", ".join(only_meme) or "(none)")
    print("JASPAR-format only:", ", ".join(only_jaspar) or "(none)")
    print("MISSING from both:", ", ".join(missing_both) or "(none)")

    print("\n=== Motif row counts per target ===")
    for sym in TARGETS:
        pat = re.compile(sym, re.IGNORECASE)
        cm = sum(1 for n in meme_names if pat.search(n))
        cj = sum(1 for n in jaspar_names if pat.search(n))
        if cm or cj:
            print(f"  {sym}: meme={cm}, jaspar_txt={cj}")


if __name__ == "__main__":
    main()
