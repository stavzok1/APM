#!/usr/bin/env python3
"""
Scan all HOCOMOCO v14 MEME bundles in a directory for motif name hits.

Compares full ``MOTIF <name>`` / ``><name>`` lines against a panel of TF labels
(case-insensitive substring or regex per row).

Example::

    .venv/bin/python3 scripts/motifs/hocomoco_multi_bundle_tf_report.py /home/stavz/HOCOMOCO
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


# Pipeline default panel + breast / MHC-II / chromatin / NF-Y / TGF-β / hormone extras
QUERIES: list[tuple[str, re.Pattern]] = [
    ("STAT1", re.compile(r"STAT1", re.I)),
    ("STAT2", re.compile(r"STAT2", re.I)),
    ("STAT3", re.compile(r"STAT3", re.I)),
    ("IRF1", re.compile(r"IRF1", re.I)),
    ("IRF2", re.compile(r"IRF2", re.I)),
    ("IRF3", re.compile(r"IRF3", re.I)),
    ("IRF9", re.compile(r"IRF9", re.I)),
    ("NFKB1", re.compile(r"NFKB1", re.I)),
    ("NFKB2", re.compile(r"NFKB2", re.I)),
    ("RELA", re.compile(r"RELA", re.I)),
    ("p65_RELA_alias", re.compile(r"p65", re.I)),
    ("RELB", re.compile(r"RELB", re.I)),
    ("REL", re.compile(r"\bREL\.", re.I)),  # avoid RELB / RELA false REL. token
    ("FOS", re.compile(r"\bFOS\.", re.I)),
    ("FOSL1", re.compile(r"FOSL1", re.I)),
    ("FOSL2", re.compile(r"FOSL2", re.I)),
    ("JUN", re.compile(r"\bJUN\.", re.I)),
    ("JUNB", re.compile(r"JUNB", re.I)),
    ("JUND", re.compile(r"JUND", re.I)),
    ("ATF3", re.compile(r"ATF3", re.I)),
    ("CTCF", re.compile(r"CTCF", re.I)),
    ("CTCFL", re.compile(r"CTCFL", re.I)),
    ("NLRC5", re.compile(r"NLRC5", re.I)),
    ("MYC", re.compile(r"\bMYC\.", re.I)),
    ("CIITA", re.compile(r"CIITA", re.I)),
    ("RFX5", re.compile(r"RFX5", re.I)),
    ("NFYA", re.compile(r"NFYA", re.I)),
    ("NFYB", re.compile(r"NFYB", re.I)),
    ("NFYC", re.compile(r"NFYC", re.I)),
    ("EZH2", re.compile(r"EZH2", re.I)),
    ("SUZ12", re.compile(r"SUZ12", re.I)),
    ("SMAD_family", re.compile(r"SMAD", re.I)),
    ("BRD4", re.compile(r"BRD4", re.I)),
    ("FOXA1", re.compile(r"FOXA1", re.I)),
    ("ESR1", re.compile(r"ESR1", re.I)),
    ("GATA3", re.compile(r"GATA3", re.I)),
]


def parse_meme_motif_names(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("MOTIF "):
            parts = line.split()
            if len(parts) >= 2:
                names.append(parts[1])
    return names


def parse_jaspar_txt_headers(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith(">"):
            tok = line[1:].strip().split()
            if tok:
                names.append(tok[0])
    return names


def count_hits(names: list[str], pat: re.Pattern) -> int:
    return sum(1 for n in names if pat.search(n))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dir", type=Path, nargs="?", default=Path("/home/stavz/HOCOMOCO"))
    args = ap.parse_args()
    root: Path = args.dir
    if not root.is_dir():
        print(f"Not a directory: {root}", file=sys.stderr)
        sys.exit(1)

    meme_files = sorted(p for p in root.glob("H14*_meme_format.meme") if p.is_file())
    jaspar_txts = sorted(p for p in root.glob("H14*_jaspar_format.txt") if p.is_file())

    if not meme_files:
        print(f"No H14*_meme_format.meme under {root}", file=sys.stderr)
        sys.exit(1)

    bundle_data: dict[str, list[str]] = {}
    for p in meme_files:
        bundle_data[p.name] = parse_meme_motif_names(p)
    for p in jaspar_txts:
        key = p.name + " (jaspar headers)"
        bundle_data[key] = parse_jaspar_txt_headers(p)

    print(f"# HOCOMOCO multi-bundle TF scan\n\nRoot: `{root}`\n")
    print("| Bundle | motifs / headers |")
    print("|---|---:|")
    for k, names in sorted(bundle_data.items()):
        print(f"| {k} | {len(names)} |")
    print()

    # Per-query: which bundles have >=1 hit
    print("## Coverage by TF (any hit in bundle name string)\n")
    print("| TF / family | " + " | ".join(sorted(bundle_data.keys())) + " |")
    print("|" + "---|" * (len(bundle_data) + 1))

    for label, pat in QUERIES:
        cells = []
        for bname in sorted(bundle_data.keys()):
            c = count_hits(bundle_data[bname], pat)
            cells.append(str(c) if c else "—")
        print(f"| {label} | " + " | ".join(cells) + " |")

    # Summary: missing from ALL meme files (not jaspar - same names typically)
    meme_only = {k: v for k, v in bundle_data.items() if k.endswith(".meme")}
    print("\n## Missing from **all** MEME bundles (0 in every `*_meme_format.meme`)\n")
    missing = []
    for label, pat in QUERIES:
        if all(count_hits(v, pat) == 0 for v in meme_only.values()):
            missing.append(label)
    print(", ".join(missing) if missing else "(none)")

    print(
        "\n---\n"
        "**Note:** `REL` uses pattern `REL.` to reduce false hits on `RELA`/`RELB`. "
        "`FOS` / `JUN` / `MYC` use a dot after the token to match HOCOMOCO names like `FOS.H14...`."
    )


if __name__ == "__main__":
    main()
