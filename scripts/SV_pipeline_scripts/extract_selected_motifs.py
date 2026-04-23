import re
import sys
from pathlib import Path

# -----------------------------
# Config: which TFs to keep
# -----------------------------
TARGET_TF_SYMBOLS = [
    # STATs
    "STAT1", "STAT2", "STAT3",
    # IRFs
    "IRF1", "IRF2", "IRF3",
    # NF-kB family
    "NFKB1", "NFKB2", "RELA", "p65", "RELB", "REL",
    # AP-1-ish
    "FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF3",
    # CTCF family
    "CTCF", "CTCFL",
    # NLRC5 (likely absent, but we try)
    "NLRC5", "MYC"
]

# We’ll match these as substrings of the motif name, case-insensitive
TARGET_PATTERNS = [re.compile(sym, re.IGNORECASE) for sym in TARGET_TF_SYMBOLS]


def motif_matches_any_target(motif_name: str) -> bool:
    return any(pat.search(motif_name) for pat in TARGET_PATTERNS)


def extract_selected_motifs(in_meme: Path, out_meme: Path):
    """
    Read a MEME-format motif file and write a new MEME file containing only
    motifs whose name matches any of the TARGET_TF_SYMBOLS (substring match).
    """

    in_text = in_meme.read_text()

    # Split header vs motifs. In MEME v4 files, header ends before first "MOTIF ".
    parts = in_text.split("\nMOTIF ")
    if len(parts) == 1:
        raise ValueError("Could not find any 'MOTIF ' blocks in the input MEME file.")

    header = parts[0]
    motif_blocks_raw = parts[1:]

    selected_blocks = []
    found_symbols = set()

    for block_raw in motif_blocks_raw:
        # Re-add the leading 'MOTIF ' that was removed by split
        block = "MOTIF " + block_raw

        # First line of the block should look like: "MOTIF <name> [optional stuff]"
        first_line = block.splitlines()[0]
        # Example motif line in HOCOMOCO: "MOTIF STAT1_HUMAN.H14CORE.0.A"
        # We'll take token after "MOTIF" as the motif name.
        tokens = first_line.split()
        if len(tokens) < 2:
            continue
        motif_name = tokens[1]

        if motif_matches_any_target(motif_name):
            selected_blocks.append(block)
            # Try to detect which symbol matched (for reporting)
            for sym in TARGET_TF_SYMBOLS:
                if re.search(sym, motif_name, re.IGNORECASE):
                    found_symbols.add(sym)

    # Build output text
    out_lines = [header.rstrip()]  # keep header as-is
    out_lines.append("")  # blank line

    for block in selected_blocks:
        out_lines.append(block.rstrip())
        out_lines.append("")  # blank line between motifs

    out_text = "\n".join(out_lines).rstrip() + "\n"
    out_meme.write_text(out_text)

    # Simple report to stdout
    missing = [sym for sym in TARGET_TF_SYMBOLS if sym not in found_symbols]
    print(f"Selected {len(selected_blocks)} motif blocks into {out_meme}")
    if missing:
        print("No motif found (by name substring) for:", ", ".join(missing))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_selected_motifs.py <input_meme> <output_meme>")
        sys.exit(1)

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    extract_selected_motifs(in_path, out_path)
