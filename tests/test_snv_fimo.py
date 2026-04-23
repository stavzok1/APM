"""Unit tests for SNV reference-window FIMO helpers (no ``fimo`` binary required)."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.SNV.snv_fimo import (
    build_snv_fimo_bed,
    variant_reference_window_0based,
)
from pipeline.SNV.snv_fimo import aggregate_fimo_hits_per_variant


class TestVariantWindow(unittest.TestCase):
    def test_centered_span(self) -> None:
        s0, e = variant_reference_window_0based(1000, flank_bp=2)
        self.assertEqual(s0, 997)
        self.assertEqual(e, 1002)
        self.assertEqual(e - s0, 5)

    def test_clamp_start_at_zero(self) -> None:
        s0, e = variant_reference_window_0based(2, flank_bp=10)
        self.assertEqual(s0, 0)


class TestBuildBed(unittest.TestCase):
    def test_names_are_row_indices(self) -> None:
        df = pd.DataFrame(
            {
                "chrom": ["chr1", "chrom2"],
                "pos": [100, 200],
                "ref": ["A", "C"],
                "alt": ["G", "T"],
            }
        )
        with tempfile.TemporaryDirectory() as td:
            bed = Path(td) / "w.bed"
            build_snv_fimo_bed(df, bed, flank_bp=1)
            lines = bed.read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].endswith("\t0"))
        self.assertTrue(lines[1].endswith("\t1"))


class TestHitsFromFimoDf(unittest.TestCase):
    def test_genomic_lift_and_cap(self) -> None:
        fimo_df = pd.DataFrame(
            [
                {
                    "sequence_name": "0",
                    "motif_id": "CTCF.MA0139.1",
                    "start": 1,
                    "stop": 10,
                    "strand": "+",
                    "p-value": 1e-5,
                    "q-value": 0.01,
                    "score": 12.3,
                    "matched_sequence": "ACGTACGTAC",
                },
                {
                    "sequence_name": "0",
                    "motif_id": "JUN.MA0489.1",
                    "start": 2,
                    "stop": 8,
                    "strand": "-",
                    "p-value": 1e-3,
                    "q-value": 0.02,
                    "score": 5.0,
                    "matched_sequence": "AAAAAA",
                },
            ]
        )
        # Window for POS=100, flank=10 -> start0 = 89
        window_starts0 = [89]
        variant_positions = [100]
        hits = aggregate_fimo_hits_per_variant(
            fimo_df, window_starts0, variant_positions, max_hits_per_variant=1
        )
        self.assertEqual(len(hits), 1)
        self.assertEqual(len(hits[0]), 1)
        self.assertEqual(hits[0][0]["TF"], "CTCF")
        self.assertEqual(hits[0][0]["start"], 89)
        self.assertEqual(hits[0][0]["end"], 99)
        self.assertEqual(hits[0][0]["variant_pos"], 100)


if __name__ == "__main__":
    unittest.main()
