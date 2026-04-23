"""Tests for SNV strict-overlap ChIP annotation."""

from __future__ import annotations

import unittest
from pathlib import Path

import pandas as pd

from pipeline.SNV.snv_chip import annotate_snvs_with_chip


class TestSnvChipStrictOverlap(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = Path(__file__).resolve().parent / "_tmp_snv_chip.parquet"
        chip = pd.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr1"],
                "start": [100, 103, 200],
                "end": [103, 200, 250],
                "tf": ["CTCF", "CTCF", "STAT3"],
                "cell_type": ["MCF7", "T47D", "MCF7"],
                "source": ["ENCODE", "ENCODE", "CHIP_ATLAS"],
                "score_norm": [1.0, 2.0, 3.0],
                "sample_id": ["CTCF_a", "CTCF_b", "STAT3_x"],
                "cell_subtype": ["", "Basal_like", ""],
            }
        )
        chip.to_parquet(self.tmp, index=False)

    def tearDown(self) -> None:
        if self.tmp.is_file():
            self.tmp.unlink()

    def test_pos104_overlaps_second_peak_span(self) -> None:
        """POS 104 (0-based point 103) lies in [103,200) only, not in [100,103)."""
        df = pd.DataFrame({"chrom": ["chr1"], "pos": [104]})
        out = annotate_snvs_with_chip(df, self.tmp)
        hits = out.iloc[0]["snv_chip_hits"]
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]["peak_start"], 103)
        self.assertEqual(hits[0]["tf"], "CTCF")

    def test_pos103_overlaps_first_peak(self) -> None:
        """POS 103 overlaps [100,103): last base 0-based 102; peak end exclusive 103 covers 102."""
        df = pd.DataFrame({"chrom": ["chr1"], "pos": [103]})
        out = annotate_snvs_with_chip(df, self.tmp)
        hits = out.iloc[0]["snv_chip_hits"]
        tfs = {h["peak_start"] for h in hits}
        self.assertIn(100, tfs)

    def test_aggregate_nonempty(self) -> None:
        df = pd.DataFrame({"chrom": ["chr1"], "pos": [104]})
        out = annotate_snvs_with_chip(df, self.tmp)
        agg = out.iloc[0]["snv_chip_aggregate"]
        self.assertIn("by_tf_source_stratum", agg)
        self.assertTrue(len(agg["by_tf_source_stratum"]) >= 1)


if __name__ == "__main__":
    unittest.main()
