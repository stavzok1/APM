import unittest

import pandas as pd

from pipeline.tad_annotation.mirroring import mirror_hits_into_boundaries_by_overlap


class TestTadMirroringBoundariesNonePayload(unittest.TestCase):
    def test_ccre_hit_extractor_handles_none_payload(self):
        boundaries = pd.DataFrame(
            {
                "boundary_id": ["b1"],
                "chrom": ["chr1"],
                "start": [100],
                "end": [200],
            }
        )
        ccres = pd.DataFrame(
            {
                "cCRE_id": ["c1"],
                "chrom": ["chr1"],
                "start": [150],
                "end": [160],
            }
        )

        out = mirror_hits_into_boundaries_by_overlap(
            boundaries,
            ccres,
            biosample="X",
            feature_kind="ccre",
        )
        self.assertIn("cCRE_hits", out.columns)
        cell = out.loc[0, "cCRE_hits"]
        self.assertIsInstance(cell, dict)
        self.assertIn("X", cell)
        self.assertEqual(len(cell["X"]), 1)
        self.assertEqual(cell["X"][0].get("cCRE_id"), "c1")
        # rel should exist and be None when payload missing
        self.assertTrue("rel" in cell["X"][0])
        self.assertIsNone(cell["X"][0].get("rel"))

