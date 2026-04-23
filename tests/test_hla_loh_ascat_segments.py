import unittest

import pandas as pd

from pipeline.qc.hla_loh_ascat_segments import hla_loh_from_segments_for_genes


class TestHlaLohAscatSegments(unittest.TestCase):
    def test_copy_neutral_loh_minor_zero(self):
        segs = pd.DataFrame(
            {
                "chrom": ["chr6"],
                "start": pd.array([29941200], dtype="Int64"),
                "end": pd.array([29946000], dtype="Int64"),
                "cn_total": [2.0],
                "cn_major": [2.0],
                "cn_minor": [0.0],
            }
        )
        intervals = {"HLA-A": ("chr6", 29941260, 29945884)}
        rec = hla_loh_from_segments_for_genes(segs, intervals, min_overlap_frac=0.30)
        self.assertTrue(rec["HLA-A_seg_loh"])
        self.assertTrue(rec["hla_loh_seg_any"])
