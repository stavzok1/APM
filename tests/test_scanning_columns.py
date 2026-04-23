import sys
import unittest
from pathlib import Path

import pandas as pd


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


class TestScanningColumns(unittest.TestCase):
    def test_elem_focus_scanning_columns_basic(self):
        from pipeline.scanning_columns import derive_elem_focus_scanning_columns

        df = pd.DataFrame(
            [
                {
                    "cCRE_id": "EH38E0000001.1",
                    "gene_links": {"HLA-A": {}, "B2M": {}},
                    "screen_exp": {
                        "per_biosample": {
                            "MCF7": {"Intact-HiC": {"score": 0.2, "strength": "strong"}},
                            "MCF10A": {"Intact-HiC": {"score": 0.1, "strength": "weak"}},
                        },
                        "conservation_breast": {"Intact-HiC": {"n_strong": 2}},
                    },
                    "ABC_enhancers": [
                        {
                            "start": 1,
                            "end": 2,
                            "ABC_full": {
                                "MCF7_ENCODE": {"ABC_score": 0.12, "is_strong": True, "is_self_promoter": False},
                                "MCF10A-Ji2017": {"ABC_score": 0.05, "is_strong": False, "is_self_promoter": True},
                            },
                        }
                    ],
                    "hichip": {"MCF7": {"n_loops": 3, "max_score": 1.2}},
                    "chip_hits": {"STAT1": {"MCF7": {"ENCODE": {"n_peaks": 1}}}, "CTCF": {}},
                }
            ]
        )
        out = derive_elem_focus_scanning_columns(df)
        self.assertEqual(out.loc[0, "scan_gene_links_n_genes"], 2)
        self.assertTrue(int(out.loc[0, "scan_screen_exp_n_biosamples_any"]) >= 1)
        self.assertEqual(out.loc[0, "scan_screen_exp_MCF7_max_score"], 0.2)
        self.assertEqual(out.loc[0, "scan_ABC_max_score_MCF7_ENCODE"], 0.12)
        self.assertEqual(out.loc[0, "scan_hichip_n_loops_MCF7"], 3)
        self.assertTrue(out.loc[0, "scan_chip_hits_has_STAT1"])
        self.assertTrue(out.loc[0, "scan_chip_hits_has_CTCF"])

    def test_atac_peak_scanning_columns_basic(self):
        from pipeline.scanning_columns import derive_atac_peak_scanning_columns

        df = pd.DataFrame(
            [
                {
                    "peak_id": "chr1:1-2",
                    "gene_links": {"HLA-A": {"dist_to_tss": 10, "tier": "0–100kb"}},
                    "lncrna_links": {"NEAT1": {"dist_to_tss": 20, "tier": "0–100kb"}},
                    "ccre_links": [{"cCRE_id": "EH38E0000001.1"}],
                    "TAD_domains": {"Kim_T47D": {"domains": {}, "primary": {}}},
                    "TAD_boundary_overlaps": {"Kim_T47D": {"overlaps_boundary": True}},
                }
            ]
        )
        out = derive_atac_peak_scanning_columns(df)
        self.assertEqual(out.loc[0, "scan_gene_links_n_genes"], 1)
        self.assertEqual(out.loc[0, "scan_ccre_links_n_ccres"], 1)
        self.assertEqual(out.loc[0, "scan_TAD_domains_n_biosamples"], 1)
        self.assertTrue(out.loc[0, "scan_TAD_boundary_overlaps_any"])


if __name__ == "__main__":
    unittest.main()

