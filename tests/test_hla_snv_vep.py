import unittest

import pandas as pd

from pipeline.qc.hla_snv_vep import summarize_hla_snv_from_df


class TestHlaSnvVep(unittest.TestCase):
    def test_lof_on_hla_transcript_detected(self):
        df = pd.DataFrame(
            {
                "gene_hits": [
                    [
                        {
                            "Feature_type": "Transcript",
                            "SYMBOL": "HLA-A",
                            "Consequence": "stop_gained",
                            "IMPACT": "HIGH",
                            "CLIN_SIG": None,
                            "gnomADe_AF": "0.0001",
                        }
                    ]
                ],
                "tumor_vaf": [0.4],
            }
        )
        s = summarize_hla_snv_from_df(df)
        self.assertTrue(s["hla_snv_lof_any"])
        self.assertTrue(s["hla_snv_lof_rare_any"])
        self.assertEqual(s["hla_snv_n_lof_variants"], 1)
        self.assertEqual(s["hla_snv_genes_hit"], "HLA-A")

    def test_high_missense_separate_from_lof(self):
        df = pd.DataFrame(
            {
                "gene_hits": [
                    [
                        {
                            "Feature_type": "Transcript",
                            "SYMBOL": "B2M",
                            "Consequence": "missense_variant",
                            "IMPACT": "HIGH",
                        }
                    ]
                ],
                "tumor_vaf": [0.1],
            }
        )
        s = summarize_hla_snv_from_df(df)
        self.assertFalse(s["hla_snv_lof_any"])
        self.assertTrue(s["hla_snv_high_missense_any"])
        # No pathogenicity scores on hit ⇒ HIGH missense still counts toward damage.
        self.assertTrue(s["hla_snv_high_missense_damage_any"])

    def test_lof_common_pop_excluded_from_damage(self):
        df = pd.DataFrame(
            {
                "gene_hits": [
                    [
                        {
                            "Feature_type": "Transcript",
                            "SYMBOL": "HLA-B",
                            "Consequence": "stop_gained",
                            "IMPACT": "HIGH",
                            "gnomADe_AF": "0.5",
                        }
                    ]
                ],
                "tumor_vaf": [0.2],
            }
        )
        s = summarize_hla_snv_from_df(df)
        self.assertTrue(s["hla_snv_lof_any"])
        self.assertFalse(s["hla_snv_lof_rare_any"])

    def test_high_missense_with_low_cadd_excluded(self):
        df = pd.DataFrame(
            {
                "gene_hits": [
                    [
                        {
                            "Feature_type": "Transcript",
                            "SYMBOL": "B2M",
                            "Consequence": "missense_variant",
                            "IMPACT": "HIGH",
                            "CADD_PHRED": "8",
                            "gnomADe_AF": "0.0001",
                        }
                    ]
                ],
                "tumor_vaf": [0.1],
            }
        )
        s = summarize_hla_snv_from_df(df)
        self.assertTrue(s["hla_snv_high_missense_any"])
        self.assertFalse(s["hla_snv_high_missense_damage_any"])

    def test_high_missense_cadd_passes(self):
        df = pd.DataFrame(
            {
                "gene_hits": [
                    [
                        {
                            "Feature_type": "Transcript",
                            "SYMBOL": "B2M",
                            "Consequence": "missense_variant",
                            "IMPACT": "HIGH",
                            "CADD_PHRED": "28",
                            "gnomADe_AF": "0.0001",
                        }
                    ]
                ],
                "tumor_vaf": [0.1],
            }
        )
        s = summarize_hla_snv_from_df(df)
        self.assertTrue(s["hla_snv_high_missense_damage_any"])
