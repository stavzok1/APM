import sys
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


class TestCovariatesModule(unittest.TestCase):
    def test_default_gene_sets_contains_expected(self):
        from pipeline.RNA_exp.signatures import default_gene_sets

        gs = default_gene_sets()
        self.assertIn("HYPOXIA_BUFFA_mean", gs)
        self.assertIn("AUTOPHAGY_core_mean", gs)
        self.assertIn("MHCII_readout_mean", gs)

    def test_tp53_classifier_empty_maf(self):
        from pipeline.covariates.genomics import classify_tp53_from_maf

        maf = pd.DataFrame(columns=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"])
        out = classify_tp53_from_maf(maf)
        self.assertIsInstance(out, pd.DataFrame)
        self.assertEqual(len(out), 0)

    def test_sample_to_sample_vial_lift(self):
        from pipeline.covariates.build_covariates import _lift_sample_to_sample_vial

        df_sample = pd.DataFrame(
            {"value": [1.0, 2.0]},
            index=pd.Index(["TCGA-AA-0001-01", "TCGA-AA-0002-01"], name="sample"),
        )
        target_vials = pd.Index(["TCGA-AA-0001-01A", "TCGA-AA-0002-01A"], name="sample_vial")
        lifted = _lift_sample_to_sample_vial(df_sample, target_vials=target_vials)
        self.assertEqual(list(lifted.index), list(target_vials))
        self.assertEqual(lifted.loc["TCGA-AA-0001-01A", "value"], 1.0)
        self.assertEqual(lifted.loc["TCGA-AA-0002-01A", "value"], 2.0)


if __name__ == "__main__":
    unittest.main()

