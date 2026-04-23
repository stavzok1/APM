import unittest

import pandas as pd

from pipeline.config import PIPELINE_GENE_PANEL, TIER1_LNCRNA_GENES
from pipeline.scanning_columns import derive_elem_focus_scanning_columns


class TestScanningColumnsSpec(unittest.TestCase):
    def test_scan_gene_links_n_genes_matches_len_gene_links(self):
        panel = list(PIPELINE_GENE_PANEL)[:3]
        tier1 = list(TIER1_LNCRNA_GENES)[:1]
        gene_links = {g: {"k": 1} for g in panel + tier1}
        df = pd.DataFrame(
            {
                "gene_links": [gene_links],
                "screen_exp": [{}],
                "ABC_enhancers": [[]],
                "hichip": [{}],
                "chip_hits": [{}],
            }
        )
        out = derive_elem_focus_scanning_columns(df)
        self.assertEqual(int(out["scan_gene_links_n_genes"].iloc[0]), len(gene_links))
        self.assertEqual(
            bool(out["scan_gene_links_any_panel_gene"].iloc[0]),
            any(g in set(PIPELINE_GENE_PANEL) for g in gene_links),
        )
