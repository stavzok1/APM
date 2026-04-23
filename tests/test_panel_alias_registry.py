"""Unit tests for canonical alias registries (genes + miRNA)."""

import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


class TestPanelAliasRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from pipeline.genes.panel_alias_registry import clear_gene_alias_caches

        clear_gene_alias_caches()

    def test_manual_gene_alias_pvrl2(self):
        from pipeline.genes.panel_alias_registry import get_gene_flat_alias_map

        m = get_gene_flat_alias_map(400_000)
        self.assertEqual(m.get("PVRL2"), "NECTIN2")

    def test_canonical_inverse_contains_self(self):
        from pipeline.genes.panel_alias_registry import get_gene_canonical_to_aliases

        inv = get_gene_canonical_to_aliases(400_000)
        self.assertIn("NECTIN2", inv)
        self.assertIn("NECTIN2", inv["NECTIN2"])

    def test_mirna_normalize(self):
        from pipeline.genes.panel_alias_registry import normalize_mirna_symbol

        self.assertEqual(normalize_mirna_symbol("miR-155-5p"), "hsa-miR-155-5p")

    def test_resolve_symbol_to_panel(self):
        from pipeline.genes.symbol_normalization import resolve_symbol_to_panel

        m = {"PVRL2": "NECTIN2"}
        self.assertEqual(resolve_symbol_to_panel("PVRL2", {"NECTIN2"}, m), "NECTIN2")

    def test_resolve_ensembl_id_to_panel(self):
        from pipeline.genes.symbol_normalization import resolve_symbol_to_panel

        m = {"ENSG00000166710": "B2M", "ENSG00000166710.24": "B2M"}
        self.assertEqual(resolve_symbol_to_panel("ENSG00000166710", {"B2M"}, m), "B2M")
        self.assertEqual(resolve_symbol_to_panel("ENSG00000166710.24", {"B2M"}, m), "B2M")
        m2 = {"ENSG00000166710": "B2M"}
        self.assertEqual(resolve_symbol_to_panel("ENSG00000166710.24", {"B2M"}, m2), "B2M")

    def test_b2m_ensembl_in_flat_map(self):
        from pipeline.genes.panel_alias_registry import get_gene_flat_alias_map

        m = get_gene_flat_alias_map(50_000)
        self.assertEqual(m.get("ENSG00000166710"), "B2M")

    def test_mirna_mimat_resolves_to_arm(self):
        from pipeline.genes.panel_alias_registry import normalize_mirna_symbol

        out = normalize_mirna_symbol("MIMAT0004571")
        self.assertTrue(out.startswith("hsa-miR-"))


class TestRnaMappingGate(unittest.TestCase):
    def test_ucsc_only_mapping_contains_pvrl2(self):
        from pipeline.config import UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE

        self.assertEqual(UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE.get("PVRL2"), "NECTIN2")


if __name__ == "__main__":
    unittest.main()
