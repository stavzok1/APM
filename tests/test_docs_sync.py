import re
import sys
import unittest
from pathlib import Path
from typing import Dict, List, Set


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


DOC = REPO / "research_plan" / "01_gene_panel_extended.md"


def _extract_generated_block(text: str, key: str) -> str:
    begin = f"<!-- GENERATED:{key} BEGIN -->"
    end = f"<!-- GENERATED:{key} END -->"
    if begin not in text or end not in text:
        raise AssertionError(f"Missing generated block markers for {key!r}")
    m = re.search(re.escape(begin) + r"(.*?)" + re.escape(end), text, flags=re.DOTALL)
    if not m:
        raise AssertionError(f"Could not extract generated block {key!r}")
    return m.group(1).strip()


def _extract_backticked_symbols(text: str) -> Set[str]:
    # Extremely permissive; we filter for all-caps genes + digits + hyphen.
    toks = set(re.findall(r"`([A-Za-z0-9][A-Za-z0-9\\-]*)`", text))
    # keep gene-ish tokens only (avoid counts etc.)
    out = {t for t in toks if re.fullmatch(r"[A-Z0-9][A-Z0-9\\-]*", t)}
    return out


class TestGenePanelDocSync(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        from pipeline import config as cfg

        cls.cfg = cfg
        cls.doc_text = DOC.read_text(encoding="utf-8")

    def test_generated_blocks_present(self):
        for key in [
            "panel_counts",
            "extended_primary_config_order",
            "tier2_medium_genes_config_order",
            "tier3_cnv_only_genes_config_order",
            "stratifier_genes_config_order",
            "tier4_readout_genes_config_order",
        ]:
            _extract_generated_block(self.doc_text, key)

    def test_counts_table_matches_config(self):
        from scripts.panel.regen_gene_panel_doc_fragments import build_blocks

        expected = build_blocks()["panel_counts"].strip()
        got = _extract_generated_block(self.doc_text, "panel_counts")
        self.assertEqual(got, expected)

    def test_config_order_blocks_match_config(self):
        from scripts.panel.regen_gene_panel_doc_fragments import build_blocks

        expected_blocks: Dict[str, str] = build_blocks()
        for key in [
            "extended_primary_config_order",
            "tier2_medium_genes_config_order",
            "tier3_cnv_only_genes_config_order",
            "stratifier_genes_config_order",
            "tier4_readout_genes_config_order",
        ]:
            got = _extract_generated_block(self.doc_text, key)
            self.assertEqual(got, expected_blocks[key].strip())

    def test_all_gene_symbols_mentioned_are_in_panel_alias_seed(self):
        mentioned = _extract_backticked_symbols(self.doc_text)

        # Allowlist: non-gene tokens we intentionally backtick in prose.
        allow = {
            "PIPELINE_GENE_PANEL",
            "PRIMARY_GENES",
            "EXTENDED_PRIMARY_GENES",
            "TIER1_LNCRNA_GENES",
            "TIER2_MEDIUM_GENES",
            "TIER3_CNV_ONLY_GENES",
            "TIER_SNV_CNV_STRATIFIER_GENES",
            "TIER4_READOUT_GENES",
            "FULL_INTEGRATION_GENES",
            "CNV_GENES",
            "RNA_EXPRESSION_ALIAS_SEED_SYMBOLS",
            "PANEL_ALIAS_SEED_SYMBOLS",
            # Legacy aliases that are intentionally mentioned in the doc but
            # are not panel symbols (they are remapped via LEGACY_DATASET_SYMBOL_RENAMES).
            "TMEM173",
            "MB21D1",
            "STING",
        }
        mentioned = {m for m in mentioned if m not in allow}

        panel = set(self.cfg.RNA_EXPRESSION_ALIAS_SEED_SYMBOLS)
        missing = sorted(mentioned - panel)
        self.assertFalse(
            missing,
            "Gene-like symbols mentioned in doc but not in panel alias seed set: "
            + ", ".join(missing),
        )


if __name__ == "__main__":
    unittest.main()

