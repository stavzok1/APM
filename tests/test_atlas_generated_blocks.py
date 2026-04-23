import re
import unittest
from pathlib import Path

from pipeline.rppa.rppa_schemas import empty_rppa_panel_scores


REPO = Path(__file__).resolve().parents[1]
ATLAS = REPO / "pipeline" / "md" / "DATA_STRUCTURE_ATLAS.md"


def _extract_generated_block(text: str, key: str) -> str:
    begin = f"<!-- GENERATED:{key} BEGIN -->"
    end = f"<!-- GENERATED:{key} END -->"
    m = re.search(re.escape(begin) + r"(.*?)" + re.escape(end), text, flags=re.DOTALL)
    if not m:
        raise AssertionError(f"Missing generated block {key!r}")
    return m.group(1).strip()


class TestAtlasGeneratedBlocks(unittest.TestCase):
    def test_rppa_panel_names_match_schema(self):
        txt = ATLAS.read_text(encoding="utf-8")
        body = _extract_generated_block(txt, "rppa_panel_score_columns")
        expected = ", ".join(empty_rppa_panel_scores().keys())
        self.assertEqual(body, expected)
