import json
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


SNAPSHOT_PATH = REPO / "tests" / "fixtures" / "pipeline_schema_snapshot.json"


class TestPipelineSchemaSnapshot(unittest.TestCase):
    def test_schema_snapshot_matches(self):
        if not SNAPSHOT_PATH.exists():
            self.fail(
                f"Missing schema snapshot fixture at {SNAPSHOT_PATH}. "
                f"Generate it via: .venv/bin/python3 scripts/schema/dump_pipeline_schema_snapshot.py"
            )

        from scripts.schema.dump_pipeline_schema_snapshot import build_snapshot

        expected = build_snapshot()
        got = json.loads(SNAPSHOT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(got, expected)


if __name__ == "__main__":
    unittest.main()

