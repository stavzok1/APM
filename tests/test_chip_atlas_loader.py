"""Unit tests for ChIP-Atlas narrow BED loading (UCSC conversion in ``chip_atlas_convert``)."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pipeline.CHIP.chip_atlas_convert import (
    chip_atlas_gfftags_cell_line,
    convert_chip_atlas_file_to_narrow,
    needs_ucsc_conversion,
)
from pipeline.CHIP.chip_loader import (
    CHIP_ATLAS_UNSPECIFIED_CELL_LINE,
    load_chip_atlas_bed,
)


class TestChipAtlasGffTags(unittest.TestCase):
    def test_cell_line_from_encoded_tags(self) -> None:
        meta = (
            "ID=SRX5213646;Name=BRD4%20(@%20HCC1937);Title=...;"
            "cell%20line=HCC1937;chip%20antibody=BRD4"
        )
        self.assertEqual(chip_atlas_gfftags_cell_line(meta), "HCC1937")

    def test_cell_from_name_only(self) -> None:
        meta = (
            "ID=ERX371713;Name=BRD4%20(@%20Mammary%20epithelial%20cells);"
            "Title=Illumina%20HiSeq%202000%20sequencing;Cell%20group=Breast;"
        )
        self.assertEqual(
            chip_atlas_gfftags_cell_line(meta),
            "Mammary epithelial cells",
        )

    def test_source_name_fallback(self) -> None:
        meta = (
            "ID=SRX14612984;Name=JUN%20(@%20T-47D);Cell%20group=Breast;"
            "<br>source_name=T47D;cell%20type=Mammary%20gland%20ductal%20carcinoma"
        )
        self.assertEqual(chip_atlas_gfftags_cell_line(meta), "T47D")


class TestChipAtlasConvertThenLoad(unittest.TestCase):
    def test_ucsc_bed6_gfftags_roundtrip(self) -> None:
        content = (
            'track name="BRD4 (@ Breast) 50"\n'
            "chr1\t9905\t10457\t"
            "ID=SRX5213646;Name=BRD4%20(@%20HCC1937);cell%20line=HCC1937;\t581\t.\n"
        )
        with tempfile.TemporaryDirectory() as td:
            src = Path(td) / "src_ucsc.bed"
            narrow = Path(td) / "BRD4.bed"
            src.write_text(content, encoding="utf-8")
            self.assertTrue(needs_ucsc_conversion(src))
            convert_chip_atlas_file_to_narrow(src, narrow, tf_from_filename="BRD4")
            self.assertFalse(needs_ucsc_conversion(narrow))
            df = load_chip_atlas_bed(narrow)
            self.assertEqual(len(df), 1)
            self.assertEqual(df.iloc[0]["tf"], "BRD4")
            self.assertEqual(df.iloc[0]["cell_type"], "HCC1937")


class TestLoadChipAtlasBed(unittest.TestCase):
    def test_legacy_six_column(self) -> None:
        content = (
            "chrom\tstart\tend\texperiment\ttf\tcell_type\n"
            "chr1\t100\t200\tSRX1\tSTAT1\tMCF-7\n"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
            tmp.write(content)
            p = Path(tmp.name)
        try:
            df = load_chip_atlas_bed(p)
            self.assertEqual(len(df), 1)
            self.assertEqual(df.iloc[0]["tf"], "STAT1")
            self.assertEqual(df.iloc[0]["cell_type"], "MCF-7")
        finally:
            p.unlink(missing_ok=True)

    def test_legacy_empty_cell_line_gets_sentinel(self) -> None:
        content = "chr1\t100\t200\tSRX9\tSTAT3\t\n"
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
            tmp.write(content)
            p = Path(tmp.name)
        try:
            df = load_chip_atlas_bed(p)
            self.assertEqual(len(df), 1)
            self.assertEqual(df.iloc[0]["cell_type"], CHIP_ATLAS_UNSPECIFIED_CELL_LINE)
        finally:
            p.unlink(missing_ok=True)

    def test_ucsc_without_preprocess_raises(self) -> None:
        content = (
            "chr1\t9905\t10457\t"
            "ID=SRX1;Name=STAT3%20peakset;\t581\t.\n"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as tmp:
            tmp.write(content)
            p = Path(tmp.name)
        try:
            with self.assertRaises(ValueError):
                load_chip_atlas_bed(p)
        finally:
            p.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
