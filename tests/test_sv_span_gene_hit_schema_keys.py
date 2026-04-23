import unittest

import pandas as pd

from pipeline.SV.spatial_mapping import classify_span_sv_gene_hit


class TestSvSpanGeneHitSchemaKeys(unittest.TestCase):
    def test_classify_span_sv_gene_hit_keys_stable(self):
        """
        Guardrail: ``classify_span_sv_gene_hit`` output keys are part of the on-disk SV
        ``gene_hits`` schema contract. If you add/remove keys, update this test + Atlas.
        """
        gf = pd.DataFrame()
        d = classify_span_sv_gene_hit(
            sv_start=1000,
            sv_end=2000,
            gene_start=1500,
            gene_end=2500,
            strand="+",
            signed_dist=0,
            gene_features=gf,
        )
        expected = frozenset(
            {
                "overlap_bp",
                "overlap_percent",
                "promoter_flag",
                "gene_body_flag",
                "mane_cds_flag",
                "cds_flag",
                "utr_flag",
                "exon_flag",
                "mane_transcript_flag",
                "intron_only_flag",
                "upstream_5kb_flag",
                "downstream_5kb_flag",
                "region_hit",
                "hit_side",
                "stop_codon_flag",
                "start_codon_flag",
                "transcript_id",
                "transcript_type",
                "exon_interval_ids",
                "exon_ids",
                "cds_interval_ids",
                "utr_interval_ids",
                "start_codon_interval_ids",
                "stop_codon_interval_ids",
            }
        )
        self.assertEqual(frozenset(d.keys()), expected)
