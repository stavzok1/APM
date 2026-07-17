"""Evidence layer → the PMID-deduped, METHOD-centric ledger (Design §Decision D/E). IMPLEMENTED (Phase 2).

Refactor of build_edges + evidence_scoring (BUILD_PLAN §1b). `ledger.py` builds one ledger keyed by
(arm, gene, pmid, assay_class) from miRTarBase + TarBase v9, **de-duplicated globally by (edge × PMID ×
assay_class)** (physical AND functional), and fuses a **method-centric** per-edge weight
(Σ_class CLASS_WEIGHT · log1p(distinct PMIDs)). `methods.py` is the assay taxonomy + directness weights.

Built result (Hallmark scope): 1.92M raw (m,g,pmid,class) rows → 1.11M deduped (**42% collapsed**, the
same paper in both DBs); 6,313 distinct PMIDs; top fused edges are canonical (miR-21→PDCD4/PTEN,
miR-200c→ZEB1, miR-221→CDKN1B). Wired into learned.data.assemble_gene(w_prior_source='ledger').

Not yet: ENCORI/POSTAR3/Manakov CLIP folded into the same classes (union, PMID-deduped); learned method
weights vs an external repression label (replacing the hand-set CLASS_WEIGHT); K_D as the separate
sequence channel.
"""
from mirna_hallmark.learned.evidence.ledger import build_ledger, edge_weights  # noqa: F401
from mirna_hallmark.learned.evidence.methods import CLASS_WEIGHT, classify      # noqa: F401
