# Gene symbol normalization & alias registry

This repo uses **canonical HGNC/GENCODE-style symbols** (as listed in `pipeline/config.py`) as the join keys for
gene-centric evidence. External sources frequently use alternate spellings, aliases, or Ensembl IDs; this document
describes the **single, pipeline-wide** normalization mechanism.

---

## What is normalized (and where)

The same alias registry is used across modules that ingest external gene tokens, including (non-exhaustive):
- RNA expression normalization (`pipeline/main.py` Step 11)
- SNV / SV / CNV annotation helpers (VEP-derived symbols, external tables)
- miRTarBase / TargetScan gene tokens
- methylation probe gene lists

If a module accepts `APM_USE_GENE_SYMBOL_MAPPING=0`, it will skip the heavy alias registry and only apply
small manual fixes (see below).

---

## Source of truth: the alias registry

**Code**: `pipeline/genes/panel_alias_registry.py`

It produces two maps:
- `get_gene_flat_alias_map(max_lines=None)`: **`alias → canonical`**
- `get_gene_canonical_to_aliases(max_lines=None)`: **`canonical → {aliases}`**

You usually won’t call these directly; modules tend to use:

**Helpers**: `pipeline/genes/symbol_normalization.py`
- `default_symbol_mapping()` → full registry (slow once; cached)
- `get_panel_symbol_mapping(max_lines)` → optionally capped for dev speed

---

## Where aliases come from (merged layers)

The registry merges these sources **in order**, later layers overriding earlier keys when needed:

1. **HGNC alias table**: `annotations/Alias_v5.22.xls`
   - Header is tab-separated; body lines are comma-separated CSV rows.
   - Only rows whose **approved symbol** is in `PANEL_ALIAS_SEED_SYMBOLS` are kept.
2. **UCSC fixes**: `UCSC_RNA_SEQ_GENE_SYMBOL_CHANGE` in `pipeline/config.py`
   - Small curated remaps seen in UCSC/TCGA builds.
3. **Legacy dataset remaps**: `LEGACY_DATASET_SYMBOL_RENAMES` in `pipeline/config.py`
   - Example: `TMEM173 → STING1`, `MB21D1 → CGAS`.
4. **Hand-curated canonical alias lists**: `MANUAL_CANONICAL_TO_ALIASES` in `panel_alias_registry.py`
   - Useful for “nickname” biology tokens (e.g. STING synonyms) that sometimes aren’t stable in HGNC exports.
5. **Ensembl gene ids**: `ENSG…` (+ version-stripped) → canonical `gene_name`
   - Derived from `PATHS.gencode_gtf_pq` for every symbol in `PANEL_ALIAS_SEED_SYMBOLS`.
   - This is why Ensembl-keyed matrices can still join to canonical symbols.

---

## Performance knob: `APM_HGNC_ALIAS_MAX_LINES`

`pipeline/main.py` (RNA Step 11) supports:
- `APM_HGNC_ALIAS_MAX_LINES=<N>`: limit the number of HGNC body lines scanned.

This is **only** a dev-speed feature. It can **miss** aliases that appear later in the HGNC table.
If you see unexpected unmapped tokens (e.g. `NF-YA`, `p65`, `GCN5`, …), rerun with the variable unset so the
registry performs a full scan (cached).

---

## Minimal example (interactive sanity check)

From repo root:

```bash
PYTHONPATH=/home/stavz/masters/gdc/APM .venv/bin/python3 -c "
from pipeline.genes.panel_alias_registry import get_gene_flat_alias_map
m = get_gene_flat_alias_map(None)
print('TMEM173 ->', m.get('TMEM173'))
print('MB21D1  ->', m.get('MB21D1'))
print('p65    ->', m.get('p65'))
"
```

---

## Disabling mapping

Many loaders respect:
- `APM_USE_GENE_SYMBOL_MAPPING=0`

When disabled, pipelines should not be expected to reconcile broad HGNC alias space; only small manual maps may apply.

