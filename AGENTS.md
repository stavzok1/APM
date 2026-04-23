# Agent onboarding: APM regulatory pipeline

Read this **before** editing code or answering project-specific questions. It lists the documentation to ingest, the execution spine of the repo, and where common changes land.

---

## 1. What this repository is

**APM** is a Python **pipeline** under `pipeline/` that builds **regulatory element–centric** tables (ENCODE cCREs) with **multi-source evidence** (SCREEN, ABC, HiChIP, ChIP, optional omics modules), **gene / lncRNA** context, and optional **TAD** annotation. Outputs are mostly **pandas** tables and **nested dicts** (serialized in CSV/JSON columns).

Data roots and defaults live in **`pipeline/config.py`** (`PATHS`, `BIOSAMPLES`, `THRESHOLDS`, gene lists).

---

## 2. Documentation to read (in order)

### Tier A — always read for any pipeline task

| Order | Path | Why |
|------:|------|-----|
| 1 | `pipeline/md/PIPELINE_STRUCTURE_AND_USAGE.md` | Directory layout, **step-by-step flow**, CLI / function entry points, output locations. |
| 2 | `pipeline/md/DATA_STRUCTURE_ATLAS.md` | **Visual contract** for nested columns (`gene_links`, `chip_hits`, SCREEN blocks, ABC, etc.). Use the TOC; jump to the section you touch. |
| 3 | `pipeline/config.py` | **Single source** for paths, biosample panels, thresholds. Prefer extending here over hard-coding paths in modules. |

### Tier B — read when the task touches IDs, biosamples, or cutoffs

| Path | Why |
|------|-----|
| `pipeline/md/SAMPLE_ID_MATCHING_GUIDE.md` | TCGA-style **join keys** (`sample_vial`, etc.) and **canonical biosample** rules (`biosample_names`, subtype map pointer). |
| `pipeline/md/THRESHOLDS_AND_SIGNAL_STRENGTH.md` | Where **weak / strong** and other numeric cutoffs are defined and how they behave. |
| `pipeline/biosample_names.py` | Canonical **cell line / biosample** strings vs **on-disk** names; SCREEN/ABC slug maps. |
| `pipeline/cell_line_subtype_map.py` | **PAM50-style** tumor groups (LumA, LumB, HER2, Basal) and **healthy / Normal_like** labels for canonical keys. |

### Tier C — downstream analysis

| Path | Why |
|------|-----|
| `analysis/README.md` | **Sample × module coverage** script: counts and intersections across RNA/SNV/SV/CNV/ATAC/methylation/RPPA/HLA/immune using normalized TCGA keys. Regenerate RNA/ATAC/HLA column manifests with `scripts/annotations/build_annotation_sample_manifests.py`. |

### Tier D — module-deep dives (read only if your task names that area)

All under `pipeline/md/module_specific_processing_md/` unless noted (Tier D):

| Topic | Doc |
|-------|-----|
| SNV | `snv_documentation.md` |
| SV + motifs | `sv_→_motif_pipeline_vep_fimo.md` (canonical; `pipeline/SV/SV_PIPELINE_README.md` is a short pointer only) |
| Methylation | `METHYLATION_MODULE.md` |
| CNV | `pipeline/CNV/CNV_MODULE.md` |
| ATAC peaks | `ATAC_PEAKS_USAGE.md` |
| miRTarBase | `mirtarbase_module.md` |

**Historical / superseded:** `pipeline/md/prev_md/` — do not treat as current unless comparing behavior.

**Ignore for onboarding:** anything under `.venv/`, `old_pipelines/`, and duplicate `test/pipeline/` trees unless you are explicitly migrating legacy code.

---

## 3. Core Python spine (minimal mental model)

```
pipeline/main.py          → run_full_pipeline(), run_evidence_only(), …
pipeline/config.py        → PATHS, BIOSAMPLES, THRESHOLDS, PRIMARY_GENES (core 66), PIPELINE_GENE_PANEL, CNV_GENES
pipeline/schemas.py       → empty_*, ensure_* for nested evidence shapes
pipeline/utils.py         → shared geometry / chr / strength helpers

pipeline/genes/           → GENCODE load, promoters, lncRNA, miRNA / miRTarBase
pipeline/regulatory_elements/  → cCRE load, cell-line signals, distance matching, elem_focus
pipeline/evidence/        → SCREEN, ABC, HiChIP, merge → gene_links
pipeline/CHIP/            → unified ChIP, chip_hits on elements / SV
pipeline/tad_annotation/  → TAD domains on genes / elements (separate biosample registry)
```

**Evidence merge attachment:** `pipeline/evidence/evidence_merger.py` builds the per–cCRE `gene_links` column and attaches it in `main.py`.

**Do not** assume `pipeline/__init__.py` re-exports the main API; import from `pipeline.main` and `pipeline.config` explicitly (see `PIPELINE_STRUCTURE_AND_USAGE.md`).

---

## 4. Before you start a concrete task

1. **Confirm scope** — regulatory element pipeline vs SNV/SV/CNV/methylation/RPPA/ATAC submodule; they share `config` patterns but different entry scripts.
2. **Check `DATA_STRUCTURE_ATLAS.md`** for the exact nested shape you are modifying so serializers and smoke expectations stay aligned.
3. **Check `THRESHOLDS_AND_SIGNAL_STRENGTH.md`** if you change behavior of “weak/strong”, filters, or tiers.
4. **If renaming biosamples or adding panels** — update **`biosample_names.py`** maps and **`BiosampleConfig`** comments; keep **disk paths** and **upstream ENCODE strings** distinct (see `SAMPLE_ID_MATCHING_GUIDE.md`).
5. **Run or extend smoke tests** — `tests/run_smoke_tests.py` is the primary automated sanity check for this repo layout (run with **`.venv/bin/python3 -m tests.run_smoke_tests`** from the repo root on WSL).

---

## 5. Running commands reliably under WSL from Windows-hosted tools

**Cursor:** the always-on project rule **`.cursor/rules/wsl-shell-bridge.mdc`** restates the agent-specific `working_directory` + `wsl ... bash -lc` pattern below.

This repo lives on a **WSL filesystem** (`/home/stavz/...`). If the command runner is a **Windows PowerShell** shell, bash idioms like `&&`, `set -euo pipefail`, `$(...)`, and paths like `/home/...` can fail or silently do the wrong thing.

### Recommended pattern (PowerShell → WSL bash)

1. **Move the shell to a native Windows directory first** (avoids empty stdout from some runners when cwd is a WSL UNC path):

```powershell
Set-Location $env:USERPROFILE
```

2. **Run one-shot work in Linux** with `cd` inside the bash script:

```powershell
wsl.exe -d Ubuntu -- bash -lc 'cd /home/stavz/masters/gdc/APM && .venv/bin/python3 scripts/annotations/build_annotation_sample_manifests.py'
```

Notes:
- Use **single quotes** around the entire `bash -lc '…'` script so PowerShell does not expand `$VAR` or `$(...)`.
- Inside the single-quoted `-lc` argument, use normal bash (`&&`, `|`, `$(...)` inside bash if you wrap inner parts carefully, or use a here-doc—see below).
- **Python on WSL:** bare `python3` is often **system** Python and may **not** have project deps (e.g. `pandas`). From `…/APM` use **`.venv/bin/python3`** (or `source .venv/bin/activate` first) for `tests/run_smoke_tests.py`, `scripts/*.py`, and `python -m pipeline.main`.
- **Why step 1 matters:** if the tool’s persisted cwd is **`\\wsl.localhost\...` or `\\wsl$\...`**, captured **stdout/stderr may look empty** even when the command exits `0`. Setting cwd to e.g. `C:\Users\<you>` before `wsl.exe` fixes that for Cursor’s command runner. **Alternative:** open the workspace with **Remote - WSL** so paths are `/home/...` and the default shell is bash—then you can run commands directly with no bridge.
- **Heavy quoting:** for multi-line Python or nested quotes, avoid inline `python -c "..."` through PowerShell. Prefer `wsl … bash -lc 'cd …/APM && .venv/bin/python3 some_script.py'` or a bash **here-doc** inside the single-quoted `-lc` string so bash parses the script, not PowerShell.

To confirm the bridge is working, create a probe file and verify it exists:

```powershell
Set-Location $env:USERPROFILE
wsl.exe -d Ubuntu -- bash -lc 'cd /home/stavz/masters/gdc/APM; mkdir -p analysis/output; date -u +%Y%m%d_%H%M%S > analysis/output/_bridge_test_agent.txt'
```

---
## 6. Related work outside this tree

- **TAD calling / Hi-C heavy processing** may live in sibling repo **`TADs/`** (e.g. `TADs/config.py` for PAM50 cell-line metadata). The APM pipeline **consumes processed TADs** paths via `PATHS.tads_processed` and `tad_annotation/`; keep cross-repo facts consistent when editing subtype or TAD biosample names.

---

## 7. Quick checklist (copy for a PR / task comment)

- [ ] Read Tier A docs above  
- [ ] Opened relevant Atlas section + `config.py`  
- [ ] Grep’d for callers before renaming keys (`gene_links`, `screen_exp`, etc.)  
- [ ] Ran or noted `tests/run_smoke_tests.py` for touched behavior  
- [ ] Updated **Tier A/B** markdown only if behavior or contracts actually changed  

---

*Maintainers: when you add a new top-level output or module, link it from `PIPELINE_STRUCTURE_AND_USAGE.md` and (if user-facing structure) add a section to `DATA_STRUCTURE_ATLAS.md`, then extend this file’s Tier B/C tables if a new deep-dive doc exists.*
