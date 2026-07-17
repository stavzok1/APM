# miRNA-Hallmark WIP talk — build

Self-contained deck for the work-in-progress lab seminar. Narrative + figure spec live in the
plan; this folder is the runnable artifact.

## Files
- `talk.qmd` — the Quarto deck (reveal.js + pptx), LaTeX math, references `figures/`.
- `make_figures.py` — regenerates every data-driven figure from `mirna_hallmark/output/*` and
  copies the pre-rendered pipeline figures (grant heatmap, CNV boxplot, NMF heatmap,
  edge-transition, Buffa) into `figures/`. Pure pandas+matplotlib, one figure per function,
  failure-isolated.
- `figures/` — generated PNGs (F*.png).
- `talk.html` / `talk.pptx` — rendered outputs.

## Rebuild
```bash
# 1. figures (after any pipeline refresh)
.venv/bin/python3 mirna_hallmark/presentation/make_figures.py

# 2. render — needs Quarto on PATH
export PATH="/cs/usr/stavzok/opt/quarto/bin:$PATH"
cd mirna_hallmark/presentation
quarto render talk.qmd              # both formats
quarto render talk.qmd --to pptx    # PowerPoint object only
quarto render talk.qmd --to revealjs
quarto preview talk.qmd             # live preview while editing
```

## Notes
- Quarto 1.9.38 installed at `/cs/usr/stavzok/opt/quarto/`.
- `.pptx` is produced by Quarto→Pandoc directly (LaTeX `$…$` → native PowerPoint equations,
  figures embedded). No PowerPoint/LibreOffice/Chromium needed to build it.
- reveal.js extras (fragments, two-column layouts) render in HTML; in PowerPoint they flatten to
  simple one-image-per-slide layouts — keep that in mind if `.pptx` is the primary object.
- To brand the PowerPoint: save a themed blank deck as `template.pptx` and set
  `pptx: {reference-doc: template.pptx}` in the YAML.
- Provenance: figures are from the 2026-06-23/24 S1+MIMAT refresh; the dosage figure (F13) /
  slide 15 number follow the latest `mirna_locus_cnv` run.
