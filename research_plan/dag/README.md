# DAG sources

Graphviz DOT files. DOT was chosen over mermaid for these because:

- Layout engines (`dot`, `neato`, `fdp`, `twopi`) handle dense graphs far better.
- Text-based and version-controlled.
- Render to SVG / PNG for papers and slides without a live browser.
- Online viewers work with raw DOT: [https://dreampuf.github.io/GraphvizOnline](https://dreampuf.github.io/GraphvizOnline), [https://edotor.net](https://edotor.net).

## Files


| File                            | What it is                                                                                               | Recommended layout                  |
| ------------------------------- | -------------------------------------------------------------------------------------------------------- | ----------------------------------- |
| `01_biological_circuit.dot`     | Closed immune-visibility circuit: antigen → APM → T-cell → IFN-γ loopback → NK axis + counter-regulators | `dot` top-to-bottom                 |
| `02_pipeline_architecture.dot`  | Layers L0 (pipeline outputs) → L1 (builders) → L2 (scores) → L3 (engine) → L4 (papers)                   | `dot` left-to-right                 |
| `03_within_panel_relations.dot` | Intra-panel gene–gene wiring: bottleneck, paralogs, antagonist ratios, feedback, coactivators            | `neato` or `fdp` for force-directed |
| `04_noncoding_space.dot`        | Methylation × miRNA × lncRNA integrated hypothesis space                                                 | `dot` top-to-bottom                 |
| `05_hypothesis_dependency.dot`  | H1–H58 as nodes; edges = shared feature builders; clusters = families A–L                                | `fdp` or `sfdp`                     |


## Render commands

Install Graphviz (`apt install graphviz` or `brew install graphviz`), then:

```bash
# SVG (scalable, recommended for papers)
dot -Tsvg 01_biological_circuit.dot -o 01_biological_circuit.svg

# High-DPI PNG (for slides)
dot -Tpng -Gdpi=200 02_pipeline_architecture.dot -o 02_pipeline_architecture.png

# Force-directed for dense relations
fdp -Tsvg 03_within_panel_relations.dot -o 03_within_panel_relations.svg
fdp -Tsvg 05_hypothesis_dependency.dot -o 05_hypothesis_dependency.svg

# One-shot all
for f in *.dot; do dot -Tsvg "$f" -o "${f%.dot}.svg"; done
```

## Editing guidelines

- Keep cluster labels short; detail lives in `02_hypothesis_catalog.md`.
- Use `shape=box` for data / tables, `shape=ellipse` for biological entities, `shape=diamond` for composite scores.
- Edge style: solid = direct mechanism, dashed = statistical / inferred, dotted = negative / inhibitory.
- Color palette kept small: `lightblue` = data/tables, `lightyellow` = features/scores, `lightgreen` = biological, `mistyrose` = counter-regulators, `lightgrey` = pipeline infrastructure.

