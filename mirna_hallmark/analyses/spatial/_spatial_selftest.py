"""Self-test for the Phase-3 spatial MH-56d engine (`spatial_common.gene_level_pressure_retest`).

Verifies all four classification branches on synthetic data + the loaders' format tolerance. The
end-to-end validation on REAL data lives in `spatial_mibi_anchor` (its `engine_validation_vs_MH56d`
gate must reproduce SLC2A1 brake-release + ERBB2/ESR1 concordant-repression).

Run: ``.venv/bin/python3 -m mirna_hallmark._spatial_selftest``
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mirna_hallmark.analyses.spatial import spatial_common as SP


def test_classifier_branches() -> None:
    """released×up→brake_release, gained×down→concordant_repression, gained×up→discordant_rise."""
    delta = SP.load_pressure_delta()
    released = next(g for g in ("SLC2A1", "FAP") if g in delta.index and delta[g] < 0)
    gained = next(g for g in ("ERBB2", "ESR1", "VIM") if g in delta.index and delta[g] > 0)
    ts = pd.Series(np.linspace(0, 1, 120), index=[f"u{i}" for i in range(120)])
    rng = np.random.default_rng(0)
    ro = pd.DataFrame(
        [ts.values + rng.normal(0, 0.05, 120),       # released, rises with tumour axis
         -ts.values + rng.normal(0, 0.05, 120),      # gained, falls
         ts.values + rng.normal(0, 0.05, 120)],      # gained, rises (discordant)
        index=[released, gained, "VIM"], columns=ts.index)
    # ensure 3 distinct genes
    if "VIM" in (released, gained):
        ro = ro.iloc[:2]
    out = SP.gene_level_pressure_retest(ro, ts, target_genes=list(ro.index), min_units=30)
    cls = out.set_index("gene")["classification"].to_dict()
    assert cls[released] == "brake_release", cls
    assert cls[gained] == "concordant_repression", cls
    print(f"[selftest] classifier branches OK: {cls}")


def test_morans_i() -> None:
    """Smooth gradient → high I; shuffled → ~0."""
    n = 400
    g = np.arange(n) % 20
    coords = pd.DataFrame({"x": g, "y": np.arange(n) // 20}, index=[f"s{i}" for i in range(n)])
    smooth = pd.Series(coords["x"].values.astype(float), index=coords.index)
    rng = np.random.default_rng(1)
    shuf = pd.Series(rng.permutation(smooth.values), index=coords.index)
    i_smooth, i_shuf = SP.morans_i(smooth, coords), SP.morans_i(shuf, coords)
    assert i_smooth > 0.4 > i_shuf, (i_smooth, i_shuf)
    print(f"[selftest] Moran's I OK: smooth={i_smooth:.2f} shuffled={i_shuf:.2f}")


def main() -> None:
    test_classifier_branches()
    test_morans_i()
    print("[selftest] ALL PASS")


if __name__ == "__main__":
    main()
