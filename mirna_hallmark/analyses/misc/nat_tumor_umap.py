"""UMAP of TCGA-BRCA miRNA-seq: NAT / tumor / optional GTEx healthy breast.

Comparison modes:

  matched   — 104 NAT + matched primary tumors (n=103)
  unmatched — 104 NAT + 104 primary tumors from participants without NAT

With ``--include-gtex`` (default): adds 346 GTEx v10 healthy breast donors.
Axis labels show post-hoc linear variance captured per UMAP dimension; subtitle
includes PCA PC1/PC2 reference (UMAP is nonlinear — see manifest).

Run:
  .venv/bin/python3 -m mirna_hallmark.nat_tumor_umap --comparison unmatched
  .venv/bin/python3 -m mirna_hallmark.nat_tumor_umap --comparison matched --no-gtex
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Literal, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.analyses.cross_state.cross_state_landscape import _gtex_all_arm_matrix
from mirna_hallmark.family_normal_reference import (
    _gtex_available,
    _load_full_mirna,
    _participant,
    _split_types,
)

OUT_DIR = C.TISSUE_REFERENCE_DIR / "nat_tumor_umap"
Comparison = Literal["matched", "unmatched"]

TISSUE_COLORS = {
    "GTEx_healthy": "#59A14F",
    "NAT": "#EDC948",
    "tumor": "#E15759",
    "tumor_unmatched": "#E15759",
}

GTEX_BARCODE_PREFIX = "GTEx:"


def _tumor_by_participant(tumor: pd.DataFrame) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for c in tumor.columns:
        out.setdefault(_participant(c), []).append(c)
    return out


def _gtex_sample_barcodes() -> List[str]:
    gtex = _load_gtex_matrix()
    return [f"{GTEX_BARCODE_PREFIX}{c}" for c in gtex.columns.astype(str)]


def _load_gtex_matrix() -> pd.DataFrame:
    if not _gtex_available():
        raise FileNotFoundError(
            "GTEx v10 miRNA files missing under data/GTEx/ — needed for --include-gtex"
        )
    return _gtex_all_arm_matrix()


def _append_gtex_meta(meta: pd.DataFrame, sample_order: List[str], comparison: Comparison) -> None:
    for col in _load_gtex_matrix().columns.astype(str):
        bid = f"{GTEX_BARCODE_PREFIX}{col}"
        sample_order.append(bid)
        meta.loc[len(meta)] = {
            "barcode": bid,
            "participant": col,
            "tissue": "GTEx_healthy",
            "comparison": comparison,
            "is_matched_pair": False,
        }


def select_cohort(
    nat: pd.DataFrame,
    tumor: pd.DataFrame,
    comparison: Comparison,
    *,
    n_unmatched_tumors: int = 104,
    include_gtex: bool = True,
) -> Tuple[List[str], pd.DataFrame]:
    nat_by_p = {_participant(c): c for c in nat.columns}
    tumor_by_p = _tumor_by_participant(tumor)
    nat_participants = set(nat_by_p)

    rows: List[dict] = []
    sample_order: List[str] = []

    for p in sorted(nat_by_p):
        nc = nat_by_p[p]
        sample_order.append(nc)
        rows.append(
            {
                "barcode": nc,
                "participant": p,
                "tissue": "NAT",
                "comparison": comparison,
                "is_matched_pair": False,
            }
        )

    if comparison == "matched":
        for p in sorted(nat_by_p):
            if p not in tumor_by_p:
                continue
            tc = sorted(tumor_by_p[p])[0]
            sample_order.append(tc)
            rows.append(
                {
                    "barcode": tc,
                    "participant": p,
                    "tissue": "tumor",
                    "comparison": comparison,
                    "is_matched_pair": True,
                }
            )
    elif comparison == "unmatched":
        candidates: List[str] = []
        for p in sorted(tumor_by_p):
            if p in nat_participants:
                continue
            candidates.append(sorted(tumor_by_p[p])[0])
        if len(candidates) < n_unmatched_tumors:
            raise ValueError(
                f"Only {len(candidates)} non-NAT-participant tumors available; "
                f"requested {n_unmatched_tumors}"
            )
        for tc in candidates[:n_unmatched_tumors]:
            sample_order.append(tc)
            rows.append(
                {
                    "barcode": tc,
                    "participant": _participant(tc),
                    "tissue": "tumor_unmatched",
                    "comparison": comparison,
                    "is_matched_pair": False,
                }
            )
    else:
        raise ValueError(f"Unknown comparison {comparison!r}")

    meta = pd.DataFrame(rows)
    if include_gtex:
        _append_gtex_meta(meta, sample_order, comparison)
    return sample_order, meta


def _combined_expression(
    tcga_full: pd.DataFrame,
    sample_ids: Sequence[str],
    *,
    include_gtex: bool,
) -> pd.DataFrame:
    """arm x samples log2(PM+1); GTEx columns prefixed ``GTEx:``."""
    tcga_ids = [s for s in sample_ids if not str(s).startswith(GTEX_BARCODE_PREFIX)]
    parts = [tcga_full.loc[:, tcga_ids]]
    if include_gtex:
        gtex = _load_gtex_matrix()
        gtex_ids = [s for s in sample_ids if str(s).startswith(GTEX_BARCODE_PREFIX)]
        donors = [s.removeprefix(GTEX_BARCODE_PREFIX) for s in gtex_ids]
        g = gtex.loc[:, donors].copy()
        g.columns = [f"{GTEX_BARCODE_PREFIX}{c}" for c in g.columns]
        parts.append(g)
    combined = pd.concat(parts, axis=1)
    return combined.groupby(level=0).mean()


def _feature_matrix(
    combined: pd.DataFrame,
    samples: Sequence[str],
    *,
    hsa_only: bool = True,
    top_var: int = 500,
    min_var: float = 1e-6,
) -> Tuple[pd.DataFrame, pd.Series]:
    sub = combined.loc[:, list(samples)]
    if hsa_only:
        sub = sub.loc[sub.index.astype(str).str.startswith("hsa-")]
    mat = sub.T.copy()
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    var = mat.var(axis=0)
    keep = var[var > min_var].sort_values(ascending=False)
    if top_var and len(keep) > top_var:
        keep = keep.head(top_var)
    return mat.loc[:, keep.index], keep


def _run_umap(
    X: np.ndarray,
    *,
    n_neighbors: int = 15,
    min_dist: float = 0.25,
    metric: str = "euclidean",
    random_state: int = 0,
) -> np.ndarray:
    import umap

    n = X.shape[0]
    nn = min(n_neighbors, max(2, n - 1))
    reducer = umap.UMAP(
        n_neighbors=nn,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
        n_jobs=1,
    )
    return reducer.fit_transform(X)


def _variance_summary(Xz: np.ndarray, emb: np.ndarray) -> Dict[str, float]:
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LinearRegression

    pca = PCA(n_components=min(5, Xz.shape[0], Xz.shape[1]))
    pca.fit(Xz)
    evr = pca.explained_variance_ratio_

    def _r2(y: np.ndarray, x: np.ndarray) -> float:
        return float(LinearRegression().fit(x, y).score(x, y))

    u1 = _r2(Xz, emb[:, [0]]) * 100
    u2 = _r2(Xz, emb[:, [1]]) * 100
    u12 = _r2(Xz, emb[:, :2]) * 100
    spread = emb.var(axis=0)
    spread1 = float(spread[0] / spread.sum() * 100) if spread.sum() > 0 else float("nan")

    return {
        "pc1_variance_pct": float(evr[0] * 100),
        "pc2_variance_pct": float(evr[1] * 100) if len(evr) > 1 else float("nan"),
        "pc1_pc2_cumulative_pct": float(evr[:2].sum() * 100) if len(evr) > 1 else float("nan"),
        "umap1_linear_r2_pct": u1,
        "umap2_linear_r2_pct": u2,
        "umap1_umap2_linear_r2_pct": u12,
        "umap1_embedding_spread_pct": spread1,
    }


def plot_umap(
    coords: pd.DataFrame,
    meta: pd.DataFrame,
    out_path: Path,
    *,
    comparison: Comparison,
    var_stats: Dict[str, float],
    include_gtex: bool,
    title: Optional[str] = None,
) -> None:
    import matplotlib.pyplot as plt

    df = coords.merge(meta, on="barcode", how="left")
    fig, ax = plt.subplots(figsize=(9.5, 7.5))

    plot_groups: List[Tuple[str, str, str, int, float]] = []
    if include_gtex:
        plot_groups.append(("GTEx_healthy", "GTEx healthy", TISSUE_COLORS["GTEx_healthy"], 18, 0.45))
    plot_groups.append(("NAT", "NAT", TISSUE_COLORS["NAT"], 42, 0.82))
    plot_groups.append(
        (
            "tumor_unmatched" if comparison == "unmatched" else "tumor",
            "tumor (no NAT)" if comparison == "unmatched" else "tumor (matched)",
            TISSUE_COLORS["tumor"],
            42,
            0.82,
        )
    )

    for tissue_key, label, color, size, alpha in plot_groups:
        sub = df.loc[df["tissue"].eq(tissue_key)]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP1"],
            sub["UMAP2"],
            c=color,
            label=f"{label} (n={len(sub)})",
            s=size,
            alpha=alpha,
            edgecolors="white",
            linewidths=0.35,
            zorder=2 if tissue_key == "GTEx_healthy" else 4,
        )

    if comparison == "matched":
        for p, grp in df.groupby("participant"):
            if grp["tissue"].nunique() < 2:
                continue
            nat = grp.loc[grp["tissue"].eq("NAT"), ["UMAP1", "UMAP2"]]
            tum = grp.loc[grp["tissue"].eq("tumor"), ["UMAP1", "UMAP2"]]
            if nat.empty or tum.empty:
                continue
            ax.plot(
                [nat["UMAP1"].iloc[0], tum["UMAP1"].iloc[0]],
                [nat["UMAP2"].iloc[0], tum["UMAP2"].iloc[0]],
                color="#cccccc",
                linewidth=0.6,
                alpha=0.45,
                zorder=1,
            )

    if title is None:
        title = "miRNA-seq UMAP: NAT vs primary tumor"
        if include_gtex:
            title += " + GTEx healthy"

    n_nat = int((meta["tissue"] == "NAT").sum())
    n_tum = int(meta["tissue"].isin(["tumor", "tumor_unmatched"]).sum())
    n_gtex = int((meta["tissue"] == "GTEx_healthy").sum())
    parts = [f"{n_nat} NAT", f"{n_tum} tumor"]
    if include_gtex:
        parts.append(f"{n_gtex} GTEx healthy")
    subtitle = " + ".join(parts)
    if comparison == "unmatched":
        subtitle += " (tumors: participants without NAT)"

    var_line = (
        f"Axis variance (linear R² vs {int(var_stats.get('n_features', 500))} miRNA features): "
        f"UMAP1={var_stats['umap1_linear_r2_pct']:.1f}%, "
        f"UMAP2={var_stats['umap2_linear_r2_pct']:.1f}% "
        f"(joint={var_stats['umap1_umap2_linear_r2_pct']:.1f}%)  ·  "
        f"PCA ref: PC1={var_stats['pc1_variance_pct']:.1f}%, "
        f"PC2={var_stats['pc2_variance_pct']:.1f}%"
    )

    ax.set_title(f"{title}\n{subtitle}\n{var_line}", fontsize=10)
    ax.set_xlabel(f"UMAP1  ({var_stats['umap1_linear_r2_pct']:.1f}% profile variance)")
    ax.set_ylabel(f"UMAP2  ({var_stats['umap2_linear_r2_pct']:.1f}% profile variance)")
    ax.legend(frameon=True, loc="best", fontsize=9)
    ax.grid(alpha=0.2, linestyle="--")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def run(
    *,
    comparison: Comparison = "matched",
    out_dir: Path = OUT_DIR,
    include_gtex: bool = True,
    top_var: int = 500,
    n_unmatched_tumors: int = 104,
    n_neighbors: int = 15,
    min_dist: float = 0.25,
    random_state: int = 0,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)

    gtex_tag = "_gtex" if include_gtex else ""
    prefix = f"nat_{comparison}_tumor{gtex_tag}_umap"

    full = _load_full_mirna()
    split = _split_types(full)
    samples, meta = select_cohort(
        split["nat"],
        split["tumor"],
        comparison,
        n_unmatched_tumors=n_unmatched_tumors,
        include_gtex=include_gtex,
    )
    combined = _combined_expression(full, samples, include_gtex=include_gtex)
    feat, var = _feature_matrix(combined, samples, top_var=top_var)

    X = feat.to_numpy(dtype=float)
    Xz = (X - X.mean(axis=0)) / (X.std(axis=0, ddof=0) + 1e-8)
    nn = n_neighbors if not include_gtex else max(n_neighbors, 30)
    emb = _run_umap(Xz, n_neighbors=nn, min_dist=min_dist, random_state=random_state)
    var_stats = _variance_summary(Xz, emb)
    var_stats["n_features"] = float(feat.shape[1])

    coords = pd.DataFrame({"barcode": feat.index, "UMAP1": emb[:, 0], "UMAP2": emb[:, 1]})
    coords = coords.merge(meta, on="barcode", how="left")
    coords.to_csv(out_dir / f"{prefix}_coordinates.tsv", sep="\t", index=False)
    meta.to_csv(out_dir / f"{prefix}_samples.tsv", sep="\t", index=False)
    var.rename("variance").to_csv(out_dir / f"{prefix}_features.tsv", sep="\t", header=True)

    png = out_dir / "figures" / f"{prefix}.png"
    plot_umap(
        coords[["barcode", "UMAP1", "UMAP2"]],
        meta,
        png,
        comparison=comparison,
        var_stats=var_stats,
        include_gtex=include_gtex,
    )

    manifest = {
        "module": "mirna_hallmark.nat_tumor_umap",
        "comparison": comparison,
        "include_gtex": include_gtex,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_nat": int((meta["tissue"] == "NAT").sum()),
        "n_tumor": int(meta["tissue"].isin(["tumor", "tumor_unmatched"]).sum()),
        "n_gtex": int((meta["tissue"] == "GTEx_healthy").sum()),
        "n_features": int(feat.shape[1]),
        "top_var": top_var,
        "umap": {"n_neighbors": nn, "min_dist": min_dist, "random_state": random_state},
        "variance_summary": var_stats,
        "variance_note": (
            "UMAP axes are nonlinear; axis % = post-hoc linear R² of z-scored miRNA "
            "features regressed on each UMAP coordinate. PCA % is the linear reference."
        ),
        "cross_platform_note": "GTEx TPM/URS vs TCGA RPM/MIMAT — joint UMAP is shape exploration only",
        "figure": str(png.relative_to(C.OUTPUT_ROOT)),
    }
    (out_dir / f"{prefix}_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(
        f"[nat_tumor_umap] {comparison}{'+GTEx' if include_gtex else ''}: "
        f"samples={len(samples)} ({manifest['n_nat']} NAT, {manifest['n_tumor']} tumor, "
        f"{manifest['n_gtex']} GTEx); features={manifest['n_features']}"
    )
    print(
        f"[nat_tumor_umap] variance: UMAP1={var_stats['umap1_linear_r2_pct']:.1f}% "
        f"UMAP2={var_stats['umap2_linear_r2_pct']:.1f}% "
        f"(PCA PC1+PC2={var_stats['pc1_pc2_cumulative_pct']:.1f}%)"
    )
    print(f"[nat_tumor_umap] figure -> {png}")
    return coords


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison", choices=["matched", "unmatched"], default="unmatched")
    ap.add_argument("--include-gtex", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--top-var", type=int, default=500)
    ap.add_argument("--n-unmatched", type=int, default=104)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--min-dist", type=float, default=0.25)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    C.ensure_output_dirs()
    run(
        comparison=args.comparison,
        out_dir=args.out_dir,
        include_gtex=args.include_gtex,
        top_var=args.top_var,
        n_unmatched_tumors=args.n_unmatched,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        random_state=args.seed,
    )


if __name__ == "__main__":
    main()
