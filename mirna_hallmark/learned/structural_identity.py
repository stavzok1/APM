"""Structural identity — the LOSS lens (refined precursor identity; Design §Decision I loss variant).

"Who is DESIGNED to specifically repress g", **expression-free**, so it stays high for a *silenced specialist*
that the coupling-based lenses (M / budget / Shapley identity) all miss (M is gated by current coupling → 0).

Candidate regulators = POOLED-HE edges = miRTarBase-HE ∪ TarBase-v9 low-throughput functional
(reporter/western/proteomics, non-weak) ≈ 5,940 (HE 5,300 + ~640 TarBase) — the multi-source analog of
miRTarBase-HE; excludes prediction-only, binding-only CLIP, and high-throughput RNA-seq/microarray. Restricted
to arms with a biochemical potential = scanMiR-modelled, with a TargetScan-context++ FALLBACK for arms scanMiR
lacks (e.g. miR-137). Expression-free (no arm floor — expression re-enters at lost_specialist).

    potential(f,g)   = repression POTENTIAL — biochemical (scanMiR K_D / TargetScan context++) FUSED with
                       evidence (ledger existence·directness, ×(1+log1p)) and direction-filtered (drop TarBase-v9
                       validated-activating edges). "How strong (and validated) a repressor f WOULD be if expressed."
    specificity(f,g) = ev-mass fraction  Σ_members ev(·,g) / Σ_members Σ_g' ev(·,g')  — the family's POOLED
                       evidence-mass fraction aimed at g (specialist = concentrated; hub = spread). Evidence-only;
                       there is NO biochemical specialist to floor to (every arm hits ~700 targets ~uniformly).
    confidence(f)    = min(1, ev_depth(f)/CONF_K)   (ev_depth = # genes with ANY family evidence) — trust in the
                       specificity claim; enters as a MULTIPLIER (no evidence ⇒ no ownership), not a blend.
    structural(f,g)  = potential(f,g) · specificity(f,g) · confidence(f)  → gene-local normalised = **structural_share**.

Family = the identified unit (§F). Overlaps M only on the PRIOR; diverges exactly on silenced arms.
`lost_specialist(gene)` adds the loss layer (methylation-gain + baseline-active + tumour-silenced) — separate.

CLI: `python -m mirna_hallmark.learned.structural_identity PTEN ESR1 [--source scanmir|targetscan|fused]`
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from mirna_hallmark import data_loaders as D
from mirna_hallmark.learned import families as FAM

_CACHE: dict = {}
CONF_K = 10   # pooled-family study depth (# genes) at which specificity is fully trusted; below → structural shrinks ×conf
def _functional_regulators(gene: str) -> set:
    """Arms on `gene` in the POOLED-HE set (`ledger.pooled_he_edges`: miRTarBase-HE ∪ TarBase-v9 low-throughput
    functional, ~5,940). SAME set the learned model now uses (`assemble_gene(he_only=True)`), so the lenses share
    one inclusion policy."""
    if "func_edges" not in _CACHE:
        from mirna_hallmark.learned.evidence import ledger as LG
        _CACHE["func_edges"] = LG.pooled_he_edges().rename(columns={"miRNA": "arm"})[["arm", "gene"]]
    fe = _CACHE["func_edges"]
    return set(fe.loc[fe["gene"] == gene, "arm"])


def _family_members() -> dict:
    """seed-family → ALL its member arms (over the full scanMiR arm universe), cached. Used to pool evidence
    over the WHOLE family (co-seed arms share the targetome) rather than only a gene's candidate members."""
    if "fam_members" not in _CACHE:
        strength, _, _ = _potential_matrix()
        fm = FAM.family_of(pd.Index(list(strength.index)))
        _CACHE["fam_members"] = {f: list(v) for f, v in fm.groupby(fm).groups.items()}
    return _CACHE["fam_members"]


def _potential_matrix(source: str = "fused"):
    """(strength arm×gene, biochem_promiscuity arm, evidence arm×gene). strength = biochemical (scanMiR |K_D| /
    TargetScan) FUSED with evidence (ledger, direction-filtered). biochem_promisc = genome-wide strong-site
    (8mer/7mer) target count. evidence = ledger weight matrix (for the evidence-targetome specificity + depth)."""
    key = f"pot_{source}"
    if key in _CACHE:
        return _CACHE[key]
    from mirna_hallmark.learned import kd as KD
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.evidence import ledger as LG
    aff = KD.affinity()
    strong = (aff["n_8mer"].fillna(0) + aff["n_7mer"].fillna(0)) >= 1     # canonical strong site
    promisc = aff[strong].groupby("arm").size()                          # # strong-site target genes (biochemical)
    sc = aff.assign(mag=(-aff["repression"]).clip(lower=0)).pivot_table(
        index="arm", columns="gene", values="mag", aggfunc="max", fill_value=0.0)
    tsdf = LD._targetscan_context()                                      # TargetScan context++ magnitude (arm×gene)
    ts = (tsdf.pivot_table(index="arm", columns="gene", values="ts_mag", aggfunc="max", fill_value=0.0)
          if tsdf is not None else pd.DataFrame())
    genes = sc.columns.union(ts.columns)
    sc = sc.reindex(columns=genes, fill_value=0.0); ts = ts.reindex(columns=genes, fill_value=0.0)
    if source == "targetscan":
        base = ts
    elif source == "scanmir":
        base = sc
    else:  # fused = scanMiR base, with a TargetScan FALLBACK for arms scanMiR doesn't model (e.g. miR-137),
        missing = ts.index.difference(sc.index)                          # rescaled to scanMiR's magnitude scale so
        scpos, tspos = sc.values[sc.values > 0], ts.values[ts.values > 0]  # a fallback arm is comparable, not dominant
        scale = float(np.median(scpos) / np.median(tspos)) if len(tspos) and len(scpos) else 1.0
        base = pd.concat([sc, ts.loc[missing] * scale]) if len(missing) else sc
    ev = LG.edge_weights().pivot_table(index="arm", columns="gene", values="ledger_weight", aggfunc="max",
                                       fill_value=0.0).reindex(index=base.index, columns=base.columns, fill_value=0.0)
    strength = base if source in ("scanmir", "targetscan") else base * (1.0 + np.log1p(ev))  # fused × evidence
    try:                                                                 # direction filter: drop validated-activating
        from mirna_hallmark.learned.evidence.direction import tarbase_direction
        for r in tarbase_direction().query("direction=='activating'").itertuples():
            if r.arm in strength.index and r.gene in strength.columns:
                strength.at[r.arm, r.gene] = 0.0
    except Exception:
        pass
    _CACHE[key] = (strength, promisc, ev)
    return _CACHE[key]


def structural_identity(gene: str, *, source: str = "fused") -> pd.DataFrame:
    """Baseline-free structural-identity shares for gene g's regulator seed-families. Candidate arms = the
    pooled-ledger FUNCTIONAL edges (multi-source analog of miRTarBase-HE; no prediction/binding-only),
    ∩ scanMiR-modelled arms; expression-free (arm floor re-enters at `lost_specialist`). Evidence is POOLED
    across each family's members (co-seed arms are one curated body). Columns: potential (family's strongest
    member's fused strength on g), ev_depth/confidence (family-pooled study depth), specificity (pooled
    evidence-mass fraction on g), structural = potential·specificity·confidence, structural_share (gene-local)."""
    strength, promisc, ev = _potential_matrix(source)
    if gene not in strength.columns:
        print(f"{gene}: not in potential matrix"); return pd.DataFrame()
    arms = sorted(_functional_regulators(gene) & set(strength.index))  # pooled-ledger FUNCTIONAL edges ∩ scanMiR
    if not arms:                                                        # (expression-free; arm floor re-enters at lost_specialist)
        print(f"{gene}: no functional regulator arms with scanMiR potential"); return pd.DataFrame()
    ev_mass = ev.sum(axis=1)                                            # per-arm total curated-evidence mass (for pooling)
    cand_fams = set(FAM.family_of(pd.Index(arms)))                      # seed families with a functional regulator of g
    full = _family_members()                                            # family → ALL its seed-member arms (cached)
    rows = {}
    for f in cand_fams:
        mem = [m for m in full.get(f, []) if m in strength.index]       # FULL seed family — co-seed arms share the
        if not mem:                                                     # targetome, so pool the WHOLE family's evidence
            continue                                                    # (else a thinly-studied member fakes concentration)
        top = strength.loc[mem, gene].idxmax()                         # POTENTIAL = family's biochem capacity = strongest member
        pot = float(strength.loc[top, gene])
        # SPECIFICITY = evidence-mass fraction on g, POOLED across the full family. CONFIDENCE (pooled study depth)
        # enters as a MULTIPLIER (no evidence ⇒ no ownership) — there is no biochemical specialist to floor to.
        evm = ev.loc[mem]                                              # members × gene
        ev_fam_g, ev_fam_mass = float(evm[gene].sum()), float(ev_mass.loc[mem].sum())
        ev_depth_fam = int((evm > 0).any(axis=0).sum())                # # genes with ANY family evidence (pooled breadth)
        spec = ev_fam_g / ev_fam_mass if ev_fam_mass > 0 else 0.0      # fraction of family's evidence-mass aimed at g
        conf = min(1.0, ev_depth_fam / CONF_K)
        rows[f] = {"potential": round(pot, 2), "ev_depth": ev_depth_fam, "confidence": round(conf, 2),
                   "specificity": round(spec, 4), "structural": pot * spec * conf}
    df = pd.DataFrame(rows).T
    tot = df["structural"].sum()
    df["structural_share"] = (df["structural"] / tot).round(3) if tot > 0 else 0.0
    df = df.sort_values("structural_share", ascending=False)
    with pd.option_context("display.width", 150):
        print(f"\n=== {gene} — STRUCTURAL identity (potential × pooled-evidence specificity × confidence, source={source}) ===")
        print(df.round(3).to_string())
    return df


def compare_lenses(gene: str, *, source: str = "fused") -> pd.DataFrame:
    """Structural (designed owner) vs Shapley identity (current coupling ownership) vs budget (current force).
    A silenced specialist = HIGH structural, LOW identity, LOW budget."""
    from mirna_hallmark.learned.attribution import identity_vs_magnitude
    import contextlib, io
    s = structural_identity(gene, source=source)
    with contextlib.redirect_stdout(io.StringIO()):
        im = identity_vs_magnitude(gene)
    if s.empty or im.empty:
        return pd.DataFrame()
    m = s[["structural_share"]].join(im[["indiv_rho", "identity_share", "magnitude_share"]], how="outer").fillna(0.0)
    m = m.sort_values("structural_share", ascending=False)
    m["struct_minus_current"] = (m["structural_share"] - m[["identity_share", "magnitude_share"]].max(axis=1)).round(3)
    print(f"\n=== {gene} — STRUCTURAL (designed) vs IDENTITY (couples) vs BUDGET (force) ===")
    print(m.round(3).to_string())
    cand = m[(m["structural_share"] > 0.15) & (m["identity_share"] < 0.05) & (m["magnitude_share"] < 0.05)]
    if len(cand):
        print(f"  ⚑ silenced-specialist candidates (high structural, ~no current coupling/force): {list(cand.index)}")
    return m


def lost_specialist(gene: str, *, source: str = "fused", floor: float = np.log2(11), top: int = 20) -> pd.DataFrame:
    """LOSS layer (Design §H loss): which STRUCTURAL specialists for g were **baseline-active** (expressed in
    GTEx and/or NAT) but are **tumour-silenced** (abundance dropped below the functional floor)? These are the
    'lost repressor' candidates. The structural score nominates (who's designed to own g); the baseline
    proves it was present; the tumour drop shows it's gone. **Methylation gate** — tumour-vs-normal promoter
    Δβ over CpG probes that DIRECTLY overlap the miRNA promoter window (`methylation.locus_methylation`; no
    cCRE hop) — is the mechanism confirmation (`meth_confirmed`), priced only for actual loss calls. Requires
    the beta cache (`python -m mirna_hallmark.learned.methylation --build`); shown as ⏳ when absent."""
    from mirna_hallmark.learned import state as STA
    from mirna_hallmark.learned import states as ST
    from mirna_hallmark.learned import methylation as METH
    s = structural_identity(gene, source=source)
    if s.empty:
        return s
    fam = FAM.family_of(pd.Index(sorted(_functional_regulators(gene))))  # same candidate set as structural_identity
    fam_arms = {f: [a for a in members] for f, members in fam.groupby(fam).groups.items()}
    gtex = STA._gtex_mirna(); Xn, _ = ST.state_matrices("11"); Xt, _ = ST.state_matrices("01")
    have_meth = METH._BETA_CACHE.exists()                                # methylation mechanism leg (else ⏳)

    def _mean(mat, arms):
        arms = [a for a in arms if a in mat.index]
        return float(mat.loc[arms].mean(axis=1).max()) if arms else np.nan   # strongest member's mean level

    rows = []
    for f in s.index:
        arms = fam_arms.get(f, [f])
        g_ex, n_ex, t_ex = _mean(gtex, arms), _mean(Xn, arms), _mean(Xt, arms)
        baseline_active = ((g_ex == g_ex and g_ex > floor) or (n_ex == n_ex and n_ex > floor))
        tumour_silenced = (t_ex == t_ex and t_ex < floor)
        lost = bool(s.loc[f, "structural_share"] > 0.02 and baseline_active and tumour_silenced)
        dbeta, meth_conf = np.nan, False
        if have_meth and lost:                                          # only price methylation for actual loss calls
            mr = METH.locus_methylation(arms)
            if mr.get("n_probes", 0) > 0:
                dbeta, meth_conf = mr["delta_beta"], mr["hypermethylated"]
        rows.append({"family": f, "structural_share": s.loc[f, "structural_share"],
                     "confidence": s.loc[f, "confidence"], "gtex_abund": round(g_ex, 1) if g_ex == g_ex else np.nan,
                     "nat_abund": round(n_ex, 1) if n_ex == n_ex else np.nan,
                     "tum_abund": round(t_ex, 1) if t_ex == t_ex else np.nan,
                     "baseline_active": baseline_active, "tumour_silenced": tumour_silenced,
                     "delta_beta": dbeta, "meth_confirmed": meth_conf, "LOST": lost})
    df = pd.DataFrame(rows).sort_values(["LOST", "meth_confirmed", "structural_share"], ascending=False)
    with pd.option_context("display.width", 170):
        print(f"\n=== {gene} — LOST-SPECIALIST candidates (structural + baseline-active + tumour-silenced"
              f"{' + promoter Δβ' if have_meth else ''}) ===")
        print(df.head(top).to_string(index=False))
    lost = df[df["LOST"]]
    tag = ("promoter Δβ shown; meth_confirmed = hyper-methylated (mechanism)" if have_meth
           else "⏳ build methylation cache (`-m ...learned.methylation --build`) for the mechanism leg")
    print(f"  ⚑ LOST candidates: {list(lost['family'])}  ({tag})")
    return df


if __name__ == "__main__":
    src = next((a.split("=")[1] for a in sys.argv if a.startswith("--source=")), "fused")
    genes = [a for a in sys.argv[1:] if not a.startswith("-")] or ["PTEN", "ESR1"]
    if "--lost" in sys.argv:
        for g in genes:
            lost_specialist(g, source=src)
    else:
        for g in genes:
            compare_lenses(g, source=src)
