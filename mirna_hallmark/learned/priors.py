"""Per-edge priors π (inclusion), μ (magnitude), τ (confidence) — Design §Decision D/E, now BUILT.

Non-circular priors for the spike-and-slab (`spike_slab.py`): everything here comes from **curation +
sequence/biochemistry only, never TCGA (X,Y)**. Three objects, per candidate regulator of a gene:

    π_{g,m}  INCLUSION prob   — evidence-graded, monotone in the PMID-deduped ledger weight (§2). The
                               spike-and-slab's Bernoulli(π) inclusion prior; sets who is a-priori likely real.
    μ_{g,m}  MAGNITUDE        — biochemical affinity from scanMiR predicted repression (family = strongest
                               member); normalised within-gene to O(1). Widens the slab for high-affinity
                               edges (the Bayesian analog of the adaptive-lasso inverse penalty) — location
                               stays at 0 (no bias), only the allowed SCALE grows, so a real-but-weak arm is
                               kept small-but-nonzero instead of being zeroed. Design "prior sets ordering,
                               data sets magnitude" (§3) is preserved: μ never pulls the coefficient up.
    τ_{g,m}  CONFIDENCE       — in [0,1], increasing in evidence depth (HE deep → 1, orphan/thin → →0). A
                               multiplier on the slab scale: a thinly-evidenced edge gets a tighter slab.

`edge_priors(gene)` returns these for a gene's POOLED-HE candidates at arm or seed-family resolution.
The slab scale the sampler should use is `slab_scale = 1 + mu_gain·μ·τ` (see `slab_scale`).

CLI: `.venv/bin/python3 -m mirna_hallmark.learned.priors PTEN`
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy.special import expit

from mirna_hallmark.learned import data as LD
from mirna_hallmark.learned import families as FAM


def inclusion_prior(w: np.ndarray, *, p0: float = 0.3, kappa: float = 1.5, rarity_gain: float = 0.0,
                    rarity: Optional[np.ndarray] = None, spec_gain: float = 0.0,
                    spec: Optional[np.ndarray] = None, spec_gate: Optional[np.ndarray] = None,
                    lo: float = 0.02, hi: float = 0.98) -> np.ndarray:
    """Evidence-graded Bernoulli inclusion prob π, monotone in the prior weight w AND (optionally) seed rarity.
    π = expit( logit(p0) + κ·z_w + rarity_gain·z_rarity ). z_w = robust-standardized log1p(evidence);
    z_rarity = standardized seed rarity. rarity NaN (uncovered) → 0 contribution. p0 = base rate at median evidence.

    ⚠ rarity_gain DEFAULTS TO 0 (rarity OFF) — the term is NOT validated (`rarity_bench.py` + edge-coupling audit,
    2026-07-09). The idea: a rare seed = a-priori more likely a genuine SPECIFIC regulator, tiebreaking collinear
    co-seed arms the likelihood can't. BUT as a global high-gain inclusion channel it FAILS: it lifts the
    globally-rarest seeds, which are disproportionately UNEXPRESSED / low-confidence orphans with tiny (often
    artifactual) targetomes — not coupling specialists. Audit: Spearman(rarity, edge-ρ) = +0.10/+0.11/−0.01 on
    PTEN/CCND1/VEGFA (wrong sign); the lifted "specialists" are unexpressed (miR-105 RPM 0.48) or non-coupling
    (miR-99/100 q1.0, miR-718 q0.45). The legit case (miR-18a in miR-17~92) is MODERATE-rarity + expressed + in a
    coupling family. Re-enable only behind an expression gate + within-cooperative-family scope + coupling
    validation (see LEARNED_MODEL_DISCOVERY_SYNTHESIS §6a). Rarity is INCLUSION-only regardless — see slab_scale.

    `spec` is the VALIDATED replacement (`seq_specificity.py`): per-edge affinity CONCENTRATION (z-scored),
    added as spec_gain·gate·z(spec). `spec_gate` ∈(0,1] is the EXPRESSION gate — the fix rarity lacked, so an
    undetectable arm gets no boost however specific its sequence. Gated seq-spec lifts real expressed+coupling
    specialists (PTEN miR-25/92 PIP 0.52→0.86) and suppresses unexpressed ones (miR-494 gate 0.17). Default OFF
    (spec_gain=0): the lift is modest/uneven and can still promote expressed-non-coupling edges; enable when the
    discovery-precision validation (SYNTHESIS §6a role scoping) is done."""
    s = np.log1p(np.clip(np.asarray(w, dtype=float), 0, None))
    med = np.median(s)
    mad = np.median(np.abs(s - med)) * 1.4826
    z = (s - med) / mad if mad > 1e-9 else np.zeros_like(s)
    logit = np.log(p0 / (1 - p0)) + kappa * z
    if rarity is not None:                                                # rarity ORDERING channel (RETIRED; see docstring)
        r = np.asarray(rarity, dtype=float)
        rc = r[~np.isnan(r)]
        if len(rc) > 1 and rc.std() > 1e-9:
            zr = np.where(np.isnan(r), 0.0, (r - rc.mean()) / rc.std())   # NaN → neutral (0)
            logit = logit + rarity_gain * zr
    if spec is not None and spec_gain:                                    # SEQ-SPEC ordering channel (validated; §6a)
        s2 = np.asarray(spec, dtype=float); sc = s2[~np.isnan(s2)]        # affinity-concentration, EXPRESSION-GATED
        if len(sc) > 1 and sc.std() > 1e-9:
            zs = np.where(np.isnan(s2), 0.0, (s2 - sc.mean()) / sc.std())
            gate = 1.0 if spec_gate is None else np.clip(np.asarray(spec_gate, dtype=float), 0.0, 1.0)
            logit = logit + spec_gain * gate * zs                        # gate kills undetectable arms (rarity's failure)
    return np.clip(expit(logit), lo, hi)


def _biochem_magnitude(gene: str, candidates: pd.Index, *, by: str) -> pd.Series:
    """scanMiR predicted-repression magnitude per candidate (≥0), normalised within-gene to O(1).
    Family mode uses the strongest member (kd._family_affinity). Unmatched candidate → 0 (no biochem prior)."""
    from mirna_hallmark.learned import kd as KD
    try:
        if by == "family":
            fa = KD._family_affinity()
            rg = fa.loc[fa["gene"] == gene].set_index("arm")["repression"]     # 'arm' col = family label
        else:
            aff = KD.affinity()
            rg = aff.loc[aff["gene"] == gene].set_index("arm")["repression"]
    except Exception:                                                          # scanMiR cache absent → no biochem prior
        rg = pd.Series(dtype=float)
    mag = (-rg.reindex(candidates).fillna(0.0)).clip(lower=0.0)                # ≤0 logFC → ≥0 magnitude
    med = mag[mag > 0].median()
    if med and med > 0:
        mag = mag / med                                                       # within-gene O(1) scale
    return pd.Series(mag.to_numpy(dtype=float), index=candidates, name="mu")


def _confidence(w: np.ndarray, *, half: float = 3.0) -> np.ndarray:
    """Evidence-depth confidence τ ∈ (0,1], saturating: τ = w/(w+half). `half` = weight at τ=0.5."""
    w = np.clip(np.asarray(w, dtype=float), 0, None)
    return w / (w + half)


def edge_priors(gene: str, *, by: str = "family", w_prior_source: str = "ledger",
                candidates: Optional[pd.Index] = None, w_evidence: Optional[pd.Series] = None,
                p0: float = 0.3, kappa: float = 1.5, rarity_gain: float = 0.0,
                spec_gain: float = 0.0, expr: Optional[pd.Series] = None, expr_floor: float = 3.46,
                expr_scale: float = 1.5) -> pd.DataFrame:
    """π / μ / τ for a gene's candidate regulators. If `candidates`+`w_evidence` are supplied (e.g. the
    family predictors + weights already assembled by the caller) they are used directly; otherwise the
    gene's POOLED-HE candidates are assembled here (arm mode, optionally collapsed to families).
    `rarity_gain` defaults 0 — seed-rarity channel OFF (unvalidated). `spec_gain` (>0) turns on the VALIDATED
    seq-specificity channel (affinity concentration, `seq_specificity.py`); `expr` = per-candidate median log2
    expression → the EXPRESSION GATE (floor=log2(11) detectability, soft) so undetectable arms get no boost."""
    if candidates is None or w_evidence is None:
        _, X, _, w = LD.assemble_gene(gene, w_prior_source=w_prior_source)
        if by == "family":
            X, w, _ = FAM.collapse_by_family(X, w, FAM.family_of(X.columns))
        candidates, w_evidence = X.columns, w
        if expr is None:
            expr = X.median()                                             # per-arm median expr for the gate
    w = w_evidence.reindex(candidates).fillna(0.0).to_numpy(dtype=float)
    out = pd.DataFrame(index=candidates)
    out["rarity"] = _seed_rarity(candidates)                              # RETIRED channel (kept for provenance)
    out["seq_spec"] = _seq_specificity(gene, candidates)                 # VALIDATED specificity (affinity concentration)
    gate = None
    if expr is not None:
        e = expr.reindex(candidates).to_numpy(dtype=float)
        gate = expit((e - expr_floor) / expr_scale)                      # soft detectability gate ∈(0,1]
        out["expr_gate"] = gate
    out["pi"] = inclusion_prior(w, p0=p0, kappa=kappa, rarity_gain=rarity_gain, rarity=out["rarity"].to_numpy(),
                                spec_gain=spec_gain, spec=out["seq_spec"].to_numpy(), spec_gate=gate)
    out["mu"] = _biochem_magnitude(gene, candidates, by=by).to_numpy()
    out["tau"] = _confidence(w)
    out["w"] = w
    return out


def _seq_specificity(gene: str, candidates: pd.Index) -> np.ndarray:
    """Per-edge affinity-concentration specificity (HE universe). NaN where uncovered. See seq_specificity."""
    try:
        from mirna_hallmark.learned import seq_specificity as SQ
        return SQ.seq_spec(gene, candidates, universe="HE").reindex(candidates).to_numpy(dtype=float)
    except Exception:
        return np.full(len(candidates), np.nan)


def _seed_rarity(candidates: pd.Index) -> np.ndarray:
    """Seed-rarity ∈ (0,1] per candidate (rare seed → 1). Sequence-only ⇒ resolves the collinear direction the
    likelihood cannot. NaN where uncovered (STAGE-1 non-conserved) → caller treats as neutral. See seed_rarity."""
    try:
        from mirna_hallmark.learned import seed_rarity as SR
        return SR.rarity(candidates).reindex(candidates).to_numpy(dtype=float)
    except Exception:
        return np.full(len(candidates), np.nan)


def slab_scale(priors: pd.DataFrame, *, mu_gain: float = 1.0, floor: float = 0.25) -> np.ndarray:
    """Per-edge slab-scale multiplier for the sampler: 1 + mu_gain·μ·τ, floored. High biochemical affinity
    AND deep evidence ⇒ a looser slab (larger allowed magnitude); thin/weak ⇒ tighter (shrink toward 0).
    Location is untouched — this is a scale prior, not a magnitude clamp (Design 'prior sets ordering').

    Seed rarity acts on π/inclusion, NOT here — and this is empirically settled, not just a design choice
    (`rarity_bench.py`, 2026-07-09): adding rarity to the slab (evidence-gated ·τ OR standalone) moves the DENSE
    ridge's OOF coupling by ΔρOOF = +0.000 at every gain/gene and leaves attribution magnitudes identical to 3
    decimals — because at n≫p the likelihood dominates the slab prior (a loosened prior variance the data already
    pins down does nothing; miR-105 slab 1.36→2.96 ⇒ ΔM 3e-5). The slab has leverage only on selection, and there
    inclusion π (a discrete decision the prior can tip) is the correct locus. So rarity stays out of the slab."""
    m = (1.0 + mu_gain * priors["mu"].to_numpy(dtype=float) * priors["tau"].to_numpy(dtype=float))
    return np.clip(m, floor, None)


def _report(gene: str = "PTEN") -> None:
    p = edge_priors(gene)
    print(f"[priors] {gene}: {len(p)} family candidates")
    print(f"[priors]  pi  range {p['pi'].min():.2f}-{p['pi'].max():.2f} (median {p['pi'].median():.2f})")
    print(f"[priors]  mu  nonzero {int((p['mu']>0).sum())}/{len(p)} (max {p['mu'].max():.2f})")
    print(f"[priors]  tau range {p['tau'].min():.2f}-{p['tau'].max():.2f}")
    p = p.assign(slab=slab_scale(p))
    print(p.sort_values("pi", ascending=False).head(10).to_string(
        float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    import sys
    _report(sys.argv[1] if len(sys.argv) > 1 else "PTEN")
