"""Phase-2 GENERATOR for the canonical architecture map (ARCHITECTURE_PLAN.md).

Joins the three dimensions BY REFERENCE — never copying content — and emits `docs/ARCHITECTURE.md`:
  * RESULTS  : `docs/derived/axis_assignment.tsv` (reviewed MH→axis) × the registry headlines.
  * ANALYSES : every module, auto-tagged to an axis by its docstring (first pass; reviewable later).
  * MODELS   : the per-axis model note (hand-curated below — the only prose the map holds, one line/axis).
The map is a VIEW: every finding is an `MH-##` link, every module a path — regenerate, never hand-edit.

Run: .venv/bin/python3 -m mirna_hallmark.analyses.ops.gen_architecture
"""
from __future__ import annotations

import glob
import re
from collections import defaultdict

from mirna_hallmark import config as C
from mirna_hallmark.analyses.ops.axis_bootstrap import classify

ROOT = C.SUBPROJECT_DIR
REG = ROOT / "docs" / "DISCOVERY_REGISTRY.md"
ASSIGN = ROOT / "docs" / "derived" / "axis_assignment.tsv"
OUT = ROOT / "docs" / "ARCHITECTURE.md"

# The closed axis taxonomy (Phase 0): tag → (title, scientific question, one-line MODEL note). Ordered.
AXES = [
    ("model", "The learned model", "what the estimator IS",
     "ONE dense learned-τ² non-negative Gibbs posterior/gene; two readouts (π≡1 dense→coupling; evidence-π→discovery). Lasso RETIRED."),
    ("edge-existence", "Edge existence & coupling", "do the curated/predicted edges exist AND act (coupling|C)?",
     "partial-Spearman(arm, mRNA | C) via the Gibbs posterior; C = CPE+target_cn+mal_prolif (+8 deconv fractions)."),
    ("cn-causal", "CN-dose causal", "CN-dose causal identification (instrument, exclusion)",
     "CN-locus 2SLS instrument + Hansen-J over-ID + T1 screen. ⛔ instruments RETRACTED; DESIGN live (highest-value open item)."),
    ("attribution", "Attribution / identity", "WHO owns a gene's regulation — identity vs magnitude?",
     "Shapley/LMG on R² (identity, collinearity-fair) vs Gibbs β (magnitude). `share`=β_f/Σβ is NOT identity."),
    ("decoy", "Decoy / specificity", "does curation beat an abundance/variance-matched fake?",
     "site-free / abundance+variance-matched decoy; Hungarian on signed loadings+dose+variance. Gap ≈ −0.012 (two designs agree)."),
    ("discovery", "Discovery", "novel edges beyond curation (its own lane)",
     "site-free empirical null (heavy-tailed) → Simes-within-family → BH; per-edge EMPTY, signal is set-level + convergent-evidence."),
    ("protein", "CPTAC / protein", "does it hold at the protein layer?",
     "CPTAC prospective (n=101); β_source = a·β_target gauge; composition confound reframe (epithelial survive, stromal collapse)."),
    ("progression", "Progression / state", "GTEx→NAT→tumor trajectory / state",
     "nested state model M_t = a·M_H + Δ; GTEx-healthy anchor (NAT is field-effect-contaminated)."),
    ("subtype", "Subtype", "PAM50-stratified coupling / who-is-pressured",
     "per-PAM50 coupling test (learned β). ⬜ full subtype contrasts + subtype-stratified discovery = board §Y.2/#1."),
    ("outcome", "Outcome / prognostic", "does it predict survival?",
     "survival on pressure/β. Pressure magnitude non-prognostic (data-bound: TCGA=OS only). ⬜ learned-β prognostic untested."),
    ("external", "External cohorts", "independent-cohort replication",
     "Buffa / METABRIC / SCAN-B / CCLE. Cross-cohort loss ~80% (cohort boundary), ~0% (mRNA→protein layer)."),
    ("dcis-ev", "DCIS / EV / pre-malignant", "the DCIS / extracellular-vesicle lane",
     "in-situ→invasive coupling; selective EV export; MH-55 corroborates MH-114 16 days early."),
]
AXIS_TITLE = {t: (title, q, m) for t, title, q, m in AXES}


def _mh_headlines() -> dict:
    reg = REG.read_text(encoding="utf-8")
    rows = re.findall(r'^\| (MH-[\w/]+) \|(.+?)(?=\n\| MH-|\n\Z)', reg, re.S | re.M)
    out = {}
    for mh, body in rows:
        strength = ""
        m = re.search(r'\| ([SAPR](?:/[SAPR])*) \|', body)          # the strength-tag column
        if m:
            strength = m.group(1)
        h = re.sub(r'[*`⭐⛔⚠✅🔨]', '', body).strip()
        h = re.sub(r'\s+', ' ', h)
        # trim to the first sentence-ish
        h = re.split(r'(?<=[.—])\s', h)[0][:150]
        out[mh] = (h, strength)
    return out


MODULE_AXIS = ROOT / "docs" / "derived" / "module_axis.tsv"


def _modules_by_axis() -> dict:
    """Module→axis from the persisted, reviewable `module_axis.tsv` (built by axis_bootstrap logic + curated
    `shared`/`model` overrides; every module homed, incl. a `shared` infrastructure bucket). Edit that file to
    correct a tag, then regenerate — the map never re-guesses."""
    by = defaultdict(list)
    for line in MODULE_AXIS.read_text().splitlines()[1:]:
        mod, ax, *_ = line.split("\t")
        by[ax].append(mod)
    return by


def run() -> None:
    assign = defaultdict(list)
    for line in ASSIGN.read_text().splitlines()[1:]:
        mh, ax, *_ = line.split("\t")
        assign[ax].append(mh)
    heads = _mh_headlines()
    mods = _modules_by_axis()

    L = ["# Canonical architecture — axes × models × analyses × results",
         "",
         "> **GENERATED** by `analyses/ops/gen_architecture.py` — do NOT hand-edit. A materialized VIEW over the",
         "> one-home docs (registry, catalog, STATE_OF_PLAY): every finding is an `MH-##` (→ registry), every",
         "> module a path. Regenerate after tagging changes. Join key: `docs/derived/axis_assignment.tsv`.",
         f"> Axes: **{len(AXES)}** · findings mapped: **{sum(len(v) for v in assign.values())}** · "
         f"modules tagged: **{sum(len(v) for v in mods.values())}**.", ""]
    for tag, title, q, model in AXES:
        mhs = assign.get(tag, [])
        ms = mods.get(tag, [])
        L.append(f"## {title}  `{tag}`")
        L.append(f"*{q}*")
        L.append("")
        L.append(f"- **model:** {model}")
        L.append(f"- **status → verdict:** see [`STATE_OF_PLAY.md`](STATE_OF_PLAY.md) · **{len(mhs)}** findings · **{len(ms)}** modules")
        if mhs:
            L.append(f"- **results ({len(mhs)}):** " + " · ".join(f"[{mh}](DISCOVERY_REGISTRY.md)" for mh in mhs))
        if ms:
            show = ms[:14]
            L.append(f"- **analyses ({len(ms)}):** " + ", ".join(f"`{m}`" for m in show) + (" …" if len(ms) > 14 else ""))
        L.append("")
    sh = mods.get("shared", [])
    if sh:
        L.append("## Shared / infrastructure  `shared`")
        L.append("*axis-agnostic loaders, builders, batch/CN helpers, and the retired-pressure builders — serve every axis*")
        L.append("")
        L.append(f"- **modules ({len(sh)}):** " + ", ".join(f"`{m}`" for m in sh[:20]) + (" …" if len(sh) > 20 else ""))
        L.append("")
    OUT.write_text("\n".join(L))
    _emit_json(assign, heads, mods)                          # data for the artifact view (rendered separately)
    print(f"[gen-architecture] wrote {OUT}: {len(AXES)} axes, "
          f"{sum(len(v) for v in assign.values())} findings, {sum(len(v) for v in mods.values())} modules tagged")


# per-axis status (curated from STATE_OF_PLAY — the ONE line of prose not derivable by join)
STATUS = {"model": "strong", "edge-existence": "mixed", "cn-causal": "design-only", "attribution": "strong",
          "decoy": "strong", "discovery": "strong", "protein": "mixed", "progression": "strong",
          "subtype": "open", "outcome": "data-bound", "external": "strong", "dcis-ev": "strong"}


def _emit_json(assign, heads, mods) -> None:
    """Emit docs/derived/architecture_data.json — the data the HTML artifact renders (regenerate → re-publish
    architecture.html, same file path keeps the URL)."""
    import json
    axes = []
    for tag, title, q, model in AXES:
        fs = [{"mh": mh, "h": heads.get(mh, ("", ""))[0][:120], "s": heads.get(mh, ("", ""))[1]}
              for mh in assign.get(tag, [])]
        axes.append({"tag": tag, "title": title, "q": q, "model": model, "status": STATUS.get(tag, ""),
                     "findings": fs, "modules": sorted(mods.get(tag, []))})
    if mods.get("shared"):
        axes.append({"tag": "shared", "title": "Shared / infrastructure", "status": "shared",
                     "q": "axis-agnostic loaders, builders, batch/CN helpers, retired-pressure builders",
                     "model": "not an axis — cross-cutting infrastructure every axis depends on.",
                     "findings": [], "modules": sorted(mods["shared"])})
    out = {"axes": axes, "n_ax": len([a for a in axes if a["tag"] != "shared"]),
           "n_find": sum(len(a["findings"]) for a in axes),
           "n_mod": sum(len(a["modules"]) for a in axes)}
    (ROOT / "docs" / "derived" / "architecture_data.json").write_text(json.dumps(out))


if __name__ == "__main__":
    run()
