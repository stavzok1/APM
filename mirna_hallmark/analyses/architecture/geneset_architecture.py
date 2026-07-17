"""Stage 3 (extension 2b, pilot): gene-set architecture-weighted miRNA pressure.

The Hallmark/gene-set pressure used elsewhere treats a set as a FLAT bag of genes
(mean / sum over members). This pilot overlays the set's regulatory **topology** so
that pressure on a node counts according to its role:

  - **master-regulator weight** = reverse-PageRank on the set's induced directed
    signed network (a node that influences many downstream nodes scores high);
  - **sign-aware effect** = a node whose out-edges are mostly inhibitory is an
    *inhibitor*; a miRNA repressing an inhibitor *de-represses* the program, so its
    effective pressure sign on the program output is flipped (double-negative);
  - **redundancy discount** = nodes sharing most of their downstream targets
    (Jaccard >= 0.5) split their weight (one of several redundant effectors matters
    less than a non-redundant bottleneck).

Network source: **OmniPath** REST (`https://omnipathdb.org/interactions`, directed +
signed; aggregates SIGNOR/SignaLink/... ) fetched via urllib and cached locally (no
heavy package / no networkx dependency; PageRank is a small numpy power iteration).

Pilot systems (both are MSigDB Hallmark sets, so the Stage-1 per-edge contributions
already exist): EMT (`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`, the miR-200/ZEB
axis) and PI3K/AKT/PTEN (`HALLMARK_PI3K_AKT_MTOR_SIGNALING`). For each, flat-mean vs
architecture-weighted miRNA pressure ranking + the rank shift.

Run:
  .venv/bin/python3 -m mirna_hallmark.geneset_architecture
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from mirna_hallmark import config as C
from mirna_hallmark.hallmark_sets import HallmarkSets
from mirna_hallmark.mirna_state_class import OUT_DIR as STATE_DIR

OUT_DIR = C.TISSUE_REFERENCE_DIR / "geneset_architecture"
CACHE_DIR = OUT_DIR / "network_cache"
PILOT_SETS = ["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_PI3K_AKT_MTOR_SIGNALING"]
OMNIPATH_URL = "https://omnipathdb.org/interactions"
REDUNDANCY_JACCARD = 0.5
PAGERANK_D = 0.85
# signed-influence propagation: finite-horizon (Katz-style) damped sum over signed paths
PROP_ALPHA = 0.25        # per-hop damping (< 1/spectral_radius keeps the series convergent)
PROP_HORIZON = 6         # max path length summed
# curated-TF hierarchy boost: w_hier = w_arch · (1 + TF_HIER_BOOST · is_tf). A MANUAL, tunable
# up-weight of genes that are curated human TFs (Lambert 2018), layered on the topology w_arch — the
# topology only weakly tracks TF identity (≈48% of curated TFs are below median w_arch). 0 = no boost.
TF_HIER_BOOST = 0.5

# --------------------------------------------------------------------------- #
# Tumor-direction prior (coarse BC-context literature polarity per Hallmark)
#   +1 = program whose UP-regulation is pro-tumorigenic in breast cancer
#   -1 = program whose UP-regulation is tumor-SUPPRESSIVE / anti-tumor (immune, death)
#    0 = ambiguous / context-dependent / metabolic-differentiation (excluded from the
#        pro-tumor call, NOT from the topology score)
# These are deliberately conservative; ambiguous sets are 0 rather than guessed.
# --------------------------------------------------------------------------- #
TUMOR_DIRECTION_PRIOR: Dict[str, int] = {
    # pro-tumorigenic (+1)
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": +1,
    "HALLMARK_HYPOXIA": +1,
    "HALLMARK_GLYCOLYSIS": +1,
    "HALLMARK_MTORC1_SIGNALING": +1,
    "HALLMARK_PI3K_AKT_MTOR_SIGNALING": +1,
    "HALLMARK_E2F_TARGETS": +1,
    "HALLMARK_G2M_CHECKPOINT": +1,
    "HALLMARK_MYC_TARGETS_V1": +1,
    "HALLMARK_MYC_TARGETS_V2": +1,
    "HALLMARK_MITOTIC_SPINDLE": +1,
    "HALLMARK_ANGIOGENESIS": +1,
    "HALLMARK_TGF_BETA_SIGNALING": +1,
    "HALLMARK_WNT_BETA_CATENIN_SIGNALING": +1,
    "HALLMARK_NOTCH_SIGNALING": +1,
    "HALLMARK_KRAS_SIGNALING_UP": +1,
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE": +1,
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": +1,
    # tumor-suppressive / anti-tumor (-1)
    "HALLMARK_APOPTOSIS": -1,
    "HALLMARK_P53_PATHWAY": -1,
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": -1,
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": -1,
    "HALLMARK_INFLAMMATORY_RESPONSE": -1,
    "HALLMARK_ALLOGRAFT_REJECTION": -1,
}


def _prior(sname: str) -> int:
    return int(TUMOR_DIRECTION_PRIOR.get(sname, 0))


# --------------------------------------------------------------------------- #
# OmniPath signed directed network (cached)
# --------------------------------------------------------------------------- #
def _fetch_omnipath(genes: Sequence[str], *, cache: Path, force: bool = False) -> pd.DataFrame:
    if cache.is_file() and not force:
        return pd.read_csv(cache, sep="\t")
    cache.parent.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    # chunk the partners list to keep the GET URL within limits
    glist = sorted(set(genes))
    for i in range(0, len(glist), 60):
        chunk = glist[i:i + 60]
        q = OMNIPATH_URL + "?" + urllib.parse.urlencode({
            "partners": ",".join(chunk), "genesymbols": "yes",
            "fields": "is_directed,is_stimulation,is_inhibition", "format": "tsv",
        })
        try:
            txt = urllib.request.urlopen(q, timeout=90).read().decode()
        except Exception as exc:  # transient network / 5xx: skip chunk, keep going
            print(f"[architecture]   WARN OmniPath chunk {i // 60} failed ({exc}); skipping")
            continue
        rows.extend(csv.DictReader(io.StringIO(txt), delimiter="\t"))
    df = pd.DataFrame(rows)
    if df.empty:
        df = pd.DataFrame(columns=["source_genesymbol", "target_genesymbol", "is_directed",
                                   "is_stimulation", "is_inhibition"])
    df.to_csv(cache, sep="\t", index=False)
    return df


_EMPTY_EDGES = pd.DataFrame(columns=["source_genesymbol", "target_genesymbol", "stim", "inh", "net_sign"])


def _induced_edges(net: pd.DataFrame, members: set) -> pd.DataFrame:
    # tolerate empty / malformed / partial OmniPath responses: a set with no usable
    # signed-directed edges simply gets a flat (uniform-weight) architecture.
    if net is None or net.empty or not {"source_genesymbol", "target_genesymbol"}.issubset(net.columns):
        return _EMPTY_EDGES.copy()
    def _b(v) -> bool:
        return str(v).strip().lower() in ("true", "1", "yes")
    e = net.copy()
    e = e[e["source_genesymbol"].isin(members) & e["target_genesymbol"].isin(members)]
    if "is_directed" in e.columns:
        e = e[e["is_directed"].map(_b)]
    e["stim"] = e["is_stimulation"].map(_b) if "is_stimulation" in e.columns else False
    e["inh"] = e["is_inhibition"].map(_b) if "is_inhibition" in e.columns else False
    if e.empty:
        return _EMPTY_EDGES.copy()
    # CALIBRATED edge sign: instead of collapsing each (source,target) to a hard +1/-1/0,
    # aggregate ALL OmniPath records for the pair into a continuous net_sign in [-1, +1]:
    #   net_sign = (n_stimulation - n_inhibition) / (n_stimulation + n_inhibition)
    # A pair that is purely stimulatory -> +1, purely inhibitory -> -1, dual / conflicting
    # (the "ambiguous" case that used to be flattened to 0) -> a graded value reflecting how
    # the evidence leans, so repressing such a node carries a calibrated (not zeroed) effect.
    e["_st"] = e["stim"].astype(int)
    e["_ih"] = e["inh"].astype(int)
    agg = (e.groupby(["source_genesymbol", "target_genesymbol"], as_index=False)
           .agg(n_stim=("_st", "sum"), n_inh=("_ih", "sum")))
    denom = (agg["n_stim"] + agg["n_inh"]).replace(0, np.nan)
    agg["net_sign"] = ((agg["n_stim"] - agg["n_inh"]) / denom).fillna(0.0)
    agg["stim"] = agg["net_sign"] > 0
    agg["inh"] = agg["net_sign"] < 0
    return agg[["source_genesymbol", "target_genesymbol", "stim", "inh", "net_sign", "n_stim", "n_inh"]]


# --------------------------------------------------------------------------- #
# Topology -> per-gene architecture weights
# --------------------------------------------------------------------------- #
def _pagerank(A: np.ndarray, *, d: float = PAGERANK_D, iters: int = 200) -> np.ndarray:
    n = A.shape[0]
    if n == 0:
        return np.array([])
    out = A.sum(axis=1)
    M = np.where(out[:, None] > 0, A / np.where(out[:, None] == 0, 1.0, out[:, None]), 1.0 / n)
    Mt = M.T
    pr = np.ones(n) / n
    for _ in range(iters):
        pr = (1 - d) / n + d * Mt.dot(pr)
        s = pr.sum()
        if s > 0:
            pr /= s
    return pr


def _signed_influence(nodes: Sequence[str], edges: pd.DataFrame, *,
                      alpha: float = PROP_ALPHA, horizon: int = PROP_HORIZON) -> np.ndarray:
    """Net signed downstream influence of *increasing* each node on total program output.

    Builds a signed adjacency S (S[i,j] = +1 activating, -1 inhibiting; unsigned/ambiguous
    edges contribute 0 sign) and sums damped signed path products up to `horizon` hops:
        T = sum_{k=1..H} (alpha*S)^k
    output_influence(g) = 1 + sum_j T[g, j]   (self + net effect on the rest of the set).
    A central activator scores >0; an inhibitor hub (represses many) scores <0. Repressing
    such a node (a miRNA) therefore flips the program-output direction accordingly.
    """
    idx = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)
    if n == 0:
        return np.array([])
    S = np.zeros((n, n))
    has_net = "net_sign" in edges.columns
    for _, e in edges.iterrows():
        s, t = e["source_genesymbol"], e["target_genesymbol"]
        if s not in idx or t not in idx:
            continue
        if has_net:
            # calibrated continuous sign in [-1,+1]: ambiguous/dual edges contribute a
            # graded value (not a hard 0), so a node's repression carries a calibrated effect
            S[idx[s], idx[t]] += float(e["net_sign"])
        else:
            st, ih = bool(e["stim"]), bool(e["inh"])
            if st and not ih:
                S[idx[s], idx[t]] += 1.0
            elif ih and not st:
                S[idx[s], idx[t]] -= 1.0
    T = np.zeros((n, n))
    P = np.eye(n)
    for _ in range(horizon):
        P = alpha * P.dot(S)
        T += P
    return 1.0 + T.sum(axis=1)


def _gene_weights(members: Sequence[str], edges: pd.DataFrame) -> pd.DataFrame:
    nodes = sorted(members)
    idx = {g: i for i, g in enumerate(nodes)}
    n = len(nodes)
    A = np.zeros((n, n))
    out_targets: Dict[str, set] = {g: set() for g in nodes}
    out_inh: Dict[str, int] = {g: 0 for g in nodes}
    out_tot: Dict[str, int] = {g: 0 for g in nodes}
    out_netsign: Dict[str, float] = {g: 0.0 for g in nodes}
    has_net = "net_sign" in edges.columns
    for _, e in edges.iterrows():
        s, t = e["source_genesymbol"], e["target_genesymbol"]
        if s in idx and t in idx:
            A[idx[s], idx[t]] = 1.0
            out_targets[s].add(t)
            out_tot[s] += 1
            if bool(e["inh"]) and not bool(e["stim"]):
                out_inh[s] += 1
            if has_net:
                out_netsign[s] += float(e["net_sign"])
    # master-regulator weight = reverse-PageRank (influence / upstreamness)
    rev_pr = _pagerank(A.T)
    pr = _pagerank(A)
    rev_pr = rev_pr if rev_pr.size else np.ones(n) / n
    # effect sign: predominantly-inhibitory out-edges => inhibitor (-1). The BINARY call is a
    # coarse majority vote; `effect_sign_soft` is the calibrated mean net_sign over the node's
    # out-edges in [-1,+1] (ambiguous / dual-role nodes land near 0 instead of being forced to
    # +/-1), so a node with conflicting out-edges carries a graded, not flattened, sign.
    sign = []
    sign_soft = []
    for g in nodes:
        frac = out_inh[g] / out_tot[g] if out_tot[g] else 0.0
        sign.append(-1 if frac > 0.5 else 1)
        sign_soft.append(out_netsign[g] / out_tot[g] if out_tot[g] else 0.0)
    # redundancy groups by downstream-target Jaccard
    group = list(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            ti, tj = out_targets[nodes[i]], out_targets[nodes[j]]
            if not ti or not tj:
                continue
            jac = len(ti & tj) / len(ti | tj)
            if jac >= REDUNDANCY_JACCARD:
                gj = group[j]
                gi = group[i]
                for k in range(n):
                    if group[k] == gj:
                        group[k] = gi
    grp_size = pd.Series(group).value_counts().to_dict()
    df = pd.DataFrame({
        "gene": nodes,
        "centrality_revpagerank": rev_pr,
        "pagerank_downstream": pr if pr.size else np.ones(n) / n,
        "out_degree": [out_tot[g] for g in nodes],
        "in_degree": [int(A[:, idx[g]].sum()) for g in nodes],
        "effect_sign": sign,
        "effect_sign_soft": sign_soft,
        "redundancy_group": group,
    })
    df["redundancy_size"] = df["redundancy_group"].map(grp_size)
    # architecture weight: centrality damped by redundancy, normalized so mean = 1
    raw = (df["centrality_revpagerank"] + 1.0 / n) / df["redundancy_size"]
    df["w_arch"] = raw / raw.mean()
    # net signed downstream influence on program output (topology-propagated)
    df["output_influence"] = _signed_influence(nodes, edges)
    # per-gene MALIGNANCY role (oncogene +1 / TSG -1 / dual-unknown 0) - lets repression of a
    # node be scored for malignancy direction INDEPENDENTLY of the gene-gene topology (which is
    # blind to PTEN's lipid phosphatase, etc.). repressing a TSG is pro-tumor, an oncogene anti-.
    try:
        from mirna_hallmark.gene_roles import load_gene_roles
        roles = load_gene_roles(nodes).set_index("gene")["malignancy_sign"]
        df["malignancy_sign"] = df["gene"].map(roles).fillna(0).astype(int)
    except Exception:
        df["malignancy_sign"] = 0
    # CONTINUOUS, ALL-GENE role weight (DepMap CRISPR Chronos gene-effect, breast panel) — the graded
    # counterpart to the binary malignancy_sign, defined for ~18.5k genes not just the 232 drivers.
    # dep_role_weight is already oriented like -malignancy_sign (>0 = repressing g is pro-tumor / TSG-like,
    # <0 = repressing g hits a dependency / anti-tumor). One-sided: resolves the dependency (anti-tumor)
    # tail well, the TSG (pro-tumor) tail weakly — so it is the all-gene graded onco-collateral axis.
    try:
        from mirna_hallmark.gene_roles import load_gene_dependency
        dep = load_gene_dependency(nodes).set_index("gene")["dep_role_weight"]
        df["dep_role_weight"] = df["gene"].map(dep).fillna(0.0)
    except Exception:
        df["dep_role_weight"] = 0.0
    # CURATED TF / master-regulator hierarchy (Lambert 2018) — an identity-based annotation that the
    # topology w_arch only weakly captures. `is_tf` flags curated human TFs; `w_hier` is w_arch with a
    # MANUAL, tunable up-weight for TFs (TF_HIER_BOOST). w_arch is left untouched so the malignancy_arch
    # scores stay pure-topology; w_hier is an opt-in overlay for "regulator-level" weighting.
    try:
        from mirna_hallmark.gene_roles import load_tf_census
        tfc = load_tf_census(nodes).set_index("gene")["is_tf"]
        df["is_tf"] = df["gene"].map(tfc).fillna(False).astype(bool)
    except Exception:
        df["is_tf"] = False
    df["w_hier"] = df["w_arch"] * (1.0 + TF_HIER_BOOST * df["is_tf"].astype(float))
    return df


# --------------------------------------------------------------------------- #
# Flat vs architecture-weighted miRNA pressure
# --------------------------------------------------------------------------- #
def _mirna_pressure(edge: pd.DataFrame, members: set, weights: pd.DataFrame,
                    prior: int = 0) -> pd.DataFrame:
    sub = edge[edge["gene"].isin(members)].copy()
    if sub.empty or "c_tumor" not in sub.columns:
        return pd.DataFrame()
    w = weights.set_index("gene")
    sub["w_arch"] = sub["gene"].map(w["w_arch"]).fillna(1.0)
    sub["effect_sign"] = sub["gene"].map(w["effect_sign"]).fillna(1)
    sub["output_influence"] = sub["gene"].map(w["output_influence"]).fillna(1.0)
    sub["malignancy_sign"] = sub["gene"].map(w["malignancy_sign"]).fillna(0) if "malignancy_sign" in w.columns else 0
    sub["dep_role_weight"] = sub["gene"].map(w["dep_role_weight"]).fillna(0.0) if "dep_role_weight" in w.columns else 0.0
    sub["w_hier"] = sub["gene"].map(w["w_hier"]).fillna(sub["w_arch"]) if "w_hier" in w.columns else sub["w_arch"]
    sub["is_tf"] = sub["gene"].map(w["is_tf"]).fillna(False).astype(bool) if "is_tf" in w.columns else False
    sub["c_arch"] = sub["c_tumor"] * sub["w_arch"]
    sub["c_arch_signed"] = sub["c_arch"] * sub["effect_sign"]
    # MALIGNANCY-weighted: repressing a TSG (sign -1) credits +c toward malignancy; an oncogene
    # (+1) credits -c. Independent of the gene-gene topology, so it captures miR-21->PTEN/PDCD4.
    sub["mal_change"] = -sub["malignancy_sign"] * sub["c_tumor"]
    # ARCHITECTURE-calibrated malignancy: gene-role DIRECTION (the only sign) × edge strength
    # (c_tumor) × program-importance weight (w_arch = master-regulator/centrality, ALWAYS > 0).
    # We deliberately do NOT multiply by the signed output_influence: that double-counts direction
    # (a TSG with negative output_influence would flip the role credit) and lets a single hub edge
    # (e.g. TP53, w_arch=21 x output_influence=10) explode the score. w_arch alone answers the
    # "weight of the gene relative to the program (master TF/inhibitor)" requirement cleanly.
    sub["mal_arch_change"] = -sub["malignancy_sign"] * sub["c_tumor"] * sub["w_arch"]
    sub["tsg_credit"] = np.where(sub["malignancy_sign"] < 0,
                                 sub["c_tumor"] * sub["w_arch"], 0.0)
    sub["onco_collateral"] = np.where(sub["malignancy_sign"] > 0,
                                      sub["c_tumor"] * sub["w_arch"], 0.0)
    # topology-propagated: repressing g (miRNA, c>0) changes program output by -c * influence(g)
    sub["prog_change"] = -sub["c_tumor"] * sub["output_influence"]
    # ACQUIRED variant (parallel to abundance): weight by the per-edge acquired pressure
    # gain max(c_tumor - c_gtex, 0) (collapse-repaired, magnitude-aware) instead of total
    # tumor abundance -> isolates the pro-tumor pressure that is TUMOR-ACQUIRED vs healthy.
    if "c_gtex" in sub.columns:
        sub["acq_gain"] = (pd.to_numeric(sub["c_tumor"], errors="coerce")
                           - pd.to_numeric(sub["c_gtex"], errors="coerce")).clip(lower=0.0).fillna(0.0)
    else:
        sub["acq_gain"] = sub["c_tumor"]
    sub["prog_change_acq"] = -sub["acq_gain"] * sub["output_influence"]
    sub["mal_change_acq"] = -sub["malignancy_sign"] * sub["acq_gain"]
    sub["mal_arch_change_acq"] = -sub["malignancy_sign"] * sub["acq_gain"] * sub["w_arch"]
    # CONTINUOUS all-gene malignancy (dep_role_weight already = repression-credit orientation, so no
    # extra minus): the graded analogue of mal_change over EVERY gene in the set, not just curated drivers.
    sub["mal_change_cont"] = sub["dep_role_weight"] * sub["c_tumor"]
    sub["mal_change_cont_acq"] = sub["dep_role_weight"] * sub["acq_gain"]
    # HIERARCHY-weighted malignancy: same gene-role direction × edge strength, but the program-importance
    # weight is w_hier (= w_arch up-weighted for curated TFs) instead of w_arch, so repressing a curated
    # master-regulator counts MORE than an equally-central non-TF effector. Parallel to mal_arch_change.
    sub["mal_change_hier"] = -sub["malignancy_sign"] * sub["c_tumor"] * sub["w_hier"]
    sub["mal_change_hier_acq"] = -sub["malignancy_sign"] * sub["acq_gain"] * sub["w_hier"]
    # HIERARCHY-weighted program-output change: topology damage with curated-TF master regulators
    # up-weighted (w_hier/w_arch ratio applied to the propagated influence). Isolates pressure that
    # routes through transcriptional masters vs generic effectors.
    sub["hier_boost"] = (sub["w_hier"] / sub["w_arch"]).replace([np.inf, -np.inf], 1.0).fillna(1.0)
    sub["prog_change_hier"] = -sub["c_tumor"] * sub["output_influence"] * sub["hier_boost"]
    sub["tf_pressure"] = np.where(sub["is_tf"], sub["c_tumor"], 0.0)  # raw pressure landing on curated TFs
    agg = sub.groupby("miRNA").agg(
        n_targets=("gene", "nunique"),
        flat_pressure=("c_tumor", "sum"),
        arch_pressure=("c_arch", "sum"),
        arch_signed_pressure=("c_arch_signed", "sum"),
        program_output_change=("prog_change", "sum"),
        acquired_pressure=("acq_gain", "sum"),
        program_output_change_acquired=("prog_change_acq", "sum"),
        mal_pro_tumor=("mal_change", "sum"),
        mal_pro_tumor_acquired=("mal_change_acq", "sum"),
        mal_pro_tumor_arch=("mal_arch_change", "sum"),
        mal_pro_tumor_arch_acquired=("mal_arch_change_acq", "sum"),
        mal_pro_tumor_cont=("mal_change_cont", "sum"),
        mal_pro_tumor_cont_acquired=("mal_change_cont_acq", "sum"),
        mal_pro_tumor_hier=("mal_change_hier", "sum"),
        mal_pro_tumor_hier_acquired=("mal_change_hier_acq", "sum"),
        program_output_change_hier=("prog_change_hier", "sum"),
        tf_pressure=("tf_pressure", "sum"),
        tsg_credit=("tsg_credit", "sum"),
        onco_collateral=("onco_collateral", "sum"),
        n_tsg_targets=("malignancy_sign", lambda s: int((s < 0).sum())),
        n_onco_targets=("malignancy_sign", lambda s: int((s > 0).sum())),
        n_dep_targets=("dep_role_weight", lambda s: int((s < -0.5).sum())),
        n_tf_targets=("is_tf", lambda s: int(s.sum())),
    ).reset_index()
    # combine topology-propagated program change with the tumor-direction prior:
    #   pro_tumor_pressure > 0  -> the arm's pressure pushes this program toward malignancy
    #   pro_tumor_pressure < 0  -> toward suppression (damages a pro-tumor program /
    #                              releases an anti-tumor one)
    agg["tumor_direction_prior"] = int(prior)
    agg["pro_tumor_pressure"] = prior * agg["program_output_change"]                  # abundance-weighted
    agg["pro_tumor_acquired"] = prior * agg["program_output_change_acquired"]         # acquired-weighted (parallel)
    agg["flat_rank"] = agg["flat_pressure"].rank(ascending=False).astype(int)
    agg["arch_rank"] = agg["arch_pressure"].rank(ascending=False).astype(int)
    agg["rank_shift"] = agg["flat_rank"] - agg["arch_rank"]  # +ve = promoted by architecture
    return agg.sort_values("arch_pressure", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Family-level rollup (seed-approx) + total-layer pressure per set
# --------------------------------------------------------------------------- #
_FAM_RE = re.compile(r"(miR|let|mir)-(\d+)([a-z]*)(?:-(\d+))?(?:-(3p|5p))?$")


def _family_key(arm: str) -> str:
    """Seed-APPROX family key from the arm name: collapse the alphabetic paralog suffix and
    copy-number index, keep the numeric stem + the mature arm (5p/3p). e.g.
    hsa-miR-30a-5p / -30b-5p / -30e-5p -> miR-30-5p; hsa-let-7a-5p -> let-7-5p;
    hsa-miR-200c-3p -> miR-200-3p. NOTE: this is a NAME-stem proxy for the TargetScan seed
    family (e.g. miR-200a vs miR-200b/c actually differ by one seed nt) - documented as such."""
    a = str(arm).replace("hsa-", "")
    m = _FAM_RE.match(a)
    if not m:
        return a
    fam = f"{m.group(1)}-{m.group(2)}"
    return f"{fam}-{m.group(5)}" if m.group(5) else fam


def _family_cross_set(all_pressure: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Family-level analog of the per-miRNA cross-set summary: sum each seed-approx family's
    pressure within every Hallmark set (long form) and roll up the pro-tumor breadth/magnitude
    across polarity sets (wide form). Answers 'what is the impact of a whole miRNA FAMILY?'."""
    if all_pressure.empty:
        return pd.DataFrame(), pd.DataFrame()
    ap = all_pressure.copy()
    ap["family"] = ap["miRNA"].map(_family_key)
    agg_cols = {c: "sum" for c in ["flat_pressure", "arch_pressure", "arch_signed_pressure",
                                   "program_output_change", "acquired_pressure",
                                   "program_output_change_acquired", "pro_tumor_pressure",
                                   "pro_tumor_acquired", "mal_pro_tumor_arch",
                                   "mal_pro_tumor_hier", "tf_pressure", "tsg_credit",
                                   "onco_collateral"] if c in ap.columns}
    long = (ap.groupby(["family", "hallmark", "tumor_direction_prior"], as_index=False)
            .agg(n_arms=("miRNA", "nunique"), **{k: (k, v) for k, v in agg_cols.items()}))
    pol = long[long["tumor_direction_prior"] != 0]
    if pol.empty:
        return long, pd.DataFrame()
    gp = pol.groupby("family")
    wide = pd.DataFrame({
        "n_arms_total": ap.groupby("family")["miRNA"].nunique(),
        "n_polarity_sets": gp.size(),
        "n_sets_pro_tumor": gp["pro_tumor_pressure"].apply(lambda s: int((s > 0).sum())),
        "n_sets_anti_tumor": gp["pro_tumor_pressure"].apply(lambda s: int((s < 0).sum())),
        "sum_pro_tumor_pressure": gp["pro_tumor_pressure"].sum(),
        "sum_pro_tumor_acquired": gp["pro_tumor_acquired"].sum() if "pro_tumor_acquired" in pol.columns else np.nan,
    }).reset_index()
    top = pol.sort_values("pro_tumor_pressure", ascending=False).drop_duplicates("family", keep="first")
    wide = wide.merge(top[["family", "hallmark"]].rename(columns={"hallmark": "top_pro_tumor_set"}),
                      on="family", how="left")
    return long, wide.sort_values("sum_pro_tumor_pressure", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def _cross_set_summary(all_pressure: pd.DataFrame) -> pd.DataFrame:
    """Per miRNA, how its architecture-weighted role distributes across Hallmark sets.

    arch_signed < 0 = the arm mainly represses *inhibitor* nodes -> net pathway-ACTIVATING;
    arch_signed > 0 = net pathway-repressing. We count sets in each direction and the mean."""
    if all_pressure.empty:
        return pd.DataFrame()
    g = all_pressure.groupby("miRNA")
    out = pd.DataFrame({
        "n_sets_present": g.size(),
        "n_sets_net_activating": g["arch_signed_pressure"].apply(lambda s: int((s < 0).sum())),
        "n_sets_net_repressing": g["arch_signed_pressure"].apply(lambda s: int((s > 0).sum())),
        "mean_arch_signed": g["arch_signed_pressure"].mean(),
        "min_arch_signed": g["arch_signed_pressure"].min(),
        "mean_arch_pressure": g["arch_pressure"].mean(),
        "mean_rank_shift": g["rank_shift"].mean(),
    }).reset_index()
    # ---- MALIGNANCY-role aggregation (gene-role based, prior-INDEPENDENT, ALL sets) ---- #
    # captures pro-tumor effect via repressing TSGs even where the gene-gene topology is blind
    # (e.g. miR-21 -> PTEN/PDCD4). Summed over every set the arm touches.
    if "mal_pro_tumor" in all_pressure.columns:
        gm = pd.DataFrame({
            "sum_mal_pro_tumor": g["mal_pro_tumor"].sum(),
            "mean_mal_pro_tumor": g["mal_pro_tumor"].mean(),
            "sum_mal_pro_tumor_acquired": g["mal_pro_tumor_acquired"].sum()
            if "mal_pro_tumor_acquired" in all_pressure.columns else np.nan,
            "sum_mal_pro_tumor_arch": g["mal_pro_tumor_arch"].sum()
            if "mal_pro_tumor_arch" in all_pressure.columns else np.nan,
            "sum_mal_pro_tumor_cont": g["mal_pro_tumor_cont"].sum()
            if "mal_pro_tumor_cont" in all_pressure.columns else np.nan,
            "sum_mal_pro_tumor_cont_acquired": g["mal_pro_tumor_cont_acquired"].sum()
            if "mal_pro_tumor_cont_acquired" in all_pressure.columns else np.nan,
            "sum_mal_pro_tumor_hier": g["mal_pro_tumor_hier"].sum()
            if "mal_pro_tumor_hier" in all_pressure.columns else np.nan,
            "sum_mal_pro_tumor_hier_acquired": g["mal_pro_tumor_hier_acquired"].sum()
            if "mal_pro_tumor_hier_acquired" in all_pressure.columns else np.nan,
            "sum_tf_pressure": g["tf_pressure"].sum() if "tf_pressure" in all_pressure.columns else np.nan,
            "n_tf_targets": g["n_tf_targets"].max() if "n_tf_targets" in all_pressure.columns else np.nan,
            "sum_tsg_credit": g["tsg_credit"].sum() if "tsg_credit" in all_pressure.columns else np.nan,
            "sum_onco_collateral": g["onco_collateral"].sum()
            if "onco_collateral" in all_pressure.columns else np.nan,
            "n_tsg_targets": g["n_tsg_targets"].max(),
            "n_onco_targets": g["n_onco_targets"].max(),
        }).reset_index()
        out = out.merge(gm, on="miRNA", how="left")
    # set where the arm is most net-activating
    idx = all_pressure.sort_values("arch_signed_pressure").drop_duplicates("miRNA", keep="first")
    out = out.merge(idx[["miRNA", "hallmark", "arch_signed_pressure"]]
                    .rename(columns={"hallmark": "top_activating_set",
                                     "arch_signed_pressure": "top_activating_signed"}),
                    on="miRNA", how="left")
    # ---- prior-weighted pro-tumor aggregation (only polarity-assigned sets) ---- #
    if "pro_tumor_pressure" in all_pressure.columns:
        pol = all_pressure[all_pressure["tumor_direction_prior"] != 0]
        if not pol.empty:
            gp = pol.groupby("miRNA")
            pt = pd.DataFrame({
                "n_polarity_sets": gp.size(),
                "n_sets_pro_tumor": gp["pro_tumor_pressure"].apply(lambda s: int((s > 0).sum())),
                "n_sets_anti_tumor": gp["pro_tumor_pressure"].apply(lambda s: int((s < 0).sum())),
                "sum_pro_tumor_pressure": gp["pro_tumor_pressure"].sum(),
                "mean_pro_tumor_pressure": gp["pro_tumor_pressure"].mean(),
            }).reset_index()
            tp = pol.sort_values("pro_tumor_pressure", ascending=False).drop_duplicates("miRNA", keep="first")
            pt = pt.merge(tp[["miRNA", "hallmark", "pro_tumor_pressure"]]
                          .rename(columns={"hallmark": "top_pro_tumor_set",
                                           "pro_tumor_pressure": "top_pro_tumor_pressure"}),
                          on="miRNA", how="left")
            # ACQUIRED-weighted parallel roll-up
            if "pro_tumor_acquired" in pol.columns:
                ga = pol.groupby("miRNA")
                acq = pd.DataFrame({
                    "n_sets_pro_tumor_acquired": ga["pro_tumor_acquired"].apply(lambda s: int((s > 0).sum())),
                    "sum_pro_tumor_acquired": ga["pro_tumor_acquired"].sum(),
                    "mean_pro_tumor_acquired": ga["pro_tumor_acquired"].mean(),
                }).reset_index()
                ta = pol.sort_values("pro_tumor_acquired", ascending=False).drop_duplicates("miRNA", keep="first")
                acq = acq.merge(ta[["miRNA", "hallmark"]].rename(columns={"hallmark": "top_pro_tumor_acquired_set"}),
                                on="miRNA", how="left")
                pt = pt.merge(acq, on="miRNA", how="left")
            out = out.merge(pt, on="miRNA", how="left")
    if "sum_pro_tumor_pressure" in out.columns:
        return out.sort_values("sum_pro_tumor_pressure", ascending=False,
                               na_position="last").reset_index(drop=True)
    return out.sort_values(["n_sets_net_activating", "mean_arch_signed"],
                           ascending=[False, True]).reset_index(drop=True)


def run(*, out_dir: Path = OUT_DIR, state_dir: Path = STATE_DIR,
        sets: Sequence[str] = PILOT_SETS, force_fetch: bool = False,
        all_hallmarks: bool = False, write_per_set: Optional[bool] = None) -> Dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    hs = HallmarkSets.load()
    edge = pd.read_csv(state_dir / "mirna_gene_edge_class.tsv", sep="\t")

    if all_hallmarks:
        sets = sorted(hs.sets.keys())
    # per-set files would be 3xN; for the full collection default to combined-only
    if write_per_set is None:
        write_per_set = not all_hallmarks

    results: Dict[str, pd.DataFrame] = {}
    manifest_sets = {}
    all_pressure_rows: List[pd.DataFrame] = []
    all_weight_rows: List[pd.DataFrame] = []
    set_summary_rows: List[dict] = []
    for n, sname in enumerate(sets, 1):
        members = set(hs.sets.get(sname, []))
        if not members:
            print(f"[architecture] WARN set {sname} not found; skipping")
            continue
        print(f"[architecture] ({n}/{len(sets)}) {sname}: {len(members)} genes")
        cache = CACHE_DIR / f"omnipath_{sname}.tsv"
        net = _fetch_omnipath(sorted(members), cache=cache, force=force_fetch)
        edges = _induced_edges(net, members)
        weights = _gene_weights(sorted(members), edges)
        prior = _prior(sname)
        pressure = _mirna_pressure(edge, members, weights, prior=prior)

        tag = sname.replace("HALLMARK_", "")
        if write_per_set:
            weights.to_csv(out_dir / f"architecture_{tag}_gene_weights.tsv", sep="\t", index=False)
            edges.to_csv(out_dir / f"architecture_{tag}_induced_edges.tsv", sep="\t", index=False)
            if not pressure.empty:
                pressure.to_csv(out_dir / f"architecture_{tag}_mirna_pressure.tsv", sep="\t", index=False)
        results[sname] = pressure
        if not pressure.empty:
            all_pressure_rows.append(pressure.assign(hallmark=tag))
        all_weight_rows.append(weights.assign(hallmark=tag))

        n_inh = int((weights["effect_sign"] < 0).sum())
        manifest_sets[sname] = {"n_members": len(members), "n_induced_edges": int(len(edges)),
                                "n_inhibitor_nodes": n_inh,
                                "n_redundancy_groups": int(weights["redundancy_group"].nunique())}
        top_mr = weights.sort_values("centrality_revpagerank", ascending=False).head(3)["gene"].tolist()
        top_act = (pressure.sort_values("arch_signed_pressure").head(1)
                   if not pressure.empty else pd.DataFrame())
        top_pt = (pressure.sort_values("pro_tumor_pressure", ascending=False).head(1)
                  if (not pressure.empty and prior != 0) else pd.DataFrame())
        # TOTAL miRNA-layer effect on the program (sum over ALL arms): the net program-output
        # change the whole miRNA layer imposes. <0 = program net DAMAGED (repressed), >0 = net
        # activated. For AMBIGUOUS (prior=0) sets this is the only directional readout we keep
        # (no malignancy sign), so it answers 'where/which way does pressure land inside it'.
        tot_poc = float(pressure["program_output_change"].sum()) if not pressure.empty else np.nan
        tot_pt = float(pressure["pro_tumor_pressure"].sum()) if (not pressure.empty and prior != 0) else np.nan
        tot_pta = float(pressure["pro_tumor_acquired"].sum()) if (not pressure.empty and prior != 0) else np.nan
        # GENE-ROLE malignancy aggregate (the prior-FREE per-set total): sums the per-arm
        # mal_pro_tumor = Σ_g −malignancy_sign·c_tumor over the whole set. Unlike total_pro_tumor
        # (which needs a hand-assigned set polarity and is `ambig` for prior=0 sets), this is defined
        # for ALL 50 sets from the curated TSG/oncogene roles. >0 = net pro-tumor (TSG repression
        # dominates oncogene collateral) regardless of the coarse prior.
        def _sum(col):
            return float(pressure[col].sum()) if (not pressure.empty and col in pressure.columns) else np.nan
        tot_mal = _sum("mal_pro_tumor")
        set_summary_rows.append({
            "hallmark": tag, "n_members": len(members), "n_induced_edges": int(len(edges)),
            "tumor_direction_prior": prior,
            "n_inhibitor_nodes": n_inh, "master_regulators": ";".join(top_mr),
            "top_net_activating_mirna": top_act["miRNA"].iloc[0] if not top_act.empty else None,
            "top_net_activating_signed": float(top_act["arch_signed_pressure"].iloc[0]) if not top_act.empty else np.nan,
            "top_pro_tumor_mirna": top_pt["miRNA"].iloc[0] if not top_pt.empty else None,
            "top_pro_tumor_pressure": float(top_pt["pro_tumor_pressure"].iloc[0]) if not top_pt.empty else np.nan,
            "total_program_output_change": tot_poc,
            "net_program_direction": ("damaged" if (pd.notna(tot_poc) and tot_poc < 0)
                                      else ("activated" if pd.notna(tot_poc) else None)),
            "total_pro_tumor_pressure": tot_pt,
            "total_pro_tumor_acquired": tot_pta,
            "total_mal_pro_tumor": tot_mal,
            "total_mal_pro_tumor_acquired": _sum("mal_pro_tumor_acquired"),
            "total_mal_pro_tumor_hier": _sum("mal_pro_tumor_hier"),
            "total_mal_pro_tumor_hier_acquired": _sum("mal_pro_tumor_hier_acquired"),
            "total_program_output_change_hier": _sum("program_output_change_hier"),
            "total_tf_pressure": _sum("tf_pressure"),
            "total_tsg_credit": _sum("tsg_credit"),
            "total_onco_collateral": _sum("onco_collateral"),
            "malignancy_direction": ("pro_tumor" if (pd.notna(tot_mal) and tot_mal > 0)
                                     else ("anti_tumor" if pd.notna(tot_mal) else None)),
            "total_mal_pro_tumor_cont": _sum("mal_pro_tumor_cont"),
            "total_mal_pro_tumor_cont_acquired": _sum("mal_pro_tumor_cont_acquired"),
        })
        if write_per_set:
            _print_set(sname, weights, pressure)

    # combined long-form tables + cross-set miRNA summary
    set_summary = pd.DataFrame(set_summary_rows)
    all_pressure = pd.concat(all_pressure_rows, ignore_index=True) if all_pressure_rows else pd.DataFrame()
    all_weights = pd.concat(all_weight_rows, ignore_index=True) if all_weight_rows else pd.DataFrame()
    cross = _cross_set_summary(all_pressure)
    # ---- Q6 brake-release annotation: join each arm's OWN cross-state trajectory so the
    # anti-tumor (suppressive) leaders can be read together with whether the brake itself is
    # being LOST in tumor (dHT/log2fc/pressure_move from Stage-1 state class). ------------- #
    if not cross.empty:
        sc_path = state_dir / "mirna_state_class.tsv"
        if sc_path.is_file():
            sc = pd.read_csv(sc_path, sep="\t")
            keep = [c for c in ["miRNA", "dHT", "log2fc_tumor_vs_healthy", "pressure_move",
                                "acquired_axis", "acquired_gainer"] if c in sc.columns]
            cross = cross.merge(
                sc[keep].rename(columns={"dHT": "own_dHT", "log2fc_tumor_vs_healthy": "own_log2fc",
                                         "pressure_move": "own_pressure_move",
                                         "acquired_axis": "own_acquired_axis",
                                         "acquired_gainer": "own_acquired_gainer"}),
                on="miRNA", how="left")
    # ---- Q7 family-level + total-layer rollup --------------------------------- #
    fam_long, fam_cross = _family_cross_set(all_pressure)
    if all_hallmarks or not write_per_set:
        if not set_summary.empty:
            set_summary.sort_values("n_inhibitor_nodes", ascending=False).to_csv(
                out_dir / "architecture_all_set_summary.tsv", sep="\t", index=False)
        if not all_pressure.empty:
            all_pressure.to_csv(out_dir / "architecture_all_mirna_pressure.tsv", sep="\t", index=False)
        if not all_weights.empty:
            all_weights.to_csv(out_dir / "architecture_all_gene_weights.tsv", sep="\t", index=False)
        if not cross.empty:
            cross.to_csv(out_dir / "architecture_mirna_cross_set.tsv", sep="\t", index=False)
        if not fam_long.empty:
            fam_long.to_csv(out_dir / "architecture_family_set_pressure.tsv", sep="\t", index=False)
        if not fam_cross.empty:
            fam_cross.to_csv(out_dir / "architecture_family_cross_set.tsv", sep="\t", index=False)
        _print_cross(set_summary, cross)
        _print_family(fam_cross)

    manifest = {
        "module": "mirna_hallmark.geneset_architecture",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "network_source": "OmniPath REST interactions (directed+signed)",
        "state_class_dir": str(state_dir),
        "redundancy_jaccard": REDUNDANCY_JACCARD, "pagerank_d": PAGERANK_D,
        "propagation": {"alpha": PROP_ALPHA, "horizon": PROP_HORIZON,
                        "method": "finite-horizon signed-path (Katz) influence"},
        "tumor_direction_prior": {k.replace("HALLMARK_", ""): v for k, v in TUMOR_DIRECTION_PRIOR.items()},
        "n_prior_sets": len(TUMOR_DIRECTION_PRIOR),
        "all_hallmarks": all_hallmarks, "n_sets": len(manifest_sets),
        "sets": manifest_sets,
    }
    (out_dir / "method_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    results["_set_summary"] = set_summary
    results["_cross_set"] = cross
    results["_family_cross_set"] = fam_cross
    results["_family_set_pressure"] = fam_long
    return results


def _print_family(fam_cross: pd.DataFrame) -> None:
    if fam_cross is None or fam_cross.empty or "sum_pro_tumor_pressure" not in fam_cross.columns:
        return
    print("\n[architecture] FAMILY-level net pro-tumor pressure (seed-approx; sum over family arms):")
    for _, r in fam_cross.head(8).iterrows():
        print(f"  {r.family} ({int(r.n_arms_total)} arms): pro-tumor in {int(r.n_sets_pro_tumor)}/"
              f"{int(r.n_polarity_sets)} sets, sum {r.sum_pro_tumor_pressure:+.2f} (top {r.top_pro_tumor_set})")
    print("[architecture] most ANTI-TUMOR families:")
    for _, r in fam_cross.sort_values("sum_pro_tumor_pressure").head(5).iterrows():
        print(f"  {r.family}: anti-tumor in {int(r.n_sets_anti_tumor)}/{int(r.n_polarity_sets)} "
              f"sets, sum {r.sum_pro_tumor_pressure:+.2f}")


def _print_cross(set_summary: pd.DataFrame, cross: pd.DataFrame) -> None:
    if not set_summary.empty:
        print(f"\n[architecture] {len(set_summary)} sets; richest internal topology (most inhibitor nodes):")
        for _, r in set_summary.sort_values("n_inhibitor_nodes", ascending=False).head(6).iterrows():
            print(f"  {r.hallmark}: {r.n_inhibitor_nodes} inhibitor nodes / {r.n_induced_edges} edges; "
                  f"top net-activating {r.top_net_activating_mirna} ({r.top_net_activating_signed:+.2f})")
    if not cross.empty and "sum_pro_tumor_pressure" in cross.columns:
        print("\n[architecture] miRNAs with highest NET PRO-TUMOR pressure "
              "(topology x tumor-direction prior, polarity-assigned sets only):")
        for _, r in cross.dropna(subset=["sum_pro_tumor_pressure"]).head(8).iterrows():
            print(f"  {r.miRNA}: pro-tumor in {int(r.n_sets_pro_tumor)}/{int(r.n_polarity_sets)} "
                  f"polarity sets, sum {r.sum_pro_tumor_pressure:+.2f} (top {r.top_pro_tumor_set})")
        print("\n[architecture] miRNAs with most ANTI-TUMOR (suppressive) net pressure:")
        for _, r in cross.dropna(subset=["sum_pro_tumor_pressure"]).sort_values(
                "sum_pro_tumor_pressure").head(5).iterrows():
            print(f"  {r.miRNA}: anti-tumor in {int(r.n_sets_anti_tumor)}/{int(r.n_polarity_sets)} "
                  f"polarity sets, sum {r.sum_pro_tumor_pressure:+.2f}")
        if "sum_pro_tumor_acquired" in cross.columns:
            print("\n[architecture] highest net pro-tumor pressure that is TUMOR-ACQUIRED "
                  "(parallel: weighted by c_tumor - c_gtex, not total abundance):")
            for _, r in cross.dropna(subset=["sum_pro_tumor_acquired"]).sort_values(
                    "sum_pro_tumor_acquired", ascending=False).head(8).iterrows():
                print(f"  {r.miRNA}: acquired pro-tumor sum {r.sum_pro_tumor_acquired:+.2f} "
                      f"(top {r.top_pro_tumor_acquired_set})")
    elif not cross.empty:
        print("\n[architecture] miRNAs most consistently net pathway-ACTIVATING across Hallmark sets:")
        for _, r in cross.head(8).iterrows():
            print(f"  {r.miRNA}: net-activating in {int(r.n_sets_net_activating)}/{int(r.n_sets_present)} sets, "
                  f"mean signed {r.mean_arch_signed:+.3f} (top {r.top_activating_set})")


def _print_set(sname: str, weights: pd.DataFrame, pressure: pd.DataFrame) -> None:
    print(f"\n[architecture] {sname}: {len(weights)} genes, "
          f"{int((weights['effect_sign'] < 0).sum())} inhibitor nodes")
    top_c = weights.sort_values("centrality_revpagerank", ascending=False).head(6)
    print("  master regulators (reverse-PageRank):",
          "; ".join(f"{r.gene}({r.centrality_revpagerank:.3f},sign{r.effect_sign:+d})" for _, r in top_c.iterrows()))
    if not pressure.empty:
        promoted = pressure.sort_values("rank_shift", ascending=False).head(5)
        print("  miRNAs PROMOTED by architecture (flat_rank->arch_rank):",
              "; ".join(f"{r.miRNA}({r.flat_rank}->{r.arch_rank})" for _, r in promoted.iterrows()))
        neg = pressure.sort_values("arch_signed_pressure").head(4)
        print("  net pathway-ACTIVATING pressure (represses inhibitors, signed<0):",
              "; ".join(f"{r.miRNA}({r.arch_signed_pressure:+.2f})" for _, r in neg.iterrows()))


def set_gene_hits(*, out_dir: Path = OUT_DIR, state_dir: Path = STATE_DIR) -> pd.DataFrame:
    """Per-(set, gene) 'what is hit' rollup, derived from the already-written gene weights + edge
    table (no OmniPath re-fetch). For each gene in each Hallmark set it reports the incoming tumour
    miRNA pressure that lands on it, the program-output damage that routes through it (pressure ×
    `output_influence` = the topology hub view), and the malignancy-weighted hit (pressure on a TSG =
    pro-tumor, the gene-role view). This localises a set's miRNA pressure to specific genes / hubs and
    cross-checks the two directional lenses (topology vs gene-role) at gene resolution."""
    gw = pd.read_csv(out_dir / "architecture_all_gene_weights.tsv", sep="\t")
    edge = pd.read_csv(state_dir / "mirna_gene_edge_class.tsv", sep="\t")
    ct = pd.to_numeric(edge["c_tumor"], errors="coerce")
    inc = ct.groupby(edge["gene"]).sum().rename("incoming_pressure_tumor")
    top_arm = (edge.assign(_c=ct).dropna(subset=["_c"]).sort_values("_c", ascending=False)
               .groupby("gene").first()["miRNA"].rename("top_arm"))
    g = (gw.merge(inc, left_on="gene", right_index=True, how="left")
           .merge(top_arm, left_on="gene", right_index=True, how="left"))
    g["incoming_pressure_tumor"] = g["incoming_pressure_tumor"].fillna(0.0)
    g["prog_change_gene"] = g["incoming_pressure_tumor"] * pd.to_numeric(g["output_influence"], errors="coerce")
    g["mal_hit"] = -pd.to_numeric(g["malignancy_sign"], errors="coerce") * g["incoming_pressure_tumor"]
    g["is_hub"] = False
    for _, sub in g.groupby("hallmark"):
        thr_oi = sub["output_influence"].quantile(0.9)
        thr_wa = sub["w_arch"].quantile(0.9)
        g.loc[sub.index, "is_hub"] = ((sub["output_influence"] >= thr_oi) | (sub["w_arch"] >= thr_wa))
    if "is_tf" not in g.columns:
        g["is_tf"] = False
    # a curated TF that is ALSO a topology hub = high-confidence master regulator (identity ∧ topology)
    g["hub_tf"] = g["is_hub"].astype(bool) & g["is_tf"].astype(bool)
    cols = ["hallmark", "gene", "incoming_pressure_tumor", "output_influence", "w_arch",
            "prog_change_gene", "malignancy_sign", "mal_hit", "is_hub", "is_tf", "hub_tf",
            "top_arm", "effect_sign_soft"]
    out = g[cols].sort_values(["hallmark", "prog_change_gene"], ascending=[True, False]).reset_index(drop=True)
    out.to_csv(out_dir / "architecture_set_gene_hits.tsv", sep="\t", index=False)
    print(f"[architecture] set_gene_hits: {len(out)} (set,gene) rows over "
          f"{out['hallmark'].nunique()} sets → architecture_set_gene_hits.tsv")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state-dir", type=Path, default=STATE_DIR)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--sets", nargs="*", default=PILOT_SETS)
    ap.add_argument("--all-hallmarks", action="store_true",
                    help="run every MSigDB Hallmark set (combined long-form outputs + cross-set summary)")
    ap.add_argument("--write-per-set", action="store_true",
                    help="also write per-set files when running --all-hallmarks")
    ap.add_argument("--force-fetch", action="store_true", help="re-download the OmniPath network")
    ap.add_argument("--hits-only", action="store_true",
                    help="only (re)build the per-(set,gene) hit rollup from existing outputs (no re-fetch)")
    args = ap.parse_args()
    if args.hits_only:
        set_gene_hits(out_dir=args.out_dir, state_dir=args.state_dir)
        return
    run(out_dir=args.out_dir, state_dir=args.state_dir, sets=args.sets, force_fetch=args.force_fetch,
        all_hallmarks=args.all_hallmarks, write_per_set=True if args.write_per_set else None)
    set_gene_hits(out_dir=args.out_dir, state_dir=args.state_dir)


if __name__ == "__main__":
    main()
