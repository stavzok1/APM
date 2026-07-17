"""Within-family structure — classify what a collapsed seed family actually IS inside (§8 QC + taxonomy).

The model collapses same-seed arms into one family predictor (β_f), assuming §8 (same seed → identical targeting).
This module looks INSIDE the collapse along family-intrinsic axes and flags where the assumption breaks:

  • EXPRESSION   — per member median RPM: expressed (≥ floor) vs silent (a subset carries the family).
  • SEED-isomiR  — per member 5'-isomiR fraction from RAW TCGA/CPTAC isoform-seq (`seed_shift`): a 5' shift moves
                   the seed → the member's true seed ≠ the family's → the collapse is partly WRONG (mis-collapse).
  • LOCUS-origin — per member `n_active_source` (`_arm_focal_locus`): from an active transcriptional locus vs silent.

CLASS (family-intrinsic, gene-independent):
  homogeneous-collapsed | dose-partitioned (some silent) | seed-heterogeneous (mis-collapse) | all-silent/low.

Gene-conditional refinement (the "within-family confidence") — `family_card(gene)` folds in the Bayesian between-
family PIP-entry (`channel_cn.between_family_bayes`, MH-94) + `resolution_report` member-resolution. The graded
δ-pooling width is the future upgrade of the member-resolution axis.
"""
from __future__ import annotations
import glob, os
from collections import defaultdict
import numpy as np
import pandas as pd

ROOT = "/sci/labs/michall/stavzok/APM"
_ISO_TCGA = f"{ROOT}/data/miRNA/GDCdata_isoform"
_ISO_TEST = f"{ROOT}/data/miRNA/GDCdata_isoform_test"
_ISO_CPTAC = f"{ROOT}/data/CPTAC/cptac2_mirna_isoform"
_CACHE_SEED = f"{ROOT}/mirna_hallmark/output/learned/within_family/seed_shift.tsv"
_WIN = 6          # real 5'-isomiRs are within ±WIN nt of the dominant 5'; farther = a different locus (paralog)
FLOOR = 10.0      # RPM detectability floor (matches the arm-expression floor)


_MANIFEST = f"{ROOT}/data/miRNA/gdc_mirna_isoform_download_manifest.tsv"


def _file_meta() -> dict:
    """file UUID → (participant 12-char TCGA id, sample_type). From the GDC manifest. The isoform files are named
    `<uuid>.isoforms.quantification.txt`; the GDCdata_isoform_test files use the GDC data uuid as the dir name."""
    m = pd.read_csv(_MANIFEST, sep="\t")
    part = m["cases.submitter_id"].astype(str) if "cases.submitter_id" in m else m["cases"].astype(str).str[:12]
    return {str(u): (p[:12], str(st)) for u, p, st in zip(m["id"], part, m["sample_type"])}


def _iso_files(sample_types=("Primary Tumor",), include_cptac: bool = False) -> list[tuple]:
    """Raw isoform files as (path, participant, sample_type), FILTERED to `sample_types` (default Primary Tumor —
    the model's cohort; the manifest also carries NAT / Metastatic which the model does NOT use). TCGA-only by
    default (CPTAC is a different platform whose 5'-calling would contaminate the seed reference). Normal-like PAM50
    subtype is NOT filtered here — it is inherited downstream by matching participants to the model's `X` columns
    (which are already normal-like-excluded). `sample_types=None` ⇒ all types."""
    meta = _file_meta()
    out = []
    for f in glob.glob(f"{_ISO_TCGA}/*.isoforms.quantification.txt"):
        uuid = os.path.basename(f).split(".")[0]
        part, st = meta.get(uuid, (None, None))
        if sample_types is None or st in sample_types:
            out.append((f, part, st))
    for f in glob.glob(f"{_ISO_TEST}/**/*.isoforms.quantification.txt", recursive=True):
        uuid = os.path.basename(os.path.dirname(f))                            # test dir named by uuid
        part, st = meta.get(uuid, (None, None))
        if sample_types is None or st in sample_types:
            out.append((f, part, st))
    if include_cptac:
        out += [(f, None, "CPTAC") for f in glob.glob(f"{_ISO_CPTAC}/*.isoforms.quantification.txt")]
    seen = {}                                                                  # ONE file per participant (matches the model's
    for f, p, st in out:                                                       # multi-vial mean-collapse; CPTAC has no participant)
        seen.setdefault(p if p is not None else f, (f, p, st))
    return sorted(seen.values())


_ISO_COLS = ["miRNA_region", "isoform_coords", "read_count", "reads_per_million_miRNA_mapped"]


def _parse_mature(f: str):
    """One isoform file → ({mimat: {5'_start: read_count}}, {mimat: total_rpm}) for mature isoforms — VECTORISED
    (usecols + groupby, no per-row Python loop). Module-level so it's picklable for the parallel parse."""
    try:
        df = pd.read_csv(f, sep="\t", usecols=_ISO_COLS)
    except Exception:
        return {}, {}
    df = df[df["miRNA_region"].astype(str).str.startswith("mature")]
    if not len(df):
        return {}, {}
    cc = df["isoform_coords"].astype(str).str.split(":", expand=True)
    se = cc[2].str.split("-", expand=True)
    start = pd.to_numeric(se[0], errors="coerce"); end = pd.to_numeric(se[1], errors="coerce")
    g = pd.DataFrame({"mi": df["miRNA_region"].astype(str).str.split(",").str[1].to_numpy(),
                      "fp": np.where(cc[3].to_numpy() == "+", start.to_numpy(), end.to_numpy()),
                      "rc": pd.to_numeric(df["read_count"], errors="coerce").fillna(0.0).to_numpy(),
                      "rpm": pd.to_numeric(df["reads_per_million_miRNA_mapped"], errors="coerce").fillna(0.0).to_numpy()})
    g = g.dropna(subset=["fp"]); g["fp"] = g["fp"].astype(int)
    dist: dict = {}
    for (mi, fp), r in g.groupby(["mi", "fp"])["rc"].sum().items():
        dist.setdefault(mi, {})[int(fp)] = float(r)
    return dist, g.groupby("mi")["rpm"].sum().to_dict()


def seed_shift(force: bool = False, min_reads: int = 20, sample_types=("Primary Tumor",)) -> pd.DataFrame:
    """Per mature arm: 5'-isomiR seed-shift fraction, computed **PER SAMPLE** then averaged (NOT pooled — pooling
    conflates inter-sample 5'-dominant disagreement into false heterogeneity: miR-181b/374a read ~0.48 pooled but
    ~0.01 per-sample). Two passes: (1) pooled reads → per-arm GLOBAL-mode 5' (the canonical seed reference, robust to
    miRBase-version offset); (2) per file, fraction of that sample's within-locus (±6nt) mature reads NOT at the mode.
    Arm summary = mean ± sd of per-sample fractions over samples with ≥`min_reads`. Columns: mimat, arm, total_rpm,
    seed_shift_frac (per-sample MEAN), seed_shift_sd, n_samples, top_shift_nt, n_files. Cached."""
    if os.path.exists(_CACHE_SEED) and not force:
        return pd.read_csv(_CACHE_SEED, sep="\t")
    loci = pd.read_csv(f"{ROOT}/data/miRNA/mirna_mature_loci.csv")
    mimat2arm = dict(zip(loci["mature_accession"], loci["mirbase_mature_id"]))
    fileinfo = _iso_files(sample_types=sample_types)                        # (path, participant, sample_type), model cohort
    paths = [fi[0] for fi in fileinfo]
    from concurrent.futures import ProcessPoolExecutor                      # parse files in PARALLEL across cores
    nproc = min(8, (os.cpu_count() or 4))
    with ProcessPoolExecutor(max_workers=nproc) as ex:
        parsed = list(ex.map(_parse_mature, paths, chunksize=8))           # ~n_files/nproc parse time
    # pass 1: pooled reads → per-arm GLOBAL mode 5' (reference) + cohort RPM + pooled alt-shift direction
    pool: dict = defaultdict(lambda: defaultdict(float)); rpmtot: dict = defaultdict(float)
    for dist, rpm in parsed:
        for mi, d in dist.items():
            for fp, c in d.items():
                pool[mi][fp] += c
        for mi, p in rpm.items():
            rpmtot[mi] += p
    mode = {mi: max(d, key=d.get) for mi, d in pool.items()}
    # pass 2: per-sample fraction relative to the fixed mode
    persamp: dict = defaultdict(list)
    for dist, _ in parsed:
        for mi, d in dist.items():
            m0 = mode[mi]
            win = {fp: c for fp, c in d.items() if abs(fp - m0) <= _WIN}
            tot = sum(win.values())
            if tot >= min_reads:
                persamp[mi].append(1.0 - win.get(m0, 0.0) / tot)
    rows = []
    for mi, fr in persamp.items():
        if rpmtot[mi] < 50 or len(fr) < 3:
            continue
        fr = np.array(fr)
        m0 = mode[mi]
        alts = sorted([(int(fp - m0), c) for fp, c in pool[mi].items() if fp != m0 and abs(fp - m0) <= _WIN],
                      key=lambda x: -x[1])
        rows.append({"mimat": mi, "arm": mimat2arm.get(mi, mi), "total_rpm": round(rpmtot[mi], 1),
                     "seed_shift_frac": round(float(fr.mean()), 3), "seed_shift_sd": round(float(fr.std()), 3),
                     "n_samples": len(fr), "top_shift_nt": alts[0][0] if alts else 0})
    out = pd.DataFrame(rows)
    out["n_files"] = len(paths)
    os.makedirs(os.path.dirname(_CACHE_SEED), exist_ok=True)
    out.to_csv(_CACHE_SEED, sep="\t", index=False)
    return out


def _seed_family_maps():
    """(seed→family, mimat→mature_seq) from miR_Family_Info (human) — for resolving each 5'-offset's seed to a family."""
    fi = pd.read_csv(f"{ROOT}/data/miRNA/miR_Family_Info.txt", sep="\t")
    fi = fi[fi["MiRBase ID"].astype(str).str.startswith("hsa-")]
    seed2fam = dict(zip(fi["Seed+m8"].astype(str).str.upper().str.replace("T", "U"), fi["miR family"]))
    seq = {a: str(s).upper().replace("T", "U")
           for a, s in zip(fi["MiRBase Accession"], fi["Mature sequence"]) if isinstance(s, str)}
    return seed2fam, seq


def seed_composition(sample_types=("Primary Tumor",), shrink_k: float = 20.0) -> pd.DataFrame:
    """THE model-facing per-sample tensor (Phase 1, WITH the JUMP leg) — NOT averaged. Per (participant × arm), the
    fraction of arm m's mature reads in participant s carrying **each family's seed**: `c_{s,m,f}`. An arm's reads
    split over destination families by their 5'-offset seed — the CANONICAL family (offset 0), any JUMP family (a
    5'-isomiR whose shifted seed = another family's, e.g. miR-29b +1 → miR-767-5p), and `orphan` (shifted seed matches
    no known family → dropped from every family). This is what corrects `x(m,s)` into `X_fam,f[s]` (Phase 2:
    `X_fam,f[s] = log2(1+Σ_{m} c_{s,m,f}·RPM_{s,m})` — a family GAINS an arm's jumped reads and LOSES the reads that
    jump away). Per-sample variation kept (a constant is absorbed by z-scoring); low-read cells Dirichlet-shrunk toward
    the arm's cohort family-distribution (pseudo-count `shrink_k`). Keyed by participant → intersect with the model's
    `X` columns. Columns: participant, arm, target_family, frac, n_reads (`is_canonical` marks the offset-0 family)."""
    fileinfo = _iso_files(sample_types=sample_types)
    paths = [fi[0] for fi in fileinfo]; parts = [fi[1] for fi in fileinfo]
    loci = pd.read_csv(f"{ROOT}/data/miRNA/mirna_mature_loci.csv")
    mimat2arm = dict(zip(loci["mature_accession"], loci["mirbase_mature_id"]))
    seed2fam, seq = _seed_family_maps()
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as ex:
        parsed = list(ex.map(_parse_mature, paths, chunksize=8))
    pool: dict = defaultdict(lambda: defaultdict(float))
    for dist, _ in parsed:
        for mi, d in dist.items():
            for fp, c in d.items():
                pool[mi][fp] += c
    mode = {mi: max(d, key=d.get) for mi, d in pool.items()}

    def fam_at(mi, offset):                                                  # the family whose seed sits at this 5'-offset
        s = seq.get(mi); i = 1 + offset
        if not s or i < 0 or i + 7 > len(s):
            return "orphan"
        return seed2fam.get(s[i:i + 7], "orphan")

    # per-arm offset→family (only within the window) + cohort family distribution (Dirichlet shrink target)
    off_fam: dict = {}; coh: dict = defaultdict(lambda: defaultdict(float))
    for mi, d in pool.items():
        m0 = mode[mi]; of = {}
        for fp, c in d.items():
            if abs(fp - m0) <= _WIN:
                f = fam_at(mi, int(fp - m0)); of[fp] = f; coh[mi][f] += c
        off_fam[mi] = of
        tot = sum(coh[mi].values()) or 1.0
        coh[mi] = {f: c / tot for f, c in coh[mi].items()}
    rows = []
    for (dist, _), part in zip(parsed, parts):
        if part is None:
            continue
        for mi, d in dist.items():
            of = off_fam.get(mi)
            if not of:
                continue
            m0 = mode[mi]; fam_reads: dict = defaultdict(float)
            for fp, c in d.items():
                if fp in of:
                    fam_reads[of[fp]] += c
            n = sum(fam_reads.values())
            if n < 1:
                continue
            canon = of.get(m0, "orphan")
            union = set(coh[mi]) | set(fam_reads)
            raw = {f: fam_reads.get(f, 0.0) + shrink_k * coh[mi].get(f, 0.0) for f in union}   # obs + EB prior
            tot = sum(raw.values()) or 1.0
            for f, v in raw.items():
                frac = v / tot                                                # renormalise per cell → Σ_f frac = 1 exactly
                if frac > 1e-3:
                    rows.append((part, mimat2arm.get(mi, mi), f, round(float(frac), 4), int(n), f == canon))
    return pd.DataFrame(rows, columns=["participant", "arm", "target_family", "frac", "n_reads", "is_canonical"])


def _classify(n, n_expr, n_silent, max_shift):
    if max_shift == max_shift and max_shift > 0.30:
        return "seed-heterogeneous"                                # a member's 5'-isomiR moves its seed off the family
    if n_expr == 0:
        return "all-silent/low"                                    # no expressed member (collapse moot in BRCA)
    if n_silent > 0:
        return "dose-partitioned"                                  # expressed subset carries it, silent passengers
    return "homogeneous-collapsed"                                 # all expressed + seed-homogeneous = one clean number


def family_structure(min_members: int = 2, seed_het_thr: float = 0.30) -> pd.DataFrame:
    """Family-intrinsic within-family taxonomy over all multi-member seed families (gene-independent). One row/family:
    n members, n_expr / n_silent, max_seed_shift (+ which member), n_active_loci, and the class."""
    from mirna_hallmark.learned import data as LD, families as FAM, instrument as IN
    sh = dict(zip(seed_shift()["arm"], seed_shift()["seed_shift_frac"]))
    X = LD._load()["X"]                                            # arm × participant, log2(RPM+1)
    med = (2.0 ** X - 1.0).median(1)                               # per-arm median RPM
    fam = FAM.family_of(pd.Index(X.index))
    d = pd.DataFrame({"arm": X.index, "family": fam.reindex(X.index).values, "med_rpm": med.reindex(X.index).values})
    d = d.dropna(subset=["family"])
    d["seed_shift"] = d["arm"].map(sh)
    d["n_active"] = d["arm"].map(lambda a: IN._arm_focal_locus(a).get("n_active_source", 0))
    rows = []
    for f, g in d.groupby("family"):
        if len(g) < min_members:
            continue
        expr = g[g.med_rpm >= FLOOR]; sil = g[g.med_rpm < FLOOR]
        cov = g.dropna(subset=["seed_shift"])
        maxsh = cov["seed_shift"].max() if len(cov) else np.nan
        het_arm = cov.loc[cov["seed_shift"].idxmax(), "arm"] if len(cov) and maxsh > seed_het_thr else ""
        rows.append({"family": f, "n": len(g), "n_expr": len(expr), "n_silent": len(sil),
                     "max_seed_shift": round(maxsh, 3) if maxsh == maxsh else np.nan,
                     "seed_het_arm": het_arm.replace("hsa-", ""), "n_active_loci": int((g["n_active"] > 0).sum()),
                     "class": _classify(len(g), len(expr), len(sil), maxsh)})
    return pd.DataFrame(rows).sort_values(["class", "max_seed_shift"], ascending=[True, False])


def shifted_seed_family(seed_het_thr: float = 0.30) -> pd.DataFrame:
    """DISPOSITIVE mis-collapse test: for each seed-heterogeneous arm, does its dominant 5'-isomiR SEED (nt 2-8,
    shifted by the observed 5' offset) belong to a DIFFERENT known seed family? If yes, that fraction of the arm's
    reads targets another family's genome-wide targetome ⇒ §8 is literally wrong for them. Self-validates by
    recovering miRBase's own annotated seed variants (miR-411-5p.1↔.2 etc.). Handles ±1/+2 shifts from the mature seq.
    Columns: arm, shift_nt, seed_shift_frac, canon_seed/canon_fam, shifted_seed/shifted_fam, jumps_family."""
    fi = pd.read_csv(f"{ROOT}/data/miRNA/miR_Family_Info.txt", sep="\t")
    fi = fi[fi["MiRBase ID"].astype(str).str.startswith("hsa-")].copy()
    fi["seed"] = fi["Seed+m8"].astype(str).str.upper().str.replace("T", "U")
    seed2fam = dict(zip(fi["seed"], fi["miR family"]))
    rec = fi.dropna(subset=["Mature sequence"]).drop_duplicates("MiRBase Accession")
    mimat = {r["MiRBase Accession"]: (str(r["Mature sequence"]).upper().replace("T", "U"), r["seed"], r["miR family"])
             for _, r in rec.iterrows()}
    q = seed_shift()
    q = q[(q.seed_shift_frac > seed_het_thr) & (q.top_shift_nt.abs() <= 2)]
    rows = []
    for _, r in q.iterrows():
        m = mimat.get(r["mimat"])
        if not m:
            continue
        seq, cseed, cfam = m
        i = 1 + int(r["top_shift_nt"])
        sseed = seq[i:i + 7] if 0 <= i and i + 7 <= len(seq) else None       # seed of the 5'-shifted isoform
        if not sseed or len(sseed) < 7:
            continue
        sfam = seed2fam.get(sseed, "(no known family)")
        rows.append({"arm": r["arm"].replace("hsa-", ""), "shift_nt": int(r["top_shift_nt"]),
                     "seed_shift_frac": r["seed_shift_frac"], "canon_seed": cseed, "canon_fam": cfam,
                     "shifted_seed": sseed, "shifted_fam": sfam,
                     "jumps_family": sfam not in ("(no known family)", cfam)})
    return pd.DataFrame(rows).sort_values("seed_shift_frac", ascending=False)


def phantom_scan(floor: float = 10.0, min_jump: float = 0.05, min_donor_rpm: float = 100.0) -> pd.DataFrame:
    """Systematic PHANTOM-REGULATOR scan (generalises the miR-767 case, MH-97). A phantom = an **HE regulator arm
    that is unexpressed in BRCA** (median RPM < `floor`) **whose seed is manufactured as a 5'-isomiR by an ABUNDANT
    arm** (a jump into the phantom's family) — so its curated edges can only be that abundant arm's isomiRs here, and
    the model credits an absent regulator. Joins: (1) HE arms + median RPM + n curated targets; (2) each arm's seed
    family; (3) `seed_composition` jumps (donor arm → target_family). Reports each low-RPM HE arm whose family is a
    jump target from an expressed donor, with the donor, jump fraction, and shared curated targets."""
    from mirna_hallmark.learned import data as LD
    from mirna_hallmark.learned.evidence import ledger as LG
    he = LG.pooled_he_edges()
    tgt = he.groupby("miRNA")["gene"].apply(lambda s: set(s))                  # arm -> curated target set
    X = LD._load()["X"]; med = (2.0 ** X - 1.0).median(1)                      # per-arm median RPM
    # arm -> seed family (canonical), from miR_Family_Info
    fi = pd.read_csv(f"{ROOT}/data/miRNA/miR_Family_Info.txt", sep="\t")
    fi = fi[fi["MiRBase ID"].astype(str).str.startswith("hsa-")]
    arm2fam = dict(zip(fi["MiRBase ID"], fi["miR family"]))
    comp = seed_composition()
    jumps = (comp[(~comp["is_canonical"]) & (comp["target_family"] != "orphan")]
             .groupby(["arm", "target_family"])["frac"].mean().reset_index())  # donor arm -> target family -> mean frac
    jumps = jumps[jumps["frac"] >= min_jump]
    rows = []
    for arm in tgt.index:                                                      # HE regulator arms only
        if med.get(arm, 0.0) >= floor:                                         # expressed → not a phantom
            continue
        fam = arm2fam.get(arm)
        if fam is None:
            continue
        donors = jumps[(jumps["target_family"] == fam) & (jumps["arm"] != arm)]
        donors = donors[donors["arm"].map(lambda a: med.get(a, 0.0)) >= min_donor_rpm]  # abundant donor
        for _, dr in donors.iterrows():
            shared = tgt[arm] & tgt.get(dr["arm"], set())
            rows.append({"phantom_arm": arm.replace("hsa-", ""), "phantom_rpm": round(float(med.get(arm, 0)), 1),
                         "phantom_family": fam, "n_curated_targets": len(tgt[arm]),
                         "donor_arm": dr["arm"].replace("hsa-", ""), "donor_rpm": round(float(med.get(dr["arm"], 0)), 0),
                         "jump_frac": round(float(dr["frac"]), 3), "n_shared_targets": len(shared),
                         "shared_targets": ",".join(sorted(shared)[:6])})
    return pd.DataFrame(rows).sort_values(["n_shared_targets", "jump_frac"], ascending=False)


def delta_pooling(gene: str, family: str, use_cn: bool = True, n_iter: int = 1500, burn: int = 500, seed: int = 0):
    """δ-pooling — the GRADED within-family member resolution (`M_member = M_fam + δ_m`, Decision F). Same-seed members
    repress identically (§8), so this grades **dose-delivery** resolution — which member(s) the family's realized effect
    is attributed to — NOT repression identity. **Fuses THREE signals** into one per-member share:
      (1) mRNA — members as SEPARATE predictors (dense posterior); collinearity ⇒ wide `sd`, low resolution `z=|M/sd|`;
      (2) per-segment CN — `instrument.family_causal_attribution` (Ring-1): resolves members on INDEPENDENT loci
          (BETWEEN segments), weighted by first-stage `seg_F`;
      (3) chimeric — Ring-1's within-segment L2 (`l2_src`), resolves same-locus members the CN can't.
    **Fusion = INVERSE-VARIANCE combination of two INDEPENDENT share estimators** — the mRNA share (precision `z²`) and
    the exogenous-CN share (precision `seg_F`): `share ∝ (z²·mrna_share + seg_F·cn_share)/(z²+seg_F)`. This is the
    RIGHT model here (a collinear observational estimator + an exogenous instrument), NOT a joint *likelihood*: a
    fully-joint hierarchical Gibbs was tested and is WORSE — the mRNA (n≈1000) dominates any CN prior, and a global τ²
    lets the resolved members' mRNA↔CN disagreement loosen the prior on the *unresolved* member, burying exactly the
    member the CN should rescue (+ unphysical negative M). `use_cn=False` = mRNA-only. Returns (per-member DataFrame,
    confidence=top share). CN+chimeric sharpen members the collinear mRNA cannot (e.g. a member on its own CN segment)."""
    from mirna_hallmark.learned import data as LD, families as FAM, attribution_eb as AE, spike_slab as SS, instrument as IN
    Y, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger", deconv=True)
    fam = FAM.family_of(pd.Index(X.columns))
    members = [a for a in X.columns if str(fam.get(a)) == family]
    if len(members) < 2:
        return pd.DataFrame({"arm": members, "z": [np.nan]}), np.nan
    yr, Xz, cols = AE._prep(Y, X[members], C)                                  # MEMBER-level (not family-collapsed)
    m, s, _ = SS._gibbs_posterior(Xz, yr, np.ones(Xz.shape[1]), n_iter=n_iter, burn=burn, seed=seed)
    z = np.abs(m) / (s + 1e-9)
    mrna_share = np.clip(m, 0, None) / (np.clip(m, 0, None).sum() or 1.0)
    df = pd.DataFrame({"arm": [c.replace("hsa-", "") for c in cols], "M": m.round(4), "sd": s.round(4),
                       "z": z.round(2), "mrna_share": mrna_share.round(3)})
    df["fused_share"] = df["mrna_share"]
    if use_cn:
        try:
            r1 = IN.family_causal_attribution(gene, members)                   # per-member CN+chimeric portion + seg_F + l2_src
        except Exception:
            r1 = pd.DataFrame()
        if len(r1):
            r1 = r1.assign(arm=r1["arm"].str.replace("hsa-", "", regex=False))
            df = df.merge(r1[["arm", "portion", "seg_F", "l2_src"]], on="arm", how="left")
            prec_m = df["z"].fillna(0) ** 2                                    # mRNA precision (resolution z²)
            prec_c = df["seg_F"].fillna(0).clip(upper=200)                     # exogenous-CN precision (first-stage F)
            cn_share = df["portion"].fillna(0.0)
            fused = (prec_m * df["mrna_share"] + prec_c * cn_share) / (prec_m + prec_c).replace(0, np.nan)
            df["fused_share"] = (fused.fillna(df["mrna_share"]) / (fused.fillna(df["mrna_share"]).sum() or 1.0)).round(3)
    df = df.sort_values("fused_share", ascending=False)
    return df, float(df["fused_share"].max())


def family_card(gene: str, min_members: int = 2) -> pd.DataFrame:
    """Gene-conditional within-family card — the FULL four-axis view per family targeting `gene`:
      1. intrinsic CLASS (expression + seed-isomiR + locus);   2. LOCUS: `n_active_loci` (active-source members);
      3. gene-conditional ENTRY: the Bayesian between-family PIP (`channel_cn.between_family_bayes`, MH-94);
      4. within-family CONFIDENCE (graded): `delta_conf` = the δ-pooling max member `|M/sd|` (≥2 ⇒ a member is
         resolvable; <2 ⇒ family-only) — supersedes the discrete `resolution_report` axis (kept for reference)."""
    from mirna_hallmark.learned import families as FAM, channel_cn as CH
    struct = family_structure(min_members=min_members).set_index("family")
    res = FAM.resolution_report(gene).set_index("family") if gene else pd.DataFrame()
    try:
        bf = CH.between_family_bayes(gene)
        pip = bf.groupby("family")["pip_cn"].max() if len(bf) else pd.Series(dtype=float)
    except Exception:
        pip = pd.Series(dtype=float)
    fams = [f for f in res.index if f in struct.index] if len(res) else []
    rows = []
    for f in fams:
        r = struct.loc[f]
        try:
            _, dconf = delta_pooling(gene, f)                                  # graded within-family confidence
        except Exception:
            dconf = np.nan
        rows.append({"family": f, "class": r["class"], "n_expr": int(r["n_expr"]), "n_silent": int(r["n_silent"]),
                     "n_active_loci": int(r["n_active_loci"]), "max_seed_shift": r["max_seed_shift"],
                     "pip_entry": round(float(pip.get(f, np.nan)), 3),
                     "delta_conf": round(float(dconf), 2) if dconf == dconf else np.nan,
                     "resolution": res.loc[f, "resolution"]})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--seed":
        print(seed_shift(force=True).sort_values("seed_shift_frac", ascending=False).head(20).to_string(index=False))
    else:
        s = family_structure()
        print("=== within-family taxonomy (family-intrinsic) ===")
        print(s["class"].value_counts().to_string())
        print(s[s["class"] == "seed-heterogeneous"].head(15).to_string(index=False))
