"""Phase 0 of the miRNA-locus GENOMIC-CONTEXT axis: classify every HE arm's precursor locus by its relation to
the host gene it sits in — STRAND-AWARE, because only a SAME-STRAND (sense) intronic miRNA is co-transcribed from
the host pre-mRNA (⇒ its abundance carries the host's transcriptional program; host mRNA is a valid co-transcription
proxy). An ANTISENSE overlap shares position but NOT transcription — it has its own/antisense promoter, so it is
autonomous-like. Classes (host-TYPE primary; `in_exon` is a secondary flag, over-calls):
  sense_coding_host   -- co-transcribed from a protein-coding host  (the MCM7/miR-93 case; AL host-proxy available + confound)
  sense_lncRNA_host   -- co-transcribed from a lncRNA host (often the pri-miRNA's OWN unit: MIR17HG, DLEU2 -> proxy, weak confound)
  antisense_overlap   -- overlaps a gene body on the OPPOSITE strand -> NOT co-transcribed (autonomous-like)
  intergenic          -- no gene-body overlap -> autonomous (own promoter)
plus flags: `clustered` (within CLUSTER_KB of another miRNA = polycistronic), `mirgenedb` (MirGeneDB-validated
bona-fide, quality layer), `n_loci`/`context_mixed` (multi-locus arms — every copy classified, representative reported).

COORDINATE AUTHORITY = GENCODE `gene_type==miRNA` (100% HE coverage, single provenance; see DATA_SOURCES).
CLI: `.venv/bin/python3 -m mirna_hallmark.learned.genomic_context`
"""
from __future__ import annotations
import os, re, numpy as np, pandas as pd
from mirna_hallmark import config as C
from mirna_hallmark.learned import data as LD, instrument as IN
from mirna_hallmark.learned.evidence import ledger as LG

GENCODE = "data/gencode.v49.annotation.gtf.csv"
OUT = "mirna_hallmark/output/learned/genomic_context.tsv"
CLUSTER_KB = 10.0
_CACHE: dict = {}

def _precursor_loci() -> pd.DataFrame:
    if "loci" not in _CACHE:
        d = pd.read_csv(C.PATHS.mirna_precursor_loci_csv,
                        usecols=["chrom", "start", "end", "strand", "gene_id", "mature_accessions"])
        d = d.dropna(subset=["chrom", "start", "end"])
        _CACHE["loci"] = d.set_index("gene_id")
        mm = {}                                                          # MIMAT -> MI precursor (coverage fallback)
        for mi, acc in zip(d["gene_id"], d["mature_accessions"].fillna("")):
            for a in str(acc).split(","):
                if a:
                    mm.setdefault(a.strip(), mi)
        _CACHE["mimat2mi"] = mm
    return _CACHE["loci"]

def _arm_mimat() -> dict:
    if "arm2mimat" not in _CACHE:
        el = IN._entity_long()
        arm = el[el["entity_level"] == "arm"][["mirbase_mature_id", "mature_accession"]].dropna().drop_duplicates()
        _CACHE["arm2mimat"] = dict(zip(arm["mirbase_mature_id"], arm["mature_accession"]))
    return _CACHE["arm2mimat"]

def _norm_name(arm: str) -> str:
    """hsa-miR-93-5p -> MIR93 ; hsa-let-7a-5p -> MIRLET7A ; hsa-miR-1185-1-3p -> MIR1185-1 (keep internal copy hyphen;
    strip ONLY the hsa- prefix + the -5p/-3p arm suffix). Do NOT strip all hyphens — MIR1185-1 must not become MIR11851."""
    b = re.sub(r"-[35]p\*?$", "", arm.replace("hsa-", ""))
    return b.replace("miR-", "MIR").replace("let-", "MIRLET").upper()

def _gencode_mir():
    """GENCODE `gene_type==miRNA` loci — THE authoritative + EXHAUSTIVE coordinate source for this axis (1879 loci,
    strand, same coordinate system as the host genes). Returns {upper name -> [(chrom,start,end,strand), ...]} keeping
    ALL loci per name (a name can recur), plus a df for clustering."""
    if "gcmir" not in _CACHE:
        genes, _ = _gencode()
        rows = [(str(r.gene_name).upper(), c, int(r.start), int(r.end), r.strand)
                for c, gg in genes.items() for r in gg[gg["gene_type"] == "miRNA"].itertuples()]
        nm: dict = {}
        for name, c, s, e, st in rows:
            nm.setdefault(name, []).append((c, s, e, st))
        _CACHE["gcmir"] = nm
        _CACHE["gcmir_names"] = list(nm)
        _CACHE["gcmir_df"] = pd.DataFrame(rows, columns=["name", "chrom", "start", "end", "strand"])
    return _CACHE["gcmir"]

def _arm_gencode_loci(arm: str):
    """ALL GENCODE loci for an arm — handles multi-COPY families GENCODE numbers/letters (MIRLET7A1/2/3, MIR548A1/2,
    MIR517A/B/C). Base = normalized name; match exact, then a copy suffix. Numeric bases (MIR93) take a family LETTER
    (never a digit — so MIR93 never falsely grabs MIR935); letter bases (MIRLET7A) take a copy DIGIT."""
    gc = _gencode_mir(); base = _norm_name(arm); names = _CACHE["gcmir_names"]; b = re.escape(base)
    if base in gc:                                                        # exact hit (single locus, incl. MIR1185-1)
        return gc[base], "gencode"
    m = re.fullmatch(r"(.+)-([0-9]+)", base)                             # copy-specific arm MIR125B-1 -> GENCODE MIR125B1
    if m and (m.group(1) + m.group(2)) in gc:
        return gc[m.group(1) + m.group(2)], "gencode"
    hits = [n for n in names if re.fullmatch(b + r"-[0-9]+", n)]          # hyphenated copies: MIR1 -> MIR1-1/MIR1-2
    if not hits and base[-1:].isalpha():                                  # concatenated copies: MIRLET7A -> MIRLET7A1/2/3
        hits = [n for n in names if re.fullmatch(b + r"[0-9]+", n)]
    if not hits and base[-1:].isdigit():                                  # family letter (numeric base): MIR517 -> MIR517A(-1)
        hits = [n for n in names if re.fullmatch(b + r"[A-Z](-?[0-9]+)?", n)]
    loci = [t for n in hits for t in gc[n]]
    return (loci, "gencode") if loci else (None, "unmapped")

def _mirgenedb_validated() -> set:
    """arms whose mature accession IS in the MirGeneDB CN-precursor set (cnv_miRNA.csv) — a QUALITY/confidence flag
    (bona-fide miRNA per MirGeneDB curation), NOT a coordinate source."""
    if "mgdb" not in _CACHE:
        loci = _precursor_loci()                                         # populates _CACHE['mimat2mi']
        arm2mimat = _arm_mimat(); mimat_set = set(_CACHE["mimat2mi"])
        _CACHE["mgdb"] = {a for a, mm in arm2mimat.items() if mm in mimat_set}
    return _CACHE["mgdb"]

def _gencode():
    if "genes" not in _CACHE:
        g = pd.read_csv(GENCODE, usecols=["seqname", "feature", "start", "end", "strand", "gene_type", "gene_name"])
        genes = g[g["feature"] == "gene"].rename(columns={"seqname": "chrom"})
        exons = g[g["feature"] == "exon"].rename(columns={"seqname": "chrom"})[["chrom", "start", "end", "gene_name"]]
        _CACHE["genes"] = {c: gg.sort_values("start") for c, gg in genes.groupby("chrom")}
        _CACHE["exons"] = {c: ee.sort_values("start") for c, ee in exons.groupby("chrom")}
    return _CACHE["genes"], _CACHE["exons"]

_CODING = {"protein_coding"}
_LNC = {"lncRNA", "lincRNA", "antisense", "processed_transcript", "sense_intronic", "TEC"}

def _classify_locus(chrom, s, e, strand):
    """Host-TYPE is primary (coding vs lncRNA) since that governs the co-transcription/AL-proxy logic; intron-vs-exon
    is a secondary `in_exon` flag (over-calls, because ANY minor-isoform exon counts — not a canonical-transcript call)."""
    genes, exons = _gencode()
    gg = genes.get(chrom)
    base = dict(host=None, host_type=None, strand_rel=None, in_exon=False)
    if gg is None:
        return {**base, "mir_class": "intergenic"}
    ov = gg[(gg["start"] <= e) & (gg["end"] >= s)]                        # gene bodies overlapping the precursor
    ov = ov[ov["gene_type"] != "miRNA"]                                   # ignore the miRNA's own gene entry
    if ov.empty:
        return {**base, "mir_class": "intergenic"}
    sense = ov[ov["strand"] == strand]
    if sense.empty:                                                       # only antisense overlap -> autonomous-like
        h = ov.iloc[(ov["end"] - ov["start"]).values.argmin()]           # smallest antisense host
        return dict(mir_class="antisense_overlap", host=h["gene_name"], host_type=h["gene_type"],
                    strand_rel="antisense", in_exon=False)
    coding = sense[sense["gene_type"].isin(_CODING)]
    pick = coding if not coding.empty else sense                         # prefer a protein-coding sense host
    h = pick.iloc[(pick["end"] - pick["start"]).values.argmin()]         # smallest containing sense host
    ee = exons.get(chrom)
    in_exon = bool(ee is not None and not ee[(ee["gene_name"] == h["gene_name"]) &
                                             (ee["start"] <= e) & (ee["end"] >= s)].empty)
    cls = "sense_coding_host" if h["gene_type"] in _CODING else "sense_lncRNA_host"
    return dict(mir_class=cls, host=h["gene_name"], host_type=h["gene_type"], strand_rel="sense", in_exon=in_exon)

def _is_clustered(chrom, s, e, strand) -> bool:
    """within CLUSTER_KB of ANOTHER GENCODE miRNA locus on the same chrom/strand (polycistronic)."""
    _gencode_mir(); allp = _CACHE["gcmir_df"]
    near = allp[(allp["chrom"] == chrom) & (allp["strand"] == strand) &
                (allp["start"] <= e + CLUSTER_KB * 1e3) & (allp["end"] >= s - CLUSTER_KB * 1e3)]
    return len(near) > 1                                                  # >1 because the arm's own locus is counted

# representative-class priority when a multi-locus arm's loci DISAGREE: coding-host first, so a host-confound
# locus is never hidden behind a majority of autonomous copies (the abundance is the SUM over loci).
_PRIORITY = ["sense_coding_host", "sense_lncRNA_host", "antisense_overlap", "intergenic"]

LOCUS_CTX_OUT = "mirna_hallmark/output/learned/locus_context.tsv"


def locus_host_map(persist: bool = False) -> dict:
    """PER-LOCUS host classification keyed on the CN instrument's `locus_id` (`MI*` precursor accession = `locus_cn`
    columns / `arm_loci_map` locus_id). This is the per-locus RESOLUTION of `classify_he_arms`'s single representative
    `host`: each hairpin locus is classified at ITS OWN hg38 coordinate (bridge = `_precursor_loci` gene_id `MI*` →
    (chrom,start,end,strand), which are exactly the CN ids — verified 506/506 — then `_classify_locus`). So the CN
    exclusion can condition each ACTIVE locus on the host it ACTUALLY sits in, instead of the arm's coding-first
    representative — which for a multi-locus arm can belong to a DIFFERENT, possibly SILENT, locus (fixes miR-16@SMC4-vs-
    DLEU2, miR-26a@CTDSP2-vs-CTDSPL, and miR-194/BPNT1 & miR-30c/NFYC whose coding host is at a silent locus). Returns
    {locus_id → dict(mir_class, host, host_type, strand_rel, in_exon, chrom, start, end, strand)} for every CN locus."""
    if "locus_host" in _CACHE:
        return _CACHE["locus_host"]
    pl = _precursor_loci()
    m: dict = {}
    for lid in [l for l in IN.locus_cn().columns if l in pl.index]:
        r = pl.loc[lid]
        r = r.iloc[0] if getattr(r, "ndim", 1) > 1 else r                     # dup gene_id guard
        try:
            c, s, e, st = str(r["chrom"]), int(r["start"]), int(r["end"]), r["strand"]
            m[lid] = {**_classify_locus(c, s, e, st), "chrom": c, "start": s, "end": e, "strand": st}
        except Exception:
            continue
    _CACHE["locus_host"] = m
    if persist:
        df = pd.DataFrame([{"locus_id": k, **v} for k, v in m.items()])[
            ["locus_id", "chrom", "start", "end", "strand", "mir_class", "host", "host_type", "strand_rel", "in_exon"]]
        os.makedirs(os.path.dirname(LOCUS_CTX_OUT), exist_ok=True)
        df.to_csv(LOCUS_CTX_OUT, sep="\t", index=False)
        print(f"wrote {LOCUS_CTX_OUT} ({len(df)} CN loci)  class dist: {df['mir_class'].value_counts().to_dict()}")
    return m

_FANTOM5 = "data/external_cache/fantom5_human_miRNA_promoters.tsv"


def _fantom5_promoters() -> dict:
    """mature arm name → list of FANTOM5 CAGE promoters as (chrom, tss, promoter_string) (de Rie 2017). Key = `miRNA`
    with the trailing `*` (old star-arm nomenclature) stripped; multi-locus arms accumulate a promoter per locus."""
    if "fantom5" not in _CACHE:
        import os
        m: dict = {}
        if os.path.exists(_FANTOM5):
            f = pd.read_csv(_FANTOM5, sep="\t")
            f["arm_norm"] = f["miRNA"].astype(str).str.rstrip("*")
            for a, g in f.groupby("arm_norm"):
                m[a] = [(str(c), t, str(p)) for c, t, p in
                        zip(g["chromosome"], pd.to_numeric(g["tss"], errors="coerce"), g["promoter"])
                        if pd.notna(p)]
        _CACHE["fantom5"] = m
    return _CACHE["fantom5"]


_CHAIN = "data/external_cache/hg19ToHg38.over.chain.gz"


def _gene_body() -> dict:
    """gene_name → (chrom, lo, hi) GENCODE (hg38) gene body span, for the accession-label within-host-body test."""
    if "gene_body" not in _CACHE:
        genes, _ = _gencode()
        sp: dict = {}
        for chrom, gg in genes.items():
            for name, s, e in zip(gg["gene_name"], gg["start"], gg["end"]):
                lo, hi = (s, e) if s <= e else (e, s)
                if name in sp and sp[name][0] == chrom:
                    sp[name] = (chrom, min(sp[name][1], lo), max(sp[name][2], hi))
                else:
                    sp.setdefault(name, (chrom, lo, hi))
        _CACHE["gene_body"] = sp
    return _CACHE["gene_body"]


def _transcript_gene() -> dict:
    """Ensembl transcript_id (versionless) → gene_name (GENCODE v49), to resolve `pN@ENST…` promoter labels to a gene."""
    if "tx_gene" not in _CACHE:
        g = pd.read_csv(GENCODE, usecols=["transcript_id", "gene_name"])
        g = g[g["transcript_id"].notna()]
        _CACHE["tx_gene"] = {str(t).split(".")[0]: n for t, n in zip(g["transcript_id"], g["gene_name"])}
    return _CACHE["tx_gene"]


def _liftover():
    """Cached hg19→hg38 LiftOver (FANTOM5/de Rie is hg19; GENCODE v49 is hg38). None if pyliftover/chain unavailable."""
    if "liftover" not in _CACHE:
        import os
        lo = None
        if os.path.exists(_CHAIN):
            try:
                from pyliftover import LiftOver
                lo = LiftOver(_CHAIN)
            except Exception:
                lo = None
        _CACHE["liftover"] = lo
    return _CACHE["liftover"]


def _prom_gene(p: str):
    """The GENE named by a FANTOM5 promoter string, or None for a coordinate-only (`p@chr…`) = own promoter."""
    if "@" not in p:
        return None
    tail = p.split("@", 1)[1]
    if tail.startswith("chr"):                                           # coordinate promoter ⇒ the miRNA's OWN
        return None
    return tail.split(",")[0].split(":")[0]


def _feature_number(chrom, pos, gene, strand):
    """Which numbered feature of `gene` a position falls in: ('exon'|'intron', N) numbered 5′→3′ (TRANSCRIPTION order,
    strand-aware), or (None, None) if outside the gene body. Uses the UNION of the gene's GENCODE exons (merged), so N
    is the collapsed-transcript feature index (not a single-isoform call — same convention as `in_exon`)."""
    _, exons = _gencode()
    ee = exons.get(chrom)
    if ee is None:
        return (None, None)
    g = ee[ee["gene_name"] == gene]
    if g.empty:
        return (None, None)
    iv = sorted({(int(s), int(e)) for s, e in zip(g["start"], g["end"])})
    merged = []
    for s, e in iv:                                                      # merge overlapping/adjacent exons (genomic order)
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    n = len(merged)
    for i, (s, e) in enumerate(merged):
        if s <= pos <= e:
            return ("exon", i + 1 if strand == "+" else n - i)          # transcription-order exon number
    for i in range(n - 1):
        if merged[i][1] < pos < merged[i + 1][0]:
            return ("intron", i + 1 if strand == "+" else (n - 1) - i)   # transcription-order intron number
    return (None, None)


def _host_strand(gene):
    """gene_name → genomic strand from GENCODE (for exon/intron numbering)."""
    if "host_strand" not in _CACHE:
        genes, _ = _gencode()
        m = {}
        for chrom, gg in genes.items():
            for name, st in zip(gg["gene_name"], gg["strand"]):
                m.setdefault(name, st)
        _CACHE["host_strand"] = m
    return _CACHE["host_strand"].get(gene)


def _locus_at(chrom, pos, prefer=None):
    """INDEPENDENT of any de Rie label — the GENCODE (hg38) annotation at a genomic position: (gene_name, gene_type,
    feature) with feature ∈ {exon, intron, intergenic}. Among the non-miRNA gene bodies containing `pos`, returns the
    `prefer` gene if it is one of them (nested/overlapping genes — e.g. MCM7 & AP4M1 share the MCM7-promoter region),
    else the SMALLEST containing gene."""
    genes, exons = _gencode()
    gg = genes.get(chrom)
    if gg is None:
        return (None, None, "intergenic")
    ov = gg[(gg["start"] <= pos) & (gg["end"] >= pos) & (gg["gene_type"] != "miRNA")]
    if ov.empty:
        return (None, None, "intergenic")
    h = ov[ov["gene_name"] == prefer].iloc[0] if prefer is not None and (ov["gene_name"] == prefer).any() \
        else ov.iloc[(ov["end"] - ov["start"]).values.argmin()]
    feat, num = _feature_number(chrom, pos, h["gene_name"], h["strand"])     # NUMBERED feature (5′→3′)
    return (h["gene_name"], h["gene_type"], f"{feat} {num}" if feat else "exon")


def _host_tss() -> dict:
    """gene_name → (chrom, tss) EXACT hg38 TSS from GENCODE (5′ start; `start` on +, `end` on −). No liftover — the host
    side of the co-transcription comparison is native hg38."""
    if "host_tss" not in _CACHE:
        genes, _ = _gencode()
        m: dict = {}
        for chrom, gg in genes.items():
            for name, s, e, st in zip(gg["gene_name"], gg["start"], gg["end"], gg["strand"]):
                m.setdefault(name, (chrom, s if st == "+" else e))
        _CACHE["host_tss"] = m
    return _CACHE["host_tss"]


def _promoter_coord(rows, host):
    """COORDINATE annotation of the arm's promoter, INDEPENDENT of the de Rie label. Lift the FANTOM5 hg19 TSS → hg38
    (`_liftover`; 0% unmapped/multi-map on this data) and (1) read the GENCODE gene / biotype / exon-vs-intron there
    (`_locus_at`, host-preferring for nested genes); (2) measure its distance to the HOST's EXACT hg38 TSS (`_host_tss`,
    no liftover) — the direct 'does the miRNA fire from the host's OWN promoter' signal (small ⇒ co-transcribed from the
    host promoter; large ⇒ an alternative/independent promoter). Returns (gene, gene_type, feature, coord_class,
    host_tss_kb). coord_class ∈ {host_shared, shared_other, intergenic, unknown}. This is POSITIONAL and can disagree
    with the LABEL/transcriptional `promoter_class` (de Rie's CAGE call) — an intronic-OWN promoter is positionally in
    the host body but fires independently (miR-1205: own by label, inside PVT1 by coordinate)."""
    lo = _liftover()
    if not rows or lo is None:
        return (None, None, None, "unknown", None)
    ht = _host_tss().get(host) if host else None
    fallback = None; hit = None; dist = None
    for c, t, p in rows:
        if t != t:
            continue
        conv = lo.convert_coordinate(c, int(t))
        if not conv:
            continue
        cc, pos = conv[0][0], conv[0][1]
        if ht and cc == ht[0]:                                           # distance to the host's exact hg38 TSS
            d = abs(pos - ht[1]); dist = d if dist is None else min(dist, d)
        g, gt, feat = _locus_at(cc, pos, prefer=host)
        if g and host and g == host and hit is None:
            hit = (g, gt, feat)
        if fallback is None:
            fallback = (g, gt, feat)
    dkb = round(dist / 1000.0, 1) if dist is not None else None
    if hit is not None:
        return (*hit, "host_shared", dkb)
    if fallback is None:
        return (None, None, None, "unknown", dkb)
    return (*fallback, "shared_other" if fallback[0] else "intergenic", dkb)


def _classify_promoter(rows, host, pad: int = 2000):
    """FANTOM5 promoter → (chosen promoter string, `promoter_class`). de Rie's LABEL is the discriminator (a gene/transcript
    label ⇒ the promoter drives THAT transcript; a coordinate label ⇒ an independent CAGE peak). ACCESSION RESOLVER (fixes
    the under-calls): a `pN@ENST…` label is resolved to its gene via GENCODE; a non-symbol transcript accession (ENST-not-in-
    v49 / GenBank `AK…`) is confirmed as the host when its FANTOM5 hg19 TSS, LIFTED to hg38, lands INSIDE the host body — that
    within-body test is applied ONLY to accession labels (they name a host transcript, e.g. miR-21→p1@AK310806=VMP1), NOT to
    coordinate labels (an intronic-OWN peak is also within the body, so it stays `own`). Classes:
      host_shared  — label names the host (symbol or ENST→host) OR an accession-label TSS lifts inside the host body;
      shared_other — label names a DIFFERENT gene / transcript;
      own          — coordinate-only or the miRNA's OWN `MIR*`/`LET*` gene promoter (independent);
      unknown      — precursor absent from the FANTOM5 CAGE atlas.
    Position (GENCODE class) and transcription (this) can DISAGREE — an intronic arm may be `own`/`shared_other`."""
    if not rows:
        return (None, "unknown")
    tx = _transcript_gene()
    resolved = []                                                        # (chrom, tss, promoter, gene, is_accession)
    for c, t, p in rows:
        g = _prom_gene(p); acc = False
        if g and g.startswith("ENST"):
            g2 = tx.get(g.split(".")[0])
            if g2:
                g = g2                                                   # ENST resolved to its gene
            else:
                acc = True                                               # ENST not in v49 → treat as an accession
        elif g and re.match(r"^[A-Z]{2}\d", g):                          # GenBank accession (AK…, BC…, …)
            acc = True
        resolved.append((c, t, p, g, acc))
    for c, t, p, g, acc in resolved:                                     # (1) label names the host (symbol or ENST→host)
        if g and host and g == host:
            return (p, "host_shared")
    body = _gene_body().get(host) if host else None; lo = _liftover()    # (2) ACCESSION label → within host body (lifted)?
    if body and lo is not None:
        for c, t, p, g, acc in resolved:
            if acc and t == t:
                conv = lo.convert_coordinate(c, int(t))
                if conv and conv[0][0] == body[0] and body[1] - pad <= conv[0][1] <= body[2] + pad:
                    return (p, "host_shared")
    named = [(p, g) for c, t, p, g, acc in resolved if g and not g.upper().startswith(("MIR", "LET"))]
    if named:                                                            # (3) names a DIFFERENT gene / transcript
        return (named[0][0], "shared_other")
    return (rows[0][2], "own")                                           # (4) coordinate-only / own MIR* promoter


def classify_he_arms():
    """One row per HE arm. MULTI-LOCUS arms (let-7a = 3 GENCODE loci) are handled explicitly: every locus is
    classified; `n_loci` = GENCODE copy count; `context_mixed` = the loci DISAGREE on class; `mir_class` = the
    representative (majority, ties broken toward `sense_coding_host` since the abundance is the sum and a single
    host-coupled locus already injects the host program); `classes` = the full per-locus breakdown. Dose-weighted
    host-coupled fraction (using per-locus expression / w_norm) is the Phase-1 refinement, noted not built here."""
    from collections import Counter
    ed = LG.pooled_he_edges(); X = LD._load()["X"]
    he_arms = sorted(set(ed["miRNA"]) & set(X.index))
    mgdb = _mirgenedb_validated(); fant = _fantom5_promoters()
    rows = []
    for arm in he_arms:
        loci, src = _arm_gencode_loci(arm)
        base = dict(arm=arm, mirgenedb=arm in mgdb, coord_source=src)
        if loci is None:
            fr = fant.get(arm, [])
            prom, pc = _classify_promoter(fr, None)
            pg, pgt, pfeat, pcc, pdkb = _promoter_coord(fr, None)
            rows.append({**base, "n_loci": 0, "mir_class": "unmapped", "host": None, "host_type": None,
                         "strand_rel": None, "in_exon": False, "host_region": None, "context_mixed": False, "classes": "",
                         "clustered": False, "promoter": prom, "promoter_class": pc, "promoter_gene": pg,
                         "promoter_gene_type": pgt, "promoter_feature": pfeat, "promoter_coord_class": pcc,
                         "promoter_host_tss_kb": pdkb})
            continue
        cls = [_classify_locus(c, s, e, st) for (c, s, e, st) in loci]
        cnt = Counter(x["mir_class"] for x in cls)
        top = max(cnt, key=lambda k: (cnt[k], -_PRIORITY.index(k) if k in _PRIORITY else -99))
        rep_idx = next(i for i, x in enumerate(cls) if x["mir_class"] == top)
        rep = cls[rep_idx]; rc, rs, re_, rst = loci[rep_idx]             # representative locus coordinates
        host_region = None                                              # which NUMBERED intron/exon of the host the miRNA sits in
        if rep["host"]:
            feat, num = _feature_number(rc, (rs + re_) // 2, rep["host"], _host_strand(rep["host"]) or rst)
            if feat:
                host_region = f"{feat} {num}"
        fr = fant.get(arm, [])
        prom, pc = _classify_promoter(fr, rep["host"])                   # LABEL/resolver = TRANSCRIPTIONAL co-transcription call
        pg, pgt, pfeat, pcc, pdkb = _promoter_coord(fr, rep["host"])     # COORDINATE = independent positional call + host-TSS distance
        rows.append({**base, "n_loci": len(loci), "mir_class": top, "host": rep["host"], "host_type": rep["host_type"],
                     "strand_rel": rep["strand_rel"], "in_exon": rep["in_exon"], "host_region": host_region,
                     "context_mixed": len(cnt) > 1, "classes": ",".join(f"{k}:{v}" for k, v in cnt.items()),
                     "clustered": any(_is_clustered(*t) for t in loci), "promoter": prom, "promoter_class": pc,
                     "promoter_gene": pg, "promoter_gene_type": pgt, "promoter_feature": pfeat,
                     "promoter_coord_class": pcc, "promoter_host_tss_kb": pdkb})
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(OUT), exist_ok=True); df.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT} ({len(df)} HE arms)  coord source: {df['coord_source'].value_counts().to_dict()}")
    print(f"multi-locus arms (n_loci>1): {int((df['n_loci'] > 1).sum())}   of which context-MIXED: {int(df['context_mixed'].sum())}")
    print("\n=== genomic-context distribution (representative class) ===")
    print(df["mir_class"].value_counts().to_string())
    print(f"\nclustered (polycistronic): {int(df['clustered'].sum())}/{len(df)}   MirGeneDB-validated: {int(df['mirgenedb'].sum())}")
    hc = df[df["mir_class"] == "sense_coding_host"]
    print(f"\n=== sense_coding_host (positionally intronic in a protein-coding host): {len(hc)} arms ===")
    print(hc["host"].value_counts().head(10).to_string())
    print(f"\n=== promoter_class (FANTOM5 CAGE co-transcription call) ===")
    print(df["promoter_class"].value_counts().to_string())
    print(f"of sense_coding_host: {hc['promoter_class'].value_counts().to_dict()} "
          f"(host_shared={100 * (hc['promoter_class'] == 'host_shared').mean():.0f}% co-transcribed; "
          f"own/shared_other = positionally-intronic but transcriptionally INDEPENDENT)")
    return df

def _resid1(v, c):
    ok = np.isfinite(v) & np.isfinite(c)
    if ok.sum() < 30:
        return None, ok
    A = np.column_stack([np.ones(ok.sum()), c[ok]]); b, *_ = np.linalg.lstsq(A, v[ok], rcond=None)
    r = np.full(len(v), np.nan); r[ok] = v[ok] - A @ b
    return r, ok

def _residN(v, Z):
    """Residualize v on a covariate MATRIX Z (n×k), NaN-safe. Returns (resid, ok-mask)."""
    Z = np.atleast_2d(np.asarray(Z, float))
    if Z.shape[0] != len(v):
        Z = Z.T
    ok = np.isfinite(v) & np.isfinite(Z).all(1)
    if ok.sum() < 30:
        return None, ok
    A = np.column_stack([np.ones(ok.sum()), Z[ok]]); b, *_ = np.linalg.lstsq(A, v[ok], rcond=None)
    r = np.full(len(v), np.nan); r[ok] = v[ok] - A @ b
    return r, ok

def host_target_relatedness(min_r: float = 0.3):
    """#4 — the CN-instrument host pleiotropy (and the coupling host-confound) is a threat ONLY where the host gene
    is actually RELATED to the target. For every sense_coding_host edge (arm intronic in host H, target T), score
    partial Spearman(H, T | CPE, **miRNA**) — host-target co-expression beyond purity AND beyond the arm itself. Partialling
    the miRNA (upgrade 2026-07-11, user-driven) isolates the host's OWN target-association (host protein → T) from the host
    merely PROXYING the miR's transcription (sense co-transcription, ρ≈0.29) — since the miR is the thing we study, removing
    it sharpens the gate for BOTH consumers (CN exclusion pleiotropy + coupling confound). Empirically near-identical to
    CPE-only (mean |ρ| 0.39→0.36; 16/19 down-weighted segments survive) but correctly DEMOTES the co-transcription-entangled
    MCM7 tail (E2F5/BTG2@MCM7 drop <0.3). Proliferation kept IN (it mediates for prolif hosts). Columns: `r_host_target`
    (CPE-only, reference) + `r_host_target_mir` (CPE+miR); **`related` uses the miRNA-partialled `r_host_target_mir`**. Arms
    absent from X (unexpressed) fall back to CPE-only. `related` = |r|>=min_r ⇒ a real host→T path ⇒ genuine pleiotropy/confound."""
    from scipy.stats import spearmanr
    gc = pd.read_csv(OUT, sep="\t")
    rep_host = dict(gc[gc["mir_class"] == "sense_coding_host"][["arm", "host"]].dropna().itertuples(index=False))
    lhm = locus_host_map()                                             # MI* → per-locus host (the more-correct enumeration)
    # PER-LOCUS coding hosts (fix, 2026-07-11): a multi-locus arm has a DISTINCT coding host per locus (miR-26a: CTDSP2 +
    # CTDSPL; miR-218: SLIT2 + SLIT3; miR-103a: PANK2 + PANK3), and the coding host of miR-194/miR-30c sits on a SILENT
    # locus — so scoring only the coding-first REPRESENTATIVE both misses the second host and can mis-locate. Enumerate
    # EVERY distinct coding host across an arm's CN loci; keep `host_locus` for traceability; representative = fallback.
    arm_hosts: dict = {}                                              # arm → [(host, host_locus), ...]
    for arm in rep_host:
        lm = IN.arm_loci_map().get(arm)
        seen: dict = {}
        if lm is not None:
            for lid in lm["locus_id"].astype(str):
                h = lhm.get(lid)
                if h and h["mir_class"] == "sense_coding_host" and isinstance(h["host"], str):
                    seen.setdefault(h["host"], lid)                   # first locus per distinct coding host
        if not seen and isinstance(rep_host.get(arm), str):
            seen[rep_host[arm]] = ""                                  # arm without CN-locus coverage → representative fallback
        if seen:
            arm_hosts[arm] = list(seen.items())
    ed = LG.pooled_he_edges(); Y = LD._load()["Y"]; clin = LD._load()["C"]; X = LD._load()["X"]
    cpe = pd.to_numeric(clin.get("CPE"), errors="coerce") if "CPE" in clin.columns else None
    nmap = LD._arm_name_map(X) if hasattr(LD, "_arm_name_map") else {}
    sub = ed[ed["miRNA"].isin(arm_hosts)][["miRNA", "gene"]].drop_duplicates()
    rows = []
    for arm, gene in sub.itertuples(index=False):
        for host, host_locus in arm_hosts[arm]:
            if host not in Y.index or gene not in Y.index or host == gene:
                continue
            xname = nmap.get(arm) or (arm if arm in X.index else None)
            parts = Y.columns if cpe is None else Y.columns.intersection(cpe.index)
            if xname is not None:
                parts = parts.intersection(X.columns)                # common set so r_cpe & r_mir are comparable
            h = Y.loc[host, parts].to_numpy(float); t = Y.loc[gene, parts].to_numpy(float)
            if cpe is None:
                r_cpe = r_mir = spearmanr(h, t).correlation
            else:
                c = cpe.reindex(parts).to_numpy(float)
                hr, ok = _resid1(h, c); tr, _ = _resid1(t, c)
                if hr is None:
                    continue
                r_cpe = spearmanr(hr[ok], tr[ok]).correlation
                if xname is not None:                                # partial the miRNA too (arm expression)
                    mir = X.loc[xname, parts].to_numpy(float)
                    hr2, ok2 = _residN(h, np.column_stack([c, mir])); tr2, _ = _residN(t, np.column_stack([c, mir]))
                    r_mir = spearmanr(hr2[ok2], tr2[ok2]).correlation if hr2 is not None else r_cpe
                else:
                    r_mir = r_cpe                                     # arm unexpressed → can't partial it, fall back
            rows.append(dict(arm=arm, host=host, host_locus=host_locus, target=gene, r_host_target=round(float(r_cpe), 3),
                             r_host_target_mir=round(float(r_mir), 3), related=abs(r_mir) >= min_r))
    df = pd.DataFrame(rows)
    out = "mirna_hallmark/output/learned/host_target_relatedness.tsv"
    df.to_csv(out, sep="\t", index=False)
    n_rel = int(df["related"].sum())
    print(f"wrote {out}: {len(df)} coding-host edges | RELATED (|r_mir|>={min_r}): {n_rel} ({100*n_rel//max(1,len(df))}%) "
          f"= genuine pleiotropy/confound; the rest = benign shared locus")
    print(f"mean |r| CPE-only {df.r_host_target.abs().mean():.3f} -> +miR {df.r_host_target_mir.abs().mean():.3f} "
          f"(drop = miR-shared share); related under CPE-only would be {int((df.r_host_target.abs()>=min_r).sum())}")
    print("per-arm: fraction of an arm's targets that are host-related (top pleiotropy-risk arms):")
    ar = df.groupby(["arm", "host"]).agg(n=("target", "size"), rel_frac=("related", "mean")).reset_index()
    print(ar[ar["n"] >= 3].sort_values("rel_frac", ascending=False).head(10).to_string(index=False))
    return df

HOST_REL_OUT = "mirna_hallmark/output/learned/host_target_relatedness.tsv"
HOST_LENS_OUT = "mirna_hallmark/output/learned/host_confound_lens.tsv"

def _host_edge_verdict(gene, host, arm, folds, n_iter, Y):
    """Generalized host-gene confounder verdict for ONE coding-host edge — MH-100's OOF 2×2, but conditioning on the
    SPECIFIC host gene H (not the broad proliferation metagene). Does controlling H de-confound the arm's coupling
    to `gene` (confound) or over-control it (mechanism)? Reuses `prolif_verdict`'s OOF/classify/edge helpers."""
    from mirna_hallmark.learned import prolif_verdict as PV, spike_slab as SS, families as FAM
    from scipy.stats import spearmanr
    Yg, X, C, w = LD.assemble_gene(gene, w_prior_source="ledger")
    fam = FAM.family_of(X.columns); afam = fam.get(arm)
    Xf, wf, _ = FAM.collapse_by_family(X, w, fam)
    if host not in Y.index or afam not in Xf.columns or host == gene:
        return None
    parts = Yg.index
    Hr = PV._resid(Y.loc[host, parts].to_numpy(float), C.to_numpy(float))       # host axis orthogonal to C
    Caug = C.copy(); Caug["host_axis"] = Hr
    yv = Yg.to_numpy(float); Cc, Ca = C.to_numpy(float), Caug.to_numpy(float)
    pC = PV._oof_pred(Yg, Xf, C, wf, folds, n_iter)
    pCR = PV._oof_pred(Yg, Xf, Caug, wf, folds, n_iter)
    # REPRESSION-DIRECTED coupling = −ρ (repression ⇒ pred anti-correlates with target ⇒ ρ<0 ⇒ −ρ>0). NOT abs():
    # abs conflates real repression-strengthening (−0.35→−0.47) with sign-FLIPS to anti-repression (−0.04→+0.16).
    def rho(pred, S): return -spearmanr(PV._resid(pred, S), PV._resid(yv, S)).correlation
    dC = rho(pCR, Cc) - rho(pC, Cc); dCR = rho(pCR, Ca) - rho(pC, Ca)           # >0 = controlling H STRENGTHENS repression
    gclass = PV._classify(dC, dCR)
    Mb, _ = SS.fit_gene_ss(Yg, Xf, C, wf, gene=gene, n_iter=n_iter, burn=n_iter // 3)
    Ma, _ = SS.fit_gene_ss(Yg, Xf, Caug, wf, gene=gene, n_iter=n_iter, burn=n_iter // 3)
    b0, b1 = float(Mb.get(afam, 0.0)), float(Ma.get(afam, 0.0))
    rel = (b1 - b0) / (b0 + 1e-6)
    return dict(arm=arm, host=host, target=gene, family=afam, beta_base=round(b0, 4), beta_host=round(b1, 4),
                rel=round(rel, 3), dC=round(float(dC), 3), dCR=round(float(dCR), 3),
                gene_verdict=gclass, edge_verdict=PV._edge_verdict(gclass, rel))

_HL = {"folds": 4, "n_iter": 700}
def _hl_worker(t):
    Y = LD._load()["Y"]
    try:
        return _host_edge_verdict(t[0], t[1], t[2], _HL["folds"], _HL["n_iter"], Y)
    except Exception as e:
        return ("err", t, str(e))

def host_confound_lens(*, workers=8, folds=4, n_iter=700, limit=None):
    """The MH-99/100 UNIFICATION: run the host-gene OOF 2×2 over every RELATED coding-host edge (host_target_relatedness),
    classifying per-edge {confound (control the host) / over_control (host is the mechanism) / fragile / neutral}."""
    import time
    from concurrent.futures import ProcessPoolExecutor
    _HL["folds"], _HL["n_iter"] = folds, n_iter
    rel = pd.read_csv(HOST_REL_OUT, sep="\t")
    rel = rel[rel["related"]][["target", "host", "arm"]]
    if limit:
        rel = rel.head(limit)
    tasks = [(g, h, a) for g, h, a in rel.itertuples(index=False)]
    LD._load()                                                                 # warm parent caches pre-fork
    print(f"host-gene lens: {len(tasks)} related coding-host edges · {workers} workers", flush=True)
    rows, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, r in enumerate(ex.map(_hl_worker, tasks, chunksize=4)):
            if isinstance(r, dict):
                rows.append(r)
            if i % 100 == 0:
                print(f"  [{i}/{len(tasks)}] {(time.time()-t0)/60:.1f}min", flush=True)
    df = pd.DataFrame(rows)
    # NOVEL-host flag: is the confounding host a proliferation gene (already caught by MH-100) or a NOVEL non-prolif
    # confounder? host_prolif_corr = |Spearman(host, broad cell-cycle metagene)|; host_novel = <0.4 (MH-100-invisible).
    from scipy.stats import spearmanr
    from mirna_hallmark.learned import prolif_verdict as PV
    Y = LD._load()["Y"]; pg = [g for g in PV._PROLIF if g in Y.index]
    P = Y.loc[pg]; P = P.sub(P.mean(1), axis=0).div(P.std(1) + 1e-9, axis=0).mean(0).reindex(Y.columns)
    hpc = {h: abs(spearmanr(Y.loc[h], P).correlation) for h in df["host"].unique() if h in Y.index}
    df["host_prolif_corr"] = df["host"].map(hpc).round(2)
    df["host_novel"] = df["host_prolif_corr"] < 0.4
    df.to_csv(HOST_LENS_OUT, sep="\t", index=False)
    print(f"\nwrote {HOST_LENS_OUT} ({len(df)} edges) in {(time.time()-t0)/60:.1f}min")
    print("edge verdicts:", df["edge_verdict"].value_counts().to_dict())
    print("gene verdicts:", df["gene_verdict"].value_counts().to_dict())
    conf = df[df["gene_verdict"] == "confound"]
    print(f"confound: {len(conf)} | NOVEL non-prolif hosts (|host_prolif_corr|<0.4): {int(conf['host_novel'].sum())}")
    print(f"mean dC={df.dC.mean():+.3f} (>0 = host de-confounds; <0 = over-control/mechanism)")
    return df

FAM_CTX_OUT = "mirna_hallmark/output/learned/family_context.tsv"


def family_context():
    """MH-101 (c) — per seed-family CONTEXT MIX. A seed family pools arms of possibly-DIFFERENT genomic contexts (some
    coding-host, some intergenic/antisense; some co-transcribed, some independent), and the coupling unit is the pooled
    `X_fam`, so the family's host-program confound / CN-instrument validity is a DOSE-WEIGHTED MIX of its members, not a
    single class. This crosses the per-arm annotation (`genomic_context.tsv`) with the family structure + member abundance.
    Per family: member context composition; **dose-weighted coding-host & co-transcribed fractions** (Σ RPM-share over the
    sense_coding_host / host_shared members — the fraction of the family dose that carries a host-program confound); the
    dominant (dose-max) arm + its context; the distinct coding hosts (single-host ⇒ uniform confound; multi ⇒ mixed).
    Generalises MH-100's arm-source finding (miR-93-5p carries the miR-17~92 family's confound) to every family."""
    from mirna_hallmark.learned import families as FAM
    X = LD._load()["X"]                                                                  # arms × participants (arms = index)
    gc = pd.read_csv(OUT, sep="\t").set_index("arm")
    fam = FAM.family_of(list(X.index))
    rpm = pd.Series((np.power(2.0, X.to_numpy(float)) - 1.0).mean(1), index=X.index)      # mean RPM per arm (over participants) = dose proxy
    rows = []
    for f, arms in fam.groupby(fam).groups.items():
        arms = [a for a in arms if a in gc.index]                                        # HE-annotated members only
        if not arms:
            continue
        d = rpm.reindex(arms).fillna(0.0); share = d / (d.sum() or 1.0)
        cls, prom, hosts = gc.loc[arms, "mir_class"], gc.loc[arms, "promoter_class"], gc.loc[arms, "host"]
        comp = cls.value_counts().to_dict()
        is_ch = (cls == "sense_coding_host")
        dom = share.idxmax()
        ch_hosts = sorted({h for h in hosts[is_ch] if isinstance(h, str)})
        rows.append(dict(
            family=f, n_members=len(arms), composition=",".join(f"{k}:{v}" for k, v in comp.items()),
            context_homogeneous=(len(comp) == 1),
            dose_coding_host_frac=round(float(share[is_ch].sum()), 3),
            dose_cotranscribed_frac=round(float(share[is_ch & (prom == "host_shared")].sum()), 3),
            dominant_arm=dom.replace("hsa-", ""), dominant_class=cls.get(dom),
            dominant_host=hosts.get(dom) if is_ch.get(dom, False) else None,
            dominant_promoter=prom.get(dom),
            n_distinct_coding_hosts=len(ch_hosts), coding_hosts=",".join(ch_hosts[:4])))
    df = pd.DataFrame(rows).sort_values("dose_cotranscribed_frac", ascending=False).reset_index(drop=True)
    os.makedirs(os.path.dirname(FAM_CTX_OUT), exist_ok=True); df.to_csv(FAM_CTX_OUT, sep="\t", index=False)
    multi = df[df["n_members"] > 1]
    print(f"wrote {FAM_CTX_OUT} ({len(df)} families; {len(multi)} multi-arm)")
    print(f"context-MIXED multi-arm families: {int((~multi['context_homogeneous']).sum())}/{len(multi)}")
    print(f"families with a host-confound dose (dose_coding_host_frac≥0.5): {int((df['dose_coding_host_frac'] >= 0.5).sum())} "
          f"| co-transcribed≥0.5: {int((df['dose_cotranscribed_frac'] >= 0.5).sum())}")
    print("\n=== top families by dose-weighted CO-TRANSCRIBED host fraction (the host-confounded ones) ===")
    cols = ["family", "n_members", "dose_cotranscribed_frac", "dominant_arm", "dominant_host", "n_distinct_coding_hosts", "coding_hosts"]
    print(df[df["dose_cotranscribed_frac"] > 0].head(15)[cols].to_string(index=False))
    return df


if __name__ == "__main__":
    import sys
    a = sys.argv[1:]
    if a and a[0] == "--relatedness":
        host_target_relatedness()
    elif a and a[0] == "--host-lens":
        host_confound_lens()
    elif a and a[0] == "--family-context":
        family_context()
    else:
        classify_he_arms()
