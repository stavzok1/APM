"""Geometric helpers for CNV distance computation."""
import pandas as pd

def compute_signed_distance_with_overlap(ds, de, gs, ge, strand):
    if de >= gs and ds <= ge:
        return 0, max(ds, gs), min(de, ge)
    if de < gs:
        d = gs - de
        return (-d if strand == "+" else +d), None, None
    if ds > ge:
        d = ds - ge
        return (+d if strand == "+" else -d), None, None
    return 0, None, None

def compute_signed_distance_interval(ss, se, es, ee):
    if se >= es and ss <= ee: return 0
    if se < es: return -(es - se)
    if ss > ee: return ss - ee
    return 0
def classify_gene_hit(ds, de, gs, ge, strand, signed_dist, gf,
                      pu=2000, pd_=500, pw=5000):
    ibs, ibe = max(ds, gs), min(de, ge)
    if ibe >= ibs:
        ob = ibe - ibs + 1
        op = (ob / (ge - gs + 1)) * 100
        bd = 1
    else:
        ob = op = 0
        bd = 0
    if strand == "+":
        ps, pe = gs - pu, gs + pd_
    else:
        ps, pe = ge - pd_, ge + pu
    pm = 1 if min(de, pe) >= max(ds, ps) else 0
    he = hu = hs = hst = False
    ti = tt = None
    if gf is not None and not gf.empty:
        for _, f in gf.iterrows():
            fs, fe = int(f["start"]), int(f["end"])
            if de < fs or ds > fe: continue
            ft = f.get("feature", None)
            if ft in ("exon","CDS","start_codon","stop_codon","transcript"):
                he = True
                if ft == "stop_codon": hs = True
                if ft == "start_codon": hst = True
                elif ft == "transcript":
                    ti = f["transcript_id"]
                    tt = f["transcript_type"]
            elif ft == "UTR": hu = True
    ef = 1 if (he or hu) else 0
    io = 1 if (bd and not ef) else 0
    rh = ""
    if pm or ob > 0:
        if pm: rh += "promoter,"
        if ob > 0: rh += ("exon_or_UTR," if ef else "intron,")
    else:
        if abs(signed_dist) <= pw:
            rh += ("upstream_5kb," if signed_dist < 0 else "downstream_5kb,")
        else: rh += "intergenic"
    return {"overlap_bp": ob, "overlap_percent": op,
            "promoter_flag": pm, "gene_body_flag": bd,
            "exon_flag": ef, "intron_only_flag": io,
            "upstream_5kb_flag": int("upstream_5kb" in rh),
            "downstream_5kb_flag": int("downstream_5kb" in rh),
            "region_hit": rh,
            "stop_codon_flag": int(hs), "start_codon_flag": int(hst),
            "transcript_id": ti, "transcript_type": tt}
def classify_elem_hit(ss, se, es, ee, sd, pw=5000):
    i_s, i_e = max(ss, es), min(se, ee)
    if i_e >= i_s:
        ol = i_e - i_s + 1
        pc = (ol / (ee - es + 1)) * 100
        of = 1
    else:
        ol = pc = 0
        of = 0
    if of: rh = "overlaps"
    elif abs(sd) <= pw:
        if sd < 0: rh = "proximal_upstream"
        elif sd > 0: rh = "proximal_downstream"
        else: rh = "proximal"
    else: rh = "distal"
    return {"overlap_bp": ol, "overlap_percent": pc,
            "overlaps_flag": of, "region_hit": rh,
            "proximal_flag": int("proximal" in rh),
            "distal_flag": int("distal" in rh)}
