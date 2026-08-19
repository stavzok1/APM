import pandas as pd, numpy as np, pathlib, warnings
warnings.filterwarnings("ignore")
B = pathlib.Path("mirna_hallmark/output/learned")
g = pd.read_csv(B/"realization/gene_card.tsv", sep="\t", low_memory=False)
def num(s):
    if s.dtype==object: return s.map({True:1.0,False:0.0,"True":1.0,"False":0.0}).astype(float)
    return pd.to_numeric(s, errors="coerce")
h = g[num(g["regulatory_handoff"])==1].copy()
both = g["dominant_HLY"].notna() & g["dominant_TUM"].notna()
print(f"handoff defined on {both.sum()} genes (needs HLY and TUM both present); fires on {len(h)} ({len(h)/both.sum():.1%})")
hly_leg = h[h["dominant_HLY"] != h["dominant_TUM"]]
nat_leg = h[h["dominant_NAT"] != h["dominant_TUM"]]
print(f"  HLY!=TUM leg: {len(hly_leg)}   NAT!=TUM leg: {len(nat_leg)}   both: {len(h[(h.dominant_HLY!=h.dominant_TUM)&(h.dominant_NAT!=h.dominant_TUM)])}")
print(f"  fires ONLY via the NAT leg (HLY==TUM): {len(h)-len(hly_leg)}  <-- these are NOT healthy->tumour switches")
print("\n=== GENUINE healthy -> tumour switches (HLY != TUM), n=%d ===" % len(hly_leg))
pair = hly_leg.groupby(["dominant_HLY","dominant_TUM"]).size().sort_values(ascending=False)
print(pair.head(10).to_string())
net = (hly_leg["dominant_TUM"].value_counts() - hly_leg["dominant_HLY"].value_counts().reindex(
        hly_leg["dominant_TUM"].value_counts().index).fillna(0))
allarms = sorted(set(hly_leg["dominant_HLY"])|set(hly_leg["dominant_TUM"]))
t = hly_leg["dominant_TUM"].value_counts().reindex(allarms).fillna(0)
hh = hly_leg["dominant_HLY"].value_counts().reindex(allarms).fillna(0)
net = (t-hh).sort_values()
print("\nNET dominance change across the healthy->tumour switch (tumour count - healthy count):")
print("  GAINERS:", {k: int(v) for k,v in net.tail(8).items()})
print("  LOSERS :", {k: int(v) for k,v in net.head(8).items()})
print("\nshift_class of the dominant edge, HLY!=TUM genes:")
print(hly_leg["dominant_edge_shift_class"].value_counts().head(8).to_string())
