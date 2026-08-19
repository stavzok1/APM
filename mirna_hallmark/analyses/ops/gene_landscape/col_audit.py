"""COLUMN AUDIT — block by block: is every column meaningful, and is it documented?"""
import pandas as pd, numpy as np, pathlib, re, warnings, collections, json
warnings.filterwarnings("ignore")
B = pathlib.Path("mirna_hallmark/output/learned")
from mirna_hallmark.learned import card_glossary as GL
CARDS = {"arm": B/"arm_card.tsv", "edge": B/"realization/edge_card.tsv",
         "gene": B/"realization/gene_card.tsv", "gene_family": B/"gene_family_card.tsv",
         "seed_family": B/"seed_family_card.tsv"}
reg = pd.read_csv(B/"card_registry.tsv", sep="\t", dtype=str).fillna("")
dom = {(r.card, r.column): r.domain for r in reg.itertuples()}
def blk(c):
    m = re.match(r"^([a-z]{2,12}?)_", c); return m.group(1)+"_" if m else "(bare)"

flags = collections.defaultdict(list); rows=[]
for card, path in CARDS.items():
    df = pd.read_csv(path, sep="\t", low_memory=False)
    n = len(df)
    for c in df.columns:
        s = df[c]; nn = s.notna().sum(); frac = nn/n
        nu = s.nunique(dropna=True)
        f = []
        if nn == 0: f.append("ALL-NaN")
        elif nu == 1: f.append("CONSTANT")
        elif frac < 0.05: f.append("<5% populated")
        if frac < 0.5 and not dom.get((card,c)): f.append("sparse, NO domain")
        rows.append({"card":card,"block":blk(c),"col":c,"n":nn,"frac":frac,"nuniq":nu,
                     "desc":bool(GL.describe(c,card)),"flags":";".join(f)})
        for x in f: flags[x].append(f"{card}:{c}")
A = pd.DataFrame(rows)
print(f"columns audited: {len(A)}   described: {A.desc.sum()} ({A.desc.mean():.1%})\n")

print("="*100); print("BLOCK-BY-BLOCK AUDIT (blocks with >=4 columns, worst fill first)"); print("="*100)
gb = A.groupby(["card","block"]).agg(ncol=("col","size"), med_fill=("frac","median"),
                                     min_fill=("frac","min"), flagged=("flags", lambda s:(s!="").sum()),
                                     described=("desc","sum")).reset_index()
gb = gb[gb.ncol>=4].sort_values("med_fill")
print(gb.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

print("\n" + "="*100); print("⛔ COLUMNS THAT CARRY NO INFORMATION"); print("="*100)
for k in ["ALL-NaN","CONSTANT"]:
    v = flags.get(k, [])
    print(f"\n{k}: {len(v)}")
    for x in v: print(f"   {x}")
print(f"\n<5% populated: {len(flags.get('<5% populated',[]))}")
for x in flags.get("<5% populated", []): print(f"   {x}")

print("\n" + "="*100); print("⚠ SPARSE WITH NO DOMAIN ENTRY (a blank here is unexplained)"); print("="*100)
sp = flags.get("sparse, NO domain", [])
print(f"{len(sp)} columns:")
for x in sp[:40]: print(f"   {x}")
if len(sp)>40: print(f"   … and {len(sp)-40} more")

print("\n" + "="*100); print("DUPLICATE CONTENT — columns that are identical to another on the same card"); print("="*100)
for card, path in CARDS.items():
    df = pd.read_csv(path, sep="\t", low_memory=False)
    num = df.select_dtypes(include=[np.number])
    seen = {}; dups=[]
    for c in num.columns:
        key = tuple(pd.util.hash_pandas_object(num[c].fillna(-1.23456e9), index=False))[:0] or None
        h = pd.util.hash_pandas_object(num[c].fillna(-1.23456e9), index=False).sum()
        if h in seen:
            if num[c].equals(num[seen[h]]): dups.append((seen[h], c))
        else: seen[h] = c
    if dups:
        print(f"\n[{card}] {len(dups)} identical pairs:")
        for a,b in dups[:12]: print(f"   {a}  ==  {b}")
