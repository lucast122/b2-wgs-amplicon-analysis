#!/usr/bin/env python3
"""Osmolyte / compatible-solute gene profiling. Maps DIAMOND-vs-(curated UniProt osmolyte DB)
hits -> genes -> osmotic-stress pathways, normalised to hits per million reads (HPM), and tests
the drought contrast. Hypothesis (parallel to the plant osmolyte response in Bai et al.): does
microbial osmotic-stress gene potential increase under drought?"""
import os, glob, re
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import stats

OZ = "/mnt/disk4/timo/gbi/b2/functional/osmolyte"; HITS = f"{OZ}/hits"; MAP = f"{OZ}/id2gene.map"
DOUT = "/mnt/disk4/timo/gbi/diamond_out"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

def pathway(g):
    g = g.lower()
    if g in ("otsa","otsb","tres","trey","trez"): return "Trehalose"
    if g.startswith("ect"): return "Ectoine"
    if g in ("beta","betb","bett","gbsa","gbsb"): return "Glycine-betaine"
    if g == "ggps": return "Glucosylglycerol"
    if g == "mtld": return "Mannitol"
    if g in ("prov","prow","prox","opud"): return "Solute transport"
    return "other"

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]; m=len(p)
    if m==0: return p
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(q,0,1); return out
def cliffs(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)==0 or len(b)==0: return np.nan
    return (sum((bj>a).sum() for bj in b)-sum((bj<a).sum() for bj in b))/(len(a)*len(b))

def main():
    id2g={}
    for line in open(MAP):
        a=line.rstrip("\n").split("\t")
        if len(a)>=2: id2g[a[0]]=a[1]
    files=sorted(glob.glob(f"{HITS}/*.osm.tsv"))
    print(f"parsing {len(files)} osmolyte hit files")
    gene={}; path={}
    for fp in files:
        acc=re.search(r"(SRR\d+)",os.path.basename(fp)).group(1)
        gc=defaultdict(int); pc=defaultdict(int)
        for line in open(fp):
            p=line.split("\t")
            if len(p)<2: continue
            g=id2g.get(p[1])
            if not g: continue
            gc[g]+=1; pc[pathway(g)]+=1
        if gc: gene[acc]=pd.Series(gc); path[acc]=pd.Series(pc)
    reads={acc:(sum(1 for _ in open(f"{DOUT}/{acc}.diamond_lca.tsv")) if os.path.exists(f"{DOUT}/{acc}.diamond_lca.tsv") else 5_000_000) for acc in gene}
    gmat=pd.DataFrame(gene).fillna(0).T; pmat=pd.DataFrame(path).fillna(0).T
    if "other" in pmat: pmat=pmat.drop(columns=["other"])
    hpm_g=gmat.div(pd.Series(reads).reindex(gmat.index)/1e6,axis=0)
    hpm_p=pmat.div(pd.Series(reads).reindex(pmat.index)/1e6,axis=0)
    hpm_g.to_csv(f"{RES}/osmolyte_gene_hpm.tsv",sep="\t"); hpm_p.to_csv(f"{RES}/osmolyte_pathway_hpm.tsv",sep="\t")
    meta=pd.read_csv(f"{A}/metadata.tsv",sep="\t",index_col=0)
    s=[x for x in hpm_p.index if x in meta.index]; md=meta.loc[s]
    pre=[x for x in md.index[md.condition=="pre-drought"]]; dro=[x for x in md.index[md.condition=="drought"]]
    print(f"matrix: {hpm_p.shape[0]} samples; pathways: {list(hpm_p.columns)}")
    rows=[]
    print("\nosmotic-stress pathway potential (HPM), pre vs drought:")
    for c in sorted(hpm_p.columns):
        a=hpm_p.loc[pre,c].values; b=hpm_p.loc[dro,c].values
        try:_,p=stats.mannwhitneyu(a,b,alternative="two-sided")
        except ValueError:p=np.nan
        d=cliffs(a,b); fc=np.log2((b.mean()+1e-3)/(a.mean()+1e-3))
        arrow="UP" if b.mean()>a.mean() else "down"
        print(f"  {c:18s} pre={a.mean():6.2f} drought={b.mean():6.2f}  {arrow}  log2FC={fc:+.2f} p={p:.3f} d={d:+.2f}")
        rows.append({"pathway":c,"pre":a.mean(),"drought":b.mean(),"log2FC":fc,"p":p,"cliffs_delta":d})
    r=pd.DataFrame(rows)
    if len(r): r["q_BH"]=bh(r["p"].values); r.round(4).to_csv(f"{RES}/osmolyte_drought.tsv",sep="\t",index=False)
    tot=hpm_p.sum(1)
    _,pt=stats.mannwhitneyu(tot.loc[pre],tot.loc[dro],alternative="two-sided")
    print(f"\nTOTAL osmolyte gene potential: pre={tot.loc[pre].mean():.1f} drought={tot.loc[dro].mean():.1f} HPM, p={pt:.3f}, d={cliffs(tot.loc[pre],tot.loc[dro]):+.2f}")
    nsig=int((r['q_BH']<0.05).sum()) if len(r) else 0
    open(f"{RES}/osmolyte_checks.txt","w").write(
        f"n_samples\t{hpm_p.shape[0]}\npathways\t{hpm_p.shape[1]}\ntotal_pre\t{tot.loc[pre].mean():.2f}\ntotal_drought\t{tot.loc[dro].mean():.2f}\ntotal_p\t{pt:.4f}\nn_pathway_signif_q05\t{nsig}\n")
    print(f"\n=== {nsig} pathways significant at q<0.05 ===")

if __name__=="__main__":
    main()
