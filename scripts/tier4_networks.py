#!/usr/bin/env python3
"""
TIER 4a — microbial co-occurrence networks, pre-drought vs drought (WGS genus).
Tests a NETWORK-level version of the manuscript's convergence hypothesis: does
drought simplify / fragment the co-occurrence structure (fewer edges, lower
connectance, higher modularity)? Independent of the PCoA-dispersion test.

Method: prevalent genera (present in >=25% of samples, mean >=0.1%), CLR-transformed,
Spearman correlations within each condition; edges = |rho|>=0.6 AND BH-FDR q<0.05.
Network metrics compared between conditions; significance of the edge-count difference
by node-label permutation. Honest caveat: n=17/19 -> correlation networks are noisy,
so this is exploratory and reported as such.

Out: results/network_metrics.tsv, results/network_checks.txt ; figs/11_cooccurrence_networks.png
"""
import os, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import networkx as nx
from scipy import stats
from skbio.stats.composition import clr, multi_replace

A="/mnt/disk4/timo/gbi/analysis"; FIG=f"{A}/figs"; RES=f"{A}/results"
meta=pd.read_csv(f"{A}/metadata.tsv",sep="\t",index_col=0)
genus=pd.read_csv(f"{A}/genus_relabund.tsv",sep="\t",index_col=0)
samples=[s for s in genus.index if s in meta.index]; genus=genus.loc[samples]; md=meta.loc[samples]
checks=[]

# prevalent genera only (stabilises correlations)
prev=(genus>0).mean(0); keep=genus.columns[(prev>=0.25)&(genus.mean(0)>=0.1)]
G=genus[keep]
print(f"{len(samples)} samples, {len(keep)}/{genus.shape[1]} prevalent genera kept")

def clr_t(df):
    M=df.values.astype(float); M=M/M.sum(1,keepdims=True); M=multi_replace(M)
    return pd.DataFrame(clr(M),index=df.index,columns=df.columns)

def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]; m=len(p)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(q,0,1); return out

def build_net(ids, rho_thr=0.6, q_thr=0.05):
    sub=clr_t(G.loc[ids]); cols=list(sub.columns); n=len(cols)
    rhos=[]; ps=[]; pairs=[]
    for i in range(n):
        for j in range(i+1,n):
            rho,p=stats.spearmanr(sub.iloc[:,i],sub.iloc[:,j])
            if np.isnan(rho): rho,p=0,1
            rhos.append(rho); ps.append(p); pairs.append((cols[i],cols[j]))
    q=bh(ps)
    g=nx.Graph(); g.add_nodes_from(cols)
    for (a,b),rho,qq in zip(pairs,rhos,q):
        if abs(rho)>=rho_thr and qq<q_thr: g.add_edge(a,b,weight=rho)
    return g

def metrics(g):
    deg=dict(g.degree()); ne=g.number_of_edges(); nn=g.number_of_nodes()
    nconn=sum(1 for _,d in deg.items() if d>0)
    dens=nx.density(g)
    try:
        from networkx.algorithms.community import greedy_modularity_communities, modularity
        comms=list(greedy_modularity_communities(g)) if ne>0 else []
        mod=modularity(g,comms) if comms else 0
    except Exception: mod=np.nan
    pos=sum(1 for *_ ,w in g.edges(data="weight") if w>0)
    return {"nodes_total":nn,"nodes_connected":nconn,"edges":ne,"density":round(dens,4),
            "mean_degree":round(np.mean(list(deg.values())),3),"modularity":round(mod,3) if mod==mod else np.nan,
            "avg_clustering":round(nx.average_clustering(g),3),
            "frac_positive":round(pos/ne,3) if ne else np.nan}

PRE=[s for s in samples if md.loc[s,"condition"]=="pre-drought"]
DRO=[s for s in samples if md.loc[s,"condition"]=="drought"]
gp=build_net(PRE); gd=build_net(DRO)
mp=metrics(gp); mdd=metrics(gd)
res=pd.DataFrame({"pre_drought":mp,"drought":mdd}).T
res.to_csv(f"{RES}/network_metrics.tsv",sep="\t")
print("\nNetwork metrics:\n",res.to_string())

# permutation test on edge-count difference (shuffle condition labels)
obs=abs(gp.number_of_edges()-gd.number_of_edges()); allids=PRE+DRO; rng=np.random.default_rng(42)
perm=[]
for _ in range(199):
    sh=list(allids); rng.shuffle(sh)
    a=build_net(sh[:len(PRE)]); b=build_net(sh[len(PRE):])
    perm.append(abs(a.number_of_edges()-b.number_of_edges()))
p_edge=(1+sum(x>=obs for x in perm))/(1+len(perm))
print(f"\nEdge-count diff pre={gp.number_of_edges()} drought={gd.number_of_edges()} | perm p={p_edge:.3f}")
open(f"{RES}/network_metrics.tsv","a").write(f"\n# edge_diff_perm_p\t{p_edge:.4f}\n# obs_edge_diff\t{obs}\n")

# figure: two networks side by side
fig,axes=plt.subplots(1,2,figsize=(11,5))
for ax,(g,title,m) in zip(axes,[(gp,"pre-drought",mp),(gd,"drought",mdd)]):
    h=g.subgraph([n for n,d in g.degree() if d>0]) if g.number_of_edges() else g
    pos=nx.spring_layout(h,seed=1,k=0.5)
    ecol=["#2166ac" if w>0 else "#b2182b" for *_,w in h.edges(data="weight")]
    nx.draw_networkx_edges(h,pos,ax=ax,edge_color=ecol,alpha=.5,width=1.2)
    nx.draw_networkx_nodes(h,pos,ax=ax,node_size=70,node_color="#444",alpha=.8)
    ax.set_title(f"{title}\n{m['edges']} edges, density {m['density']}, modularity {m['modularity']}",fontsize=10)
    ax.axis("off")
fig.suptitle(f"Genus co-occurrence networks (|ρ|≥0.6, q<0.05)  ·  edge-diff perm p={p_edge:.2f}",fontsize=11)
plt.tight_layout(); plt.savefig(f"{FIG}/11_cooccurrence_networks.png",dpi=130); plt.close()

# checks
checks.append(("network_metrics.tsv", os.path.exists(f"{RES}/network_metrics.tsv")))
checks.append(("fig11 written", os.path.exists(f"{FIG}/11_cooccurrence_networks.png")))
checks.append(("perm p in (0,1]", 0<p_edge<=1))
checks.append(("both nets built", gp.number_of_nodes()>0 and gd.number_of_nodes()>0))
open(f"{RES}/network_checks.txt","w").write("\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks))
print("\n=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k,v in checks]
print(f"TIER4a {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")
