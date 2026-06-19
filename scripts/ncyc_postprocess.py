#!/usr/bin/env python3
"""
Tier-2 N-cycle functional profiling. Maps DIAMOND-vs-NCyc hits -> gene families -> N-cycle
processes, normalised to hits per million subsample reads (HPM), then tests the drought contrast.
Headline link: nxrA/nxrB (nitrite oxidation, the Nitrospira gene) vs the taxonomic Nitrospira decline.

Inputs:
  functional/ncyc/hits/SRR*.ncyc.tsv   (qseqid sseqid pident length evalue bitscore)
  functional/ncyc/id2gene.map          (sseqid -> gene family)
  diamond_out/SRR*.diamond_lca.tsv      (line count = subsample read total, for normalisation)
Outputs to b2/analysis/results/:
  ncyc_genefamily_hpm.tsv  ncyc_process_hpm.tsv  ncyc_drought.tsv  ncyc_process_drought.tsv  ncyc_checks.txt
Figure: figs/15_ncyc_drought.png
"""
import os, glob, re
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

NC = "/mnt/disk4/timo/gbi/b2/functional/ncyc"
HITS = f"{NC}/hits"; MAP = f"{NC}/id2gene.map"
DOUT = "/mnt/disk4/timo/gbi/diamond_out"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# gene family -> N-cycle process (NCyc; Tu et al. 2019)
def proc_of(g):
    if re.match(r"amo[ABC]|hao|nxr[AB]", g): return "Nitrification"
    if re.match(r"nar[GHI]|nap[AB]|nir[KS]|nor[BC]|nosZ", g): return "Denitrification"
    if re.match(r"nrf[ABCD]|nir[BD]", g): return "DNRA"
    if re.match(r"nif[DHKW]|anfG|vnf", g): return "N-fixation"
    if re.match(r"hzs[ABC]|hdh", g): return "Anammox"
    if re.match(r"narB|nas[AB]|nirA|NR|NasA", g): return "Assim. nitrate red."
    return "N-assimilation/other"

def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); r = p[o]; m = len(p)
    if m == 0: return p
    q = np.minimum.accumulate((r*m/np.arange(1, m+1))[::-1])[::-1]
    out = np.empty(m); out[o] = np.clip(q, 0, 1); return out

def cliffs(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0: return np.nan
    gt = sum((bj > a).sum() for bj in b); lt = sum((bj < a).sum() for bj in b)
    return (gt - lt) / (len(a) * len(b))

def drought_test(mat, meta, label):
    s = [x for x in mat.index if x in meta.index]; md = meta.loc[s]; mat = mat.loc[s]
    pre = [x for x in md.index[md.condition=="pre-drought"]]; dro = [x for x in md.index[md.condition=="drought"]]
    rows = []
    for c in mat.columns:
        a = mat.loc[pre, c].values; b = mat.loc[dro, c].values
        if a.mean() < 1e-6 and b.mean() < 1e-6: continue
        try: _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError: p = np.nan
        rows.append({"feature":c,"pre":a.mean(),"drought":b.mean(),
                     "log2FC":np.log2((b.mean()+1e-3)/(a.mean()+1e-3)),"p":p,"cliffs_delta":cliffs(a,b)})
    r = pd.DataFrame(rows)
    if len(r): r["q_BH"]=bh(r["p"].values); r=r.sort_values("p"); r.round(4).to_csv(f"{RES}/{label}.tsv",sep="\t",index=False)
    return r

def main():
    id2g = {}
    for line in open(MAP):
        a = line.rstrip("\n").split("\t")
        if len(a) >= 2: id2g[a[0]] = a[1]
    print(f"id2gene: {len(id2g)} seqids -> {len(set(id2g.values()))} families")

    files = sorted(glob.glob(f"{HITS}/*.ncyc.tsv"))
    print(f"parsing {len(files)} hit files")
    fam = {}; tot_hits = {}
    for fp in files:
        acc = re.search(r"(SRR\d+)", os.path.basename(fp)).group(1)
        c = defaultdict(int); n = 0
        for line in open(fp):
            p = line.split("\t")
            if len(p) < 2: continue
            g = id2g.get(p[1])
            if g: c[g] += 1; n += 1
        if n: fam[acc] = pd.Series(c); tot_hits[acc] = n
    # normalise to hits per million subsample reads
    reads = {}
    for acc in fam:
        lp = f"{DOUT}/{acc}.diamond_lca.tsv"
        reads[acc] = sum(1 for _ in open(lp)) if os.path.exists(lp) else 5_000_000
    gmat = pd.DataFrame(fam).fillna(0).T
    hpm = gmat.div(pd.Series(reads).reindex(gmat.index)/1e6, axis=0)   # hits per million reads
    hpm.to_csv(f"{RES}/ncyc_genefamily_hpm.tsv", sep="\t")
    # process-level
    pmap = {g: proc_of(g) for g in hpm.columns}
    pmat = hpm.T.groupby(pmap).sum().T
    pmat.to_csv(f"{RES}/ncyc_process_hpm.tsv", sep="\t")
    print(f"matrix: {hpm.shape[0]} samples x {hpm.shape[1]} gene families; {pmat.shape[1]} processes")
    print("mean process HPM:")
    for k, v in pmat.mean().sort_values(ascending=False).items(): print(f"  {k:24s} {v:8.1f}")

    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    gd = drought_test(hpm, meta, "ncyc_drought")
    pd_ = drought_test(pmat, meta, "ncyc_process_drought")

    # nxr (Nitrospira nitrite oxidation) focus
    s=[x for x in hpm.index if x in meta.index]; md=meta.loc[s]
    pre=[x for x in md.index[md.condition=="pre-drought"]]; dro=[x for x in md.index[md.condition=="drought"]]
    nxr_cols=[c for c in hpm.columns if c.startswith("nxr")]
    nxr_line="nxr not detected"
    if nxr_cols:
        nxr=hpm[nxr_cols].sum(1); a=nxr.loc[pre].mean(); b=nxr.loc[dro].mean()
        nxr_line=f"nxrAB (nitrite oxidation, Nitrospira): pre={a:.1f} drought={b:.1f} HPM (log2FC={np.log2((b+1e-3)/(a+1e-3)):+.2f}, Cliff d={cliffs(nxr.loc[pre],nxr.loc[dro]):+.2f})"
    print("  " + nxr_line)
    amo_cols=[c for c in hpm.columns if c.startswith("amo")]
    if amo_cols:
        amo=hpm[amo_cols].sum(1); print(f"  amoABC (ammonia ox.): pre={amo.loc[pre].mean():.1f} drought={amo.loc[dro].mean():.1f} HPM")

    # figure: process HPM pre vs drought
    fig_ok=True
    try:
        procs=pmat.mean().sort_values(ascending=False).index
        yp=pmat.loc[pre,procs].mean(); yd=pmat.loc[dro,procs].mean(); y=np.arange(len(procs))
        plt.figure(figsize=(7.2,4.4))
        plt.barh(y-0.2,yp,0.4,label="pre-drought",color="#1f77b4")
        plt.barh(y+0.2,yd,0.4,label="drought",color="#d62728")
        plt.yticks(y,procs,fontsize=8); plt.xlabel("hits per million reads (NCyc)")
        plt.title("N-cycle process potential — pre vs drought"); plt.legend(frameon=False,fontsize=8)
        plt.tight_layout(); plt.savefig(f"{FIG}/15_ncyc_drought.png",dpi=130); plt.close()
    except Exception as e:
        print("fig failed:",e); fig_ok=False

    checks=[
        ("n_samples >= 36", hpm.shape[0]>=36),
        ("gene families detected > 20", hpm.shape[1]>20),
        ("nxr detected", len(nxr_cols)>0),
        ("nitrification process present", "Nitrification" in pmat.columns),
        ("drought tables written", len(gd)>0 and len(pd_)>0),
        ("figure written", fig_ok),
    ]
    open(f"{RES}/ncyc_checks.txt","w").write(
        "\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks)+f"\n{nxr_line}\nn_families_signif_q05\t{int((gd['q_BH']<0.05).sum()) if len(gd) else 0}\n")
    print("\n=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k,v in checks]
    print(f"  {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")

if __name__ == "__main__":
    main()
