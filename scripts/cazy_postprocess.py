#!/usr/bin/env python3
"""
Tier-2 CAZyme profiling. Parses DIAMOND-vs-CAZyDB (dbCAN) hits -> CAZy families & classes,
normalised to hits per million subsample reads (HPM), then tests the drought contrast.
CAZy class is in the subject seqid (e.g. '...|GH5|CBM6'): GH=glycoside hydrolases (degradation),
GT=transferases (synthesis), PL=lyases, CE=esterases, AA=auxiliary/redox (ligninolytic),
CBM=binding modules. Carbohydrate-degradation potential ties to the drought C-allocation story.
"""
import os, glob, re
from collections import defaultdict
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CZ = "/mnt/disk4/timo/gbi/b2/functional/cazy"; HITS = f"{CZ}/hits"
DOUT = "/mnt/disk4/timo/gbi/diamond_out"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
FAM_RE = re.compile(r"(GH|GT|PL|CE|AA|CBM)(\d+)")

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
    files = sorted(glob.glob(f"{HITS}/*.cazy.tsv"))
    print(f"parsing {len(files)} CAZy hit files")
    fam = {}; cls = {}
    for fp in files:
        acc = re.search(r"(SRR\d+)", os.path.basename(fp)).group(1)
        cf = defaultdict(int); cc = defaultdict(int)
        for line in open(fp):
            p = line.split("\t")
            if len(p) < 2: continue
            seen = set()
            for m in FAM_RE.finditer(p[1]):
                f = m.group(1) + m.group(2)           # e.g. GH5
                if f in seen: continue
                seen.add(f); cf[f] += 1; cc[m.group(1)] += 1
        if cf: fam[acc] = pd.Series(cf); cls[acc] = pd.Series(cc)
    reads = {acc: (sum(1 for _ in open(f"{DOUT}/{acc}.diamond_lca.tsv")) if os.path.exists(f"{DOUT}/{acc}.diamond_lca.tsv") else 5_000_000) for acc in fam}
    fmat = pd.DataFrame(fam).fillna(0).T
    cmat = pd.DataFrame(cls).fillna(0).T
    hpm_f = fmat.div(pd.Series(reads).reindex(fmat.index)/1e6, axis=0)
    hpm_c = cmat.div(pd.Series(reads).reindex(cmat.index)/1e6, axis=0)
    hpm_f.to_csv(f"{RES}/cazy_family_hpm.tsv", sep="\t"); hpm_c.to_csv(f"{RES}/cazy_class_hpm.tsv", sep="\t")
    print(f"matrix: {hpm_f.shape[0]} samples x {hpm_f.shape[1]} families ; classes: {list(hpm_c.columns)}")
    print("mean class HPM:")
    for k, v in hpm_c.mean().sort_values(ascending=False).items(): print(f"  {k:5s} {v:9.1f}")

    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    cd = drought_test(hpm_c, meta, "cazy_class_drought")
    fd = drought_test(hpm_f, meta, "cazy_family_drought")
    s=[x for x in hpm_c.index if x in meta.index]; md=meta.loc[s]
    pre=[x for x in md.index[md.condition=="pre-drought"]]; dro=[x for x in md.index[md.condition=="drought"]]
    print("\nclass pre vs drought (HPM):")
    for c in hpm_c.columns:
        print(f"  {c:5s} pre={hpm_c.loc[pre,c].mean():8.1f} drought={hpm_c.loc[dro,c].mean():8.1f}")
    nfsig = int((fd["q_BH"]<0.05).sum()) if len(fd) else 0
    ncsig = int((cd["q_BH"]<0.05).sum()) if len(cd) else 0
    print(f"  significant (q<0.05): {ncsig} classes, {nfsig} families")

    fig_ok=True
    try:
        cl = hpm_c.mean().sort_values(ascending=False).index
        yp=hpm_c.loc[pre,cl].mean(); yd=hpm_c.loc[dro,cl].mean(); y=np.arange(len(cl))
        plt.figure(figsize=(7,4))
        plt.barh(y-0.2,yp,0.4,label="pre-drought",color="#1f77b4")
        plt.barh(y+0.2,yd,0.4,label="drought",color="#d62728")
        names={"GH":"GH glycoside hydrolases","GT":"GT transferases","PL":"PL lyases","CE":"CE esterases","AA":"AA auxiliary/redox","CBM":"CBM binding"}
        plt.yticks(y,[names.get(c,c) for c in cl],fontsize=8); plt.xlabel("hits per million reads")
        plt.title("CAZyme classes — pre vs drought"); plt.legend(frameon=False,fontsize=8,loc="lower right")
        plt.tight_layout(); plt.savefig(f"{FIG}/16_cazy_classes.png",dpi=130); plt.close()
    except Exception as e:
        print("fig failed:",e); fig_ok=False

    checks=[
        ("n_samples >= 36", hpm_f.shape[0]>=36),
        ("families detected > 100", hpm_f.shape[1]>100),
        ("GH class present", "GH" in hpm_c.columns),
        ("class+family drought tables", len(cd)>0 and len(fd)>0),
        ("figure written", fig_ok),
    ]
    open(f"{RES}/cazy_checks.txt","w").write(
        "\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks)+
        f"\nn_class_signif_q05\t{ncsig}\nn_family_signif_q05\t{nfsig}\n")
    print("\n=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k,v in checks]
    print(f"  {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")

if __name__ == "__main__":
    main()
