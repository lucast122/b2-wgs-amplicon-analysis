#!/usr/bin/env python3
"""
Post-process DIAMOND -f102 (taxonomic LCA) per-sample outputs -> genus community matrix,
then CROSS-CHECK against the full-nr Kaiju genus profile (two independent protein-alignment
methods on the SAME cleaned nr_euk DB; difference = aligner/algorithm, not DB).

Inputs:
  diamond_out/SRR*.diamond_lca.tsv   (qseqid  taxid  evalue ; taxid 0 = unclassified)
  kaiju_db/{nodes,names}.dmp         (taxid -> rank/parent ; taxid -> name)
  b2/analysis/results/fullnr_genus_relabund.tsv   (kaiju full-nr genus matrix, % classified)
Outputs to b2/analysis/results/:
  diamond_genus_relabund.tsv         genus x sample, % of classified reads
  diamond_top_genera.tsv
  diamond_vs_kaiju.tsv               per-sample Spearman + top-genera overlap
  diamond_checks.txt
Figure: figs/14_diamond_vs_kaiju.png
"""
import os, glob, re, sys
from collections import Counter
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

DOUT = "/mnt/disk4/timo/gbi/diamond_out"
TAX  = "/mnt/disk4/timo/gbi/kaiju_db"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

def load_taxonomy():
    parent = {}; rank = {}
    with open(f"{TAX}/nodes.dmp") as f:
        for line in f:
            p = [x.strip() for x in line.split("|")]
            tid = int(p[0]); parent[tid] = int(p[1]); rank[tid] = p[2]
    name = {}
    with open(f"{TAX}/names.dmp") as f:
        for line in f:
            p = [x.strip() for x in line.split("|")]
            if p[3] == "scientific name":
                name[int(p[0])] = p[1]
    return parent, rank, name

def main():
    print("loading taxonomy (nodes/names.dmp) ...")
    parent, rank, name = load_taxonomy()
    cache = {}
    def to_genus(tid):
        if tid in cache: return cache[tid]
        seen = []; t = tid
        g = None
        for _ in range(100):
            seen.append(t)
            if rank.get(t) == "genus": g = name.get(t); break
            pt = parent.get(t)
            if pt is None or pt == t: break
            t = pt
        for s in seen: cache[s] = g
        return g

    files = sorted(glob.glob(f"{DOUT}/SRR*.diamond_lca.tsv"))
    print(f"parsing {len(files)} diamond outputs")
    cols = {}
    for fp in files:
        acc = re.search(r"(SRR\d+)", os.path.basename(fp)).group(1)
        tc = Counter()
        with open(fp) as f:
            for line in f:
                a = line.split("\t")
                if len(a) < 2: continue
                try: tid = int(a[1])
                except ValueError: continue
                if tid > 0: tc[tid] += 1
        gc = Counter()
        for tid, n in tc.items():
            g = to_genus(tid)
            if g: gc[g] += n
        tot = sum(gc.values())
        if tot > 0:
            cols[acc] = pd.Series({g: v/tot*100 for g, v in gc.items()})
        print(f"  {acc}: {sum(tc.values())} classified reads -> {len(gc)} genera")
    mat = pd.DataFrame(cols).fillna(0.0).T
    mat.to_csv(f"{RES}/diamond_genus_relabund.tsv", sep="\t")
    top = mat.mean().sort_values(ascending=False).head(20)
    top.round(3).to_frame("mean_pct_classified").to_csv(f"{RES}/diamond_top_genera.tsv", sep="\t")
    print(f"\ndiamond matrix: {mat.shape[0]} samples x {mat.shape[1]} genera")
    print("top 8 diamond genera:")
    for g, v in top.head(8).items(): print(f"  {g:24s} {v:5.2f}%")

    # cross-check vs kaiju full-nr
    kpath = f"{RES}/fullnr_genus_relabund.tsv"
    xrows = []; rho_overall = np.nan; top_overlap = np.nan
    if os.path.exists(kpath):
        kaj = pd.read_csv(kpath, sep="\t", index_col=0)
        shared_s = [s for s in mat.index if s in kaj.index]
        shared_g = [g for g in mat.columns if g in kaj.columns]
        print(f"\nshared with kaiju: {len(shared_s)} samples, {len(shared_g)} genera")
        for s in shared_s:
            a = mat.loc[s, shared_g].values; b = kaj.loc[s, shared_g].values
            if a.sum() > 0 and b.sum() > 0:
                rho, _ = stats.spearmanr(a, b)
                xrows.append({"sample": s, "spearman_rho": rho})
        xdf = pd.DataFrame(xrows)
        if len(xdf):
            xdf.round(4).to_csv(f"{RES}/diamond_vs_kaiju.tsv", sep="\t", index=False)
            rho_overall = xdf["spearman_rho"].median()
        # top-genera overlap
        dt = set(mat.mean().sort_values(ascending=False).head(15).index)
        kt = set(kaj.mean().sort_values(ascending=False).head(15).index)
        top_overlap = len(dt & kt) / 15.0
        print(f"  median per-sample Spearman(diamond,kaiju) = {rho_overall:.3f}")
        print(f"  top-15 genera overlap = {len(dt & kt)}/15 ({top_overlap*100:.0f}%)")
        print(f"  shared top genera: {sorted(dt & kt)[:10]}")

    # Nitrospira drought direction (descriptive; pseudoreplication caveat applies)
    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    nit_line = "Nitrospira not detected"
    nf = [g for g in mat.columns if "Nitrospira" in g]
    if nf:
        s = [x for x in mat.index if x in meta.index]; md = meta.loc[s]
        pre = md.index[md.condition=="pre-drought"]; dro = md.index[md.condition=="drought"]
        a = mat.loc[[x for x in pre], nf[0]].mean(); b = mat.loc[[x for x in dro], nf[0]].mean()
        nit_line = f"Nitrospira (diamond): pre={a:.2f}% drought={b:.2f}% (log2FC={np.log2((b+1e-3)/(a+1e-3)):+.2f})"
    print("  " + nit_line)

    # figure: diamond vs kaiju top-genera scatter (mean abundances)
    fig_ok = True
    try:
        if os.path.exists(kpath):
            common = [g for g in mat.columns if g in kaj.columns]
            dm = mat[common].mean(); km = kaj[common].mean()
            sel = (dm + km).sort_values(ascending=False).head(30).index
            plt.figure(figsize=(5.6,5.4))
            plt.scatter(km[sel], dm[sel], s=18, alpha=.7)
            mx = max(km[sel].max(), dm[sel].max())*1.05
            plt.plot([0,mx],[0,mx],'k--',lw=.8,alpha=.5)
            for g in sel[:8]: plt.annotate(g, (km[g], dm[g]), fontsize=6)
            plt.xlabel("Kaiju full-nr (% classified)"); plt.ylabel("DIAMOND (% classified)")
            plt.title(f"Top genera: DIAMOND vs Kaiju (same nr_euk DB)\nmedian per-sample Spearman={rho_overall:.2f}")
            plt.tight_layout(); plt.savefig(f"{FIG}/14_diamond_vs_kaiju.png", dpi=130); plt.close()
    except Exception as e:
        print("figure failed:", e); fig_ok = False

    checks = [
        ("n_samples == 39", mat.shape[0] == 39),
        ("genera > 100", mat.shape[1] > 100),
        ("relabund rows ~100%", bool((mat.sum(1).sub(100).abs() < 1.0).all())),
        ("kaiju cross-check ran", len(xrows) > 0),
        ("median Spearman > 0.5", (not np.isnan(rho_overall)) and rho_overall > 0.5),
        ("Nitrospira detected", len(nf) > 0),
        ("figure written", fig_ok),
    ]
    open(f"{RES}/diamond_checks.txt","w").write(
        "\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks) +
        f"\nmedian_spearman_vs_kaiju\t{rho_overall:.4f}\ntop15_overlap\t{top_overlap:.3f}\n{nit_line}\n")
    print("\n=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k,v in checks]
    print(f"  {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")

if __name__ == "__main__":
    main()
