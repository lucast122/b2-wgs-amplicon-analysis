#!/usr/bin/env python3
"""
Full-nr Kaiju genus-level community profile + drought test, as an INDEPENDENT protein-based
cross-validation of the sylph (genome-containment) and 16S community profiles.

Reads kaiju_out_full/*.genus.tsv (kaiju2table; percent = % of TOTAL reads incl unclassified).
Renormalises each sample to % of CLASSIFIED reads for community comparison.

Outputs to analysis/results/:
  fullnr_genus_relabund.tsv      genus x sample, % of classified reads
  fullnr_genus_drought.tsv       pre vs drought MWU + BH-FDR + Cliff's delta
  fullnr_top_genera.tsv          top genera overall (mean % classified)
  fullnr_vs_sylph.tsv            rank comparison of top genera (kaiju vs sylph)
  fullnr_checks.txt
Figure: figs/13_fullnr_top_genera.png
"""
import os, glob, re, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

KOUT = "/mnt/disk4/timo/gbi/b2/classifications/kaiju_out_full"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

def clean_genus(lineage):
    toks = [t for t in str(lineage).split(";") if t and t != "NA"]
    g = toks[-1] if toks else None
    # drop kaiju catch-all bins (reads classified above genus rank are not a genus)
    if g and ("cannot be assigned" in g or "unclassified" in g.lower()): return None
    return g

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

def parse_one(fp):
    g = {}
    with open(fp) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        ci_reads = idx.get("reads", 2); ci_name = idx.get("taxon_name", 4)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= ci_name: continue
            name = p[ci_name]
            if name == "unclassified": continue
            gen = clean_genus(name)
            if not gen: continue
            try: reads = float(p[ci_reads])
            except (ValueError, IndexError): continue
            g[gen] = g.get(gen, 0.0) + reads
    return pd.Series(g)

def main():
    files = sorted(glob.glob(f"{KOUT}/*.genus.tsv"))
    if not files:
        print("No full-nr genus tables."); return
    print(f"parsing {len(files)} full-nr genus tables")
    cols = {}
    for fp in files:
        m = re.search(r"(SRR\d+)", os.path.basename(fp))
        if not m: continue
        s = parse_one(fp)
        if s.sum() > 0: cols[m.group(1)] = s / s.sum() * 100.0   # % of classified reads
    mat = pd.DataFrame(cols).fillna(0.0).T   # samples x genera
    mat.to_csv(f"{RES}/fullnr_genus_relabund.tsv", sep="\t")
    print(f"  matrix: {mat.shape[0]} samples x {mat.shape[1]} genera")

    top = mat.mean().sort_values(ascending=False).head(20)
    top.round(3).to_frame("mean_pct_classified").to_csv(f"{RES}/fullnr_top_genera.tsv", sep="\t")
    print("  top 8 genera (mean % classified):")
    for gname, v in top.head(8).items(): print(f"    {gname:24s} {v:5.2f}%")

    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    s = [x for x in mat.index if x in meta.index]; md = meta.loc[s]; mat = mat.loc[s]
    pre = list(md.index[md.condition == "pre-drought"]); dro = list(md.index[md.condition == "drought"])
    rows = []
    for gname in mat.columns:
        a = mat.loc[pre, gname].values; b = mat.loc[dro, gname].values
        if a.mean() < 0.01 and b.mean() < 0.01: continue
        try: _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError: p = np.nan
        rows.append({"genus": gname, "pre": a.mean(), "drought": b.mean(),
                     "log2FC": np.log2((b.mean()+1e-3)/(a.mean()+1e-3)),
                     "p": p, "cliffs_delta": cliffs(a, b)})
    dr = pd.DataFrame(rows)
    if len(dr):
        dr["q_BH"] = bh(dr["p"].values); dr = dr.sort_values("p")
        dr.round(5).to_csv(f"{RES}/fullnr_genus_drought.tsv", sep="\t", index=False)
    nsig = int((dr["q_BH"] < 0.05).sum()) if len(dr) else 0
    print(f"  drought test: {len(dr)} genera tested, {nsig} significant at q<0.05")

    # Cross-validate vs sylph genus_drought (rank of shared top genera)
    xcheck = pd.DataFrame()
    try:
        syl = pd.read_csv(f"{RES}/genus_drought.tsv", sep="\t")
        gcol = "genus" if "genus" in syl.columns else syl.columns[0]
        syl_top = syl.head(40)[gcol].astype(str).tolist()
        shared = [g for g in top.index if any(g in s2 or s2.endswith(g) for s2 in syl_top)]
        xcheck = pd.DataFrame({"kaiju_top_genus": top.index[:20],
                               "in_sylph_top40": [g in shared for g in top.index[:20]]})
        xcheck.to_csv(f"{RES}/fullnr_vs_sylph.tsv", sep="\t", index=False)
        print(f"  cross-check: {len(shared)}/{min(20,len(top))} kaiju-top genera also in sylph top-40")
    except Exception as e:
        print("  sylph cross-check skipped:", e)

    # Nitrospira specifically (the headline sylph/16S finding)
    nitro = [g for g in mat.columns if "Nitrospira" in g]
    nitro_line = "Nitrospira not detected"
    if nitro:
        g0 = nitro[0]; a = mat.loc[pre, g0].mean(); b = mat.loc[dro, g0].mean()
        nitro_line = f"Nitrospira: pre={a:.2f}% drought={b:.2f}% (log2FC={np.log2((b+1e-3)/(a+1e-3)):+.2f})"
    print("  " + nitro_line)

    # figure: top genera pre vs drought
    fig_ok = True
    try:
        gtop = top.head(14).index[::-1]
        yp = mat.loc[pre, gtop].mean(); yd = mat.loc[dro, gtop].mean(); y = np.arange(len(gtop))
        plt.figure(figsize=(7.5, 5.2))
        plt.barh(y-0.2, yp, 0.4, label="pre-drought", color="#1f77b4")
        plt.barh(y+0.2, yd, 0.4, label="drought", color="#d62728")
        plt.yticks(y, gtop, fontsize=8); plt.xlabel("% of classified reads (full-nr Kaiju)")
        plt.title("Top genera, protein-based (Kaiju nr_euk) — pre vs drought")
        plt.legend(frameon=False, fontsize=8); plt.tight_layout()
        plt.savefig(f"{FIG}/13_fullnr_top_genera.png", dpi=130); plt.close()
    except Exception as e:
        print("  figure failed:", e); fig_ok = False

    checks = [
        ("n_samples == 39", mat.shape[0] == 39),
        ("genera detected > 100", mat.shape[1] > 100),
        ("relabund rows ~100%", bool((mat.sum(1).sub(100).abs() < 1.0).all())),
        ("Nitrospira detected", len(nitro) > 0),
        ("drought table written", len(dr) > 0),
        ("figure written", fig_ok),
    ]
    open(f"{RES}/fullnr_checks.txt", "w").write(
        "\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k, v in checks) +
        f"\nn_significant_q05\t{nsig}\n{nitro_line}\n")
    print("=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k, v in checks]
    print(f"  {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")

if __name__ == "__main__":
    main()
