#!/usr/bin/env python3
"""
Post-process kaiju fungi-DB genus tables for the GBI samples -> AMF/EMF community.
Run: python kaiju_postprocess.py /mnt/disk4/timo/gbi/b2/classifications/kaiju_out_fungi

Outputs to /mnt/disk4/timo/gbi/b2/analysis/results/:
  kaiju_domain_fractions.tsv            (% of ALL reads by group, incl unclassified)
  kaiju_amf_matrix.tsv / kaiju_emf_matrix.tsv         (genus x sample, % of TOTAL reads)
  kaiju_amf_fungalnorm.tsv / kaiju_emf_fungalnorm.tsv (% of FUNGAL = AMF+EMF+Fungi_other reads)
  kaiju_amf_drought.tsv / kaiju_emf_drought.tsv        (pre vs drought, MWU + BH-FDR + Cliff's d)
  kaiju_checks.txt
Figure: figs/12_amf_emf_genera.png   (idx >=12 to avoid clobbering 09/10/11)
Excel : /mnt/disk4/timo/gbi/GBI_AMF_EMF_kaiju.xlsx  (all sheets)

AMF/EMF expressed as % of FUNGAL reads (not total) because ~96% of soil reads are
non-fungal and unclassified by the fungi DB. Drought tests on the fungal-normalised matrices.
"""
import os, glob, re, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

KOUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/disk4/timo/gbi/b2/classifications/kaiju_out_fungi"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
XLSX = "/mnt/disk4/timo/gbi/GBI_AMF_EMF_kaiju.xlsx"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

AMF_GENERA = ["Rhizophagus","Funneliformis","Glomus","Claroideoglomus","Diversispora",
    "Gigaspora","Scutellospora","Acaulospora","Paraglomus","Ambispora","Archaeospora",
    "Geosiphon","Racocetra","Dentiscutata","Cetraspora","Entrophospora","Oehlia",
    "Septoglomus","Sclerocystis","Redeckera","Pacispora","Dominikia","Rhizoglomus"]
EMF_GENERA = ["Cortinarius","Inocybe","Tricholoma","Russula","Lactarius","Amanita",
    "Boletus","Suillus","Laccaria","Pisolithus","Scleroderma","Hebeloma","Paxillus",
    "Cenococcum","Tuber","Thelephora","Tomentella","Hydnum","Cantharellus","Ramaria",
    "Sebacina","Clavulina","Hygrophorus","Entoloma","Craterellus","Hysterangium"]
AMF_CLADE = re.compile(r"Glomeromycot|Glomerale|Diversisporale|Archaeosporale|Paraglomerale", re.I)

def clean_genus(lineage):
    toks = [t for t in str(lineage).split(";") if t and t != "NA"]
    return toks[-1] if toks else str(lineage)

def classify(blob):
    for g in AMF_GENERA:
        if re.search(rf"\b{g}\b", blob): return "AMF"
    if AMF_CLADE.search(blob): return "AMF"
    for g in EMF_GENERA:
        if re.search(rf"\b{g}\b", blob): return "EMF"
    if re.search(r"Fungi|Basidiomycot|Ascomycot|Mucoromycot|Mortierellomycot|mycota|mycetes", blob): return "Fungi_other"
    if re.search(r"Bacteria", blob): return "Bacteria"
    if re.search(r"Archaea", blob): return "Archaea"
    if re.search(r"Viridiplantae|Streptophyta|Embryophyta", blob): return "Plant"
    return "Other"

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
    rows = []
    with open(fp) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        ci_reads = idx.get("reads", 2); ci_pct = idx.get("percent", 1); ci_name = idx.get("taxon_name", 4)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= ci_name: continue
            lineage = p[ci_name]
            try: reads = float(p[ci_reads]); pct = float(p[ci_pct])
            except (ValueError, IndexError): continue
            rows.append((clean_genus(lineage), reads, pct, classify(lineage)))
    return pd.DataFrame(rows, columns=["genus", "reads", "percent", "group"])

def main():
    files = sorted(glob.glob(f"{KOUT}/*.genus.tsv"))
    if not files:
        print("No *.genus.tsv yet."); return
    print(f"parsing {len(files)} genus tables")
    amf = {}; emf = {}; domain = {}; fungal = {}
    for fp in files:
        m = re.search(r"(SRR\d+)", os.path.basename(fp))
        if not m: continue
        acc = m.group(1); df = parse_one(fp)
        if df.empty: continue
        gsum = df.groupby("group")["percent"].sum(); domain[acc] = gsum
        fungal[acc] = float(gsum.get("AMF",0)+gsum.get("EMF",0)+gsum.get("Fungi_other",0))
        amf[acc] = df[df.group=="AMF"].groupby("genus")["percent"].sum()
        emf[acc] = df[df.group=="EMF"].groupby("genus")["percent"].sum()
    dom = pd.DataFrame(domain).T.fillna(0); dom.to_csv(f"{RES}/kaiju_domain_fractions.tsv", sep="\t")
    amf_m = pd.DataFrame(amf).T.fillna(0); emf_m = pd.DataFrame(emf).T.fillna(0)
    amf_m.to_csv(f"{RES}/kaiju_amf_matrix.tsv", sep="\t"); emf_m.to_csv(f"{RES}/kaiju_emf_matrix.tsv", sep="\t")
    ft = pd.Series(fungal).replace(0, np.nan)
    amf_fn = (amf_m.div(ft, axis=0)*100).fillna(0); emf_fn = (emf_m.div(ft, axis=0)*100).fillna(0)
    amf_fn.to_csv(f"{RES}/kaiju_amf_fungalnorm.tsv", sep="\t"); emf_fn.to_csv(f"{RES}/kaiju_emf_fungalnorm.tsv", sep="\t")
    print(f"AMF genera: {list(amf_m.columns)}")
    print(f"EMF genera: {list(emf_m.columns)}")
    print(f"mean fungal fraction of reads: {ft.mean():.2f}% ; AMF as % of fungi (mean): {amf_fn.sum(1).mean():.1f}")

    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    def drought_test(mat, label):
        s = [x for x in mat.index if x in meta.index]
        if not s: return pd.DataFrame()
        mat = mat.loc[s]; md = meta.loc[s]
        pre = [x for x in md.index[md.condition=="pre-drought"]]; dro = [x for x in md.index[md.condition=="drought"]]
        rows = []
        for g in mat.columns:
            a = mat.loc[pre, g].values; b = mat.loc[dro, g].values
            if a.mean()<1e-4 and b.mean()<1e-4: continue
            try: _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            except ValueError: p = np.nan
            rows.append({"genus":g,"pre":a.mean(),"drought":b.mean(),
                         "log2FC":np.log2((b.mean()+1e-3)/(a.mean()+1e-3)),"p":p,"cliffs_delta":cliffs(a,b)})
        r = pd.DataFrame(rows)
        if len(r): r["q_BH"] = bh(r["p"].values); r = r.sort_values("p")
        if len(r): r.to_csv(f"{RES}/kaiju_{label}_drought.tsv", sep="\t", index=False)
        print(f"\n{label.upper()} drought (fungal-norm):\n{r.head(8).to_string(index=False) if len(r) else '  (no data)'}")
        return r
    amf_dr = drought_test(amf_fn, "amf"); emf_dr = drought_test(emf_fn, "emf")

    # figure: top AMF + EMF genera (mean % of fungi), pre vs drought
    try:
        s = [x for x in amf_fn.index if x in meta.index]; md = meta.loc[s]
        pre = md.index[md.condition=="pre-drought"]; dro = md.index[md.condition=="drought"]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
        for ax, mat, title in [(axes[0], amf_fn.loc[s], "AMF genera"), (axes[1], emf_fn.loc[s], "EMF genera")]:
            top = mat.mean().sort_values(ascending=False).head(10).index[::-1]
            yp = mat.loc[[x for x in pre if x in mat.index], top].mean()
            yd = mat.loc[[x for x in dro if x in mat.index], top].mean()
            y = np.arange(len(top))
            ax.barh(y-0.2, yp, 0.4, label="pre-drought", color="#1f77b4")
            ax.barh(y+0.2, yd, 0.4, label="drought", color="#d62728")
            ax.set_yticks(y); ax.set_yticklabels(top, fontsize=8); ax.set_xlabel("% of fungal reads")
            ax.set_title(title); ax.legend(fontsize=8, frameon=False)
        plt.tight_layout(); plt.savefig(f"{FIG}/12_amf_emf_genera.png", dpi=130); plt.close()
        fig_ok = True
    except Exception as e:
        print("figure failed:", e); fig_ok = False

    # Excel deliverable
    xlsx_ok = True
    try:
        with pd.ExcelWriter(XLSX) as xl:
            dom.round(4).to_excel(xl, "domain_fractions")
            amf_fn.round(4).to_excel(xl, "AMF_pct_of_fungi"); emf_fn.round(4).to_excel(xl, "EMF_pct_of_fungi")
            amf_m.round(5).to_excel(xl, "AMF_pct_of_total"); emf_m.round(5).to_excel(xl, "EMF_pct_of_total")
            if len(amf_dr): amf_dr.round(5).to_excel(xl, "AMF_drought", index=False)
            if len(emf_dr): emf_dr.round(5).to_excel(xl, "EMF_drought", index=False)
        print(f"\nExcel written: {XLSX}")
    except Exception as e:
        print("Excel failed:", e); xlsx_ok = False

    checks = [
        ("n_samples_parsed", len(dom) > 0),
        ("AMF detected (Rhizophagus)", "Rhizophagus" in amf_m.columns),
        ("EMF detected", emf_m.shape[1] > 0),
        ("fungalnorm sums sane (<=100.1)", bool((amf_fn.sum(1) <= 100.1).all())),
        ("domain fractions row ~100", bool(((dom.sum(1)-100).abs() < 1).all())),
        ("figure 12 written", fig_ok),
        ("excel written", xlsx_ok),
    ]
    open(f"{RES}/kaiju_checks.txt","w").write("\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks))
    print("\n=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k,v in checks]
    print(f"samples parsed: {len(dom)}/39  |  {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")

if __name__ == "__main__":
    main()
