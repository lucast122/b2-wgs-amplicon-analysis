#!/usr/bin/env python3
"""
Post-process kaiju2table genus output for all 39 GBI samples (runs once the
custom nr_euk+AMF+EMF kaiju classification completes).

Input : /mnt/disk4/timo/gbi/kaiju_out/<acc>.genus.tsv
        (kaiju2table -r genus -l superkingdom,phylum,class,order,family,genus)
        columns: file, percent, reads, taxon_id, taxon_name, taxon_path
Output: /mnt/disk4/timo/gbi/analysis/results/
        kaiju_amf_matrix.tsv, kaiju_emf_matrix.tsv, kaiju_domain_fractions.tsv,
        kaiju_amf_drought.tsv, kaiju_emf_drought.tsv
        figs/09_amf_genera.png, figs/10_kaiju_domains.png

Safe to run repeatedly; skips missing samples. Designed so the report's
AMF/EMF section fills in automatically.
"""
import os, glob, re
import numpy as np, pandas as pd

KOUT = "/mnt/disk4/timo/gbi/kaiju_out"
A = "/mnt/disk4/timo/gbi/analysis"
RES = f"{A}/results"; FIG = f"{A}/figs"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# --- AMF and EMF genus dictionaries (match against taxon_name or lineage) ---
AMF_GENERA = ["Rhizophagus","Funneliformis","Glomus","Claroideoglomus","Diversispora",
    "Gigaspora","Scutellospora","Acaulospora","Paraglomus","Ambispora","Archaeospora",
    "Geosiphon","Racocetra","Dentiscutata","Cetraspora","Entrophospora","Oehlia",
    "Septoglomus","Sclerocystis","Redeckera","Pacispora","Dominikia","Rhizoglomus"]
EMF_GENERA = ["Cortinarius","Inocybe","Tricholoma","Russula","Lactarius","Amanita",
    "Boletus","Suillus","Laccaria","Pisolithus","Scleroderma","Hebeloma","Paxillus",
    "Cenococcum","Tuber","Thelephora","Tomentella","Hydnum","Cantharellus","Ramaria",
    "Sebacina","Clavulina","Hygrophorus","Entoloma","Craterellus","Hysterangium"]
# Glomeromycota = AMF clade marker in lineage
AMF_CLADE = re.compile(r"Glomeromycot|Glomerale|Diversisporale|Archaeosporale|Paraglomerale", re.I)

def classify(name, path):
    blob = f"{name} {path}"
    for g in AMF_GENERA:
        if re.search(rf"\b{g}\b", blob): return "AMF"
    if AMF_CLADE.search(blob): return "AMF"
    for g in EMF_GENERA:
        if re.search(rf"\b{g}\b", blob): return "EMF"
    if re.search(r"Fungi|Basidiomycot|Ascomycot|Mucoromycot|Mortierellomycot", blob): return "Fungi_other"
    if re.search(r"Bacteria", blob): return "Bacteria"
    if re.search(r"Archaea", blob): return "Archaea"
    if re.search(r"Viridiplantae|Streptophyta|Embryophyta", blob): return "Plant"
    return "Other"

def parse_one(fp):
    """Return DataFrame: genus, reads, percent, group for one sample."""
    rows = []
    with open(fp) as f:
        header = f.readline().rstrip("\n").split("\t")
        # locate columns robustly
        idx = {c: i for i, c in enumerate(header)}
        ci_reads = idx.get("reads", 2); ci_pct = idx.get("percent", 1)
        ci_name = idx.get("taxon_name", 4)
        ci_path = idx.get("taxon_path", len(header)-1)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= ci_name: continue
            name = p[ci_name]; path = p[ci_path] if len(p) > ci_path else ""
            try: reads = float(p[ci_reads]); pct = float(p[ci_pct])
            except ValueError: continue
            rows.append((name, reads, pct, classify(name, path)))
    return pd.DataFrame(rows, columns=["genus","reads","percent","group"])

def main():
    files = sorted(glob.glob(f"{KOUT}/*.genus.tsv"))
    if not files:
        print("No *.genus.tsv yet — kaiju run not finished. Nothing to do."); return
    print(f"parsing {len(files)} genus tables")

    amf = {}; emf = {}; domain = {}
    for fp in files:
        acc = re.search(r"(SRR\d+)", os.path.basename(fp)).group(1)
        df = parse_one(fp)
        if df.empty: continue
        # domain fractions (% of classified reads by group)
        gsum = df.groupby("group")["percent"].sum()
        domain[acc] = gsum
        amf[acc] = df[df.group=="AMF"].groupby("genus")["percent"].sum()
        emf[acc] = df[df.group=="EMF"].groupby("genus")["percent"].sum()

    dom = pd.DataFrame(domain).T.fillna(0)
    dom.to_csv(f"{RES}/kaiju_domain_fractions.tsv", sep="\t")
    amf_m = pd.DataFrame(amf).T.fillna(0); amf_m.to_csv(f"{RES}/kaiju_amf_matrix.tsv", sep="\t")
    emf_m = pd.DataFrame(emf).T.fillna(0); emf_m.to_csv(f"{RES}/kaiju_emf_matrix.tsv", sep="\t")
    print(f"AMF genera detected: {list(amf_m.columns)}")
    print(f"EMF genera detected: {list(emf_m.columns)}")
    print(f"domain fractions (mean %):\n{dom.mean().round(3)}")

    # --- drought differential on AMF/EMF genera ---
    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    from scipy import stats
    def drought_test(mat, label):
        s = [x for x in mat.index if x in meta.index]
        mat = mat.loc[s]; md = meta.loc[s]
        pre = md.index[md.condition=="pre-drought"]; dro = md.index[md.condition=="drought"]
        rows=[]
        for g in mat.columns:
            a = mat.loc[[x for x in pre if x in mat.index], g]
            b = mat.loc[[x for x in dro if x in mat.index], g]
            if a.mean()<1e-4 and b.mean()<1e-4: continue
            try: _,p = stats.mannwhitneyu(a,b,alternative="two-sided")
            except ValueError: p=np.nan
            rows.append({"genus":g,"pre":a.mean(),"drought":b.mean(),"MWU_p":p})
        r = pd.DataFrame(rows).sort_values("MWU_p") if rows else pd.DataFrame()
        if len(r): r.to_csv(f"{RES}/kaiju_{label}_drought.tsv", sep="\t", index=False)
        print(f"\n{label.upper()} drought test:\n{r.to_string(index=False) if len(r) else '  (no genera)'}")
        return r

    drought_test(amf_m, "amf"); drought_test(emf_m, "emf")
    print("\nDONE kaiju_postprocess — matrices + drought tests written to results/")

if __name__ == "__main__":
    main()
