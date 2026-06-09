#!/usr/bin/env python3
"""
Independent 16S-amplicon test of the WGS drought signal.
Uses the GG2 phylum (L2) relative-abundance table from the 2026 amplicon
reanalysis, subset to the drought-labelled samples (Drought = drought vs pre),
and asks whether Actinomycetota / Nitrospirota move the same way as in the WGS.
"""
import os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

AMP = "/mnt/disk4/timo/gbi/amplicon/reanalysis_2026"
L2 = f"{AMP}/07_wgs_compare/amplicon_tables/gg2_L2_rel.tsv"
META = f"{AMP}/metadata/sample-metadata.tsv"
OUT = "/mnt/disk4/timo/gbi/analysis"; FIG = f"{OUT}/figs"; RES = f"{OUT}/results"

# ---- load L2 (skip '# Constructed' line; 2nd line header) ----
tab = pd.read_csv(L2, sep="\t", skiprows=1, index_col=0)
# rows = 'd__..;p__..' ; collapse to phylum name
def phy(x):
    parts = x.split(";")
    for p in parts:
        if p.startswith("p__"):
            return p[3:] or "Unassigned"
    return "Unassigned"
tab.index = [phy(i) for i in tab.index]
tab = tab.groupby(level=0).sum()              # taxa x samples (fractions 0-1)
tab = tab * 100.0                              # to percent

# ---- metadata ----
md = pd.read_csv(META, sep="\t")
md = md.rename(columns={"#SampleID": "sid"}).set_index("sid")
dro_ids = md.index[md["Drought"] == "drought"].tolist()
pre_ids = md.index[md["Drought"] == "pre"].tolist()
dro_ids = [s for s in dro_ids if s in tab.columns]
pre_ids = [s for s in pre_ids if s in tab.columns]
print(f"amplicon drought-labelled: {len(dro_ids)} drought, {len(pre_ids)} pre")

# ---- test each phylum drought vs pre ----
focus = ["Actinomycetota", "Nitrospirota", "Pseudomonadota", "Acidobacteriota",
         "Chloroflexota", "Methylomirabilota", "Thermoproteota", "Planctomycetota"]
rows = []
for ph in tab.index:
    if ph not in tab.index: continue
    a = tab.loc[ph, pre_ids].astype(float)
    b = tab.loc[ph, dro_ids].astype(float)
    if a.mean() < 0.1 and b.mean() < 0.1: continue
    try: _, pv = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError: pv = np.nan
    rows.append({"phylum": ph, "pre_mean": a.mean(), "drought_mean": b.mean(),
                 "log2FC": np.log2((b.mean()+1e-3)/(a.mean()+1e-3)), "MWU_p": pv})
res = pd.DataFrame(rows).sort_values("MWU_p")
res.to_csv(f"{RES}/amplicon_phylum_drought.tsv", sep="\t", index=False)
print("\n16S amplicon drought response (top):")
print(res.head(10).to_string(index=False))

# ---- focused figure: same phyla as WGS headline ----
# GG2 names carry GTDB suffixes (e.g. Nitrospirota_A_437815); map to clean labels.
ri = res.set_index("phylum")
wanted = [("Actinomycetota","Actinomycetota"),
          ("Nitrospirota_A_437815","Nitrospirota"),
          ("Pseudomonadota","Pseudomonadota"),
          ("Acidobacteriota","Acidobacteriota"),
          ("Chloroflexota","Chloroflexota"),
          ("Methylomirabilota","Methylomirabilota")]
wanted = [(k, lab) for k, lab in wanted if k in ri.index]
show = ri.loc[[k for k, _ in wanted]].reset_index()
show["label"] = [lab for _, lab in wanted]
fig, ax = plt.subplots(figsize=(6.8, 4))
x = np.arange(len(show)); w = 0.38
ax.bar(x - w/2, show["pre_mean"], w, label="pre-drought", color="#1f77b4")
ax.bar(x + w/2, show["drought_mean"], w, label="drought", color="#d62728")
for i, (_, r) in enumerate(show.iterrows()):
    star = "*" if (pd.notna(r.MWU_p) and r.MWU_p < 0.05) else ""
    ax.text(i, max(r.pre_mean, r.drought_mean)+0.3,
            f"p={r.MWU_p:.2g}{star}" if pd.notna(r.MWU_p) else "", ha="center", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(show["label"], rotation=35, ha="right", fontsize=8)
ax.set_ylabel("mean relative abundance (%)")
ax.set_title(f"16S amplicon drought response\n({len(pre_ids)} pre vs {len(dro_ids)} drought, GG2; * p<0.05)")
ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(f"{FIG}/06_amplicon_drought.png", dpi=130); plt.close()

# ---- concordance note: direction agreement with WGS ----
wgs = pd.read_csv(f"{RES}/phylum_shift_pre_vs_drought.tsv", sep="\t").set_index("phylum")
print("\nDirection concordance (WGS vs 16S), shared phyla:")
for ph in ["Actinomycetota","Nitrospirota","Acidobacteriota","Methylomirabilota","Planctomycetota"]:
    if ph in wgs.index and ph in res.set_index("phylum").index:
        wd = "UP" if wgs.loc[ph,"drought_mean"]>wgs.loc[ph,"pre_mean"] else "DOWN"
        rr = res.set_index("phylum").loc[ph]
        ad = "UP" if rr["drought_mean"]>rr["pre_mean"] else "DOWN"
        print(f"  {ph:18s} WGS {wd:4s} | 16S {ad:4s} | {'AGREE' if wd==ad else 'differ'}")
print("\nDONE amplicon_drought")
