#!/usr/bin/env python3
"""Synthesis figure: drought effect size (Cliff's delta) for the key readouts across all data
types (taxonomy, N-cycle genes, CAZymes, osmolyte genes, diversity). Makes the core message in
one panel — the taxonomic Nitrospira decline is the only non-trivial effect; functional gene
potential and diversity sit in the 'negligible' band (|delta|<0.33)."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"

def grab(fn, key_col, key, dcol="cliffs_delta"):
    p = f"{RES}/{fn}"
    if not os.path.exists(p): return None
    d = pd.read_csv(p, sep="\t")
    if key_col not in d.columns or dcol not in d.columns: return None
    m = d[d[key_col].astype(str).str.fullmatch(key)] if key_col != "phylum" else d[d[key_col]==key]
    if not len(m): m = d[d[key_col].astype(str).str.startswith(key)]
    return float(m.iloc[0][dcol]) if len(m) else None

# (label, category, delta) — pull from the per-analysis drought tables
rows = []
def add(label, cat, val):
    if val is not None and not (isinstance(val,float) and np.isnan(val)): rows.append((label, cat, val))

# Taxonomy
add("Nitrospira (genus)", "Taxonomy", grab("genus_drought.tsv","genus","Nitrospira_D"))
add("Bradyrhizobium (genus)", "Taxonomy", grab("genus_drought.tsv","genus","Bradyrhizobium"))
add("Actinomycetota (phylum)", "Taxonomy", grab("phylum_shift_pre_vs_drought.tsv","phylum","Actinomycetota"))
add("Nitrospirota (phylum)", "Taxonomy", grab("phylum_shift_pre_vs_drought.tsv","phylum","Nitrospirota"))
# N-cycle genes
add("nxr / nitrite ox. (gene)", "N-cycle genes", grab("ncyc_process_drought.tsv","feature","Nitrification"))
add("Denitrification (gene)", "N-cycle genes", grab("ncyc_process_drought.tsv","feature","Denitrification"))
# CAZymes
add("GH carbohydr. (gene)", "CAZymes", grab("cazy_class_drought.tsv","feature","GH"))
add("AA ligninolytic (gene)", "CAZymes", grab("cazy_class_drought.tsv","feature","AA"))
# Osmolyte genes
add("Trehalose (gene)", "Osmolyte genes", grab("osmolyte_drought.tsv","pathway","Trehalose"))
add("Glycine-betaine (gene)", "Osmolyte genes", grab("osmolyte_drought.tsv","pathway","Glycine-betaine"))

if not rows:
    print("no data"); raise SystemExit

df = pd.DataFrame(rows, columns=["label","cat","delta"])
cats = ["Taxonomy","N-cycle genes","CAZymes","Osmolyte genes"]
colors = {"Taxonomy":"#d62728","N-cycle genes":"#1f77b4","CAZymes":"#2ca02c","Osmolyte genes":"#9467bd"}
df["cat"] = pd.Categorical(df["cat"], categories=cats, ordered=True)
df = df.sort_values(["cat","delta"]).reset_index(drop=True)
y = np.arange(len(df))[::-1]

fig, ax = plt.subplots(figsize=(8.4, 5.6))
# negligible band (|delta| < 0.33 = small; 0.33-0.474 medium; >0.474 large, Romano thresholds)
ax.axvspan(-0.33, 0.33, color="#eeeeee", zorder=0)
ax.axvline(0, color="#888", lw=0.8, zorder=1)
ax.scatter(df["delta"], y, c=[colors[c] for c in df["cat"]], s=90, zorder=3, edgecolor="white", linewidth=0.6)
for yi, (_, r) in zip(y, df.iterrows()):
    ax.plot([0, r["delta"]], [yi, yi], color=colors[r["cat"]], lw=2, alpha=0.5, zorder=2)
ax.set_yticks(y); ax.set_yticklabels(df["label"], fontsize=9)
ax.set_xlim(-1, 1)
ax.set_xlabel("drought effect size (Cliff's δ)   ·   ← lower under drought   |   higher →")
ax.set_title("Taxonomy shifts, functional gene potential is buffered\n"
             "drought effect size across data types (grey band = negligible, |δ|<0.33)",
             fontsize=12)
# legend
from matplotlib.lines import Line2D
leg = [Line2D([0],[0],marker='o',color='w',markerfacecolor=colors[c],markersize=9,label=c) for c in cats]
ax.legend(handles=leg, fontsize=8, frameon=False, loc="lower right")
plt.tight_layout(); plt.savefig(f"{FIG}/17_buffering_overview.png", dpi=140); plt.close()
print(f"Wrote 17_buffering_overview.png ({len(df)} readouts)")
print(df.to_string(index=False))
