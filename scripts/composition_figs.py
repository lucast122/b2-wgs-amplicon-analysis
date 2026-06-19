#!/usr/bin/env python3
"""Community composition overview + drought-response figure for the report."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

A = "/mnt/disk4/timo/gbi/b2/analysis"; FIG = f"{A}/figs"
meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
phylum = pd.read_csv(f"{A}/phylum_relabund.tsv", sep="\t", index_col=0)
samples = [s for s in phylum.index if s in meta.index]
phylum = phylum.loc[samples]; md = meta.loc[samples]

# ---- mean composition by condition (top phyla) ----
top = phylum.mean().sort_values(ascending=False).head(10).index.tolist()
comp = []
for c in ["pre-drought", "drought"]:
    m = phylum.loc[md.condition == c, top].mean()
    comp.append(m)
comp = pd.DataFrame(comp, index=["pre-drought", "drought"])
comp["Other"] = 100 - comp.sum(axis=1)

colors = plt.cm.tab20(np.linspace(0, 1, len(comp.columns)))
fig, ax = plt.subplots(figsize=(6, 4.2))
left = np.zeros(2)
for j, ph in enumerate(comp.columns):
    ax.barh(["pre-drought", "drought"], comp[ph], left=left, color=colors[j], label=ph)
    left += comp[ph].values
ax.set_xlabel("mean relative abundance (%)")
ax.set_title("Soil community composition (phylum, sylph)")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7, frameon=False)
ax.set_xlim(0, 100)
plt.tight_layout(); plt.savefig(f"{FIG}/04_composition.png", dpi=130, bbox_inches="tight"); plt.close()

# ---- drought-response: significant phyla ----
from scipy import stats
sig = []
for ph in phylum.columns:
    a = phylum.loc[md.condition == "pre-drought", ph]
    b = phylum.loc[md.condition == "drought", ph]
    if a.mean() < 0.1 and b.mean() < 0.1: continue
    try: _, pv = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError: pv = 1.0
    sig.append((ph, a.mean(), b.mean(), pv))
sig = pd.DataFrame(sig, columns=["phylum","pre","drought","p"]).sort_values("p")
show = sig.head(6)

fig, ax = plt.subplots(figsize=(6, 4))
x = np.arange(len(show)); w = 0.38
ax.bar(x - w/2, show["pre"], w, label="pre-drought", color="#1f77b4")
ax.bar(x + w/2, show["drought"], w, label="drought", color="#d62728")
for i, (_, r) in enumerate(show.iterrows()):
    star = "*" if r["p"] < 0.05 else ""
    ax.text(i, max(r["pre"], r["drought"]) + 0.2, f"p={r['p']:.2g}{star}",
            ha="center", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(show["phylum"], rotation=35, ha="right", fontsize=8)
ax.set_ylabel("mean relative abundance (%)")
ax.set_title("Phylum-level drought response (* p<0.05)")
ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(f"{FIG}/05_drought_response.png", dpi=130); plt.close()
print("Wrote 04_composition.png, 05_drought_response.png")
print("\nMean composition (%):")
print(comp.round(2).to_string())
