#!/usr/bin/env python3
"""Regenerate the genus-level (08) and Nitrospira-across-ranks (09) drought figures with
proper axis headroom so p-value labels no longer overflow the plot border. These two figures
had lost their generator in an earlier refactor; this is now their canonical source."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"

# ---- 08: genus-level drought response (named genera) ----
g = pd.read_csv(f"{RES}/genus_drought.tsv", sep="\t")
if "named" in g.columns:
    g = g[g["named"] == True]
g = g.sort_values("p").head(8).reset_index(drop=True)
y = np.arange(len(g))[::-1]                       # smallest-p at top
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.barh(y + 0.2, g["pre"], 0.4, label="pre-drought", color="#1f77b4")
ax.barh(y - 0.2, g["drought"], 0.4, label="drought", color="#d62728")
xmax = float(max(g["pre"].max(), g["drought"].max()))
ax.set_xlim(0, xmax * 1.34)                       # room for p-labels on the right
for yi, (_, r) in zip(y, g.iterrows()):
    star = "*" if r["p"] < 0.05 else ""
    ax.text(max(r["pre"], r["drought"]) + xmax * 0.015, yi,
            f"p={r['p']:.2g}{star}", va="center", fontsize=7.5)
ax.set_yticks(y); ax.set_yticklabels(g["genus"], fontsize=8)
ax.set_xlabel("mean relative abundance (%)")
ax.set_title("Genus-level drought response (named genera, * p<0.05)")
ax.legend(frameon=False, loc="lower right")
plt.tight_layout(); plt.savefig(f"{FIG}/08_genus_drought.png", dpi=130); plt.close()

# ---- 09: Nitrospira consistent decline across ranks ----
rows = []
try:
    ph = pd.read_csv(f"{RES}/phylum_shift_pre_vs_drought.tsv", sep="\t")
    r = ph[ph.phylum.str.contains("Nitrospirota", case=False, na=False)]
    if len(r): rows.append(("Phylum  Nitrospirota", r.iloc[0]["pre_mean"], r.iloc[0]["drought_mean"]))
except Exception: pass
try:
    gd = pd.read_csv(f"{RES}/genus_drought.tsv", sep="\t")
    for _, r in gd[gd.genus.str.startswith("Nitrospira")].iterrows():
        rows.append((f"Genus  {r.genus}", r["pre"], r["drought"]))
except Exception: pass
try:
    sp = pd.read_csv(f"{RES}/species_drought.tsv", sep="\t")
    sp = sp[sp.species.str.contains("Nitrospira", na=False)].sort_values("pre", ascending=False).head(2)
    for _, r in sp.iterrows():
        nm = r.species if len(r.species) < 26 else r.species[:24] + "…"
        rows.append((f"Species  {nm}", r["pre"], r["drought"]))
except Exception: pass

if rows:
    lab = [x[0] for x in rows]; pre = [x[1] for x in rows]; dro = [x[2] for x in rows]
    y = np.arange(len(rows))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.barh(y + 0.2, pre, 0.4, label="pre-drought", color="#1f77b4")
    ax.barh(y - 0.2, dro, 0.4, label="drought", color="#d62728")
    xmax = float(max(max(pre), max(dro)))
    ax.set_xlim(0, xmax * 1.18)
    for yi, p, d in zip(y, pre, dro):
        ax.text(max(p, d) + xmax * 0.015, yi, f"{(d-p)/ (p+1e-9)*100:+.0f}%", va="center", fontsize=7.5, color="#555")
    ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=8)
    ax.set_xlabel("mean relative abundance (%)")
    ax.set_title("Nitrospira declines consistently across taxonomic ranks")
    ax.legend(frameon=False, loc="lower right")
    plt.tight_layout(); plt.savefig(f"{FIG}/09_nitrospira_ranks.png", dpi=130); plt.close()

print(f"Wrote 08_genus_drought.png ({len(g)} genera) and 09_nitrospira_ranks.png ({len(rows)} rank entries)")
