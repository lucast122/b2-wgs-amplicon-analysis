#!/usr/bin/env python3
"""
13C integration (Tier-4 capstone): pair each microbial group's DROUGHT RESPONSE measured two ways
 - ACTIVITY: delta-13C uptake of recent plant carbon (PLFA-SIP, companion B2 data), dry vs ambient
 - ABUNDANCE: standing relative abundance (this WGS metagenome), drought vs pre-drought
A per-plot join is not possible (the WGS and PLFA samplings use different plot/rep schemes &
timepoint scales), so groups are compared at the group level — which is exactly the axis the
report's thesis turns on: standing abundance vs carbon-uptake activity.

Outputs: b2/analysis/results/c13_group_comparison.tsv, c13_checks.txt ; figs/18_c13_vs_abundance.png
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"; FIG = f"{A}/figs"
XLS = "/mnt/disk4/timo/gbi/b2/c13/B2_plot_envmicrobe_data.xlsx"
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# ---------- 13C activity response (dry - amb, post-label) per PLFA group ----------
sm = pd.read_excel(XLS, sheet_name="special microbe")
post = sm[sm.tim > 0]
def c13_resp(g):
    a = post[post.tre == "amb"][g].dropna(); d = post[post.tre == "dry"][g].dropna()
    return d.mean() - a.mean()
GROUPS = {"AC": "Actinobacteria", "AMF": "Arbuscular mycorrhiza", "GN": "Gram-negative",
          "GP": "Gram-positive", "SF": "Saprotrophic fungi", "bac": "Total bacteria"}

# ---------- WGS abundance response (log2 drought/pre) per matching group ----------
m = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
phy = pd.read_csv(f"{A}/phylum_relabund.tsv", sep="\t", index_col=0)
s = [x for x in phy.index if x in m.index]; phy = phy.loc[s]; md = m.loc[s]
pre = md.index[md.condition == "pre-drought"]; dro = md.index[md.condition == "drought"]
def phyl(prefixes):
    cc = [c for c in phy.columns if any(c.startswith(p) for p in prefixes)]
    v = phy[cc].sum(1); return v.loc[pre].mean(), v.loc[dro].mean()
GPOS = ["Actinomycetota","Actinobacteriota","Firmicutes","Bacillota"]
GNEG = ["Pseudomonadota","Proteobacteria","Bacteroidota","Acidobacteriota","Verrucomicrobiota","Gemmatimonadota","Myxococcota"]
amf = pd.read_csv(f"{RES}/kaiju_amf_fungalnorm.tsv", sep="\t", index_col=0)
dom = pd.read_csv(f"{RES}/kaiju_domain_fractions.tsv", sep="\t", index_col=0)
sa = [x for x in amf.index if x in m.index]; mda = m.loc[sa]
pa = mda.index[mda.condition=="pre-drought"]; da = mda.index[mda.condition=="drought"]
def wgs_abund(code):
    if code == "AC":  return phyl(["Actinomycetota","Actinobacteriota"])
    if code == "GP":  return phyl(GPOS)
    if code == "GN":  return phyl(GNEG)
    if code == "bac": p,d = phyl(GPOS+GNEG+["Chloroflexota","Planctomycetota","Nitrospirota"]); return p,d
    if code == "AMF": return amf.loc[pa].sum(1).mean(), amf.loc[da].sum(1).mean()
    if code == "SF":  return dom.loc[pa,"Fungi_other"].mean(), dom.loc[da,"Fungi_other"].mean()
    return (np.nan, np.nan)

rows = []
for code, name in GROUPS.items():
    c13 = c13_resp(code)
    ap, ad = wgs_abund(code)
    log2fc = np.log2((ad + 1e-3) / (ap + 1e-3))
    rows.append({"group": code, "name": name, "c13_dry_minus_amb": c13,
                 "wgs_pre": ap, "wgs_drought": ad, "wgs_log2FC": log2fc})
df = pd.DataFrame(rows)
df.round(3).to_csv(f"{RES}/c13_group_comparison.tsv", sep="\t", index=False)
print(df.round(2).to_string(index=False))

# correlation across groups (do activity and abundance responses track each other?)
from scipy import stats
rho, p = stats.spearmanr(df["c13_dry_minus_amb"], df["wgs_log2FC"])
print(f"\nSpearman(13C-activity response, WGS-abundance response) across {len(df)} groups: rho={rho:.2f}, p={p:.2f}")

# ---------- figure: activity vs abundance response, per group ----------
fig_ok = True
try:
    plt.figure(figsize=(7.6, 6))
    plt.axhline(0, color="#bbb", lw=.8); plt.axvline(0, color="#bbb", lw=.8)
    # quadrant shading for the decoupled corners
    plt.gca().add_patch(plt.Rectangle((-3,0),3,3,color="#fde6e6",zorder=0))   # activity-down, abundance-up
    plt.gca().add_patch(plt.Rectangle((0,-3),3,3,color="#fde6e6",zorder=0))   # activity-up, abundance-down
    cols = {"AMF":"#9467bd","SF":"#2ca02c","AC":"#d62728","GP":"#1f77b4","GN":"#17becf","bac":"#7f7f7f"}
    for _, r in df.iterrows():
        plt.scatter(r["c13_dry_minus_amb"], r["wgs_log2FC"], s=130, color=cols.get(r["group"],"#333"),
                    edgecolor="white", zorder=3)
        plt.annotate(f" {r['name']}", (r["c13_dry_minus_amb"], r["wgs_log2FC"]), fontsize=8.5, va="center")
    plt.xlabel("← less   ¹³C carbon-uptake response to drought (Δδ¹³C, dry−amb)   more →")
    plt.ylabel("← lower   WGS abundance response (log₂ drought/pre)   higher →")
    plt.title("Activity vs abundance under drought — two different axes\n"
              "(pink corners = decoupled: abundance and carbon-uptake move oppositely)", fontsize=11)
    plt.xlim(-2.2, 3); plt.ylim(-0.9, 1.1)
    plt.tight_layout(); plt.savefig(f"{FIG}/18_c13_vs_abundance.png", dpi=140); plt.close()
except Exception as e:
    print("fig failed:", e); fig_ok = False

checks = [
    ("6 groups paired", len(df) == 6),
    ("AMF 13C-activity up", float(df[df.group=="AMF"]["c13_dry_minus_amb"]) > 0),
    ("Actino 13C-activity down but abundance up (decoupled)",
        float(df[df.group=="AC"]["c13_dry_minus_amb"]) < 0 and float(df[df.group=="AC"]["wgs_log2FC"]) > 0),
    ("figure written", fig_ok),
]
open(f"{RES}/c13_checks.txt","w").write(
    "\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks) +
    f"\nspearman_activity_vs_abundance\t{rho:.3f}\nspearman_p\t{p:.3f}\n")
print("\n=== CHECKS ==="); [print(f"  [{'OK' if v else 'FAIL'}] {k}") for k,v in checks]
print(f"  {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL'}")
