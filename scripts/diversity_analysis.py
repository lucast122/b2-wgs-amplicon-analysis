#!/usr/bin/env python3
"""
Core ecology analysis for the GBI Biosphere 2 drought metagenomes (sylph).
Ties directly to the paper's central claim: drought REDUCES site-to-site
variation in microbial community structure (convergence).

Analyses:
  1. Alpha diversity (Shannon) per sample; pre-drought vs drought
  2. Bray-Curtis PCoA ordination (condition, site, timepoint)
  3. PERMANOVA: does condition structure the community?
  4. PERMDISP (betadisper): is multivariate dispersion LOWER under drought?
     -> direct test of the paper's "convergence" hypothesis
  5. Bacteria:Fungi and phylum-level shifts pre-drought -> drought

Outputs to /mnt/disk4/timo/gbi/analysis/figs and results/*.tsv
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from skbio.stats.distance import permanova, permdisp, DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity

A = "/mnt/disk4/timo/gbi/analysis"
FIG = f"{A}/figs"; RES = f"{A}/results"
os.makedirs(FIG, exist_ok=True); os.makedirs(RES, exist_ok=True)

meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
genus = pd.read_csv(f"{A}/genus_relabund.tsv", sep="\t", index_col=0)
phylum = pd.read_csv(f"{A}/phylum_relabund.tsv", sep="\t", index_col=0)

# align
samples = [s for s in genus.index if s in meta.index]
genus = genus.loc[samples]; phylum = phylum.loc[samples]; md = meta.loc[samples]
cond = md["condition"].values
print(f"{len(samples)} samples | conditions: {pd.Series(cond).value_counts().to_dict()}")

# palette
pal = {"pre-drought": "#1f77b4", "drought": "#d62728"}

# ---- 1. Alpha diversity (Shannon on genus relabund) -----------------------
def shannon(row):
    p = row[row > 0] / 100.0
    return float(-(p * np.log(p)).sum())
md = md.copy()
md["shannon"] = [shannon(genus.loc[s]) for s in samples]
md["richness"] = (genus > 0).sum(axis=1).values

pre = md.loc[md.condition == "pre-drought", "shannon"]
dro = md.loc[md.condition == "drought", "shannon"]
u, p_alpha = stats.mannwhitneyu(pre, dro, alternative="two-sided")
print(f"Shannon  pre={pre.mean():.3f}  drought={dro.mean():.3f}  MWU p={p_alpha:.3g}")

fig, ax = plt.subplots(figsize=(4, 4))
for i, c in enumerate(["pre-drought", "drought"]):
    vals = md.loc[md.condition == c, "shannon"]
    ax.scatter(np.random.normal(i, 0.06, len(vals)), vals, color=pal[c], alpha=.7, s=30)
    ax.bar(i, vals.mean(), width=.5, color=pal[c], alpha=.25)
ax.set_xticks([0, 1]); ax.set_xticklabels(["pre-drought", "drought"])
ax.set_ylabel("Shannon diversity (genus)")
ax.set_title(f"Alpha diversity  (MWU p={p_alpha:.2g})")
plt.tight_layout(); plt.savefig(f"{FIG}/01_alpha_shannon.png", dpi=130); plt.close()

# ---- 2-4. Beta diversity, PCoA, PERMANOVA, PERMDISP -----------------------
# Bray-Curtis on genus relative abundances
bc = beta_diversity("braycurtis", genus.values, ids=samples)
dm = DistanceMatrix(bc.data, ids=samples)

pn = permanova(dm, grouping=list(cond), permutations=999)
pd_disp = permdisp(dm, grouping=list(cond), permutations=999)
print(f"PERMANOVA  F={pn['test statistic']:.2f}  p={pn['p-value']:.3g}")
print(f"PERMDISP   F={pd_disp['test statistic']:.2f}  p={pd_disp['p-value']:.3g}")

ord_res = pcoa(dm)
pc = ord_res.samples.iloc[:, :2].copy()
pc.columns = [0, 1]; pc.index = samples
ve = ord_res.proportion_explained[:2].values * 100

fig, ax = plt.subplots(figsize=(5, 4.4))
for c in ["pre-drought", "drought"]:
    idx = [s for s in samples if md.loc[s, "condition"] == c]
    ax.scatter(pc.loc[idx, 0], pc.loc[idx, 1], color=pal[c], label=c, alpha=.8, s=45,
               edgecolor="k", linewidth=.4)
ax.set_xlabel(f"PCo1 ({ve[0]:.1f}%)"); ax.set_ylabel(f"PCo2 ({ve[1]:.1f}%)")
ax.set_title(f"Bray-Curtis PCoA\nPERMANOVA p={pn['p-value']:.3g} | PERMDISP p={pd_disp['p-value']:.3g}")
ax.legend(frameon=False)
plt.tight_layout(); plt.savefig(f"{FIG}/02_pcoa_braycurtis.png", dpi=130); plt.close()

# ---- Convergence: mean distance-to-centroid per condition -----------------
# (skbio permdisp uses median; also compute mean within-group BC for clarity)
def within_group_dist(ids):
    sub = dm.filter(ids).data
    iu = np.triu_indices(len(ids), 1)
    return sub[iu]
wd_pre = within_group_dist([s for s in samples if md.loc[s, "condition"] == "pre-drought"])
wd_dro = within_group_dist([s for s in samples if md.loc[s, "condition"] == "drought"])
u2, p_disp_mwu = stats.mannwhitneyu(wd_pre, wd_dro, alternative="greater")  # pre>drought => convergence
print(f"Within-group BC dissimilarity  pre={wd_pre.mean():.3f}  drought={wd_dro.mean():.3f}")
print(f"  (convergence test, pre>drought) MWU p={p_disp_mwu:.3g}")

fig, ax = plt.subplots(figsize=(4, 4))
ax.boxplot([wd_pre, wd_dro], labels=["pre-drought", "drought"])
ax.set_ylabel("within-condition Bray-Curtis dissimilarity")
ax.set_title(f"Community heterogeneity\n(lower = convergence)  MWU p={p_disp_mwu:.2g}")
plt.tight_layout(); plt.savefig(f"{FIG}/03_convergence.png", dpi=130); plt.close()

# ---- 5. Phylum shifts pre-drought -> drought ------------------------------
rows = []
for ph in phylum.columns:
    a = phylum.loc[md.condition == "pre-drought", ph]
    b = phylum.loc[md.condition == "drought", ph]
    if a.mean() < 0.1 and b.mean() < 0.1:
        continue
    try:
        _, pv = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        pv = np.nan
    rows.append({"phylum": ph, "pre_mean": a.mean(), "drought_mean": b.mean(),
                 "log2FC": np.log2((b.mean()+1e-3)/(a.mean()+1e-3)), "MWU_p": pv})
shift = pd.DataFrame(rows).sort_values("MWU_p")
shift.to_csv(f"{RES}/phylum_shift_pre_vs_drought.tsv", sep="\t", index=False)
print("\nTop phylum shifts (pre -> drought):")
print(shift.head(8).to_string(index=False))

md[["condition","site","timepoint","shannon","richness"]].to_csv(f"{RES}/alpha_diversity.tsv", sep="\t")
with open(f"{RES}/stats_summary.txt", "w") as f:
    f.write(f"n_samples\t{len(samples)}\n")
    f.write(f"shannon_pre\t{pre.mean():.4f}\nshannon_drought\t{dro.mean():.4f}\nshannon_MWU_p\t{p_alpha:.4g}\n")
    f.write(f"PERMANOVA_F\t{pn['test statistic']:.4f}\nPERMANOVA_p\t{pn['p-value']:.4g}\n")
    f.write(f"PERMDISP_F\t{pd_disp['test statistic']:.4f}\nPERMDISP_p\t{pd_disp['p-value']:.4g}\n")
    f.write(f"withinBC_pre\t{wd_pre.mean():.4f}\nwithinBC_drought\t{wd_dro.mean():.4f}\nconvergence_MWU_p\t{p_disp_mwu:.4g}\n")
print("\nDONE diversity_analysis — figs in figs/, results in results/")
