#!/usr/bin/env python3
"""
TIER 1 — compositional & robustness stats (report review #2,#3,#4).
  (3) compositional artifact -> CLR/Aitchison ordination+PERMANOVA; CLR Welch+BH;
                                skbio dirmult_ttest (Dirichlet-multinomial) and ancombc.
  (2) pseudoreplication      -> skbio dirmult_lme with PLOT as random effect, plus a
                                plot-collapsed MWU sensitivity check.
  (4) power                  -> min detectable Cohen's d at n=17 vs 19.
Outputs results/*.tsv|txt + figs/10_aitchison_pcoa.png ; self-checks at end.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
from skbio.stats.composition import clr, multi_replace, ancombc, dirmult_ttest, dirmult_lme
from skbio.stats.distance import permanova, permdisp, DistanceMatrix
from skbio.stats.ordination import pcoa

A = "/mnt/disk4/timo/gbi/analysis"; FIG = f"{A}/figs"; RES = f"{A}/results"
meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
genus = pd.read_csv(f"{A}/genus_relabund.tsv", sep="\t", index_col=0)
phylum = pd.read_csv(f"{A}/phylum_relabund.tsv", sep="\t", index_col=0)
samples = [s for s in genus.index if s in meta.index]
genus = genus.loc[samples]; phylum = phylum.loc[samples]; md = meta.loc[samples].copy()
md["timepoint"] = md["timepoint"].astype(str); md["plot"] = md["plot"].astype(str)
cond = md["condition"].values
PRE = md.index[md.condition == "pre-drought"]; DRO = md.index[md.condition == "drought"]
pal = {"pre-drought": "#1f77b4", "drought": "#d62728"}
checks = []
def bh(p):
    p=np.asarray(p,float); o=np.argsort(p); r=p[o]; m=len(p)
    q=np.minimum.accumulate((r*m/np.arange(1,m+1))[::-1])[::-1]
    out=np.empty(m); out[o]=np.clip(q,0,1); return out

# ---------- CLR + Aitchison ordination ----------
def clr_frame(df):
    M = df.values.astype(float); M = M/M.sum(1,keepdims=True); M = multi_replace(M)
    return pd.DataFrame(clr(M), index=df.index, columns=df.columns)
clr_g = clr_frame(genus)
ait = np.sqrt(((clr_g.values[:,None,:]-clr_g.values[None,:,:])**2).sum(-1))
dm = DistanceMatrix(ait, ids=list(samples))
pn = permanova(dm, grouping=list(cond), permutations=999)
pdsp = permdisp(dm, grouping=list(cond), permutations=999)
ordr = pcoa(dm); pc = ordr.samples.iloc[:,:2].copy(); pc.columns = [0,1]; pc.index = samples
ve = ordr.proportion_explained[:2].values*100
with open(f"{RES}/aitchison_stats.txt","w") as f:
    f.write(f"PERMANOVA_Aitchison_F\t{pn['test statistic']:.4f}\nPERMANOVA_Aitchison_p\t{pn['p-value']:.4g}\n")
    f.write(f"PERMDISP_Aitchison_F\t{pdsp['test statistic']:.4f}\nPERMDISP_Aitchison_p\t{pdsp['p-value']:.4g}\n")
print(f"Aitchison PERMANOVA p={pn['p-value']:.3g} | PERMDISP p={pdsp['p-value']:.3g}")
tpm={"0h":"o","6h":"^","48h":"s"}
fig,ax=plt.subplots(figsize=(5.6,4.6))
for c in ["pre-drought","drought"]:
    for tp in ["0h","6h","48h"]:
        idx=[s for s in samples if md.loc[s,"condition"]==c and md.loc[s,"timepoint"]==tp]
        if idx: ax.scatter(pc.loc[idx,0],pc.loc[idx,1],color=pal[c],marker=tpm[tp],s=50,alpha=.85,edgecolor="k",linewidth=.4)
ax.set_xlabel(f"PCo1 ({ve[0]:.1f}%)"); ax.set_ylabel(f"PCo2 ({ve[1]:.1f}%)")
ax.set_title(f"Aitchison (CLR) PCoA\nPERMANOVA p={pn['p-value']:.3g} | PERMDISP p={pdsp['p-value']:.3g}")
plt.tight_layout(); plt.savefig(f"{FIG}/10_aitchison_pcoa.png",dpi=130); plt.close()

# ---------- CLR Welch + BH ----------
rows=[]
for g in clr_g.columns:
    a=clr_g.loc[PRE,g]; b=clr_g.loc[DRO,g]; t,p=stats.ttest_ind(a,b,equal_var=False)
    rows.append({"genus":g,"clr_diff":b.mean()-a.mean(),"welch_p":p})
clrda=pd.DataFrame(rows); clrda["q_BH"]=bh(clrda["welch_p"].values); clrda=clrda.sort_values("welch_p")
clrda.to_csv(f"{RES}/clr_diff_abundance.tsv",sep="\t",index=False)
nit_clr=clrda[clrda.genus.str.startswith("Nitrospira")]
print("\nCLR Welch+BH top6:"); print(clrda.head(6).to_string(index=False))
print("Nitrospira (CLR):"); print(nit_clr.to_string(index=False))

# ---------- dirmult_ttest (Dirichlet-multinomial compositional) ----------
grp = pd.Series(cond, index=samples)
def save_df(obj, name):
    try:
        df = obj if isinstance(obj, pd.DataFrame) else pd.DataFrame(obj)
        df.to_csv(f"{RES}/{name}", sep="\t"); return df
    except Exception as e:
        print(f"  could not save {name}: {e}"); return None
try:
    dmt = dirmult_ttest(genus, grp, treatment="drought", reference="pre-drought", seed=42)
    save_df(dmt, "dirmult_ttest_genus.tsv")
    sig_dmt = dmt[dmt["Signif"]] if "Signif" in dmt.columns else dmt[dmt.get("qvalue",1)<0.05]
    print(f"\ndirmult_ttest: {len(sig_dmt)} significant genera")
    print(sig_dmt.head(10).to_string()); checks.append(("dirmult_ttest", True))
except Exception as e:
    print("dirmult_ttest failed:", e); checks.append(("dirmult_ttest", False))

# ---------- dirmult_lme: PLOT random effect (pseudoreplication fix) ----------
try:
    md_l = md.copy()
    lme = dirmult_lme(genus, md_l, formula="condition", grouping="plot",
                      draws=64, seed=42, re_formula="1")
    save_df(lme, "dirmult_lme_plot_genus.tsv")
    sig_lme = lme[lme["Signif"]] if "Signif" in lme.columns else lme[lme.get("qvalue",1)<0.05]
    print(f"\ndirmult_lme (plot RE): {len(sig_lme)} significant genera")
    print(sig_lme.head(10).to_string()); checks.append(("dirmult_lme", True))
except Exception as e:
    print("dirmult_lme failed:", e); checks.append(("dirmult_lme", False))

# ---------- ancombc ----------
try:
    tbl_pos = pd.DataFrame(multi_replace(genus.values/genus.values.sum(1,keepdims=True)),
                           index=samples, columns=genus.columns)   # zero-free for ANCOM-BC
    abc = ancombc(tbl_pos, md, formula="condition")
    save_df(abc if isinstance(abc,pd.DataFrame) else pd.DataFrame(abc), "ancombc_genus.tsv")
    print("\nancombc ran ->", type(abc).__name__); checks.append(("ancombc", True))
except Exception as e:
    print("ancombc failed:", e); checks.append(("ancombc", False))

# ---------- plot-collapsed sensitivity ----------
gp = genus.groupby(md["plot"]).mean(); pc_cond = md.groupby("plot")["condition"].first()
pP=gp.index[pc_cond[gp.index]=="pre-drought"]; dP=gp.index[pc_cond[gp.index]=="drought"]
rows=[]
for g in ["Nitrospira_D","Bradyrhizobium","Methyloceanibacter"]:
    if g not in gp.columns: continue
    a=gp.loc[pP,g]; b=gp.loc[dP,g]
    try: _,p=stats.mannwhitneyu(a,b,alternative="two-sided")
    except ValueError: p=np.nan
    rows.append({"genus":g,"n_pre_plots":len(pP),"n_dro_plots":len(dP),"pre":a.mean(),"drought":b.mean(),"MWU_p_plot":p})
plotdf=pd.DataFrame(rows); plotdf.to_csv(f"{RES}/plot_collapsed_drought.tsv",sep="\t",index=False)
print(f"\nPlot-collapsed ({len(pP)} pre vs {len(dP)} drought plots):"); print(plotdf.to_string(index=False))

# ---------- power ----------
try:
    from statsmodels.stats.power import TTestIndPower
    an=TTestIndPower(); d80=an.solve_power(effect_size=None,nobs1=17,ratio=19/17,alpha=0.05,power=0.8)
    p08=an.power(effect_size=0.8,nobs1=17,ratio=19/17,alpha=0.05)
    open(f"{RES}/power_analysis.txt","w").write(
        f"n_pre\t17\nn_drought\t19\nmin_detectable_Cohens_d_80pct\t{d80:.3f}\npower_at_d0.8\t{p08:.3f}\n")
    print(f"\nPower: min detectable d={d80:.2f} (80%); power@d=0.8 = {p08:.2f}"); checks.append(("power",True))
except Exception as e:
    print("power failed:",e); checks.append(("power",False))

# ---------- CHECKS ----------
print("\n=== CHECKS ===")
for o in ["aitchison_stats.txt","clr_diff_abundance.tsv","plot_collapsed_drought.tsv","power_analysis.txt"]:
    checks.append((o, os.path.exists(f"{RES}/{o}") and os.path.getsize(f"{RES}/{o}")>0))
checks.append(("aitchison_p in (0,1]", 0<pn['p-value']<=1))
checks.append(("CLR Nitrospira down", bool(len(nit_clr) and (nit_clr["clr_diff"]<0).any())))
checks.append(("fig10 written", os.path.exists(f"{FIG}/10_aitchison_pcoa.png")))
for k,v in checks: print(f"  [{'OK' if v else 'FAIL'}] {k}")
open(f"{RES}/tier1_checks.txt","w").write("\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks))
print(f"\nTIER1 {'ALL PASS' if all(v for _,v in checks) else 'SOME FAIL (see above)'}")
