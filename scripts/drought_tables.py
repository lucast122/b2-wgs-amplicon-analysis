#!/usr/bin/env python3
"""
Regenerate the per-taxon drought differential tables from the (n=39) sylph matrices,
replacing lost ad-hoc producers. Outputs columns matching what build_report.py reads;
add_fdr_effectsizes.py then appends q_BH + Cliff's delta.
  genus_drought.tsv   : genus, pre, drought, log2FC, p, named
  species_drought.tsv : species, pre, drought, log2FC, p
  uncultured_fraction.txt : mean_uncultured_pct, min, max  (over s__ species per sample)
Self-checks -> results/drought_tables_checks.txt
"""
import re, numpy as np, pandas as pd
from scipy import stats
A="/mnt/disk4/timo/gbi/analysis"; RES=f"{A}/results"
meta=pd.read_csv(f"{A}/metadata.tsv",sep="\t",index_col=0)
PRE=set(meta.index[meta.condition=="pre-drought"]); DRO=set(meta.index[meta.condition=="drought"])

def diff_table(mat, namecol):
    pre=[s for s in mat.index if s in PRE]; dro=[s for s in mat.index if s in DRO]
    rows=[]
    for t in mat.columns:
        a=mat.loc[pre,t].values; b=mat.loc[dro,t].values
        if a.mean()<1e-4 and b.mean()<1e-4: continue
        try: _,p=stats.mannwhitneyu(a,b,alternative="two-sided")
        except ValueError: p=np.nan
        rows.append({namecol:t,"pre":a.mean(),"drought":b.mean(),
                     "log2FC":np.log2((b.mean()+1e-3)/(a.mean()+1e-3)),"p":p})
    return pd.DataFrame(rows).sort_values("p")

# genus
g=pd.read_csv(f"{A}/genus_relabund.tsv",sep="\t",index_col=0)
gd=diff_table(g,"genus")
gd["named"]=gd["genus"].apply(lambda s: bool(re.search(r"[a-z]", str(s))))  # GTDB codes are caps+digits
gd.to_csv(f"{RES}/genus_drought.tsv",sep="\t",index=False)

# species
sp=pd.read_csv(f"{A}/species_relabund.tsv",sep="\t",index_col=0)
sd=diff_table(sp,"species")
sd.to_csv(f"{RES}/species_drought.tsv",sep="\t",index=False)

# uncultured fraction (per sample, over species)
def is_uncultured(s):
    m=re.match(r"^[A-Z][a-z]+(?:_[A-Z])? ([a-z]+)", str(s))   # Genus epithet
    return not (m and not re.match(r"sp\d", m.group(1)))      # cultured if real epithet (not sp123)
unc_cols=[c for c in sp.columns if is_uncultured(c)]
tot=sp.sum(1).replace(0,np.nan)
unc_pct=(sp[unc_cols].sum(1)/tot*100).dropna()
with open(f"{RES}/uncultured_fraction.txt","w") as f:
    f.write(f"mean_uncultured_pct\t{unc_pct.mean():.4f}\nmin\t{unc_pct.min():.4f}\nmax\t{unc_pct.max():.4f}\n")

# checks
nit_g = gd[gd.genus.str.startswith("Nitrospira")]
nit_s = sd[sd.species.str.startswith("Nitrospira_D sp")]
checks=[
 ("genus_drought rows>0", len(gd)>0),
 ("species_drought rows>0", len(sd)>0),
 ("Nitrospira genus present", len(nit_g)>0),
 ("Nitrospira species present", len(nit_s)>0),
 ("Nitrospira_D down in drought", bool(len(nit_g) and (nit_g[nit_g.genus=='Nitrospira_D']['drought'].values[0] < nit_g[nit_g.genus=='Nitrospira_D']['pre'].values[0]) if 'Nitrospira_D' in set(gd.genus) else False)),
 ("uncultured mean in (0,100]", 0<unc_pct.mean()<=100),
]
open(f"{RES}/drought_tables_checks.txt","w").write("\n".join(f"{'OK' if v else 'FAIL'}\t{k}" for k,v in checks))
print(f"genus_drought: {len(gd)} genera; species_drought: {len(sd)}; uncultured mean {unc_pct.mean():.1f}% (range {unc_pct.min():.0f}-{unc_pct.max():.0f})")
print("Nitrospira_D genus:", nit_g[nit_g.genus=='Nitrospira_D'][['pre','drought','p']].to_string(index=False) if 'Nitrospira_D' in set(gd.genus) else "absent")
for k,v in checks: print(f"  [{'OK' if v else 'FAIL'}] {k}")
print("ALL PASS" if all(v for _,v in checks) else "SOME FAIL")
