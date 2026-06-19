#!/usr/bin/env python3
"""
Post-hoc rigor pass (report review #1, #6): add Benjamini-Hochberg FDR q-values
and Cliff's delta effect sizes to every differential-abundance table, in place.

Reads the existing *_drought.tsv tables + their relabund matrices, adds:
  q_BH         - BH-FDR over the p-values WITHIN that table
  cliffs_delta - non-parametric effect size, drought vs pre (+ = higher in drought)
  cliffs_mag   - negligible/small/medium/large (Romano thresholds)
Rewrites each TSV with the new columns appended. Safe to re-run.
"""
import numpy as np, pandas as pd, os

A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"
meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
PRE  = set(meta.index[meta.condition == "pre-drought"])
DRO  = set(meta.index[meta.condition == "drought"])

def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.full(n, np.nan)
    ok = ~np.isnan(p); idx = np.where(ok)[0]
    if len(idx) == 0: return q
    pv = p[idx]; order = np.argsort(pv); ranked = pv[order]; m = len(pv)
    qv = ranked * m / (np.arange(1, m+1))
    qv = np.minimum.accumulate(qv[::-1])[::-1]      # enforce monotonicity
    out = np.empty(m); out[order] = np.clip(qv, 0, 1)
    q[idx] = out; return q

def cliffs(a, b):
    """delta in [-1,1]; + means b (drought) tends higher than a (pre)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0: return np.nan
    gt = sum((bj > a).sum() for bj in b); lt = sum((bj < a).sum() for bj in b)
    return (gt - lt) / (len(a) * len(b))

def mag(d):
    d = abs(d)
    return "negligible" if d < .147 else "small" if d < .33 else "medium" if d < .474 else "large"

def annotate(tsv, matrix, name_col, p_col):
    fp = f"{RES}/{tsv}"
    if not os.path.exists(fp): print(f"  skip {tsv} (absent)"); return
    df = pd.read_csv(fp, sep="\t")
    df["q_BH"] = bh(df[p_col].values)
    # effect size from the matrix if available
    if matrix and os.path.exists(f"{A}/{matrix}"):
        M = pd.read_csv(f"{A}/{matrix}", sep="\t", index_col=0)
        pre = [s for s in M.index if s in PRE]; dro = [s for s in M.index if s in DRO]
        ds = []
        for nm in df[name_col]:
            if nm in M.columns:
                ds.append(cliffs(M.loc[pre, nm].values, M.loc[dro, nm].values))
            else:
                ds.append(np.nan)
        df["cliffs_delta"] = ds
        df["cliffs_mag"] = [mag(x) if pd.notna(x) else "" for x in ds]
    df.to_csv(fp, sep="\t", index=False)
    nsig = int((df["q_BH"] < 0.05).sum())
    print(f"  {tsv}: {len(df)} taxa, {nsig} pass q<0.05 (was p<0.05: {int((df[p_col]<0.05).sum())})")
    keep = [name_col, p_col, "q_BH"] + (["cliffs_delta","cliffs_mag"] if "cliffs_delta" in df else [])
    print(df.sort_values(p_col).head(6)[keep].to_string(index=False))

print("=== FDR + effect sizes ===")
annotate("phylum_shift_pre_vs_drought.tsv", "phylum_relabund.tsv", "phylum", "MWU_p")
annotate("genus_drought.tsv",   "genus_relabund.tsv",   "genus",   "p")
annotate("species_drought.tsv", "species_relabund.tsv", "species", "p")
annotate("amplicon_phylum_drought.tsv", None, "phylum", "MWU_p")   # different sample set; q only
print("DONE — q_BH + cliffs_delta written into the TSVs")
