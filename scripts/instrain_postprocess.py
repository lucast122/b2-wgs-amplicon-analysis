#!/usr/bin/env python3
"""
Strain-level drought test from inStrain profiles. For a target MAG (e.g. Nitrospira),
collect per-sample nucleotide diversity (pi), SNV/SNS density, coverage and breadth from each
sample's inStrain genome_info.tsv, then test the drought contrast. This probes WITHIN-population
microdiversity — the axis that read-abundance and gene-potential tests cannot see.

Usage: instrain_postprocess.py <label>   (label = subdir under assembly/instrain/)
Outputs to b2/analysis/results/: instrain_<label>.tsv (+ appended to a checks file).
"""
import os, glob, sys
import numpy as np, pandas as pd
from scipy import stats

LABEL = sys.argv[1] if len(sys.argv) > 1 else "Nitrospira"
ISROOT = f"/mnt/disk4/timo/gbi/b2/assembly/instrain/{LABEL}"
A = "/mnt/disk4/timo/gbi/b2/analysis"; RES = f"{A}/results"

def cliffs(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0: return np.nan
    return (sum((bj > a).sum() for bj in b) - sum((bj < a).sum() for bj in b)) / (len(a) * len(b))

def main():
    rows = []
    for d in sorted(glob.glob(f"{ISROOT}/*.IS")):
        acc = os.path.basename(d).replace(".IS", "")
        gi = glob.glob(f"{d}/output/*genome_info.tsv")
        if not gi: continue
        try: g = pd.read_csv(gi[0], sep="\t")
        except Exception: continue
        if not len(g): continue
        r = g.iloc[0]
        rows.append({"sample": acc,
                     "coverage": r.get("coverage", np.nan),
                     "breadth": r.get("breadth", np.nan),
                     "nucl_diversity": r.get("nucl_diversity", np.nan),
                     "SNV_count": r.get("SNV_count", np.nan),
                     "SNS_count": r.get("SNS_count", np.nan),
                     "divergent_sites": r.get("divergent_site_count", np.nan)})
    df = pd.DataFrame(rows)
    if not len(df):
        print(f"no inStrain genome_info for {LABEL} yet"); return
    # require reasonable breadth so pi is meaningful
    df = df[df["breadth"].astype(float) >= 0.5]
    meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
    df["condition"] = [meta.loc[s, "condition"] if s in meta.index else None for s in df["sample"]]
    df.to_csv(f"{RES}/instrain_{LABEL}.tsv", sep="\t", index=False)
    pre = df[df.condition == "pre-drought"]; dro = df[df.condition == "drought"]
    print(f"=== inStrain {LABEL}: {len(df)} samples with breadth>=0.5 (pre={len(pre)} drought={len(dro)}) ===")
    out = [f"label\t{LABEL}", f"n_samples_breadth50\t{len(df)}"]
    for metric in ["nucl_diversity", "coverage", "breadth", "SNV_count"]:
        a = pre[metric].astype(float).dropna(); b = dro[metric].astype(float).dropna()
        if len(a) < 3 or len(b) < 3:
            print(f"  {metric}: insufficient n"); continue
        try: _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError: p = np.nan
        d = cliffs(a, b)
        arrow = "UP" if b.mean() > a.mean() else "down"
        print(f"  {metric:16s} pre={a.mean():.4g} drought={b.mean():.4g} {arrow}  p={p:.3f}  Cliff δ={d:+.2f}")
        out.append(f"{metric}_pre\t{a.mean():.5g}"); out.append(f"{metric}_drought\t{b.mean():.5g}")
        out.append(f"{metric}_p\t{p:.4f}"); out.append(f"{metric}_cliffs\t{d:.3f}")
    open(f"{RES}/instrain_{LABEL}_checks.txt", "w").write("\n".join(out) + "\n")
    print(f"  -> results/instrain_{LABEL}.tsv")

if __name__ == "__main__":
    main()
