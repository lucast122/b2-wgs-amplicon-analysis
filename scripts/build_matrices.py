#!/usr/bin/env python3
"""
Parse sylph .sylphmpa taxonomic profiles for all 39 GBI Biosphere 2 drought
samples into tidy abundance matrices, merged with the experimental design map.

Outputs (to /mnt/disk4/timo/gbi/b2/analysis/):
  - phylum_relabund.tsv   : samples x phylum  (relative abundance %)
  - genus_relabund.tsv    : samples x genus
  - class_relabund.tsv    : samples x class
  - metadata.tsv          : cleaned design map keyed by SRR
  - merged_long.tsv       : long-format (sample, clade, rank, relabund, +metadata)
"""
import os, glob, re
import pandas as pd

TAXPROF = "/mnt/disk4/timo/gbi/b2/sylph_output/taxprof"
DESIGN  = "/mnt/disk4/timo/gbi/b2/sylph_output/sample_design_map.tsv"
OUT     = "/mnt/disk4/timo/gbi/b2/analysis"
os.makedirs(OUT, exist_ok=True)

# ---- metadata -------------------------------------------------------------
meta = pd.read_csv(DESIGN, sep="\t")
meta = meta.rename(columns={"Run": "sample"})
meta = meta.set_index("sample")
meta.to_csv(f"{OUT}/metadata.tsv", sep="\t")

# ---- parse each sylphmpa --------------------------------------------------
# Format: clade_name \t relative_abundance \t sequence_abundance \t ANI \t Cov
# clade_name is a |-delimited lineage: d__..|p__..|c__..|o__..|f__..|g__..|s__..
rank_prefix = {"d": "domain", "p": "phylum", "c": "class",
               "o": "order", "f": "family", "g": "genus", "s": "species"}

records = []
for fp in sorted(glob.glob(f"{TAXPROF}/*.sylphmpa")):
    # sample id from the #SampleID header line (path to fastq)
    srr = None
    with open(fp) as fh:
        for line in fh:
            if line.startswith("#SampleID"):
                m = re.search(r'(SRR\d+)', line)
                if m: srr = m.group(1)
                continue
            if line.startswith("clade_name") or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            clade, relab = parts[0], parts[1]
            try:
                relab = float(relab)
            except ValueError:
                continue
            # deepest rank in this lineage
            levels = clade.split("|")
            deepest = levels[-1]
            rp = deepest.split("__")[0]
            rank = rank_prefix.get(rp, "unknown")
            name = deepest.split("__", 1)[-1]
            records.append({"sample": srr, "clade": clade, "name": name,
                            "rank": rank, "relabund": relab})

df = pd.DataFrame(records)
df = df.dropna(subset=["sample"])
print(f"Parsed {df['sample'].nunique()} samples, {len(df)} clade records")

# ---- pivot per rank -------------------------------------------------------
def pivot_rank(rank):
    sub = df[df["rank"] == rank]
    # use the deepest-name as column; sum duplicates (shouldn't be many)
    mat = (sub.pivot_table(index="sample", columns="name",
                           values="relabund", aggfunc="sum")
              .fillna(0.0))
    return mat

for rank in ["phylum", "class", "genus", "species"]:
    mat = pivot_rank(rank)
    mat.to_csv(f"{OUT}/{rank}_relabund.tsv", sep="\t")
    print(f"  {rank}: {mat.shape[0]} samples x {mat.shape[1]} taxa")

# ---- merged long with metadata -------------------------------------------
merged = df.merge(meta, left_on="sample", right_index=True, how="left")
merged.to_csv(f"{OUT}/merged_long.tsv", sep="\t", index=False)
print(f"Wrote merged_long.tsv ({len(merged)} rows)")
print("DONE build_matrices")
