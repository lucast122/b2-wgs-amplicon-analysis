# B2 WGS + amplicon analysis

Companion shotgun-metagenomic (WGS) and 16S-amplicon analysis extending the
Biosphere 2 WALD ¹³CO₂ drought study — adding a taxonomic "who is present" axis to
the ¹³C carbon-flux measurements in the manuscript.

**📊 Report:** [`index.html`](index.html) — self-contained HTML (open in a browser, or
view via GitHub Pages once enabled).

## Key findings
- Bulk soil bacterial/archaeal community is **compositionally resilient** to drought
  (no significant shift in diversity or overall structure).
- Two phyla move coherently: drought-tolerant **Actinomycetota ↑**, nitrite-oxidising
  **Nitrospirota ↓** — reproduced in an independent 328-sample 16S amplicon survey.
- Genus level sharpens this to a **nitrogen-cycle signal**: *Nitrospira* (nitrite
  oxidation, p=0.0009) and *Bradyrhizobium* (N-fixation / legume symbiont) both decline.
- Composition is stable across the 0/6/48 h ¹³C pulse-chase window → the ¹³C-activity
  dynamics reflect metabolic shifts within a stable community, not community turnover.
- AMF/EMF mycorrhizal section pending custom Kaiju (nr_euk + AMF + EMF) classification.

## Layout
```
index.html          self-contained report (figures embedded)
scripts/            analysis pipeline (Python)
  build_matrices.py     parse sylph .sylphmpa -> abundance matrices
  diversity_analysis.py alpha/beta diversity, PERMANOVA, PERMDISP, phylum shifts
  composition_figs.py   composition + drought-response figures
  amplicon_drought.py   independent 16S drought test (QIIME2 GG2 L2 table)
  build_report.py       assemble the HTML report
results/            relative-abundance matrices, stats tables, summary
figs/               generated figures (PNG)
```

## Methods (brief)
- WGS: 39 paired-end Illumina soil metagenomes; profiling with `sylph` vs GTDB r232,
  FungiRefSeq-2025, and a dereplicated 32-genome AMF database; tables via `sylph-tax`.
  36/39 passed sketch QC.
- 16S: V4 515F/806R EMP amplicons; QIIME 2 2026.4 + DADA2 1.38 + Greengenes2 2024.09 /
  SILVA 138.1 (328 samples, 34,594 ASVs).
- Stats: Mann–Whitney U (taxa), PERMANOVA & PERMDISP on Bray–Curtis (999 permutations).

To regenerate the report: `python3 scripts/build_report.py`.
