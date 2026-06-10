#!/usr/bin/env python3
"""Build a standalone detailed METHODS page (methods.html) for the GitHub Pages
site, linked from the main report (index.html). Pulls a few live numbers from
results/ so it stays in sync. Self-contained styling matching the report."""
import os, pandas as pd
A = "/mnt/disk4/timo/gbi/analysis"; RES = f"{A}/results"
OUT = f"{A}/GBI_methods.html"

def stat(d, k, default="—"):
    return d.get(k, default)
stats = {}
if os.path.exists(f"{RES}/stats_summary.txt"):
    for line in open(f"{RES}/stats_summary.txt"):
        if "\t" in line: k, v = line.strip().split("\t"); stats[k] = v
ait = {}
if os.path.exists(f"{RES}/aitchison_stats.txt"):
    for line in open(f"{RES}/aitchison_stats.txt"):
        if "\t" in line: k, v = line.strip().split("\t"); ait[k] = v
pw = {}
if os.path.exists(f"{RES}/power_analysis.txt"):
    for line in open(f"{RES}/power_analysis.txt"):
        if "\t" in line: k, v = line.strip().split("\t"); pw[k] = v

CSS = """body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
max-width:820px;margin:32px auto;padding:0 18px;color:#1a1a1a;line-height:1.6}
h1{font-size:1.55em;border-bottom:2px solid #444;padding-bottom:6px}
h2{font-size:1.2em;margin-top:1.7em;color:#222;border-bottom:1px solid #ddd;padding-bottom:3px}
h3{font-size:1.02em;margin-top:1.2em;color:#333}
code{background:#f0f0f0;padding:1px 4px;border-radius:3px;font-size:.88em}
.muted{color:#666;font-size:.9em}.back{display:inline-block;margin:8px 0 18px;font-size:.95em}
table{border-collapse:collapse;margin:10px 0;font-size:.88em}td,th{border:1px solid #ccc;padding:4px 9px;text-align:left}
.note{background:#f3f7fb;border-left:4px solid #1f77b4;padding:8px 13px;margin:12px 0;border-radius:3px;font-size:.92em}"""

html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Detailed Methods — GBI Biosphere 2 Drought Microbiome</title>
<style>{CSS}</style></head><body>
<a class="back" href="index.html">&larr; Back to results report</a>
<h1>Detailed Methods</h1>
<p class="muted">Companion methods for <a href="index.html">“Soil microbiome response to drought in
the Biosphere 2 tropical rainforest.”</a> Shotgun + 16S amplicon profiling of Biosphere 2 WALD
soil, extending the ¹³CO₂ drought manuscript. All code and intermediate tables:
<code>github.com/lucast122/b2-wgs-amplicon-analysis</code>.</p>

<h2>1. Samples and experimental design</h2>
<p>39 paired-end Illumina shotgun metagenomes of Biosphere 2 tropical-rainforest soil
(WALD drought campaign, 2019), retrieved from SRA (accessions SRR24887495–SRR24888648).
The design is a 2 × 3 × 3 layout: <b>condition</b> (pre-drought, drought) ×
<b>timepoint</b> (0 h, 6 h, 48 h of a ¹³C pulse-chase) × <b>site</b> (Site1, Site2, and a
CTRL plot), with plot-level replication. {stat(stats,'n_samples','36')}/39 samples passed
sketch quality control; the 3 excluded (2 CTRL, 1 Site1) leave the CTRL plot with a single
shotgun sample, so site contrasts rest on Site1 vs Site2.</p>
<div class="note"><b>Design caveat (carried through all statistics):</b> the three timepoints
are repeated measures of the same physical plots and are therefore <i>not</i> independent.
Tests that pool samples (Mann–Whitney, PERMANOVA) treat them as independent and are
anticonservative; we therefore anchor every taxon claim on effect size and cross-assay /
cross-rank reproduction, and add an explicit plot-level mixed model and a plot-collapsed
sensitivity test (§5).</p>

<h2>2. Shotgun taxonomic profiling (sylph)</h2>
<p>Reads were profiled with <code>sylph</code> (k-mer genome containment) against
GTDB r232, FungiRefSeq-2025, and a dereplicated 32-genome arbuscular-mycorrhizal-fungi (AMF)
database; taxonomic relative-abundance tables (phylum→species) were produced with
<code>sylph-tax</code>. sylph estimates sequence-abundance via genome containment, so values
are relative abundances of detected genomes.</p>

<h2>3. 16S rRNA amplicon re-analysis (independent validation)</h2>
<p>A separate V4 515F/806R 16S dataset from the same Biosphere 2 soil was re-analysed with
QIIME 2 2026.4 (DADA2 1.38 denoising; Greengenes2 2024.09 + SILVA 138.1 taxonomy):
<b>328 samples, 34,594 ASVs</b>. Phylum (L2) relative abundances were compared to the
shotgun profile, and the drought-labelled amplicon subset was tested for the same headline
phyla. This is an independent assay and sample set, so agreement is at the level of
biological signal, not paired replicates.</p>

<h2>4. Read-level protein classification for mycorrhizal fungi (custom Kaiju)</h2>
<p>sylph’s genome-containment approach does not recover AMF/ectomycorrhizal (EMF) fungi
(too low-coverage in bulk soil shotgun). To resolve them we built custom <b>Kaiju</b>
protein databases and classified all 39 samples at the read level.</p>
<h3>4.1 Root-cause fix for index corruption</h3>
<p>Initial full <code>nr_euk</code> Kaiju indices classified 100% of reads as unclassified.
The cause was the <code>*</code> character (Prodigal stop codons) embedded in the protein
FASTA: Kaiju uses <code>*</code> as its internal sequence terminator, so any <code>*</code>
inside a sequence shatters the FM-index. The fix is to strip <code>*</code> and all
non-letter characters from the FASTA before building
(<code>awk '/^&gt;/{{...}}{{gsub(/[^A-Z]/,"",seq)}}'</code>), retaining X/B/Z/U/O/J as
wildcards. Indices were validated by back-translating the first database protein to DNA and
confirming it classifies (<code>C</code>) to its own taxon.</p>
<h3>4.2 Three databases</h3>
<table><tr><th>Database</th><th>Source</th><th>Build</th><th>Purpose</th></tr>
<tr><td>Fungi-only</td><td>cleaned fungi + AMF + EMF FASTA (17 Gbp)</td><td><code>kaiju-mkbwt</code> 16 threads</td><td>AMF/EMF genus matrices</td></tr>
<tr><td>Full nr_euk</td><td>cleaned nr_euk (150 Gbp)</td><td><code>kaiju-mkbwt</code> 16 threads, no soft RAM cap</td><td>specificity + domain fractions</td></tr>
<tr><td>DIAMOND</td><td>same cleaned nr_euk FASTA</td><td><code>diamond makedb</code></td><td>protein-level cross-check + functional</td></tr></table>
<p>Builds run on a 503 GB / 64-core host; <code>kaiju-mkbwt</code> thread count was tuned to
peak memory. Classification used Kaiju MEM mode; per-rank tables via
<code>kaiju2table -l superkingdom,phylum,class,order,family,genus</code>. AMF/EMF genera were
mapped against curated Glomeromycota/EMF genus lists with drought differential tests.</p>

<h2>5. Statistical analysis</h2>
<h3>5.1 Diversity and ordination</h3>
<p>Alpha diversity: genus-level Shannon index; pre-drought vs drought by Mann–Whitney U.
Beta diversity: Bray–Curtis on genus relative abundances, ordinated by Principal Coordinates
Analysis (PCoA); PERMANOVA and PERMDISP (999 permutations) tested condition structure and
multivariate dispersion. Because Bray–Curtis on proportions is compositional, we additionally
computed a <b>centred-log-ratio (CLR) / Aitchison</b> ordination (zeros handled by
multiplicative replacement) and PERMANOVA
(Aitchison PERMANOVA p={stat(ait,'PERMANOVA_Aitchison_p','—')},
PERMDISP p={stat(ait,'PERMDISP_Aitchison_p','—')}) as a compositionally-robust check.</p>
<h3>5.2 Differential abundance</h3>
<p>Per-taxon pre-vs-drought tests were run at phylum, genus and species rank with
Mann–Whitney U, then corrected for multiple testing with the Benjamini–Hochberg false
discovery rate (q-values) and reported with <b>Cliff’s δ</b> effect sizes
(negligible/small/medium/large by Romano thresholds). To guard against compositional
artifacts we cross-checked with three compositional methods from scikit-bio: a CLR Welch
t-test, the Dirichlet-multinomial t-test (<code>dirmult_ttest</code>), and ANCOM-BC
(<code>ancombc</code>).</p>
<h3>5.3 Pseudoreplication controls</h3>
<p>To respect the repeated-measures design we fit a Dirichlet-multinomial <b>linear mixed
model</b> (<code>dirmult_lme</code>) with condition as fixed effect and <b>plot as a random
effect</b>, and separately collapsed samples to one mean profile per plot and re-tested the
headline taxa (Mann–Whitney). Concordance across these design-respecting analyses, rather
than a single pooled p-value, is the basis for the reported signals.</p>
<h3>5.4 Convergence tests</h3>
<p>The manuscript’s “drought reduces site-to-site variation” hypothesis was tested two ways
on taxonomic composition: (i) within-condition Bray–Curtis dissimilarity (pre vs drought),
and (ii) between-site Bray–Curtis distance (Site1 vs Site2 pairs, pre vs drought), each by
one-sided Mann–Whitney (pre &gt; drought = convergence).</p>
<h3>5.5 Power</h3>
<p>With n={pw.get('n_pre','17')} pre-drought vs n={pw.get('n_drought','19')} drought samples,
the minimum effect detectable at 80% power is Cohen’s
d&nbsp;&asymp;&nbsp;{pw.get('min_detectable_Cohens_d_80pct','—')}; the design therefore has
limited power for small per-taxon shifts, and a non-significant result is read as
“no shift detectable at this sample size,” not as evidence of no effect.</p>

<h2>6. Functional and genome-resolved analyses</h2>
<p class="muted">The following extend the taxonomic survey toward mechanism and are added as
they complete (see report for current status).</p>
<h3>6.1 Functional gene profiling</h3>
<p>Reads are aligned (DIAMOND blastx) against curated functional databases to quantify
nitrogen-cycle marker genes (<i>amoA, nxrA/B, nifH, nirK/S, nosZ</i>; NCyc), carbohydrate-active
enzymes (dbCAN/CAZy), and desiccation/osmotic-stress pathways (<i>ectABC, otsAB</i>, glycine-betaine,
sporulation), tested pre vs drought.</p>
<h3>6.2 Assembly and MAG recovery</h3>
<p>Co-assembly (MEGAHIT), binning (MetaBAT2 / SemiBin2), quality (CheckM2), taxonomy
(GTDB-Tk), and dereplication (dRep) recover metagenome-assembled genomes — targeting the
large uncultured fraction of this enclosed ecosystem. Recovered genomes are screened for
stress pathways; strain dynamics of the dominant <i>Nitrospira</i> are assessed with inStrain
(popANI/SNV), and Actinomycetota biosynthetic gene clusters with antiSMASH.</p>
<h3>6.3 Ecological networks and viruses</h3>
<p>Per-condition genus co-occurrence networks (Spearman on CLR-transformed prevalent genera;
edges at |ρ|≥0.6 and BH-FDR q&lt;0.05; connectance compared by label permutation) test
network-level convergence. Viral contigs via geNomad; MAGs via MEGAHIT → MetaBAT2 → CheckM2 →
GTDB-Tk (databases staged locally under <code>ebt-metagenomics/resources</code>).</p>

<h2>7. Reproducibility</h2>
<p>Profiling/QC in dedicated conda environments; statistics in Python 3.12 with scikit-bio,
scipy, statsmodels and pandas (env <code>qiime2-amplicon-2026.4</code>). The report pipeline
runs <code>diversity_analysis.py → add_fdr_effectsizes.py → build_report.py</code>
(wrapper <code>run_report_pipeline.sh</code>). All scripts, result tables and figures are in
the repository.</p>
<a class="back" href="index.html">&larr; Back to results report</a>
<p class="muted">Auto-generated methods page · GBI metagenomics pipeline.</p>
</body></html>"""
open(OUT, "w").write(html)
print("Wrote", OUT, f"({os.path.getsize(OUT)//1024} KB)")
