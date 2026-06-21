# What to ship with the manuscript — dissemination & citation plan

Goal: maximise citations and community reuse of the B2 WALD soil-metagenomics companion work.
Principle: the artifacts that get independently cited are **deposited, DOI'd, standards-compliant
genomes and reusable tools** — not figures in a PDF. Below, in priority order, what to ship, where,
and why it earns citations.

## Tier A — the artifacts that get cited on their own

### 1. MAGs → ENA/NCBI BioProject (MIMAG-compliant)  ★ highest-impact
- **What:** the dereplicated genome-resolved MAGs (currently 88 HQ+MQ from 3/39 samples; expect
  several hundred across all 39, dereplicated to a non-redundant set).
- **Where:** ENA or NCBI under a dedicated **BioProject** (link to the existing SRA reads
  PRJ). Each MAG gets its own MIMAG BioSample + accession; >99.9% of MAGs get distinct BioSample
  IDs from the read sets, so they are independently findable/citable.
- **Why it cites:** genome-resolved MAGs from **Biosphere 2 tropical-rainforest soil** — an
  enclosed, understudied, largely *uncultured* community (≈80% of our reads are uncultured GTDB
  lineages). These genomes will surface in GTDB, GlobDB, and every future query for these taxa,
  and are reused/iterated (each reuse = a citation). This is the single best citation lever.
- **Metadata to attach:** CheckM2 completeness/contamination, GTDB-Tk taxonomy, assembler
  (MEGAHIT), binner (MetaBAT2), source sample + condition (pre/drought) + depth, N50.

### 2. Custom AMF/EMF Kaiju database → Zenodo (DOI)  ★ most novel/reusable tool
- **What:** the fungal-enriched Kaiju DB (nr_euk fungi + 32 dereplicated AMF genomes + EMF refs)
  + the build recipe + the `kaiju_postprocess.py` AMF/EMF classifier.
- **Where:** Zenodo (free, permanent DOI, keyword-indexed). FMI is large; if >50 GB, ship the
  **cleaned source FAA + build script + AMF genome manifest** instead of the index, so anyone can
  rebuild — smaller, more reusable, version-independent.
- **Why it cites:** detecting arbuscular/ectomycorrhizal fungi *from bulk shotgun* is something
  most pipelines (sylph, SSU, standard Kaiju) cannot do — we showed it recovers Rhizophagus-
  dominated AMF that those methods miss. A turnkey resource for the large mycorrhizal-ecology
  community = a tool paper-style citation magnet.

## Tier B — reproducibility & supplementary (credit + trust, supports reuse)

### 3. Code → GitHub release + Zenodo DOI
- The analysis repo (`github.com/lucast122/b2-wgs-amplicon-analysis`) already hosts the report;
  add the 33 analysis scripts + a README + environment spec, tag a release, and switch on the
  **GitHub–Zenodo integration** so the release mints a DOI. Cite that DOI in the paper.

### 4. Standardised supplementary tables → Zenodo / journal SI
- Taxonomic relative-abundance matrices (phylum/genus/species; sylph, Kaiju, DIAMOND).
- Functional HPM matrices (NCyc, dbCAN/CAZy, osmolyte) + drought-test tables.
- AMF/EMF Excel workbooks (already built).
- One machine-readable `samples.tsv` (SRR ↔ condition/site/plot/depth/timepoint) — the key to
  anyone re-using the reads.

### 5. Methods note (small but real): the nr_euk `*` corruption fix
- nr_euk FASTA contains `*` stop characters that silently corrupt a Kaiju FMI (the index builds
  but mis-classifies). Document the strip-and-rebuild fix; others hitting it will cite the
  methods/repo. Ship the one-line cleaning step in the README.

## How (concrete steps, in order)
1. Finish MAGs: all-39 assembly → MetaBAT2 → **CheckM2** filter (keep ≥50% / <10%) →
   **GTDB-Tk** taxonomy → **dRep** dereplicate → MIMAG metadata table.
2. Register a BioProject; bulk-submit MAGs (tooling: `subMG` / ENA webin / METAGENOTE automate
   the MIMAG submission).
3. Package the AMF/EMF DB bundle (cleaned FAA + build + manifest + classifier) → Zenodo, rich
   keywords ("arbuscular mycorrhiza", "metagenomics", "Kaiju", "shotgun", "Glomeromycota").
4. GitHub release → Zenodo DOI for code.
5. Drop the **Data & Code Availability** statement (drafted below) into the manuscript.

## Draft Data & Code Availability statement (for the manuscript)
> Raw shotgun metagenomic reads are available in the NCBI SRA under accessions SRR24887495–
> SRR24888648 (BioProject PRJNAxxxxxx). Metagenome-assembled genomes generated in this study are
> deposited in ENA/NCBI under BioProject PRJNAyyyyyy, with per-MAG MIMAG metadata (CheckM2 quality,
> GTDB-Tk taxonomy) in Supplementary Table Sx. The custom AMF/EMF-enriched Kaiju database and its
> classification scripts are archived on Zenodo (DOI: 10.5281/zenodo.xxxxxxx). All analysis code is
> available at github.com/lucast122/b2-wgs-amplicon-analysis (archived at Zenodo DOI:
> 10.5281/zenodo.yyyyyyy). Taxonomic and functional abundance matrices are provided as
> Supplementary Data.

## Notes
- Big Kaiju/DIAMOND nr_euk indexes (210 GB / 162 GB) are **not** worth depositing — ship recipes.
- Keep the report/site as the human-readable companion; the deposited artifacts are what cite.
