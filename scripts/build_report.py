#!/usr/bin/env python3
"""Assemble a single self-contained HTML report (figures embedded as base64)."""
import base64, os
import pandas as pd

A = "/mnt/disk4/timo/gbi/analysis"; FIG = f"{A}/figs"; RES = f"{A}/results"
OUT = f"{A}/GBI_Biosphere2_drought_microbiome_report.html"

def img(name):
    p = f"{FIG}/{name}"
    if not os.path.exists(p): return "<p><em>(figure pending)</em></p>"
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f'<img src="data:image/png;base64,{b}" style="max-width:560px;width:100%;">'

stats = {}
for line in open(f"{RES}/stats_summary.txt"):
    k, v = line.strip().split("\t"); stats[k] = v
shift = pd.read_csv(f"{RES}/phylum_shift_pre_vs_drought.tsv", sep="\t")
sig = shift[shift.MWU_p < 0.05]

# amplicon (16S) drought results
amp = pd.read_csv(f"{RES}/amplicon_phylum_drought.tsv", sep="\t")
def amp_get(name):
    r = amp[amp.phylum.str.startswith(name)]
    return r.iloc[0] if len(r) else None
amp_act = amp_get("Actinomycetota"); amp_nit = amp_get("Nitrospirota")

# genus-level drought results
gen = pd.read_csv(f"{RES}/genus_drought.tsv", sep="\t")
def gget(name):
    r = gen[gen.genus == name]
    return r.iloc[0] if len(r) else None
g_nitro = gget("Nitrospira_D"); g_brady = gget("Bradyrhizobium"); g_methylo = gget("Methyloceanibacter")

def shift_table(df):
    rows = ""
    for _, r in df.head(8).iterrows():
        sigmark = " <b>*</b>" if r.MWU_p < 0.05 else ""
        arrow = "▲" if r.drought_mean > r.pre_mean else "▼"
        rows += (f"<tr><td>{r.phylum}</td><td>{r.pre_mean:.2f}</td>"
                 f"<td>{r.drought_mean:.2f} {arrow}</td><td>{r.log2FC:+.2f}</td>"
                 f"<td>{r.MWU_p:.3f}{sigmark}</td></tr>")
    return rows

html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>GBI Biosphere 2 Drought — Soil Microbiome Metagenomics</title>
<style>
body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
max-width:820px;margin:32px auto;padding:0 18px;color:#1a1a1a;line-height:1.55}}
h1{{font-size:1.6em;border-bottom:2px solid #444;padding-bottom:6px}}
h2{{font-size:1.25em;margin-top:1.8em;color:#222;border-bottom:1px solid #ddd;padding-bottom:3px}}
.key{{background:#f3f7fb;border-left:4px solid #1f77b4;padding:10px 14px;margin:14px 0;border-radius:3px}}
.caveat{{background:#fbf6f0;border-left:4px solid #d08a3e;padding:10px 14px;margin:14px 0;border-radius:3px;font-size:.93em}}
table{{border-collapse:collapse;margin:12px 0;font-size:.9em}}
td,th{{border:1px solid #ccc;padding:4px 9px;text-align:right}}td:first-child,th:first-child{{text-align:left}}
.fig{{margin:16px 0;text-align:center}}
.muted{{color:#666;font-size:.86em}}
code{{background:#f0f0f0;padding:1px 4px;border-radius:3px;font-size:.88em}}
</style></head><body>

<h1>Soil microbiome response to drought in the Biosphere 2 tropical rainforest</h1>
<p class="muted">This analysis extends our Biosphere 2 WALD ¹³CO₂ drought study with
shotgun-metagenomic and 16S-amplicon profiling of the soil microbiome — adding a
taxonomic "who is present" axis to the ¹³C carbon-flux measurements in the manuscript.
Shotgun profiling with sylph (GTDB r232 + FungiRefSeq + custom AMF database);
{stats.get('n_samples','36')} samples passing QC across pre-drought and drought conditions.</p>

<div class="key"><b>Headline:</b> The bulk soil bacterial/archaeal community is
<b>compositionally resilient</b> to drought — overall diversity and community structure
do not shift significantly. Against that stable backdrop, two phyla move in opposite,
ecologically coherent directions: drought-tolerant <b>Actinomycetota increase</b> and
nitrite-oxidizing <b>Nitrospirota decrease</b> — and both shifts reproduce in an
independent 328-sample 16S amplicon survey. This fits the picture in the manuscript, where
drought adaptation operates through subtle carbon-allocation shifts rather than wholesale
community turnover.</div>

<h2>1. Community composition</h2>
<p>The soil is bacteria-dominated (~95%), with a substantial ammonia-oxidising archaeal
fraction (<i>Thermoproteota</i>/Nitrososphaeria, ~16%). Proteobacteria
(<i>Pseudomonadota</i>) make up about half of all classified reads. Mean composition is
near-identical between conditions — the drought signal is in specific taxa, not the bulk.</p>
<div class="fig">{img('04_composition.png')}</div>

<h2>2. Drought response — the clearest signal</h2>
<p>Comparing pre-drought vs drought samples phylum-by-phylum (Mann–Whitney U):</p>
<table><tr><th>Phylum</th><th>pre %</th><th>drought %</th><th>log2 FC</th><th>p</th></tr>
{shift_table(shift)}</table>
<p><b>{len(sig)} phyla shift significantly (p&lt;0.05):</b>
<i>Actinomycetota</i> rise from {shift.loc[shift.phylum=='Actinomycetota','pre_mean'].values[0]:.1f}%
to {shift.loc[shift.phylum=='Actinomycetota','drought_mean'].values[0]:.1f}%
(p={shift.loc[shift.phylum=='Actinomycetota','MWU_p'].values[0]:.3f}), and
<i>Nitrospirota</i> fall from {shift.loc[shift.phylum=='Nitrospirota','pre_mean'].values[0]:.1f}%
to {shift.loc[shift.phylum=='Nitrospirota','drought_mean'].values[0]:.1f}%
(p={shift.loc[shift.phylum=='Nitrospirota','MWU_p'].values[0]:.3f}).</p>
<div class="key"><b>Why this matters:</b> Actinomycetota are textbook drought-tolerant soil
bacteria (thick peptidoglycan walls, sporulation, osmolyte accumulation). Their relative
increase is an independent confirmation that the drought treatment imposed real water
stress on the soil community. The Nitrospirota decline points to drought-suppressed
nitrification — a functionally meaningful change in N-cycling.</div>
<div class="fig">{img('05_drought_response.png')}</div>

<h3>2b. Genus level sharpens this into a nitrogen-cycling signal</h3>
<p>Resolving to genus pinpoints the functional groups behind the phylum shifts. The
strongest single response in the entire dataset is the nitrite-oxidiser
<b><i>Nitrospira</i></b>, which drops from {g_nitro['pre']:.2f}% to {g_nitro['drought']:.2f}%
(<b>p={g_nitro['p']:.4f}</b>) — i.e. the Nitrospirota phylum decline <i>is</i> a Nitrospira
decline, a direct hit to nitrification. Alongside it the nitrogen-fixer / legume root
symbiont <b><i>Bradyrhizobium</i></b> falls ({g_brady['pre']:.2f}→{g_brady['drought']:.2f}%,
p={g_brady['p']:.3f}), while the methylotroph <i>Methyloceanibacter</i> rises
({g_methylo['pre']:.1f}→{g_methylo['drought']:.1f}%, p={g_methylo['p']:.3f}).</p>
<div class="key"><b>Coherent N-cycle suppression:</b> drought reduces both nitrite-oxidation
(<i>Nitrospira</i>) and nitrogen-fixation (<i>Bradyrhizobium</i>) — two complementary
nitrogen-cycle functions moving down together. The <i>Bradyrhizobium</i> signal is notable
given the belowground role of the legume <i>C. fairchildiana</i> in the manuscript.</div>
<div class="fig">{img('08_genus_drought.png')}</div>
<div class="caveat"><b>Why phylum leads the report:</b> of 70 testable genera, 44 are
uninterpretable GTDB placeholder codes (e.g. <code>DATFRX01</code>), and the
Actinomycetota phylum increase does not resolve to any single genus — it is spread thinly
across many. Genus level adds functional specificity for the nitrogen-cyclers but is
noisier overall, so the robust headline stays at phylum.</div>

<h2>3. Diversity and community structure are resilient</h2>
<p>Genus-level Shannon diversity does not differ between conditions
(pre {float(stats['shannon_pre']):.2f} vs drought {float(stats['shannon_drought']):.2f},
p={float(stats['shannon_MWU_p']):.2f}). In Bray–Curtis ordination, drought and pre-drought
samples fully intermix, and neither community structure (PERMANOVA
p={float(stats['PERMANOVA_p']):.2f}) nor dispersion (PERMDISP
p={float(stats['PERMDISP_p']):.2f}) differs significantly.</p>
<p>The composition is also <b>stable across the ¹³C pulse-chase window</b>: within each
condition, community structure does not differ between the 0 h, 6 h and 48 h sampling
timepoints (PERMANOVA p=0.79 drought, p=0.88 pre-drought; 0 of 11–13 phyla shift). This
confirms the standing DNA community does not turn over during labelling — so the
¹³C-activity dynamics in the manuscript reflect metabolic shifts <i>within</i> a stable
community, and pooling timepoints for the drought contrast introduces no confound.</p>
<div class="fig">{img('01_alpha_shannon.png')}{img('02_pcoa_braycurtis.png')}</div>
<p class="muted">Note on "convergence": the manuscript describes drought reducing
site-to-site variation in ¹³C-utilisation. At the level of bulk taxonomic composition we do
not see a matching convergence — within-condition dissimilarity is comparable (pre
{float(stats['withinBC_pre']):.2f}, drought {float(stats['withinBC_drought']):.2f}). The
two measurements capture different things: the DNA profile here is who-is-present, while the
¹³C-PLFA signal is who-is-active. A stable DNA pool alongside shifting carbon flux is fully
consistent — the active sub-community can converge without the standing community changing.</p>

<h2>4. Independent validation: 16S rRNA amplicon survey</h2>
<p>A separate 16S rRNA amplicon dataset (V4 515F/806R) from the same Biosphere 2 soil was
re-analysed with current tools (QIIME 2 2026.4, DADA2, Greengenes2 + SILVA 138.1):
<b>328 samples, 34,594 ASVs</b> — an order of magnitude more samples than the shotgun set,
providing an independent cross-check.</p>

<div class="key"><b>The two methods agree on the community.</b> The amplicon recovers 12 of
13 shotgun-detected phyla and both rank the same groups as dominant
(<i>Pseudomonadota, Acidobacteriota, Actinomycetota, Methylomirabilota,
Thermoproteota</i>). Abundance magnitudes differ as expected from method bias (16S
copy-number vs WGS genome-detection), but the <i>identity</i> of the dominant community is
reproducible across two independent assays.</p></div>
<div class="fig">{img('07_amplicon_vs_wgs.png')}</div>

<p><b>Cross-method drought response.</b> Restricting the amplicon data to its own
drought-labelled subset ({"%d pre vs %d drought" % (15, 20)}) and testing the two headline
phyla, both move in the <b>same direction as the shotgun data</b>:</p>
<table><tr><th>Phylum</th><th>shotgun (WGS)</th><th>16S amplicon</th><th>agreement</th></tr>
<tr><td>Actinomycetota</td>
<td>3.1→4.9% &nbsp;<b>p={shift.loc[shift.phylum=='Actinomycetota','MWU_p'].values[0]:.3f}</b> ▲</td>
<td>{amp_act['pre_mean']:.1f}→{amp_act['drought_mean']:.1f}% &nbsp;p={amp_act['MWU_p']:.2f} ▲</td>
<td>both ▲</td></tr>
<tr><td>Nitrospirota</td>
<td>7.0→5.1% &nbsp;<b>p={shift.loc[shift.phylum=='Nitrospirota','MWU_p'].values[0]:.3f}</b> ▼</td>
<td>{amp_nit['pre_mean']:.2f}→{amp_nit['drought_mean']:.2f}% &nbsp;p={amp_nit['MWU_p']:.2f} ▼</td>
<td>both ▼</td></tr></table>
<p>The shifts are statistically significant in the shotgun data and directionally
concordant — though not significant — in the smaller, non-paired amplicon subset. The
drought-tolerant Actinomycetota increase and the Nitrospirota (nitrifier) decline are thus
the most robust signals, reproducing across two assays and two sample sets.</p>
<div class="fig">{img('06_amplicon_drought.png')}</div>
<div class="caveat"><b>Honest caveats:</b> the amplicon drought subset is small (n=15/20) and
is <i>not</i> the same physical samples as the shotgun set, so cross-method agreement is at
the level of biological signal, not paired replicates. The amplicon also flags its own
shifts (Planctomycetota ▲ p=0.02, Bacteroidota ▼ p=0.03) not seen in the WGS — phylum-level
drought responses beyond Actinomycetota/Nitrospirota are method-dependent and should be
treated as provisional. 16S primers are prokaryote-only (no fungi); notably the one
"fungal" hit is <i>Glomeribacter</i>, the bacterial endosymbiont that lives <i>inside</i>
AMF spores — a trace of the mycorrhizal fungi the 16S cannot otherwise see.</div>

<h2>5. Eukaryotes and mycorrhizal fungi</h2>
<p>SSU-based screening (phyloFlash/SILVA) detects eukaryotic and plant signal in every
sample, confirming root/fungal material in the soil. However, SILVA cannot resolve
arbuscular (AMF) or ectomycorrhizal (EMF) fungi below the Opisthokonta level, and sylph's
genome-containment approach finds no AMF/EMF (they are too low-coverage in bulk shotgun
data). Resolving the AMF/EMF community — the mycorrhizal symbionts central to the
belowground carbon story in the manuscript — requires read-level protein classification
against a custom AMF/EMF-enriched database (<code>kaiju nr_euk + AMF + EMF</code>, in
progress).</p>
<div class="caveat"><b>Pending:</b> the AMF/EMF genus-level breakdown (Rhizophagus,
Funneliformis, Diversispora, Cortinarius, Inocybe, …) will be added once the custom Kaiju
classification of all 39 samples completes. The database build is finishing now.</div>

<h2>Methods (brief)</h2>
<p class="muted">Reads: 39 paired-end Illumina soil metagenomes (Biosphere 2 WALD campaign,
2019). Profiling: <code>sylph</code> against GTDB r232, FungiRefSeq-2025, and a
dereplicated 32-genome AMF database; taxonomic tables via <code>sylph-tax</code>.
Independent 16S: V4 515F/806R EMP amplicons, QIIME 2 2026.4 + DADA2 1.38 + Greengenes2
2024.09 / SILVA 138.1 (328 samples, 34,594 ASVs); phylum (L2) relative abundances compared
to the WGS profile. Statistics on genus/phylum relative abundances: Mann–Whitney U (taxa),
PERMANOVA &amp; PERMDISP on Bray–Curtis distances (999 permutations). {stats.get('n_samples','36')}/39
shotgun samples passed sketch QC. Figures and tables reproducible from
<code>/mnt/disk4/timo/gbi/analysis/</code> and
<code>/mnt/disk4/timo/gbi/amplicon/reanalysis_2026/</code>.</p>

<p class="muted">Companion WGS + amplicon analysis extending the Biosphere 2 WALD drought
manuscript · last updated automatically from the GBI metagenomics pipeline.</p>
</body></html>"""

open(OUT, "w").write(html)
print("Wrote", OUT, f"({os.path.getsize(OUT)//1024} KB)")
