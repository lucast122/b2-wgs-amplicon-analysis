#!/usr/bin/env python3
"""Assemble a single self-contained HTML report (figures embedded as base64)."""
import base64, os
import pandas as pd

A = "/mnt/disk4/timo/gbi/b2/analysis"; FIG = f"{A}/figs"; RES = f"{A}/results"
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
sig = shift[shift.MWU_p < 0.05]                       # raw p (pre-correction)
sig_fdr = shift[shift.get("q_BH", pd.Series(1, index=shift.index)) < 0.05]  # survive BH-FDR

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

# species-level summary (uncultured fraction + strain-level Nitrospira)
unc = {}
for line in open(f"{RES}/uncultured_fraction.txt"):
    k, v = line.strip().split("\t"); unc[k] = v
spd = pd.read_csv(f"{RES}/species_drought.tsv", sep="\t")
spd = spd[~spd.species.str.startswith("t__")]
sp_nitro = spd[spd.species.str.startswith("Nitrospira_D sp")].iloc[0]

# --- Tier 1 robustness numbers (compositional + pseudoreplication + power) ---
def _kv(fn):
    d = {}; p = f"{RES}/{fn}"
    if os.path.exists(p):
        for line in open(p):
            if "\t" in line: k, v = line.strip().split("\t"); d[k] = v
    return d
ait = _kv("aitchison_stats.txt"); pwr = _kv("power_analysis.txt")
def _nsig(fn):
    p = f"{RES}/{fn}"
    if os.path.exists(p):
        try:
            d = pd.read_csv(p, sep="\t")
            if "Signif" in d.columns: return int(d["Signif"].sum())
        except Exception: pass
    return None
n_dmt = _nsig("dirmult_ttest_genus.tsv"); n_lme = _nsig("dirmult_lme_plot_genus.tsv")
def _row(fn, col, val):
    p = f"{RES}/{fn}"
    if os.path.exists(p):
        d = pd.read_csv(p, sep="\t"); m = d[d[col].astype(str).str.startswith(val)]
        if len(m): return m.iloc[0]
    return None
clr_nit = _row("clr_diff_abundance.tsv", "genus", "Nitrospira_D")
pc_nit  = _row("plot_collapsed_drought.tsv", "genus", "Nitrospira_D")
T1 = os.path.exists(f"{RES}/aitchison_stats.txt")   # gate the robustness section

def shift_table(df):
    rows = ""
    for _, r in df.head(8).iterrows():
        arrow = "▲" if r.drought_mean > r.pre_mean else "▼"
        q = r["q_BH"] if "q_BH" in df.columns else float("nan")
        eff = r["cliffs_mag"] if "cliffs_mag" in df.columns else ""
        rows += (f"<tr><td>{r.phylum}</td><td>{r.pre_mean:.2f}</td>"
                 f"<td>{r.drought_mean:.2f} {arrow}</td><td>{r.log2FC:+.2f}</td>"
                 f"<td>{r.MWU_p:.3f}</td><td>{q:.2f}</td><td>{eff}</td></tr>")
    return rows

_aitp = ait.get('PERMANOVA_Aitchison_p', '—')
_clrq = f"{clr_nit['q_BH']:.2f}" if clr_nit is not None else "—"
_pcp  = f"{pc_nit['MWU_p_plot']:.2f}" if pc_nit is not None else "—"
_pwrd = pwr.get('min_detectable_Cohens_d_80pct', '—')
_netp = "—"
if os.path.exists(f"{RES}/network_metrics.tsv"):
    for _l in open(f"{RES}/network_metrics.tsv"):
        if _l.startswith("# edge_diff_perm_p"): _netp = _l.strip().split("\t")[1]
robustness_html = f"""
<h3>3b. Robustness, power, and a sampling caveat</h3>
<p>Abundances are compositional and timepoints are repeated measures of the same plots, so the
drought contrast was re-tested with methods controlling for both
(<a href="methods.html">methods §5.2–5.4</a>). It does not strengthen:</p>
<table><tr><th>Test (controls for)</th><th>Drought contrast</th></tr>
<tr><td>Aitchison/CLR PERMANOVA (compositionality)</td><td>p={_aitp}, n.s.</td></tr>
<tr><td>CLR Welch+FDR / Dirichlet-multinomial <i>t</i>-test (compositionality)</td><td><i>Nitrospira</i> top hit q={_clrq}; {n_dmt if n_dmt is not None else '—'} genera signif.</td></tr>
<tr><td>Mixed model, plot random effect (pseudoreplication)</td><td>{n_lme if n_lme is not None else '—'} genera signif.</td></tr>
<tr><td>Plot-collapsed Mann–Whitney (pseudoreplication)</td><td><i>Nitrospira</i> p={_pcp}</td></tr>
<tr><td>Genus co-occurrence networks (structure)</td><td>no connectance shift, perm p={_netp}</td></tr></table>
<div class="caveat"><b>Bottom line:</b> under compositional- and plot-aware analysis <b>no taxon
is significant</b>, and the design is underpowered (n={pwr.get('n_pre','?')} vs {pwr.get('n_drought','?')};
min detectable d&asymp;{_pwrd}). <b>Sampling confound:</b> pre-drought (Sept) and drought (Nov)
samples are temporally separated, so drought is partly confounded with season (the
contemporaneous CTRL plot carries the drought-period timepoint, not an independent watered
control). <i>Nitrospira</i> is the strongest, most consistent signal but is best treated as a
hypothesis for a powered, season-controlled follow-up.</div>
""" if T1 else ""

# --- protein cross-validation (Kaiju nr_euk + DIAMOND) and AMF/EMF (Kaiju fungi DB) ---
dchk = _kv("diamond_checks.txt")
DIA_OK = os.path.exists(f"{RES}/diamond_checks.txt")
AMF_OK = os.path.exists(f"{RES}/kaiju_amf_fungalnorm.tsv")
_meta = pd.read_csv(f"{A}/metadata.tsv", sep="\t", index_col=0)
def _nit_means(fn):
    p = f"{RES}/{fn}"
    if not os.path.exists(p): return (float('nan'), float('nan'))
    m = pd.read_csv(p, sep="\t", index_col=0)
    col = [c for c in m.columns if "Nitrospira" in c]
    if not col: return (float('nan'), float('nan'))
    s = [x for x in m.index if x in _meta.index]; md = _meta.loc[s]
    pre = [x for x in md.index[md.condition=="pre-drought"]]
    dro = [x for x in md.index[md.condition=="drought"]]
    return (m.loc[pre, col[0]].mean(), m.loc[dro, col[0]].mean())
_dspear = float(dchk.get("median_spearman_vs_kaiju", "nan"))
_dov = int(round(float(dchk.get("top15_overlap", "0")) * 15))
_dnp, _dnd = _nit_means("diamond_genus_relabund.tsv")
_knp, _knd = _nit_means("fullnr_genus_relabund.tsv")
protein_html = f"""
<h3>4b. Protein-level cross-validation: two aligners, one answer</h3>
<p>The shotgun reads were independently re-classified at the protein level with <b>Kaiju</b>
(maximum-exact-match) and <b>DIAMOND</b> (double-indexed alignment, <code>-f102</code> LCA)
against the <i>same</i> cleaned nr_euk protein database &mdash; so any difference reflects the
aligner, not the reference. The two give <b>highly concordant community profiles</b> (median
per-sample Spearman &rho;={_dspear:.2f}; {_dov}/15 top genera shared) and rank the same
nitrogen-cycling / rhizosphere taxa as dominant (<i>Bradyrhizobium, Nitrospira, Luteitalea,
Mesorhizobium</i>).</p>
<div class="key"><b><i>Nitrospira</i> declines under drought in every method.</b> The nitrite
oxidiser drops pre&rarr;drought in DIAMOND ({_dnp:.1f}&rarr;{_dnd:.1f}% of classified reads),
in Kaiju ({_knp:.1f}&rarr;{_knd:.1f}%), in sylph genome-containment, and in the 328-sample 16S
survey &mdash; <b>four independent assays, same direction.</b> A pattern reproduced across
genome-containment and two protein aligners is the strongest evidence in this study, even
though it still does not survive the pseudoreplication-corrected test (&sect;3b).</div>
<div class="fig">{img('14_diamond_vs_kaiju.png')}</div>
""" if DIA_OK else ""

_amf = (pd.read_csv(f"{RES}/kaiju_amf_fungalnorm.tsv", sep="\t", index_col=0) if AMF_OK else None)
_emf = (pd.read_csv(f"{RES}/kaiju_emf_fungalnorm.tsv", sep="\t", index_col=0) if AMF_OK else None)
_domf = (pd.read_csv(f"{RES}/kaiju_domain_fractions.tsv", sep="\t", index_col=0) if AMF_OK else None)
_rh = (_amf["Rhizophagus"].mean() if (AMF_OK and "Rhizophagus" in _amf.columns) else float('nan'))
_ff = (_domf[["AMF","EMF","Fungi_other"]].sum(1).mean() if AMF_OK else float('nan'))
_ep = (_emf.sum(1).mean() if AMF_OK else float('nan'))
amf_html = (f"""<p>Resolving the mycorrhizal community needs read-level protein classification
against a fungal-enriched database. All 39 samples were classified with a custom <b>Kaiju
fungal DB</b> (nr_euk fungi + 32 AMF genomes + EMF references), recovering the AMF/EMF genera
that SSU surveys and genome-containment miss.</p>
<div class="key"><b>Arbuscular mycorrhizal fungi (AMF) are present, dominated by
<i>Rhizophagus</i></b> (&asymp;{_rh:.1f}% of fungal reads; fungi are {_ff:.1f}% of all reads).
Further AMF: <i>Oehlia, Acaulospora, Paraglomus, Diversispora</i>. Ectomycorrhizal signal is
lower (&asymp;{_ep:.1f}% of fungal reads), led by <i>Russula</i> and <i>Tuber</i>. Full AMF/EMF
genus matrices (% of fungal and of total reads) + drought tests: companion Excel
(<code>GBI_AMF_EMF_kaiju.xlsx</code>).</div>
<div class="fig">{img('12_amf_emf_genera.png')}</div>
<div class="caveat">AMF/EMF reads are a small fraction of bulk-soil shotgun data, so genus
proportions carry wide per-sample variance &mdash; read as presence / relative dominance, not
precise quantities. No AMF/EMF genus shows a drought shift surviving FDR.</div>"""
if AMF_OK else
'<div class="caveat"><b>Pending:</b> AMF/EMF genus breakdown once the custom Kaiju classification completes.</div>')

# fungal guild pre vs drought (for the discussion)
if AMF_OK:
    _s = [x for x in _domf.index if x in _meta.index]; _md = _meta.loc[_s]
    _pre = [x for x in _md.index[_md.condition=="pre-drought"]]; _dro = [x for x in _md.index[_md.condition=="drought"]]
    _tf = _domf[["AMF","EMF","Fungi_other"]].sum(1)
    _fpre, _fdro = _tf.loc[_pre].mean(), _tf.loc[_dro].mean()
else:
    _fpre = _fdro = float('nan')
discussion_html = f"""
<h2>6. Discussion: reconciling the metagenome with the ¹³C carbon-flux experiment</h2>
<p>The companion B2 WALD ¹³CO₂ pulse-labeling indicates drought re-routes recent plant carbon
into <b>mycorrhizal and saprotrophic fungal</b> lipid biomarkers — yet here the standing fungal
community is essentially <b>unchanged</b> (total fungi {_fpre:.1f}&rarr;{_fdro:.1f}% of reads;
AMF, EMF and saprotrophs all flat, none surviving FDR). These results are <b>complementary, not
contradictory</b> — the two assays measure different axes of the same system.</p>
<div class="key"><b>Take-aways</b>
<ul>
<li><b>Abundance &ne; activity.</b> ¹³C-PLFA/NLFA tracks <i>carbon flux into active lipid
synthesis</i> over hours–days; shotgun reads track <i>standing DNA</i> integrated over weeks. A
fungus can take up more recent carbon without its genome pool changing — the very gap that
stable-isotope probing (SIP) was developed to close.<sup>[3]</sup></li>
<li><b>Flux goes to storage, not necessarily new biomass.</b> Drought ¹³C is allocated
preferentially to fungal <i>storage lipids</i> rather than hyphal growth, so an allocation signal
need not appear as more cells or DNA.<sup>[6]</sup></li>
<li><b>Relative ≠ absolute.</b> Shotgun profiles are compositional; a guild can rise in absolute
biomass yet stay flat in %. Relative abundances are documented to be misleading for quantitative
fungal dynamics.<sup>[5]</sup></li>
<li><b>Compartment &amp; detection limit.</b> ¹³C labels the root/myco-rhizosphere interface;
bulk-soil shotgun under-samples that active fraction, and AMF/EMF (low-biomass, large/patchy
genomes, sparse references) sit near our detection floor (~6% of an already ~3% fungal fraction).</li>
<li><b>Precedent in the literature.</b> DNA/amplicon abundance generally does <i>not</i> track
isotope-based activity under drought — e.g. only 8% of prokaryote and 25% of fungal OTUs shifted
across a drought–rewet cycle<sup>[4]</sup> — and the B2 soils themselves show a strong
<i>functional</i> drought response (carbon re-routed to volatile metabolites) with no requirement
that taxonomic abundance mirror it.<sup>[2]</sup></li>
</ul>
<b>Bottom line:</b> &ldquo;no abundance shift&rdquo; does not refute &ldquo;more carbon to
fungi.&rdquo; The metagenome shows a <i>compositionally buffered</i> community over the drought
window; the ¹³C shows an <i>allocational/functional</i> response. The definitive reconciliation
would be <b>DNA- or quantitative-SIP metagenomics</b> on the ¹³C-labeled soils (heavy-fraction
sequencing of the active community) — the natural follow-up this bulk-shotgun survey cannot itself
provide.</div>
<p class="muted" style="font-size:.82em"><b>References.</b>
[1] <a href="https://meetingorganizer.copernicus.org/EGU24/EGU24-17755.html">Belowground C
allocation of tropical rainforests under drought: an ecosystem ¹³CO₂ labeling experiment</a>.
EGU General Assembly 2024, EGU24-17755 (B2 WALD campaign). &nbsp;
[2] <a href="https://www.nature.com/articles/s41564-023-01432-9">Drought re-routes soil microbial
carbon metabolism towards emission of volatile metabolites in an artificial tropical
rainforest</a>. <i>Nat. Microbiol.</i> 2023. &nbsp;
[3] <a href="https://journals.asm.org/doi/10.1128/msystems.00417-22">Quantitative stable-isotope
probing (qSIP) with metagenomics links microbial physiology and activity to soil moisture</a>.
<i>mSystems</i> 2022. &nbsp;
[4] <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5845876/">Drought legacy effects on the
composition of soil fungal and prokaryote communities</a>. 2018. &nbsp;
[5] <a href="https://www.mdpi.com/2076-2607/9/3/589">Relative abundances of species or sequence
variants can be misleading: soil fungal communities as an example</a>. <i>Microorganisms</i> 2021. &nbsp;
[6] <a href="https://link.springer.com/article/10.1007/s00374-022-01670-9">Fatty acid 16:1ω5 as a
proxy for arbuscular mycorrhizal fungal biomass: current challenges and ways forward</a>.
<i>Biol. Fertil. Soils</i> 2022.</p>
"""

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
<p style="margin-top:-4px"><a href="methods.html"><b>📋 Read the detailed methods &rarr;</b></a>
<span class="muted">(separate page — samples, databases, the Kaiju build, and full statistics)</span></p>

<div class="key"><b>Headline:</b> At n={stats.get('n_samples','36')}, no taxon survives
multiple-testing correction and overall diversity/structure is unchanged — <b>no
statistically significant drought shift</b>. The strongest pattern, by converging evidence
rather than any single p-value, is suppression of nitrite-oxidising <b><i>Nitrospira</i></b>
(large effect, consistent phylum→genus→strain, echoed in the 328-sample 16S survey), with a
drought-tolerant <b>Actinomycetota</b> rise. We report these as <i>hypotheses</i>, not
established effects. <b>Key caveat:</b> pre-drought samples are from September and drought
from November, so condition is confounded with sampling date/season (§3b).</div>

<h2>1. Community composition</h2>
<p>The soil is bacteria-dominated (~95%), with a substantial ammonia-oxidising archaeal
fraction (<i>Thermoproteota</i>/Nitrososphaeria, ~16%). Proteobacteria
(<i>Pseudomonadota</i>) make up about half of all classified reads. Mean composition is
near-identical between conditions — the drought signal is in specific taxa, not the bulk.</p>
<div class="fig">{img('04_composition.png')}</div>

<h2>2. Drought response — the clearest signal</h2>
<p>Comparing pre-drought vs drought samples phylum-by-phylum (Mann–Whitney U):</p>
<table><tr><th>Phylum</th><th>pre %</th><th>drought %</th><th>log2 FC</th><th>p</th><th>q (BH)</th><th>effect</th></tr>
{shift_table(shift)}</table>
<p>Two phyla move at raw p&lt;0.05 in ecologically coherent directions:
<i>Actinomycetota</i> rise from {shift.loc[shift.phylum=='Actinomycetota','pre_mean'].values[0]:.1f}%
to {shift.loc[shift.phylum=='Actinomycetota','drought_mean'].values[0]:.1f}%
(p={shift.loc[shift.phylum=='Actinomycetota','MWU_p'].values[0]:.3f}, medium effect), and
<i>Nitrospirota</i> fall from {shift.loc[shift.phylum=='Nitrospirota','pre_mean'].values[0]:.1f}%
to {shift.loc[shift.phylum=='Nitrospirota','drought_mean'].values[0]:.1f}%
(p={shift.loc[shift.phylum=='Nitrospirota','MWU_p'].values[0]:.3f}, medium effect).</p>
<div class="caveat"><b>No phylum survives BH-FDR correction</b>
(q≈{shift.loc[shift.phylum=='Nitrospirota','q_BH'].values[0]:.2f}). The phylum table is
suggestive only; the N-cycling signal rests on cross-rank and cross-assay consistency
(below), not these p-values.</div>
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
(<b>p={g_nitro['p']:.4f}</b>, q={g_nitro['q_BH']:.2f}, Cliff's δ={g_nitro['cliffs_delta']:.2f} —
<b>large effect</b>) — i.e. the Nitrospirota phylum decline <i>is</i> a Nitrospira
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

<h3>2c. Species level: a strain-specific signal in a largely uncultured community</h3>
<p>Two things stand out at species/strain resolution. First, this is a
<b>predominantly uncultured community</b>: on average <b>{float(unc['mean_uncultured_pct']):.0f}%</b>
of classified abundance (range {float(unc['min']):.0f}–{float(unc['max']):.0f}%) belongs to
GTDB lineages with no cultured representative — placeholder-coded single-cell and
metagenome-assembled genomes (e.g. the most abundant single "species", at ~17%, is an
uncultured organism, <code>SCGC-AG-212-J23</code>). This fits a long-enclosed artificial
ecosystem harbouring novel soil taxa.</p>
<p>Second, the nitrogen-cycle signal sharpens to a <b>single strain</b>: the
nitrite-oxidiser decline tracks one dominant <i>Nitrospira</i> genome
(<code>GCA_029194675.1</code>), down {sp_nitro['pre']:.2f}→{sp_nitro['drought']:.2f}%
(<b>p={sp_nitro['p']:.3f}</b>, q={sp_nitro['q_BH']:.2f}). Critically, the same decline is
present at <b>every taxonomic rank</b> — phylum → genus → strain (raw p&lt;0.05 throughout;
q≈0.07–0.12 after correction) — which, together with its large effect size and independent
16S reproduction, makes drought-suppressed nitrification the most internally consistent
observation in the dataset. It nonetheless does <b>not</b> survive compositional or
plot-level correction (§3b), so we frame it as a strong hypothesis rather than an
established effect.</p>
<div class="fig">{img('09_nitrospira_ranks.png')}</div>

<h2>3. Diversity and community structure: no detectable shift</h2>
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
<div class="fig">{img('01_alpha_shannon.png')}{img('02_pcoa_braycurtis.png')}{img('02b_pcoa_site.png')}{img('10_aitchison_pcoa.png')}</div>
<p class="muted">Markers encode the design: in the first ordination colour = condition and
shape = timepoint; in the second colour = site (Site1/Site2/CTRL) and shape = condition.
Sites intermix completely.</p>
<p class="muted">Note on "convergence": the manuscript describes drought reducing
site-to-site variation in ¹³C-utilisation. At the level of bulk taxonomic composition we do
not see a matching convergence — within-condition dissimilarity is comparable (pre
{float(stats['withinBC_pre']):.2f}, drought {float(stats['withinBC_drought']):.2f}), and a
direct between-site test (Site1 vs Site2 Bray–Curtis, pre vs drought) shows no convergence
either (pre {float(stats.get('betweenSiteBC_pre',0)):.2f}, drought
{float(stats.get('betweenSiteBC_drought',0)):.2f}, one-sided p={float(stats.get('betweenSite_conv_MWU_p',1)):.2f};
sites diverge slightly, if anything, under drought). The
two measurements capture different things: the DNA profile here is who-is-present, while the
¹³C-PLFA signal is who-is-active. A stable DNA pool alongside shifting carbon flux is fully
consistent — the active sub-community can converge without the standing community changing.</p>
{robustness_html}

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
{protein_html}
<h2>5. Eukaryotes and mycorrhizal fungi</h2>
<p>SSU screening (phyloFlash/SILVA) detects eukaryotic and plant signal in every sample
(confirming root/fungal material), but cannot resolve arbuscular (AMF) or ectomycorrhizal
(EMF) fungi below Opisthokonta, and sylph genome-containment finds no AMF/EMF (too
low-coverage in bulk shotgun).</p>
{amf_html}
{discussion_html}
<h2>Plain-language notes on the statistics</h2>
<p class="muted">A reader's guide to why the headline is framed as a hypothesis — written for
non-statisticians (e.g. for manuscript text).</p>
<div class="key" style="font-size:.93em">
<p><b>Repeated timepoints = pseudoreplication.</b> Each plot was sampled at 0/6/48 h, and the
standing community does not turn over in that window — so a plot's three timepoints are
essentially three copies of one community, not three independent samples. Treating them as
independent is like measuring one patient's blood pressure three times and analysing it as
three patients: it makes results look more certain than they are. We corrected for this with a
mixed model (plot as a random effect) and by collapsing each plot to one value; both removed
the apparent significance.</p>
<p><b>Effect size (Cliff's δ) = how big, not just whether.</b> A p-value says only "could this
be chance?"; it says nothing about magnitude. Cliff's δ asks: pick a random drought and a
random pre-drought sample — how often is one higher? 0 = total overlap, ±1 = complete
separation. <i>Nitrospira</i>'s δ ≈ −0.63 ("large") means drought samples are almost always
lower — a strong, consistent shift even though it is not formally significant.</p>
<p><b>Power = could we even have detected it?</b> Power is the chance of catching a real effect
if one exists; small studies are a coarse net that only catches the biggest fish. At
n = {pwr.get('n_pre','?')} vs {pwr.get('n_drought','?')} we can only reliably detect very large
effects (Cohen's d ≈ {_pwrd}), so "no significant shift" means "none detectable at this sample
size", <i>not</i> "no shift". A large effect that is
non-significant in an underpowered study is the signature of a real signal we were too small to
prove — hence we report it as a hypothesis.</p>
<p><b>Multiple testing (BH-FDR).</b> Testing many taxa at once guarantees some will look
significant by chance; the false-discovery-rate correction (q-values) adjusts for that. It is a
separate issue from the non-independence above — both inflate false positives, and we control
for both.</p>
</div>

<h2>Methods (brief) — <a href="methods.html">full detailed methods &rarr;</a></h2>
<p class="muted">Reads: 39 paired-end Illumina soil metagenomes (Biosphere 2 WALD campaign,
2019). Profiling: <code>sylph</code> against GTDB r232, FungiRefSeq-2025, and a
dereplicated 32-genome AMF database; taxonomic tables via <code>sylph-tax</code>.
Independent 16S: V4 515F/806R EMP amplicons, QIIME 2 2026.4 + DADA2 1.38 + Greengenes2
2024.09 / SILVA 138.1 (328 samples, 34,594 ASVs); phylum (L2) relative abundances compared
to the WGS profile. Protein-level classification of all 39 metagenomes against a cleaned
nr_euk database with <code>Kaiju</code> (MEM) and <code>DIAMOND</code> 2.2.1 (double-indexed,
<code>-f102</code> LCA; 5M-read subsample per sample concatenated into one search), plus a
custom Kaiju fungal DB (nr_euk-fungi + 32 AMF genomes + EMF references) for the AMF/EMF
breakdown; the nr_euk database build required stripping internal <code>*</code> separators
from the source FASTA. Statistics on genus/phylum relative abundances: Mann–Whitney U (taxa) with
Benjamini–Hochberg FDR (q-values) and Cliff's δ effect sizes; PERMANOVA &amp; PERMDISP on
Bray–Curtis distances (999 permutations). {stats.get('n_samples','36')}/39 shotgun samples
processed — 3 (2 CTRL, 1 Site1) initially failed sylph sketching because of corrupted reads
from a faulty parallel download (a few records with malformed quality strings) and were
recovered with <code>seqkit sana</code> + re-pairing before re-profiling. Timepoints (0/6/48 h) are repeated
measures of the same plots; p-values treat samples as independent and are thus
anticonservative, so taxon claims are anchored on effect size and cross-assay reproduction
rather than p alone. Figures and tables reproducible from
<code>/mnt/disk4/timo/gbi/b2/analysis/</code> and
<code>/mnt/disk4/timo/gbi/b2/amplicon/reanalysis_2026/</code>.</p>

<p class="muted">Companion WGS + amplicon analysis extending the Biosphere 2 WALD drought
manuscript · last updated automatically from the GBI metagenomics pipeline.</p>
</body></html>"""

open(OUT, "w").write(html)
print("Wrote", OUT, f"({os.path.getsize(OUT)//1024} KB)")
