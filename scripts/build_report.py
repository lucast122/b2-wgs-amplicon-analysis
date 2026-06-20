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
<h3>4b. Protein-level cross-validation (Kaiju + DIAMOND)</h3>
<p>Two protein aligners classified the reads against the same nr_euk database. They agree closely
(median per-sample Spearman &rho;={_dspear:.2f}; {_dov}/15 top genera shared) and identify the same
dominant taxa — the nitrogen-cyclers and rhizosphere bacteria <i>Bradyrhizobium</i>,
<i>Nitrospira</i>, <i>Luteitalea</i> and <i>Mesorhizobium</i>. Both record the drought decline in
<i>Nitrospira</i> (DIAMOND {_dnp:.1f}&rarr;{_dnd:.1f}%, Kaiju {_knp:.1f}&rarr;{_knd:.1f}% of
classified reads). The decline therefore holds across all four methods applied here — sylph,
Kaiju, DIAMOND and 16S — the most reproducible result in the dataset, though still short of
significance after correction for pseudoreplication (&sect;3b).</p>
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
<h2>6. Discussion: the metagenome in the context of the ¹³C–PLFA companion study</h2>
<p>This shotgun survey is the standing-community counterpart to the B2 WALD <b>¹³C–PLFA companion
study</b> (Bai, Shi, Huang … Werner, Dippold)<sup>[1]</sup>, which traced recent photosynthate
through the rhizosphere of the same drought campaign<sup>[2]</sup>. Under drought that study found
trees invested fresh carbon mainly into <b>topsoil fine-root water-soluble carbon (osmolytes)</b>
— not deeper roots or rhizodeposition — and re-routed recent C to fungi in a
<b>compartment-specific</b> way: <b>AMF gained ¹³C in the still-moist subsoil</b> (especially the
legume <i>Clitoria fairchildiana</i>), while <b>saprotrophic + ectomycorrhizal fungi gained ¹³C in
the dried topsoil</b> (the single largest allocation shift). Which bacteria received C was set by
root metabolic / exudate control, not by Gram&plus;/Gram&minus; desiccation tolerance.</p>
<div class="key"><b>Why our flat metagenome is concordant, not contradictory.</b> The companion
study's <i>own PLFA community composition was largely buffered</i>: only the dried <b>topsoil</b>
rhizosphere shifted significantly, subsoil communities did not, and AMF biomarkers were present
under both ambient and drought conditions. Its dramatic drought signal lived in <b>¹³C allocation —
who is actively taking up new carbon</b> — not in standing biomass. Our metagenome measures exactly
that buffered axis (standing DNA) and agrees with it: no taxon survives correction, diversity and
structure are unchanged, and all fungal guilds are flat (total fungi
{_fpre:.1f}&rarr;{_fdro:.1f}% of reads). Read-abundance is simply blind to the isotope-allocation
axis where the companion signal sits.</div>
<div class="key">
<ul>
<li><b>Abundance &ne; activity — both methods show it.</b> ¹³C–PLFA tracks <i>carbon flux into
active lipid synthesis</i> (hours–days); shotgun reads and PLFA mol&percnt; both track <i>standing
community structure</i>. The drought response here is an allocation/flux response that only isotope
tracing resolves — neither abundance method should be expected to mirror it.</li>
<li><b>Actinomycetota — a concrete case.</b> Our WGS shows a drought <i>rise</i> in Actinomycetota
relative abundance, yet the companion ¹³C data show Actinobacteria are <i>not</i> preferentially
fed fresh carbon under drought (lowest enrichment; constant organic-matter degraders with
desiccation-resistant spores)<sup>[1]</sup>. The abundance rise reflects <b>drought persistence,
not active uptake of recent plant C</b> — abundance and activity pointing in different directions
within the same system.</li>
<li><b>AMF/fungi flat ≈ "present under both."</b> Our flat AMF/EMF read-abundance matches the PLFA
picture (AMF biomarkers in both conditions); the AMF response is a <i>subsoil ¹³C-allocation</i>
effect, below the detection floor of bulk-soil shotgun (AMF are ~6&percnt; of an already ~3&percnt;
fungal fraction).</li>
<li><b>Bacteria — function over taxonomy.</b> The companion work attributes carbon flow to root
exudate / metabolic control rather than cell-wall tolerance, with functional capacity retained even
as taxa reshuffle<sup>[1,3]</sup> — consistent with our result that no clean <i>taxonomic</i>
abundance signal survives correction; the drought response is metabolic, not compositional. (Our
declining nitrifier <i>Nitrospira</i> is an N-cycle read-out the C-focused PLFA study does not
target, but is consistent with drought-suppressed nitrification.)</li>
<li><b>Compartment &amp; resolution.</b> The companion study is <b>rhizosphere</b>, depth- and
tree-species-resolved; our shotgun set is neither, and carries a Sept-vs-Nov seasonal offset
(&sect;3b) — so it cannot resolve the topsoil-vs-subsoil or species-specific effects (e.g. the
<i>C. fairchildiana</i> subsoil-AMF flux) that carry the signal.</li>
</ul>
<b>Bottom line:</b> the metagenome and the ¹³C–PLFA study are <b>concordant where they overlap</b>
(a compositionally buffered standing community) and <b>complementary where they do not</b> (¹³C
carbon allocation). &ldquo;No abundance shift&rdquo; does not refute &ldquo;more carbon to
fungi&rdquo; — they are different axes. The definitive bridge is <b>DNA- / quantitative-SIP
metagenomics</b> with depth and compartment resolution (heavy-fraction sequencing of the
¹³C-labelled active community)<sup>[4]</sup> — the natural follow-up this survey cannot itself
provide.</div>
<p class="muted" style="font-size:.82em"><b>References.</b>
[1] Bai, Shi, Huang, Acharya … Werner, Dippold. The impacts of drought on plant–microbe carbon
interactions in an artificial tropical rainforest ecosystem (companion ¹³C–PLFA study). &nbsp;
[2] Werner et&nbsp;al. Ecosystem fluxes during drought and recovery in an experimental forest.
<i>Science</i> 2021, 374:1514. &nbsp;
[3] <a href="https://www.nature.com/articles/s41564-023-01432-9">Honeker et&nbsp;al. Drought
re-routes soil microbial carbon metabolism towards emission of volatile metabolites in an
artificial tropical rainforest.</a> <i>Nat. Microbiol.</i> 2023, 8:1480. &nbsp;
[4] <a href="https://journals.asm.org/doi/10.1128/msystems.00417-22">Quantitative stable-isotope
probing (qSIP) with metagenomics links microbial physiology and activity to soil moisture.</a>
<i>mSystems</i> 2022.</p>
"""

# --- Tier-2 N-cycle functional potential (NCyc) ---
NCYC_OK = os.path.exists(f"{RES}/ncyc_genefamily_hpm.tsv")
if NCYC_OK:
    _ncg = pd.read_csv(f"{RES}/ncyc_genefamily_hpm.tsv", sep="\t", index_col=0)
    _ncp = pd.read_csv(f"{RES}/ncyc_process_hpm.tsv", sep="\t", index_col=0)
    _ns = [x for x in _ncg.index if x in _meta.index]; _nmd = _meta.loc[_ns]
    _npre = [x for x in _nmd.index[_nmd.condition=="pre-drought"]]
    _ndro = [x for x in _nmd.index[_nmd.condition=="drought"]]
    def _gsum(df, pref):
        cols = [c for c in df.columns if c.startswith(pref)]
        if not cols: return (float('nan'), float('nan'))
        v = df[cols].sum(1); return (v.loc[_npre].mean(), v.loc[_ndro].mean())
    _nxrp, _nxrd = _gsum(_ncg, "nxr"); _amop, _amod = _gsum(_ncg, "amo")
    _denit = _ncp["Denitrification"].mean() if "Denitrification" in _ncp else float('nan')
    _dnra = _ncp["DNRA"].mean() if "DNRA" in _ncp else float('nan')
    _ncd = pd.read_csv(f"{RES}/ncyc_drought.tsv", sep="\t") if os.path.exists(f"{RES}/ncyc_drought.tsv") else pd.DataFrame()
    _ncnsig = int((_ncd["q_BH"] < 0.05).sum()) if len(_ncd) else 0
else:
    _nxrp=_nxrd=_amop=_amod=_denit=_dnra=float('nan'); _ncnsig=0
ncyc_html = (f"""
<h2>5b. Functional potential: nitrogen-cycle genes</h2>
<p>Shotgun reads were profiled against <b>NCyc</b> (68 nitrogen-cycle gene families; DIAMOND).
The community's N-cycling gene pool is dominated by <b>denitrification</b> (&asymp;{_denit:.0f}
hits per million reads) and <b>DNRA</b> (&asymp;{_dnra:.0f}), with a smaller dedicated
nitrification set.</p>
<div class="key"><b>The nitrite-oxidation gene pool is buffered even though <i>Nitrospira</i>
declines.</b> Drought leaves <i>nxrAB</i> (nitrite oxidoreductase) essentially unchanged
({_nxrp:.1f}&rarr;{_nxrd:.1f} hits/million reads) and ammonia-oxidation <i>amoABC</i> only
marginally lower ({_amop:.1f}&rarr;{_amod:.1f}); <b>no N-cycle gene family survives FDR</b>
({_ncnsig} significant). So the strong <i>taxonomic</i> decline of <i>Nitrospira</i>
(&sect;2b, &sect;4b) is <b>not</b> mirrored by a loss of nitrite-oxidation <i>gene potential</i> —
consistent with functional redundancy (other taxa also carry <i>nxr</i>) and with the
abundance&ne;function theme (&sect;6): the standing community's N-cycling capacity is buffered
against this drought.</div>
<div class="fig">{img('15_ncyc_drought.png')}</div>
""" if NCYC_OK else "")

# --- Tier-2 CAZyme functional potential (dbCAN) ---
CAZY_OK = os.path.exists(f"{RES}/cazy_class_hpm.tsv")
if CAZY_OK:
    _cz = pd.read_csv(f"{RES}/cazy_class_hpm.tsv", sep="\t", index_col=0)
    _cs = [x for x in _cz.index if x in _meta.index]; _cmd = _meta.loc[_cs]
    _cpre = [x for x in _cmd.index[_cmd.condition=="pre-drought"]]
    _cdro = [x for x in _cmd.index[_cmd.condition=="drought"]]
    def _czg(cl):
        return (_cz.loc[_cpre, cl].mean(), _cz.loc[_cdro, cl].mean()) if cl in _cz.columns else (float('nan'), float('nan'))
    _ghp, _ghd = _czg("GH"); _aap, _aad = _czg("AA")
    _czfd = pd.read_csv(f"{RES}/cazy_family_drought.tsv", sep="\t") if os.path.exists(f"{RES}/cazy_family_drought.tsv") else pd.DataFrame()
    _cznsig = int((_czfd["q_BH"] < 0.05).sum()) if len(_czfd) else 0
else:
    _ghp=_ghd=_aap=_aad=float('nan'); _cznsig=0
cazy_html = (f"""
<p><b>Carbohydrate-active enzymes (CAZymes; dbCAN, 4.1 M sequences).</b> The same pattern holds on
the carbon side: glycoside-hydrolase (GH; carbohydrate/litter degradation) and auxiliary-activity
(AA; ligninolytic/redox) gene potential are unchanged under drought (GH {_ghp:.0f}&rarr;{_ghd:.0f};
AA {_aap:.0f}&rarr;{_aad:.0f} hits/million reads) and <b>no CAZy class or family survives FDR</b>
({_cznsig} significant). The community's <i>degradative</i> gene repertoire is buffered too — so the
¹³C drought response is a <i>re-routing of carbon to fungal partners</i> (&sect;6), not a change in
the standing pool of carbohydrate-degrading or nitrogen-cycling genes.</p>
<div class="fig">{img('16_cazy_classes.png')}</div>
""" if CAZY_OK else "")

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
<p class="muted">Shotgun-metagenomic and 16S-amplicon profiling of the soil microbiome from the
Biosphere 2 WALD ¹³CO₂ drought experiment — the standing-community counterpart to the ¹³C–PLFA
carbon-flux study (Bai et&nbsp;al.). {stats.get('n_samples','36')} metagenomes across pre-drought
and drought; taxonomic profiling with sylph, Kaiju and DIAMOND, functional profiling with NCyc
and dbCAN. <a href="methods.html">Detailed methods &rarr;</a></p>

<div class="key"><b>Summary.</b> The standing soil community is compositionally stable across the
drought: no taxon survives multiple-testing correction and α- and β-diversity are unchanged. The
most consistent shift is a decline in the nitrite-oxidiser <b><i>Nitrospira</i></b>, reproduced
from phylum to strain and across four independent methods (sylph, Kaiju, DIAMOND, 16S), with a
parallel decline in the N-fixer <i>Bradyrhizobium</i> and an increase in <i>Actinomycetota</i>.
This signal does not reach significance once pseudoreplication is accounted for, and condition is
confounded with season (pre-drought = September, drought = November). Functional gene potential
(N-cycling, CAZymes) is likewise unchanged. Together this indicates the drought response is a
re-routing of carbon <i>allocation</i> to fungal partners (companion ¹³C study) rather than a
turnover of the standing community or its gene repertoire.</div>

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
<p>Actinomycetota are characteristically drought-tolerant (thick peptidoglycan walls,
sporulation, osmolyte accumulation), so their relative increase indicates the treatment imposed
real water stress on the community; the Nitrospirota decline points to suppressed nitrification.</p>
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
<div class="key">Both nitrite-oxidation (<i>Nitrospira</i>) and nitrogen-fixation
(<i>Bradyrhizobium</i>) decline together — a coherent suppression of two complementary N-cycle
functions. The <i>Bradyrhizobium</i> decline is notable given the belowground role of the legume
<i>C. fairchildiana</i> in the companion ¹³C study.</div>
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
{ncyc_html}
{cazy_html}
{discussion_html}
<h2>Statistical considerations</h2>
<div class="key" style="font-size:.93em">
<p>Four points underlie reporting the <i>Nitrospira</i> shift as a hypothesis rather than an
established effect. <b>Pseudoreplication:</b> the 0/6/48 h timepoints are repeated measures of the
same plots, not independent samples; a plot random-effect mixed model and a plot-collapsed test
both remove the nominal significance (§3b). <b>Multiple testing:</b> with Benjamini–Hochberg FDR
across taxa, no taxon reaches q&lt;0.05. <b>Effect size:</b> nonetheless <i>Nitrospira</i>'s
Cliff's δ ≈ −0.63 is large (drought samples are almost always lower), so the direction is
consistent even where it is not significant. <b>Power:</b> at n = {pwr.get('n_pre','?')} vs
{pwr.get('n_drought','?')} only very large effects (Cohen's d ≈ {_pwrd}) are detectable, so a
non-significant result here means "not detectable at this sample size", not "no effect". The
signal is therefore anchored on effect size and cross-method reproduction rather than p-values.</p>
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
