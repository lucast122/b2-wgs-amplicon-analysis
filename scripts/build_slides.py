#!/usr/bin/env python3
"""Build a self-contained 10-slide HTML deck (figures embedded as base64; arrow-key / click
navigation) summarising the metagenomics report for a lab/colleague presentation."""
import base64, os
import pandas as pd

A = "/mnt/disk4/timo/gbi/b2/analysis"; FIG = f"{A}/figs"; RES = f"{A}/results"
OUT = f"{A}/GBI_slides.html"

def img(name, h="62vh"):
    p = f"{FIG}/{name}"
    if not os.path.exists(p): return ""
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f'<img src="data:image/png;base64,{b}" style="max-height:{h};max-width:88%;object-fit:contain;">'

def kv(fn):
    d = {}; p = f"{RES}/{fn}"
    if os.path.exists(p):
        for line in open(p):
            if "\t" in line:
                k, v = line.strip().split("\t", 1); d[k] = v
    return d

# --- numbers ---
stats = kv("stats_summary.txt"); n = stats.get("n_samples", "39")
dchk = kv("diamond_checks.txt"); rho = float(dchk.get("median_spearman_vs_kaiju", "nan"))
gd = pd.read_csv(f"{RES}/genus_drought.tsv", sep="\t")
def g(name, col):
    r = gd[gd.genus == name]
    return float(r.iloc[0][col]) if len(r) else float("nan")
nit_pre, nit_dro, nit_p, nit_d = g("Nitrospira_D","pre"), g("Nitrospira_D","drought"), g("Nitrospira_D","p"), g("Nitrospira_D","cliffs_delta")
pwr = kv("power_analysis.txt"); dmin = pwr.get("min_detectable_Cohens_d_80pct","large")

SL = []
def slide(html): SL.append(html)

# 1 title
slide(f"""<div class="center">
<div class="kicker">Biosphere 2 · WALD drought campaign</div>
<h1>Soil microbiome response to drought<br>in a tropical rainforest</h1>
<p class="sub">Shotgun-metagenomic & 16S profiling — the standing-community counterpart to the
¹³C–PLFA carbon-flux study (Bai et&nbsp;al.)</p>
<p class="meta">{n} metagenomes · sylph · Kaiju · DIAMOND · NCyc · dbCAN · GTDB r232</p></div>""")

# 2 question / design
slide(f"""<h2>The question</h2>
<div class="cols"><div>
<p class="lead">The ¹³C–PLFA study shows drought <b>re-routes carbon</b> to mycorrhizal & saprotrophic fungi.
Does that leave a fingerprint in <b>who is present</b> and <b>what genes they carry</b>?</p>
<ul>
<li>{n} soil metagenomes, pre-drought (Sept) vs drought (Nov)</li>
<li>Taxonomy: sylph (genome containment) + Kaiju & DIAMOND (protein)</li>
<li>Independent 16S amplicon survey (328 samples)</li>
<li>Function: N-cycle (NCyc) and CAZymes (dbCAN)</li>
</ul></div>
<div class="big">🌧️→🍄</div></div>""")

# 3 composition
slide(f"""<h2>1 · A stable, bacteria-dominated, largely uncultured community</h2>
<div class="fig">{img('04_composition.png')}</div>
<p class="cap">~95% bacteria (Pseudomonadota ≈ half), ~16% ammonia-oxidising archaea. Mean composition
is near-identical pre- vs drought — any signal is in <b>specific taxa</b>, not the bulk.</p>""")

# 4 drought response
slide(f"""<h2>2 · The clearest signal: nitrite-oxidising <i>Nitrospira</i> declines</h2>
<div class="row">{img('09_nitrospira_ranks.png','52vh')}{img('05_drought_response.png','52vh')}</div>
<p class="cap"><i>Nitrospira</i> drops consistently from phylum → genus → strain
({nit_pre:.1f}→{nit_dro:.1f}%, p={nit_p:.3f}, Cliff's δ={nit_d:.2f}, large); drought-tolerant
<b>Actinomycetota</b> rise. A coherent suppression of nitrification.</p>""")

# 5 honesty / robustness
slide(f"""<h2>3 · …but the standing community is statistically buffered</h2>
<div class="cols"><div>
<ul class="big-ul">
<li><b>No taxon survives</b> multiple-testing correction (FDR)</li>
<li>α- and β-diversity <b>unchanged</b> (PERMANOVA n.s.)</li>
<li>Timepoints are pseudoreplicates → plot-level models <b>remove</b> significance</li>
<li>Underpowered (only d ≈ {dmin} detectable) and <b>confounded with season</b> (Sept vs Nov)</li>
</ul></div><div class="quote">
The shift is reported as a <b>hypothesis</b>, anchored on effect size and cross-method
reproduction — not on a single p-value.</div></div>""")

# 6 cross-validation
slide(f"""<h2>4 · Four independent methods agree</h2>
<div class="fig">{img('14_diamond_vs_kaiju.png','56vh')}</div>
<p class="cap">Two protein aligners (Kaiju, DIAMOND) give concordant profiles (Spearman ρ={rho:.2f})
and the same dominant taxa — <i>Bradyrhizobium, Nitrospira, Luteitalea, Mesorhizobium</i>. The
<i>Nitrospira</i> decline reproduces across <b>sylph, Kaiju, DIAMOND and 16S</b>.</p>""")

# 7 mycorrhizal
slide(f"""<h2>5 · Mycorrhizal fungi recovered (custom Kaiju fungal DB)</h2>
<div class="fig">{img('12_amf_emf_genera.png','56vh')}</div>
<p class="cap">Arbuscular mycorrhizal fungi are present and <i>Rhizophagus</i>-dominated — the
symbionts central to the manuscript's belowground carbon story, which SSU surveys and
genome-containment miss. No AMF/EMF genus shifts significantly.</p>""")

# 8 functional
slide(f"""<h2>6 · Functional gene potential is buffered too</h2>
<div class="row">{img('15_ncyc_drought.png','50vh')}{img('16_cazy_classes.png','50vh')}</div>
<p class="cap">Nitrite-oxidation genes (<i>nxr</i>) are flat <b>despite</b> the <i>Nitrospira</i>
decline (functional redundancy); carbohydrate-degradation (CAZy GH/AA) is flat too. Neither the
N-cycling nor the C-degrading repertoire turns over.</p>""")

# 9 reconciliation
slide(f"""<h2>7 · Reconciling with the ¹³C–PLFA study</h2>
<div class="cols"><div>
<p class="lead">Concordant where they overlap, complementary where they don't.</p>
<ul>
<li>PLFA <b>composition</b> is also buffered (only dried topsoil shifts) — the ¹³C signal is in
<b>carbon allocation</b>, not standing biomass</li>
<li>Our DNA abundance measures exactly that buffered axis</li>
<li>Example: Actinomycetota <b>abundance ↑</b> here, yet ¹³C shows they are <b>not</b> fed fresh
carbon — abundance ≠ activity</li>
</ul></div><div class="quote">
"No abundance shift" does not refute "more carbon to fungi" — they are <b>different axes</b> of
the same system.</div></div>""")

# 10 conclusions
slide(f"""<h2>Conclusions & next steps</h2>
<ul class="big-ul">
<li>Drought leaves the <b>standing soil community and its gene repertoire essentially unchanged</b></li>
<li>Most consistent signal: a <b>nitrification (Nitrospira) decline</b> — reproducible across 4 methods, large effect, but not significant after correction / season-confounded</li>
<li>Consistent with the ¹³C story: drought acts on <b>carbon allocation & activity</b>, not community turnover</li>
<li><b>In progress:</b> genome-resolved MAGs (assembly + binning of all 39) → strain-level & functional confirmation; ¹³C per-plot integration</li>
</ul>
<p class="meta" style="margin-top:30px">Full report + figures + tables: GitHub Pages · companion to Bai et&nbsp;al.</p>""")

slides_html = "".join(f'<section class="slide">{s}<div class="pg">{i+1} / {len(SL)}</div></section>' for i, s in enumerate(SL))
doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>B2 drought metagenomics — slides</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#1a2233;background:#0b1f33}}
.slide{{display:none;width:100vw;height:100vh;padding:5vh 6vw;background:#fff;flex-direction:column;
align-items:center;justify-content:flex-start;position:fixed;top:0;left:0;overflow:hidden}}
.slide.on{{display:flex}}
h1{{font-size:3.0em;line-height:1.12;color:#10233b;text-align:center}}
h2{{font-size:1.7em;color:#10233b;border-bottom:3px solid #1f77b4;padding-bottom:8px;margin-bottom:2.2vh;align-self:stretch}}
.center{{height:100%;display:flex;flex-direction:column;align-items:center;justify-content:center;text-align:center}}
.kicker{{color:#1f77b4;font-weight:700;letter-spacing:.12em;text-transform:uppercase;font-size:1em;margin-bottom:18px}}
.sub{{font-size:1.35em;color:#33485f;margin-top:22px;max-width:60%}}
.meta{{color:#7088a0;margin-top:26px;font-size:1.05em}}
.lead{{font-size:1.5em;line-height:1.4;margin-bottom:18px}}
ul{{font-size:1.35em;line-height:1.7;margin-left:1.1em;max-width:46vw}} .big-ul li{{font-size:1.05em;margin-bottom:10px}}
.cols{{display:flex;gap:5vw;width:100%;align-items:center;justify-content:space-between;flex:1}}
.row{{display:flex;gap:2vw;width:100%;align-items:center;justify-content:center;flex:1}}
.fig{{flex:1;display:flex;align-items:center;justify-content:center;width:100%}}
.cap{{font-size:1.2em;color:#33485f;text-align:center;max-width:84%;margin-top:2vh;line-height:1.45}}
.big{{font-size:9em}} .big-ul{{list-style:none;margin:0}} .big-ul li{{padding-left:1.4em;position:relative}}
.big-ul li:before{{content:"▸";color:#1f77b4;position:absolute;left:0}}
.quote{{font-size:1.5em;font-style:italic;color:#10233b;border-left:5px solid #1f77b4;padding:14px 22px;background:#f3f7fb;max-width:34vw;line-height:1.4}}
.pg{{position:absolute;bottom:2.5vh;right:3vw;color:#9fb0c4;font-size:.95em}}
i{{font-style:italic}}
</style></head><body>
{slides_html}
<script>
let i=0; const s=document.querySelectorAll('.slide');
function show(n){{i=(n+s.length)%s.length; s.forEach((x,k)=>x.classList.toggle('on',k===i));}}
document.addEventListener('keydown',e=>{{if(['ArrowRight',' ','PageDown'].includes(e.key))show(i+1);
else if(['ArrowLeft','PageUp'].includes(e.key))show(i-1);}});
document.addEventListener('click',e=>show(i+1));
show(0);
</script></body></html>"""
open(OUT, "w").write(doc)
print(f"Wrote {OUT} ({len(SL)} slides, {os.path.getsize(OUT)//1024} KB)")
