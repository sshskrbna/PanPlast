"""Self-contained HTML report."""
from __future__ import annotations

import base64
import html
import os
from pathlib import Path

import pandas as pd

from .analysis import CATEGORY_ORDER

CSS = """
:root{--bg:#f7f8fa;--card:#fff;--ink:#1f2329;--muted:#646b76;--line:#e2e5ea;--accent:#1f5f75;
--core:#1f5f75;--soft:#5fa8bf;--shell:#d9a441;--cloud:#b9bec6;--warn:#b5651d;--err:#b3261e}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#15181c;--card:#1d2127;--ink:#e7e9ec;
--muted:#9aa1ab;--line:#2e333a;--accent:#7cc0d6}}
:root[data-theme="dark"]{--bg:#15181c;--card:#1d2127;--ink:#e7e9ec;--muted:#9aa1ab;--line:#2e333a;--accent:#7cc0d6}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.5 -apple-system,
"Segoe UI",Inter,Roboto,sans-serif}main{max-width:1280px;margin:0 auto;padding:28px 16px 80px}
h1{font-size:26px;margin:0 0 4px}h2{font-size:19px;margin:36px 0 10px;padding-top:8px}
.sub{color:var(--muted);margin:0 0 20px}nav{display:flex;flex-wrap:wrap;gap:8px;margin:12px 0 4px}
nav a{color:var(--accent);text-decoration:none;border:1px solid var(--line);border-radius:999px;padding:3px 12px;
background:var(--card);font-size:13px}.cards{display:grid;grid-template-columns:repeat(auto-fill,minmax(170px,1fr));gap:12px}
.card{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.card .v{font-size:26px;font-weight:650;line-height:1.1}.card .v.sm{font-size:19px;white-space:nowrap}.card .k{color:var(--muted);font-size:12.5px;margin-top:4px}
.card .bar{height:4px;border-radius:2px;margin-bottom:10px}
figure{margin:0;background:#fff;border:1px solid var(--line);border-radius:12px;padding:10px;overflow:auto}
figure img{display:block;max-width:100%;height:auto}figcaption{color:var(--muted);font-size:12.5px;margin-top:8px}
figure img.tall{width:1800px;max-width:none}.dl{font-size:12.5px;color:var(--muted)}.dl a{color:var(--accent)}
.tbl{background:var(--card);border:1px solid var(--line);border-radius:12px;overflow:auto;max-height:560px}
table{border-collapse:collapse;width:100%;font-size:12.5px}th,td{padding:5px 9px;border-bottom:1px solid var(--line);
text-align:left;white-space:nowrap}th{position:sticky;top:0;background:var(--card);cursor:pointer;user-select:none;
font-weight:600}td.num{text-align:right;font-variant-numeric:tabular-nums}tr:hover td{background:rgba(127,140,160,.08)}
.filter{margin:6px 0 8px;padding:6px 10px;border:1px solid var(--line);border-radius:8px;background:var(--card);
color:var(--ink);width:min(360px,100%)}.pill{display:inline-block;padding:1px 8px;border-radius:999px;font-size:11.5px;
color:#fff}.core{background:var(--core)}.soft-core{background:var(--soft)}.shell{background:var(--shell)}
.cloud{background:var(--cloud);color:#222}.genes{display:flex;flex-wrap:wrap;gap:4px}.genes span{font-style:italic;
background:var(--card);border:1px solid var(--line);border-radius:6px;padding:0 6px;font-size:12.5px}
.gallery{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:10px}
.gallery a{display:block;background:#fff;border:1px solid var(--line);border-radius:10px;padding:6px;color:#333;
text-decoration:none;font-size:12px}.gallery img{width:100%;height:auto}.note{color:var(--muted);font-size:13px}
i{font-style:italic}
"""

JS = """
document.querySelectorAll('table.sortable').forEach(t=>{const th=t.querySelectorAll('th');
th.forEach((h,i)=>h.addEventListener('click',()=>{const tb=t.tBodies[0];const rows=[...tb.rows];
const asc=h.dataset.asc!=='1';th.forEach(x=>delete x.dataset.asc);h.dataset.asc=asc?'1':'0';
rows.sort((a,b)=>{const x=a.cells[i].dataset.v??a.cells[i].innerText,y=b.cells[i].dataset.v??b.cells[i].innerText;
const nx=parseFloat(x),ny=parseFloat(y);const r=(!isNaN(nx)&&!isNaN(ny))?nx-ny:x.localeCompare(y,'ru');
return asc?r:-r});rows.forEach(r=>tb.appendChild(r))}))});
document.querySelectorAll('input.filter').forEach(inp=>inp.addEventListener('input',()=>{
const t=document.getElementById(inp.dataset.t);const q=inp.value.toLowerCase();
[...t.tBodies[0].rows].forEach(r=>r.style.display=r.innerText.toLowerCase().includes(q)?'':'none')}));
"""


def _img(path, cls=""):
    data = base64.b64encode(Path(path).read_bytes()).decode()
    return f'<img class="{cls}" src="data:image/png;base64,{data}" alt="">'


def _links(fig, root):
    return " · ".join(f'<a href="{html.escape(os.path.relpath(p, root))}">{ext.upper()}</a>' for ext, p in fig.items())


def _table(df: pd.DataFrame, tid: str, italic_cols=("species", "gene"), pill_col="category"):
    if df is None or df.empty:
        return f'<p class="note">{_r("nodata")}</p>'
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for _, r in df.iterrows():
        cells = []
        for c in df.columns:
            v = r[c]
            if isinstance(v, float):
                txt = "" if pd.isna(v) else (f"{v:,.2f}".replace(",", " ") if abs(v) < 1000 else f"{v:,.0f}".replace(",", " "))
                cells.append(f'<td class="num" data-v="{"" if pd.isna(v) else v}">{txt}</td>')
            elif isinstance(v, (int,)) or (hasattr(v, "dtype") and str(getattr(v, "dtype", "")).startswith("int")):
                cells.append(f'<td class="num" data-v="{v}">{int(v):,}</td>'.replace(",", " "))
            elif c == pill_col:
                cells.append(f'<td><span class="pill {html.escape(str(v))}">{html.escape(str(v))}</span></td>')
            elif c in italic_cols:
                cells.append(f"<td><i>{html.escape(str(v))}</i></td>")
            else:
                cells.append(f"<td>{html.escape(str(v))}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<input class="filter" data-t="{tid}" placeholder="{_r("filter")}">'
            f'<div class="tbl"><table id="{tid}" class="sortable"><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div>')


R = {
    'nodata': 'No data.',
    'filter': 'Filter…',
    'title': 'PanPlast report',
    'h1': 'Chloroplast genome annotation',
    'sub': 'LSC/IR/SSC boundaries, junction genes, core and soft-core genes · {n} genomes · tables:',
    'nav_sum': 'Summary',
    'nav_junc': 'IR junctions',
    'nav_reg': 'Regions',
    'nav_matrix': 'Gene matrix',
    'nav_maps': 'Maps',
    'c_genomes': 'genomes',
    'c_len': 'genome length, bp',
    'c_lsc': 'LSC, bp',
    'c_ir': 'IR, bp',
    'c_ssc': 'SSC, bp',
    'c_pan': 'unique genes in total (pan)',
    'c_core_all': 'core: in all {n} genomes',
    'c_soft': 'soft core: {a}–{b} of genomes',
    'c_soft_total': 'soft core incl. core (≥{a})',
    'download': 'Download',
    'h_cats': 'Core and soft-core genes',
    'cap_pan': 'Left: genes per category, by type. Middle: how many genes occur in exactly k genomes. Right: growth of the pan and core sets with the number of genomes.',
    'h_cat_table': 'Category table',
    'h_junc': 'IR junctions (JLB, JSB, JSA, JLA)',
    'junc_note': 'JLB: LSC/IRb, JSB: IRb/SSC, JSA: SSC/IRa, JLA: IRa/LSC. Boundaries come from self-BLAST (the inverted repeat itself), not from the annotation.',
    'cap_junc': 'Genes at the junctions: forward strand above the bar, reverse strand below. A gene crossing the junction shows «bp left | bp right»; the others show their distance to the junction. Scroll sideways.',
    'h_junc_table': 'Junction genes',
    'h_reg': 'Region lengths',
    'cap_reg': 'LSC, IRb, SSC and IRa of every genome.',
    'h_sum_table': 'Per-genome summary',
    'h_matrix': 'Gene matrix',
    'cap_var': 'Variable genes only: absent somewhere or with a different copy number.',
    'full_matrix': 'Full matrix (all genes)',
    'cap_full': 'Colour: copy number; left band: gene category.',
    'h_maps': 'Circular maps',
    'maps_note': 'Click opens the PDF. Outer ring: forward-strand genes, inner ring: reverse strand; grey ring: GC content.',
    'h_qc': 'Quality control',
    'qc_note': 'Worth checking by hand: genomes missing soft-core genes (assembly or annotation problem), boundaries that differ from the GenBank annotation, unresolved tRNA names.',
    'h_ren': 'tRNAs renamed by homology',
    'ren_note': 'ARAGORN misreads the anticodon of some intron-containing tRNAs; the name is taken from the best homolog in the other genomes.',
}


def _r(key, **kw):
    s = R[key]
    return s.format(**kw) if kw else s


def write_report(path, summary, cats, junctions, qc, renames, figs, maps, labels, thr, root):
    root = Path(root)
    n = len(summary)
    cnt = cats.category.astype(str).value_counts().reindex(CATEGORY_ORDER, fill_value=0)
    core_t, soft_t, shell_t = thr
    soft_total = int((cats.percent >= 100 * soft_t).sum())

    def card(v, k, color=None):
        bar = f'<div class="bar" style="background:{color}"></div>' if color else ""
        cls = "v sm" if "–" in str(v) else "v"
        return f'<div class="card">{bar}<div class="{cls}">{v}</div><div class="k">{k}</div></div>'

    def rng(col):
        s = summary[col].dropna()
        return "—" if s.empty else (f"{int(s.min()):,}–{int(s.max()):,}".replace(",", " "))

    cards = "".join([
        card(n, _r("c_genomes")),
        card(rng("length_bp"), _r("c_len")),
        card(rng("LSC_bp"), _r("c_lsc")), card(rng("IRb_bp"), _r("c_ir")), card(rng("SSC_bp"), _r("c_ssc")),
        card(len(cats), _r("c_pan")),
        card(cnt["core"], _r("c_core_all", n=n) if core_t >= 1 else f"core (≥{core_t:.0%})", "var(--core)"),
        card(cnt["soft-core"], _r("c_soft", a=f"{soft_t:.0%}", b=f"{core_t:.0%}"), "var(--soft)"),
        card(soft_total, _r("c_soft_total", a=f"{soft_t:.0%}"), "linear-gradient(90deg,var(--core),var(--soft))"),
        card(cnt["shell"], f"shell ({shell_t:.0%}–{soft_t:.0%})", "var(--shell)"),
        card(cnt["cloud"], f"cloud (<{shell_t:.0%})", "var(--cloud)"),
    ])

    def genelist(cat):
        sub = cats[cats.category.astype(str) == cat]
        if sub.empty:
            return '<p class="note">—</p>'
        parts = []
        for tp in ("CDS", "tRNA", "rRNA"):
            g = sub[sub.type == tp]
            if len(g):
                spans = "".join(f'<span title="{r.n_genomes}/{n}">{html.escape(r.gene)}</span>' for r in g.itertuples())
                parts.append(f'<p class="note">{tp}: {len(g)}</p><div class="genes">{spans}</div>')
        return "".join(parts)

    cat_cols = ["gene", "type", "functional_group", "n_genomes", "percent", "category", "copies_min", "copies_max", "absent_in"]
    sum_cols = [c for c in summary.columns if c not in ("source",)]
    jcols = ["species", "genome", "junction", "gene", "strand", "relation", "distance_bp", "bp_left", "bp_right",
             "gene_start", "gene_end", "fragment"]
    jt = junctions[jcols] if len(junctions) else junctions

    def fig_block(key, caption, cls=""):
        f = figs.get(key)
        if not f:
            return ""
        return (f'<figure>{_img(f["png"], cls)}<figcaption>{caption}<br><span class="dl">{_r("download")}: '
                f'{_links(f, root)}</span></figcaption></figure>')

    gallery = ""
    if maps:
        items = "".join(
            f'<a href="{html.escape(os.path.relpath(m["pdf"], root))}"><img loading="lazy" '
            f'src="{html.escape(os.path.relpath(m["png"], root))}"><i>{html.escape(labels.get(g, g))}</i> · {g}</a>'
            for g, m in maps)
        gallery = f'<div class="gallery">{items}</div>'

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{_r("title")}</title>
<style>{CSS}</style></head><body><main>
<h1>{_r("h1")}</h1>
<p class="sub">{_r("sub", n=n)} <a href="panplast_results.xlsx">panplast_results.xlsx</a>, <a href="tables/">tables/*.tsv</a></p>
<nav><a href="#sum">{_r("nav_sum")}</a><a href="#cats">Core / soft core</a><a href="#junc">{_r("nav_junc")}</a>
<a href="#reg">{_r("nav_reg")}</a><a href="#matrix">{_r("nav_matrix")}</a><a href="#maps">{_r("nav_maps")}</a><a href="#qc">QC</a></nav>
<section id="sum"><h2>{_r("nav_sum")}</h2><div class="cards">{cards}</div></section>
<section id="cats"><h2>{_r("h_cats")}</h2>
{fig_block("pangenome", _r("cap_pan"))}
<h3>Core ({cnt["core"]})</h3>{genelist("core")}
<h3>Soft core ({cnt["soft-core"]})</h3>{genelist("soft-core")}
<h3>Shell ({cnt["shell"]})</h3>{genelist("shell")}
<h3>Cloud ({cnt["cloud"]})</h3>{genelist("cloud")}
<h3>{_r("h_cat_table")}</h3>{_table(cats[cat_cols], "t_cats")}</section>
<section id="junc"><h2>{_r("h_junc")}</h2>
<p class="note">{_r("junc_note")}</p>
{fig_block("junctions", _r("cap_junc"), "tall")}
<h3>{_r("h_junc_table")}</h3>{_table(jt, "t_junc")}</section>
<section id="reg"><h2>{_r("h_reg")}</h2>{fig_block("regions", _r("cap_reg"))}
<h3>{_r("h_sum_table")}</h3>{_table(summary[sum_cols], "t_sum")}</section>
<section id="matrix"><h2>{_r("h_matrix")}</h2>
{fig_block("heatmap_var", _r("cap_var"))}
<details><summary>{_r("full_matrix")}</summary>{fig_block("heatmap", _r("cap_full"))}</details></section>
<section id="maps"><h2>{_r("h_maps")}</h2><p class="note">{_r("maps_note")}</p>{gallery}</section>
<section id="qc"><h2>{_r("h_qc")}</h2>
<p class="note">{_r("qc_note")}</p>{_table(qc, "t_qc")}
<h3>{_r("h_ren")}</h3>
<p class="note">{_r("ren_note")}</p>
{_table(renames, "t_ren", italic_cols=())}</section>
</main><script>{JS}</script></body></html>"""
    Path(path).write_text(doc)
