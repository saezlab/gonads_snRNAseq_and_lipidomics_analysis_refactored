"""Build the panel trace and readable expanded-marker zonation update report."""
import argparse
import csv
import html
import json
import os
import re
from collections import Counter
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import resolve_recorded_path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether


def read_csv(path):
    return list(csv.DictReader(path.open(encoding="utf-8-sig"))) if path.exists() else []


def ascii_text(x):
    return str(x or "").translate(str.maketrans({"—":"-", "–":"-", "−":"-", "→":"->", "×":"x", "≤":"<=", "≥":">=", "±":"+/-", "’":"'", "“":'"', "”":'"', "⁻":"-", "¹":"1", "²":"2", "³":"3"}))


def status_text(value):
    labels={"Deferred—James":"Deferred to James", "Deferred - James":"Deferred to James",
            "PENDING_NEW_RESULTS":"Pending new results", "PENDING_MARKER_RULE_CONFIRMATION":"Methods detail for James to confirm",
            "PRESERVE_INDEPENDENT_MEASUREMENT":"Preserve independent measurement", "PRESERVE_ANALYSIS_PARAMETER":"Keep analysis setting",
            "UPDATE_R_PRESERVE_EXPERIMENTAL_KD_T":"Update correlation; preserve KD effect", "CORRECT_SEX_SPECIFIC_SAMPLE_UNIT":"Correct the sex-specific sample count",
            "REBUILD_FROM_EXPLICIT_CLASS_NUMERATORS_AND_DENOMINATORS":"Rebuild from explicit class fractions"}
    return labels.get(value,str(value or "").replace("_"," ").capitalize())


def panel_trace(root):
    rows = []
    def add(panel, status, workflow, files, changes, manuscript):
        existing = [str((root / f).resolve()) for f in files if (root / f).exists()]
        rows.append(dict(panel=panel, status=status, workflow=workflow, generated_files="; ".join(existing),source_tables="",
                         expected_files="; ".join(files), change_to_review=changes, manuscript_locations=manuscript))
    add("Figure4A", "Deferred - James", "Nucleus-level assignment/composition", [], "Nucleus counts, PP/PC proportions and composition p-values deferred by user; no pseudobulk substitution.", "P91; Figure4A legend/image")
    add("Figure4B", "Recomputed", "WF05 paired FCG model", ["fcg05/zonation_vs_sex_bias_scatter.pdf"], "Both zonation and sex axes; FDR counts and PP/PC splits; correct positive-sex-axis direction.", "P93; Figure4B legend; TableS10")
    add("Figure4C-E", "Preserved independent measurement", "Xenium", [], "No computational refit from these sci-RNA pseudobulks; retain independent counts/images.", "P93; TableS11")
    add("Figure4F-I", "Preserved independent measurement", "MALDI / RNA FISH", [], "No change to spatial-lipid assays or their numerical results.", "P95; TableS12")
    add("Figure5A", "Recomputed", "WF05 selection", ["fcg05/sex_specific_zonation_quadrant.pdf"], "Reoptimized threshold, selected genes, female/male coordinates and examples.", "P99/P101; Figure5A legend; TableS13")
    add("Figure5B", "Counts regenerated; schematic assembly", "WF05 selection", ["fcg05/sex_specific_zonation_summary_bar.pdf"], "Replace four programme counts; preserve schematic. Clarify assignment/overlap definition.", "P101/P103; Figure5B legend/image")
    add("Figure5C", "Recomputed; legend clarification", "WF05 network/class summaries", ["fcg05/lipid_class_zonation_jitter_vertical.pdf", "fcg05/lipid_class_zonation_proportions.pdf", "fcg05/Fig5C_four_context_associations.pdf"], "Separate legacy lipid-abundance-sex butterfly fractions from gene-programme sex/zone memberships; source denominators explicitly.", "P103; Figure5C legend; TableS14")
    add("Figure5D", "Recomputed network sources", "WF05 networks", ["fcg05/Fig5D_representative_networks.pdf"], "Refresh representative partner coefficients and eligibility. SIF/ATT files are the authoritative network source.", "P105; Figure5D legend; TableS14")
    add("Figure5E", "Preserved independent measurement", "Mtmr11 knockdown RNA efficacy", [], "KD efficacy values, four replicate design and original experimental p-value retained.", "Figure5E legend/image")
    add("Figure5F", "Recomputed", "WF07 zonation validation", ["fcg07/mtmr11_zonation_validation_scatter.pdf", "fcg07/Fig5F_Mtmr11_male_periportal.pdf"], "Update in-vivo r, ULM/p and regression against unchanged KD lipid effects; correct male-only animal count.", "P107; Figure5F legend/image; TableS15")
    add("Figure6A", "Recomputed", "WF14 models -> WF15 -> WF16 -> layouts", ["cf453_rewiring/panels/Fig6A_anchor_gene_shifts.pdf", "cf453_layouts/panels/Fig6A_option1_heatmap.pdf", "cf453_layouts/panels/Fig6A_option2_paired_dots.pdf"], "Recompute anchors, raw effects and direct-test marks; both layouts use the new set and observed counts.", "Working Figure6 draft; absent from inspected five-figure manuscript")
    add("Figure6B", "Recomputed", "WF16 class associations", ["cf453_rewiring/panels/Fig6B_lipid_class_shift.pdf"], "New selected programmes, same paired SD/WD lipid universe, fresh 999-draw uncertainty.", "Working Figure6 draft")
    add("Figure6C", "Recomputed", "WF16 gained/lost examples", ["cf453_rewiring/panels/Fig6C_network_gain_loss.pdf"], "New correlations, selected partners and omission/bootstrap checks. Exploratory threshold-defined examples.", "Working Figure6 draft; final assembly letters provisional")
    add("FigureS4B-C", "Deferred - James", "Nucleus-level score/marker QC", [], "Update old two-marker gating illustration and marker-expression/count displays later.", "FigureS4B-C legends/images")
    add("FigureS4A/D/E", "Preserved independent measurement", "Xenium / MALDI / FISH", [], "Independent spatial images and measurements retained.", "FigureS4 legends; TablesS11/S12")
    add("FigureS5A", "Recomputed", "WF05 threshold optimization", ["fcg05/optimal_threshold_detection.pdf"], "Replace operating-point label and curve; retain fixed-historical-threshold sensitivity separately.", "P101; FigureS5A legend/image")
    add("FigureS5B", "Image preserved; interpretation checked", "Independent Aox3 Xenium", [], "Image unchanged; confirm gene pattern and context against refreshed RNA.", "P101; FigureS5B")
    add("FigureS5C", "Preserved independent measurement", "Cultured hepatocyte RNA", [], "Zonation-marker expression after culture/KD is independent of these inputs.", "P107; FigureS5C")
    add("TablesS10/S13/S14/S15", "Recomputed result sources", "WF05/WF07", ["fcg05/ttop_full_long.csv", "fcg05/sex_specific_zonation_genes_summary.csv", "fcg07/zonation_network_validation_scores.csv"], "Replace DE, selected gene, gene-lipid-TF network and validation tables from the new versioned sources.", "Supplementary table list")
    add("GDX supporting analysis", "New paired support; no main-panel replacement", "Repaired four-group GDX branch", ["gdx_zonation/GDX_zonation_interaction_comparison.pdf"], "Old/new inputs fitted with the same explicit paired model; not the incompatible legacy hormone-replacement analysis.", "No verified current main-figure consumer")
    add("Figures1-3 / S1-S3", "Preserved", "Whole-hepatocyte / lipid / MOFA / KO", [], "No zone-input dependency. Figure3 editorial/KO network-version decisions remain separate.", "Unchanged manuscript sections")
    add("WF12 zonal comparison", "Historical snapshot retained", "Conditional cross-cohort rebase", [], "Must refresh FCG reference and zone RNA together before reuse with the new definitions. Not required for standalone Figure6.", "Optional supporting comparison, not a refreshed main panel")
    add("WF13 CF459 comparison", "Held pending equivalent CF459 zones", "Three-cohort robustness", [], "No expanded-marker CF459 input supplied; do not mix compartment definitions.", "Historical supporting comparison retained")
    # Merge any production-agent panel manifests so the trace follows actual exports.
    manifests = [p for p in root.rglob("*manifest*.csv") if "inputs" not in p.parts and "report" not in p.parts]
    for path in manifests:
        for r in read_csv(path):
            panel = str(r.get("panel_id") or r.get("panel") or r.get("figure_panel") or "").replace(" ", "")
            canonical=lambda x:re.sub(r"^(Figure|Fig)","",str(x).replace(" ","")).split("_",1)[0]
            for item in rows:
                if canonical(panel)==canonical(item["panel"]):
                    f = r.get("generated_file") or r.get("file") or r.get("path") or r.get("output_file")
                    if f:
                        q = Path(f)
                        if q.is_absolute(): q = resolve_recorded_path(q)
                        if not q.is_absolute():
                            q = root / q if (root / q).exists() else path.parent / q
                        if q.exists() and str(q.resolve()) not in item["generated_files"]:
                            item["generated_files"] = "; ".join(filter(None,[item["generated_files"],str(q.resolve())]))
                    source=r.get("source_table")
                    if source:
                        q=Path(source)
                        if q.is_absolute():q=resolve_recorded_path(q)
                        if not q.is_absolute():q=root/q if (root/q).exists() else path.parent/q
                        if q.exists() and str(q.resolve()) not in item["source_tables"]:
                            item["source_tables"]="; ".join(filter(None,[item["source_tables"],str(q.resolve())]))
    table_rows=read_csv(root/'tables/table_manifest.csv')
    for item in rows:
        if item['panel']=='TablesS10/S13/S14/S15':
            item['source_tables']='; '.join(str(resolve_recorded_path(r['file'])) for r in table_rows if resolve_recorded_path(r['file']).exists())
    target = root / "report/panel_trace.csv"
    with target.open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows


def build(root, draft=False):
    root = Path(root).resolve(); report = root / "report"; report.mkdir(exist_ok=True)
    claims = read_csv(report / "manuscript_numeric_inventory.csv")
    if not claims:
        raise RuntimeError("Manuscript numerical inventory is required")
    pending = [r for r in claims if r.get("status") == "PENDING_NEW_RESULTS"]
    if pending and not draft:
        raise RuntimeError("Unresolved production values: " + ", ".join(r["metric_key"] for r in pending))
    panels = panel_trace(root)
    cm = read_csv(root / "comparisons/cf453_model_summary.csv")
    cg = read_csv(root / "comparisons/fcg_selected_gene_counts.csv")
    f4 = read_csv(root / "comparisons/Fig4B_FDR_counts_old_vs_new.csv")
    bykey = {r["metric_key"]:r for r in claims}
    def current(key,default="See numerical ledger"):
        r=bykey.get(key,{})
        return r.get("new_value") or r.get("proposed_replacement") or default
    headline = [
      "The supplied inputs contain paired periportal/pericentral RNA profiles from 43 Four Core Genotypes mice, 22 gonadectomy/sham mice and 28 wild-type standard/Western-diet mice. Both zones remain paired by animal. Whole-liver lipid measurements are unchanged.",
      "Figure4A and S4B-C are deferred to James by explicit user instruction. Their nucleus counts, proportions and composition p-values are not estimated from pseudobulks.",
      "FCG gene-level models, selected programmes, gene-lipid networks, class associations and zonation-network validation were recomputed from the supplied expanded-marker pseudobulks.",
      "The diet-study chain was rebuilt through the direct-interaction models, selected-programme associations and both alternative Figure6A layouts. The lipid measurements and independent experimental results remain fixed.",
      "The report distinguishes marker-input changes, the consequence of reoptimizing the operating threshold, and pre-existing legend/denominator/sample-unit issues. The manuscript itself has not been edited."
    ]
    if f4:
        old=next(r for r in f4 if r["version"]=="historical");new=next(r for r in f4 if r["version"]=="revised")
        headline.append(f"Figure4B: {old['dual_bias']} of {old['sex_biased']} -> {new['dual_bias']} of {new['sex_biased']} genes have both gonadal and zonation FDR<=0.05; PC {old['pericentral']} -> {new['pericentral']}, PP {old['periportal']} -> {new['periportal']}.")
    if cg:
        counts={v:sum(int(r["n_genes"]) for r in cg if r["version"]==v) for v in {r["version"] for r in cg}}
        headline.append(f"Figure5A: {counts.get('historical')} selected genes at old threshold 5.4; {counts.get('revised_fixed_historical_threshold')} with updated inputs at that same cut; {counts.get('revised_reoptimized')} after reoptimization to threshold {current('interaction_threshold')}. Update the threshold and gene counts together.")
    if cm:
        old=next(r for r in cm if r["version"]=="old");new=next(r for r in cm if r["version"]!="old")
        headline.append(f"Diet-study eligible genes {old['eligible_genes']} -> {new['eligible_genes']}; paired-mouse correlation {float(old['rho']):.4f} -> {float(new['rho']):.4f}. Whole-liver lipid interaction discoveries remain {new['lipid_interaction_q05']}.")
    headline.append("The refreshed Figure6 uses 35 unique genes from 40 anchor selections, with ten direct sex-dependent diet-interaction marks. Its four network representative genes remain the same; six of the 24 selected partners change. Association-level conclusions remain exploratory: no correlation passes the declared BH family and all 120 class bootstrap ranges include zero.")
    headline.append("Mtmr11 male-periportal ULM: " + bykey.get("Mtmr11_ULM",{}).get("old_value","-7.21") + " -> " + current("Mtmr11_ULM") + "; p: " + current("Mtmr11_ULM_p") + ". Use the explicitly prespecified male-PP context, separately from adaptive context selection.")
    headline.append("Figure5 associated lipid species: "+bykey.get("unique_associated_lipids",{}).get("old_value","442")+" -> "+current("unique_associated_lipids")+". The class and network sections distinguish primary reoptimization from the fixed-historical-threshold sensitivity.")
    decisions = [
      "Use the reoptimized gene set with its new threshold and retain the fixed 5.4 sensitivity. Do not attribute the whole difference in selected-gene counts to large gene-level changes.",
      "Figure5C has two distinct encodings: the faithful butterfly uses whole-liver lipid sex bias, whereas the explicit four-context view uses the sex/zone of the correlated gene programme. Match each legend to its actual source definition.",
      "Gene-lipid edge fractions, unique lipid-species counts, shared-category counts and mean partners per gene have different denominators. The numerical ledger states which quantity each replacement represents.",
      "Correct Figure4B's sex-axis direction to match the plotted Testes-minus-Ovaries contrast. Correct Figure5F's male-only n using the matched animal roster; the pooled n=43 is not its sampling unit.",
      "Replace the old two-marker Methods description with the supplied multi-gene UCell scoring procedure. Exact per-cohort retained markers and nucleus composition outputs remain James's provenance/figure tasks.",
      "Figure5E, independent Xenium/MALDI/FISH measurements, whole-hepatocyte models and the Figure3 KO analysis are preserved. Historical WF12 zonal/WF13 results must not be mixed with the new compartment definitions."
    ]
    source_paths = ["inputs/import_manifest.csv","inputs/cohort_qc.csv","comparisons/fcg_historical_model_exact_replay.csv",
                    "comparisons/Fig4B_claim_denominator_definition_audit.csv","comparisons/fcg_contrast_summary_old_vs_new.csv",
                    "comparisons/fcg_selected_gene_counts.csv","comparisons/cf453_model_summary.csv","comparisons/cf453_chain_checks.csv",
                    "comparisons/gdx_summary.csv","QA/historical_sha256.csv"]
    # Editable self-contained HTML with full literal claim text and result links.
    esc=lambda s:html.escape(str(s or ""))
    def htable(headers, records):
        return '<table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in records)+'</tbody></table>'
    def artifact_links(item):
        paths=[Path(x) for x in item['generated_files'].split('; ') if x and Path(x).suffix.lower()=='.pdf']
        sources=[Path(x) for x in item['source_tables'].split('; ') if x]
        return '<br>'.join('<a href="'+esc(os.path.relpath(q,report))+'">'+esc(q.name)+'</a>' for q in paths+sources)
    panel_html='<table><thead><tr>'+''.join('<th>'+x+'</th>' for x in ['Panel','Status / analysis','Change and manuscript location','Panel and source files'])+'</tr></thead><tbody>'+''.join('<tr><td>'+esc(r['panel'])+'</td><td>'+esc(r['status']+' / '+r['workflow'])+'</td><td>'+esc(r['change_to_review']+' '+r['manuscript_locations'])+'</td><td>'+artifact_links(r)+'</td></tr>' for r in panels)+'</tbody></table>'
    h=['<!doctype html><html><head><meta charset="utf-8"><title>Expanded marker liver zonation update</title><style>body{font:15px/1.5 -apple-system,BlinkMacSystemFont,Arial,sans-serif;color:#202631;max-width:1250px;margin:35px auto;padding:0 25px}h1,h2{line-height:1.2}h1{font-size:30px}h2{margin-top:36px;font-size:23px}table{border-collapse:collapse;width:100%;margin:20px 0;font-size:13px}td,th{border:1px solid #d7dce2;padding:9px;vertical-align:top;text-align:left;overflow-wrap:anywhere}th{background:#e9eef4}tr:nth-child(even){background:#f8f9fb}.summary{background:#edf4f8;padding:16px 23px;border-left:4px solid #32647e}.claim{border-top:1px solid #ddd;padding:15px 0}code{font-size:12px;overflow-wrap:anywhere}a{color:#175f85}.status{font-weight:bold}</style></head><body>',
       '<h1>Expanded marker liver zonation update</h1><p>16 November 2026. Pseudobulk-dependent rerun and manuscript replacement audit.</p><div class="summary">'+''.join('<p>'+esc(x)+'</p>' for x in headline)+'</div>',
       '<h2>Figure panel trace</h2>',panel_html,
       '<h2>Numerical replacement ledger</h2><p>Identifiers such as body.p0093 mean top-level paragraph 93 in the frozen local manuscript, using one-based numbering. The complete CSV retains source paths, literal text, replacement values and notes. Unchanged experimental values and deferred James items are included deliberately.</p>',
       '<p><a href="manuscript_numeric_inventory.csv">Complete claim inventory CSV</a> | <a href="new_replacements.csv">Old-to-new replacements CSV</a> | <a href="panel_trace.csv">Panel trace CSV</a></p>']
    for r in claims:
        h.append('<div class="claim"><b>'+esc(r['metric_key'])+'</b> <span class="status">'+esc(status_text(r.get('status')))+'</span><p>'+esc(r.get('paragraph_id'))+' / '+esc(r.get('section_or_panel'))+'</p><p><b>Old:</b> '+esc(r.get('old_value'))+'<br><b>New:</b> '+esc(r.get('new_value') or r.get('proposed_replacement') or 'No replacement generated; see status')+'</p><p>'+esc(r.get('note'))+'</p><details><summary>Original manuscript wording and result source</summary><p>'+esc(r.get('exact_text'))+'</p><code>'+esc(r.get('new_result_source') or r.get('source'))+'</code></details></div>')
    h.extend(['<h2>Interpretation and author actions</h2>',''.join('<p>'+esc(x)+'</p>' for x in decisions),'<h2>Reproducibility files</h2><ul>'])
    for f in source_paths:
        if (root/f).exists(): h.append('<li><a href="../'+esc(f)+'">'+esc(f)+'</a></li>')
    h.append('</ul><p>Full versioned outputs and analysis-specific QA remain under this run directory. Historical input/result SHA256 preservation is checked separately.</p></body></html>')
    (report/'zonation_update_report.html').write_text('\n'.join(h))
    # Paginated PDF: readable panel trace and complete concise numerical ledger.
    styles=getSampleStyleSheet()
    styles.add(ParagraphStyle(name='SmallCellZU',fontName='Helvetica',fontSize=7.8,leading=10,spaceAfter=2))
    styles.add(ParagraphStyle(name='BodyZU',fontName='Helvetica',fontSize=9.2,leading=13,spaceAfter=7))
    styles.add(ParagraphStyle(name='SubZU',fontName='Helvetica-Bold',fontSize=12,leading=15,spaceBefore=12,spaceAfter=7))
    def p(x,style='BodyZU'):
        return Paragraph(html.escape(ascii_text(x)).replace('\n','<br/>'),styles[style])
    story=[p('Expanded marker liver zonation update','Title'),p('16 November 2026 | Pseudobulk rerun and manuscript replacement audit'+(' | DRAFT' if draft else ''))]
    for x in headline: story.append(p(x))
    story.append(p('Figure panel trace','SubZU'))
    def table(headers,values,widths):
        data=[[p(x,'SmallCellZU') for x in headers]]+[[p(x,'SmallCellZU') for x in row] for row in values]
        t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dfe9f0')),('VALIGN',(0,0),(-1,-1),'TOP'),('GRID',(0,0),(-1,-1),.35,colors.HexColor('#cbd3da')),('LEFTPADDING',(0,0),(-1,-1),5),('RIGHTPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f6f8fa')])]))
        return t
    story.append(table(['Panel','Status / analysis','Changes to review'],[[r['panel'],r['status']+'\n'+r['workflow'],r['change_to_review']+'\n'+r['manuscript_locations']+'\nFiles: '+', '.join(Path(f).name for f in r['generated_files'].split('; ') if f and Path(f).suffix.lower()=='.pdf')] for r in panels],[65,145,305]))
    story.append(p('Direct panel links','SubZU'))
    for item in panels:
        links=[Path(f) for f in item['generated_files'].split('; ') if f and Path(f).suffix.lower()=='.pdf']
        if links:
            label=html.escape(item['panel'])+': '+ ' | '.join('<link href="'+html.escape(q.as_uri(),quote=True)+'" color="#175f85">'+html.escape(q.name)+'</link>' for q in links)
            story.append(Paragraph(label,styles['SmallCellZU']))
    story.extend([PageBreak(),p('Manuscript numerical replacements','SubZU'),p(f'All {len(claims)} inventoried claims are listed below. body.p0093 means top-level manuscript paragraph 93 (one-based numbering in the frozen source). Exact text and machine-readable result paths are retained in manuscript_numeric_inventory.csv and the HTML report. Blank replacements are intentional for preserved or deferred items.')])
    for r in claims:
        label=r['metric_key'].replace('_',' ')
        status_raw=r.get('status','');status=status_text(status_raw)
        value=r.get('new_value') or r.get('proposed_replacement') or ('Retain existing value' if status_raw.startswith('PRESERVE') else 'See action / status')
        story.append(KeepTogether([p(r.get('paragraph_id','')+' | '+r.get('section_or_panel','')+' | '+label,'SubZU'),
          table(['Original manuscript','Updated result or action','Status'],[[r.get('old_value',''),value,status]],[145,250,120]),
          p(r.get('note',''),'SmallCellZU'),Spacer(1,5)]))
    story.extend([PageBreak(),p('Interpretation and author actions','SubZU')])
    story.extend(p(x) for x in decisions)
    story.append(p('Reproducibility and files','SubZU'))
    tm=root/'tables/table_manifest.csv'
    if tm.exists():story.append(Paragraph('<link href="'+html.escape(tm.as_uri(),quote=True)+'" color="#175f85">Updated supplementary source table manifest</link>',styles['BodyZU']))
    story.extend(p(x,'SmallCellZU') for x in source_paths if (root/x).exists())
    story.append(p('Full artifact paths and version-specific source tables are in panel_trace.csv. The report does not alter the manuscript or the preserved historical analyses.'))
    pdf=report/'zonation_update_report.pdf'
    def footer(canvas,doc):
        canvas.saveState();canvas.setFont('Helvetica',8);canvas.setFillColor(colors.HexColor('#606b76'))
        canvas.drawString(14*mm,10*mm,'Expanded marker zonation update | 16 November 2026')
        canvas.drawRightString(A4[0]-14*mm,10*mm,str(doc.page));canvas.restoreState()
    SimpleDocTemplate(str(pdf),pagesize=A4,leftMargin=14*mm,rightMargin=14*mm,topMargin=14*mm,bottomMargin=17*mm).build(story,onFirstPage=footer,onLaterPages=footer)
    summary=dict(report_pdf=str(pdf),report_html=str(report/'zonation_update_report.html'),claims=len(claims),panels=len(panels),statuses=dict(Counter(r.get('status','') for r in claims)),draft=draft)
    (report/'report_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='results/zonation_marker_update_20261116')
    parser.add_argument('--draft',action='store_true')
    args=parser.parse_args();build(args.output,args.draft)
