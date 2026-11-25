"""Verify original-file preservation and package the versioned report and panels."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'maintenance'))
from repository_paths import resolve_recorded_path


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig'))) if path.exists() else []


def package(root):
    root=Path(root).resolve()
    before=read(root/'QA/historical_sha256.csv')
    assert before, 'Historical snapshot required'
    preserved=[]
    for row in before:
        p=resolve_recorded_path(row['path']);ok=p.exists() and sha(p)==row['sha256']
        preserved.append(dict(path=row['path'],resolved_path=str(p),sha256_before=row['sha256'],unchanged=ok))
    assert all(x['unchanged'] for x in preserved), 'Historical artifact content changed'
    with (root/'QA/historical_preservation_checks.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(preserved[0]));w.writeheader();w.writerows(preserved)
    report_manifest=json.loads((root/'report/report_manifest.json').read_text())
    assert not report_manifest['draft'], 'Do not package a draft report'
    selected=set()
    for item in read(root/'report/panel_trace.csv'):
        for field in ('generated_files','source_tables'):
            for name in item.get(field,'').split('; '):
                if name:
                    p=Path(name)
                    assert p.is_file(), str(p)
                    selected.add(p)
                    if p.suffix.lower()=='.pdf':
                        for suffix in ('.svg','.png'):
                            if p.with_suffix(suffix).exists():selected.add(p.with_suffix(suffix))
    for folder in ('report','comparisons','tables'):
        for p in (root/folder).glob('*'):
            if p.is_file() and p.suffix.lower() in ('.pdf','.html','.csv','.json','.txt'):
                selected.add(p)
    for p in root.rglob('*manifest*.csv'):
        if not any(x in p.parts for x in ('cli_replay','QA')):
            selected.add(p)
    for p in (root/'inputs').glob('*.csv'):selected.add(p)
    if (root/'completion_manifest.json').exists():selected.add(root/'completion_manifest.json')
    for p in (root/'inputs/preprocessing_reference').glob('*'):
        if p.is_file() and p.suffix.lower() in ('.R','.csv'):selected.add(p)
    for p in (root/'QA').glob('*.csv'):selected.add(p)
    for p in (root/'report/qa').glob('*.csv'):selected.add(p)
    for p in (root/'report/qa/layout_final').rglob('*'):
        if p.is_file() and p.suffix.lower() in ('.csv','.json'):selected.add(p)
    for p in (root/'QA/fcg_visual').glob('*manifest*.csv'):selected.add(p)
    for p in (root/'runtime').glob('*.R'):selected.add(p)
    manuscript=root/'cf453_manuscript_panels'
    for pattern in ('panels/*.pdf','panels/*.png','source_data/*.csv','README.md','input_sha256.csv'):
        selected.update(p for p in manuscript.glob(pattern) if p.is_file())
    for d in ('cf453_models','cf453_direct_diet','cf453_rewiring','cf453_layouts','fcg05','fcg07','gdx_zonation'):
        for name in ('README.md','README.txt','METHODS.md','PANEL_LEGENDS.md','sessionInfo.txt','provenance_sha256.csv','completion_manifest.json'):
            p=root/d/name
            if p.exists():selected.add(p)
    payload={str(p.relative_to(root)):p.read_bytes() for p in selected}
    code=[]
    candidates=[Path('scripts/core/zonation_update_fcg.R'),Path('scripts/core/zonation_update_gdx.R'),
      Path('scripts/reports/run_zonation_update_fcg.R'),Path('scripts/reports/run_zonation_update_gdx.R'),
      Path('scripts/reports/stage_zonation_marker_update.R'),Path('scripts/reports/zonation_update_manuscript_audit.py'),
      Path('scripts/reports/build_zonation_update_report.py'),Path('scripts/reports/package_zonation_marker_update.py'),
      Path('scripts/reports/prepare_zonation_update_report_sources.py'),
      Path('scripts/reports/compare_cf453_zonation_update.R'),Path('tests/qa_zonation_marker_update.R'),
      Path('scripts/workflows/17_expanded_marker_zonation_update.R')]
    candidates += list(Path('scripts/core').glob('cf453*.R'))
    candidates += list(Path('scripts/workflows').glob('1[456]_cf453*.R'))
    candidates += [Path('scripts/core/wdsd_intake.R'),Path('scripts/reports/cf453_figure6a_layouts.R'),
                   Path('scripts/reports/cf453_manuscript_panels.R'),
                   Path('docs/workflow_plans/17_expanded_marker_zonation_update.md')]
    for p in candidates:
        if p.exists():
            name='code/'+str(p)
            payload[name]=p.read_bytes()
            code.append(dict(source=str(p.resolve()),bundle_path=name,sha256=sha(p)))
    b=io.StringIO();w=csv.DictWriter(b,fieldnames=['source','bundle_path','sha256']);w.writeheader();w.writerows(code)
    payload['code/source_manifest.csv']=b.getvalue().encode()
    payload['README.txt']=b'Expanded marker zonation update, 16 November 2026.\nOpen report/zonation_update_report.html for linked panel/source navigation or the PDF report for a printable audit.\nreport/manuscript_numeric_inventory.csv retains all 89 paragraph/caption claims and new_replacements.csv also tracks the 24 historical Figure5D edges.\nFigure4A and S4B-C are deferred to James. Independent measurements and historical analyses are preserved.\nThe bundle contains reports, panels and analysis sources; full versioned model objects and raw/staged RNA inputs remain in the analysis workspace and are identified by recorded hashes.\n'
    b=io.StringIO();w=csv.writer(b);w.writerow(['path','sha256','bytes'])
    for name,data in sorted(payload.items()):w.writerow([name,hashlib.sha256(data).hexdigest(),len(data)])
    payload['SHA256SUMS.csv']=b.getvalue().encode()
    dest=root/'Zonation_update_report_panels_and_source_data.zip'
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,data in sorted(payload.items()):z.writestr(name,data)
    with zipfile.ZipFile(dest) as z:
        assert z.testzip() is None
        assert all(hashlib.sha256(z.read(name)).digest()==hashlib.sha256(data).digest() for name,data in payload.items())
    summary=dict(bundle=str(dest),bytes=dest.stat().st_size,files=len(payload),historical_files_unchanged=len(preserved),
                 csv_sha256_and_zip_crc='passed',report_claims=report_manifest['claims'],panel_trace_rows=report_manifest['panels'])
    (root/'bundle_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', required=True)
    package(p.parse_args().output)
