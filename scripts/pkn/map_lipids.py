"""Build a resolution-aware SwissLipids crosswalk; no causal edges are created.

Run with Python >=3.9 and pygoslin==2.2.5, scipy==1.13.1, numpy==2.0.2:
  python scripts/pkn/map_lipids.py
The gzip snapshot is downloaded with --download only when missing. The script
keeps original feature IDs, reports nomenclature resolution (not experimental
identification confidence), and never treats class mappings as exact matches.
"""
import argparse
import collections
import csv
import datetime
import functools
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re
import sys
import urllib.request

from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/causal_modelling/cosmos_mouse_pkn'
RAW = ROOT / 'support/cosmos_mouse/raw'
URL = 'https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv'
PARSER = LipidParser()
LEVELS = [x for x in LipidLevel if x.value >= LipidLevel.CLASS.value]
HEADGROUP_EQUIVALENCE = {'GlcCer': 'HexCer', 'GalCer': 'HexCer',
                         'SulfoGalCer': 'SHexCer', 'SulfoGlcCer': 'SHexCer'}
NOMENCLATURE_ALIASES = {'TAG': 'TG', 'DAG': 'DG', 'PE-Cer': 'EPC',
                        'SulfoHexCer': 'SHexCer', 'GM3': 'Hex2NeuAcCer'}
FIELDS = {'Lipid ID': 'swisslipids_id', 'Level': 'swisslipids_level',
          'Name': 'name', 'Abbreviation*': 'abbreviation', 'CHEBI': 'chebi_ids',
          'Lipid class*': 'class_ids', 'Parent': 'parent_ids',
          'LIPID MAPS': 'lipidmaps_ids', 'HMDB': 'hmdb_ids'}
# These are deliberately broad classes verified against the downloaded names.
# Class annotations cannot imply structural specificity or signed regulation.
CLASS_IDS = {
    'PC': 'SLM:000055211', 'PE': 'SLM:000055213', 'PA': 'SLM:000055209',
    'PG': 'SLM:000055215', 'PI': 'SLM:000055216', 'PS': 'SLM:000055218',
    'LPC': 'SLM:000055200', 'LPE': 'SLM:000055202', 'LPA': 'SLM:000055199',
    'LPG': 'SLM:000055204', 'LPI': 'SLM:000055205', 'LPS': 'SLM:000055207',
    'Cer': 'SLM:000399814', 'HexCer': 'SLM:000399707', 'SM': 'SLM:000001000',
    'BMP': 'SLM:000477293', 'CL': 'SLM:000000102', 'MLCL': 'SLM:000000362',
    'GM3': 'SLM:000489935', 'FA': 'SLM:000000984',
    'TG': 'SLM:000000400', 'DG': 'SLM:000000401', 'EPC': 'SLM:000399706',
    'Hex2Cer': 'SLM:000399708', 'SHexCer': 'SLM:000399709',
    'Hex2NeuAcCer': 'SLM:000489935', 'SE 27:1': 'SLM:000000470',
}
ETHER_BROAD_IDS = {'PC': 'SLM:000001188', 'PE': 'SLM:000001189',
                   'PA': 'SLM:000001187', 'PG': 'SLM:000001190',
                   'PI': 'SLM:000001191', 'PS': 'SLM:000001192',
                   'LPC': 'SLM:000001188', 'LPE': 'SLM:000001189',
                   'LPA': 'SLM:000001187', 'LPG': 'SLM:000001190',
                   'LPI': 'SLM:000001191', 'LPS': 'SLM:000001192',
                   'TG': 'SLM:000117141', 'DG': 'SLM:000117140'}


def checksum(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def input_form(name):
    """Modern oxygen notation requires a space before chains in Goslin 2.2.5."""
    match = re.fullmatch(r'([^() ]+)\((.*)\)', name.strip())
    if match and ';' in match.group(2):
        return match.group(1) + ' ' + match.group(2)
    return name.strip()


def canonical_unordered(name):
    # Goslin preserves input ordering of '_' chains by default. Normalize that
    # ordering only; '/' remains untouched, including sphingoid-base/acyl roles.
    if ' ' in name:
        cls, body = name.split(' ', 1)
        if '_' in body and '/' not in body:
            return cls + ' ' + '_'.join(sorted(body.split('_')))
    return name


@functools.lru_cache(maxsize=None)
def parse(name):
    try:
        lipid = PARSER.parse(input_form(name))
        original_class = name.split('(')[0].split(' ')[0]
        if original_class == 'CL' and lipid.get_lipid_string(LipidLevel.CLASS) != 'CL':
            raise ValueError('Ambiguous two-moiety CL notation: automatic conversion to lysocardiolipin rejected')
        level = lipid.lipid.info.level
        strings = {}
        for item in LEVELS:
            if item.value <= level.value:
                try:
                    strings[item.name] = canonical_unordered(lipid.get_lipid_string(item))
                except Exception:
                    pass
        return {'canonical': canonical_unordered(lipid.get_lipid_string()), 'level': level.name,
                'rank': level.value, 'strings': strings,
                'class': lipid.get_lipid_string(LipidLevel.CLASS), 'error': ''}
    except Exception as error:
        return {'canonical': '', 'level': 'UNPARSED', 'rank': 0, 'strings': {},
                'class': name.split('(')[0].split(' ')[0], 'error': str(error)}


def class_alias(value):
    return HEADGROUP_EQUIVALENCE.get(value, value)


def replace_class(value):
    parts = value.split(' ', 1)
    return class_alias(parts[0]) + (' ' + parts[1] if len(parts) > 1 else '')


def rough_key(name):
    """A permissive prefilter only; acceptance always uses Goslin structures."""
    cls = re.split(r'[ (]', name, maxsplit=1)[0]
    chains = [(int(c), int(d)) for c, d in re.findall(r'(\d+):(\d+)', name)]
    if not chains:
        return None
    nonzero = sum(c != 0 for c, _ in chains)
    if cls in {'PC', 'PE', 'PA', 'PG', 'PI', 'PS'} and len(chains) > 1 and nonzero == 1:
        cls = 'L' + cls
    cls = NOMENCLATURE_ALIASES.get(cls, cls)
    carbon = sum(c for c, _ in chains)
    double_bonds = sum(d for _, d in chains) + name.count('P-')
    if cls == 'CE':
        cls, carbon, double_bonds = 'SE', carbon + 27, double_bonds + 1
    return class_alias(cls), carbon, double_bonds


def classify(measured, candidate):
    if not measured['rank'] or not candidate['rank']:
        return None
    if class_alias(measured['class']) != class_alias(candidate['class']):
        return None
    comparison_rank = min(measured['rank'], candidate['rank'])
    level = LipidLevel(comparison_rank).name
    left = measured['strings'].get(level)
    right = candidate['strings'].get(level)
    if not left or not right or replace_class(left) != replace_class(right):
        return None
    if measured['class'] != candidate['class']:
        return 'compatible_headgroup_isomer'
    if measured['rank'] == candidate['rank']:
        return 'exact_reported_resolution'
    if measured['rank'] > candidate['rank']:
        return 'compatible_lower_resolution'
    return 'compatible_more_specific'


def compact(row):
    result = {dest: row.get(src, '') for src, dest in FIELDS.items()}
    result['chebi_ids'] = '|'.join('CHEBI:' + x for x in re.findall(r'\d+', result['chebi_ids']))
    for field in ['parent_ids', 'class_ids']:
        result[field] = '|'.join(re.findall(r'SLM:\d+', result[field]))
    return result


def read_swiss(path):
    # The public download contains non-UTF8 characters; identifiers are ASCII.
    with gzip.open(path, 'rt', encoding='latin-1', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def export_chebi_nodes(candidate_rows):
    """Export identity and pH aliases as mappings only, never causal links."""
    path = RAW / 'rhea-chebi_pH7_3_mapping.tsv'
    ph73 = {}
    if path.exists():
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = 'CHEBI:' + row['CHEBI']
                value = ('CHEBI:' + row['CHEBI_PH7_3'], row['ORIGIN'])
                if key in ph73:
                    assert ph73[key] == value
                ph73[key] = value
    result = []
    for row in candidate_rows:
        for chebi in filter(None, row['chebi_ids'].split('|')):
            common = {key: row[key] for key in ['feature_id', 'swisslipids_id', 'match_type', 'in_factor1']}
            common['original_chebi_id'] = chebi
            result.append(dict(common, node_chebi_id=chebi, mapping_relation='SwissLipids identifier',
                               protonation_origin='not_applied'))
            if chebi in ph73 and ph73[chebi][0] != chebi:
                target, origin = ph73[chebi]
                result.append(dict(common, node_chebi_id=target,
                                   mapping_relation='Rhea major microspecies at pH 7.3', protonation_origin=origin))
    assert result and all(row['node_chebi_id'].startswith('CHEBI:') for row in result)
    write_csv(OUT / 'mapping_chebi_nodes.csv', result, list(result[0]))
    return {'rows': len(result), 'pH73_crosswalk_path': str(path.relative_to(ROOT)) if path.exists() else None,
            'pH73_crosswalk_sha256': checksum(path) if path.exists() else None,
            'pH73_changed_rows': sum(r['protonation_origin'] != 'not_applied' for r in result)}


def tests():
    assert parse('Cer(18:1;O2/18:0)')['canonical'] == 'Cer 18:1;O2/18:0'
    assert classify(parse('PE(20:4_18:0)'), parse('PE(18:0_20:4)')) == 'exact_reported_resolution'
    assert classify(parse('PC(18:0_20:4)'), parse('PC(18:0/20:4)')) == 'compatible_more_specific'
    assert classify(parse('PC(18:0/20:4)'), parse('PC(20:4/18:0)')) is None
    assert classify(parse('PC(18:0_20:4)'), parse('PC(38:4)')) == 'compatible_lower_resolution'
    assert classify(parse('HexCer(18:1;O2/18:0)'), parse('GlcCer(d18:1/18:0)')) == 'compatible_headgroup_isomer'
    assert classify(parse('Cer(18:1;O2/18:0)'), parse('Cer(t18:1/18:0)')) is None
    assert classify(parse('PC(O-18:0_20:4)'), parse('PC(18:0_20:4)')) is None
    assert classify(parse('FA(18:1;O)'), parse('FA(18:1)')) is None
    assert rough_key('PC(18:0/0:0)') == rough_key('LPC(18:0)')
    assert parse('CL(38:3_38:4)')['level'] == 'UNPARSED'
    assert rough_key('CE(18:2)') == rough_key(parse('CE(18:2)')['canonical'])
    assert rough_key('SulfoHexCer(d34:2)') == rough_key('SHexCer 34:2;O2')


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('--download', action='store_true')
    args.add_argument('--test-only', action='store_true')
    opts = args.parse_args()
    tests()
    if opts.test_only:
        print('Lipid crosswalk semantic tests passed.')
        return
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / 'swisslipids_lipids.tsv.gz'
    if not path.exists() and opts.download:
        urllib.request.urlretrieve(URL, path)
    if not path.exists():
        raise FileNotFoundError(f'{path}: use --download')
    with (ROOT / 'results/lipidomic/gonade_type_ttop.csv').open() as handle:
        measured_ids = [r['ID'] for r in csv.DictReader(handle)]
    with (ROOT / 'results/mofa/networks/gene_lipid_correlations_factor1.csv').open() as handle:
        factor_ids = set(next(csv.reader(handle))[1:])
    assert len(measured_ids) == len(set(measured_ids))
    universe = list(dict.fromkeys(measured_ids + sorted(factor_ids)))
    features = {fid: parse(fid) for fid in universe}
    lookup = collections.defaultdict(list)
    for fid, value in features.items():
        if value['rank']:
            lookup[rough_key(value['strings'].get('SPECIES', value['canonical']))].append(fid)
    print(f'Parsed {sum(bool(p["rank"]) for p in features.values())}/{len(features)} feature labels.', flush=True)
    catalog, matches, levels = {}, {}, collections.Counter()
    parse_failures = collections.Counter()
    total = 0
    # All ChEBI-bearing entries are retained as resource annotations. Species and
    # molecular records are parsed only in relevant composition buckets; detailed
    # isomers without ChEBI are intentionally not exhaustively enumerated.
    for row in read_swiss(path):
        total += 1
        levels[row['Level']] += 1
        sid = row['Lipid ID']
        if row['CHEBI'] or row['Level'] in {'Class', 'Category'}:
            catalog[sid] = compact(row)
        if row['Level'] in {'Class', 'Category', ''}:
            continue
        if row['Level'] == 'Isomeric subspecies' and not row['CHEBI']:
            continue
        for abbreviation in row['Abbreviation*'].split(' | '):
            possible = lookup.get(rough_key(abbreviation), [])
            if not possible:
                continue
            parsed = parse(abbreviation)
            if not parsed['rank']:
                parse_failures[row['Level']] += 1
                continue
            for fid in possible:
                kind = classify(features[fid], parsed)
                if not kind:
                    continue
                # Avoid generic species → thousands of isomers without useful
                # identifiers. More-specific candidates are retained only when
                # ChEBI-mapped or their shorthand is an exact reported label.
                if kind in {'compatible_more_specific', 'compatible_headgroup_isomer'} and not row['CHEBI']:
                    if parsed['rank'] > features[fid]['rank']:
                        continue
                key = fid, sid, kind
                matches[key] = {'feature_id': fid, 'swisslipids_id': sid, 'match_type': kind,
                                'matched_abbreviation': abbreviation, 'candidate_canonical': parsed['canonical'],
                                'candidate_nomenclature_level': parsed['level'], 'via_swisslipids_id': ''}
                catalog[sid] = compact(row)
        if total % 100000 == 0:
            print(f'Scanned {total:,} SwissLipids rows; {len(matches):,} candidates.', flush=True)
    # Parent paths are only annotations. They never become ±1 causal edges.
    for fid in universe:
        cls = features[fid]['class']
        sid = CLASS_IDS.get(cls)
        if re.search(r'[OP]-', features[fid]['canonical']):
            sid = ETHER_BROAD_IDS.get(cls)
        if sid:
            matches[(fid, sid, 'hierarchy_class')] = {
                'feature_id': fid, 'swisslipids_id': sid, 'match_type': 'hierarchy_class',
                'matched_abbreviation': '', 'candidate_canonical': '',
                'candidate_nomenclature_level': 'CLASS', 'via_swisslipids_id': 'explicit_broad_class_annotation'}
    pending = {p for row in catalog.values() for field in ['parent_ids', 'class_ids']
               for p in row[field].split('|') if p and p not in catalog}
    while pending:
        previous = len(catalog)
        for row in read_swiss(path):
            if row['Lipid ID'] in pending:
                catalog[row['Lipid ID']] = compact(row)
        pending = {p for row in catalog.values() for field in ['parent_ids', 'class_ids']
                   for p in row[field].split('|') if p and p not in catalog}
        if len(catalog) == previous:
            break
    # Ancestors are obtained from explicit SwissLipids Parent and lipid-class
    # links. Preserve which candidate introduced each relation.
    seeds = list(matches.values())
    for match in seeds:
        fid, origin = match['feature_id'], match['swisslipids_id']
        todo = [(origin, 0)]
        visited = set()
        while todo:
            sid, depth = todo.pop()
            if sid in visited or sid not in catalog:
                continue
            visited.add(sid)
            row = catalog[sid]
            if depth:
                kind = 'hierarchy_class' if row['swisslipids_level'] in {'Class', 'Category'} else 'hierarchy_ancestor'
                key = fid, sid, kind
                matches.setdefault(key, {'feature_id': fid, 'swisslipids_id': sid, 'match_type': kind,
                                        'matched_abbreviation': '', 'candidate_canonical': '',
                                        'candidate_nomenclature_level': '', 'via_swisslipids_id': origin})
            todo.extend((p, depth + 1) for field in ['parent_ids', 'class_ids'] for p in row[field].split('|') if p)
    candidate_rows = []
    for (_, sid, _), match in sorted(matches.items()):
        fid = match['feature_id']
        candidate_rows.append(dict(match, in_lipid_de=fid in measured_ids, in_factor1=fid in factor_ids,
                                   measured_canonical=features[fid]['canonical'],
                                   measured_nomenclature_level=features[fid]['level'], **{
                                       key: value for key, value in catalog[sid].items() if key != 'swisslipids_id'}))
    by_feature = collections.defaultdict(list)
    for row in candidate_rows:
        by_feature[row['feature_id']].append(row)
    summary_rows = []
    for fid in universe:
        rows = by_feature[fid]
        kinds = collections.Counter(row['match_type'] for row in rows)
        exact = [row for row in rows if row['match_type'] == 'exact_reported_resolution']
        compatible = [row for row in rows if row['match_type'].startswith('compatible_')]
        status = 'exact_reported_resolution' if exact else 'compatible_only' if compatible else 'class_or_ancestor_only' if rows else 'unmapped'
        summary_rows.append({'feature_id': fid, 'in_lipid_de': fid in measured_ids,
                             'in_factor1': fid in factor_ids, 'canonical_name': features[fid]['canonical'],
                             'class': features[fid]['class'], 'nomenclature_level': features[fid]['level'],
                             'parse_error': features[fid]['error'], 'best_mapping_status': status,
                             'exact_candidates': len(exact), 'compatible_candidates': len(compatible),
                             'hierarchy_candidates': sum(n for k, n in kinds.items() if k.startswith('hierarchy_')),
                             'exact_chebi_ids': '|'.join(sorted({c for r in exact for c in r['chebi_ids'].split('|') if c})),
                             'compatible_chebi_ids': '|'.join(sorted({c for r in compatible for c in r['chebi_ids'].split('|') if c})),
                             'ambiguity': 'hexose_isomer_unresolved' if features[fid]['class'] == 'HexCer' else
                                          'multiple_exact_resource_records' if len(exact) > 1 else ''})
    assert {row['feature_id'] for row in summary_rows} == set(universe)
    assert all(row['match_type'] != 'exact_reported_resolution' or
               row['measured_canonical'] == row['candidate_canonical'] for row in candidate_rows)
    write_csv(OUT / 'mapping_candidates.csv', candidate_rows, list(candidate_rows[0]))
    write_csv(OUT / 'mapping_features.csv', summary_rows, list(summary_rows[0]))
    write_csv(OUT / 'mapping_swisslipids_catalog.csv', sorted(catalog.values(), key=lambda r: r['swisslipids_id']), list(FIELDS.values()))
    chebi_node_metadata = export_chebi_nodes(candidate_rows)
    def coverage(rows):
        return {'n_features': len(rows), 'statuses': dict(collections.Counter(r['best_mapping_status'] for r in rows)),
                'exact_chebi_features': sum(bool(r['exact_chebi_ids']) for r in rows),
                'compatible_chebi_features': sum(bool(r['compatible_chebi_ids']) for r in rows),
                'parsed_features': sum(r['nomenclature_level'] != 'UNPARSED' for r in rows)}
    metadata = {'built_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'source': {'url': URL, 'path': str(path.relative_to(ROOT)), 'sha256_download': checksum(path),
                           'encoding': 'latin-1', 'rows': total, 'levels': dict(levels)},
                'dependencies': {name: importlib.metadata.version(name) for name in ['pygoslin', 'scipy', 'numpy']},
                'inputs': {str(p.relative_to(ROOT)): checksum(p) for p in [ROOT / 'results/lipidomic/gonade_type_ttop.csv', ROOT / 'results/mofa/networks/gene_lipid_correlations_factor1.csv']},
                'all_features': coverage(summary_rows), 'factor1': coverage([r for r in summary_rows if r['in_factor1']]),
                'candidate_rows': len(candidate_rows), 'catalog_rows': len(catalog),
                'chebi_node_mapping': chebi_node_metadata,
                'candidate_parse_failures': dict(parse_failures), 'unresolved_parent_ids': sorted(pending),
                'policy': ['Nomenclature resolution does not establish experimental structural identification confidence.',
                           'Exact means identical Goslin canonical shorthand at the reported resolution; multiple resource records remain alternatives.',
                           'Underscore and slash are distinct; matching to higher/lower resolution is explicitly labeled.',
                           'HexCer does not establish glucose versus galactose; such candidates are labeled headgroup isomers.',
                           'Two-moiety CL labels automatically reinterpreted as lysocardiolipins by Goslin are rejected; class annotation only is retained.',
                           'Class and ancestor links are annotations, never signed causal edges.',
                           'Detailed isomers without ChEBI are not exhaustively enumerated; compatibility candidate counts are not a complete chemical-isomer enumeration.',
                           'SwissLipids structures are not organism-filtered. Mouse specificity belongs to the catalytic/regulatory evidence, not lipid chemical identity.']}
    download_metadata = RAW / 'swisslipids_downloads.json'
    if download_metadata.exists():
        metadata['source']['metadata_file'] = str(download_metadata.relative_to(ROOT))
        metadata['source']['metadata_sha256'] = checksum(download_metadata)
        metadata['source']['release_metadata'] = next(
            row for row in json.loads(download_metadata.read_text()) if row['file'] == 'lipids.tsv')
    (OUT / 'mapping_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps({k: metadata[k] for k in ['all_features', 'factor1', 'candidate_rows', 'catalog_rows']}, indent=2))


if __name__ == '__main__':
    main()
