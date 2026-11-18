"""Reproduce the 2026-09-08 PKN feasibility audit from frozen public downloads.

Run from the project root: python3 scripts/exploration/audit_cosmos_pkn.py
No package dependencies, downloads, omics scoring, or existing-output overwrites.
The protein backbone export is a candidate component, not a complete lipid PKN.
"""
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/cosmos_pkn_audit'


def read_csv(path, delimiter=','):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(name, rows, fields):
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def flag(value):
    value = str(value).lower()
    if value not in {'true', 'false', '1', '0'}:
        raise ValueError(f'Unexpected boolean: {value!r}')
    return value in {'true', '1'}


def main():
    manifest = json.loads((OUT / 'manifest.json').read_text())
    for source in manifest['sources']:
        data = gzip.decompress((ROOT / source['file']).read_bytes())
        assert hashlib.sha256(data).hexdigest() == source['sha256_uncompressed']
    with gzip.open(OUT / 'raw/om-mouse-pkn.json.gz', 'rt') as handle:
        cached = json.load(handle)
    network = cached['network']
    assert len(network) == cached['meta']['total_edges']
    resources = []
    for resource in sorted({e['resource'] for e in network}):
        counts = collections.Counter(e['mor'] for e in network if e['resource'] == resource)
        resources.append(dict(resource=resource, records=sum(counts.values()),
                              positive=counts[1], negative=counts[-1], unsigned=counts[0]))
    write_csv('cached_resource_audit.csv', resources,
              ['resource', 'records', 'positive', 'negative', 'unsigned'])

    # Keep explicit signed, directed evidence; quarantine conflicts across BOTH layers.
    signed = collections.defaultdict(set)
    provenance = collections.defaultdict(set)
    symbols = collections.defaultdict(set)
    fresh_stats = {}
    for layer in ['grn', 'ppi']:
        rows = read_csv(OUT / f'raw/om-fresh-{layer}.tsv.gz', '\t')
        counts = collections.Counter()
        counts['records'] = len(rows)
        for row in rows:
            for endpoint in ['source', 'target']:
                symbols[row[endpoint + '_genesymbol']].add(row[endpoint])
            if not flag(row['is_directed']):
                counts['undirected'] += 1
                continue
            signs = {s for key, s in [('is_stimulation', 1), ('is_inhibition', -1)]
                     if flag(row[key])}
            counts['dual_sign' if len(signs) == 2 else 'single_sign' if signs else 'unsigned'] += 1
            pair = row['source'], row['target']
            if not all(pair) or pair[0] == pair[1]:
                counts['missing_or_self'] += 1
                continue
            signed[pair].update(signs)
            provenance[pair].add((layer, row['sources'], row['references']))
        fresh_stats[layer] = dict(counts)
    backbone = [dict(source=a, interaction=next(iter(signs)), target=b)
                for (a, b), signs in sorted(signed.items()) if len(signs) == 1]
    conflicts = [dict(source=a, target=b) for (a, b), signs in sorted(signed.items())
                 if len(signs) == 2]
    assert backbone and all(e['interaction'] in {-1, 1} for e in backbone)
    assert len({(e['source'], e['target']) for e in backbone}) == len(backbone)
    write_csv('candidate_protein_backbone.csv', backbone, ['source', 'interaction', 'target'])
    write_csv('protein_sign_conflicts.csv', conflicts, ['source', 'target'])
    evidence = [dict(source=a, target=b, layer=layer, resources=res, references=refs)
                for (a, b), values in sorted(provenance.items())
                for layer, res, refs in sorted(values)]
    write_csv('protein_edge_evidence.csv', evidence,
              ['source', 'target', 'layer', 'resources', 'references'])

    lipids = read_csv(ROOT / 'results/lipidomic/gonade_type_ttop.csv')
    ids = [r['ID'] for r in lipids]
    assert len(ids) == len(set(ids))
    with (ROOT / 'support/linex2_adjacency_matrix_OT.csv').open() as handle:
        reader = csv.reader(handle)
        columns = next(reader)[1:]
        matrix = list(reader)
    assert all(len(row) == len(columns) + 1 for row in matrix)
    assert [row[0] for row in matrix] == columns
    pairs = {(row[0], columns[i]) for row in matrix
             for i, value in enumerate(row[1:]) if float(value) != 0}
    connected = {n for pair in pairs for n in pair}
    with (ROOT / 'results/mofa/networks/gene_lipid_correlations_factor1.csv').open() as handle:
        reader = csv.reader(handle)
        mofa_lipids = set(next(reader)[1:])
        mofa_genes = [row[0] for row in reader]
    inventory = [dict(feature_id=identifier, class_label=identifier.split('(')[0],
                      in_linex_nodes=identifier in columns,
                      connected_in_linex=identifier in connected,
                      in_mofa_factor1_matrix=identifier in mofa_lipids,
                      pkn_mapping_status='not_yet_mapped') for identifier in ids]
    write_csv('lipid_feature_inventory.csv', inventory, list(inventory[0]))
    enzymes = read_csv(OUT / 'raw/swisslipids-enzymes.tsv.gz', '\t')
    mouse_enzymes = [row for row in enzymes if row['Protein taxon'] == '10090']
    write_csv('swisslipids_mouse_enzymes.csv', mouse_enzymes, list(enzymes[0]))
    uniprot = read_csv(OUT / 'raw/focus-uniprot.tsv.gz', '\t')
    panel = ['Cyp2b9', 'Cyp2b13', 'Cyp2c38', 'Cyp2d9', 'Mtmr11', 'Sult3a1', 'Abcb1a']
    panel_rows = []
    for gene in panel:
        # Union of exact gene-name hits is an audit, not a final protein mapping.
        accessions = {r['Entry'] for r in uniprot if r['Gene Names (primary)'] == gene}
        counts = collections.Counter(e['resource'] for e in network
                                    if e['source'] in accessions or e['target'] in accessions)
        panel_rows.append(dict(gene=gene, uniprot_exact_gene_hits=';'.join(sorted(accessions)),
                               uniprot_in_fresh_protein_exports=';'.join(sorted(symbols[gene])),
                               cached_incident_records=json.dumps(dict(counts), sort_keys=True),
                               swisslipids_mouse_enzyme_rows=sum(
                                   r['Gene name'] == gene for r in mouse_enzymes)))
    write_csv('focus_gene_audit.csv', panel_rows, list(panel_rows[0]))
    summary = dict(
        cached_records=len(network), cached_resources=resources, fresh_protein_exports=fresh_stats,
        candidate_protein_edges=len(backbone), quarantined_protein_pairs=len(conflicts),
        current_lipid_DE_features=len(ids), linex_nodes=len(columns),
        linex_nonzero_matrix_entries=len(pairs), linex_symmetric=all((b, a) in pairs for a, b in pairs),
        DE_features_in_linex_nodes=len(set(ids) & set(columns)),
        DE_features_connected_in_linex=len(set(ids) & connected),
        mofa_factor1_matrix_genes=len(mofa_genes), mofa_factor1_matrix_lipids=len(mofa_lipids),
        mofa_lipids_connected_in_linex=len(mofa_lipids & connected),
        swisslipids_enzyme_rows=len(enzymes), swisslipids_mouse_enzyme_rows=len(mouse_enzymes),
        coverage_note='LINEX overlap is exact label overlap, not curated directed PKN coverage. '
                      'ChEBI/SwissLipids mapping and a complete lipid PKN remain unbuilt.')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
