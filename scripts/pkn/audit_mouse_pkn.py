"""Chemical-normalization sensitivity graph and coverage (never input scores).

The optional pH7.3 graph uses ChEBI alt_ids plus Rhea target forms only when
the forms belong to an explicit ChEBI conjugate acid/base component, or the
Rhea mapping itself is curated. No ontology class relation becomes an edge.
"""
import collections
import csv
import gzip
import json
from pathlib import Path
from assemble_mouse_pkn import OUT, table, jsonlines


def write_table(name, rows, fields=None):
    fields = fields or list(rows[0])
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class ChemicalMap:
    def __init__(self):
        self.alias = {r['alias']: r['primary_chebi']
                      for r in table(OUT / 'regulation_chebi_aliases.csv')}
        self.parent = {}
        for row in table(OUT / 'regulation_chebi_acid_base_relations.csv'):
            a, b = self.find(row['acid']), self.find(row['base'])
            self.parent[a] = b
        targets = collections.defaultdict(set)
        reasons = collections.defaultdict(set)
        self.ph_origins = collections.defaultdict(set)
        for row in table(OUT / 'rhea_chebi_ph73.tsv', '\t'):
            source, target = self.primary(row['source']), self.primary(row['target'])
            if source == target:
                continue
            explicit = self.find(source) == self.find(target)
            curated = row['origin'].lower() in {'curation', 'curated'}
            if explicit or curated:
                targets[source].add(target)
                reasons[source].add('ChEBI_conjugate_acid_base' if explicit else 'Rhea_curated_ph73')
                self.ph_origins[source].add(row['origin'])
        self.targets = {s: next(iter(ts)) for s, ts in targets.items() if len(ts) == 1}
        self.reasons = {s: '|'.join(sorted(rs)) for s, rs in reasons.items()}

    def primary(self, identifier):
        seen = set()
        while identifier in self.alias and identifier not in seen:
            seen.add(identifier)
            identifier = self.alias[identifier]
        return identifier

    def find(self, identifier):
        identifier = self.primary(identifier)
        if identifier not in self.parent:
            self.parent[identifier] = identifier
        if self.parent[identifier] != identifier:
            self.parent[identifier] = self.find(self.parent[identifier])
        return self.parent[identifier]

    def map(self, identifier):
        primary = self.primary(identifier)
        return self.targets.get(primary, primary)


def main():
    nodes = {r['node_id']: r for r in jsonlines(OUT / 'nodes.jsonl')}
    chem = ChemicalMap()
    aliases = []
    node_map = {}
    for old, row in nodes.items():
        new = old
        if row['node_type'] == 'metabolite':
            mapped = chem.map(row['chemical_id'])
            new = 'Metab__' + mapped + '_' + row['compartment']
            if new != old:
                aliases.append(dict(original_node=old, normalized_node=new,
                    original_chemical=row['chemical_id'], normalized_chemical=mapped,
                    compartment=row['compartment'],
                    basis=chem.reasons.get(chem.primary(row['chemical_id']), 'ChEBI_alt_id'),
                    ph_target_selection_origin='|'.join(sorted(chem.ph_origins.get(chem.primary(row['chemical_id']), [])))))
        node_map[old] = new
    write_table('chemical_node_aliases.csv', aliases,
                ['original_node', 'normalized_node', 'original_chemical', 'normalized_chemical', 'compartment', 'basis', 'ph_target_selection_origin'])
    normalized_nodes = {}
    for old, row in nodes.items():
        new = node_map[old]
        if new not in normalized_nodes:
            item = dict(row)
            item['node_id'] = new
            if row['node_type'] == 'metabolite':
                item['chemical_id'] = chem.map(row['chemical_id'])
            normalized_nodes[new] = item
            normalized_nodes[new]['original_nodes'] = []
        normalized_nodes[new]['original_nodes'].append(old)
    with (OUT / 'nodes_ph73.jsonl').open('w') as f:
        for key, value in sorted(normalized_nodes.items()):
            f.write(json.dumps(value, sort_keys=True) + '\n')

    normalized = collections.defaultdict(set)
    # Normalize ALL eligible evidence before sign resolution, including conflicts
    # excluded from the original graph. An alias must not resurrect a conflict.
    with gzip.open(OUT / 'edge_evidence.jsonl.gz', 'rt') as f:
        for line in f:
            edge = json.loads(line)
            if edge['first_variant'] not in {'core', 'lipid_curated'}:
                continue
            a, b = node_map[edge['source']], node_map[edge['target']]
            if a != b:
                normalized[a, b].add(int(edge['interaction']))
    ph_edges = [dict(source=a, interaction=next(iter(signs)), target=b)
                for (a, b), signs in sorted(normalized.items()) if len(signs) == 1]
    write_table('pkn_lipid_curated_ph73.csv', ph_edges, ['source', 'interaction', 'target'])
    write_table('conflicts_lipid_curated_ph73.csv',
                [dict(source=a, target=b) for (a,b), signs in sorted(normalized.items()) if len(signs)>1],
                ['source','target'])

    candidates = table(OUT / 'mapping_candidates.csv')
    mapping = collections.defaultdict(lambda: collections.defaultdict(set))
    for r in candidates:
        category = ('exact' if r['match_type'] == 'exact_reported_resolution' else
                    'compatible' if r['match_type'].startswith('compatible_') else 'hierarchy')
        mapping[r['feature_id']][category].update(x for x in r['chebi_ids'].split('|') if x)
    features = table(OUT / 'mapping_features.csv')
    variants = ['core', 'lipid_curated', 'lipid_curated_ph73', 'lipid_rules', 'exploratory']
    coverage, gene_rows, variant_stats, topology = [], [], {}, []
    panel = ['Cyp2b9', 'Cyp2b13', 'Cyp2c38', 'Cyp2d9', 'Mtmr11', 'Sult3a1', 'Abcb1a']
    all_lipids = {r['id'] for r in table(OUT / 'rhea_all_lipid_chebi.tsv', '\t')}
    all_lipids.update(r['source'] for r in table(OUT / 'regulation_chebi_source_crosswalk.csv')
                      if r['lipid_is_a'] == 'True')
    all_lipids.update(chem.map(x) for x in list(all_lipids))
    examples = []
    for variant in variants:
        nt = normalized_nodes if variant.endswith('ph73') else nodes
        edges = table(OUT / f'pkn_{variant}.csv')
        present = {r[k] for r in edges for k in ['source', 'target']}
        chemicals = {nt[n]['chemical_id'] for n in present if nt[n]['node_type'] == 'metabolite'}
        outgoing, incoming = collections.defaultdict(set), collections.defaultdict(set)
        for e in edges:
            outgoing[e['source']].add(e['target'])
            incoming[e['target']].add(e['source'])
        if variant == 'lipid_rules':
            factor_features = {r['feature_id'] for r in features if r['in_factor1'] == 'True'}
            for seed in ['Esr1', 'Ar', 'Stat5b', 'Nfic', 'Dbp']:
                seen = {'Gene__' + seed}
                queue = collections.deque([('Gene__' + seed, 0)])
                while queue:
                    node, depth = queue.popleft()
                    if depth >= 6:
                        continue
                    for child in outgoing[node]:
                        if child not in seen:
                            seen.add(child)
                            queue.append((child, depth + 1))
                hits = sorted(f for f in factor_features if
                              'Metab__MeasuredLipid:' + f + '_unknown' in seen)
                topology.append(dict(query_gene=seed, maximum_edges=6,
                    rule_factor1_lipids_reachable=len(hits), feature_ids='|'.join(hits),
                    interpretation='Topology query only; no TF activity, causal validation or MOON input selection.'))
        for row in features:
            fid = row['feature_id']
            hits = {}
            for tier in ['exact', 'compatible', 'hierarchy']:
                ids = mapping[fid][tier]
                if variant.endswith('ph73'):
                    ids = {chem.map(x) for x in ids}
                hits[tier] = sorted(ids & chemicals)
            coverage.append(dict(feature_id=fid, in_factor1=row['in_factor1'], variant=variant,
                exact_chemical_ids='|'.join(hits['exact']), compatible_chemical_ids='|'.join(hits['compatible']),
                hierarchy_chemical_ids='|'.join(hits['hierarchy']),
                species_rule_connected=('MeasuredLipid:' + fid) in chemicals,
                interpretation='Mapping candidates only; compatible/hierarchy matches do not establish species specificity.'))
        for gene in panel:
            gid = 'Gene__' + gene
            instances = {n for n in outgoing[gid] if nt[n]['node_type'] == 'enzyme_reaction_instance'}
            lipid_inputs, lipid_outputs = set(), set()
            def is_lipid(n):
                return nt[n]['node_type'] == 'metabolite' and (nt[n]['chemical_id'] in all_lipids or
                        nt[n]['chemical_id'].startswith('MeasuredLipid:'))
            for n in sorted(instances):
                inputs = {x for x in incoming[n] if is_lipid(x)}
                outputs = {x for x in outgoing[n] if is_lipid(x)}
                lipid_inputs.update(inputs)
                lipid_outputs.update(outputs)
                if variant == 'lipid_curated' and inputs and outputs and sum(x['gene'] == gene for x in examples) < 5:
                    examples.append(dict(gene=gene, reaction_id=nt[n]['reaction_id'],
                        resource=nt[n]['resource'], enzyme_node=n,
                        lipid_substrates='|'.join(sorted(inputs)), lipid_products='|'.join(sorted(outputs)),
                        note='Curated reaction traversal; lipid nodes may be chemical classes, not measured species.'))
            gene_rows.append(dict(gene=gene, variant=variant, present=gid in present,
                enzyme_instances=len(instances), lipid_substrate_nodes=len(lipid_inputs), lipid_product_nodes=len(lipid_outputs),
                incoming_lipid_regulators=sum(is_lipid(n) for n in incoming[gid])))
        subset = [r for r in coverage if r['variant'] == variant]
        counts = {}
        for label, rows in [('all_features', subset), ('factor1', [r for r in subset if r['in_factor1'] == 'True'])]:
            counts[label] = dict(total=len(rows),
                exact=sum(bool(r['exact_chemical_ids']) for r in rows),
                compatible=sum(bool(r['compatible_chemical_ids']) for r in rows),
                hierarchy=sum(bool(r['hierarchy_chemical_ids']) for r in rows))
            counts[label]['species_rule_connected'] = sum(r['species_rule_connected'] for r in rows)
        variant_stats[variant] = dict(edges=len(edges), nodes=len(present), chemical_ids=len(chemicals), coverage=counts)
    write_table('feature_coverage.csv', coverage)
    write_table('focus_gene_coverage.csv', gene_rows)
    write_table('topology_queries.csv', topology)
    write_table('focus_reaction_examples.csv', examples,
                ['gene','reaction_id','resource','enzyme_node','lipid_substrates','lipid_products','note'])
    (OUT / 'coverage_summary.json').write_text(json.dumps(variant_stats, indent=2) + '\n')
    print(json.dumps(variant_stats, indent=2))


if __name__ == '__main__':
    main()
