"""Assemble auditable COSMOS projections from independently built mouse layers.

Node IDs are explicit; use the node table, not legacy COSMOS name regexes.
Metabolism is a positive reaction traversal (substrate -> enzyme instance ->
product), not an inferred concentration/flux relationship. Complex requirements
and uncertain directions are allowed only in the exploratory sensitivity PKN.
"""
import collections
import csv
import gzip
import hashlib
import itertools
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'results/causal_modelling/cosmos_mouse_pkn'
RAW = ROOT / 'support/cosmos_mouse/raw'
AUDIT = ROOT / 'results/causal_modelling/cosmos_pkn_audit'
VARIANTS = ['core', 'lipid_curated', 'lipid_rules', 'exploratory']


def table(path, delimiter=','):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def jsonlines(path):
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def flag(x):
    x = str(x).lower()
    if x not in {'true', 'false', '1', '0'}:
        raise ValueError(f'Unknown boolean {x}')
    return x in {'true', '1'}


def chemical_id(identifier):
    value = str(identifier).strip()
    if re.fullmatch(r'(CHEBI:)?\d+', value):
        return 'CHEBI:' + value.split(':')[-1]
    return value


def dnf(tree, limit=256):
    """Preserve AND complexes and OR alternatives; fail on combinatorial blowup."""
    if not tree:
        return []
    if 'gene' in tree:
        return [frozenset([tree['gene']])]
    if 'or' in tree:
        groups = set().union(*(set(dnf(t, limit)) for t in tree['or']))
    elif 'and' in tree:
        groups = {frozenset()}
        for t in tree['and']:
            branch = dnf(t, limit)
            if len(groups) * len(branch) > limit:
                raise ValueError('GPR expansion exceeds limit')
            groups = {a | b for a in groups for b in branch}
    else:
        raise ValueError('Unrecognized GPR tree')
    if len(groups) > limit:
        raise ValueError('GPR expansion exceeds limit')
    return sorted(groups, key=lambda s: sorted(s))


class Builder:
    def __init__(self):
        self.nodes = {}
        self.evidence = []
        self.exclusions = collections.Counter()
        self.gene_map = {}
        self.mouse_accessions = set()
        for row in table(RAW / 'uniprot_mouse.tsv.gz', '\t'):
            self.mouse_accessions.add(row['Entry'])
            symbol = row['Gene Names (primary)']
            if re.fullmatch(r'[A-Za-z0-9_.-]+', symbol):
                self.gene_map[row['Entry']] = symbol

    def node(self, identifier, kind, label='', **attributes):
        if identifier not in self.nodes:
            self.nodes[identifier] = dict(node_id=identifier, node_type=kind,
                                          label=label or identifier, **attributes)
        else:
            assert self.nodes[identifier]['node_type'] == kind
        return identifier

    def gene(self, symbol):
        return self.node('Gene__' + symbol, 'gene_product', symbol,
                         gene_symbol=symbol, organism=10090)

    def protein(self, accession):
        base = re.sub(r'-\d+$', '', accession)
        if accession in self.gene_map or base in self.gene_map:
            return self.gene(self.gene_map.get(accession, self.gene_map.get(base)))
        # Unmapped/multi-locus accessions keep their identity, never split arbitrarily.
        return self.node('UniProt__' + accession, 'protein_unmapped', accession,
                         uniprot=accession, organism=10090 if base in self.mouse_accessions else None,
                         mouse_accession_verified=base in self.mouse_accessions)

    def metabolite(self, identifier, compartment='unknown', label=''):
        identifier = chemical_id(identifier)
        compartment = compartment or 'unknown'
        return self.node('Metab__' + identifier + '_' + compartment, 'metabolite',
                         label or identifier, chemical_id=identifier, compartment=compartment)

    def emit(self, source, target, sign, variant, resource, mechanism, **metadata):
        assert sign in {-1, 1} and variant in VARIANTS
        if any(self.nodes[n]['node_type'] == 'protein_unmapped' and
               not self.nodes[n].get('mouse_accession_verified') for n in [source, target]):
            variant = 'exploratory'
        if source == target:
            self.exclusions['self_edges'] += 1
            return
        self.evidence.append(dict(source=source, target=target, interaction=sign,
                                  first_variant=variant, resource=resource,
                                  mechanism=mechanism, **metadata))

    def proteins(self):
        for layer in ['grn', 'ppi']:
            for row in table(AUDIT / f'raw/om-fresh-{layer}.tsv.gz', '\t'):
                if not flag(row['is_directed']):
                    self.exclusions['protein_undirected'] += 1
                    continue
                signs = [sign for key, sign in [('is_stimulation', 1), ('is_inhibition', -1)]
                         if flag(row[key])]
                if not signs:
                    self.exclusions['protein_unsigned'] += 1
                for sign in signs:
                    self.emit(self.protein(row['source']), self.protein(row['target']), sign,
                              'core', 'OmniPath:' + layer, 'transcription' if layer == 'grn' else 'signaling',
                              references=row['references'], resources=row['sources'],
                              source_accession=row['source'], target_accession=row['target'],
                              species_evidence='OmniPath mouse export; orthology may underlie source annotations')

    def reaction(self, record, groups, directions, first_variant):
        if not groups or not record['substrates'] or not record['products']:
            self.exclusions['reaction_no_catalyst_or_boundary'] += 1
            return
        for reverse, uncertain in directions:
            substrates, products = record['substrates'], record['products']
            if reverse:
                substrates, products = products, substrates
            for genes in groups:
                unverified = any(self.nodes[g]['node_type'] == 'protein_unmapped' and
                                 not self.nodes[g].get('mouse_accession_verified') for g in genes)
                variant = 'exploratory' if len(genes) > 1 or uncertain or unverified else first_variant
                signature = json.dumps([record['resource'], record['reaction_id'], reverse, sorted(genes)])
                digest = hashlib.sha256(signature.encode()).hexdigest()[:16]
                node = self.node('Enzyme__' + digest, 'enzyme_reaction_instance',
                                 record['reaction_id'] + (' reverse' if reverse else ' forward'),
                                 reaction_id=record['reaction_id'], resource=record['resource'],
                                 reverse=reverse, catalysts=sorted(genes),
                                 complex_required=len(genes) > 1, uncertain_direction=uncertain)
                metadata = dict(reaction_id=record['reaction_id'], reverse=reverse,
                                uncertain_direction=uncertain, complex_required=len(genes) > 1,
                                source_version=record.get('source_version', ''),
                                annotated_direction=record.get('direction', 'model_flux_bounds'),
                                direction_basis=record.get('direction_basis', 'Mouse-GEM allowed flux bounds'),
                                master_record='mouse_gem_reactions.jsonl' if record['resource'] == 'Mouse-GEM'
                                else 'lipid_rules_reactions.jsonl' if record['resource'] == 'LINEX2_rule'
                                else 'rhea_reactions.jsonl')
                for gene in genes:
                    self.emit(gene, node, 1, variant, record['resource'], 'catalyst_enables_reaction', **metadata)
                for participant in substrates:
                    metab = self.metabolite(participant['id'], participant.get('compartment', 'unknown'), participant.get('name', ''))
                    self.emit(metab, node, 1, variant, record['resource'], 'substrate_participation',
                              stoichiometry=participant['stoichiometry'], **metadata)
                for participant in products:
                    metab = self.metabolite(participant['id'], participant.get('compartment', 'unknown'), participant.get('name', ''))
                    self.emit(node, metab, 1, variant, record['resource'], 'product_formation',
                              stoichiometry=participant['stoichiometry'], **metadata)

    def gem(self):
        for m in jsonlines(OUT / 'mouse_gem_metabolites.jsonl'):
            self.metabolite(m['id'], m['compartment'], m['name'])
        for row in jsonlines(OUT / 'mouse_gem_reactions.jsonl'):
            try:
                groups = [[self.gene(g) for g in group] for group in dnf(row['gpr'])]
            except ValueError:
                self.exclusions['gpr_over_limit'] += 1
                continue
            directions = [(False, False)] if row['forward_allowed'] else []
            if row['reverse_allowed']:
                directions.append((True, False))
            self.reaction(row, groups, directions, 'core')

    def rhea(self):
        rows = jsonlines(OUT / 'rhea_reactions.jsonl')
        directional = {(r['master_reaction_id'], p) for r in rows
                       if r['reversible'] is not None for p in r['catalysts']}
        for row in rows:
            if row.get('is_transport'):
                self.exclusions['rhea_transport_unresolved_compartments'] += 1
                continue
            catalysts = [p for p in row['catalysts'] if row['reversible'] is not None or
                         (row['master_reaction_id'], p) not in directional]
            if not catalysts:
                self.exclusions['rhea_redundant_UN_annotation'] += 1
                continue
            groups = [[self.protein(p)] for p in catalysts]
            if row['reversible'] is None:
                directions = [(False, True), (True, True)]
            else:
                directions = [(False, False)] + ([(True, False)] if row['reversible'] else [])
            generic = any(not p['id'].startswith('CHEBI:') or
                          not isinstance(p['stoichiometry'], (float, int))
                          for p in row['substrates'] + row['products'])
            self.reaction(row, groups, directions, 'exploratory' if generic else 'lipid_curated')

    def regulatory(self):
        for row in table(OUT / 'regulation_edges.csv'):
            # The cache-normalized STITCH export lost raw actor orientation.
            # A direction-safe STITCH layer is assembled separately from raw
            # source records by build_stitch_directional.py.
            if row['resource'] == 'STITCH':
                self.exclusions['cached_STITCH_raw_actor_orientation_unavailable'] += 1
                continue
            sign = int(row['interaction'])
            # Fresh native evidence in curated; cached fallbacks only in sensitivity network.
            tier = row['evidence_tier']
            variant = 'lipid_curated' if tier == 'curated_mouse_current' else 'exploratory'
            self.emit(self.metabolite(row['source'], row.get('compartment', 'unknown'), row.get('source_name', '')),
                      self.protein(row['target']), sign, variant, row['resource'], row['mechanism'],
                      evidence_tier=tier, references=row.get('references', ''),
                      source_version=row.get('source_version', ''),
                      original_edge_id=row['edge_id'],
                      source_accession=row['source'], target_accession=row['target'])

    def lipid_rules(self):
        for row in jsonlines(OUT / 'lipid_rules_reactions.jsonl'):
            assert row['reversible'] is False and row['direction'] in {'LR', 'RL'}
            groups = [[self.protein(p)] for p in row['catalysts']]
            self.reaction(row, groups, [(False, False)], 'lipid_rules')

    def export(self):
        summary = {'variants': {}, 'excluded_records': dict(self.exclusions),
                   'semantics': 'Positive metabolic traversal; no concentration or flux inference',
                   'protein_namespace': 'Gene__mouse_primary_symbol; unmapped accessions preserved',
                   'compartments': 'unknown is distinct; no synthetic compartment bridges'}
        with gzip.open(OUT / 'edge_evidence.jsonl.gz', 'wt') as f:
            for edge in self.evidence:
                f.write(json.dumps(edge, sort_keys=True) + '\n')
        for v, variant in enumerate(VARIANTS):
            pairs = collections.defaultdict(set)
            for e in self.evidence:
                if VARIANTS.index(e['first_variant']) <= v:
                    pairs[e['source'], e['target']].add(e['interaction'])
            edges = [(a, next(iter(signs)), b) for (a, b), signs in sorted(pairs.items()) if len(signs) == 1]
            conflicts = [(a, b) for (a, b), signs in sorted(pairs.items()) if len(signs) > 1]
            nodes = {n for a, sign, b in edges for n in (a, b)}
            connected_pairs = {(a, b) for a, sign, b in edges}
            for node in nodes:
                if self.nodes[node]['node_type'] == 'enzyme_reaction_instance':
                    assert all((g, node) in connected_pairs for g in self.nodes[node]['catalysts']), (
                        variant, 'Reaction lost a required catalyst connector', node)
            with (OUT / f'pkn_{variant}.csv').open('w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['source', 'interaction', 'target'])
                writer.writerows(edges)
            with (OUT / f'conflicts_{variant}.csv').open('w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['source', 'target'])
                writer.writerows(conflicts)
            summary['variants'][variant] = dict(edges=len(edges), nodes=len(nodes),
                conflicting_pairs_excluded=len(conflicts),
                node_types=dict(collections.Counter(self.nodes[n]['node_type'] for n in nodes)))
        with (OUT / 'nodes.jsonl').open('w') as f:
            for key, node in sorted(self.nodes.items()):
                f.write(json.dumps(node, sort_keys=True) + '\n')
        fields = ['node_id','label','node_type','gene_symbol','chemical_id','compartment',
                  'reaction_id','resource','catalysts','complex_required','uncertain_direction',
                  'organism','mouse_accession_verified']
        with (OUT / 'nodes.csv').open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for node in self.nodes.values():
                writer.writerow({k: json.dumps(node[k]) if isinstance(node.get(k), (list,dict))
                                 else node.get(k, '') for k in fields})
        (OUT / 'build_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(json.dumps(summary, indent=2))


def main():
    builder = Builder()
    builder.proteins()
    builder.gem()
    builder.rhea()
    builder.regulatory()
    builder.lipid_rules()
    builder.export()


if __name__ == '__main__':
    main()
