"""Lossless Mouse-GEM SBML extraction; no fitted data or inferred signs.

Run from any directory: python3 scripts/pkn/build_mouse_gem.py
The master reactions retain the SBML GPR tree, full stoichiometry, bounds and
all chemical cross-references. Ambiguous ChEBI mappings remain model entities.
"""
import csv
import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / 'support/cosmos_mouse/raw'
OUT = ROOT / 'results/cosmos_mouse_pkn'
COMMIT = 'fc769ac1da421c4e4e9e1cf2a9aec3fe4f02ef6a'
CORE = '{http://www.sbml.org/sbml/level3/version1/core}'
FBC = '{http://www.sbml.org/sbml/level3/version1/fbc/version2}'
RDF = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'


def xrefs(element):
    return sorted({e.attrib[RDF + 'resource'] for e in element.iter()
                   if RDF + 'resource' in e.attrib})


def gpr_tree(element):
    if element is None:
        return None
    kind = element.tag.rsplit('}', 1)[-1]
    if kind == 'geneProductRef':
        return {'gene': element.attrib[FBC + 'geneProduct']}
    children = [gpr_tree(e) for e in element]
    if kind == 'geneProductAssociation':
        assert len(children) == 1
        return children[0]
    if kind not in {'and', 'or'}:
        raise ValueError(f'Unknown GPR operator {kind}')
    return {kind: children}


def parse_model(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rb') as handle:
        root = ET.parse(handle).getroot()
    model = root.find(CORE + 'model')
    parameters = {e.attrib['id']: float(e.attrib['value'])
                  for e in model.findall(f'{CORE}listOfParameters/{CORE}parameter')}
    metabolites = {}
    for e in model.findall(f'{CORE}listOfSpecies/{CORE}species'):
        refs = xrefs(e)
        chebi = sorted({m.group(0) for ref in refs for m in re.finditer(r'CHEBI:\d+', ref)})
        model_id = e.attrib['id']
        # No arbitrary first-ID selection for ambiguous chemical mappings.
        identifier = chebi[0] if len(chebi) == 1 else 'MouseGEM:' + model_id
        metabolites[model_id] = dict(
            id=identifier, model_id=model_id, name=e.attrib.get('name', model_id),
            compartment=e.attrib['compartment'], chebi_ids=chebi, xrefs=refs,
            formula=e.attrib.get(FBC + 'chemicalFormula'),
            charge=e.attrib.get(FBC + 'charge'), mapping_status=(
                'unique_chebi' if len(chebi) == 1 else 'ambiguous_chebi' if chebi else 'model_only'))
    reactions = []
    for e in model.findall(f'{CORE}listOfReactions/{CORE}reaction'):
        sides = {}
        for side, tag in [('substrates', 'listOfReactants'), ('products', 'listOfProducts')]:
            sides[side] = []
            for participant in e.findall(f'{CORE}{tag}/{CORE}speciesReference'):
                entry = metabolites[participant.attrib['species']]
                coeff = float(participant.attrib.get('stoichiometry', '1'))
                assert coeff > 0
                sides[side].append(dict(id=entry['id'], model_id=entry['model_id'],
                                        compartment=entry['compartment'], stoichiometry=coeff))
        lower = parameters[e.attrib[FBC + 'lowerFluxBound']]
        upper = parameters[e.attrib[FBC + 'upperFluxBound']]
        assert lower <= upper
        reactions.append(dict(
            reaction_id=e.attrib['id'], name=e.attrib.get('name'), resource='Mouse-GEM',
            organism=10090, source_version=COMMIT, **sides,
            gpr=gpr_tree(e.find(FBC + 'geneProductAssociation')),
            lower_bound=lower, upper_bound=upper, forward_allowed=upper > 0,
            reverse_allowed=lower < 0, reversible=lower < 0 < upper,
            sbml_reversible=e.attrib.get('reversible'), xrefs=xrefs(e),
            notes=' '.join(e.find(CORE + 'notes').itertext()).strip()
            if e.find(CORE + 'notes') is not None else ''))
    assert metabolites and reactions
    return list(metabolites.values()), reactions


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    path = RAW / 'Mouse-GEM.xml'
    metabolites, reactions = parse_model(path)
    for filename, rows in [('mouse_gem_metabolites.jsonl', metabolites),
                           ('mouse_gem_reactions.jsonl', reactions)]:
        with (OUT / filename).open('w') as handle:
            for row in rows:
                handle.write(json.dumps(row, sort_keys=True) + '\n')
    summary = dict(model_commit=COMMIT, metabolites=len(metabolites), reactions=len(reactions),
                   reactions_with_gpr=sum(r['gpr'] is not None for r in reactions),
                   blocked_reactions=sum(not r['forward_allowed'] and not r['reverse_allowed']
                                         for r in reactions),
                   reversible_reactions=sum(r['reversible'] for r in reactions),
                   all_participants_preserved=True)
    sources = []
    for name, upstream in [('Mouse-GEM.xml', 'Mouse-GEM.xml'),
                           ('mouse_gem_metabolites.tsv', 'metabolites.tsv'),
                           ('mouse_gem_reactions.tsv', 'reactions.tsv')]:
        file = RAW / name
        sources.append(dict(file=str(file.relative_to(ROOT)),
                            sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
                            url=f'https://raw.githubusercontent.com/SysBioChalmers/Mouse-GEM/{COMMIT}/model/{upstream}'))
    (OUT / 'mouse_gem_manifest.json').write_text(json.dumps(dict(summary=summary, sources=sources), indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
