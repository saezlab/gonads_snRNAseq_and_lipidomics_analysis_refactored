"""Semantic regression tests for reaction extraction and COSMOS projection."""
import json
import tempfile
import unittest
from pathlib import Path

from build_mouse_gem import parse_model
from assemble_mouse_pkn import Builder, dnf, chemical_id, flag


def empty_builder():
    b = Builder.__new__(Builder)
    import collections
    b.nodes, b.evidence, b.gene_map = {}, [], {}
    b.mouse_accessions = set()
    b.exclusions = collections.Counter()
    return b


class ReactionTests(unittest.TestCase):
    def test_unbalanced_participant_counts_and_reverse_only_bounds(self):
        # Unequal numbers of reactants/products must never be zipped/truncated.
        xml = '''<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
          xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"><model>
          <listOfParameters><parameter id="lo" value="-1000"/><parameter id="hi" value="0"/></listOfParameters>
          <listOfSpecies><species id="A" compartment="c"/><species id="B" compartment="c"/>
          <species id="C" compartment="m"/></listOfSpecies>
          <listOfReactions><reaction id="R" reversible="true" fbc:lowerFluxBound="lo" fbc:upperFluxBound="hi">
          <listOfReactants><speciesReference species="A" stoichiometry="2"/></listOfReactants>
          <listOfProducts><speciesReference species="B" stoichiometry="1"/>
          <speciesReference species="C" stoichiometry="3"/></listOfProducts>
          <fbc:geneProductAssociation><fbc:and><fbc:geneProductRef fbc:geneProduct="G1"/>
          <fbc:geneProductRef fbc:geneProduct="G2"/></fbc:and></fbc:geneProductAssociation>
          </reaction></listOfReactions></model></sbml>'''
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'model.xml'
            path.write_text(xml)
            molecules, reactions = parse_model(path)
        r = reactions[0]
        self.assertEqual(len(r['substrates']), 1)
        self.assertEqual(len(r['products']), 2)
        self.assertEqual(r['products'][1]['stoichiometry'], 3)
        self.assertEqual(r['products'][1]['compartment'], 'm')
        self.assertFalse(r['forward_allowed'])
        self.assertTrue(r['reverse_allowed'])
        self.assertEqual(dnf(r['gpr']), [frozenset(['G1', 'G2'])])

    def test_gpr_preserves_and_or(self):
        rule = {'and': [{'gene': 'A'}, {'or': [{'gene': 'B'}, {'gene': 'C'}]}]}
        self.assertEqual(set(dnf(rule)), {frozenset(['A', 'B']), frozenset(['A', 'C'])})
        self.assertRaises(ValueError, dnf, rule, 1)

    def test_projection_direction_and_complex_tiers(self):
        b = empty_builder()
        genes = [b.gene('G1'), b.gene('G2')]
        r = dict(resource='fixture', reaction_id='R',
                 substrates=[dict(id='CHEBI:1', compartment='c', stoichiometry=2)],
                 products=[dict(id='CHEBI:2', compartment='m', stoichiometry=1),
                           dict(id='CHEBI:3', compartment='m', stoichiometry=3)])
        b.reaction(r, [genes], [(True, False)], 'core')
        substrate_edges = [e for e in b.evidence if e['mechanism'] == 'substrate_participation']
        self.assertEqual({e['source'] for e in substrate_edges}, {'Metab__CHEBI:2_m', 'Metab__CHEBI:3_m'})
        self.assertTrue(all(e['first_variant'] == 'exploratory' for e in b.evidence))
        self.assertEqual(len(b.evidence), 5)  # two catalysts + two inputs + one output
        b = empty_builder()
        b.reaction(r, [[b.gene('G1')]], [(False, True)], 'lipid_curated')
        self.assertTrue(all(e['first_variant'] == 'exploratory' for e in b.evidence))

    def test_unknown_compartment_never_collapses_to_cytoplasm(self):
        b = empty_builder()
        self.assertNotEqual(b.metabolite('CHEBI:1', 'unknown'), b.metabolite('CHEBI:1', 'c'))
        self.assertEqual(chemical_id('15377'), 'CHEBI:15377')
        self.assertTrue(flag('True'))
        self.assertFalse(flag('0'))
        self.assertRaises(ValueError, flag, 'unknown')

    def test_unverified_catalyst_downgrades_entire_reaction(self):
        b = empty_builder()
        protein = b.protein('not_validated_mouse_accession')
        r = dict(resource='Rhea', reaction_id='R',
                 substrates=[dict(id='CHEBI:1', compartment='unknown', stoichiometry=1)],
                 products=[dict(id='CHEBI:2', compartment='unknown', stoichiometry=1)])
        b.reaction(r, [[protein]], [(False, False)], 'lipid_curated')
        self.assertEqual(len(b.evidence), 3)
        self.assertTrue(all(e['first_variant'] == 'exploratory' for e in b.evidence))

    def test_sign_conflicts_are_quarantined_after_gene_id_collapse(self):
        import assemble_mouse_pkn as module
        b = empty_builder()
        b.gene_map = {'P1': 'SameGene', 'P2': 'SameGene', 'P3': 'Target'}
        b.emit(b.protein('P1'), b.protein('P3'), 1, 'core', 'A', 'signaling')
        b.emit(b.protein('P2'), b.protein('P3'), -1, 'core', 'B', 'transcription')
        original = module.OUT
        with tempfile.TemporaryDirectory() as directory:
            module.OUT = Path(directory)
            try:
                b.export()
                self.assertEqual(len((module.OUT / 'pkn_core.csv').read_text().splitlines()), 1)
                self.assertEqual(json.loads((module.OUT / 'build_summary.json').read_text())
                                 ['variants']['core']['conflicting_pairs_excluded'], 1)
            finally:
                module.OUT = original

    def test_cached_stitch_is_never_reassembled_from_a_stale_regulation_table(self):
        import assemble_mouse_pkn as module
        b = empty_builder()
        original = module.OUT
        with tempfile.TemporaryDirectory() as directory:
            module.OUT = Path(directory)
            try:
                (module.OUT / 'regulation_edges.csv').write_text('resource\nSTITCH\n')
                b.regulatory()
                self.assertEqual(b.evidence, [])
                self.assertEqual(b.exclusions['cached_STITCH_raw_actor_orientation_unavailable'], 1)
            finally:
                module.OUT = original


if __name__ == '__main__':
    unittest.main()
