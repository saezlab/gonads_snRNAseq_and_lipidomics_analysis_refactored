"""Regression tests for raw STITCH direction and conflict admission rules."""

import unittest

from build_stitch_directional import (
    candidate_from_action,
    detailed_link_reason,
    preliminary_action_reason,
    project_signed_pairs,
)


def action(**overrides):
    record = {
        "item_id_a": "CIDm00000271", "item_id_b": "10090.ENSMUSP00000100212",
        "mode": "binding", "action": "activation", "a_is_acting": "t", "score": "900",
    }
    record.update(overrides)
    return record


class StitchDirectionalTests(unittest.TestCase):
    def test_chemical_actor_explicit_action_is_accepted_even_when_mode_is_binding(self):
        self.assertIsNone(preliminary_action_reason(action(), 700))

    def test_protein_actor_orientation_is_rejected(self):
        self.assertEqual(preliminary_action_reason(action(item_id_a="10090.ENSMUSP00000100212",
                                                         item_id_b="CIDm00000271"), 700),
                         "raw_actor_is_not_STITCH_chemical")

    def test_false_acting_flag_is_rejected(self):
        self.assertEqual(preliminary_action_reason(action(a_is_acting="f"), 700),
                         "raw_chemical_is_not_acting_partner")

    def test_unsigned_action_and_score_below_cutoff_are_rejected(self):
        self.assertEqual(preliminary_action_reason(action(action="binding"), 700),
                         "raw_action_is_not_explicitly_signed")
        self.assertEqual(preliminary_action_reason(action(score="699"), 700),
                         "action_score_below_threshold")

    def test_ambiguous_chemical_mapping_is_rejected(self):
        record, reason = candidate_from_action(
            action(), {"271": {"CHEBI:1", "CHEBI:2"}},
            {"10090.ENSMUSP00000100212": {"Q9TEST"}}, {"Q9TEST": "Test"},
            {"actions_uncompressed": "a", "actions_compressed": "b", "links_compressed": "c"}, 42)
        self.assertIsNone(record)
        self.assertEqual(reason, "chemical_has_ambiguous_current_CHEBI_mapping")

    def test_unmapped_protein_and_high_textmining_are_rejected(self):
        record, reason = candidate_from_action(
            action(), {"271": {"CHEBI:1"}}, {}, {},
            {"actions_uncompressed": "a", "actions_compressed": "b", "links_compressed": "c"}, 42)
        self.assertIsNone(record)
        self.assertEqual(reason, "mouse_protein_has_no_frozen_UniProt_mapping")
        self.assertEqual(detailed_link_reason({"textmining": "700"}, 700),
                         "detailed_link_textmining_at_or_above_threshold")
        self.assertEqual(detailed_link_reason({}, 700), "invalid_detailed_link_textmining_score")

    def test_opposing_signs_are_quarantined_after_node_projection(self):
        graph, conflicts, _ = project_signed_pairs([
            {"source": "Metab__CHEBI:1_unknown", "target": "Gene__Target", "interaction": 1},
            {"source": "Metab__CHEBI:1_unknown", "target": "Gene__Target", "interaction": -1},
            {"source": "Metab__CHEBI:2_unknown", "target": "Gene__Target", "interaction": 1},
        ])
        self.assertEqual(conflicts, [("Metab__CHEBI:1_unknown", "Gene__Target")])
        self.assertEqual(graph, [("Metab__CHEBI:2_unknown", 1, "Gene__Target")])


if __name__ == "__main__":
    unittest.main()
