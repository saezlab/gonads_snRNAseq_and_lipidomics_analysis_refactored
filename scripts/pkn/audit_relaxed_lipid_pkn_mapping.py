#!/usr/bin/env python3
"""Audit the duplication induced by relaxed measured-lipid-to-PKN mappings.

The frozen mouse PKN uses ChEBI-addressable metabolite nodes.  This audit asks
how many *existing* PKN metabolite nodes and reaction instances a lipidomic
feature would reach if its score were copied to mappings that are less specific
than the strict, reported-resolution policy.  It does not create edges, alter
the PKN, or fit MOON.

The tiers are intentionally separate:

* ``exact_reported_resolution``: an exact SwissLipids-to-ChEBI candidate at the
  reported nomenclature resolution;
* ``compatible_more_specific``: a SwissLipids candidate with additional
  structural detail.  This can include an sn-position assignment, double-bond
  geometry/position, or stereochemistry; it is therefore not a pure
  sn-position-only mapping;
* ``compatible_lower_resolution``: an aggregate species such as PC(42:6).
  It no longer preserves the measured acyl-chain partition;
* ``compatible_headgroup_isomer``: an unresolved headgroup isomer, for example
  glucose versus galactose; and
* ``hierarchy_class`` / ``hierarchy_ancestor``: class or ancestor candidates.
  These are deliberately broad and must never be confused with a
  molecular-species mapping.

The separate ``generic_diacyl_class_template`` tier is generated only for the
same seven explicit diacyl classes accepted by ``build_lipid_rules.py``.  It
maps a feature to a generic Rhea class participant and is reported as a broad
class-template stress test, never as an sn-relaxed molecular-species match.

Run from the repository root:

    python3 scripts/pkn/audit_relaxed_lipid_pkn_mapping.py \
      --output results/causal_modelling/relaxed_lipid_mapping_audit

All output is a new audit directory.  The default PKN is the immutable
``pkn_lipid_rules`` resource used for the Factor 1 TF-to-lipid workflow.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
from typing import Any, Iterable

# Reuse the exact diacyl eligibility and generic Rhea participant policy of the
# measured-species rule builder.  This import is standard-library-only and does
# not execute a build.  It keeps the broad-template sensitivity anchored to the
# same seven classes rather than inferring a class from a label ad hoc.
from build_lipid_rules import GENERIC as DIACYL_TEMPLATE_CHEBI
from build_lipid_rules import parse_feature as parse_diacyl_feature


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PKN = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
DEFAULT_OUTPUT = ROOT / "results/causal_modelling/relaxed_lipid_mapping_audit"

TIER_ORDER = [
    "existing_measured_species_rule",
    "exact_reported_resolution",
    "compatible_more_specific",
    "compatible_lower_resolution",
    "compatible_headgroup_isomer",
    "generic_diacyl_class_template",
    "hierarchy_class",
    "hierarchy_ancestor",
]

TIER_DESCRIPTION = {
    "existing_measured_species_rule": (
        "Existing measured-species LINEX2 rule node in pkn_lipid_rules; this is "
        "a chain-matched reaction hypothesis, not a curated ChEBI species reaction."
    ),
    "exact_reported_resolution": (
        "SwissLipids candidate with the same canonical shorthand at the reported "
        "nomenclature resolution and a ChEBI identifier."
    ),
    "compatible_more_specific": (
        "SwissLipids candidate carrying extra structural detail. It may resolve sn "
        "position, double-bond position/geometry, or stereochemistry, so it is not "
        "restricted to positional isomers."
    ),
    "compatible_lower_resolution": (
        "SwissLipids aggregate species candidate. It can preserve class and total "
        "carbon:double-bond composition, but not the measured acyl-chain partition."
    ),
    "compatible_headgroup_isomer": (
        "Compatible SwissLipids headgroup isomer, such as a glucose/galactose "
        "alternative for an unresolved hexosylceramide."
    ),
    "generic_diacyl_class_template": (
        "Generic Rhea 1,2-diacyl class participant used by the measured-species rule "
        "builder. This is a broad class-template stress test, not a species or "
        "sn-position mapping."
    ),
    "hierarchy_class": (
        "SwissLipids class-level annotation. It intentionally discards molecular "
        "species identity."
    ),
    "hierarchy_ancestor": (
        "SwissLipids non-class ancestor annotation. It is broader than the measured "
        "species and may overlap a compatible candidate."
    ),
}

# The cumulative policy columns expose different decisions rather than silently
# treating every candidate as equally valid.  The generic template branch stays
# separate from position/structure compatibility because generic Rhea terms do
# not encode the observed acyl-chain composition.
POLICIES = collections.OrderedDict([
    (
        "native_position_or_structure_relaxed",
        ["exact_reported_resolution", "compatible_more_specific"],
    ),
    (
        "strict_existing",
        ["existing_measured_species_rule", "exact_reported_resolution"],
    ),
    (
        "position_or_structure_relaxed",
        [
            "existing_measured_species_rule",
            "exact_reported_resolution",
            "compatible_more_specific",
        ],
    ),
    (
        "aggregate_composition_relaxed",
        [
            "existing_measured_species_rule",
            "exact_reported_resolution",
            "compatible_more_specific",
            "compatible_lower_resolution",
        ],
    ),
    (
        "headgroup_isomer_sensitivity",
        [
            "existing_measured_species_rule",
            "exact_reported_resolution",
            "compatible_more_specific",
            "compatible_lower_resolution",
            "compatible_headgroup_isomer",
        ],
    ),
    (
        "generic_diacyl_template_stress_test",
        [
            "existing_measured_species_rule",
            "exact_reported_resolution",
            "generic_diacyl_class_template",
        ],
    ),
    (
        "class_hierarchy_stress_test",
        [
            "existing_measured_species_rule",
            "exact_reported_resolution",
            "hierarchy_class",
            "hierarchy_ancestor",
        ],
    ),
    (
        "all_broad_candidates_stress_test",
        list(TIER_ORDER),
    ),
])


def read_csv(path: Path, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_lines(path: Path) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    with path.open() as handle:
        for line in handle:
            if line.strip():
                row = json.loads(line)
                result[row["node_id"]] = row
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truth(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1"}


def pipe(values: Iterable[str]) -> str:
    return "|".join(sorted({value for value in values if value}))


def percentile(values: list[int], probability: float) -> float:
    """Linear percentile without a numerical dependency."""
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return float(ordered[0])
    location = (len(ordered) - 1) * probability
    lower = math.floor(location)
    upper = math.ceil(location)
    if lower == upper:
        return float(ordered[lower])
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (location - lower)


def summary_statistics(values: list[int]) -> dict[str, float | int]:
    if not values:
        return {"min": 0, "p25": 0.0, "median": 0.0, "p75": 0.0,
                "p90": 0.0, "p95": 0.0, "p99": 0.0, "max": 0,
                "mean": 0.0}
    return {
        "min": min(values),
        "p25": round(percentile(values, 0.25), 3),
        "median": round(statistics.median(values), 3),
        "p75": round(percentile(values, 0.75), 3),
        "p90": round(percentile(values, 0.90), 3),
        "p95": round(percentile(values, 0.95), 3),
        "p99": round(percentile(values, 0.99), 3),
        "max": max(values),
        "mean": round(statistics.mean(values), 3),
    }


def candidate_rows(features: list[dict[str, str]], mapping_rows: list[dict[str, str]],
                   active_metabolite_nodes: set[str]) -> list[dict[str, Any]]:
    """Make one provenance-preserving candidate identity row per tier and ChEBI ID."""
    feature_lookup = {row["feature_id"]: row for row in features}
    if len(feature_lookup) != len(features):
        raise AssertionError("mapping_features.csv must contain one row per feature")

    records: dict[tuple[str, str, str], dict[str, Any]] = {}

    def add(feature_id: str, tier: str, chemical_id: str, origin: str,
            swisslipids_id: str = "", matched_abbreviation: str = "",
            candidate_canonical: str = "", candidate_level: str = "",
            match_type: str = "") -> None:
        if tier not in TIER_ORDER:
            raise AssertionError("Unknown mapping tier: " + tier)
        key = (feature_id, tier, chemical_id)
        record = records.setdefault(key, {
            "feature_id": feature_id,
            "mapping_tier": tier,
            "mapping_tier_description": TIER_DESCRIPTION[tier],
            "candidate_chemical_id": chemical_id,
            "candidate_origin": set(),
            "match_types": set(),
            "swisslipids_ids": set(),
            "matched_abbreviations": set(),
            "candidate_canonical_names": set(),
            "candidate_nomenclature_levels": set(),
        })
        record["candidate_origin"].add(origin)
        record["match_types"].add(match_type)
        record["swisslipids_ids"].add(swisslipids_id)
        record["matched_abbreviations"].add(matched_abbreviation)
        record["candidate_canonical_names"].add(candidate_canonical)
        record["candidate_nomenclature_levels"].add(candidate_level)

    # Only direct SwissLipids identifiers are considered for the unnormalised
    # pkn_lipid_rules. pH 7.3 aliases belong to the separate pH-normalized PKN
    # sensitivity and are intentionally not imported here.
    accepted_match_types = {
        "exact_reported_resolution",
        "compatible_more_specific",
        "compatible_lower_resolution",
        "compatible_headgroup_isomer",
        "hierarchy_class",
        "hierarchy_ancestor",
    }
    for row in mapping_rows:
        match_type = row["match_type"]
        if match_type not in accepted_match_types:
            raise AssertionError("Unexpected crosswalk match type: " + match_type)
        for chemical_id in filter(None, row["chebi_ids"].split("|")):
            add(
                row["feature_id"], match_type, chemical_id, "SwissLipids_to_ChEBI",
                swisslipids_id=row["swisslipids_id"],
                matched_abbreviation=row["matched_abbreviation"],
                candidate_canonical=row["candidate_canonical"],
                candidate_level=row["candidate_nomenclature_level"],
                match_type=match_type,
            )

    # Existing measured-species rule nodes provide the strict rule-hypothesis
    # baseline. A generic template mapping is admitted only for a feature that
    # passes the same explicit two-chain eligibility predicate as the rule build.
    for feature in features:
        feature_id = feature["feature_id"]
        measured_id = "MeasuredLipid:" + feature_id
        if measured_id in active_metabolite_nodes:
            add(feature_id, "existing_measured_species_rule", measured_id,
                "pkn_lipid_rules_existing_measured_species_node",
                matched_abbreviation=feature["canonical_name"],
                candidate_canonical=feature["canonical_name"],
                candidate_level=feature["nomenclature_level"],
                match_type="existing_measured_species_rule")
        parsed = parse_diacyl_feature(feature)
        if parsed is None:
            continue
        for chemical_id, lipid_class in DIACYL_TEMPLATE_CHEBI.items():
            if lipid_class == parsed["class"]:
                add(feature_id, "generic_diacyl_class_template", chemical_id,
                    "build_lipid_rules_GENERIC_Rhea_template",
                    matched_abbreviation=feature["canonical_name"],
                    candidate_canonical=feature["canonical_name"],
                    candidate_level=feature["nomenclature_level"],
                    match_type="generic_diacyl_class_template")

    result = []
    for record in records.values():
        result.append({
            "feature_id": record["feature_id"],
            "mapping_tier": record["mapping_tier"],
            "mapping_tier_description": record["mapping_tier_description"],
            "candidate_chemical_id": record["candidate_chemical_id"],
            "candidate_origin": pipe(record["candidate_origin"]),
            "match_types": pipe(record["match_types"]),
            "swisslipids_ids": pipe(record["swisslipids_ids"]),
            "matched_abbreviations": pipe(record["matched_abbreviations"]),
            "candidate_canonical_names": pipe(record["candidate_canonical_names"]),
            "candidate_nomenclature_levels": pipe(record["candidate_nomenclature_levels"]),
        })
    return sorted(result, key=lambda row: (row["feature_id"], TIER_ORDER.index(row["mapping_tier"]),
                                            row["candidate_chemical_id"]))


def node_neighbourhoods(edges: list[dict[str, str]], nodes: dict[str, dict[str, Any]],
                        active_metabolites: set[str]) -> dict[str, dict[str, set[str]]]:
    """Summarise reaction and non-reaction adjacency for each active metabolite."""
    result = {
        node: {
            "reaction_instances": set(),
            "product_reaction_instances": set(),
            "substrate_reaction_instances": set(),
            "nonreaction_neighbours": set(),
            "incident_edges": set(),
        }
        for node in active_metabolites
    }
    for edge in edges:
        source, target = edge["source"], edge["target"]
        signed = (source, edge["interaction"], target)
        for metabolite, other, direction in ((source, target, "outgoing"), (target, source, "incoming")):
            if metabolite not in result:
                continue
            result[metabolite]["incident_edges"].add(signed)
            if nodes[other]["node_type"] == "enzyme_reaction_instance":
                result[metabolite]["reaction_instances"].add(other)
                if direction == "incoming":
                    result[metabolite]["product_reaction_instances"].add(other)
                else:
                    result[metabolite]["substrate_reaction_instances"].add(other)
            else:
                result[metabolite]["nonreaction_neighbours"].add(other)
    return result


def union_metrics(node_ids: set[str], nodes: dict[str, dict[str, Any]],
                  neighbourhoods: dict[str, dict[str, set[str]]]) -> dict[str, Any]:
    metrics: dict[str, set[str]] = {
        "chemical_ids": set(),
        "reaction_instances": set(),
        "product_reaction_instances": set(),
        "substrate_reaction_instances": set(),
        "nonreaction_neighbours": set(),
        "incident_edges": set(),
    }
    for node in node_ids:
        metrics["chemical_ids"].add(nodes[node]["chemical_id"])
        for key in metrics:
            if key != "chemical_ids":
                metrics[key].update(neighbourhoods[node][key])
    product_nodes = {
        node for node in node_ids if neighbourhoods[node]["product_reaction_instances"]
    }
    substrate_nodes = {
        node for node in node_ids if neighbourhoods[node]["substrate_reaction_instances"]
    }
    return {
        "n_candidate_metabolite_nodes": len(node_ids),
        "n_candidate_product_metabolite_nodes": len(product_nodes),
        "n_candidate_substrate_metabolite_nodes": len(substrate_nodes),
        "n_candidate_both_role_metabolite_nodes": len(product_nodes & substrate_nodes),
        "n_candidate_chemical_ids": len(metrics["chemical_ids"]),
        "n_candidate_reaction_instances": len(metrics["reaction_instances"]),
        "n_candidate_product_reaction_instances": len(metrics["product_reaction_instances"]),
        "n_candidate_substrate_reaction_instances": len(metrics["substrate_reaction_instances"]),
        "n_candidate_nonreaction_neighbours": len(metrics["nonreaction_neighbours"]),
        "n_candidate_incident_pkn_edges": len(metrics["incident_edges"]),
        "candidate_metabolite_nodes": pipe(node_ids),
        "candidate_product_metabolite_nodes": pipe(product_nodes),
        "candidate_substrate_metabolite_nodes": pipe(substrate_nodes),
        "candidate_chemical_ids": pipe(metrics["chemical_ids"]),
        "candidate_reaction_instance_nodes": pipe(metrics["reaction_instances"]),
    }


def validate_inputs(pkn: Path) -> list[Path]:
    paths = [
        pkn / "pkn_lipid_rules.csv",
        pkn / "nodes.jsonl",
        pkn / "mapping_features.csv",
        pkn / "mapping_candidates.csv",
        ROOT / "scripts/pkn/build_lipid_rules.py",
        Path(__file__).resolve(),
    ]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing audit dependency: " + ", ".join(missing))
    return paths


def self_test() -> None:
    assert POLICIES["strict_existing"] == [
        "existing_measured_species_rule", "exact_reported_resolution"
    ]
    assert POLICIES["native_position_or_structure_relaxed"] == [
        "exact_reported_resolution", "compatible_more_specific"
    ]
    assert "compatible_more_specific" in POLICIES["position_or_structure_relaxed"]
    assert "compatible_lower_resolution" not in POLICIES["position_or_structure_relaxed"]
    assert "generic_diacyl_class_template" not in POLICIES["position_or_structure_relaxed"]
    assert "generic_diacyl_class_template" in POLICIES["generic_diacyl_template_stress_test"]
    assert set(POLICIES["all_broad_candidates_stress_test"]) == set(TIER_ORDER)
    assert percentile([0, 10], 0.5) == 5.0
    test_nodes = {
        "M": {"chemical_id": "CHEBI:1", "node_type": "metabolite"},
        "E": {"node_type": "enzyme_reaction_instance"},
        "G": {"node_type": "gene_product"},
    }
    test_edges = [
        {"source": "M", "interaction": "1", "target": "E"},
        {"source": "E", "interaction": "1", "target": "M"},
        {"source": "M", "interaction": "1", "target": "G"},
    ]
    neighbourhood = node_neighbourhoods(test_edges, test_nodes, {"M"})["M"]
    assert neighbourhood["product_reaction_instances"] == {"E"}
    assert neighbourhood["substrate_reaction_instances"] == {"E"}
    assert neighbourhood["nonreaction_neighbours"] == {"G"}
    metrics = union_metrics({"M"}, test_nodes, {"M": neighbourhood})
    assert metrics["n_candidate_product_metabolite_nodes"] == 1
    assert metrics["n_candidate_substrate_metabolite_nodes"] == 1
    print("Relaxed lipid mapping audit self-tests passed.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pkn", type=Path, default=DEFAULT_PKN,
                        help="Frozen PKN directory (default: causal_modelling/cosmos_mouse_pkn).")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="New or empty audit directory.")
    parser.add_argument("--self-test", action="store_true",
                        help="Run semantic unit checks without reading the PKN.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite audit files in an existing output directory.")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return

    pkn = args.pkn.resolve()
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or (any(output.iterdir()) and not args.overwrite)):
        parser.error("Output directory is occupied: " + str(output))
    dependencies = validate_inputs(pkn)

    edges = read_csv(pkn / "pkn_lipid_rules.csv")
    if not edges or set(edges[0]) != {"source", "interaction", "target"}:
        raise AssertionError("pkn_lipid_rules.csv must contain source, interaction and target only")
    signed_edges = {(row["source"], row["interaction"], row["target"]) for row in edges}
    if len(signed_edges) != len(edges):
        raise AssertionError("pkn_lipid_rules.csv contains duplicate signed triples")
    if {row["interaction"] for row in edges} - {"-1", "1"}:
        raise AssertionError("pkn_lipid_rules.csv has an invalid sign")

    nodes = json_lines(pkn / "nodes.jsonl")
    active_node_ids = {endpoint for edge in edges for endpoint in (edge["source"], edge["target"])}
    if not active_node_ids <= set(nodes):
        raise AssertionError("Active PKN node absent from nodes.jsonl")
    active_metabolites = {
        node for node in active_node_ids if nodes[node]["node_type"] == "metabolite"
    }
    metabolite_nodes_by_chemical: dict[str, set[str]] = collections.defaultdict(set)
    for node in active_metabolites:
        metabolite_nodes_by_chemical[nodes[node]["chemical_id"]].add(node)
    neighbourhoods = node_neighbourhoods(edges, nodes, active_metabolites)

    features = read_csv(pkn / "mapping_features.csv")
    if len({row["feature_id"] for row in features}) != len(features):
        raise AssertionError("mapping_features.csv must have unique feature IDs")
    mapping = read_csv(pkn / "mapping_candidates.csv")
    if not {row["feature_id"] for row in mapping} <= {row["feature_id"] for row in features}:
        raise AssertionError("mapping_candidates.csv has a feature outside mapping_features.csv")
    raw_candidates = candidate_rows(features, mapping, set(metabolite_nodes_by_chemical))

    # Candidate chemical identity inventory: include zero-node identities so the
    # audit separates an unsupported chemical identity from a low-duplication
    # identity.  In particular, many SwissLipids positional records have no
    # ChEBI handle and thus cannot be PKN nodes without inventing an identifier.
    chemical_inventory: list[dict[str, Any]] = []
    node_inventory: list[dict[str, Any]] = []
    reaction_inventory: list[dict[str, Any]] = []
    tier_nodes: dict[str, dict[str, set[str]]] = collections.defaultdict(lambda: collections.defaultdict(set))
    feature_tier_identities: dict[str, dict[str, set[str]]] = collections.defaultdict(lambda: collections.defaultdict(set))
    for candidate in raw_candidates:
        feature_id = candidate["feature_id"]
        tier = candidate["mapping_tier"]
        chemical_id = candidate["candidate_chemical_id"]
        pkn_nodes = sorted(metabolite_nodes_by_chemical.get(chemical_id, set()))
        feature_tier_identities[feature_id][tier].add(chemical_id)
        tier_nodes[feature_id][tier].update(pkn_nodes)
        local_metrics = union_metrics(set(pkn_nodes), nodes, neighbourhoods)
        chemical_inventory.append({
            **candidate,
            "n_active_pkn_metabolite_nodes": len(pkn_nodes),
            "active_pkn_metabolite_nodes": pipe(pkn_nodes),
            "n_active_pkn_reaction_instances": local_metrics["n_candidate_reaction_instances"],
            "n_active_pkn_product_metabolite_nodes": local_metrics["n_candidate_product_metabolite_nodes"],
            "n_active_pkn_substrate_metabolite_nodes": local_metrics["n_candidate_substrate_metabolite_nodes"],
            "n_active_pkn_product_reaction_instances": local_metrics["n_candidate_product_reaction_instances"],
            "n_active_pkn_substrate_reaction_instances": local_metrics["n_candidate_substrate_reaction_instances"],
            "n_active_pkn_nonreaction_neighbours": local_metrics["n_candidate_nonreaction_neighbours"],
            "n_active_pkn_incident_edges": local_metrics["n_candidate_incident_pkn_edges"],
        })
        for node in pkn_nodes:
            local = union_metrics({node}, nodes, neighbourhoods)
            node_inventory.append({
                **candidate,
                "pkn_metabolite_node": node,
                "pkn_node_label": nodes[node].get("label", ""),
                "compartment": nodes[node].get("compartment", ""),
                **local,
            })
            roles: dict[str, set[str]] = collections.defaultdict(set)
            for reaction in neighbourhoods[node]["product_reaction_instances"]:
                roles[reaction].add("product")
            for reaction in neighbourhoods[node]["substrate_reaction_instances"]:
                roles[reaction].add("substrate")
            for reaction, role_set in sorted(roles.items()):
                reaction_node = nodes[reaction]
                reaction_inventory.append({
                    "feature_id": feature_id,
                    "mapping_tier": tier,
                    "candidate_chemical_id": chemical_id,
                    "pkn_metabolite_node": node,
                    "reaction_instance_node": reaction,
                    "reaction_role": pipe(role_set),
                    "reaction_id": reaction_node.get("reaction_id", ""),
                    "reaction_resource": reaction_node.get("resource", ""),
                    "reaction_reverse": reaction_node.get("reverse", ""),
                    "reaction_catalysts": pipe(reaction_node.get("catalysts", [])),
                })

    # One long row per feature/policy makes the duplicate score-assignment
    # burden directly queryable without asking users to reconstruct a union of
    # mapping tiers.  The node-level inventory retains the provenance needed to
    # audit every member of that union.
    feature_policy_rows: list[dict[str, Any]] = []
    policy_to_node_features: dict[str, dict[str, set[str]]] = {
        policy: collections.defaultdict(set) for policy in POLICIES
    }
    for feature in features:
        feature_id = feature["feature_id"]
        for policy, tiers in POLICIES.items():
            node_ids = set().union(*(tier_nodes[feature_id][tier] for tier in tiers))
            identity_ids = set().union(*(feature_tier_identities[feature_id][tier] for tier in tiers))
            metrics = union_metrics(node_ids, nodes, neighbourhoods)
            for node in node_ids:
                policy_to_node_features[policy][node].add(feature_id)
            feature_policy_rows.append({
                "feature_id": feature_id,
                "in_lipid_de": feature["in_lipid_de"],
                "in_historical_factor1_64": feature["in_factor1"],
                "canonical_name": feature["canonical_name"],
                "lipid_class": feature["class"],
                "nomenclature_level": feature["nomenclature_level"],
                "best_mapping_status": feature["best_mapping_status"],
                "policy": policy,
                "included_mapping_tiers": pipe(tiers),
                "n_candidate_chemical_id_records": len(identity_ids),
                **metrics,
            })

    policy_row = {(row["feature_id"], row["policy"]): row for row in feature_policy_rows}
    checks = {
        "every_feature_has_every_policy": len(policy_row) == len(features) * len(POLICIES),
        "candidate_node_inventory_only_contains_active_metabolites": all(
            row["pkn_metabolite_node"] in active_metabolites for row in node_inventory
        ),
        "reaction_inventory_only_contains_active_reaction_instances": all(
            row["reaction_instance_node"] in active_node_ids and
            nodes[row["reaction_instance_node"]]["node_type"] == "enzyme_reaction_instance"
            for row in reaction_inventory
        ),
        "product_candidate_metabolite_counts_do_not_exceed_total": all(
            row["n_candidate_product_metabolite_nodes"] <= row["n_candidate_metabolite_nodes"]
            for row in feature_policy_rows
        ),
        "position_or_structure_policy_superset_of_strict_existing": all(
            policy_row[(feature["feature_id"], "position_or_structure_relaxed")]["n_candidate_metabolite_nodes"] >=
            policy_row[(feature["feature_id"], "strict_existing")]["n_candidate_metabolite_nodes"]
            for feature in features
        ),
        "aggregate_composition_policy_superset_of_position_or_structure": all(
            policy_row[(feature["feature_id"], "aggregate_composition_relaxed")]["n_candidate_metabolite_nodes"] >=
            policy_row[(feature["feature_id"], "position_or_structure_relaxed")]["n_candidate_metabolite_nodes"]
            for feature in features
        ),
        "headgroup_sensitivity_policy_superset_of_aggregate_composition": all(
            policy_row[(feature["feature_id"], "headgroup_isomer_sensitivity")]["n_candidate_metabolite_nodes"] >=
            policy_row[(feature["feature_id"], "aggregate_composition_relaxed")]["n_candidate_metabolite_nodes"]
            for feature in features
        ),
    }
    if not all(checks.values()):
        raise AssertionError("Relaxed mapping audit check failed: " + json.dumps(checks, sort_keys=True))

    # Tier-specific summary makes the important distinction requested by the
    # user visible before any cumulative policy is chosen.
    tier_summary_rows: list[dict[str, Any]] = []
    for tier in TIER_ORDER:
        rows = [row for row in chemical_inventory if row["mapping_tier"] == tier]
        feature_ids = {row["feature_id"] for row in rows}
        active_rows = [row for row in rows if row["n_active_pkn_metabolite_nodes"]]
        tier_nodes_union = {node for row in node_inventory if row["mapping_tier"] == tier
                            for node in row["candidate_metabolite_nodes"].split("|") if node}
        reaction_union = {node for row in reaction_inventory if row["mapping_tier"] == tier
                          for node in [row["reaction_instance_node"]]}
        tier_summary_rows.append({
            "mapping_tier": tier,
            "mapping_tier_description": TIER_DESCRIPTION[tier],
            "candidate_identity_records": len(rows),
            "features_with_candidate_identity": len(feature_ids),
            "candidate_identity_records_active_in_pkn": len(active_rows),
            "features_with_active_pkn_node": len({row["feature_id"] for row in active_rows}),
            "unique_active_pkn_metabolite_nodes": len(tier_nodes_union),
            "unique_active_pkn_reaction_instances": len(reaction_union),
        })

    # Aggregate all-panel and the historical 64-feature subset separately.  The
    # latter retains the legacy mapping_features flag and is not a claim that
    # the remaining Factor 1 weights are absent.
    population_summary: list[dict[str, Any]] = []
    explicit_chain_levels = {
        "MOLECULAR_SPECIES", "SN_POSITION", "FULL_STRUCTURE", "COMPLETE_STRUCTURE"
    }
    for population, predicate in [
        ("all_1254_measured_lipids", lambda row: True),
        ("explicit_acyl_chain_composition_features", lambda row:
         row["nomenclature_level"] in explicit_chain_levels),
        ("underscore_molecular_species_features", lambda row:
         row["nomenclature_level"] == "MOLECULAR_SPECIES" and
         "_" in row["canonical_name"] and "/" not in row["canonical_name"]),
        ("historical_factor1_64_subset", lambda row: truth(row["in_historical_factor1_64"])),
    ]:
        for policy in POLICIES:
            rows = [row for row in feature_policy_rows if row["policy"] == policy and predicate(row)]
            node_assignments = sum(row["n_candidate_metabolite_nodes"] for row in rows)
            reaction_assignments = sum(row["n_candidate_reaction_instances"] for row in rows)
            used_nodes = policy_to_node_features[policy]
            if population != "all_1254_measured_lipids":
                selected_features = {row["feature_id"] for row in rows}
                used_nodes = {
                    node: feature_ids & selected_features
                    for node, feature_ids in used_nodes.items()
                }
                used_nodes = {node: feature_ids for node, feature_ids in used_nodes.items() if feature_ids}
            reused_counts = [len(feature_ids) for feature_ids in used_nodes.values()]
            metrics = summary_statistics([row["n_candidate_reaction_instances"] for row in rows])
            node_metrics = summary_statistics([row["n_candidate_metabolite_nodes"] for row in rows])
            population_summary.append({
                "population": population,
                "policy": policy,
                "features": len(rows),
                "features_with_any_active_pkn_node": sum(row["n_candidate_metabolite_nodes"] > 0 for row in rows),
                "features_with_any_product_side_pkn_node": sum(
                    row["n_candidate_product_metabolite_nodes"] > 0 for row in rows
                ),
                "features_with_any_candidate_reaction": sum(row["n_candidate_reaction_instances"] > 0 for row in rows),
                "feature_to_metabolite_node_assignments": node_assignments,
                "feature_to_product_metabolite_node_assignments": sum(
                    row["n_candidate_product_metabolite_nodes"] for row in rows
                ),
                "unique_pkn_metabolite_nodes": len(used_nodes),
                "metabolite_assignment_to_unique_node_ratio": round(node_assignments / len(used_nodes), 3) if used_nodes else 0.0,
                "maximum_features_collapsed_to_one_pkn_node": max(reused_counts, default=0),
                "feature_to_reaction_instance_assignments": reaction_assignments,
                **{f"reaction_instances_{key}": value for key, value in metrics.items()},
                **{f"metabolite_nodes_{key}": value for key, value in node_metrics.items()},
                **{f"product_metabolite_nodes_{key}": value for key, value in summary_statistics([
                    row["n_candidate_product_metabolite_nodes"] for row in rows
                ]).items()},
            })

    reuse_rows: list[dict[str, Any]] = []
    for policy, mapping_by_node in policy_to_node_features.items():
        for node, feature_ids in sorted(mapping_by_node.items()):
            reuse_rows.append({
                "policy": policy,
                "pkn_metabolite_node": node,
                "chemical_id": nodes[node]["chemical_id"],
                "node_label": nodes[node].get("label", ""),
                "compartment": nodes[node].get("compartment", ""),
                "n_measured_features_assigned": len(feature_ids),
                "feature_ids": pipe(feature_ids),
                **union_metrics({node}, nodes, neighbourhoods),
            })

    output.mkdir(parents=True, exist_ok=True)
    write_csv(
        output / "candidate_chemical_identity_inventory.csv", chemical_inventory,
        [
            "feature_id", "mapping_tier", "mapping_tier_description", "candidate_chemical_id",
            "candidate_origin", "match_types", "swisslipids_ids", "matched_abbreviations",
            "candidate_canonical_names", "candidate_nomenclature_levels",
            "n_active_pkn_metabolite_nodes", "active_pkn_metabolite_nodes",
            "n_active_pkn_reaction_instances", "n_active_pkn_product_reaction_instances",
            "n_active_pkn_substrate_reaction_instances", "n_active_pkn_product_metabolite_nodes",
            "n_active_pkn_substrate_metabolite_nodes", "n_active_pkn_nonreaction_neighbours",
            "n_active_pkn_incident_edges",
        ],
    )
    write_csv(
        output / "candidate_pkn_node_inventory.csv", node_inventory,
        [
            "feature_id", "mapping_tier", "mapping_tier_description", "candidate_chemical_id",
            "candidate_origin", "match_types", "swisslipids_ids", "matched_abbreviations",
            "candidate_canonical_names", "candidate_nomenclature_levels", "pkn_metabolite_node",
            "pkn_node_label", "compartment", "n_candidate_metabolite_nodes",
            "n_candidate_product_metabolite_nodes", "n_candidate_substrate_metabolite_nodes",
            "n_candidate_both_role_metabolite_nodes",
            "n_candidate_chemical_ids", "n_candidate_reaction_instances",
            "n_candidate_product_reaction_instances", "n_candidate_substrate_reaction_instances",
            "n_candidate_nonreaction_neighbours", "n_candidate_incident_pkn_edges",
            "candidate_metabolite_nodes", "candidate_product_metabolite_nodes",
            "candidate_substrate_metabolite_nodes", "candidate_chemical_ids",
            "candidate_reaction_instance_nodes",
        ],
    )
    write_csv(
        output / "candidate_reaction_instance_inventory.csv", reaction_inventory,
        [
            "feature_id", "mapping_tier", "candidate_chemical_id", "pkn_metabolite_node",
            "reaction_instance_node", "reaction_role", "reaction_id", "reaction_resource",
            "reaction_reverse", "reaction_catalysts",
        ],
    )
    write_csv(
        output / "per_feature_policy_summary.csv", feature_policy_rows,
        [
            "feature_id", "in_lipid_de", "in_historical_factor1_64", "canonical_name",
            "lipid_class", "nomenclature_level", "best_mapping_status", "policy",
            "included_mapping_tiers", "n_candidate_chemical_id_records",
            "n_candidate_metabolite_nodes", "n_candidate_product_metabolite_nodes",
            "n_candidate_substrate_metabolite_nodes", "n_candidate_both_role_metabolite_nodes",
            "n_candidate_chemical_ids",
            "n_candidate_reaction_instances", "n_candidate_product_reaction_instances",
            "n_candidate_substrate_reaction_instances", "n_candidate_nonreaction_neighbours",
            "n_candidate_incident_pkn_edges", "candidate_metabolite_nodes",
            "candidate_product_metabolite_nodes", "candidate_substrate_metabolite_nodes",
            "candidate_chemical_ids", "candidate_reaction_instance_nodes",
        ],
    )
    write_csv(
        output / "mapping_tier_summary.csv", tier_summary_rows,
        [
            "mapping_tier", "mapping_tier_description", "candidate_identity_records",
            "features_with_candidate_identity", "candidate_identity_records_active_in_pkn",
            "features_with_active_pkn_node", "unique_active_pkn_metabolite_nodes",
            "unique_active_pkn_reaction_instances",
        ],
    )
    write_csv(
        output / "population_duplication_summary.csv", population_summary,
        list(population_summary[0]),
    )
    write_csv(
        output / "pkn_node_reuse_by_policy.csv", reuse_rows,
        [
            "policy", "pkn_metabolite_node", "chemical_id", "node_label", "compartment",
            "n_measured_features_assigned", "feature_ids", "n_candidate_metabolite_nodes",
            "n_candidate_product_metabolite_nodes", "n_candidate_substrate_metabolite_nodes",
            "n_candidate_both_role_metabolite_nodes",
            "n_candidate_chemical_ids", "n_candidate_reaction_instances",
            "n_candidate_product_reaction_instances", "n_candidate_substrate_reaction_instances",
            "n_candidate_nonreaction_neighbours", "n_candidate_incident_pkn_edges",
            "candidate_metabolite_nodes", "candidate_product_metabolite_nodes",
            "candidate_substrate_metabolite_nodes", "candidate_chemical_ids",
            "candidate_reaction_instance_nodes",
        ],
    )
    write_csv(
        output / "dependency_hashes.csv",
        [{"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for path in dependencies],
        ["path", "sha256"],
    )

    summary = {
        "scope": (
            "Potential score duplication across existing pkn_lipid_rules metabolite nodes and "
            "their adjacent reaction instances. No PKN edge or node was created, and no MOON "
            "analysis was fitted."
        ),
        "pkn": str(pkn.relative_to(ROOT)),
        "pkn_edges": len(edges),
        "pkn_active_nodes": len(active_node_ids),
        "pkn_active_metabolite_nodes": len(active_metabolites),
        "feature_count": len(features),
        "tier_definitions": TIER_DESCRIPTION,
        "cumulative_policies": {
            policy: {
                "tiers": tiers,
                "interpretation": (
                    "A score would be copied to every distinct active PKN metabolite node "
                    "supported by these tiers; counts are a sensitivity burden, not evidence "
                    "that all assignments are true."
                ),
            }
            for policy, tiers in POLICIES.items()
        },
        "important_limitations": [
            "compatible_more_specific includes more than sn-position ambiguity; it can add stereochemical or double-bond-detail candidates.",
            "compatible_lower_resolution discards the observed acyl-chain partition and is not a positional-isomer policy.",
            "generic_diacyl_class_template and hierarchy tiers are broad class-level stress tests, not species-resolved mappings.",
            "The frozen mapping intentionally does not exhaustively retain detailed SwissLipids isomers without a ChEBI identifier; such records cannot be assigned to the current ChEBI-addressed PKN without inventing nodes.",
            "This audit uses the unnormalised pkn_lipid_rules; pH 7.3 aliases belong to the separate pkn_lipid_curated_ph73 sensitivity resource.",
            "Adjacent reaction instances describe qualitative PKN traversal. They do not establish substrate specificity, tissue flux, or biochemical activity in the samples.",
        ],
        "checks": checks,
        "tier_summary": tier_summary_rows,
        "population_duplication_summary": population_summary,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Relaxed lipid-to-PKN mapping audit\n\n"
        "This audit quantifies the noise introduced if a measured lipid score were copied to "
        "multiple existing PKN metabolite nodes. It reads the frozen `pkn_lipid_rules` graph; "
        "it does not change the PKN or run MOON.\n\n"
        "## Main files\n\n"
        "- `per_feature_policy_summary.csv`: one row per measured feature and explicit mapping policy; "
        "use this for candidate metabolite nodes, product-side candidate metabolite nodes, reaction "
        "instances, product-side reaction instances, and substrate-side reaction instances.\n"
        "- `candidate_chemical_identity_inventory.csv`: every source candidate identity, including "
        "identities that have no active PKN node.\n"
        "- `candidate_pkn_node_inventory.csv`: each supported feature-to-metabolite-node assignment.\n"
        "- `candidate_reaction_instance_inventory.csv`: each feature-to-candidate-reaction link and "
        "whether the lipid node is a reaction product or substrate.\n"
        "- `pkn_node_reuse_by_policy.csv`: identifies generic nodes to which many measured features "
        "would collapse.\n"
        "- `population_duplication_summary.csv` and `mapping_tier_summary.csv`: panel-level burden.\n\n"
        "## Mapping tiers\n\n"
        "The positional/structural tier (`compatible_more_specific`) is kept apart from aggregate "
        "composition, generic diacyl template, and hierarchy mappings. It can include more than "
        "sn position, such as double-bond or stereochemical detail. A generic Rhea class term is "
        "reported only in `generic_diacyl_template_stress_test`; it is not treated as a "
        "species-specific or sn-relaxed assignment.\n\n"
        "`native_position_or_structure_relaxed` reports only direct ChEBI-addressable candidates "
        "and deliberately excludes the existing measured-species rule nodes. The "
        "`underscore_molecular_species_features` population is the narrowest relevant subset for "
        "the question of unreported sn position: it includes only molecular-species measurements "
        "with unordered acyl chains.\n\n"
        "The audit is exhaustive over the candidate chemical identities addressable in the frozen "
        "SwissLipids-to-ChEBI mapping and the active `pkn_lipid_rules` nodes. It is not an exhaustive "
        "enumeration of all chemically possible positional isomers: detailed SwissLipids records "
        "without a ChEBI ID cannot be mapped to this PKN without adding new, uncurated nodes.\n\n"
        "## MOON scoring implication\n\n"
        "Candidate mappings must not be passed to MOON as independent copies of the same lipid score: "
        "the ULM score treats duplicated rows as separate observations and would over-weight an "
        "ambiguous feature. A relaxed sensitivity should instead retain the feature-to-candidate "
        "mapping table and evaluate one supported node per feature in repeated mapping draws, or use "
        "an explicitly weighted aggregation whose interpretation is validated separately.\n"
    )

    headline = [row for row in population_summary if row["population"] == "all_1254_measured_lipids"]
    print(json.dumps({
        "output": str(output),
        "features": len(features),
        "candidate_identity_rows": len(chemical_inventory),
        "candidate_pkn_node_rows": len(node_inventory),
        "candidate_reaction_rows": len(reaction_inventory),
        "all_panel_policy_summary": headline,
    }, indent=2))


if __name__ == "__main__":
    main()
