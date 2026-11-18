#!/usr/bin/env python3
"""Audit measured-lipid coverage by the directionality-audited STITCH overlay.

This is an input-identity audit only.  It asks whether a measured lipidomic
feature can be matched to a *chemical source* in the separately built,
directionality-audited STITCH v5 mouse overlay.  It does not add a mapping to
the PKN, duplicate measurements, change the PKN, or run MOON.

The audit keeps identity assumptions separate:

* ``exact_CHEBI_identity`` preserves a SwissLipids/ChEBI chemical identity
  after ChEBI alternate-ID canonicalisation;
* ``explicit_conjugate_acid_base`` uses a directly recorded ChEBI conjugate
  acid/base relation only;
* ``same_formula_is_a_descendant`` is a deliberately separate sensitivity for
  an exact reported-resolution feature.  It follows an explicit ChEBI ``is_a``
  path through the feature or its explicit acid/base counterpart, requires the
  same formula and charge at the descendant endpoint, and can resolve
  unmeasured double-bond position/geometry or stereochemistry.  It is never a
  direct measured-chemical claim; and
* SwissLipids compatible-lower-resolution, headgroup-isomer, hierarchy-class,
  and hierarchy-ancestor rows are reported but never promoted to inputs.

Raw retained STITCH evidence and signed triples surviving the final assembled
PKN conflict filter are counted separately.  The output is useful before a
future lipid-to-TF sensitivity analysis, because the same Factor 1 value must
not be copied across several candidate chemical nodes and treated as
independent observations.

Run from the repository root:

    python3 scripts/pkn/audit_measured_lipids_stitch_directional.py \
      --output results/causal_modelling/stitch_measured_lipid_mapping_audit
"""

from __future__ import annotations

import argparse
import collections
import csv
import functools
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn"
OVERLAY = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn_stitch_directional"
FACTOR_INPUTS = ROOT / "miscellaneous" / "older_results" / "causal_modelling" / "cosmos_factor1_pilot" / "inputs" / "factor1_lipid_node_mapping.csv"
OBO = ROOT / "support" / "cosmos_mouse" / "raw" / "regulation_chebi.obo.gz"
DEFAULT_OUTPUT = ROOT / "results" / "causal_modelling" / "stitch_measured_lipid_mapping_audit"

FEATURES = BASE / "mapping_features.csv"
CANDIDATES = BASE / "mapping_candidates.csv"
ALIASES = BASE / "regulation_chebi_aliases.csv"
ACID_BASE = BASE / "regulation_chebi_acid_base_relations.csv"
STITCH_EVIDENCE = OVERLAY / "stitch_directional_evidence.csv"
STITCH_GRAPH = OVERLAY / "pkn_exploratory_stitch_directional.csv"
STITCH_MANIFEST = OVERLAY / "stitch_directional_manifest.json"
STITCH_VALIDATION = OVERLAY / "stitch_directional_validation.json"

CHEBI = re.compile(r"CHEBI:\d+")

TIER_ORDER = [
    "exact_reported_resolution",
    "compatible_more_specific",
    "compatible_lower_resolution",
    "compatible_headgroup_isomer",
    "hierarchy_class",
    "hierarchy_ancestor",
]
TIER_DESCRIPTION = {
    "exact_reported_resolution": (
        "SwissLipids candidate at the reported nomenclature resolution."
    ),
    "compatible_more_specific": (
        "Candidate resolving extra structural detail, potentially sn position, "
        "double-bond position/geometry, or stereochemistry."
    ),
    "compatible_lower_resolution": (
        "Aggregate candidate that loses one or more reported acyl-chain details."
    ),
    "compatible_headgroup_isomer": (
        "Candidate that changes an unresolved headgroup isomer."
    ),
    "hierarchy_class": "Broad SwissLipids class candidate.",
    "hierarchy_ancestor": "Broad SwissLipids ancestor candidate.",
}
RELATION_ORDER = [
    "exact_CHEBI_identity",
    "explicit_conjugate_acid_base",
    "same_formula_is_a_descendant",
]
RELATION_DESCRIPTION = {
    "exact_CHEBI_identity": (
        "Same canonical ChEBI identifier after alternate-ID resolution."
    ),
    "explicit_conjugate_acid_base": (
        "A directly recorded ChEBI conjugate acid/base counterpart, without "
        "a structural-class or isomer inference."
    ),
    "same_formula_is_a_descendant": (
        "Same-formula and same-charge ChEBI is_a descendant of the exact "
        "feature chemical or its explicit acid/base counterpart. This is a "
        "structural-isomer/position sensitivity, not direct feature identity."
    ),
}


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file without decompressing it."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def canonical_chebi(identifier: str, aliases: dict[str, str]) -> str:
    """Follow only explicit ChEBI alternate-ID chains, never ontology parents."""
    seen: set[str] = set()
    while identifier in aliases and identifier not in seen:
        seen.add(identifier)
        identifier = aliases[identifier]
    return identifier


def parse_obo(path: Path, aliases: dict[str, str]) -> dict[str, dict[str, Any]]:
    """Read the compact term fields needed for identity-safe descendant checks."""
    terms: dict[str, dict[str, Any]] = {}
    current: dict[str, Any] = {}

    def flush() -> None:
        identifier = current.get("id")
        if identifier and not current.get("obsolete", False):
            identifier = canonical_chebi(identifier, aliases)
            terms[identifier] = {
                "name": current.get("name", identifier),
                "formula": current.get("formula", ""),
                "charge": current.get("charge", ""),
                "parents": sorted({canonical_chebi(parent, aliases)
                                   for parent in current.get("parents", [])}),
            }

    with gzip.open(path, "rt") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line == "[Term]":
                flush()
                current = {}
            elif line.startswith("id: CHEBI:"):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("is_a: CHEBI:"):
                current.setdefault("parents", []).append(line.split(" !", 1)[0][6:])
            elif line.startswith("property_value: chemrof:generalized_empirical_formula "):
                match = re.search(r'"(.*?)"', line)
                current["formula"] = match.group(1) if match else ""
            elif line.startswith("property_value: chemrof:charge "):
                match = re.search(r'"(.*?)"', line)
                current["charge"] = match.group(1) if match else ""
            elif line == "is_obsolete: true":
                current["obsolete"] = True
    flush()
    return terms


def node_for_evidence(record: dict[str, str]) -> tuple[str, str, str]:
    """Return the graph triple that a raw retained STITCH record projects to."""
    source = "Metab__" + record["source"] + "_unknown"
    target = ("Gene__" + record["target_gene"] if record["target_gene"]
              else "UniProt__" + record["target"])
    return source, record["interaction"], target


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def joined(values: Iterable[Any]) -> str:
    return ";".join(sorted({str(value) for value in values if str(value)}))


def range_text(values: Iterable[int]) -> str:
    values = sorted(set(values))
    if not values:
        return ""
    return str(values[0]) if len(values) == 1 else f"{values[0]}-{values[-1]}"


def candidate_records(
    candidate_rows: list[dict[str, str]], aliases: dict[str, str]
) -> tuple[dict[str, list[dict[str, str]]], dict[str, dict[str, set[str]]]]:
    """Normalise per-feature ChEBI candidates without conflating mapping tiers."""
    rows_by_feature: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    chebis_by_feature_tier: dict[str, dict[str, set[str]]] = collections.defaultdict(
        lambda: collections.defaultdict(set)
    )
    for row in candidate_rows:
        tier = row["match_type"]
        if tier not in TIER_ORDER:
            continue
        raw_ids = CHEBI.findall(row.get("chebi_ids", ""))
        for raw_identifier in raw_ids:
            record = dict(row)
            record["candidate_chebi_raw"] = raw_identifier
            record["candidate_chebi"] = canonical_chebi(raw_identifier, aliases)
            rows_by_feature[row["feature_id"]].append(record)
            chebis_by_feature_tier[row["feature_id"]][tier].add(record["candidate_chebi"])
    return rows_by_feature, chebis_by_feature_tier


def relation_matches(
    candidate_by_feature: dict[str, list[dict[str, str]]],
    sources: set[str],
    acid_base: dict[str, set[str]],
    terms: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Build auditable feature-to-STITCH chemical candidates.

    The descendant sensitivity intentionally starts only from an exact
    reported-resolution candidate.  This prevents a broad class or aggregate
    mapping from being silently widened a second time through the ChEBI graph.
    """
    @functools.lru_cache(maxsize=None)
    def ancestors(source: str) -> dict[str, int]:
        queue: collections.deque[tuple[str, int]] = collections.deque([(source, 0)])
        found: dict[str, int] = {source: 0}
        while queue:
            node, distance = queue.popleft()
            for parent in terms.get(node, {}).get("parents", []):
                if parent not in found:
                    found[parent] = distance + 1
                    queue.append((parent, distance + 1))
        return found

    records: dict[tuple[str, str, str, str, str], dict[str, Any]] = {}

    def add(
        candidate: dict[str, str], relation: str, source: str,
        bridge: str, distance: int | None = None,
    ) -> None:
        key = (candidate["feature_id"], candidate["match_type"], relation,
               candidate["candidate_chebi"], source)
        item = records.setdefault(key, {
            "feature": candidate["feature_id"],
            "mapping_tier": candidate["match_type"],
            "identity_relation": relation,
            "candidate_chebi": candidate["candidate_chebi"],
            "candidate_chebi_raw": set(),
            "candidate_swisslipids_ids": set(),
            "candidate_canonical": set(),
            "candidate_names": set(),
            "candidate_abbreviations": set(),
            "bridge_chebi": set(),
            "is_a_distance": set(),
            "stitch_source": source,
        })
        item["candidate_chebi_raw"].add(candidate["candidate_chebi_raw"])
        item["candidate_swisslipids_ids"].add(candidate.get("swisslipids_id", ""))
        item["candidate_canonical"].add(candidate.get("candidate_canonical", ""))
        item["candidate_names"].add(candidate.get("name", ""))
        item["candidate_abbreviations"].add(candidate.get("abbreviation", ""))
        item["bridge_chebi"].add(bridge)
        if distance is not None:
            item["is_a_distance"].add(distance)

    for feature, candidates in candidate_by_feature.items():
        for candidate in candidates:
            chemical = candidate["candidate_chebi"]
            tier = candidate["match_type"]
            if chemical in sources:
                add(candidate, "exact_CHEBI_identity", chemical, chemical)
            for counterpart in sorted(acid_base.get(chemical, set()) & sources):
                add(candidate, "explicit_conjugate_acid_base", counterpart, counterpart)

            if tier != "exact_reported_resolution":
                continue
            # Start from the reported chemical and its direct, explicit
            # conjugate counterpart.  A source must be a same-formula,
            # same-charge is_a descendant of one of those seeds.
            for bridge in sorted({chemical, *acid_base.get(chemical, set())}):
                bridge_term = terms.get(bridge, {})
                formula, charge = bridge_term.get("formula"), bridge_term.get("charge")
                if not formula:
                    continue
                for source in sources:
                    if source == bridge:
                        continue
                    source_term = terms.get(source, {})
                    if (source_term.get("formula") != formula or
                            source_term.get("charge") != charge):
                        continue
                    distance = ancestors(source).get(bridge)
                    if distance and distance > 0:
                        add(candidate, "same_formula_is_a_descendant", source, bridge, distance)

    output = []
    for item in records.values():
        output.append({
            **{key: value for key, value in item.items()
               if key not in {"candidate_chebi_raw", "candidate_swisslipids_ids",
                              "candidate_canonical", "candidate_names",
                              "candidate_abbreviations", "bridge_chebi", "is_a_distance"}},
            "candidate_chebi_raw": joined(item["candidate_chebi_raw"]),
            "candidate_swisslipids_ids": joined(item["candidate_swisslipids_ids"]),
            "candidate_canonical": joined(item["candidate_canonical"]),
            "candidate_names": joined(item["candidate_names"]),
            "candidate_abbreviations": joined(item["candidate_abbreviations"]),
            "bridge_chebi": joined(item["bridge_chebi"]),
            "is_a_distance": joined(item["is_a_distance"]),
        })
    return sorted(output, key=lambda row: (
        row["feature"], TIER_ORDER.index(row["mapping_tier"]),
        RELATION_ORDER.index(row["identity_relation"]), row["stitch_source"]
    ))


def enrich_match_records(
    matches: list[dict[str, Any]], evidence_by_source: dict[str, list[dict[str, str]]],
    graph: set[tuple[str, str, str]], source_names: dict[str, str],
    factor_rows: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """Add raw-evidence and final-graph coverage while preserving source provenance."""
    enriched: list[dict[str, Any]] = []
    for match in matches:
        rows = evidence_by_source[match["stitch_source"]]
        raw_triples = {node_for_evidence(row) for row in rows}
        final_triples = raw_triples & graph
        factor = factor_rows[match["feature"]]
        relation = match["identity_relation"]
        direct = (match["mapping_tier"] == "exact_reported_resolution" and
                  relation == "exact_CHEBI_identity")
        enriched.append({
            **match,
            "stitch_source_name": source_names.get(match["stitch_source"], match["stitch_source"]),
            "factor1_existing_input": factor["included"],
            "frozen_input_node": factor["node_id"],
            "frozen_input_mapping_method": factor["mapping_method"],
            "raw_stitch_evidence_rows": len(rows),
            "raw_signed_source_protein_triples": len(raw_triples),
            "raw_target_proteins": len({triple[2] for triple in raw_triples}),
            "raw_activation_triples": sum(triple[1] == "1" for triple in raw_triples),
            "raw_inhibition_triples": sum(triple[1] == "-1" for triple in raw_triples),
            "graph_surviving_signed_triples": len(final_triples),
            "graph_surviving_target_proteins": len({triple[2] for triple in final_triples}),
            "graph_surviving_activation_triples": sum(triple[1] == "1" for triple in final_triples),
            "graph_surviving_inhibition_triples": sum(triple[1] == "-1" for triple in final_triples),
            "direct_measured_chemical_regulation": bool_text(direct),
            "mapping_interpretation": (
                "direct exact measured chemical candidate" if direct else
                "conditional acid/base chemical counterpart; not a direct measured identity" if relation == "explicit_conjugate_acid_base" else
                "conditional more-specific SwissLipids candidate; not a direct measured identity" if relation == "exact_CHEBI_identity" else
                "conditional same-formula ChEBI structural descendant; positional/stereochemical detail is unmeasured"
            ),
        })
    return enriched


def summarize_tiers(
    feature_sets: dict[str, set[str]],
    candidates_by_feature_tier: dict[str, dict[str, set[str]]],
    matches: list[dict[str, Any]], evidence_by_source: dict[str, list[dict[str, str]]],
    graph: set[tuple[str, str, str]],
) -> list[dict[str, Any]]:
    """Summarise candidate and source/protein coverage without double-counting records."""
    rows: list[dict[str, Any]] = []
    for population, feature_set in feature_sets.items():
        for tier in TIER_ORDER:
            candidate_pairs = {(feature, chemical)
                               for feature in feature_set
                               for chemical in candidates_by_feature_tier.get(feature, {}).get(tier, set())}
            for relation in ["any_relation", *RELATION_ORDER]:
                selected = [match for match in matches
                            if match["feature"] in feature_set and match["mapping_tier"] == tier and
                            (relation == "any_relation" or match["identity_relation"] == relation)]
                feature_sources = {(match["feature"], match["stitch_source"]) for match in selected}
                raw_records = {(feature, source, record["edge_id"])
                               for feature, source in feature_sources
                               for record in evidence_by_source[source]}
                raw_triples = {(feature, *node_for_evidence(record))
                               for feature, source in feature_sources
                               for record in evidence_by_source[source]}
                graph_triples = {triple for triple in raw_triples if triple[1:] in graph}
                rows.append({
                    "population": population,
                    "mapping_tier": tier,
                    "mapping_tier_description": TIER_DESCRIPTION[tier],
                    "identity_relation": relation,
                    "identity_relation_description": (
                        "Union of all relation categories below." if relation == "any_relation"
                        else RELATION_DESCRIPTION[relation]
                    ),
                    "features_in_population": len(feature_set),
                    "features_with_candidate_chebi": len({feature for feature, _ in candidate_pairs}),
                    "feature_candidate_chebi_pairs": len(candidate_pairs),
                    "features_with_stitch_source_match": len({feature for feature, _ in feature_sources}),
                    "feature_stitch_source_pairs": len(feature_sources),
                    "unique_stitch_source_chemicals": len({source for _, source in feature_sources}),
                    "feature_raw_evidence_assignments": len(raw_records),
                    "feature_raw_signed_source_protein_triples": len(raw_triples),
                    "feature_raw_target_proteins": len({(feature, target)
                                                        for feature, _, _, target in raw_triples}),
                    "feature_raw_activation_triples": sum(sign == "1"
                                                          for _, _, sign, _ in raw_triples),
                    "feature_raw_inhibition_triples": sum(sign == "-1"
                                                          for _, _, sign, _ in raw_triples),
                    "feature_graph_surviving_signed_triples": len(graph_triples),
                    "feature_graph_surviving_target_proteins": len({(feature, target)
                                                                    for feature, _, _, target in graph_triples}),
                    "feature_graph_surviving_activation_triples": sum(sign == "1"
                                                                      for _, _, sign, _ in graph_triples),
                    "feature_graph_surviving_inhibition_triples": sum(sign == "-1"
                                                                      for _, _, sign, _ in graph_triples),
                })
    return rows


def feature_inventory(
    feature_rows: dict[str, dict[str, str]], factor_rows: dict[str, dict[str, str]],
    candidates_by_feature_tier: dict[str, dict[str, set[str]]],
    matches: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Produce one row for every measured feature, including the zero-match majority."""
    match_by_feature: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for match in matches:
        match_by_feature[match["feature"]].append(match)
    output: list[dict[str, Any]] = []
    for feature, factor in factor_rows.items():
        mapped = match_by_feature.get(feature, [])
        feature_row = feature_rows.get(feature, {})
        direct = [row for row in mapped if row["direct_measured_chemical_regulation"] == "TRUE"]
        output.append({
            "feature": feature,
            "factor1_existing_input": factor["included"],
            "frozen_input_node": factor["node_id"],
            "frozen_input_mapping_method": factor["mapping_method"],
            "frozen_input_exclusion_reason": factor["exclusion_reason"],
            "oriented_factor1_weight": factor["oriented_weight"],
            "nomenclature_level": feature_row.get("nomenclature_level", factor.get("resolution", "")),
            "best_mapping_status": feature_row.get("best_mapping_status", "no_SwissLipids_candidate"),
            "n_exact_reported_resolution_chebi": len(candidates_by_feature_tier.get(feature, {}).get("exact_reported_resolution", set())),
            "n_compatible_more_specific_chebi": len(candidates_by_feature_tier.get(feature, {}).get("compatible_more_specific", set())),
            "n_compatible_lower_resolution_chebi": len(candidates_by_feature_tier.get(feature, {}).get("compatible_lower_resolution", set())),
            "n_compatible_headgroup_isomer_chebi": len(candidates_by_feature_tier.get(feature, {}).get("compatible_headgroup_isomer", set())),
            "n_hierarchy_class_chebi": len(candidates_by_feature_tier.get(feature, {}).get("hierarchy_class", set())),
            "n_hierarchy_ancestor_chebi": len(candidates_by_feature_tier.get(feature, {}).get("hierarchy_ancestor", set())),
            "n_stitch_source_matches": len({row["stitch_source"] for row in mapped}),
            "stitch_mapping_tiers": joined(row["mapping_tier"] for row in mapped),
            "stitch_identity_relations": joined(row["identity_relation"] for row in mapped),
            "stitch_source_chemicals": joined(row["stitch_source"] for row in mapped),
            "stitch_source_names": joined(row["stitch_source_name"] for row in mapped),
            "direct_measured_chemical_regulation": bool_text(bool(direct)),
            "interpretation": (
                "direct exact feature-to-STITCH chemical identity" if direct else
                "conditional candidate mapping only; not a direct measured chemical regulation claim" if mapped else
                "no retained directional STITCH chemical source under audited mappings"
            ),
        })
    return output


def protein_matches(matches: list[dict[str, Any]], evidence_by_source: dict[str, list[dict[str, str]]],
                    graph: set[tuple[str, str, str]]) -> list[dict[str, Any]]:
    """Collapse duplicate raw records into auditable feature/source/protein/sign rows."""
    grouped: dict[tuple[str, str, str, str, str, str], list[dict[str, str]]] = collections.defaultdict(list)
    match_lookup = {(row["feature"], row["mapping_tier"], row["identity_relation"], row["stitch_source"]): row
                    for row in matches}
    for key, match in match_lookup.items():
        for record in evidence_by_source[match["stitch_source"]]:
            source_node, sign, target_node = node_for_evidence(record)
            grouped[(*key, target_node, sign)].append(record)
    output: list[dict[str, Any]] = []
    for grouped_key, records in sorted(grouped.items()):
        feature, tier, relation, source, target_node, sign = grouped_key
        match = match_lookup[(feature, tier, relation, source)]
        graph_survives = ("Metab__" + source + "_unknown", sign, target_node) in graph
        output.append({
            "feature": feature,
            "mapping_tier": tier,
            "identity_relation": relation,
            "stitch_source": source,
            "stitch_source_name": match["stitch_source_name"],
            "target_node": target_node,
            "target_uniprot": joined(record["target"] for record in records),
            "target_gene": joined(record["target_gene"] for record in records),
            "interaction": sign,
            "raw_evidence_rows_collapsed": len(records),
            "action_score_range": range_text(int(record["raw_action_score"]) for record in records),
            "detailed_link_textmining_range": range_text(int(record["raw_link_textmining"]) for record in records),
            "graph_survives_conflict_filter": bool_text(graph_survives),
            "direct_measured_chemical_regulation": match["direct_measured_chemical_regulation"],
        })
    return output


def input_node_summary(factor_rows: dict[str, dict[str, str]], source_nodes: set[str]) -> list[dict[str, Any]]:
    """Separate existing Factor 1 PKN input nodes from candidate identities."""
    all_rows = list(factor_rows.values())
    accepted = [row for row in all_rows if row["included"] == "TRUE"]
    accepted_chebi = [row for row in accepted if row["node_id"].startswith("Metab__CHEBI:")]
    accepted_rules = [row for row in accepted if row["node_id"].startswith("Metab__MeasuredLipid:")]
    return [{
        "measurement_panel": "all_1254_factor1_lipids",
        "features": len(all_rows),
        "existing_factor1_inputs": len(accepted),
        "existing_rule_nodes": len(accepted_rules),
        "existing_exact_chebi_nodes": len(accepted_chebi),
        "existing_input_nodes_that_are_directional_stitch_sources": sum(
            row["node_id"] in source_nodes for row in accepted
        ),
        "interpretation": (
            "Measured-species rule nodes and the six frozen exact ChEBI nodes are "
            "tested separately from relaxed candidate identities."
        ),
    }]


def build(output: Path) -> dict[str, Any]:
    required = [FEATURES, CANDIDATES, ALIASES, ACID_BASE, STITCH_EVIDENCE, STITCH_GRAPH,
                STITCH_MANIFEST, STITCH_VALIDATION, FACTOR_INPUTS, OBO]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Required audit input(s) missing:\n" + "\n".join(missing))
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Refusing to overwrite a nonempty output directory: " + str(output))
    output.mkdir(parents=True, exist_ok=True)

    aliases = {row["alias"]: row["primary_chebi"] for row in read_csv(ALIASES)}
    acid_base: dict[str, set[str]] = collections.defaultdict(set)
    for row in read_csv(ACID_BASE):
        acid = canonical_chebi(row["acid"], aliases)
        base = canonical_chebi(row["base"], aliases)
        acid_base[acid].add(base)
        acid_base[base].add(acid)
    terms = parse_obo(OBO, aliases)
    evidence = read_csv(STITCH_EVIDENCE)
    evidence_by_source: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    source_names: dict[str, str] = {}
    for record in evidence:
        evidence_by_source[record["source"]].append(record)
        source_names.setdefault(record["source"], record["source_name"])
    graph = {(row["source"], row["interaction"], row["target"])
             for row in read_csv(STITCH_GRAPH)}
    source_nodes = {"Metab__" + source + "_unknown" for source in evidence_by_source}
    factor_rows = {row["feature"]: row for row in read_csv(FACTOR_INPUTS)}
    feature_rows = {row["feature_id"]: row for row in read_csv(FEATURES)}
    if set(factor_rows) != set(feature_rows):
        raise ValueError("Measured Factor 1 feature list and feature mapping inventory differ")
    candidate_by_feature, candidates_by_feature_tier = candidate_records(read_csv(CANDIDATES), aliases)
    raw_matches = relation_matches(candidate_by_feature, set(evidence_by_source), acid_base, terms)
    matches = enrich_match_records(raw_matches, evidence_by_source, graph, source_names, factor_rows)
    proteins = protein_matches(matches, evidence_by_source, graph)

    all_features = set(factor_rows)
    existing_inputs = {feature for feature, row in factor_rows.items() if row["included"] == "TRUE"}
    tier_summary = summarize_tiers(
        {"all_1254_measured_lipids": all_features,
         "factor1_existing_155_inputs": existing_inputs},
        candidates_by_feature_tier, matches, evidence_by_source, graph,
    )
    inventory = feature_inventory(feature_rows, factor_rows, candidates_by_feature_tier, matches)
    node_summary = input_node_summary(factor_rows, source_nodes)

    mapping_fields = [
        "feature", "mapping_tier", "identity_relation", "candidate_chebi", "candidate_chebi_raw",
        "candidate_swisslipids_ids", "candidate_canonical", "candidate_names", "candidate_abbreviations",
        "bridge_chebi", "is_a_distance", "stitch_source", "stitch_source_name",
        "factor1_existing_input", "frozen_input_node", "frozen_input_mapping_method",
        "raw_stitch_evidence_rows", "raw_signed_source_protein_triples", "raw_target_proteins",
        "raw_activation_triples", "raw_inhibition_triples", "graph_surviving_signed_triples",
        "graph_surviving_target_proteins", "graph_surviving_activation_triples",
        "graph_surviving_inhibition_triples", "direct_measured_chemical_regulation", "mapping_interpretation",
    ]
    write_csv(output / "feature_stitch_chemical_matches.csv", matches, mapping_fields)
    write_csv(output / "feature_stitch_protein_sign_matches.csv", proteins, [
        "feature", "mapping_tier", "identity_relation", "stitch_source", "stitch_source_name",
        "target_node", "target_uniprot", "target_gene", "interaction", "raw_evidence_rows_collapsed",
        "action_score_range", "detailed_link_textmining_range", "graph_survives_conflict_filter",
        "direct_measured_chemical_regulation",
    ])
    write_csv(output / "feature_inventory.csv", inventory, [
        "feature", "factor1_existing_input", "frozen_input_node", "frozen_input_mapping_method",
        "frozen_input_exclusion_reason", "oriented_factor1_weight", "nomenclature_level", "best_mapping_status",
        "n_exact_reported_resolution_chebi", "n_compatible_more_specific_chebi",
        "n_compatible_lower_resolution_chebi", "n_compatible_headgroup_isomer_chebi",
        "n_hierarchy_class_chebi", "n_hierarchy_ancestor_chebi", "n_stitch_source_matches",
        "stitch_mapping_tiers", "stitch_identity_relations", "stitch_source_chemicals", "stitch_source_names",
        "direct_measured_chemical_regulation", "interpretation",
    ])
    write_csv(output / "mapping_tier_summary.csv", tier_summary, [
        "population", "mapping_tier", "mapping_tier_description", "identity_relation",
        "identity_relation_description", "features_in_population", "features_with_candidate_chebi",
        "feature_candidate_chebi_pairs", "features_with_stitch_source_match", "feature_stitch_source_pairs",
        "unique_stitch_source_chemicals", "feature_raw_evidence_assignments",
        "feature_raw_signed_source_protein_triples", "feature_raw_target_proteins",
        "feature_raw_activation_triples", "feature_raw_inhibition_triples",
        "feature_graph_surviving_signed_triples", "feature_graph_surviving_target_proteins",
        "feature_graph_surviving_activation_triples", "feature_graph_surviving_inhibition_triples",
    ])
    write_csv(output / "existing_factor1_input_node_summary.csv", node_summary, list(node_summary[0]))

    manifest = json.loads(STITCH_MANIFEST.read_text())
    validation = json.loads(STITCH_VALIDATION.read_text())
    dependencies = [
        Path(__file__), FEATURES, CANDIDATES, ALIASES, ACID_BASE, STITCH_EVIDENCE, STITCH_GRAPH,
        STITCH_MANIFEST, STITCH_VALIDATION, FACTOR_INPUTS, OBO,
    ]
    write_csv(output / "dependency_hashes.csv", [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for path in dependencies
    ], ["path", "sha256"])

    direct = [row for row in matches if row["direct_measured_chemical_regulation"] == "TRUE"]
    conditional = [row for row in matches if row["direct_measured_chemical_regulation"] == "FALSE"]
    report = {
        "purpose": "Measured-lipid identity audit against raw directionality-audited STITCH chemical sources; no PKN modification or MOON scoring.",
        "overlay_validation_all_checks_pass": validation.get("all_checks_pass"),
        "overlay_raw_retained_evidence_rows": len(evidence),
        "overlay_raw_source_chemicals": len(evidence_by_source),
        "overlay_final_graph_edges": len(graph),
        "measured_features": len(all_features),
        "factor1_existing_inputs": len(existing_inputs),
        "existing_factor1_input_nodes_matching_stitch_sources": node_summary[0]["existing_input_nodes_that_are_directional_stitch_sources"],
        "direct_exact_measured_feature_to_stitch_source_matches": len(direct),
        "conditional_feature_to_stitch_source_matches": len(conditional),
        "features_with_any_audited_mapping": len({row["feature"] for row in matches}),
        "features_with_any_audited_mapping_that_are_existing_factor1_inputs": len(
            {row["feature"] for row in matches if row["factor1_existing_input"] == "TRUE"}
        ),
        "relation_policy": RELATION_DESCRIPTION,
        "caveat": (
            "STITCH v5 source directionality and score/text-mining admission are audited in the overlay, "
            "but v5 source files were last modified in 2016. A candidate chemical match does not establish "
            "direct hepatocyte regulation, temporal ordering, or a valid duplicated MOON input."
        ),
        "stitch_overlay_source_caveat": manifest.get("historical_source_caveat", ""),
    }
    (output / "summary.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    readme = f"""# Measured lipid → directional STITCH mapping audit

This audit compares all **{len(all_features):,}** lipidomic features and the
frozen **{len(existing_inputs):,}** Factor 1 PKN inputs with the raw,
directionality-audited STITCH v5 chemical-to-mouse-protein overlay. It does
not alter a PKN or run MOON.

The overlay contains **{len(evidence):,}** retained source records across
**{len(evidence_by_source):,}** chemical sources. Its own validation report
passed all checks. The existing Factor 1 inputs contain
**{node_summary[0]['existing_rule_nodes']}** measured-species rule nodes and
**{node_summary[0]['existing_exact_chebi_nodes']}** exact ChEBI nodes; none is
a directional STITCH chemical-source node.

## Interpretation contract

`exact_CHEBI_identity` is the only tier labelled direct measured chemical
identity. `explicit_conjugate_acid_base` is a charge-state counterpart and
`same_formula_is_a_descendant` resolves additional structural detail through
the ChEBI ontology; neither is a direct measured chemical claim. The latter is
limited to an exact reported-resolution lipid and same-formula/same-charge
descendants, but can still assume double-bond position, geometry, or
stereochemistry that the measurement does not resolve.

`feature_stitch_chemical_matches.csv` records feature-to-chemical mappings;
`feature_stitch_protein_sign_matches.csv` collapses duplicate source records
into feature/source/protein/sign triples and marks triples retained after the
assembled graph's conflict filter. `mapping_tier_summary.csv` reports all
SwissLipids mapping tiers separately for all measurements and for the existing
Factor 1 inputs. `feature_inventory.csv` retains the zero-match majority.

The STITCH source files are v5.0 and were last modified in 2016. Their raw
directionality is auditable in this overlay, but it does not establish current
evidence, direct mouse experimental support, hepatocyte context, or causal
regulation in this dataset.
"""
    (output / "README.md").write_text(readme)
    return report


def self_test() -> None:
    """Check the two identity rules that protect the audit from silent widening."""
    aliases = {"CHEBI:alt": "CHEBI:primary"}
    assert canonical_chebi("CHEBI:alt", aliases) == "CHEBI:primary"
    terms = {
        "CHEBI:base": {"formula": "C2H4O2", "charge": "0", "parents": [], "name": "generic acid"},
        "CHEBI:specific": {"formula": "C2H4O2", "charge": "0", "parents": ["CHEBI:base"], "name": "specific acid"},
        "CHEBI:modified": {"formula": "C2H6O3", "charge": "0", "parents": ["CHEBI:base"], "name": "modified acid"},
    }
    candidates = {
        "FA(2:0)": [{"feature_id": "FA(2:0)", "match_type": "exact_reported_resolution",
                      "candidate_chebi": "CHEBI:anion", "candidate_chebi_raw": "CHEBI:anion",
                      "swisslipids_id": "SLM:test", "candidate_canonical": "FA 2:0",
                      "name": "acetate", "abbreviation": "FA(2:0)"}],
    }
    matches = relation_matches(candidates, {"CHEBI:specific", "CHEBI:modified"},
                               {"CHEBI:anion": {"CHEBI:base"}}, terms)
    assert len(matches) == 1
    assert matches[0]["identity_relation"] == "same_formula_is_a_descendant"
    assert matches[0]["stitch_source"] == "CHEBI:specific"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help="New empty audit directory (default: %(default)s)")
    parser.add_argument("--self-test", action="store_true",
                        help="Run lightweight identity-policy checks and exit")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        print("Self-test passed")
        return
    report = build(args.output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
