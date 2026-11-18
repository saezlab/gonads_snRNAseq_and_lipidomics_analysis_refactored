#!/usr/bin/env python3
"""Create an auditable relaxed Factor 1 lipid-to-regulatory-PKN mapping ledger.

Each row in the default candidate panel copies one measured Factor 1 lipid
weight to one chemical source node with a signed chemical-to-protein edge that
survives the final directional exploratory-PKN conflict filter. It is a mapping
ledger only: it does not alter a PKN, resolve node collisions, or run MOON.

Default relaxed paths are deliberately limited to exact reported identity,
explicit ChEBI acid/base counterparts, same-formula/same-charge ChEBI
structural descendants, and SwissLipids-compatible more-specific candidates.
Broad class/ancestor, lower-resolution, and headgroup-isomer matches are
retained in an annotation-only table and are not silently input to MOON.

MRCLinksDB is excluded from candidates because its upstream cached sign was an
unconditional association; the nonpeptidic CellPhoneDB adapter is human-only
and unsigned. The directional STITCH v5 layer is included, but remains
exploratory because its raw source release is historical.

Run from repository root:
  python3 scripts/pkn/audit_factor1_relaxed_regulatory_mapping.py \
    --output results/causal_modelling/factor1_relaxed_lipid_regulatory_mapping
  python3 scripts/pkn/audit_factor1_relaxed_regulatory_mapping.py --self-test
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import tempfile
from pathlib import Path
from typing import Any, Iterable

# This sibling audit owns the independently tested ChEBI candidate parsing,
# explicit acid/base logic, and constrained same-formula is_a traversal.  It
# performs no I/O at import time.
import audit_measured_lipids_stitch_directional as stitch_audit


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn"
OVERLAY = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn_stitch_directional"
INPUTS = ROOT / "miscellaneous" / "older_results" / "causal_modelling" / "cosmos_factor1_pilot" / "inputs" / "factor1_lipid_node_mapping.csv"
OUTPUT_DEFAULT = ROOT / "results" / "causal_modelling" / "factor1_relaxed_lipid_regulatory_mapping"

FEATURES = BASE / "mapping_features.csv"
CANDIDATES = BASE / "mapping_candidates.csv"
ALIASES = BASE / "regulation_chebi_aliases.csv"
ACID_BASE = BASE / "regulation_chebi_acid_base_relations.csv"
OBO = ROOT / "support" / "cosmos_mouse" / "raw" / "regulation_chebi.obo.gz"
REGULATION = OVERLAY / "regulation_edges.csv"
GRAPH = OVERLAY / "pkn_exploratory_stitch_directional.csv"
NODES = OVERLAY / "nodes.csv"
STITCH_SUMMARY = OVERLAY / "stitch_directional_summary.json"
STITCH_VALIDATION = OVERLAY / "stitch_directional_validation.json"

DEFAULT_PATHS = {
    ("exact_reported_resolution", "exact_CHEBI_identity"): "strict_reported_identity",
    ("exact_reported_resolution", "explicit_conjugate_acid_base"): "relaxed_acid_base",
    ("exact_reported_resolution", "same_formula_is_a_descendant"):
        "relaxed_same_formula_is_a_descendant",
    ("compatible_more_specific", "exact_CHEBI_identity"):
        "relaxed_compatible_more_specific",
    ("compatible_more_specific", "explicit_conjugate_acid_base"):
        "relaxed_compatible_more_specific_acid_base",
}
POLICY_ORDER = [
    "strict_reported_identity",
    "relaxed_acid_base",
    "relaxed_same_formula_is_a_descendant",
    "relaxed_compatible_more_specific",
    "relaxed_compatible_more_specific_acid_base",
    "annotation_only_compatible_lower_resolution",
    "annotation_only_compatible_headgroup_isomer",
    "annotation_only_hierarchy_class",
    "annotation_only_hierarchy_ancestor",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def truth(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def joined(values: Iterable[Any]) -> str:
    return ";".join(sorted({str(value) for value in values if str(value)}))


def annotate_factor_weight_statistics(factor: dict[str, dict[str, str]]) -> None:
    """Add stable absolute-weight rank and inclusive percentile across all inputs."""
    ordered = sorted(
        factor,
        key=lambda feature: (-abs(float(factor[feature]["oriented_weight"])), feature),
    )
    denominator = max(len(ordered) - 1, 1)
    for position, feature in enumerate(ordered, start=1):
        factor[feature]["factor1_abs_weight"] = format(
            abs(float(factor[feature]["oriented_weight"])), ".17g"
        )
        factor[feature]["factor1_abs_rank_desc"] = str(position)
        factor[feature]["factor1_abs_percentile"] = format(
            100 * (position - 1) / denominator, ".12g"
        )


def policy_for(tier: str, relation: str) -> str:
    return DEFAULT_PATHS.get((tier, relation), "annotation_only_" + tier)


def policy_rank(tier: str, relation: str) -> int:
    policy = policy_for(tier, relation)
    return POLICY_ORDER.index(policy) if policy in POLICY_ORDER else len(POLICY_ORDER)


def rationale(tier: str, relation: str) -> str:
    policy = policy_for(tier, relation)
    if policy == "strict_reported_identity":
        return "Reported-resolution chemical identity equals a final signed PKN regulatory source."
    if policy == "relaxed_acid_base":
        return "Candidate source is a direct, explicit ChEBI conjugate acid/base counterpart."
    if policy == "relaxed_same_formula_is_a_descendant":
        return (
            "Candidate source is a same-formula, same-charge ChEBI descendant; its positional, "
            "geometric, or stereochemical detail is not measured."
        )
    if policy.startswith("relaxed_compatible_more_specific"):
        return (
            "SwissLipids proposes a more-specific chemical source; the measurement does not "
            "resolve all added structural detail."
        )
    return "Broad mapping annotation only; excluded from the default relaxed candidate panel."


def pkn_triple(record: dict[str, str]) -> tuple[str, str, str]:
    source = "Metab__" + record["source"] + "_unknown"
    target = "Gene__" + record["target_gene"] if record["target_gene"] else "UniProt__" + record["target"]
    return source, record["interaction"], target


def source_contexts(
    regulation: list[dict[str, str]], graph: set[tuple[str, str, str]],
    nodes: dict[str, dict[str, str]], terms: dict[str, dict[str, Any]],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Summarise only chemical-to-protein records that survived final conflicts."""
    raw_by_source: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    kept_by_source: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    for row in regulation:
        raw_by_source[row["source"]].append(row)
        if pkn_triple(row) in graph:
            kept_by_source[row["source"]].append(row)

    contexts: dict[str, dict[str, Any]] = {}
    inventory: list[dict[str, Any]] = []
    source_edges: list[dict[str, Any]] = []
    for source in sorted(raw_by_source):
        raw = raw_by_source[source]
        kept = kept_by_source.get(source, [])
        raw_triples = {pkn_triple(row) for row in raw}
        kept_triples = {pkn_triple(row) for row in kept}
        source_node = "Metab__" + source + "_unknown"
        source_name = next((row["source_name"] for row in kept if row["source_name"]), "")
        source_label = nodes.get(source_node, {}).get("label", "")
        if not source_name:
            source_name = terms.get(source, {}).get("name", source_label or source)
        context = {
            "pkn_source_chebi": source,
            "pkn_source_node": source_node,
            "pkn_source_label": source_label,
            "pkn_source_name": source_name,
            "raw_regulation_evidence_rows": len(raw),
            "raw_signed_source_target_triples": len(raw_triples),
            "surviving_regulation_evidence_rows": len(kept),
            "surviving_signed_source_target_triples": len(kept_triples),
            "surviving_target_proteins": len({target for _, _, target in kept_triples}),
            "surviving_activation_edges": sum(sign == "1" for _, sign, _ in kept_triples),
            "surviving_inhibition_edges": sum(sign == "-1" for _, sign, _ in kept_triples),
            "surviving_resources": joined(row["resource"] for row in kept),
            "surviving_evidence_tiers": joined(row["evidence_tier"] for row in kept),
            "surviving_mechanisms": joined(row["mechanism"] for row in kept),
            "surviving_source_versions": joined(row["source_version"] for row in kept),
            "edge_fingerprint_group": "edge_" + hashlib.sha256(
                "|".join(sorted(interaction + ":" + target for _, interaction, target in kept_triples)).encode()
            ).hexdigest()[:16] if kept_triples else "",
            "source_node_exists_in_catalogue": truth(source_node in nodes),
            "has_surviving_signed_regulatory_edge": truth(bool(kept_triples)),
        }
        inventory.append(context)
        if not kept_triples:
            continue
        contexts[source] = context
        by_triple: dict[tuple[str, str, str], list[dict[str, str]]] = collections.defaultdict(list)
        for row in kept:
            by_triple[pkn_triple(row)].append(row)
        for (source_node, interaction, target_node), records in sorted(by_triple.items()):
            source_edges.append({
                "pkn_source_chebi": source,
                "pkn_source_node": source_node,
                "pkn_source_name": source_name,
                "target_node": target_node,
                "interaction": interaction,
                "target_uniprot": joined(row["target"] for row in records),
                "target_gene": joined(row["target_gene"] for row in records),
                "surviving_evidence_rows_collapsed": len(records),
                "resources": joined(row["resource"] for row in records),
                "evidence_tiers": joined(row["evidence_tier"] for row in records),
                "mechanisms": joined(row["mechanism"] for row in records),
                "source_record_ids": joined(row["source_record_id"] for row in records),
                "references": joined(row["references"] for row in records),
                "source_versions": joined(row["source_version"] for row in records),
                "sign_bases": joined(row["sign_basis"] for row in records),
                "quality_notes": joined(row["quality_notes"] for row in records),
            })
    return contexts, inventory, source_edges


def candidate_paths(
    candidates: dict[str, list[dict[str, str]]], sources: set[str],
    acid_base: dict[str, set[str]], terms: dict[str, dict[str, Any]],
) -> list[dict[str, Any]]:
    """Reuse the reviewed identity/acid-base/structural traversal, then re-label it."""
    paths = stitch_audit.relation_matches(candidates, sources, acid_base, terms)
    output: list[dict[str, Any]] = []
    for row in paths:
        tier, relation = row["mapping_tier"], row["identity_relation"]
        output.append({
            "feature_id": row["feature"],
            "mapping_tier": tier,
            "identity_relation": relation,
            "candidate_chebi": row["candidate_chebi"],
            "candidate_chebi_raw": row["candidate_chebi_raw"],
            "candidate_swisslipids_ids": row["candidate_swisslipids_ids"],
            "candidate_canonical": row["candidate_canonical"],
            "candidate_names": row["candidate_names"],
            "candidate_abbreviations": row["candidate_abbreviations"],
            "bridge_chebi": row["bridge_chebi"],
            "is_a_distance": row["is_a_distance"],
            "pkn_source_chebi": row["stitch_source"],
            "path_policy": policy_for(tier, relation),
            "path_is_default_relaxed": truth((tier, relation) in DEFAULT_PATHS),
            "path_rationale": rationale(tier, relation),
        })
    return sorted(output, key=lambda row: (
        row["feature_id"], policy_rank(row["mapping_tier"], row["identity_relation"]),
        row["pkn_source_chebi"], row["candidate_chebi"],
    ))


def collapse_assignments(
    paths: list[dict[str, Any]], contexts: dict[str, dict[str, Any]],
    features: dict[str, dict[str, str]], factor: dict[str, dict[str, str]],
) -> list[dict[str, Any]]:
    """Make one assignment per measured feature and candidate PKN source node."""
    grouped: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
    for row in paths:
        grouped[(row["feature_id"], row["pkn_source_chebi"])].append(row)
    output: list[dict[str, Any]] = []
    for (feature, source), rows in sorted(grouped.items()):
        rows.sort(key=lambda row: (
            policy_rank(row["mapping_tier"], row["identity_relation"]), row["candidate_chebi"],
        ))
        default_paths = [row for row in rows if row["path_is_default_relaxed"] == "TRUE"]
        primary = default_paths[0] if default_paths else rows[0]
        context, meta, measurement = contexts[source], features[feature], factor[feature]
        assignment_id = "map_" + hashlib.sha256((feature + "|" + source).encode()).hexdigest()[:16]
        output.append({
            "assignment_id": assignment_id,
            "measurement_feature": feature,
            "candidate_measurement_value": measurement["oriented_weight"],
            "raw_factor1_weight": measurement["raw_weight"],
            "factor1_abs_weight": measurement["factor1_abs_weight"],
            "factor1_abs_rank_desc": measurement["factor1_abs_rank_desc"],
            "factor1_abs_percentile": measurement["factor1_abs_percentile"],
            "factor1_existing_frozen_input": measurement["included"],
            "frozen_input_node": measurement["node_id"],
            "frozen_input_mapping_method": measurement["mapping_method"],
            "frozen_input_exclusion_reason": measurement["exclusion_reason"],
            "measurement_canonical_name": meta.get("canonical_name", ""),
            "measurement_lipid_class": meta.get("class", ""),
            "measurement_nomenclature_level": meta.get("nomenclature_level", ""),
            "measurement_best_mapping_status": meta.get("best_mapping_status", ""),
            "pkn_source_chebi": source,
            "pkn_source_node": context["pkn_source_node"],
            "candidate_pkn_node": context["pkn_source_node"],
            "pkn_source_label": context["pkn_source_label"],
            "pkn_source_name": context["pkn_source_name"],
            "primary_mapping_policy": primary["path_policy"],
            "primary_mapping_tier": primary["mapping_tier"],
            "primary_identity_relation": primary["identity_relation"],
            "mapping_relation": primary["identity_relation"],
            "primary_mapping_rationale": primary["path_rationale"],
            "all_mapping_policies": joined(row["path_policy"] for row in rows),
            "all_mapping_tiers": joined(row["mapping_tier"] for row in rows),
            "all_identity_relations": joined(row["identity_relation"] for row in rows),
            "candidate_chebi_ids": joined(row["candidate_chebi"] for row in rows),
            "candidate_chebi_raw_ids": joined(row["candidate_chebi_raw"] for row in rows),
            "candidate_swisslipids_ids": joined(row["candidate_swisslipids_ids"] for row in rows),
            "candidate_canonical_names": joined(row["candidate_canonical"] for row in rows),
            "candidate_names": joined(row["candidate_names"] for row in rows),
            "candidate_abbreviations": joined(row["candidate_abbreviations"] for row in rows),
            "bridge_chebi_ids": joined(row["bridge_chebi"] for row in rows),
            "is_a_distance": joined(row["is_a_distance"] for row in rows),
            "n_mapping_paths_collapsed": len(rows),
            "strict_reported_identity": truth(any(
                row["path_policy"] == "strict_reported_identity" for row in rows
            )),
            "default_relaxed_candidate": truth(bool(default_paths)),
            "assignment_role": "default_relaxed_candidate" if default_paths else "annotation_only",
            "surviving_resources": context["surviving_resources"],
            "surviving_evidence_tiers": context["surviving_evidence_tiers"],
            "surviving_mechanisms": context["surviving_mechanisms"],
            "surviving_source_versions": context["surviving_source_versions"],
            "edge_fingerprint_group": context["edge_fingerprint_group"],
            "surviving_regulation_evidence_rows": context["surviving_regulation_evidence_rows"],
            "surviving_signed_source_target_triples": context["surviving_signed_source_target_triples"],
            "surviving_target_proteins": context["surviving_target_proteins"],
            "surviving_activation_edges": context["surviving_activation_edges"],
            "surviving_inhibition_edges": context["surviving_inhibition_edges"],
        })
    return output


def annotate_collisions(assignments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Flag intended feature expansion separately from source-node collisions."""
    all_by_feature: dict[str, set[str]] = collections.defaultdict(set)
    default_by_feature: dict[str, set[str]] = collections.defaultdict(set)
    default_by_source: dict[str, set[str]] = collections.defaultdict(set)
    for row in assignments:
        all_by_feature[row["measurement_feature"]].add(row["pkn_source_node"])
        if row["default_relaxed_candidate"] == "TRUE":
            default_by_feature[row["measurement_feature"]].add(row["pkn_source_node"])
            default_by_source[row["pkn_source_node"]].add(row["measurement_feature"])
    for row in assignments:
        feature, source = row["measurement_feature"], row["pkn_source_node"]
        default_nodes, source_features = default_by_feature[feature], default_by_source[source]
        row["n_candidate_source_nodes_for_measurement_all"] = len(all_by_feature[feature])
        row["n_default_candidate_source_nodes_for_measurement"] = len(default_nodes)
        row["measurement_value_is_copied_to_multiple_default_nodes"] = truth(len(default_nodes) > 1)
        row["default_candidate_source_nodes_for_measurement"] = joined(default_nodes)
        row["n_measurements_assigned_to_source_node_default"] = len(source_features)
        row["source_node_receives_multiple_measurement_values_default"] = truth(len(source_features) > 1)
        row["other_measurements_assigned_to_source_node_default"] = joined(source_features - {feature})
        row["requires_pre_scoring_collision_resolution"] = truth(
            row["default_relaxed_candidate"] == "TRUE" and len(source_features) > 1
        )
    collisions = []
    for source, features in sorted(default_by_source.items()):
        if len(features) > 1:
            related = [row for row in assignments if row["pkn_source_node"] == source]
            collisions.append({
                "pkn_source_node": source,
                "pkn_source_chebi": related[0]["pkn_source_chebi"],
                "pkn_source_name": related[0]["pkn_source_name"],
                "n_measurement_features": len(features),
                "measurement_features": joined(features),
                "candidate_measurement_values": joined(
                    row["candidate_measurement_value"] for row in related
                    if row["default_relaxed_candidate"] == "TRUE"
                ),
                "resolution_required": "TRUE",
                "required_decision": (
                    "Score this source node once as a unique topology seed; preserve both measurement assignments "
                    "for assignment-specific sign-consistency and path-extraction review."
                ),
            })
    return collisions


def unique_topology_seeds(assignments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Provide one topology seed per source node without inventing collision values."""
    grouped: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in assignments:
        if row["default_relaxed_candidate"] == "TRUE":
            grouped[row["candidate_pkn_node"]].append(row)
    output = []
    for node, rows in sorted(grouped.items()):
        rows.sort(key=lambda row: row["measurement_feature"])
        values = {row["candidate_measurement_value"] for row in rows}
        collision = len({row["measurement_feature"] for row in rows}) > 1
        output.append({
            "candidate_pkn_node": node,
            "pkn_source_chebi": rows[0]["pkn_source_chebi"],
            "pkn_source_name": rows[0]["pkn_source_name"],
            "edge_fingerprint_group": rows[0]["edge_fingerprint_group"],
            "unique_topology_seed": "TRUE",
            "n_measurement_assignments": len(rows),
            "measurement_features": joined(row["measurement_feature"] for row in rows),
            "assignment_ids": joined(row["assignment_id"] for row in rows),
            "candidate_measurement_values": joined(values),
            "n_distinct_measurement_values": len(values),
            "source_node_collision": truth(collision),
            "unique_seed_value": next(iter(values)) if len(values) == 1 else "",
            "seed_value_status": (
                "single_assignment_value" if len(values) == 1 else
                "multiple_assignment_values_preserved_for_assignment_specific_review"
            ),
            "surviving_resources": joined(row["surviving_resources"] for row in rows),
            "surviving_signed_source_target_triples": rows[0]["surviving_signed_source_target_triples"],
        })
    return output


def assignment_edges(assignments: list[dict[str, Any]], source_edges: list[dict[str, Any]]) -> list[dict[str, Any]]:
    edges_by_source: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in source_edges:
        edges_by_source[row["pkn_source_chebi"]].append(row)
    output = []
    for assignment in assignments:
        for edge in edges_by_source[assignment["pkn_source_chebi"]]:
            output.append({
                "assignment_id": assignment["assignment_id"],
                "measurement_feature": assignment["measurement_feature"],
                "candidate_measurement_value": assignment["candidate_measurement_value"],
                "default_relaxed_candidate": assignment["default_relaxed_candidate"],
                "primary_mapping_policy": assignment["primary_mapping_policy"],
                **edge,
            })
    return sorted(output, key=lambda row: (
        row["measurement_feature"], row["pkn_source_chebi"], row["target_node"], row["interaction"],
    ))


def feature_summary(
    factor: dict[str, dict[str, str]], features: dict[str, dict[str, str]],
    candidate_ids: dict[str, dict[str, set[str]]], assignments: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_feature: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in assignments:
        by_feature[row["measurement_feature"]].append(row)
    output = []
    for feature in sorted(factor):
        rows = by_feature.get(feature, [])
        default = [row for row in rows if row["default_relaxed_candidate"] == "TRUE"]
        annotation = [row for row in rows if row["default_relaxed_candidate"] != "TRUE"]
        if default:
            interpretation = "One or more default relaxed final regulatory source candidates."
        elif annotation:
            interpretation = "Only broad annotation-only regulatory source candidates."
        elif not any(candidate_ids.get(feature, {}).values()):
            interpretation = "No ChEBI-addressable SwissLipids candidate at a retained mapping tier."
        else:
            interpretation = "No final signed regulatory PKN source under the retained mapping rules."
        meta, measurement = features[feature], factor[feature]
        output.append({
            "measurement_feature": feature,
            "candidate_measurement_value": measurement["oriented_weight"],
            "raw_factor1_weight": measurement["raw_weight"],
            "factor1_abs_weight": measurement["factor1_abs_weight"],
            "factor1_abs_rank_desc": measurement["factor1_abs_rank_desc"],
            "factor1_abs_percentile": measurement["factor1_abs_percentile"],
            "factor1_existing_frozen_input": measurement["included"],
            "frozen_input_node": measurement["node_id"],
            "frozen_input_exclusion_reason": measurement["exclusion_reason"],
            "measurement_canonical_name": meta.get("canonical_name", ""),
            "measurement_lipid_class": meta.get("class", ""),
            "measurement_nomenclature_level": meta.get("nomenclature_level", ""),
            "measurement_best_mapping_status": meta.get("best_mapping_status", ""),
            "n_exact_reported_resolution_chebi": len(candidate_ids.get(feature, {}).get("exact_reported_resolution", set())),
            "n_compatible_more_specific_chebi": len(candidate_ids.get(feature, {}).get("compatible_more_specific", set())),
            "n_compatible_lower_resolution_chebi": len(candidate_ids.get(feature, {}).get("compatible_lower_resolution", set())),
            "n_compatible_headgroup_isomer_chebi": len(candidate_ids.get(feature, {}).get("compatible_headgroup_isomer", set())),
            "n_hierarchy_class_chebi": len(candidate_ids.get(feature, {}).get("hierarchy_class", set())),
            "n_hierarchy_ancestor_chebi": len(candidate_ids.get(feature, {}).get("hierarchy_ancestor", set())),
            "n_final_regulatory_source_nodes_all": len({row["pkn_source_node"] for row in rows}),
            "n_default_relaxed_candidate_source_nodes": len({row["pkn_source_node"] for row in default}),
            "n_annotation_only_source_nodes": len({row["pkn_source_node"] for row in annotation}),
            "default_relaxed_candidate_source_nodes": joined(row["pkn_source_node"] for row in default),
            "default_relaxed_candidate_source_chemicals": joined(row["pkn_source_chebi"] for row in default),
            "default_relaxed_mapping_policies": joined(row["primary_mapping_policy"] for row in default),
            "all_candidate_source_nodes": joined(row["pkn_source_node"] for row in rows),
            "measurement_value_is_copied_to_multiple_default_nodes": (
                default[0]["measurement_value_is_copied_to_multiple_default_nodes"] if default else "FALSE"
            ),
            "has_cross_feature_source_node_collision_default": truth(any(
                row["source_node_receives_multiple_measurement_values_default"] == "TRUE" for row in default
            )),
            "interpretation": interpretation,
        })
    return output


def policy_summary(assignments: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for policy in POLICY_ORDER + ["all_candidate_assignments"]:
        selected = assignments if policy == "all_candidate_assignments" else [
            row for row in assignments if row["primary_mapping_policy"] == policy
        ]
        default = [row for row in selected if row["default_relaxed_candidate"] == "TRUE"]
        rows.append({
            "mapping_policy": policy,
            "default_relaxed_policy": truth(policy in DEFAULT_PATHS.values()),
            "candidate_assignments": len(selected),
            "measurement_features": len({row["measurement_feature"] for row in selected}),
            "unique_pkn_source_nodes": len({row["pkn_source_node"] for row in selected}),
            "default_relaxed_assignments": len(default),
            "default_relaxed_measurement_features": len({row["measurement_feature"] for row in default}),
            "feature_to_many_node_copies": len({
                row["measurement_feature"] for row in default
                if row["measurement_value_is_copied_to_multiple_default_nodes"] == "TRUE"
            }),
            "cross_feature_source_node_collisions": len({
                row["pkn_source_node"] for row in default
                if row["source_node_receives_multiple_measurement_values_default"] == "TRUE"
            }),
            "resources": joined(resource for row in selected for resource in row["surviving_resources"].split(";")),
        })
    return rows


def resource_summary(
    regulation: list[dict[str, str]], graph: set[tuple[str, str, str]], assignments: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows = []
    notes = {
        "GuideToPharmacology": "Current mouse curated chemical-protein layer; secondary chemical mappings remain annotated.",
        "BRENDA": "Cached signed activity-modulation fallback; original source chemistry and per-edge provenance are incomplete.",
        "STITCH_directional_v5": "Raw chemical-actor directional layer; source release is historical (2016) and exploratory.",
    }
    for resource in sorted({row["resource"] for row in regulation}):
        raw = [row for row in regulation if row["resource"] == resource]
        kept = [row for row in raw if pkn_triple(row) in graph]
        candidates = [row for row in assignments if resource in row["surviving_resources"].split(";")]
        default = [row for row in candidates if row["default_relaxed_candidate"] == "TRUE"]
        rows.append({
            "resource": resource,
            "resource_status": "included_final_directional_exploratory_pkn",
            "raw_regulation_evidence_rows": len(raw),
            "surviving_regulation_evidence_rows": len(kept),
            "surviving_source_chemicals": len({row["source"] for row in kept}),
            "surviving_signed_source_target_triples": len({pkn_triple(row) for row in kept}),
            "candidate_assignments_all": len(candidates),
            "candidate_measurement_features_all": len({row["measurement_feature"] for row in candidates}),
            "default_relaxed_candidate_assignments": len(default),
            "default_relaxed_measurement_features": len({row["measurement_feature"] for row in default}),
            "policy_note": notes.get(resource, "Included only when signed evidence survived final conflict filtering."),
        })
    rows += [
        {
            "resource": "MRCLinksDB", "resource_status": "excluded_unsigned_annotation",
            "raw_regulation_evidence_rows": 0, "surviving_regulation_evidence_rows": 0,
            "surviving_source_chemicals": 0, "surviving_signed_source_target_triples": 0,
            "candidate_assignments_all": 0, "candidate_measurement_features_all": 0,
            "default_relaxed_candidate_assignments": 0, "default_relaxed_measurement_features": 0,
            "policy_note": "Excluded: cached sign is unconditional, including transport; it is not signed regulatory evidence.",
        },
        {
            "resource": "CellPhoneDB_nonpeptidic_adapter", "resource_status": "excluded_human_only_unsigned_annotation",
            "raw_regulation_evidence_rows": 0, "surviving_regulation_evidence_rows": 0,
            "surviving_source_chemicals": 0, "surviving_signed_source_target_triples": 0,
            "candidate_assignments_all": 0, "candidate_measurement_features_all": 0,
            "default_relaxed_candidate_assignments": 0, "default_relaxed_measurement_features": 0,
            "policy_note": "Excluded: human-only adapter assigns a positive sign without a checked effect; it is not projected to mouse.",
        },
        {
            "resource": "STITCH_cached_service_layer", "resource_status": "replaced_not_used",
            "raw_regulation_evidence_rows": 0, "surviving_regulation_evidence_rows": 0,
            "surviving_source_chemicals": 0, "surviving_signed_source_target_triples": 0,
            "candidate_assignments_all": 0, "candidate_measurement_features_all": 0,
            "default_relaxed_candidate_assignments": 0, "default_relaxed_measurement_features": 0,
            "policy_note": "Not used: orientation-unknown cached STITCH rows were replaced by raw chemical-actor directional STITCH v5 evidence.",
        },
    ]
    return rows


def write_readme(output: Path) -> None:
    (output / "README.md").write_text("""# Relaxed Factor 1 lipid-to-regulatory-PKN mapping
This pre-scoring ledger assigns each Factor 1 lipidomic measurement to all
supported candidate **chemical source nodes** in the directional exploratory
PKN. It does not change the PKN, choose a mapping draw, aggregate measurements,
or run MOON.
`default_relaxed_measurement_assignments.csv` contains only exact identity,
direct acid/base, same-formula structural-descendant, and compatible-more-
specific mappings. Every listed source has at least one signed direct
chemical-to-protein edge surviving the final PKN conflict filter.
`all_candidate_regulatory_assignments.csv` also contains broad class/ancestor,
lower-resolution, and headgroup-isomer annotations. Those rows have
`assignment_role=annotation_only` and must not be added to an upstream MOON
vector unless a separate policy is explicitly selected.
`mapping_paths.csv` preserves every SwissLipids/ChEBI route collapsed into an
assignment. `candidate_assignment_regulatory_edges.csv` provides direct
source-to-protein evidence and resource provenance. `measurement_mapping_summary.csv`
covers all 1,254 Factor 1 measurements. `candidate_node_collisions.csv` lists
the rare cases where different measurements reach the same PKN source node.
`default_relaxed_unique_topology_seeds.csv` has one row per default source node;
its value is blank for a cross-feature collision rather than inventing an
aggregation.
`candidate_measurement_value` is the unchanged oriented Factor 1 lipid weight.
Copying it to several candidate source nodes creates candidate signalling
hypotheses; it does not establish ligand potency, receptor occupancy, flux, or
causality. Feature-to-many-node copies are intentional here. A source node that
receives values from multiple *different* measurements should be scored once as
a unique topology seed; retain its collision rows for assignment-specific
sign-consistency and path-extraction review rather than counting it twice.
MRCLinksDB and the nonpeptidic CellPhoneDB adapter are reported in
`resource_policy_summary.csv` as excluded annotations. STITCH v5 is directionally
auditable but historical, and remains exploratory.
""")


def build(output: Path) -> dict[str, Any]:
    required = [FEATURES, CANDIDATES, ALIASES, ACID_BASE, OBO, INPUTS, REGULATION, GRAPH, NODES, STITCH_SUMMARY, STITCH_VALIDATION]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Required input(s) missing:\n" + "\n".join(missing))
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Refusing to overwrite nonempty output directory: " + str(output))
    output.mkdir(parents=True, exist_ok=True)

    aliases = {row["alias"]: row["primary_chebi"] for row in read_csv(ALIASES)}
    acid_base: dict[str, set[str]] = collections.defaultdict(set)
    for row in read_csv(ACID_BASE):
        acid = stitch_audit.canonical_chebi(row["acid"], aliases)
        base = stitch_audit.canonical_chebi(row["base"], aliases)
        acid_base[acid].add(base)
        acid_base[base].add(acid)
    terms = stitch_audit.parse_obo(OBO, aliases)
    regulation = read_csv(REGULATION)
    graph = {(row["source"], row["interaction"], row["target"]) for row in read_csv(GRAPH)}
    nodes = {row["node_id"]: row for row in read_csv(NODES)}
    contexts, source_inventory, source_edges = source_contexts(regulation, graph, nodes, terms)

    factor = {row["feature"]: row for row in read_csv(INPUTS)}
    annotate_factor_weight_statistics(factor)
    features = {row["feature_id"]: row for row in read_csv(FEATURES)}
    if set(factor) != set(features):
        raise ValueError("Factor 1 measurements differ from the frozen feature-mapping inventory")
    candidates, candidate_ids = stitch_audit.candidate_records(read_csv(CANDIDATES), aliases)
    paths = candidate_paths(candidates, set(contexts), acid_base, terms)
    assignments = collapse_assignments(paths, contexts, features, factor)
    collisions = annotate_collisions(assignments)
    seeds = unique_topology_seeds(assignments)
    edges = assignment_edges(assignments, source_edges)
    summaries = feature_summary(factor, features, candidate_ids, assignments)
    policies = policy_summary(assignments)
    resources = resource_summary(regulation, graph, assignments)

    assignment_fields = [
        "assignment_id", "measurement_feature", "candidate_measurement_value", "raw_factor1_weight", "factor1_abs_weight", "factor1_abs_rank_desc", "factor1_abs_percentile",
        "factor1_existing_frozen_input", "frozen_input_node", "frozen_input_mapping_method", "frozen_input_exclusion_reason",
        "measurement_canonical_name", "measurement_lipid_class", "measurement_nomenclature_level", "measurement_best_mapping_status",
        "pkn_source_chebi", "pkn_source_node", "candidate_pkn_node", "pkn_source_label", "pkn_source_name", "primary_mapping_policy",
        "primary_mapping_tier", "primary_identity_relation", "mapping_relation", "primary_mapping_rationale", "all_mapping_policies",
        "all_mapping_tiers", "all_identity_relations", "candidate_chebi_ids", "candidate_chebi_raw_ids",
        "candidate_swisslipids_ids", "candidate_canonical_names", "candidate_names", "candidate_abbreviations",
        "bridge_chebi_ids", "is_a_distance", "n_mapping_paths_collapsed", "strict_reported_identity",
        "default_relaxed_candidate", "assignment_role", "surviving_resources", "surviving_evidence_tiers",
        "surviving_mechanisms", "surviving_source_versions", "edge_fingerprint_group", "surviving_regulation_evidence_rows",
        "surviving_signed_source_target_triples", "surviving_target_proteins", "surviving_activation_edges",
        "surviving_inhibition_edges", "n_candidate_source_nodes_for_measurement_all",
        "n_default_candidate_source_nodes_for_measurement", "measurement_value_is_copied_to_multiple_default_nodes",
        "default_candidate_source_nodes_for_measurement", "n_measurements_assigned_to_source_node_default",
        "source_node_receives_multiple_measurement_values_default", "other_measurements_assigned_to_source_node_default",
        "requires_pre_scoring_collision_resolution",
    ]
    path_fields = [
        "feature_id", "mapping_tier", "identity_relation", "candidate_chebi", "candidate_chebi_raw",
        "candidate_swisslipids_ids", "candidate_canonical", "candidate_names", "candidate_abbreviations",
        "bridge_chebi", "is_a_distance", "pkn_source_chebi", "path_policy", "path_is_default_relaxed", "path_rationale",
    ]
    source_fields = [
        "pkn_source_chebi", "pkn_source_node", "pkn_source_label", "pkn_source_name", "raw_regulation_evidence_rows",
        "raw_signed_source_target_triples", "surviving_regulation_evidence_rows", "surviving_signed_source_target_triples",
        "surviving_target_proteins", "surviving_activation_edges", "surviving_inhibition_edges", "surviving_resources",
        "surviving_evidence_tiers", "surviving_mechanisms", "surviving_source_versions", "edge_fingerprint_group", "source_node_exists_in_catalogue",
        "has_surviving_signed_regulatory_edge",
    ]
    source_edge_fields = [
        "pkn_source_chebi", "pkn_source_node", "pkn_source_name", "target_node", "interaction", "target_uniprot",
        "target_gene", "surviving_evidence_rows_collapsed", "resources", "evidence_tiers", "mechanisms",
        "source_record_ids", "references", "source_versions", "sign_bases", "quality_notes",
    ]
    assignment_edge_fields = ["assignment_id", "measurement_feature", "candidate_measurement_value", "default_relaxed_candidate", "primary_mapping_policy"] + source_edge_fields
    feature_fields = [
        "measurement_feature", "candidate_measurement_value", "raw_factor1_weight", "factor1_abs_weight", "factor1_abs_rank_desc", "factor1_abs_percentile", "factor1_existing_frozen_input",
        "frozen_input_node", "frozen_input_exclusion_reason", "measurement_canonical_name", "measurement_lipid_class",
        "measurement_nomenclature_level", "measurement_best_mapping_status", "n_exact_reported_resolution_chebi",
        "n_compatible_more_specific_chebi", "n_compatible_lower_resolution_chebi", "n_compatible_headgroup_isomer_chebi",
        "n_hierarchy_class_chebi", "n_hierarchy_ancestor_chebi", "n_final_regulatory_source_nodes_all",
        "n_default_relaxed_candidate_source_nodes", "n_annotation_only_source_nodes",
        "default_relaxed_candidate_source_nodes", "default_relaxed_candidate_source_chemicals",
        "default_relaxed_mapping_policies", "all_candidate_source_nodes",
        "measurement_value_is_copied_to_multiple_default_nodes", "has_cross_feature_source_node_collision_default", "interpretation",
    ]
    policy_fields = [
        "mapping_policy", "default_relaxed_policy", "candidate_assignments", "measurement_features",
        "unique_pkn_source_nodes", "default_relaxed_assignments", "default_relaxed_measurement_features",
        "feature_to_many_node_copies", "cross_feature_source_node_collisions", "resources",
    ]
    resource_fields = [
        "resource", "resource_status", "raw_regulation_evidence_rows", "surviving_regulation_evidence_rows",
        "surviving_source_chemicals", "surviving_signed_source_target_triples", "candidate_assignments_all",
        "candidate_measurement_features_all", "default_relaxed_candidate_assignments",
        "default_relaxed_measurement_features", "policy_note",
    ]
    collision_fields = [
        "pkn_source_node", "pkn_source_chebi", "pkn_source_name", "n_measurement_features", "measurement_features",
        "candidate_measurement_values", "resolution_required", "required_decision",
    ]
    seed_fields = [
        "candidate_pkn_node", "pkn_source_chebi", "pkn_source_name", "edge_fingerprint_group", "unique_topology_seed",
        "n_measurement_assignments", "measurement_features", "assignment_ids", "candidate_measurement_values",
        "n_distinct_measurement_values", "source_node_collision", "unique_seed_value", "seed_value_status",
        "surviving_resources", "surviving_signed_source_target_triples",
    ]
    write_csv(output / "mapping_paths.csv", paths, path_fields)
    write_csv(output / "all_candidate_regulatory_assignments.csv", assignments, assignment_fields)
    write_csv(output / "default_relaxed_measurement_assignments.csv", [row for row in assignments if row["default_relaxed_candidate"] == "TRUE"], assignment_fields)
    write_csv(output / "candidate_assignment_regulatory_edges.csv", edges, assignment_edge_fields)
    write_csv(output / "measurement_mapping_summary.csv", summaries, feature_fields)
    write_csv(output / "candidate_node_collisions.csv", collisions, collision_fields)
    write_csv(output / "default_relaxed_unique_topology_seeds.csv", seeds, seed_fields)
    write_csv(output / "mapping_policy_summary.csv", policies, policy_fields)
    write_csv(output / "resource_policy_summary.csv", resources, resource_fields)
    write_csv(output / "regulatory_source_inventory.csv", source_inventory, source_fields)
    write_csv(output / "regulatory_source_edges.csv", source_edges, source_edge_fields)
    dependencies = [FEATURES, CANDIDATES, ALIASES, ACID_BASE, OBO, INPUTS, REGULATION, GRAPH, NODES, STITCH_SUMMARY, STITCH_VALIDATION, Path(__file__), Path(stitch_audit.__file__)]
    write_csv(output / "dependency_hashes.csv", [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for path in dependencies
    ], ["path", "sha256"])
    write_readme(output)

    default = [row for row in assignments if row["default_relaxed_candidate"] == "TRUE"]
    summary = {
        "analysis": "factor1_relaxed_lipid_to_regulatory_pkn_mapping",
        "no_moon_scoring_performed": True,
        "input_measurements": len(factor),
        "final_signed_regulatory_source_chemicals": len(contexts),
        "all_candidate_assignments": len(assignments),
        "all_candidate_measurement_features": len({row["measurement_feature"] for row in assignments}),
        "default_relaxed_assignments": len(default),
        "default_relaxed_measurement_features": len({row["measurement_feature"] for row in default}),
        "default_relaxed_source_nodes": len({row["pkn_source_node"] for row in default}),
        "default_relaxed_strict_identity_assignments": sum(row["strict_reported_identity"] == "TRUE" for row in default),
        "default_feature_to_many_node_copies": len({
            row["measurement_feature"] for row in default
            if row["measurement_value_is_copied_to_multiple_default_nodes"] == "TRUE"
        }),
        "default_cross_feature_source_node_collisions": len(collisions),
        "default_mapping_policies": {
            policy: sum(row["primary_mapping_policy"] == policy for row in default)
            for policy in POLICY_ORDER if policy in DEFAULT_PATHS.values()
        },
        "resource_policy": {
            "included": ["GuideToPharmacology", "BRENDA", "STITCH_directional_v5"],
            "excluded": ["MRCLinksDB", "CellPhoneDB_nonpeptidic_adapter", "STITCH_cached_service_layer"],
        },
        "scoring_requirement": (
            "Score each unique source node once as a topology seed; retain collision rows and use "
            "assignment-specific sign consistency and path extraction without counting duplicate node inputs twice."
        ),
        "stitch_caveat": json.loads(STITCH_SUMMARY.read_text())["caveat"],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    validate(output, factor, graph, nodes)
    return summary


def validate(output: Path, factor: dict[str, dict[str, str]], graph: set[tuple[str, str, str]], nodes: dict[str, dict[str, str]]) -> None:
    all_rows = read_csv(output / "all_candidate_regulatory_assignments.csv")
    default = read_csv(output / "default_relaxed_measurement_assignments.csv")
    edges = read_csv(output / "candidate_assignment_regulatory_edges.csv")
    assert len(read_csv(output / "measurement_mapping_summary.csv")) == len(factor)
    assert all(row["candidate_measurement_value"] == factor[row["measurement_feature"]]["oriented_weight"] for row in all_rows)
    assert all(row["pkn_source_node"] in nodes for row in all_rows)
    assert all(row["default_relaxed_candidate"] == "TRUE" for row in default)
    assert all(row["primary_mapping_policy"] in DEFAULT_PATHS.values() for row in default)
    assert all("MRCLinksDB" not in row["surviving_resources"].split(";") for row in all_rows)
    assert {(row["pkn_source_node"], row["interaction"], row["target_node"]) for row in edges} <= graph
    default_ids = {row["assignment_id"] for row in default}
    assert default_ids <= {row["assignment_id"] for row in edges}
    seeds = read_csv(output / "default_relaxed_unique_topology_seeds.csv")
    assert len(seeds) == len({row["candidate_pkn_node"] for row in default})
    assert all(row["unique_topology_seed"] == "TRUE" for row in seeds)


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="factor1_relaxed_regulatory_") as tmp:
        output = Path(tmp) / "audit"
        summary = build(output)
        rows = read_csv(output / "default_relaxed_measurement_assignments.csv")
        lookup = {(row["measurement_feature"], row["pkn_source_chebi"]): row for row in rows}
        assert summary["input_measurements"] == 1254
        assert lookup[("FA(16:0)", "CHEBI:15756")]["primary_mapping_policy"] == "relaxed_acid_base"
        assert lookup[("FA(20:4)", "CHEBI:15843")]["primary_mapping_policy"] == "relaxed_same_formula_is_a_descendant"
        assert lookup[("PC(18:0_18:0)", "CHEBI:83718")]["primary_mapping_policy"] == "relaxed_compatible_more_specific"
        assert not any(row["primary_mapping_policy"].startswith("annotation_only_") for row in rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_DEFAULT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        print("Self-test passed")
    else:
        print(json.dumps(build(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
