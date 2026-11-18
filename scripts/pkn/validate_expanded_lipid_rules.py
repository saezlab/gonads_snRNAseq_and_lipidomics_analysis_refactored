#!/usr/bin/env python3
"""Independently audit additive measured-species Rhea reaction overlays.

This validator reads records and projections; it never builds measurements or
scores a MOON network. Run from any directory with --output, or --self-test.
The comparison is to frozen templates and graphs, not to the builder's counts.
"""
from __future__ import annotations

import argparse
import collections
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

from expanded_lipid_matching import (parse_feature, match_headgroup,
                                    match_lyso_diacyl, match_sphingolipid)

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
OUTPUT = ROOT / "results/causal_modelling/cosmos_mouse_pkn_expanded_species"
TIERS = {"chain_compatible", "glycan_conditional"}
SIDES = ("substrates", "products")
ROLES = {"measured_lipid", "conditional_acyl_pool", "conditional_intermediate"}
# These frozen template chemicals have a fixed d18:1 base (one R substituent).
# Generic R2 glycosylceramides are deliberately NOT assigned this fixed base.
FIXED_BASE = {
    "CHEBI:52639": "18:1;O2", "CHEBI:17636": "18:1;O2",
    "CHEBI:22801": "18:1;O2", "CHEBI:17950": "18:1;O2",
    "CHEBI:18390": "18:1;O2",
    "CHEBI:31488": "18:0;O2", "CHEBI:31998": "18:0;O3",
    "CHEBI:67090": "18:0;O2", "CHEBI:78651": "18:0;O3",
}
TEMPLATE_CLASS = {
    "CHEBI:58608": "PA", "CHEBI:17815": "DG", "CHEBI:49172": "DG",
    "CHEBI:57643": "PC", "CHEBI:64482": "PC", "CHEBI:64612": "PE",
    "CHEBI:57262": "PS", "CHEBI:57880": "PI", "CHEBI:64716": "PG", "CHEBI:60523": "PG",
    "CHEBI:58168": "LPC", "CHEBI:57875": "LPC", "CHEBI:84465": "LPC",
    "CHEBI:64381": "LPE", "CHEBI:65213": "LPE", "CHEBI:64379": "LPS", "CHEBI:65214": "LPS",
    "CHEBI:57970": "LPA", "CHEBI:64982": "LPA", "CHEBI:77589": "LPA",
    "CHEBI:64771": "LPI", "CHEBI:64872": "LPI", "CHEBI:64840": "LPG", "CHEBI:76528": "LPG",
    "CHEBI:52639": "Cer", "CHEBI:31488": "Cer", "CHEBI:31998": "Cer",
    "CHEBI:139048": "Cer", "CHEBI:139051": "Cer", "CHEBI:17636": "SM",
    "CHEBI:64583": "SM", "CHEBI:83273": "Cer",
    "CHEBI:67090": "SM", "CHEBI:78651": "SM",
    "CHEBI:22801": "HexCer", "CHEBI:83264": "HexCer", "CHEBI:36498": "HexCer",
    "CHEBI:18390": "HexCer",
    "CHEBI:143593": "HexCer", "CHEBI:17950": "Hex2Cer", "CHEBI:79208": "Hex2Cer",
}
IMMUTABLE_RECORD = ("catalysts", "catalyst_semantics", "master_reaction_id",
                    "direction", "direction_basis", "reversible",
                    "reversibility_scope", "organism", "provenance",
                    "references", "is_transport")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def jsonlines(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return [json.loads(line) for line in f if line.strip()]


def rows(path, delimiter=","):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


class Audit:
    def __init__(self):
        self.checks = collections.Counter()
        self.failed = collections.Counter()
        self.examples = []

    def check(self, condition, name, context=""):
        self.checks[name] += 1
        if not condition:
            self.failed[name] += 1
            if len(self.examples) < 100:
                self.examples.append({"check": name, "context": str(context)})
        return bool(condition)


def graph(path, audit):
    values = rows(path)
    triples = []
    for i, row in enumerate(values, 2):
        sign = row.get("interaction", "")
        if not audit.check(sign in {"1", "-1", "+1"}, "numeric_binary_graph_sign", f"{path.name}:{i}"):
            continue
        triples.append((row["source"], int(sign), row["target"]))
    audit.check(len(triples) == len(set(triples)), "no_duplicate_graph_triples", path.name)
    pair_signs = collections.defaultdict(set)
    for s, sign, t in triples:
        pair_signs[s, t].add(sign)
        audit.check(s != t, "no_graph_self_loops", f"{path.name}:{s}")
    audit.check(all(len(signs) == 1 for signs in pair_signs.values()),
                "no_conflicting_graph_pair_signs", path.name)
    return set(triples)


def measured_identity(participant, substitution):
    identifier = participant.get("id", "")
    if identifier.startswith("MeasuredLipid:"):
        return identifier[len("MeasuredLipid:"):], None
    if identifier.startswith("CandidateMeasuredLipid:") and "|" in identifier:
        feature, condition = identifier[len("CandidateMeasuredLipid:"):].rsplit("|", 1)
        return feature, condition
    return None, None


def validate_record(record, templates, features, audit):
    rid = record.get("reaction_id", "")
    tid = record.get("supporting_reaction_id")
    if not audit.check(tid in templates, "known_supporting_rhea_template", rid):
        return
    template = templates[tid]
    audit.check(rid and rid != tid, "synthetic_reaction_identity", rid)
    audit.check(record.get("overlay_tier") in TIERS, "recognized_overlay_tier", rid)
    audit.check(template.get("direction") in {"LR", "RL"} and not template.get("is_transport"),
                "explicit_direction_no_transport", rid)
    audit.check(template.get("organism") == 10090, "mouse_template", rid)
    audit.check(bool(template.get("catalysts")), "mouse_catalysts_present", rid)
    for key in IMMUTABLE_RECORD:
        audit.check(record.get(key) == template.get(key), "preserved_record_" + key, rid)
    audit.check(record.get("chemically_balanced") is None and
                record.get("template_chemically_balanced") == template.get("chemically_balanced"),
                "species_atom_balance_not_overclaimed", rid)
    audit.check(record.get("supporting_master_reaction_id", template.get("master_reaction_id")) == template.get("master_reaction_id"),
                "correct_supporting_master", rid)
    audit.check(bool(record.get("assumptions")), "explicit_hypothesis_assumptions", rid)

    declared = {}
    substitutions = record.get("substitutions", [])
    audit.check(bool(substitutions), "nonempty_substitutions", rid)
    measured = []
    for sub in substitutions:
        side, index = sub.get("side"), sub.get("index")
        if not audit.check(side in SIDES and isinstance(index, int) and not isinstance(index, bool)
                           and 0 <= index < len(template.get(side, [])), "valid_substitution_index", rid):
            continue
        audit.check((side, index) not in declared, "unique_declared_substitution", rid)
        declared[side, index] = sub
        original = template[side][index]
        if index >= len(record.get(side, [])):
            audit.check(False, "declared_participant_present", rid)
            continue
        participant = record[side][index]
        audit.check(sub.get("template_chemical_id") == original.get("id"), "declared_template_chemical", rid)
        audit.check(sub.get("new_id") == participant.get("id"), "declared_new_chemical", rid)
        audit.check(sub.get("role") in ROLES, "recognized_substitution_role", rid)
        audit.check(participant.get("id") != original.get("id"), "declared_actual_substitution", rid)
        audit.check(participant.get("template_chemical_id") == original.get("id"), "participant_template_reference", rid)
        for original_key, stored_key in (("name", "template_name"), ("formula", "template_formula"), ("rhea_compound", "template_rhea_compound")):
            audit.check(participant.get(stored_key) == original.get(original_key), "preserved_" + stored_key, rid)
        for key in ("stoichiometry", "compartment", "participant_id", "side_id", "relative_location", "polymerization_index", "reactive_parts"):
            audit.check(participant.get(key) == original.get(key), "preserved_substituted_" + key, rid)
        if sub.get("role") == "measured_lipid":
            fid, condition = measured_identity(participant, sub)
            if not audit.check(fid in features and features.get(fid) is not None,
                               "parser_supported_measured_substitution", f"{rid}:{participant.get('id')}"):
                continue
            f = features[fid]
            measured.append(f)
            audit.check(participant.get("reported_chains") == f["chains"], "reported_measured_chain_roles_preserved", rid)
            audit.check(TEMPLATE_CLASS.get(original["id"]) == f["class"],
                        "measured_class_matches_template_chemical", f"{rid}:{fid}:{original['id']}")
            if participant.get("measurement_feature_id") is not None:
                audit.check(participant["measurement_feature_id"] == fid, "measured_feature_annotation", rid)
            if f["class"] in {"HexCer", "Hex2Cer"}:
                audit.check(record.get("overlay_tier") == "glycan_conditional", "glycan_stays_conditional", rid)
                audit.check(condition in {"glucose", "galactose", "lactosyl"}, "separate_glycan_candidate_identity", rid)
                audit.check(participant.get("candidate_identity_condition") == condition, "candidate_glycan_annotation", rid)
                audit.check((f["class"] == "HexCer" and condition in {"glucose", "galactose"}) or
                            (f["class"] == "Hex2Cer" and condition == "lactosyl"), "glycan_class_condition", rid)
                name = original.get("name", "").lower()
                if condition == "glucose":
                    audit.check("glucos" in name and "galactos" not in name, "glucose_template_identity", rid)
                elif condition == "galactose":
                    audit.check("galactos" in name or "gal-" in name, "galactose_template_identity", rid)
                elif condition == "lactosyl":
                    audit.check(("gal" in name and "glc" in name) or ("galactos" in name and "glucos" in name),
                                "lactosyl_template_identity", rid)
            else:
                audit.check(condition is None and participant.get("id", "").startswith("MeasuredLipid:"),
                            "ordinary_measured_identity", rid)
            if f["family"] == "sphingolipid":
                fixed = FIXED_BASE.get(original["id"])
                if fixed:
                    audit.check(f["sphingobase"] == fixed, "fixed_template_sphingobase", f"{rid}:{fid}")
                if original["id"] in {"CHEBI:64583", "CHEBI:83273", "CHEBI:83264", "CHEBI:36498", "CHEBI:79208", "CHEBI:143593"}:
                    audit.check(f["chain_objects"][0]["oxygen"] == 2 and f["chain_objects"][1]["oxygen"] == 0,
                                "generic_template_preserves_unmodified_sphingoid_oxygen", f"{rid}:{fid}")
                if f["oxygenated_acyl_requires_template_review"]:
                    audit.check("hydroxy" in original.get("name", "").lower() and "acyl" in original.get("name", "").lower(),
                                "oxygenated_n_acyl_template_support", rid)
        elif sub.get("role") == "conditional_acyl_pool":
            audit.check(bool(re.fullmatch(r"HypotheticalAcyl:(FA|ACoA)\(\d+:\d+\)", participant.get("id", ""))),
                        "conditional_acyl_identity_namespace", rid)
            audit.check("acyl" in original.get("name", "").lower() or "fatty acid" in original.get("name", "").lower(),
                        "acyl_pool_supported_by_template", rid)
            if "ACoA(" in participant.get("id", ""):
                audit.check(original["id"] in {"CHEBI:58342", "CHEBI:83139", "CHEBI:77636"},
                            "unmodified_acyl_coa_template", rid)
            else:
                audit.check(original["id"] in {"CHEBI:28868", "CHEBI:57560", "CHEBI:59558", "CHEBI:58950"},
                            "unmodified_fatty_acid_template", rid)
            match = re.fullmatch(r"HypotheticalAcyl:(FA|ACoA)\((\d+):(\d+)\)", participant.get("id", ""))
            if match:
                carrier, carbon, double_bonds = match[1], int(match[2]), int(match[3])
                expected_formula = (f"C{carbon+21}H{2*carbon+30-2*double_bonds}N7O17P3S" if carrier == "ACoA"
                                    else f"C{carbon}H{2*carbon-1-2*double_bonds}O2")
                audit.check(participant.get("formula") == expected_formula, "conditional_pool_composition_formula", rid)
                audit.check(participant.get("conditional_acyl_chain") == f"{carbon}:{double_bonds}",
                            "conditional_pool_chain_annotation", rid)
            audit.check(participant.get("unmeasured") is True and not participant.get("measurement_feature_id"),
                        "conditional_pool_not_measurement", rid)
        elif sub.get("role") == "conditional_intermediate":
            audit.check(participant.get("id", "").startswith("HypotheticalLipid:"), "conditional_intermediate_namespace", rid)
            match = re.fullmatch(r"HypotheticalLipid:(PE-NMe2?)\((\d+:\d+)[/_](\d+:\d+)\)", participant.get("id", ""))
            audit.check(bool(match), "pemt_intermediate_explicit_chain_identity", rid)
            if match:
                audit.check(original["id"] == {"PE-NMe": "CHEBI:64573", "PE-NMe2": "CHEBI:64572"}[match[1]],
                            "pemt_intermediate_stage_matches_template", rid)
                context = record.get("match_details", {}).get("chain_multiset", [])
                audit.check(sorted([match[2], match[3]]) == sorted(context), "pemt_intermediate_preserves_context_chains", rid)
            audit.check(participant.get("unmeasured") is True and not participant.get("measurement_feature_id"),
                        "conditional_intermediate_not_measurement", rid)

    for side in SIDES:
        audit.check(len(record.get(side, [])) == len(template[side]), "full_reaction_participant_count", f"{rid}:{side}")
        for index, original in enumerate(template[side]):
            if index >= len(record.get(side, [])):
                continue
            if (side, index) not in declared:
                audit.check(record[side][index] == original, "undeclared_participant_untouched", f"{rid}:{side}:{index}")
    actual_fids = {f["feature_id"] for f in measured}
    listed = record.get("matched_feature_ids", [])
    audit.check(set(listed) == actual_fids and len(listed) == len(set(listed)), "actual_measured_participants_only", rid)
    validate_matching(record, measured, template, features, audit)


def validate_matching(record, measured, template, features, audit):
    """Recompute structural compatibility without importing the builder."""
    rid = record["reaction_id"]
    details = record.get("match_details", {})
    audit.check(isinstance(details, dict) and bool(details), "nonempty_structural_match_details", rid)
    # Pemt steps may have only one measured endpoint and conditional pathway
    # intermediates, so their structural context is recorded separately.
    if len(measured) < 2:
        audit.check(record.get("rule_family", "").startswith("pemt") or
                    bool(record.get("context_pathway_feature_ids")), "single_endpoint_has_pathway_context", rid)
        context = record.get("context_pathway_feature_ids", [])
        audit.check(len(context) >= 2 and all(fid in features and features[fid] for fid in context),
                    "explicit_pemt_pathway_pair", rid)
        anchor = details.get("primary", details)
        pair = [anchor.get("source_feature_id"), anchor.get("target_feature_id")]
        if all(fid in context and fid in features and features[fid] for fid in pair):
            measured = [features[fid] for fid in pair]
            for fid in context:
                f = features.get(fid)
                if not f:
                    continue
                audit.check(f["chain_multiset"] == measured[0]["chain_multiset"], "shared_pemt_context_exact_chains", rid)
                for a in measured:
                    if f["sn_order_reported"] and a["sn_order_reported"]:
                        audit.check(f["chains"] == a["chains"], "shared_pemt_context_reported_sn", rid)
        else:
            audit.check(False, "pemt_recorded_match_in_context", rid)
            return
    candidates = []
    for i, a in enumerate(measured):
        for b in measured[i + 1:]:
            if a["family"] == b["family"] == "glycerol_diacyl":
                hit = match_headgroup(a, b)
            elif a["family"] == b["family"] == "sphingolipid":
                hit = match_sphingolipid(a, b)
            else:
                lyso = a if a["family"] == "glycerol_monoacyl" else b
                diacyl = b if a["family"] == "glycerol_monoacyl" else a
                hit = match_lyso_diacyl(lyso, diacyl)
            if hit:
                candidates.append(hit)
    audit.check(bool(candidates), "measured_pair_structurally_compatible", rid)
    recorded = [details.get("primary", details)]
    recorded.extend(details.get("coupled_matches", []))
    recorded.extend(record.get("coupled_matches", []))
    # Match records may carry template-specific gates in addition to the helper
    # result. Require their defining chain identities and positional choices.
    for hit in candidates:
        identity_fields = [key for key in ("source_feature_id", "target_feature_id", "lyso_feature_id", "diacyl_feature_id") if key in hit]
        # Directional templates can request the same headgroup structural check
        # in either order; chain chemistry remains symmetric at this gate.
        hit_ids = {hit[key] for key in identity_fields}
        saved = [value for value in recorded if isinstance(value, dict)
                 and {value[key] for key in identity_fields if key in value} == hit_ids]
        audit.check(bool(saved), "recomputed_match_has_recorded_feature_pair", rid)
        if not saved:
            continue
        value = saved[0]
        critical = [key for key in ("rule_family", "chain_multiset", "sphingobase", "n_acyl_chain", "retained_chain", "added_chain", "removed_chain",
                                   "sn_assignment_conditional", "headgroup_isomer_conditional", "oxygenated_acyl_requires_template_review",
                                   "sphingobase_oxygen", "n_acyl_oxygen", "lyso_retained_sn_position") if key in hit]
        for key in critical:
            audit.check(value.get(key) == hit[key], "recomputed_match_detail_" + key, rid)
        if hit["rule_family"] == "lyso_diacyl_acyl_exchange":
            for key in ("retained_chain", "added_chain", "removed_chain"):
                if key in value:
                    audit.check(value[key] == hit[key], "match_detail_" + key, rid)
            retained_sns = {int(re.search(r"([12])-acyl-sn", p.get("name", "")).group(1))
                            for side in SIDES for p in template[side]
                            if re.search(r"([12])-acyl-sn", p.get("name", ""))}
            # Full diacyl names include '1,2-diacyl' rather than '1-acyl', so
            # this gate only reads explicit monoacyl template participants.
            options = hit["positional_options"]
            if retained_sns:
                audit.check(any(option["retained_sn"] in retained_sns for option in options),
                            "template_retained_sn_compatible", rid)
            expected_options = [option for option in options if not retained_sns or option["retained_sn"] in retained_sns]
            audit.check(value.get("positional_options") == expected_options,
                        "exact_template_filtered_positional_options", rid)
            acyl_ids = [sub.get("new_id", "") for sub in record.get("substitutions", [])
                        if sub.get("role") == "conditional_acyl_pool"]
            audit.check(bool(acyl_ids) and all("(" + hit["added_chain"] + ")" in ident for ident in acyl_ids),
                        "complete_exact_added_removed_acyl_pool", rid)
        elif hit["rule_family"] == "sphingolipid_headgroup":
            for key in ("sphingobase", "n_acyl_chain"):
                if key in value:
                    audit.check(value[key] == hit[key], "match_detail_" + key, rid)


def protein_map():
    path = ROOT / "support/cosmos_mouse/raw/uniprot_mouse.tsv.gz"
    mapping = {}
    for row in rows(path, "\t"):
        symbol = row["Gene Names (primary)"]
        if re.fullmatch(r"[A-Za-z0-9_.-]+", symbol):
            mapping[row["Entry"]] = "Gene__" + symbol
    return mapping


def protein_id(accession, mapping):
    base = re.sub(r"-\d+$", "", accession)
    return mapping.get(accession, mapping.get(base, "UniProt__" + accession))


def projected_record(record, mapping):
    triples, evidence = set(), []
    for catalyst in sorted({protein_id(p, mapping) for p in record["catalysts"]}):
        signature = json.dumps([record["resource"], record["reaction_id"], False, [catalyst]])
        node = "Enzyme__" + hashlib.sha256(signature.encode()).hexdigest()[:16]
        triples.add((catalyst, 1, node))
        evidence.append((catalyst, 1, node, "catalyst_enables_reaction", None))
        for side in SIDES:
            for p in record[side]:
                metab = "Metab__" + p["id"] + "_" + (p.get("compartment") or "unknown")
                triple = (metab, 1, node) if side == "substrates" else (node, 1, metab)
                triples.add(triple)
                evidence.append((*triple, "substrate_participation" if side == "substrates" else "product_formation", p["stoichiometry"]))
    return triples, evidence


def validate(output, base):
    audit = Audit()
    source_files = [base / "rhea_reactions.jsonl", base / "mapping_features.csv", base / "pkn_lipid_rules.csv",
                    base / "nodes.jsonl", output / "expanded_species_reactions.jsonl", output / "nodes.jsonl",
                    output / "extension_edge_evidence.jsonl.gz", output / "pkn_expanded_chain.csv",
                    output / "pkn_expanded_glycan_conditional.csv", Path(__file__).resolve(),
                    Path(__file__).with_name("expanded_lipid_matching.py").resolve(),
                    ROOT / "support/cosmos_mouse/raw/uniprot_mouse.tsv.gz"]
    for path in source_files:
        if not path.exists():
            raise FileNotFoundError(path)
    templates = {r["reaction_id"]: r for r in jsonlines(base / "rhea_reactions.jsonl")}
    features = {r["feature_id"]: parse_feature(r) for r in rows(base / "mapping_features.csv")}
    records = jsonlines(output / "expanded_species_reactions.jsonl")
    ids = [r.get("reaction_id") for r in records]
    audit.check(bool(records), "nonempty_extended_master_records")
    audit.check(len(ids) == len(set(ids)), "unique_extended_master_ids")
    for record in records:
        validate_record(record, templates, features, audit)
    nodes_rows = jsonlines(output / "nodes.jsonl")
    nodes = {n["node_id"]: n for n in nodes_rows}
    audit.check(len(nodes_rows) == len(nodes), "unique_node_ids")
    base_nodes = {n["node_id"]: n for n in jsonlines(base / "nodes.jsonl")}
    mapping = protein_map()
    base_edges = graph(base / "pkn_lipid_rules.csv", audit)
    active_base_ids = {n for s, sign, t in base_edges for n in (s, t)}
    audit.check(all(nodes.get(k) == base_nodes.get(k) for k in active_base_ids), "baseline_node_records_unchanged")
    by_tier = collections.defaultdict(set)
    by_rid = {}
    expected_evidence = collections.Counter()
    for record in records:
        projected, evidence = projected_record(record, mapping)
        by_tier[record.get("overlay_tier")].update(projected)
        by_rid[record["reaction_id"]] = projected
        for s, sign, t, mechanism, stoich in evidence:
            expected_evidence[record["reaction_id"], s, sign, t, mechanism, stoich] += 1
        for s, sign, t in projected:
            audit.check(s in nodes and t in nodes, "full_reaction_projection_nodes_present", record["reaction_id"])
        for side in SIDES:
            for participant in record[side]:
                mid = "Metab__" + participant["id"] + "_" + (participant.get("compartment") or "unknown")
                node = nodes.get(mid, {})
                audit.check(node.get("node_type") == "metabolite" and node.get("chemical_id") == participant["id"]
                            and node.get("compartment") == (participant.get("compartment") or "unknown"),
                            "metabolite_node_chemical_identity", mid)
                if mid not in base_nodes:
                    audit.check(all(node.get(key) == participant.get(key, False if key == "unmeasured" else None)
                                    for key in ("measurement_feature_id", "candidate_identity_condition", "unmeasured")),
                                "new_metabolite_resolution_annotations", mid)
        enzyme_ids = {t for s, sign, t in projected if s.startswith(("Gene__", "UniProt__"))}
        for eid in enzyme_ids:
            node = nodes.get(eid, {})
            audit.check(node.get("node_type") == "enzyme_reaction_instance" and node.get("reaction_id") == record["reaction_id"]
                        and node.get("resource") == record.get("resource") and node.get("reverse") is False
                        and node.get("uncertain_direction") is False, "projected_enzyme_master_reference", eid)
            expected_catalysts = sorted({s for s, sign, t in projected if t == eid and s.startswith(("Gene__", "UniProt__"))})
            audit.check(node.get("catalysts") == expected_catalysts and node.get("complex_required") is False
                        and node.get("overlay_tier") == record.get("overlay_tier"), "enzyme_catalyst_and_tier_metadata", eid)
    chain = graph(output / "pkn_expanded_chain.csv", audit)
    conditional = graph(output / "pkn_expanded_glycan_conditional.csv", audit)
    for observed, expected, label in [(chain, base_edges | by_tier["chain_compatible"], "chain"),
                                      (conditional, base_edges | by_tier["chain_compatible"] | by_tier["glycan_conditional"], "glycan")]:
        audit.check(base_edges <= observed, "all_baseline_edges_retained_" + label)
        audit.check(observed == expected, "exact_additive_supergraph_" + label,
                    f"missing={len(expected-observed)} extra={len(observed-expected)}")
        audit.check(all(s in nodes and t in nodes for s, sign, t in observed), "all_graph_nodes_in_node_table_" + label)
    ledger = jsonlines(output / "extension_edge_evidence.jsonl.gz")
    record_lookup = {r["reaction_id"]: r for r in records}
    observed_evidence = collections.Counter()
    for row in ledger:
        rid = row.get("reaction_id")
        audit.check(rid in by_rid, "evidence_references_extended_master", rid)
        audit.check(row.get("master_record") == "expanded_species_reactions.jsonl", "evidence_master_filename", rid)
        audit.check(row.get("interaction") in {-1, 1}, "evidence_binary_sign", rid)
        record = record_lookup.get(rid)
        audit.check(record is not None and row.get("annotated_direction") == record.get("direction")
                    and row.get("reverse", False) is False and row.get("uncertain_direction", False) is False,
                    "evidence_direction_matches_template", rid)
        audit.check(record is not None and all(row.get(key) == record.get(key)
                    for key in ("supporting_reaction_id", "resource", "overlay_tier", "source_version")),
                    "evidence_matches_master_provenance", rid)
        observed_evidence[rid, row.get("source"), row.get("interaction"), row.get("target"), row.get("mechanism"), row.get("stoichiometry")] += 1
    audit.check(observed_evidence == expected_evidence, "exact_full_participant_catalyst_evidence_projection",
                f"missing={sum((expected_evidence-observed_evidence).values())} extra={sum((observed_evidence-expected_evidence).values())}")
    return {
        "passed": not audit.failed, "check_counts": dict(sorted(audit.checks.items())),
        "failed_checks": dict(sorted(audit.failed.items())), "failure_examples": audit.examples,
        "counts": {"extended_reaction_records": len(records), "records_by_tier": dict(collections.Counter(r.get("overlay_tier") for r in records)),
                   "base_edges": len(base_edges), "chain_graph_edges": len(chain), "glycan_graph_edges": len(conditional),
                   "extension_evidence_rows": len(ledger)},
        "dependency_sha256": {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): digest(p) for p in source_files},
        "interpretation": "Structural and projection invariants only. Matching does not establish species-level catalyst specificity, activity signs, flux, substrate availability or independent mechanistic validation.",
    }


def self_test():
    """Exercise independent template-invariant failures with minimal fixtures."""
    original = {"id": "CHEBI:57643", "name": "a 1,2-diacyl-sn-glycero-3-phosphocholine",
                "stoichiometry": 1, "compartment": "unknown", "participant_id": "a", "side_id": "left"}
    product = dict(original, id="CHEBI:64612", name="a 1,2-diacyl-sn-glycero-3-phosphoethanolamine", side_id="right")
    template = {"reaction_id": "RHEA:1", "master_reaction_id": "RHEA:0", "direction": "LR", "is_transport": False,
                "organism": 10090, "catalysts": ["P1"], "substrates": [original], "products": [product]}
    record = copy.deepcopy(template)
    record.update(reaction_id="TEST:1", supporting_reaction_id="RHEA:1", overlay_tier="chain_compatible",
                  rule_family="diacyl_headgroup", resource="test", assumptions=["species specificity is conditional"],
                  chemically_balanced=None, template_chemically_balanced=template.get("chemically_balanced"),
                  matched_feature_ids=["PC(16:0_18:1)", "PE(16:0_18:1)"], substitutions=[])
    features = {}
    for side, cls in [("substrates", "PC"), ("products", "PE")]:
        fid = cls + "(16:0_18:1)"
        features[fid] = parse_feature({"feature_id": fid, "canonical_name": cls + " 16:0_18:1", "class": cls})
        p = record[side][0]
        old = p["id"]
        p.update(id="MeasuredLipid:" + fid, template_chemical_id=old, reported_chains=features[fid]["chains"], template_name=p.get("name"),
                 template_formula=p.get("formula"), template_rhea_compound=p.get("rhea_compound"))
        record["substitutions"].append(dict(side=side, index=0, template_chemical_id=old, new_id=p["id"], role="measured_lipid"))
    record["match_details"] = match_headgroup(features["PC(16:0_18:1)"], features["PE(16:0_18:1)"])
    a = Audit(); validate_record(record, {"RHEA:1": template}, features, a)
    assert not a.failed, a.examples
    broken = copy.deepcopy(record); broken["products"][0]["stoichiometry"] = 2
    a = Audit(); validate_record(broken, {"RHEA:1": template}, features, a)
    assert a.failed["preserved_substituted_stoichiometry"] == 1
    broken = copy.deepcopy(record); broken["substitutions"].pop()
    a = Audit(); validate_record(broken, {"RHEA:1": template}, features, a)
    assert a.failed["undeclared_participant_untouched"] == 1
    broken = copy.deepcopy(record); broken["catalysts"] = ["wrong"]
    a = Audit(); validate_record(broken, {"RHEA:1": template}, features, a)
    assert a.failed["preserved_record_catalysts"] == 1
    broken = copy.deepcopy(record); broken["matched_feature_ids"].append("unmeasured")
    a = Audit(); validate_record(broken, {"RHEA:1": template}, features, a)
    assert a.failed["actual_measured_participants_only"] == 1
    print("Independent expanded-species validator self-tests passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--report", type=Path, help="Optional JSON report destination; never alters graph inputs.")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    report = validate(args.output, args.base)
    if args.report:
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
