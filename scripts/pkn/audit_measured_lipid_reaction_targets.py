#!/usr/bin/env python3
"""Audit measured lipid features as *reaction-product* targets in the frozen PKN.

This is deliberately narrower than the existing chemical-identity and lipid
regulatory audits. A measured feature can be present in a species-rule reaction
as either a substrate or a product; only an incoming ``product_formation`` edge
makes it a downstream reaction endpoint in the directed metabolic projection.

The script is read-only with respect to the PKN and frozen MOON result. It
creates an audit directory with a feature-level inventory, rule-level direction
checks, population summaries and hashes of every dependency.

Run from the repository root:

    python3 scripts/pkn/audit_measured_lipid_reaction_targets.py \\
      --output results/causal_modelling/lipid_reaction_target_mapping_audit
"""

import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
VARIANT_RANK = {"core": 0, "lipid_curated": 1, "lipid_rules": 2, "exploratory": 3}


def read_csv(path, delimiter=","):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def json_lines(path):
    with path.open() as handle:
        return [json.loads(line) for line in handle if line.strip()]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def truth(value):
    return str(value).lower() == "true"


def node_chemical_id(node):
    """Return the stable chemical identity encoded by a metabolite node."""
    if not node.startswith("Metab__"):
        return ""
    return node[len("Metab__"):].rsplit("_", 1)[0]


def measured_node(feature):
    return "Metab__MeasuredLipid:" + feature + "_unknown"


def role(source, product):
    if source and product:
        return "both"
    if source:
        return "substrate_only"
    if product:
        return "product_only"
    return "not_in_species_rule"


def product_evidence(pkn_dir, active_edges):
    """Return active and inactive product evidence up to the lipid-rule PKN."""
    active = collections.defaultdict(list)
    inactive = []
    all_records = []
    with gzip.open(pkn_dir / "edge_evidence.jsonl.gz", "rt") as handle:
        for line in handle:
            record = json.loads(line)
            if (VARIANT_RANK[record["first_variant"]] > VARIANT_RANK["lipid_rules"] or
                    record.get("mechanism") != "product_formation"):
                continue
            key = (record["source"], int(record["interaction"]), record["target"])
            all_records.append(record)
            if key in active_edges:
                active[record["target"]].append(record)
            else:
                inactive.append(record)
    return active, inactive, all_records


def product_summary(name, feature_ids, source_features, product_features,
                    exact_product_features, accepted_by_feature, product_nodes):
    feature_ids = set(feature_ids)
    accepted = [accepted_by_feature[f] for f in feature_ids if f in accepted_by_feature]
    accepted_product = [row for row in accepted if row["node_id"] in product_nodes]
    return {
        "population": name,
        "features": len(feature_ids),
        "species_rule_connected": len(feature_ids & (source_features | product_features)),
        "species_rule_substrates": len(feature_ids & source_features),
        "species_rule_products": len(feature_ids & product_features),
        "species_rule_substrate_only": len(feature_ids & (source_features - product_features)),
        "species_rule_product_only": len(feature_ids & (product_features - source_features)),
        "species_rule_both": len(feature_ids & source_features & product_features),
        "exact_chemical_product_candidates": len(feature_ids & exact_product_features),
        "any_product_candidate": len(feature_ids & (product_features | exact_product_features)),
        "accepted_input_nodes": len(accepted),
        "accepted_input_product_nodes": len(accepted_product),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pkn", type=Path,
        default=ROOT / "results/causal_modelling/cosmos_mouse_pkn",
        help="Frozen PKN directory (default: results/causal_modelling/cosmos_mouse_pkn)",
    )
    parser.add_argument(
        "--input-mapping", type=Path,
        default=ROOT / "miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv",
        help="Frozen all-feature lipid-to-node mapping CSV",
    )
    parser.add_argument(
        "--tf-to-lipid-run", type=Path,
        default=ROOT / "results/causal_modelling/factor1_tf_to_lipid",
        help="Selected TF-to-lipid MOON result directory",
    )
    parser.add_argument("--output", type=Path, required=True,
                        help="New or empty audit directory")
    args = parser.parse_args()

    pkn = args.pkn.resolve()
    input_mapping = args.input_mapping.resolve()
    run = args.tf_to_lipid_run.resolve()
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Output directory is occupied: " + str(output))

    dependencies = [
        pkn / "pkn_lipid_rules.csv",
        pkn / "edge_evidence.jsonl.gz",
        pkn / "mapping_features.csv",
        pkn / "mapping_chebi_nodes.csv",
        pkn / "lipid_rules_pairs.csv",
        pkn / "lipid_rules_reactions.jsonl",
        input_mapping,
        run / "downstream_input.csv",
        run / "network_ATT.csv",
        Path(__file__).resolve(),
    ]
    missing = [str(path) for path in dependencies if not path.is_file()]
    if missing:
        parser.error("Missing dependency: " + ", ".join(missing))

    pkn_rows = read_csv(pkn / "pkn_lipid_rules.csv")
    active_edges = {(row["source"], int(row["interaction"]), row["target"])
                    for row in pkn_rows}
    if len(active_edges) != len(pkn_rows):
        raise AssertionError("pkn_lipid_rules.csv contains duplicate signed triples")

    active_products, inactive_products, _ = product_evidence(pkn, active_edges)
    product_nodes = set(active_products)
    if inactive_products:
        raise AssertionError("A product_formation evidence triple was lost from pkn_lipid_rules")

    feature_rows = read_csv(pkn / "mapping_features.csv")
    feature_ids = {row["feature_id"] for row in feature_rows}
    if len(feature_ids) != len(feature_rows):
        raise AssertionError("mapping_features.csv must have one row per feature")
    factor1_features = {row["feature_id"] for row in feature_rows if truth(row["in_factor1"])}

    # ChEBI candidates remain candidates unless they are exact at the reported
    # resolution. Compatible and hierarchy candidates are intentionally absent.
    exact_chebi = collections.defaultdict(set)
    for row in read_csv(pkn / "mapping_chebi_nodes.csv"):
        if (row["mapping_relation"] == "SwissLipids identifier" and
                row["match_type"] == "exact_reported_resolution"):
            exact_chebi[row["feature_id"]].add(row["node_chebi_id"])
    product_chebi_nodes = collections.defaultdict(set)
    for node in product_nodes:
        chemical = node_chemical_id(node)
        if chemical.startswith("CHEBI:"):
            product_chebi_nodes[chemical].add(node)
    exact_chemical_product_nodes = {
        feature: sorted({node for chemical in chemicals
                         for node in product_chebi_nodes.get(chemical, set())})
        for feature, chemicals in exact_chebi.items()
    }
    exact_chemical_product_nodes = {
        feature: nodes for feature, nodes in exact_chemical_product_nodes.items() if nodes
    }
    exact_product_features = set(exact_chemical_product_nodes)

    pairs = read_csv(pkn / "lipid_rules_pairs.csv")
    records = {record["reaction_id"]: record
               for record in json_lines(pkn / "lipid_rules_reactions.jsonl")}
    if len(records) != len(pairs):
        raise AssertionError("Rule-pair and rule-record counts differ")
    if set(records) != {row["reaction_id"] for row in pairs}:
        raise AssertionError("Every rule pair must have exactly one reaction record")

    source_reactions = collections.defaultdict(list)
    product_reactions = collections.defaultdict(list)
    rule_checks = []
    for pair in pairs:
        reaction = records[pair["reaction_id"]]
        source = pair["source_feature"]
        target = pair["target_feature"]
        source_reactions[source].append(pair["reaction_id"])
        product_reactions[target].append(pair["reaction_id"])
        substrate_ids = {participant["id"] for participant in reaction["substrates"]}
        product_ids = {participant["id"] for participant in reaction["products"]}
        source_on_substrate_side = "MeasuredLipid:" + source in substrate_ids
        target_on_product_side = "MeasuredLipid:" + target in product_ids
        expected_target_node = measured_node(target)
        evidence = [record for record in active_products.get(expected_target_node, [])
                    if record["resource"] == "LINEX2_rule" and
                    record["reaction_id"] == pair["reaction_id"]]
        evidence_triples = {(record["source"], int(record["interaction"]), record["target"])
                            for record in evidence}
        catalyst_count = len(reaction["catalysts"])
        rule_checks.append({
            "reaction_id": pair["reaction_id"],
            "supporting_reaction_id": pair["supporting_reaction_id"],
            "headgroup_transform": pair["headgroup_transform"],
            "source_feature": source,
            "target_feature": target,
            "declared_direction": reaction["direction"],
            "irreversible_rule_projection": not reaction["reversible"],
            "source_replaced_on_substrate_side": source_on_substrate_side,
            "target_replaced_on_product_side": target_on_product_side,
            "expected_target_node": expected_target_node,
            "catalyst_count": catalyst_count,
            "active_catalyst_specific_product_edges": len(evidence_triples),
            "all_product_edges_retained": len(evidence_triples) == catalyst_count,
        })

    # One active product triple must remain for every alternative catalyst. This
    # checks side assignment and final signed-conflict handling independently.
    if not all(row["source_replaced_on_substrate_side"] for row in rule_checks):
        raise AssertionError("A measured rule source is not on the substrate side")
    if not all(row["target_replaced_on_product_side"] for row in rule_checks):
        raise AssertionError("A measured rule target is not on the product side")
    if not all(row["irreversible_rule_projection"] and
               row["declared_direction"] in {"LR", "RL"} for row in rule_checks):
        raise AssertionError("A rule does not retain its declared directional projection")
    if not all(row["all_product_edges_retained"] for row in rule_checks):
        raise AssertionError("A rule product edge was lost from the active PKN")

    source_features = set(source_reactions)
    rule_product_features = set(product_reactions)
    expected_rule_product_nodes = {measured_node(feature) for feature in rule_product_features}
    actual_rule_product_nodes = {
        node for node, evidence in active_products.items()
        if node.startswith("Metab__MeasuredLipid:") and
        any(record["resource"] == "LINEX2_rule" for record in evidence)
    }
    if expected_rule_product_nodes != actual_rule_product_nodes:
        raise AssertionError("Rule product-node inventory differs from active evidence")

    mapping_rows = read_csv(input_mapping)
    mapping_by_feature = {row["feature"]: row for row in mapping_rows}
    if set(mapping_by_feature) != feature_ids:
        raise AssertionError("Frozen input mapping and feature inventory differ")
    if len(mapping_by_feature) != len(mapping_rows):
        raise AssertionError("Frozen input mapping must have one row per feature")
    included = {feature: row for feature, row in mapping_by_feature.items()
                if truth(row["included"])}
    if len({row["node_id"] for row in included.values()}) != len(included):
        raise AssertionError("Accepted features must have unique nodes")

    pruned_rows = read_csv(run / "downstream_input.csv")
    pruned_nodes = {row["source"] for row in pruned_rows}
    if not pruned_nodes <= {row["node_id"] for row in included.values()}:
        raise AssertionError("Pruned downstream input is not a subset of frozen accepted nodes")
    solution_rows = read_csv(run / "network_ATT.csv")
    solution_nodes = {row["source"] for row in solution_rows if float(row["level"]) == 0}
    if not solution_nodes <= pruned_nodes:
        raise AssertionError("Selected level-0 observations are not retained downstream inputs")

    inventory = []
    for feature in feature_rows:
        name = feature["feature_id"]
        accepted = included.get(name)
        exact_nodes = exact_chemical_product_nodes.get(name, [])
        rule_source = name in source_features
        rule_product = name in rule_product_features
        node = accepted["node_id"] if accepted else ""
        inventory.append({
            "feature_id": name,
            "in_lipid_de": feature["in_lipid_de"],
            "in_historical_factor1_64": feature["in_factor1"],
            "best_mapping_status": feature["best_mapping_status"],
            "nomenclature_level": feature["nomenclature_level"],
            "exact_chebi_ids": feature["exact_chebi_ids"],
            "species_rule_role": role(rule_source, rule_product),
            "rule_source_reaction_ids": "|".join(sorted(source_reactions[name])),
            "rule_product_reaction_ids": "|".join(sorted(product_reactions[name])),
            "rule_product_node": measured_node(name) if rule_product else "",
            "rule_product_active_edges": sum(
                1 for record in active_products.get(measured_node(name), [])
                if record["resource"] == "LINEX2_rule"
            ),
            "exact_chemical_product_nodes": "|".join(exact_nodes),
            "exact_chemical_product_active_edges": sum(
                len(active_products[node]) for node in exact_nodes
            ),
            "has_product_candidate": bool(rule_product or exact_nodes),
            "accepted_frozen_input": bool(accepted),
            "accepted_node": node,
            "accepted_mapping_method": accepted["mapping_method"] if accepted else "",
            "accepted_node_is_active_product": bool(accepted and node in product_nodes),
            "retained_after_selected_topology_pruning": node in pruned_nodes,
            "selected_solution_level0_endpoint": node in solution_nodes,
        })

    all_features = set(feature_ids)
    population_rows = [
        product_summary("all_measured_features", all_features, source_features,
                        rule_product_features, exact_product_features,
                        included, product_nodes),
        product_summary("historical_factor1_64", factor1_features, source_features,
                        rule_product_features, exact_product_features,
                        included, product_nodes),
    ]
    accepted_features = set(included)
    population_rows.append(product_summary("frozen_accepted_inputs", accepted_features,
                                           source_features, rule_product_features,
                                           exact_product_features, included,
                                           product_nodes))
    pruned_features = {feature for feature, row in included.items()
                       if row["node_id"] in pruned_nodes}
    population_rows.append(product_summary("selected_topology_retained_inputs", pruned_features,
                                           source_features, rule_product_features,
                                           exact_product_features, included,
                                           product_nodes))
    solution_features = {feature for feature, row in included.items()
                         if row["node_id"] in solution_nodes}
    population_rows.append(product_summary("selected_solution_level0_endpoints", solution_features,
                                           source_features, rule_product_features,
                                           exact_product_features, included,
                                           product_nodes))

    accepted_product = {feature for feature, row in included.items()
                        if row["node_id"] in product_nodes}
    accepted_nonproduct = accepted_features - accepted_product
    candidate_only = (rule_product_features | exact_product_features) - accepted_product
    exact_input_product = {feature for feature, row in included.items()
                           if row["mapping_method"] == "unambiguous_exact_reported_resolution_chebi" and
                           row["node_id"] in product_nodes}
    rule_input_product = {feature for feature, row in included.items()
                          if row["mapping_method"] == "measured_species_rule_node" and
                          row["node_id"] in product_nodes}

    checks = {
        "all_rule_records_have_declared_lr_or_rl_direction": all(
            row["declared_direction"] in {"LR", "RL"} for row in rule_checks),
        "all_rule_records_are_projected_without_synthetic_reverse": all(
            row["irreversible_rule_projection"] for row in rule_checks),
        "all_rule_sources_are_substrates": all(
            row["source_replaced_on_substrate_side"] for row in rule_checks),
        "all_rule_targets_are_products": all(
            row["target_replaced_on_product_side"] for row in rule_checks),
        "all_catalyst_specific_rule_product_edges_are_active": all(
            row["all_product_edges_retained"] for row in rule_checks),
        "no_active_product_evidence_lost_to_signed_edge_conflicts": not inactive_products,
        "all_selected_topology_downstream_inputs_are_product_nodes": pruned_nodes <= product_nodes,
        "all_selected_solution_level0_endpoints_are_product_nodes": solution_nodes <= product_nodes,
    }
    if not all(checks.values()):
        raise AssertionError("Reaction-target audit check failed: " + json.dumps(checks, sort_keys=True))

    summary = {
        "scope": (
            "Directed product-side coverage for the frozen pkn_lipid_rules and selected "
            "Factor 1 TF-to-lipid MOON replay. A product candidate requires an active "
            "Enzyme__ -> Metab__ product_formation edge; compatible and hierarchy mappings are excluded."
        ),
        "definitions": {
            "species_rule_connected": "Feature appears on either substrate or product side of a measured-species rule.",
            "product_candidate": (
                "Feature is a species-rule product or has an exact SwissLipids-to-ChEBI "
                "candidate that is an active reaction product. It can still be excluded from "
                "the frozen input vector if its exact ChEBI identity is ambiguous."
            ),
            "accepted_product_endpoint": (
                "Frozen accepted input whose selected node has an active incoming "
                "product_formation edge."
            ),
        },
        "counts": {
            "all_measured_features": len(all_features),
            "historical_factor1_64_features": len(factor1_features),
            "rule_instances": len(pairs),
            "rule_catalyst_specific_product_edges": sum(row["catalyst_count"] for row in rule_checks),
            "rule_connected_features": len(source_features | rule_product_features),
            "rule_substrate_features": len(source_features),
            "rule_product_features": len(rule_product_features),
            "rule_substrate_only_features": len(source_features - rule_product_features),
            "rule_product_only_features": len(rule_product_features - source_features),
            "rule_both_side_features": len(source_features & rule_product_features),
            "all_panel_exact_chemical_product_candidate_features": len(exact_product_features),
            "all_panel_any_product_candidate_features": len(rule_product_features | exact_product_features),
            "all_panel_any_product_candidate_percent": round(
                100 * len(rule_product_features | exact_product_features) / len(all_features), 3),
            "frozen_accepted_input_nodes": len(included),
            "frozen_accepted_product_endpoints": len(accepted_product),
            "frozen_accepted_rule_product_endpoints": len(rule_input_product),
            "frozen_accepted_exact_chebi_product_endpoints": len(exact_input_product),
            "frozen_accepted_nonproduct_nodes": len(accepted_nonproduct),
            "frozen_accepted_rule_substrate_only_nodes": len(accepted_nonproduct & (source_features - rule_product_features)),
            "frozen_accepted_nonrule_nonproduct_nodes": len(accepted_nonproduct - source_features),
            "selected_topology_retained_downstream_nodes": len(pruned_nodes),
            "selected_topology_retained_product_nodes": len(pruned_nodes & product_nodes),
            "selected_solution_level0_endpoints": len(solution_nodes),
            "selected_solution_level0_product_endpoints": len(solution_nodes & product_nodes),
        },
        "notable_feature_sets": {
            "candidate_product_but_not_accepted_as_active_product_endpoint": sorted(candidate_only),
            "accepted_nonproduct_nodes": sorted(accepted_nonproduct),
            "historical_factor1_product_candidates": sorted(
                factor1_features & (rule_product_features | exact_product_features)),
            "selected_solution_product_features": sorted(solution_features),
        },
        "checks": checks,
        "population_summary": population_rows,
    }

    output.mkdir(parents=True, exist_ok=True)
    inventory_fields = [
        "feature_id", "in_lipid_de", "in_historical_factor1_64", "best_mapping_status",
        "nomenclature_level", "exact_chebi_ids", "species_rule_role", "rule_source_reaction_ids",
        "rule_product_reaction_ids", "rule_product_node", "rule_product_active_edges",
        "exact_chemical_product_nodes", "exact_chemical_product_active_edges", "has_product_candidate",
        "accepted_frozen_input", "accepted_node", "accepted_mapping_method",
        "accepted_node_is_active_product", "retained_after_selected_topology_pruning",
        "selected_solution_level0_endpoint",
    ]
    write_csv(output / "feature_reaction_target_inventory.csv", inventory, inventory_fields)
    write_csv(output / "rule_direction_validation.csv", rule_checks, list(rule_checks[0]))
    write_csv(output / "population_summary.csv", population_rows, list(population_rows[0]))
    write_csv(output / "dependency_hashes.csv", [
        {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for path in dependencies
    ], ["path", "sha256"])
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Measured lipid reaction-target audit\n\n"
        "This directory distinguishes a feature that is merely connected to a lipid rule "
        "from a feature that is an active product-side endpoint in the directed PKN. "
        "The audit does not modify the PKN or refit MOON.\n\n"
        "`feature_reaction_target_inventory.csv` contains one row per measured feature. "
        "`rule_direction_validation.csv` confirms the source/product substitution and active "
        "catalyst-specific product edges for every rule. `population_summary.csv` compares the "
        "complete panel, historical 64-lipid subset, frozen mapped inputs, selected topology and "
        "selected solution. `summary.json` contains definitions, headline counts and checks.\n"
    )
    print(json.dumps({"output": str(output), "counts": summary["counts"], "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
