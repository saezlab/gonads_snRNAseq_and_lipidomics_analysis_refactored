#!/usr/bin/env python3
"""Compare the deposited La Manno lipid modules with our frozen mouse PKN.

This is a read-only resource/mapping audit, not a PKN builder or MOON run.
All chemistry matches remain nomenclature candidates. Summed compositions do
not resolve acyl-chain identities, enzyme specificity, or compartment.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ROOT / "support/cosmos_mouse/comparison/lamanno_2026"
PKN = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
OUT = ROOT / "results/causal_modelling/lamanno_2026_pkn_comparison"


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def checksum(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def write_csv(path, values, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(values)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def truth(value):
    return str(value).lower() == "true"


def enzyme_symbols(value):
    # Parse a literal tuple, never execute its text.
    if value.startswith("("):
        symbols = ast.literal_eval(value)
        if not isinstance(symbols, tuple) or not all(isinstance(x, str) for x in symbols):
            raise ValueError("Unexpected enzyme annotation: " + value)
        return symbols
    return (value,)


def canonical_shorthand(value):
    match = re.fullmatch(r"([^ ()]+)\((.+)\)", value)
    return f"{match[1]} {match[2]}" if match else value


def normalized_original_label(value):
    # Only notation/class synonyms; no chain reordering or oxidation erasure.
    aliases = {"TAG": "TG", "DAG": "DG", "PE-Cer": "EPC",
               "SulfoHexCer": "SHexCer", "GM3": "Hex2NeuAcCer"}
    value = canonical_shorthand(value.strip())
    lipid_class, separator, body = value.partition(" ")
    return aliases.get(lipid_class, lipid_class) + separator + body


def draft_mapping(features, graph_nodes, graph_edges):
    by_label = {normalized_original_label(row["node"]): row for row in graph_nodes}
    assert len(by_label) == len(graph_nodes), "Normalized draft labels collided"
    sources = {normalized_original_label(row["source"]) for row in graph_edges}
    products = {normalized_original_label(row["target"]) for row in graph_edges}
    assert all("_" not in label and "/" not in label for label in by_label)
    assignments = []
    for feature in features:
        original = normalized_original_label(feature["feature_id"])
        canonical = feature["canonical_name"]
        candidates = {}
        if original in by_label:
            candidates[original] = "identical_normalized_original_shorthand"
        projected = sum_composition(canonical)
        if projected in by_label and projected not in candidates:
            candidates[projected] = ("conditional_sum_composition_projection" if projected != canonical
                                     else "canonical_parser_equivalent_shorthand")
        for label, match_type in candidates.items():
            node = by_label[label]
            assignments.append(dict(
                feature_id=feature["feature_id"],
                in_historical_factor1_64=truth(feature["in_factor1"]),
                measured_original=original, measured_canonical=canonical,
                graph_node=node["node"], match_type=match_type,
                graph_node_is_source=label in sources, graph_node_is_target=label in products,
                graph_node_active=label in sources | products,
                graph_node_inferred=node["inferred"],
                representative_molecular_string_not_used_for_mapping=node["molecular_string"]
            ))
    def selected(predicate=lambda row: True):
        return {row["feature_id"] for row in assignments if predicate(row)}
    exact = selected(lambda row: row["match_type"] == "identical_normalized_original_shorthand")
    historical = {row["feature_id"] for row in features if truth(row["in_factor1"])}
    return assignments, dict(
        active_nodes=len(sources | products), source_nodes=len(sources), product_nodes=len(products),
        identical_normalized_original_shorthand_features=len(exact),
        additional_conditional_sum_projection_features=len(selected(lambda row: row["match_type"] == "conditional_sum_composition_projection")),
        additional_canonical_parser_equivalent_features=len(selected(lambda row: row["match_type"] == "canonical_parser_equivalent_shorthand")),
        possible_feature_matches=len(selected()), matched_graph_nodes=len({row["graph_node"] for row in assignments}),
        identical_feature_matches_to_graph_products=len(selected(lambda row: row["graph_node_is_target"] and row["match_type"] == "identical_normalized_original_shorthand")),
        possible_feature_matches_to_graph_products=len(selected(lambda row: row["graph_node_is_target"])),
        possible_feature_matches_to_graph_sources=len(selected(lambda row: row["graph_node_is_source"])),
        identical_historical_factor1_64_matches=len(exact & historical),
        possible_historical_factor1_64_matches=len(selected() & historical),
        maximum_measurements_per_graph_node=max(Counter(row["graph_node"] for row in assignments).values()),
        match_type_counts=dict(Counter(row["match_type"] for row in assignments))
    )


def sum_composition(value):
    """Sum only supported parsed shorthand; retain oxygen and ether chemistry."""
    value = canonical_shorthand(value)
    if " " not in value:
        return None
    lipid_class, body = value.split(" ", 1)
    carbon = db = oxygen = ethers = 0
    for chain in re.split(r"[/_]", body):
        match = re.fullmatch(r"(O-|P-)?(\d+):(\d+)(?:;O(\d*))?", chain)
        if not match:
            return None
        carbon += int(match[2])
        db += int(match[3])
        if ";O" in chain:
            oxygen += int(match[4] or 1)
        if match[1]:
            ethers += 1
            if match[1] == "P-":
                db += 1
    if ethers > 1:
        return None
    oxygen_suffix = ";O" + (str(oxygen) if oxygen > 1 else "") if oxygen else ""
    return f"{lipid_class} " + ("O-" if ethers else "") + f"{carbon}:{db}" + oxygen_suffix


def self_check():
    assert sum_composition("PC 20:0_22:6") == "PC 42:6"
    assert sum_composition("Cer 18:1;O2/18:0") == "Cer 36:1;O2"
    assert sum_composition("PC O-16:0/18:1") == "PC O-34:1"
    assert sum_composition("PC P-16:0/18:1") == "PC O-34:2"
    assert sum_composition("PC O-18:0/18:1") != "PC 36:1"
    assert sum_composition("PC unknown") is None
    assert canonical_shorthand("LPC(16:0)") == "LPC 16:0"
    assert normalized_original_label("DAG(16:0_18:1)") == "DG 16:0_18:1"
    assert "16:0" not in re.split(r"[/_]", "18:0_18:1")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, default=SOURCES)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    self_check()
    if args.self_test:
        print("PASS: summed-composition parsing preserves ether/oxygen identity and chain compatibility.")
        return
    sources = args.sources.resolve()
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Choose a new/empty --output directory: " + str(out))
    manifest = json.loads((sources / "source_manifest.json").read_text())
    for item in manifest["files"]:
        path = sources / Path(item["path"]).name
        if checksum(path) != item["sha256"]:
            raise SystemExit("Comparison source changed: " + str(path))

    table = rows(sources / "LaManno_SuppTable10.csv")
    features = rows(PKN / "mapping_features.csv")
    rules = rows(PKN / "lipid_rules_pairs.csv")
    table_nodes = {row[key] for row in table for key in ("reagent", "product")}
    assignments = []
    for feature in features:
        label = sum_composition(feature["canonical_name"])
        if label in table_nodes:
            assignments.append(dict(
                feature_id=feature["feature_id"],
                in_historical_factor1_64=truth(feature["in_factor1"]),
                nomenclature_level=feature["nomenclature_level"],
                canonical_name=feature["canonical_name"], table_node=label,
                mapping_type="identical_reported_shorthand" if label == feature["canonical_name"]
                else "conditional_sum_composition_projection"
            ))
    by_node = defaultdict(list)
    for row in assignments:
        by_node[row["table_node"]].append(row)
    assignment_summary = [dict(
        table_node=node, measured_features=len(values),
        historical_factor1_64_features=sum(row["in_historical_factor1_64"] for row in values),
        feature_ids="|".join(row["feature_id"] for row in values)
    ) for node, values in sorted(by_node.items())]

    # A table row valid for a sum species need not be valid for every chain-resolved
    # measurement with that sum. Audit the retained lyso chain explicitly.
    lyso = {"LPC": "PC", "LPE": "PE", "LPG": "PG", "LPA": "PA", "LPS": "PS"}
    chain_checks = []
    for index, reaction in enumerate(table, 2):
        for direction, side in [("hydrolysis", "reagent"), ("acylation", "product")]:
            large = reaction[side]
            small = reaction["product" if side == "reagent" else "reagent"]
            small_class, small_chain = small.split(" ", 1)
            if small_class not in lyso or large.split()[0] != lyso[small_class]:
                continue
            for assignment in by_node.get(large, []):
                body = assignment["canonical_name"].split(" ", 1)[1]
                if not re.search(r"[/_]", body):
                    continue
                chains = re.split(r"[/_]", body)
                if not all(re.fullmatch(r"(O-|P-)?\d+:\d+(;O\d*)?", item) for item in chains):
                    continue
                chain_checks.append(dict(
                    table_csv_row=index, feature_id=assignment["feature_id"],
                    measured_canonical=assignment["canonical_name"], direction=direction,
                    table_reagent=reaction["reagent"], table_product=reaction["product"],
                    table_enzyme=reaction["enzyme"], retained_lyso_chain_present=small_chain in chains
                ))

    active = {row[key] for row in rows(PKN / "pkn_lipid_rules.csv") for key in ("source", "target")}
    nodes = [json.loads(line) for line in (PKN / "nodes.jsonl").open()]
    enzymes = [node for node in nodes if node["node_type"] == "enzyme_reaction_instance"
               and node["node_id"] in active]
    genes = sorted({gene for row in table for gene in enzyme_symbols(row["enzyme"])})
    enzyme_inventory = []
    for gene in genes:
        matches = [node for node in enzymes if "Gene__" + gene in node.get("catalysts", [])]
        enzyme_inventory.append(dict(
            table_enzyme=gene, exact_gene_node_in_our_pkn="Gene__" + gene in active,
            active_enzyme_reaction_instances=len(matches),
            active_rhea_instances=sum(node.get("resource") == "Rhea" for node in matches),
            active_mouse_gem_instances=sum(node.get("resource") == "Mouse-GEM" for node in matches),
            active_measured_species_rule_instances=sum(node.get("resource") == "LINEX2_rule" for node in matches),
            alias_review="Pisd is present; reconcile Psd annotation before import" if gene == "Psd" else ""
        ))

    product_rows = rows(ROOT / "results/causal_modelling/lipid_reaction_target_mapping_audit/feature_reaction_target_inventory.csv")
    accepted_products = {row["feature_id"] for row in product_rows if truth(row["accepted_node_is_active_product"])}
    table_features = {row["feature_id"] for row in assignments}
    rule_features = {row[key] for row in rules for key in ("source_feature", "target_feature")}
    canonical = {row["feature_id"]: row["canonical_name"] for row in features}
    projected_rule_pairs = {(sum_composition(canonical[row["source_feature"]]),
                             sum_composition(canonical[row["target_feature"]])) for row in rules}
    table_pairs = {(row["reagent"], row["product"]) for row in table}
    table_transforms = Counter(row["reagent"].split()[0] + " -> " + row["product"].split()[0] for row in table)
    rule_transforms = Counter(row["headgroup_transform"] for row in rules)

    graph_nodes = rows(sources / "graph_from_annotations_nodes.csv")
    graph_edges = rows(sources / "graph_from_annotations_edges.csv")
    draft = rows(sources / "reactions_linex_now.csv")
    assert [(row["source"], row["target"]) for row in graph_edges] == [(row["reagent"], row["product"]) for row in draft]
    draft_assignments, draft_mapping_summary = draft_mapping(features, graph_nodes, graph_edges)
    endpoint_overlap = [dict(
        feature_id=feature, possible_table10_nomenclature_match=feature in table_features,
        our_accepted_active_product=feature in accepted_products
    ) for feature in sorted(table_features | accepted_products)]

    summary = dict(
        scope="Deposited Table10 lipid component and older annotation-derived lipid graph; full signed CORNETO PKN was not located in the cited deposits.",
        table10=dict(
            reaction_rows=len(table), distinct_directed_pairs=len(table_pairs), lipid_nodes=len(table_nodes),
            lipid_classes=len({node.split()[0] for node in table_nodes}), enzyme_symbols=len(genes),
            row_enzyme_triplets=sum(len(enzyme_symbols(row["enzyme"])) for row in table),
            class_transformations=len(table_transforms), sha256=checksum(sources / "LaManno_SuppTable10.csv")
        ),
        mapping=dict(
            measured_features=len(features), historical_factor1_64_features=sum(truth(row["in_factor1"]) for row in features),
            identical_reported_shorthand=sum(row["mapping_type"] == "identical_reported_shorthand" for row in assignments),
            conditional_sum_projections=sum(row["mapping_type"] != "identical_reported_shorthand" for row in assignments),
            possible_table_matches=len(assignments), matched_table_nodes=len(by_node),
            maximum_measurements_per_table_node=max(map(len, by_node.values())),
            historical_factor1_64_matches=sum(row["in_historical_factor1_64"] for row in assignments),
            chain_audit_assignments=len(chain_checks), retained_lyso_chain_absent=sum(not row["retained_lyso_chain_present"] for row in chain_checks),
            affected_resolved_features=len({row["feature_id"] for row in chain_checks if not row["retained_lyso_chain_present"]}),
            note="These are reported-nomenclature overlaps; they are not chemical identity matches, validated reaction endpoints, or a proposed input vector. Every Factor1 lipid has a weight; the 64 flag is the historical correlation subset."
        ),
        our_pkn=dict(
            build_counts=json.loads((PKN / "build_summary.json").read_text())["variants"],
            rhea_counts={key: json.loads((PKN / "rhea_summary.json").read_text())[key]
                         for key in ("reaction_variants", "reaction_masters", "catalysts")},
            measured_rule_instances=len(rules), measured_rule_features=len(rule_features),
            instantiated_lipid_classes=len({canonical[feature].split()[0] for feature in rule_features}),
            class_transformations=len(rule_transforms),
            table_enzyme_exact_symbols_present=sum(row["exact_gene_node_in_our_pkn"] for row in enzyme_inventory),
            table_enzyme_symbols_with_active_rhea_instances=sum(row["active_rhea_instances"] > 0 for row in enzyme_inventory),
            table_enzyme_symbols_with_measured_rule_instances=sum(row["active_measured_species_rule_instances"] > 0 for row in enzyme_inventory),
            table_pairs_overlapping_summed_projection_of_our_rules=len(table_pairs & projected_rule_pairs),
            accepted_product_endpoints=len(accepted_products),
            table_matches_already_accepted_product_endpoints=len(table_features & accepted_products),
            table_match_candidates_outside_accepted_product_endpoints=len(table_features - accepted_products),
            accepted_product_endpoints_outside_table_matches=len(accepted_products - table_features),
            union_table_matches_and_accepted_products=len(table_features | accepted_products),
            overlap_note="A table sum-node match and an accepted chain-resolved product endpoint have different evidence/resolution; this is descriptive set overlap only."
        ),
        deposited_draft=dict(
            lipid_nodes=len(graph_nodes), edges=len(graph_edges),
            pairs_equal_reactions_linex_csv=True,
            mapping=draft_mapping_summary,
            rhea_bearing_edges=sum("RHEA:" in row["enzyme_id"] for row in graph_edges),
            node_types=dict(Counter(row["node_molecule_type"] for row in graph_nodes)),
            note="Lipid-only annotation-derived draft; pooled enzyme/Rhea/Reactome annotations and inferred representative chain strings do not establish mouse species specificity or a signed signaling network."
        )
    )
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "summary.json", summary)
    write_csv(out / "table10_measured_mapping_candidates.csv", assignments, list(assignments[0]))
    write_csv(out / "table10_node_measurement_multiplicity.csv", assignment_summary, list(assignment_summary[0]))
    write_csv(out / "table10_chain_compatibility.csv", chain_checks, list(chain_checks[0]))
    write_csv(out / "table10_enzyme_coverage.csv", enzyme_inventory, list(enzyme_inventory[0]))
    write_csv(out / "draft_measured_mapping_candidates.csv", draft_assignments, list(draft_assignments[0]))
    write_csv(out / "table10_vs_accepted_product_endpoints.csv", endpoint_overlap, list(endpoint_overlap[0]))
    transform_rows = [dict(transform=key, table10_rows=table_transforms[key], our_species_rule_instances=rule_transforms[key])
                      for key in sorted(set(table_transforms) | set(rule_transforms))]
    write_csv(out / "reaction_family_comparison.csv", transform_rows, list(transform_rows[0]))
    dependencies = [Path(__file__), sources / "source_manifest.json", *[sources / Path(item["path"]).name for item in manifest["files"]],
                    *[PKN / name for name in ("mapping_features.csv", "lipid_rules_pairs.csv", "pkn_lipid_rules.csv", "nodes.jsonl", "build_summary.json", "rhea_summary.json")],
                    ROOT / "results/causal_modelling/lipid_reaction_target_mapping_audit/feature_reaction_target_inventory.csv"]
    write_csv(out / "dependency_hashes.csv", [dict(path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                                               sha256=checksum(path)) for path in dependencies], ["path", "sha256"])
    print(f"Comparison saved: {len(table)} Table10 rows; {len(assignments)}/{len(features)} possible nomenclature matches; "
          f"{sum(row['exact_gene_node_in_our_pkn'] for row in enzyme_inventory)}/{len(genes)} exact enzyme symbols already in our PKN. {out}")


if __name__ == "__main__":
    main()
