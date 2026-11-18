#!/usr/bin/env python3
"""Prepare one-observation/one-node inputs for the expanded Factor 1 replay.

Read-only until --output is supplied. --validate-only never creates files.
The frozen 155 accepted identities are retained; only previously unmapped
features receive a new, explicit MeasuredLipid reaction-product identity.
This script does not score MOON, copy loadings across alternatives, assign
glycan identities, parse an RDS, or modify the PKN. Downstream R code should
take full-precision values from the frozen RDS by feature ID and check their
agreement with this CSV ledger (R's original CSV exports round doubles).
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
INPUTS = ROOT / "miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs"
EXPANDED = ROOT / "results/causal_modelling/cosmos_mouse_pkn_expanded_species"
SOURCES = {
    "frozen_mapping": INPUTS / "factor1_lipid_node_mapping.csv",
    "full_weights": INPUTS / "factor1_lipid_weights_all.csv",
    "factor_audit": INPUTS / "factor1_model_audit.json",
    "frozen_inputs_rds": INPUTS / "factor1_inputs.rds",
    "candidate_inventory": EXPANDED / "measured_product_candidate_inventory.csv",
    "pkn": EXPANDED / "pkn_expanded_chain.csv",
    "nodes": EXPANDED / "nodes.jsonl",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def keyed(records, field, description):
    result = {row[field]: row for row in records}
    require(len(result) == len(records), f"Duplicate {field} in {description}")
    require(all(result), f"Empty {field} in {description}")
    return result


def checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def field_value(value):
    if value is None:
        return ""
    if isinstance(value, list):
        return "|".join(str(item) for item in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def prepare():
    full_rows = rows(SOURCES["full_weights"])
    full = keyed(full_rows, "feature", "all Factor 1 lipid weights")
    require(len(full) == 1254, "Expected the frozen 1,254-feature lipid panel")
    require(all(row["factor"] == "Factor1" and row["view"] == "lipidomic"
                for row in full_rows), "Unexpected factor or omics view")

    frozen_rows = rows(SOURCES["frozen_mapping"])
    frozen = keyed(frozen_rows, "feature", "frozen lipid mapping")
    require(set(frozen) == set(full), "Frozen mapping and full weights differ")
    accepted = {feature: row for feature, row in frozen.items()
                if row["included"] == "TRUE"}
    require(len(accepted) == 155, "Expected 155 frozen accepted observations")
    require(len({row["node_id"] for row in accepted.values()}) == 155,
            "Frozen accepted node IDs are duplicated")

    factor_audit = json.loads(SOURCES["factor_audit"].read_text())
    multiplier = Decimal(str(factor_audit["saved_to_testes_positive_multiplier"]))
    require(multiplier == 1, "This replay requires the frozen +1 factor orientation")
    means = factor_audit["saved_coordinate_group_means"]
    require(means["Testes"] > means["Ovaries"], "Factor orientation audit disagrees")
    for feature, row in full.items():
        value = Decimal(row["factor1_weight"])
        require(value.is_finite() and value != 0, f"Invalid lipid weight: {feature}")
        require(value == Decimal(frozen[feature]["raw_weight"]),
                f"Full weight differs from frozen mapping: {feature}")
        require(value * multiplier == Decimal(frozen[feature]["oriented_weight"]),
                f"Oriented weight differs from frozen mapping: {feature}")

    inventory = keyed([row for row in rows(SOURCES["candidate_inventory"])
                       if row["overlay"] == "chain_compatible"],
                      "feature_id", "chain-compatible inventory")
    require(set(inventory) == set(full), "Candidate inventory has a different panel")
    candidates = {feature: row for feature, row in inventory.items()
                  if row["has_product_candidate"] == "True"}
    require(len(candidates) == 419, "Expected 419 chain-compatible product features")

    graph = rows(SOURCES["pkn"])
    require(len(graph) == 331032, "Expanded chain PKN differs from audited build")
    require(all(row["interaction"] in {"-1", "1"} for row in graph),
            "PKN has invalid signs")
    ids = {row[field] for row in graph for field in ("source", "target")}
    products = {row["target"] for row in graph
                if row["source"].startswith("Enzyme__")
                and row["target"].startswith("Metab__")}
    require(not any("CandidateMeasuredLipid:" in node for node in ids),
            "Conditional glycan nodes were found in the primary PKN")

    node_rows = []
    with SOURCES["nodes"].open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                node_rows.append(json.loads(line))
    all_nodes = keyed(node_rows, "node_id", "expanded node catalog")
    require(ids <= set(all_nodes), "Active PKN nodes missing from label catalog")
    active_nodes = [row for row in node_rows if row["node_id"] in ids]

    mapping = []
    for feature, full_row in full.items():
        candidate_nodes = candidates.get(feature, {}).get("product_nodes", "")
        if feature in accepted:
            node = accepted[feature]["node_id"]
            method = "preserve_frozen_accepted_node"
            # Retain the exact numeric text in the frozen ledger.
            weight = accepted[feature]["oriented_weight"]
        elif feature in candidates:
            node = "Metab__MeasuredLipid:" + feature + "_unknown"
            method = "new_one_to_one_measured_lipid_product"
            weight = str(Decimal(full_row["factor1_weight"]) * multiplier)
            require(node in candidate_nodes.split("|"),
                    f"Feature has no explicit measured product candidate: {feature}")
            require(node in products, f"Chosen new node is not a product: {node}")
            require(all_nodes[node].get("measurement_feature_id") == feature,
                    f"Node metadata has a different measurement identity: {node}")
        else:
            continue
        require(node in ids and node.startswith("Metab__"),
                f"Chosen observation absent or incorrectly typed: {node}")
        require(not all_nodes[node].get("candidate_identity_condition"),
                f"Conditional identity cannot be assigned in primary mapping: {node}")
        mapping.append({
            "node_id": node, "feature_id": feature, "weight": weight,
            "mapping_method": method, "is_active_product": node in products,
            "product_candidates": candidate_nodes,
        })

    require(len(mapping) == 423, "Expected 423 one-to-one observations")
    require(len({row["node_id"] for row in mapping}) == len(mapping),
            "Chosen node IDs are not unique")
    require(len({row["feature_id"] for row in mapping}) == len(mapping),
            "A feature loading was copied across nodes")
    require(sum(row["is_active_product"] for row in mapping) == 418,
            "Chosen active-product count differs from independent audit")
    methods = Counter(row["mapping_method"] for row in mapping)
    require(methods["preserve_frozen_accepted_node"] == 155 and
            methods["new_one_to_one_measured_lipid_product"] == 268,
            "Frozen/new mapping counts differ")

    dependencies = {str(path.relative_to(ROOT)): checksum(path)
                    for path in SOURCES.values()}
    dependencies[str(Path(__file__).resolve().relative_to(ROOT))] = checksum(Path(__file__))
    audit = {
        "analysis": "factor1_tf_to_lipid_expanded_chain",
        "mapping_policy": "Preserve all 155 frozen accepted nodes and values; add one explicit MeasuredLipid product node for each newly eligible feature.",
        "counts": {
            "full_weight_features": len(full), "all_weights_finite_nonzero": len(full),
            "frozen_observations_preserved": len(accepted),
            "original_active_product_features": 128,
            "chain_candidate_product_features": len(candidates),
            "newly_mapped_features": 268, "chosen_observations": len(mapping),
            "chosen_active_product_endpoints": 418, "chosen_nonproduct_observations": 5,
            "multi_product_alias_features": sum(int(row["product_node_count"]) > 1
                                                for row in candidates.values()),
            "mapping_methods": dict(methods), "active_PKN_nodes": len(ids),
        },
        "orientation": {
            "label": "testes_positive_Factor1_coordinate", "multiplier": int(multiplier),
            "saved_coordinate_group_means": means,
            "interpretation": "Common factor-axis orientation; not a direct testes-minus-ovaries differential effect.",
        },
        "eligibility_checks": {
            "frozen_mapping_and_full_weights_identical": True,
            "unique_feature_ids": True, "unique_node_ids": True,
            "all_chosen_ids_present": True, "new_ids_are_active_product_endpoints": True,
            "finite_nonzero_weights": True, "frozen_identities_and_CSV_values_retained": True,
            "no_feature_loading_copied_to_alternative_nodes": True,
            "conditional_glycan_identities_excluded": True,
        },
        "nonproduct_observations": [row for row in mapping if not row["is_active_product"]],
        "multi_product_alias_cases": [row for row in mapping
                                      if len(row["product_candidates"].split("|")) > 1],
        "assumptions": [
            "Original TF/RNA/regulon vectors and factor orientation are retained by the downstream R workflow; this script only prepares the lipid ledger.",
            "Each loading is assigned once. Alternative ChEBI/measured identities may remain unobserved intermediates in the PKN.",
            "LPC(20:0/0:0) retains frozen CHEBI:74968_unknown, although only its new MeasuredLipid alias is an active product; primary candidate coverage419 therefore yields418 chosen product endpoints.",
            "Four other original observations are substrate-only. They remain in the comparison input ledger and must pass the usual directed topology preparation.",
            "Unreported sn positions, DG occupancy, double-bond positions, stereochemistry and enzyme molecular specificity remain conditional reaction hypotheses.",
            "Bulk-liver lipid loadings and hepatocyte TF activities retain the cell-type scope of the frozen analysis.",
            "Reaction projection encodes catalyst/substrate participation and product formation; it is not a stoichiometric flux or signed concentration model.",
            "CSV weights retain frozen exported numeric text; read full-precision values from factor1_inputs.rds by feature ID when creating the R input vector and check agreement within export precision.",
        ],
        "sha256": dependencies,
        "output_files": ["input_mapping.csv", "nodes.csv", "mapping_audit.json"],
        "no_scoring_performed": True,
    }
    return mapping, active_nodes, audit


def write_csv(path, records, fields):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in records:
            writer.writerow({key: field_value(row.get(key)) for key in fields})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        help="New or empty output directory; existing results are never overwritten")
    parser.add_argument("--validate-only", action="store_true",
                        help="Validate and print the audit without writing files")
    args = parser.parse_args()
    if args.output is None and not args.validate_only:
        parser.error("Use --output=DIR or --validate-only")
    if args.output is not None and not args.validate_only:
        require(not args.output.exists() or
                (args.output.is_dir() and not any(args.output.iterdir())),
                "Output directory must be new or empty; refusing to overwrite")
    mapping, nodes, audit = prepare()
    if not args.validate_only:
        args.output.mkdir(parents=True, exist_ok=True)
        write_csv(args.output / "input_mapping.csv", mapping,
                  ["node_id", "feature_id", "weight", "mapping_method",
                   "is_active_product", "product_candidates"])
        base_fields = ["node_id", "label", "node_type", "gene_symbol", "chemical_id",
                       "compartment", "reaction_id", "resource", "catalysts"]
        fields = base_fields + sorted({key for row in nodes for key in row} - set(base_fields))
        write_csv(args.output / "nodes.csv", nodes, fields)
        (args.output / "mapping_audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"validated": True, "validate_only": args.validate_only,
                      "output": None if args.validate_only else str(args.output),
                      "counts": audit["counts"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
