#!/usr/bin/env python3
"""Declare currency-metabolite substrate connectors excluded in MOON sensitivities.

Uses frozen node labels and exact signed projection provenance. Full reactions,
all source PKNs, and chemical-to-protein regulatory edges remain unchanged.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

from export_expanded_factor1_network import CURRENCY_LABEL, currency_node, edge_key, jsonlines, read_csv, sha256, write_csv

ROOT = Path(__file__).resolve().parents[3]
MINIMAL_WATER = "Metab__CHEBI:15377_unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pkn", type=Path, default=ROOT / "results/causal_modelling/cosmos_mouse_pkn_expanded_species")
    parser.add_argument("--baseline-pkn", type=Path, default=ROOT / "results/causal_modelling/cosmos_mouse_pkn")
    parser.add_argument("--pkn-file", default="pkn_expanded_chain.csv")
    parser.add_argument("--alias-complete", action="store_true",
                        help="Version 2: additionally recognize Mouse-GEM H+, Pi and PPi labels; default reproduces v1")
    parser.add_argument("--node-output", type=Path, default=ROOT / "support/cosmos_mouse/moon_currency_traversal_exclusions.csv")
    parser.add_argument("--edge-output", type=Path, default=ROOT / "support/cosmos_mouse/moon_currency_substrate_edges.csv")
    parser.add_argument("--manifest-output", type=Path, default=ROOT / "support/cosmos_mouse/moon_currency_traversal_policy_manifest.json")
    args = parser.parse_args()
    outputs = [args.node_output, args.edge_output, args.manifest_output]
    if any(p.exists() for p in outputs):
        parser.error("Use new output paths; policy snapshots are not overwritten")
    sources = [args.pkn / "nodes.jsonl", args.pkn / args.pkn_file,
               args.baseline_pkn / "edge_evidence.jsonl.gz", args.pkn / "extension_edge_evidence.jsonl.gz",
               Path(__file__).resolve(), Path(__file__).with_name("export_expanded_factor1_network.py")]
    for path in sources:
        if not path.is_file():
            parser.error("Missing dependency: " + str(path))
    metadata = {r["node_id"]: r for r in jsonlines(sources[0])}
    graph = {edge_key(r) for r in read_csv(sources[1])}
    active = {n for s, _, t in graph for n in (s, t)}
    currency = {n for n in active if currency_node(metadata[n])} if args.alias_complete else {
        n for n in active if metadata[n]["node_type"] == "metabolite" and
        CURRENCY_LABEL.fullmatch(metadata[n].get("label", "")) and
        metadata[n].get("label", "").lower() not in {"h+", "pi", "ppi"}}
    edge_evidence, preserved_regulation = defaultdict(list), []
    for path in sources[2:4]:
        for row in jsonlines(path):
            eligible = (row.get("first_variant") in {"core", "lipid_curated", "lipid_rules"}
                        if path == sources[2] else row.get("overlay_tier") == "chain_compatible")
            key = edge_key(row)
            if not eligible or key not in graph or row["source"] not in currency:
                continue
            if (row.get("mechanism") == "substrate_participation" and
                    metadata[row["target"]]["node_type"] == "enzyme_reaction_instance"):
                edge_evidence[key].append(row)
            elif metadata[row["target"]]["node_type"] in {"gene_product", "protein_unmapped"}:
                preserved_regulation.append(row)
    if not edge_evidence or MINIMAL_WATER not in currency:
        raise ValueError("Expected expanded-PKN currency/water connectors are absent")
    degree = Counter(key[0] for key in edge_evidence)
    rows = []
    for n in sorted(currency):
        m = metadata[n]
        rows.append(dict(node_id=n, label=m["label"], chemical_id=m.get("chemical_id", ""),
            compartment=m.get("compartment", ""), minimal_water_bridge=n == MINIMAL_WATER,
            conventional_currency=True, substrate_connector_edges=degree[n],
            reason="Shared reaction currency can connect unrelated reactions when represented as independent positive substrate edges",
            mechanism_scope="substrate_participation to enzyme_reaction_instance only",
            evidence_basis="Frozen assembled node label/chemical ID and exact source/interaction/target mechanism provenance; modeling sensitivity, not a claim of biological irrelevance"))
    edge_rows = []
    for (s, i, t), ee in sorted(edge_evidence.items()):
        edge_rows.append(dict(source=s, interaction=i, target=t, source_label=metadata[s]["label"],
            target_label=metadata[t]["label"], minimal_water_bridge=s == MINIMAL_WATER,
            conventional_currency=True, mechanism_scope="substrate_participation",
            resources="|".join(sorted({r["resource"] for r in ee})),
            evidence_records=len(ee),
            evidence_basis="Exact frozen signed PKN connector with substrate_participation provenance and enzyme_reaction_instance target"))
    if not all(edge_key(r) in graph and r["interaction"] == 1 and
               metadata[r["target"]]["node_type"] == "enzyme_reaction_instance" for r in edge_rows):
        raise ValueError("Policy includes an unsupported/non-substrate connector")
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    write_csv(args.node_output, rows, list(rows[0]))
    write_csv(args.edge_output, edge_rows, list(edge_rows[0]))
    summary = dict(policy_version=2 if args.alias_complete else 1,
        alias_complete=args.alias_complete, currency_nodes=len(rows), excluded_signed_substrate_edges=len(edge_rows),
        minimal_water_unknown_substrate_edges=sum(r["minimal_water_bridge"] for r in edge_rows),
        currency_label_counts=dict(Counter(r["label"] for r in rows)),
        excluded_edge_source_label_counts=dict(Counter(r["source_label"] for r in edge_rows)),
        substrate_evidence_records=sum(r["evidence_records"] for r in edge_rows),
        chemical_to_protein_regulatory_evidence_preserved=len(preserved_regulation),
        policy="Remove only listed exact positive currency->enzyme_reaction_instance substrate_participation connectors in a separate scoring-PKN sensitivity. Retain incoming product edges, full reaction participants and all chemical->protein regulatory edges.",
        minimal_policy="Only Metab__CHEBI:15377_unknown outgoing substrate connectors, the actual water bridge found in the initial selected expanded fit.",
        conventional_policy="All declared currency nodes across every retained compartment, selected by the explicit frozen label panel; not hypothetical acyl pools, fatty acids or measured lipids.",
        source_pkn_modified=False, source_reaction_records_modified=False,
        source_hashes={str(p.resolve()): sha256(p) for p in sources},
        output_hashes={str(p.resolve()): sha256(p) for p in outputs[:2]})
    args.manifest_output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
