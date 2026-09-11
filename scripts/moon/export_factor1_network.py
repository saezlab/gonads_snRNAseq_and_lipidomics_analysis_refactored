#!/usr/bin/env python3
"""Join MOON solution nodes/edges to PKN provenance and audit representative paths.

All joins use untouched PKN IDs. One path is reported for every reachable upstream
TF/downstream observation pair: minimum edge count, then greatest minimum absolute
inferred MOON score, with deterministic node-ID tie breaking. Level-0 weights do
not enter that minimum because they are not on the inferred-score scale.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite(value, default=None):
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (ValueError, TypeError):
        return default


def sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def edge_key(row):
    return row["source"], int(float(row["interaction"])), row["target"]


def reachable(seeds, adjacency):
    visited = set(seeds)
    queue = collections.deque(seeds)
    while queue:
        for neighbor in adjacency.get(queue.popleft(), []):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)
    return visited


def write_graphml(path, nodes, edges, description):
    """Cytoscape-ready XML; numeric measurement and inference scales stay separate."""
    ns = "http://graphml.graphdrawing.org/xmlns"
    ET.register_namespace("", ns)
    root = ET.Element("{" + ns + "}graphml")
    node_fields = list(dict.fromkeys(key for row in nodes for key in row))
    node_fields += [field for field in ["factor1_weight", "moon_score"] if field not in node_fields]
    edge_fields = list(dict.fromkeys(key for row in edges for key in row))
    keys = {}
    for domain, fields in [("node", node_fields), ("edge", edge_fields), ("graph", ["description"])]:
        for field in fields:
            ident = "k" + str(len(keys))
            typ = "double" if field in {"score", "moon_score", "factor1_weight", "upstream_activity", "RNA_input"} else "int" if field in {"level", "interaction", "evidence_records"} else "string"
            keys[domain, field] = ident, typ
            ET.SubElement(root, "{" + ns + "}key", {"id": ident, "for": domain, "attr.name": field, "attr.type": typ})
    graph = ET.SubElement(root, "{" + ns + "}graph", {"id": "Factor1_MOON", "edgedefault": "directed"})

    def add_data(element, domain, row):
        for field, value in row.items():
            if value in {None, "", "NA", "NaN"}:
                continue
            key, typ = keys[domain, field]
            if typ in {"double", "int"}:
                numeric = finite(value)
                if numeric is None:
                    continue
                value = str(int(numeric)) if typ == "int" else str(numeric)
            ET.SubElement(element, "{" + ns + "}data", {"key": key}).text = str(value)

    add_data(graph, "graph", {"description": description})
    ids = {row["node_id"]: "n" + str(index) for index, row in enumerate(nodes)}
    for row in nodes:
        item = ET.SubElement(graph, "{" + ns + "}node", {"id": ids[row["node_id"]]})
        merged = dict(row, factor1_weight=row.get("score") if finite(row.get("level")) == 0 else "",
                      moon_score=row.get("score") if finite(row.get("level")) != 0 else "")
        add_data(item, "node", merged)
    for index, row in enumerate(edges):
        item = ET.SubElement(graph, "{" + ns + "}edge", {"id": "e" + str(index), "source": ids[row["source"]], "target": ids[row["target"]]})
        add_data(item, "edge", row)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def path_tree(seed, adjacency, attributes):
    """Shortest-path DAG dynamic program, maximizing the inferred-score bottleneck."""
    distance, parents, bottleneck = {seed: 0}, {}, {}
    bottleneck[seed] = abs(finite(attributes[seed].get("score"), 0))
    frontier = [seed]
    while frontier:
        following = set()
        for source in sorted(frontier):
            for target in adjacency.get(source, []):
                candidate_distance = distance[source] + 1
                if target in distance and distance[target] < candidate_distance:
                    continue
                target_score = math.inf if finite(attributes[target].get("level")) == 0 else abs(finite(attributes[target].get("score"), 0))
                candidate_bottleneck = min(bottleneck[source], target_score)
                if target not in distance or candidate_bottleneck > bottleneck[target] or (candidate_bottleneck == bottleneck[target] and source < parents[target]):
                    distance[target] = candidate_distance
                    bottleneck[target] = candidate_bottleneck
                    parents[target] = source
                following.add(target)
        frontier = sorted(following)
    return distance, parents, bottleneck


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path("results/cosmos_factor1_pilot"))
    parser.add_argument("--pkn", type=Path, default=Path("results/cosmos_mouse_pkn"))
    parser.add_argument("--input-mapping", type=Path, default=None,
                        help="Optional exact Factor 1 lipid-to-PKN mapping CSV; defaults to the frozen pilot input mapping when available.")
    parser.add_argument("--stem", default="solution_primary")
    parser.add_argument("--max-path-edges", type=int, default=None, help="Optional cap on the reported shortest routes (does not change scores or the supplied solution)")
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()
    output, pkn = args.results, args.pkn
    edges = read_csv(output / (args.stem + "_SIF.csv"))
    raw_attributes = read_csv(output / (args.stem + "_ATT.csv"))
    attributes = {row.get("source", row.get("node_id", row.get("nodes"))): row for row in raw_attributes}
    if len(attributes) != len(raw_attributes):
        raise ValueError("Duplicate node attributes")
    with (pkn / "nodes.jsonl").open() as handle:
        metadata = {row["node_id"]: row for row in map(json.loads, handle)}
    keys = {edge_key(row) for row in edges}
    if len(keys) != len(edges):
        raise ValueError("Duplicate edges in solution")
    endpoints = {node for source, _, target in keys for node in (source, target)}
    if endpoints != set(attributes):
        raise ValueError("SIF/ATT node sets differ")
    original_keys = {edge_key(row) for row in read_csv(pkn / "pkn_lipid_rules.csv")}
    if not keys <= original_keys:
        raise ValueError("Solution contains edges absent from the requested PKN")
    observed_features = collections.defaultdict(set)
    mapping_candidates = ([args.input_mapping] if args.input_mapping is not None else []) + [
        output / "inputs/factor1_lipid_node_mapping.csv",
        output.parent / "inputs/factor1_lipid_node_mapping.csv",
        Path("results/cosmos_factor1_pilot/inputs/factor1_lipid_node_mapping.csv"),
    ]
    mapping_path = next((path for path in mapping_candidates if path is not None and path.exists()), None)
    if mapping_path.exists():
        for row in read_csv(mapping_path):
            if row.get("included", "").lower() == "true" and row.get("node_id") in endpoints:
                observed_features[row["node_id"]].add(row["feature"])
    support = collections.defaultdict(list)
    with gzip.open(pkn / "edge_evidence.jsonl.gz", "rt") as handle:
        for line in handle:
            row = json.loads(line)
            if row.get("first_variant") in {"core", "lipid_curated", "lipid_rules"} and edge_key(row) in keys:
                support[edge_key(row)].append(row)
    if set(support) != keys:
        raise ValueError("Some exact solution edges have no eligible PKN provenance")
    evidence_stem = "network_evidence" if args.stem == "solution_primary" else "network_evidence_" + args.stem
    evidence_rows = []
    with gzip.open(output / (evidence_stem + ".jsonl.gz"), "wt") as handle:
        for key in sorted(keys):
            for item in support[key]:
                handle.write(json.dumps(item, sort_keys=True) + "\n")
            source, interaction, target = key
            entries = support[key]
            evidence_rows.append({
                "source": source, "interaction": interaction, "target": target,
                "source_label": metadata.get(source, {}).get("label", source),
                "target_label": metadata.get(target, {}).get("label", target),
                "mechanisms": "|".join(sorted({x.get("mechanism", "") for x in entries})),
                "resources": "|".join(sorted({x.get("resource", "") for x in entries})),
                "reaction_ids": "|".join(sorted({x["reaction_id"] for x in entries if x.get("reaction_id")})),
                "evidence_records": len(entries),
                "contains_species_rule_hypothesis": any(x.get("resource") == "LINEX2_rule" for x in entries),
                "references": "|".join(sorted({x.get("references", "") for x in entries if isinstance(x.get("references", ""), str)})),
            })
    edge_fields = ["source", "interaction", "target", "source_label", "target_label", "mechanisms", "resources", "reaction_ids", "evidence_records", "contains_species_rule_hypothesis", "references"]
    write_csv(output / (evidence_stem + ".csv"), evidence_rows, edge_fields)

    def label(node):
        meta = metadata.get(node, {})
        if observed_features.get(node):
            return "/".join(sorted(observed_features[node]))
        if meta.get("node_type") == "enzyme_reaction_instance":
            catalysts = "/".join(x.replace("Gene__", "") for x in meta.get("catalysts", []))
            rid = meta.get("reaction_id", "reaction")
            if rid.startswith("LIPID_RULE:"):
                rid = "species rule " + rid.split(":")[1][:6]
            return catalysts + " [" + rid + "]"
        return meta.get("label", node)

    all_nodes = []
    for node, row in sorted(attributes.items()):
        merged = dict(row)
        merged.update({"node_id": node, "label": label(node), "pkn_label": metadata.get(node, {}).get("label", ""),
                       "observed_features": "|".join(sorted(observed_features.get(node, []))),
                       "node_type": metadata.get(node, {}).get("node_type", "unknown"),
                       "reaction_id": metadata.get(node, {}).get("reaction_id", ""),
                       "catalysts": "|".join(metadata.get(node, {}).get("catalysts", []))})
        all_nodes.append(merged)
    node_fields = list(dict.fromkeys(["node_id", "label", "pkn_label", "observed_features", "node_type", "reaction_id", "catalysts"] + [key for row in raw_attributes for key in row]))
    write_csv(output / (evidence_stem + "_nodes.csv"), all_nodes, node_fields)
    write_graphml(output / (args.stem + ".graphml"), all_nodes, evidence_rows,
                  "Factor 1 MOON solution; exact PKN IDs and edge provenance retained. Factor 1 weights (level 0), inferred MOON scores and upstream TF activities are separate node attributes. Scoring result directory: " + str(output) + ". " +
                  ("This is the union of <=6-edge catalyst-only routes in the relaxed solution, retaining its fixed scores; no model was refitted." if args.stem == "solution_catalyst_bounded6" else "Extraction: " + args.stem))

    reaction_ids = {metadata.get(node, {}).get("reaction_id") for node in endpoints} - {None, ""}
    reaction_records = []
    for name in ["mouse_gem_reactions.jsonl", "rhea_reactions.jsonl", "lipid_rules_reactions.jsonl"]:
        with (pkn / name).open() as handle:
            for line in handle:
                row = json.loads(line)
                if row["reaction_id"] in reaction_ids:
                    reaction_records.append(row)
    with gzip.open(output / (evidence_stem + "_reactions.jsonl.gz"), "wt") as handle:
        for row in reaction_records:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    reaction_table = []
    for row in reaction_records:
        reaction_nodes = [node for node in sorted(endpoints) if metadata.get(node, {}).get("reaction_id") == row["reaction_id"]]
        catalysts = sorted({catalyst for node in reaction_nodes for catalyst in metadata[node].get("catalysts", [])})
        reaction_table.append({
            "reaction_id": row["reaction_id"], "resource": row.get("resource", ""),
            "enzyme_instance_ids": "|".join(reaction_nodes), "catalyst_node_ids": "|".join(catalysts),
            "hypothesis_type": row.get("hypothesis_type", "curated/model reaction annotation; qualitative graph projection"),
            "substrate_names": "|".join(item.get("name", item.get("id", "")) for item in row.get("substrates", [])),
            "product_names": "|".join(item.get("name", item.get("id", "")) for item in row.get("products", [])),
            "source_equation": row.get("source_template_equation") or row.get("equation", ""),
            "matched_feature_ids": "|".join(row.get("matched_feature_ids", [])),
            "supporting_reaction_id": row.get("supporting_reaction_id", ""),
            "assumptions": "|".join(row.get("assumptions", [])),
            "reference_scope": row.get("reference_scope", "Source annotation; organism, direction and stoichiometry are retained in the companion JSONL"),
        })
    write_csv(output / (evidence_stem + "_reactions.csv"), reaction_table,
              ["reaction_id", "resource", "enzyme_instance_ids", "catalyst_node_ids", "hypothesis_type", "substrate_names", "product_names", "source_equation", "matched_feature_ids", "supporting_reaction_id", "assumptions", "reference_scope"])

    adjacency = collections.defaultdict(list)
    edge_signs = {}
    for source, interaction, target in sorted(keys):
        adjacency[source].append(target)
        edge_signs[source, target] = interaction
    upstream = sorted(node for node, row in attributes.items() if row.get("type") == "upstream_input")
    downstream = sorted(node for node, row in attributes.items() if finite(row.get("level")) == 0 or row.get("type") == "level0")
    reverse_adjacency = collections.defaultdict(list)
    for source, _, target in keys:
        reverse_adjacency[target].append(source)
    uncontrollable = endpoints - reachable(upstream, adjacency)
    unobservable = endpoints - reachable(downstream, reverse_adjacency)
    if uncontrollable or unobservable:
        raise ValueError(f"Disconnected solution nodes: {len(uncontrollable)} unreachable from intended TFs; {len(unobservable)} cannot reach a measured lipid")
    def collect_paths(path_adjacency):
        rows, details = [], []
        for seed in upstream:
            distance, parent, bottleneck = path_tree(seed, path_adjacency, attributes)
            for target in downstream:
                if target not in distance or (args.max_path_edges is not None and distance[target] > args.max_path_edges):
                    continue
                nodes = [target]
                while nodes[-1] != seed:
                    nodes.append(parent[nodes[-1]])
                nodes.reverse()
                net_sign = math.prod(edge_signs[a, b] for a, b in zip(nodes, nodes[1:]))
                tf_activity = finite(attributes[seed].get("upstream_activity"))
                weight = finite(attributes[target].get("score"), 0)
                reactions = sorted({metadata.get(node, {}).get("reaction_id") for node in nodes} - {None, ""})
                resources = sorted({entry.get("resource", "") for a, b in zip(nodes, nodes[1:]) for entry in support[a, edge_signs[a, b], b]})
                row = {"path_id": "path_" + hashlib.sha256("\n".join(nodes).encode()).hexdigest()[:16],
                       "upstream_tf": seed, "upstream_label": label(seed), "downstream_node": target, "downstream_label": label(target),
                       "edge_count": distance[target], "net_sign": net_sign, "upstream_activity": tf_activity, "downstream_weight": weight,
                       "agrees_with_upstream_activity_and_weight": (sign(tf_activity) * net_sign == sign(weight)) if tf_activity is not None else "",
                       "minimum_absolute_inferred_score": bottleneck[target], "reaction_ids": "|".join(reactions),
                       "species_rule_reactions": sum(r.startswith("LIPID_RULE:") for r in reactions),
                       "resources": "|".join(resources), "path_labels": " → ".join(label(node) for node in nodes), "path_node_ids": "|".join(nodes)}
                rows.append(row)
                details.append(dict(row, nodes=nodes, edges=[{"source": a, "interaction": edge_signs[a, b], "target": b} for a, b in zip(nodes, nodes[1:])]))
        rows.sort(key=lambda row: (-abs(row["downstream_weight"]), -row["minimum_absolute_inferred_score"], row["edge_count"], row["upstream_tf"], row["downstream_node"]))
        return rows, details

    path_rows, path_details = collect_paths(adjacency)
    path_stem = "candidate_paths" if args.stem == "solution_primary" else "candidate_paths_" + args.stem
    path_fields = ["path_id", "upstream_tf", "upstream_label", "downstream_node", "downstream_label", "edge_count", "net_sign", "upstream_activity", "downstream_weight", "agrees_with_upstream_activity_and_weight", "minimum_absolute_inferred_score", "reaction_ids", "species_rule_reactions", "resources", "path_labels", "path_node_ids"]
    write_csv(output / (path_stem + ".csv"), path_rows, path_fields)
    with (output / (path_stem + ".jsonl")).open("w") as handle:
        for row in sorted(path_details, key=lambda row: row["path_id"]):
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    # Independent route semantics audit: a protein regulatory route may end in
    # exactly one annotated catalyst -> enzyme instance -> observed lipid step.
    # Neither metabolites nor reaction instances may be intermediate shortcuts.
    catalyst_adjacency = collections.defaultdict(list)
    for source, interaction, target in sorted(keys):
        source_protein = source.startswith(("Gene__", "UniProt__"))
        target_protein = target.startswith(("Gene__", "UniProt__"))
        source_enzyme = source.startswith("Enzyme__")
        target_enzyme = target.startswith("Enzyme__")
        records = support[source, interaction, target]
        regulatory_edge = source_protein and target_protein
        catalyst_edge = source_protein and target_enzyme and source in metadata.get(target, {}).get("catalysts", []) and any(row.get("mechanism") == "catalyst_enables_reaction" for row in records)
        measured_product_edge = source_enzyme and target in downstream and any(row.get("mechanism") == "product_formation" for row in records)
        if regulatory_edge or catalyst_edge or measured_product_edge:
            catalyst_adjacency[source].append(target)
    catalyst_rows, catalyst_details = collect_paths(catalyst_adjacency)
    audits = {"bounded6": [row for row in path_rows if row["edge_count"] <= 6],
              "catalyst_only": catalyst_rows,
              "catalyst_only_bounded6": [row for row in catalyst_rows if row["edge_count"] <= 6]}
    audit_counts = []
    for name, rows in audits.items():
        write_csv(output / (path_stem + "_" + name + ".csv"), rows, path_fields)
        audit_counts.append({"audit": name, "reachable_pairs": len(rows), "upstream_tfs": len({row["upstream_tf"] for row in rows}),
                             "downstream_nodes": len({row["downstream_node"] for row in rows}),
                             "maximum_path_edges": max((row["edge_count"] for row in rows), default=0)})
    for row in catalyst_details:
        assert all(node.startswith(("Gene__", "UniProt__")) for node in row["nodes"][:-2])
        assert row["nodes"][-2].startswith("Enzyme__") and row["nodes"][-1] in downstream
    write_csv(output / (path_stem + "_route_audit.csv"), audit_counts,
              ["audit", "reachable_pairs", "upstream_tfs", "downstream_nodes", "maximum_path_edges"])
    if args.stem == "solution_relaxed":
        bounded_details = [row for row in catalyst_details if row["edge_count"] <= 6]
        union_edges = {edge_key(edge) for row in bounded_details for edge in row["edges"]}
        union_nodes = {node for row in bounded_details for node in row["nodes"]}
        write_csv(output / "solution_catalyst_bounded6_SIF.csv",
                  [{"source": a, "interaction": s, "target": b} for a, s, b in sorted(union_edges)], ["source", "interaction", "target"])
        write_csv(output / "solution_catalyst_bounded6_ATT.csv", [attributes[node] for node in sorted(union_nodes)],
                  list(raw_attributes[0]) if raw_attributes else ["source", "score", "level", "type", "RNA_input", "upstream_activity"])
        with (output / "solution_catalyst_bounded6_selection.json").open("w") as handle:
            json.dump({"source_solution": "solution_relaxed", "basis_path_file": path_stem + "_catalyst_only_bounded6.csv",
                       "basis_paths": len(bounded_details), "nodes": len(union_nodes), "edges": len(union_edges),
                       "score_policy": "Exact fixed scores and upstream annotations copied from relaxed solution; no refitting.",
                       "selection": "Union of one shortest catalyst-only path per TF/observed-lipid pair where shortest distance <=6; shortest-path ties maximize minimum absolute inferred score.",
                       "caution": "The union graph can contain additional longer routes. The 6-edge bound applies to the selected basis paths, not every possible route through the union."}, handle, indent=2)

    def select_illustration(rows):
        illustrated, selected_endpoints = [], set()
        for row in rows:
            if row["downstream_node"] not in selected_endpoints:
                illustrated.append(row)
                selected_endpoints.add(row["downstream_node"])
            if len(illustrated) == 3:
                break
        return illustrated

    write_csv(output / (path_stem + "_illustrated.csv"), select_illustration(path_rows), path_fields)
    write_csv(output / (path_stem + "_catalyst_only_bounded6_illustrated.csv"), select_illustration(audits["catalyst_only_bounded6"]), path_fields)
    inconsistent_edges = [key for key in keys if sign(finite(attributes[key[0]].get("score"), 0)) * sign(finite(attributes[key[2]].get("score"), 0)) != key[1]]
    summary = {
        "solution_stem": args.stem, "nodes": len(endpoints), "edges": len(keys), "upstream_tfs": len(upstream),
        "reported_path_length_cap": args.max_path_edges,
        "measured_downstream_nodes": len(downstream), "reachable_tf_lipid_pairs": len(path_rows),
        "upstream_tfs_with_paths": len({r["upstream_tf"] for r in path_rows}), "downstream_nodes_with_paths": len({r["downstream_node"] for r in path_rows}),
        "edge_evidence_records": sum(len(value) for value in support.values()), "reaction_records": len(reaction_records),
        "edges_with_species_rule_hypotheses": sum(row["contains_species_rule_hypothesis"] for row in evidence_rows),
        "inconsistent_scored_edges": len(inconsistent_edges),
        "nodes_unreachable_from_intended_upstream": len(uncontrollable), "nodes_unable_to_reach_observation": len(unobservable),
        "maximum_chosen_path_edges": max((row["edge_count"] for row in path_rows), default=0),
        "chosen_paths_longer_than_six_edges": sum(row["edge_count"] > 6 for row in path_rows),
        "route_audits": audit_counts,
        "catalyst_only_route_definition": "Gene__/UniProt__ regulatory intermediates, then exactly one annotated catalyst->Enzyme__->measured lipid product. No internal metabolites or reaction instances. This audit does not change scores or validate exact substrate specificity.",
        "path_selection": "One shortest path per reachable TF/observation pair; ties maximize minimum absolute inferred MOON score, then deterministic node-ID ordering. Level-0 weights excluded from this minimum.",
        "illustration_selection": "Up to three distinct downstream nodes ordered by absolute Factor 1 weight; for each choose greatest minimum absolute inferred score, then shortest path, then TF ID.",
        "limitations": ["Paths are qualitative directed hypotheses, not causal identification or flux predictions.", "LINEX2_rule species reactions transfer a generic class reaction to the reported chains; source references do not establish the exact measured substrate specificity.", "A shortest path omits alternative routes and is an auditable example, not an optimal mechanistic explanation.", "Scoring layer depth does not cap the length of every directed TF/observation route in the union graph; exact edge counts are reported.", "Protein node scores reflect activity hypotheses; they are not measurements of protein abundance.", "Level-0 Factor 1 weights and inferred MOON scores are distinct quantities and should not be compared numerically."]}
    with (output / (evidence_stem + "_summary.json")).open("w") as handle:
        json.dump(summary, handle, indent=2)
    if inconsistent_edges:
        raise ValueError("Sign-inconsistent edges detected: " + str(inconsistent_edges[:5]))
    print(json.dumps(summary, indent=2))
    if not args.skip_plot:
        plotter = Path(__file__).with_name("plot_factor1_network.R")
        subprocess.run(["Rscript", str(plotter), str(output), evidence_stem, path_stem], check=True)
        if path_rows:
            subprocess.run(["Rscript", str(plotter), str(output), evidence_stem, path_stem + "_catalyst_only_bounded6"], check=True)


if __name__ == "__main__":
    main()
