#!/usr/bin/env python3
"""Audit and export a supplied Factor 1 lipid -> TF MOON solution (no refitting).

Level-zero scores are observed RNA-derived TF activities. Inferred MOON scores,
upstream lipid weights, downstream TF activities and RNA inputs remain separate.
One shortest route is reported per reachable input pair, selecting the greatest
minimum absolute inferred MOON score with deterministic signed-edge tie breaking in the shortest-path DAG. Substrate participation is never labelled regulatory signaling.
Only eligible exact signed PKN evidence is joined; original IDs are retained.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import html
import json
import math
from pathlib import Path
import textwrap
import xml.etree.ElementTree as ET

ELIGIBLE_VARIANTS = {"core", "lipid_curated", "lipid_rules"}
REGULATORY_MECHANISMS = {
    "ligand_receptor_or_channel_activity_modulation",
    "nuclear_receptor_activity_modulation", "enzyme_activity_modulation",
}
REACTION_MECHANISMS = {"substrate_participation", "product_formation", "catalyst_enables_reaction"}
PROTEIN_TYPES = {"gene_product", "protein_unmapped"}


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value, default=None):
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (ValueError, TypeError):
        return default


def truth(value):
    return str(value).lower() in {"true", "1", "yes"}


def edge_key(row):
    value = numeric(row["interaction"])
    if value not in {-1, 1}:
        raise ValueError(f"Non-signed interaction: {row}")
    return row["source"], int(value), row["target"]


def joined(values):
    return "|".join(sorted({str(value) for value in values if value not in (None, "")}))


def references(entries):
    result = set()
    for entry in entries:
        value = entry.get("references", "")
        result.update(map(str, value if isinstance(value, list) else [value]))
    return joined(result)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def shortest_routes(seed, adjacency, attributes):
    """Shortest-path DAG dynamic programming; no exponential path enumeration."""
    distance = {seed: 0}
    best = {seed: ()}
    bottleneck = {seed: abs(numeric(attributes[seed].get("moon_score"), 0))}
    frontier = [seed]
    while frontier:
        following = set()
        for source in sorted(frontier):
            for edge in adjacency.get(source, []):
                target = edge[2]
                steps = distance[source] + 1
                if target in distance and distance[target] < steps:
                    continue
                target_score = numeric(attributes[target].get("moon_score"), math.inf)
                candidate_score = min(bottleneck[source], abs(target_score))
                candidate = best[source] + (edge,)
                if (target not in distance or candidate_score > bottleneck[target] or
                        (candidate_score == bottleneck[target] and candidate < best[target])):
                    distance[target], bottleneck[target], best[target] = steps, candidate_score, candidate
                following.add(target)
        frontier = sorted(following)
    return distance, best, bottleneck


def graphml(path, nodes, edges, description):
    ns = "http://graphml.graphdrawing.org/xmlns"
    ET.register_namespace("", ns)
    tag = lambda name: "{" + ns + "}" + name
    root = ET.Element(tag("graphml"))
    keys = {}
    numeric_fields = {"score", "moon_score", "upstream_lipid_weight", "downstream_tf_activity", "RNA_input"}
    int_fields = {"interaction", "level", "evidence_records"}
    for domain, rows in (("node", nodes), ("edge", edges), ("graph", [{"description": description}])):
        for field in dict.fromkeys(key for row in rows for key in row):
            key_id = "k" + str(len(keys))
            typ = "double" if field in numeric_fields else "int" if field in int_fields else "string"
            keys[domain, field] = key_id, typ
            ET.SubElement(root, tag("key"), {"id": key_id, "for": domain, "attr.name": field, "attr.type": typ})
    graph = ET.SubElement(root, tag("graph"), {"id": "Factor1_lipid_to_TF", "edgedefault": "directed"})

    def add_data(element, domain, row):
        for field, value in row.items():
            if value is None or value in ("", "NA", "NaN"):
                continue
            key_id, typ = keys[domain, field]
            if typ in {"double", "int"}:
                number = numeric(value)
                if number is None:
                    continue
                value = int(number) if typ == "int" else number
            ET.SubElement(element, tag("data"), {"key": key_id}).text = str(value)

    add_data(graph, "graph", {"description": description})
    identifiers = {row["node_id"]: "n" + str(index) for index, row in enumerate(nodes)}
    for row in nodes:
        add_data(ET.SubElement(graph, tag("node"), {"id": identifiers[row["node_id"]]}), "node", row)
    for index, row in enumerate(edges):
        element = ET.SubElement(graph, tag("edge"), {"id": "e" + str(index), "source": identifiers[row["source"]], "target": identifiers[row["target"]]})
        add_data(element, "edge", row)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def illustrate_paths(path, routes, node_attributes, edge_attributes, limit):
    """Compact SVG of representative signed routes, not a calibrated score plot."""
    chosen = sorted(routes, key=lambda row: (row["path_length"], row["upstream_lipid"], row["downstream_tf"]))[:limit]
    row_height = 134
    max_steps = max([row["path_length"] for row in chosen] + [1])
    box_width, spacing, margin = 126, 34, 25
    width = max(1000, 2 * margin + (max_steps + 1) * (box_width + spacing))
    height = 136 + row_height * max(1, len(chosen))
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>',
           '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#555"/></marker></defs>',
           '<g font-family="Arial,sans-serif" fill="#222">',
           '<text x="25" y="30" font-size="21">Factor 1: measured lipid to TF hypothesis routes</text>',
           '<text x="25" y="54" font-size="13">Purple edges: actual chemical–protein activity modulation. Grey: other signed PKN edges. Dashed: species-rule hypothesis.</text>',
           '<text x="25" y="76" font-size="13">Blue: measured lipid. Yellow: other chemical. Grey: reaction instance. Green: protein. Red: observed TF endpoint.</text>',
           f'<text x="25" y="98" font-size="13">Showing {len(chosen)} of {len(routes)} representative shortest routes; no scores are rescaled or refitted.</text>']
    if not chosen:
        svg.append('<text x="25" y="164" font-size="17">No measured-lipid to observed-TF route in the supplied solution.</text>')
    for index, row in enumerate(chosen):
        y = 142 + row_height * index
        route_edges = [tuple(item) for item in json.loads(row["signed_edges_json"])]
        route_nodes = [row["upstream_lipid"]] + [edge[2] for edge in route_edges]
        title = f'{row["upstream_lipid_label"]} → {row["downstream_tf_label"]}: {row["path_classification"]}; {row["path_length"]} edges'
        svg.append(f'<text x="25" y="{y}" font-size="13">{html.escape(title)}</text>')
        for j, node in enumerate(route_nodes):
            x = margin + j * (box_width + spacing)
            attr = node_attributes[node]
            fill = "#b9d8f2" if j == 0 else "#f0b4b8" if j == len(route_nodes) - 1 else "#f8e5a4" if attr["node_type"] == "metabolite" else "#ddd" if attr["node_type"] == "enzyme_reaction_instance" else "#c1dfbf"
            svg.append(f'<rect x="{x}" y="{y + 14}" width="{box_width}" height="72" rx="6" fill="{fill}" stroke="#777"/>')
            lines = textwrap.wrap(attr["label"], width=18, break_long_words=True)[:4]
            for k, line in enumerate(lines):
                svg.append(f'<text x="{x + 5}" y="{y + 30 + 15 * k}" font-size="11">{html.escape(line)}</text>')
            svg.append(f'<title>{html.escape(node + ": " + attr["label"])}</title>')
            if j < len(route_edges):
                edge = route_edges[j]
                info = edge_attributes[edge]
                color = "#823cad" if info["is_regulatory_chemical_to_protein"] else "#555"
                dash = ' stroke-dasharray="5 3"' if info["contains_species_rule_hypothesis"] else ""
                svg.append(f'<line x1="{x + box_width}" y1="{y + 49}" x2="{x + box_width + spacing - 4}" y2="{y + 49}" stroke="{color}" stroke-width="2" marker-end="url(#arrow)"{dash}/>')
                svg.append(f'<text x="{x + box_width + 12}" y="{y + 39}" font-size="15" fill="{color}">{"+" if edge[1] == 1 else "−"}</text>')
    svg.extend(['</g>', '</svg>'])
    path.write_text("\n".join(svg) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--pkn", type=Path, default=Path("results/causal_modelling/cosmos_mouse_pkn"))
    parser.add_argument("--input-mapping", type=Path, required=True, help="Accepted exact lipid mapping CSV with feature, node_id, included columns")
    parser.add_argument("--prior-csv", type=Path, default=None,
                        help="Exact signed prior CSV; defaults to --pkn/pkn_lipid_rules.csv")
    parser.add_argument("--node-catalog", type=Path, default=None,
                        help="Node metadata JSONL; defaults to --pkn/nodes.jsonl")
    parser.add_argument("--node-aliases", type=Path, default=None,
                        help="Explicit original_node/normalized_node CSV; project provenance endpoints only using these supplied aliases")
    parser.add_argument("--stem", default="network")
    parser.add_argument("--max-illustrated-paths", type=int, default=12)
    parser.add_argument("--skip-plot", action="store_true")
    args = parser.parse_args()
    if args.max_illustrated_paths < 1:
        parser.error("--max-illustrated-paths must be positive")
    if Path(args.stem).name != args.stem:
        parser.error("--stem must be a filename stem, not a path")
    output, pkn = args.results, args.pkn
    prior_path = args.prior_csv or pkn / "pkn_lipid_rules.csv"
    catalog_path = args.node_catalog or pkn / "nodes.jsonl"
    aliases = {}
    if args.node_aliases is not None:
        for alias in read_csv(args.node_aliases):
            original, normalized = alias.get("original_node"), alias.get("normalized_node")
            if not original or not normalized:
                raise ValueError("Alias rows require original_node and normalized_node")
            if original in aliases and aliases[original] != alias:
                raise ValueError(f"Conflicting or duplicate alias metadata for {original}")
            aliases[original] = alias

    def project_evidence(entry):
        """Explicit one-step mapping only; no inferred aliases or transitive closure."""
        result = dict(entry)
        for field in ("source", "target"):
            original = entry[field]
            alias = aliases.get(original)
            result["original_" + field] = original
            result[field] = alias["normalized_node"] if alias is not None else original
            result[field + "_mapping_basis"] = alias.get("basis", "explicit_node_alias") if alias is not None else "unchanged_original_node"
            if alias is not None:
                result[field + "_node_alias"] = alias
        return result
    sif, att = output / (args.stem + "_SIF.csv"), output / (args.stem + "_ATT.csv")
    edges, raw_nodes = read_csv(sif), read_csv(att)
    raw_attributes = {row.get("source", row.get("node_id")): row for row in raw_nodes}
    if len(raw_attributes) != len(raw_nodes) or None in raw_attributes:
        raise ValueError("Missing or duplicate node IDs")
    keys = {edge_key(row) for row in edges}
    if len(keys) != len(edges):
        raise ValueError("Duplicate signed solution edges")
    endpoints = {node for source, _, target in keys for node in (source, target)}
    if endpoints != set(raw_attributes):
        raise ValueError("SIF/ATT node sets differ")
    metadata = {row["node_id"]: row for row in map(json.loads, catalog_path.open())}
    if not endpoints <= set(metadata):
        raise ValueError("Some solution nodes have no PKN metadata")
    original_keys = {edge_key(row) for row in read_csv(prior_path)}
    if not keys <= original_keys:
        raise ValueError("Solution contains signed edges absent from the PKN")
    mapped_features = collections.defaultdict(set)
    for row in read_csv(args.input_mapping):
        if truth(row.get("included")):
            mapped_features[row["node_id"]].add(row["feature"])
    support = collections.defaultdict(list)
    direct_input_regulatory_edges = set()
    with gzip.open(pkn / "edge_evidence.jsonl.gz", "rt") as handle:
        for row in map(json.loads, handle):
            if row.get("first_variant") not in ELIGIBLE_VARIANTS:
                continue
            row = project_evidence(row)
            key = edge_key(row)
            regulatory = (row.get("mechanism") in REGULATORY_MECHANISMS and
                          metadata.get(key[0], {}).get("node_type") == "metabolite" and
                          metadata.get(key[2], {}).get("node_type") in PROTEIN_TYPES)
            if regulatory and key[0] in mapped_features and key in original_keys:
                direct_input_regulatory_edges.add(key)
            if key in keys:
                support[key].append(row)
    if set(support) != keys:
        raise ValueError("Some exact solution edges lack eligible PKN evidence")

    def label(node):
        info = metadata[node]
        if mapped_features.get(node):
            return "/".join(sorted(mapped_features[node]))
        if info.get("node_type") == "enzyme_reaction_instance":
            catalysts = "/".join(x.replace("Gene__", "") for x in info.get("catalysts", []))
            reaction = info.get("reaction_id", "reaction")
            return catalysts + " [" + reaction + "]"
        return info.get("label", node)

    node_rows, attributes = [], {}
    for node, raw in sorted(raw_attributes.items()):
        row = dict(raw)
        level = numeric(row.get("level"))
        if level is None or level < 0 or int(level) != level:
            raise ValueError(f"Missing or invalid scoring level: {node}")
        score = numeric(row.get("score"))
        if score is None:
            raise ValueError(f"Non-finite solution score: {node}")
        if level == 0:
            activity = numeric(row.get("downstream_tf_activity"), score)
            if not math.isclose(score, activity, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f"Level-zero score differs from observed TF activity: {node}")
            row["downstream_tf_activity"] = activity
            if metadata[node]["node_type"] not in PROTEIN_TYPES:
                raise ValueError(f"Level-zero endpoint is not a protein/TF: {node}")
        if row.get("type") == "upstream_input":
            if node not in mapped_features or numeric(row.get("upstream_lipid_weight")) is None:
                raise ValueError(f"Upstream lipid lacks exact mapping or separate measured weight: {node}")
        row.update(node_id=node, label=label(node), pkn_label=metadata[node].get("label", ""),
                   node_type=metadata[node]["node_type"], observed_features=joined(mapped_features.get(node, [])),
                   moon_score=score if level > 0 else "", score_semantics="inferred_MOON_score" if level > 0 else "observed_RNA_derived_TF_activity",
                   reaction_id=metadata[node].get("reaction_id", ""), catalysts=joined(metadata[node].get("catalysts", [])))
        node_rows.append(row)
        attributes[node] = row
    node_fields = list(dict.fromkeys(["node_id", "label", "node_type", "moon_score", "upstream_lipid_weight", "downstream_tf_activity", "RNA_input"] + [key for row in node_rows for key in row]))
    prefix = args.stem + "_evidence"
    write_csv(output / (prefix + "_nodes.csv"), node_rows, node_fields)
    edge_rows, annotated_edges = [], {}
    with gzip.open(output / (prefix + ".jsonl.gz"), "wt") as handle:
        for key in sorted(keys):
            entries = support[key]
            for entry in entries:
                handle.write(json.dumps(entry, sort_keys=True) + "\n")
            source, interaction, target = key
            mechanisms = {x.get("mechanism", "") for x in entries}
            regulatory = bool(mechanisms & REGULATORY_MECHANISMS and metadata[source]["node_type"] == "metabolite" and metadata[target]["node_type"] in PROTEIN_TYPES)
            kind = "chemical_to_protein_activity_regulation" if regulatory else "reaction_projection" if mechanisms & REACTION_MECHANISMS else "protein_regulation" if metadata[source]["node_type"] in PROTEIN_TYPES and metadata[target]["node_type"] in PROTEIN_TYPES else "other"
            row = dict(source=source, interaction=interaction, target=target, source_label=label(source), target_label=label(target),
                       edge_class=kind, is_regulatory_chemical_to_protein=regulatory,
                       is_direct_measured_lipid_regulation=regulatory and source in mapped_features,
                       mechanisms=joined(mechanisms), resources=joined(x.get("resource") for x in entries),
                       reaction_ids=joined(x.get("reaction_id") for x in entries), evidence_records=len(entries),
                       contains_species_rule_hypothesis=any(x.get("resource") == "LINEX2_rule" or str(x.get("reaction_id", "")).startswith("LIPID_RULE:") for x in entries),
                       evidence_tiers=joined(x.get("evidence_tier") for x in entries), references=references(entries),
                       original_edge_ids=joined(x.get("original_edge_id") for x in entries),
                       species_evidence=joined(x.get("species_evidence") for x in entries),
                       original_sources=joined(x["original_source"] for x in entries),
                       original_targets=joined(x["original_target"] for x in entries),
                       source_mapping_basis=joined(x["source_mapping_basis"] for x in entries),
                       target_mapping_basis=joined(x["target_mapping_basis"] for x in entries),
                       original_signed_edges_json=json.dumps(sorted({(x["original_source"], int(x["interaction"]), x["original_target"]) for x in entries})))
            edge_rows.append(row)
            annotated_edges[key] = row
    edge_fields = ["source", "interaction", "target", "source_label", "target_label", "edge_class", "is_regulatory_chemical_to_protein", "is_direct_measured_lipid_regulation", "mechanisms", "resources", "reaction_ids", "evidence_records", "contains_species_rule_hypothesis", "evidence_tiers", "references", "original_edge_ids", "species_evidence", "original_sources", "original_targets", "source_mapping_basis", "target_mapping_basis", "original_signed_edges_json"]
    write_csv(output / (prefix + ".csv"), edge_rows, edge_fields)
    regulatory_rows = [row for row in edge_rows if row["is_regulatory_chemical_to_protein"]]
    write_csv(output / (args.stem + "_regulatory_transitions.csv"), regulatory_rows, edge_fields)
    description = ("Factor 1 lipid-to-TF solution. Level-zero score is an observed RNA-derived TF activity, not a lipid weight. "
                   "moon_score, upstream_lipid_weight, downstream_tf_activity and RNA_input are separate quantities. "
                   "Substrate participation is a qualitative metabolic projection and does not establish signaling or flux causality. "
                   "Species-specific reaction rules are hypotheses; chemical-to-protein regulatory edges are explicitly identified. " +
                   ("Evidence endpoints were projected using the explicit supplied alias table; original endpoints and mapping bases remain in exported provenance. " if aliases else "") +
                   "All solution edges were verified against the exact supplied prior: " + str(prior_path) + ".")
    graphml(output / (args.stem + ".graphml"), node_rows, edge_rows, description)
    adjacency = collections.defaultdict(list)
    for key in sorted(keys):
        adjacency[key[0]].append(key)
    upstream = sorted(node for node in attributes if attributes[node].get("type") == "upstream_input")
    downstream = sorted(node for node in attributes if numeric(attributes[node].get("level")) == 0)
    paths, path_steps = [], []
    for seed in upstream:
        distance, best, bottleneck = shortest_routes(seed, adjacency, attributes)
        for target in downstream:
            if target not in distance:
                continue
            route = best[target]
            route_nodes = [seed] + [edge[2] for edge in route]
            transitions = [(i, edge) for i, edge in enumerate(route) if annotated_edges[edge]["is_regulatory_chemical_to_protein"]]
            first = transitions[0] if transitions else None
            classification = "direct_measured_lipid_signaling_hypothesis" if first and first[0] == 0 else "indirect_metabolic_availability_hypothesis" if first else "no_explicit_regulatory_transition"
            product = math.prod(edge[1] for edge in route)
            route_id = "path_" + str(len(paths) + 1).zfill(5)
            transition_objects = [dict(step=i + 1, chemical=edge[0], chemical_label=label(edge[0]), protein=edge[2], protein_label=label(edge[2]), interaction=edge[1], mechanisms=annotated_edges[edge]["mechanisms"], resources=annotated_edges[edge]["resources"], references=annotated_edges[edge]["references"], original_edge_ids=annotated_edges[edge]["original_edge_ids"]) for i, edge in transitions]
            paths.append(dict(path_id=route_id, upstream_lipid=seed, upstream_lipid_label=label(seed), upstream_lipid_weight=attributes[seed].get("upstream_lipid_weight", ""), upstream_moon_score=attributes[seed]["moon_score"],
                              downstream_tf=target, downstream_tf_label=label(target), downstream_tf_activity=attributes[target]["downstream_tf_activity"],
                              path_length=distance[target], min_abs_inferred_moon_score=bottleneck[target], path_sign=product, path_classification=classification,
                              first_regulatory_transition_step=first[0] + 1 if first else "", first_regulatory_chemical=first[1][0] if first else "", first_regulatory_chemical_label=label(first[1][0]) if first else "", first_regulatory_protein=first[1][2] if first else "",
                              regulatory_transition_count=len(transitions), regulatory_transitions_json=json.dumps(transition_objects, sort_keys=True),
                              contains_species_rule_hypothesis=any(annotated_edges[edge]["contains_species_rule_hypothesis"] for edge in route),
                              reaction_ids=joined(annotated_edges[edge]["reaction_ids"] for edge in route), resources=joined(annotated_edges[edge]["resources"] for edge in route), references=joined(annotated_edges[edge]["references"] for edge in route),
                              node_ids=" -> ".join(route_nodes), node_labels=" -> ".join(label(node) for node in route_nodes), signed_edges_json=json.dumps(route)))
            for index, edge in enumerate(route):
                path_steps.append(dict(path_id=route_id, step=index + 1, **annotated_edges[edge]))
    path_fields = ["path_id", "upstream_lipid", "upstream_lipid_label", "upstream_lipid_weight", "upstream_moon_score", "downstream_tf", "downstream_tf_label", "downstream_tf_activity", "path_length", "min_abs_inferred_moon_score", "path_sign", "path_classification", "first_regulatory_transition_step", "first_regulatory_chemical", "first_regulatory_chemical_label", "first_regulatory_protein", "regulatory_transition_count", "regulatory_transitions_json", "contains_species_rule_hypothesis", "reaction_ids", "resources", "references", "node_ids", "node_labels", "signed_edges_json"]
    write_csv(output / (args.stem + "_candidate_paths.csv"), paths, path_fields)
    write_csv(output / (args.stem + "_candidate_path_steps.csv"), path_steps, ["path_id", "step"] + edge_fields)
    if not args.skip_plot:
        illustrate_paths(output / (args.stem + "_candidate_paths.svg"), paths, attributes, annotated_edges, args.max_illustrated_paths)
    counts = collections.Counter(row["path_classification"] for row in paths)
    summary = dict(nodes=len(node_rows), edges=len(edge_rows), upstream_lipids=len(upstream), downstream_tfs=len(downstream), reachable_input_pairs=len(paths), path_class_counts=dict(counts),
                   regulatory_transitions=len(regulatory_rows), direct_measured_lipid_regulatory_edges_in_solution=sum(row["is_direct_measured_lipid_regulation"] for row in edge_rows),
                   accepted_mapping_nodes=len(mapped_features), direct_exact_mapped_lipid_regulatory_edges_in_full_pkn=len(direct_input_regulatory_edges),
                   species_rule_edges=sum(row["contains_species_rule_hypothesis"] for row in edge_rows),
                   representative_path_policy="One shortest path per reachable input pair; greatest minimum absolute inferred MOON score; deterministic signed-edge tie breaking in the shortest-path DAG. Level-zero observed TF activities excluded from score bottleneck. Not all alternative routes are enumerated.",
                   interpretation=("No exact-mapped measured lipid has an eligible direct regulatory edge in the full PKN. Any selected routes with a later chemical-to-protein transition nominate indirect metabolic availability hypotheses, not direct signaling lipids or confirmed biochemical precursors." if not direct_input_regulatory_edges else "Exact-mapped measured lipids have some eligible direct regulatory PKN edges; route-level flags distinguish direct signaling hypotheses from indirect metabolic availability hypotheses."),
                   limitations=["A directed signed route is a prior-knowledge hypothesis, not evidence of causality, flux or biochemical feasibility in hepatocytes.", "Reaction instances encode substrate participation, product formation and catalyst availability; they do not implement mass balance or enzyme kinetics.", "The actual regulatory chemical may be an unmeasured metabolite or co-product; its abundance or availability is not independently established.", "Shortest-route classifications do not summarize every alternative walk through the solution.", "Lipid weights and RNA-derived TF activities are observational measurements on different scales; MOON scores are inferred statistics."],
                   prior_csv=str(prior_path), node_catalog=str(catalog_path), node_aliases=str(args.node_aliases) if args.node_aliases else None,
                   explicit_node_alias_count=len(aliases),
                   evidence_projection="Only explicitly supplied one-step node aliases; original endpoints and mapping bases retained; no identities inferred and no solution scores changed.",
                   input_sha256={str(path): sha256(path) for path in [sif, att, args.input_mapping, catalog_path, prior_path, pkn / "edge_evidence.jsonl.gz"] + ([args.node_aliases] if args.node_aliases else [])})
    (output / (args.stem + "_signaling_audit.json")).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (output / (args.stem + "_interpretation.md")).write_text("# Lipid-to-TF solution interpretation\n\n" + description + "\n\n" + summary["interpretation"] + "\n\n" + summary["representative_path_policy"] + "\n\n" + "\n".join("- " + item for item in summary["limitations"]) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key not in {"input_sha256", "limitations"}}, indent=2))


if __name__ == "__main__":
    main()
