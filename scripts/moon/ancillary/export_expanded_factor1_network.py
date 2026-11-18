#!/usr/bin/env python3
"""Export/audit an existing expanded-PKN Factor 1 TF-to-lipid solution.

Standard-library only. Uses exact signed IDs; never fits scores, remaps inputs,
changes the supplied solution, or modifies either PKN. Full reaction records are
preserved alongside participant bookkeeping and explicit currency-node audits.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[3]
PROTEINS = {"gene_product", "protein_unmapped"}
CURRENCY_LABEL = re.compile(
    r"^(H2O|water|H\(\+\)|H\+|proton|CoA|coenzyme A|ATP|ADP|AMP|GTP|GDP|"
    r"NAD\+|NAD\(\+\)|NADH|NADP\+|NADP\(\+\)|NADPH|phosphate|"
    r"diphosphate|Pi|PPi|S-adenosyl-L-methionine|S-adenosyl-L-homocysteine)$", re.I)
CURRENCY_IDS = {"CHEBI:" + ident for ident in (
    "15377 24636 15378 1146900 57287 15422 30616 16761 456216 16027 456215 "
    "15996 37565 17552 58189 15846 57540 18009 58349 16908 57945 16474 57783 "
    "18367 43474 18361 33019 59789 57856").split()}
MECHANISMS = {"catalyst_enables_reaction", "substrate_participation", "product_formation"}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def jsonlines(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def write_jsonl(path, records):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as handle:
        for row in records:
            handle.write(encoded(row) + "\n")


def encoded(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False)


def finite(value, default=None):
    try:
        value = float(value)
        return value if math.isfinite(value) else default
    except (TypeError, ValueError):
        return default


def sign(value):
    return 1 if value > 0 else -1 if value < 0 else 0


def edge_key(row):
    value = finite(row["interaction"])
    if value not in {-1, 1}:
        raise ValueError("Invalid signed edge: " + encoded(row))
    return row["source"], int(value), row["target"]


def joined(values):
    return "|".join(sorted({str(v) for v in values if v not in (None, "")}))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def currency_node(meta):
    return meta.get("node_type") == "metabolite" and (
        meta.get("chemical_id") in CURRENCY_IDS or
        CURRENCY_LABEL.fullmatch(meta.get("label", "")))


def reachable(seeds, adjacency):
    visited, queue = set(seeds), deque(seeds)
    while queue:
        for target in adjacency.get(queue.popleft(), []):
            if target not in visited:
                visited.add(target)
                queue.append(target)
    return visited


def shortest_tree(seed, adjacency, attributes):
    """Shortest paths, greatest inferred-score bottleneck, deterministic ties."""
    distance, parent = {seed: 0}, {}
    bottleneck = {seed: abs(float(attributes[seed]["score"]))}
    frontier = [seed]
    while frontier:
        following = set()
        for source in sorted(frontier):
            for target in sorted(adjacency.get(source, [])):
                d = distance[source] + 1
                if target in distance and distance[target] < d:
                    continue
                v = (math.inf if float(attributes[target]["level"]) == 0 else
                     abs(float(attributes[target]["score"])))
                b = min(bottleneck[source], v)
                if (target not in distance or b > bottleneck[target] or
                        (b == bottleneck[target] and source < parent[target])):
                    distance[target], parent[target], bottleneck[target] = d, source, b
                following.add(target)
        frontier = following
    return distance, parent, bottleneck


def graphml(path, nodes, edges):
    ns = "http://graphml.graphdrawing.org/xmlns"
    ET.register_namespace("", ns)
    tag = lambda name: "{" + ns + "}" + name
    root = ET.Element(tag("graphml"))
    fields = {"node": list(dict.fromkeys(k for r in nodes for k in r)),
              "edge": list(dict.fromkeys(k for r in edges for k in r))}
    numbers = {"score", "moon_score", "factor1_weight", "upstream_activity", "RNA_input"}
    integers = {"level", "interaction", "evidence_records"}
    keys = {}
    for domain, names in fields.items():
        for name in names:
            ident = "k" + str(len(keys))
            typ = "double" if name in numbers else "int" if name in integers else "string"
            keys[domain, name] = ident, typ
            ET.SubElement(root, tag("key"), {"id": ident, "for": domain,
                                            "attr.name": name, "attr.type": typ})
    graph = ET.SubElement(root, tag("graph"), {"id": "Factor1_expanded_MOON", "edgedefault": "directed"})
    def data(element, domain, row):
        for name, value in row.items():
            if value in (None, "", "NA", "NaN"):
                continue
            key, typ = keys[domain, name]
            if typ != "string":
                value = finite(value)
                if value is None:
                    continue
                value = int(value) if typ == "int" else value
            ET.SubElement(element, tag("data"), {"key": key}).text = str(value)
    ids = {r["node_id"]: "n" + str(i) for i, r in enumerate(nodes)}
    for row in nodes:
        data(ET.SubElement(graph, tag("node"), {"id": ids[row["node_id"]]}), "node", row)
    for i, row in enumerate(edges):
        data(ET.SubElement(graph, tag("edge"), {"id": "e" + str(i), "source": ids[row["source"]],
                                               "target": ids[row["target"]]}), "edge", row)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def participant_node(participant):
    ident = participant["id"]
    if re.fullmatch(r"(?:CHEBI:)?\d+", ident):
        ident = "CHEBI:" + ident.split(":")[-1]
    return "Metab__" + ident + "_" + (participant.get("compartment") or "unknown")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Existing fit with network_SIF/ATT.csv")
    parser.add_argument("--pkn", type=Path, default=ROOT / "results/causal_modelling/cosmos_mouse_pkn_expanded_species")
    parser.add_argument("--baseline-pkn", type=Path, default=ROOT / "results/causal_modelling/cosmos_mouse_pkn")
    parser.add_argument("--pkn-file", default="pkn_expanded_chain.csv", help="Filename inside --pkn")
    parser.add_argument("--pruned-pkn", type=Path, help="Optional final uncompressed scoring PKN CSV")
    parser.add_argument("--input-dir", type=Path, help="Prepared input-vector directory; defaults to --results")
    parser.add_argument("--output", type=Path, help="New/empty audit directory; default RESULTS/expanded_network_audit")
    parser.add_argument("--stem", default="network", help="Input SIF/ATT stem, default network")
    args = parser.parse_args()
    results, pkn, base = args.results.resolve(), args.pkn.resolve(), args.baseline_pkn.resolve()
    out = (args.output or results / "expanded_network_audit").resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error("Use a new or empty audit output directory: " + str(out))
    sif_path, att_path = results / (args.stem + "_SIF.csv"), results / (args.stem + "_ATT.csv")
    graph_path = pkn / args.pkn_file
    node_paths = [base / "nodes.jsonl", pkn / "nodes.jsonl"]
    evidence_paths = [base / "edge_evidence.jsonl.gz", pkn / "extension_edge_evidence.jsonl.gz"]
    master_paths = [base / name for name in ["mouse_gem_reactions.jsonl", "rhea_reactions.jsonl", "lipid_rules_reactions.jsonl"]]
    master_paths.append(pkn / "expanded_species_reactions.jsonl")
    dependencies = [sif_path, att_path, graph_path, *node_paths, *evidence_paths, *master_paths, Path(__file__).resolve()]
    for path in dependencies:
        if not path.is_file():
            parser.error("Missing required input: " + str(path))
    pruned_path = args.pruned_pkn or results / "pkn_final_uncompressed.csv"
    if args.pruned_pkn and not pruned_path.is_file():
        parser.error("Explicit pruned PKN is missing: " + str(pruned_path))
    pruned = {edge_key(r) for r in read_csv(pruned_path)} if pruned_path.is_file() else None
    if pruned is not None:
        dependencies.append(pruned_path)
    edges, attributes_raw = read_csv(sif_path), read_csv(att_path)
    keys = {edge_key(r) for r in edges}
    if len(keys) != len(edges):
        raise ValueError("Duplicate exact edges in solution")
    attributes = {r.get("source", r.get("node_id", r.get("nodes"))): r for r in attributes_raw}
    endpoints = {n for s, _, t in keys for n in (s, t)}
    if len(attributes) != len(attributes_raw) or set(attributes) != endpoints:
        raise ValueError("SIF/ATT IDs differ or ATT IDs are duplicated")
    for node, row in attributes.items():
        if finite(row.get("score")) is None or finite(row.get("level")) is None:
            raise ValueError("Missing/non-finite score/level: " + str(node))
    full = {edge_key(r) for r in read_csv(graph_path)}
    full_nodes = {n for s, _, t in full for n in (s, t)}
    if not keys <= full or (pruned is not None and not keys <= pruned):
        raise ValueError("Selected signed edges are absent from supplied full/final PKN")
    metadata = {}
    for path in node_paths:
        metadata.update((r["node_id"], r) for r in jsonlines(path))
    if not endpoints <= metadata.keys():
        raise ValueError("Selected IDs lack node metadata")
    support = defaultdict(list)
    for path in evidence_paths:
        for row in jsonlines(path):
            eligible = (row.get("first_variant") in {"core", "lipid_curated", "lipid_rules"}
                        if path == evidence_paths[0] else
                        row.get("overlay_tier") in ({"chain_compatible", "glycan_conditional"}
                            if "glycan" in args.pkn_file else {"chain_compatible"}))
            if eligible and edge_key(row) in keys:
                support[edge_key(row)].append(row)
    if set(support) != keys:
        raise ValueError("Selected edges lack eligible signed provenance: " + str(sorted(keys-set(support))[:5]))
    input_values = {}
    input_dir = (args.input_dir or results).resolve()
    for name in ["upstream_input.csv", "downstream_input.csv"]:
        path = input_dir / name
        if path.is_file():
            input_values[name] = {r["source"]: finite(r.get("score")) for r in read_csv(path)}
            dependencies.append(path)
    upstream = sorted(n for n, r in attributes.items() if r.get("type") == "upstream_input" or
                      n in input_values.get("upstream_input.csv", {}))
    downstream = sorted(n for n, r in attributes.items() if float(r["level"]) == 0)
    adjacency, reverse = defaultdict(list), defaultdict(list)
    pair_sign = {}
    for s, i, t in sorted(keys):
        if (s, t) in pair_sign:
            raise ValueError("Conflicting signs retained in solution")
        adjacency[s].append(t)
        reverse[t].append(s)
        pair_sign[s, t] = i
    inconsistent = [k for k in keys if sign(float(attributes[k[0]]["score"]) *
                                          float(attributes[k[2]]["score"])) != k[1]]
    uncontrollable = endpoints - reachable(upstream, adjacency)
    unobservable = endpoints - reachable(downstream, reverse)
    if inconsistent or uncontrollable or unobservable:
        raise ValueError("Independent solution audit failed: " + encoded(dict(
            inconsistent_edges=inconsistent, unreachable_from_inputs=sorted(uncontrollable),
            cannot_reach_observations=sorted(unobservable))))
    if any(s in downstream for s, _, _ in keys):
        raise ValueError("Selected solution contains outgoing downstream-observation edges")
    def label(node):
        m = metadata[node]
        if m.get("node_type") == "enzyme_reaction_instance":
            rid = m.get("reaction_id", "")
            if rid.startswith(("LIPID_RULE:", "EXPANDED_LIPID_RULE:")):
                rid = "species rule " + rid.split(":")[-1][:8]
            return "/".join(metadata.get(c, {}).get("label", c) for c in m.get("catalysts", [])) + " [" + rid + "]"
        return m.get("label", node)
    node_rows = []
    for n in sorted(endpoints):
        m, r = metadata[n], attributes[n]
        node_rows.append(dict(r, node_id=n, label=label(n), pkn_label=m.get("label", ""),
            node_type=m["node_type"], reaction_id=m.get("reaction_id", ""),
            catalysts=joined(m.get("catalysts", [])), compartment=m.get("compartment", ""),
            chemical_id=m.get("chemical_id", ""), measurement_feature_id=m.get("measurement_feature_id", ""),
            candidate_identity_condition=m.get("candidate_identity_condition", ""),
            factor1_weight=r["score"] if n in downstream else "",
            moon_score=r["score"] if n not in downstream else ""))
    edge_rows = []
    for s, i, t in sorted(keys):
        ee = support[s, i, t]
        edge_rows.append(dict(source=s, interaction=i, target=t, source_label=label(s), target_label=label(t),
            source_score_level=int(float(attributes[s]["level"])), target_score_level=int(float(attributes[t]["level"])),
            score_layer_adjacent=int(float(attributes[s]["level"])) == int(float(attributes[t]["level"])) + 1,
            mechanisms=joined(e.get("mechanism") for e in ee), resources=joined(e.get("resource") for e in ee),
            reaction_ids=joined(e.get("reaction_id") for e in ee), evidence_records=len(ee),
            contains_species_rule_hypothesis=any(e.get("resource") in {"LINEX2_rule", "Measured_species_rule_expanded"} for e in ee),
            references_json=encoded([e.get("references") for e in ee if e.get("references")]),
            template_reaction_ids=joined(e.get("supporting_reaction_id") for e in ee)))
    selected_reactions = defaultdict(list)
    for n in sorted(endpoints):
        if metadata[n]["node_type"] == "enzyme_reaction_instance":
            selected_reactions[metadata[n]["reaction_id"]].append(n)
    records = {}
    for path in master_paths:
        for r in jsonlines(path):
            if r["reaction_id"] in selected_reactions:
                if r["reaction_id"] in records and records[r["reaction_id"]] != r:
                    raise ValueError("Conflicting reaction master records")
                records[r["reaction_id"]] = r
    if records.keys() != selected_reactions.keys():
        raise ValueError("Selected reaction is missing its complete master record")
    participants, reaction_rows = [], []
    pruned_nodes = {n for s, _, t in pruned for n in (s, t)} if pruned is not None else None
    for rid, instances in sorted(selected_reactions.items()):
        record = records[rid]
        reaction_rows.append(dict(reaction_id=rid, resource=record.get("resource", ""),
            supporting_reaction_id=record.get("supporting_reaction_id", ""),
            enzyme_instance_ids=joined(instances), catalyst_node_ids=joined(c for n in instances for c in metadata[n].get("catalysts", [])),
            direction=record.get("direction", ""), source_equation=record.get("source_template_equation") or record.get("equation", ""),
            substrate_names=joined(p.get("name", p["id"]) for p in record["substrates"]),
            product_names=joined(p.get("name", p["id"]) for p in record["products"]),
            matched_feature_ids=joined(record.get("matched_feature_ids", [])),
            assumptions_json=encoded(record.get("assumptions", [])), reference_scope=record.get("reference_scope", "")))
        for n in instances:
            for side in ["substrates", "products"]:
                original_side = ("products" if side == "substrates" else "substrates") if metadata[n].get("reverse") else side
                for index, p in enumerate(record[original_side]):
                    chemical = participant_node(p)
                    k = (chemical, 1, n) if side == "substrates" else (n, 1, chemical)
                    participants.append(dict(reaction_id=rid, enzyme_instance=n, participant_role=side,
                        master_record_side=original_side, participant_index=index, participant_node=chemical,
                        participant_name=p.get("name", p["id"]), stoichiometry=p["stoichiometry"],
                        measurement_feature_id=p.get("measurement_feature_id", ""), unmeasured=p.get("unmeasured", False),
                        node_in_full_pkn=chemical in full_nodes, connector_in_full_pkn=k in full,
                        node_in_final_pkn=chemical in pruned_nodes if pruned_nodes is not None else "not_assessed",
                        connector_in_final_pkn=k in pruned if pruned is not None else "not_assessed",
                        node_in_solution=chemical in endpoints, connector_in_solution=k in keys))
            for catalyst in metadata[n].get("catalysts", []):
                k = (catalyst, 1, n)
                participants.append(dict(reaction_id=rid, enzyme_instance=n, participant_role="catalyst",
                    master_record_side="catalysts", participant_index="", participant_node=catalyst,
                    participant_name=metadata.get(catalyst, {}).get("label", catalyst), stoichiometry="",
                    measurement_feature_id="", unmeasured=False, node_in_full_pkn=catalyst in full_nodes,
                    connector_in_full_pkn=k in full,
                    node_in_final_pkn=catalyst in pruned_nodes if pruned_nodes is not None else "not_assessed",
                    connector_in_final_pkn=k in pruned if pruned is not None else "not_assessed",
                    node_in_solution=catalyst in endpoints, connector_in_solution=k in keys))
    currency = {n for n in endpoints if currency_node(metadata[n])}
    bridges = {n for n in currency if adjacency[n] and reverse[n]}
    metabolic_currency_edges = {k for k in keys if k[0] in bridges and
        metadata[k[2]]["node_type"] == "enzyme_reaction_instance" and
        any(e.get("mechanism") == "substrate_participation" for e in support[k])}
    cofactor_rows = []
    for n in sorted(currency):
        incoming = sorted(k for k in keys if k[2] == n)
        outgoing = sorted(k for k in keys if k[0] == n)
        is_product = any(e.get("mechanism") == "product_formation" for k in incoming for e in support[k])
        is_substrate = any(e.get("mechanism") == "substrate_participation" for k in outgoing for e in support[k])
        role = ("reaction_product_to_reaction_substrate_bridge" if n in bridges and is_product and is_substrate else
                "other_directed_bridge" if n in bridges else "incidental_product_sink" if incoming and not outgoing else
                "substrate_source_only" if outgoing and not incoming else "isolated")
        cofactor_rows.append(dict(node_id=n, label=label(n), chemical_id=metadata[n].get("chemical_id", ""),
            incoming_edges=len(incoming), outgoing_edges=len(outgoing), is_directed_bridge=n in bridges,
            classification=role, incoming_signed_edges_json=encoded(incoming), outgoing_signed_edges_json=encoded(outgoing)))
    edge_lookup = {edge_key(r): r for r in edge_rows}
    def routes(graph, category):
        paths = []
        for seed in upstream:
            distance, parents, bottleneck = shortest_tree(seed, graph, attributes)
            for target in downstream:
                if target not in distance:
                    continue
                nodes = [target]
                while nodes[-1] != seed:
                    nodes.append(parents[nodes[-1]])
                nodes.reverse()
                ee = [(a, pair_sign[a, b], b) for a, b in zip(nodes, nodes[1:])]
                activity = finite(attributes[seed].get("upstream_activity"), input_values.get("upstream_input.csv", {}).get(seed))
                weight, net = float(attributes[target]["score"]), math.prod(e[1] for e in ee)
                paths.append(dict(path_id=category + "_" + hashlib.sha256(encoded(nodes).encode()).hexdigest()[:16],
                    route_category=category, upstream_tf=seed, upstream_label=label(seed), downstream_node=target,
                    downstream_label=label(target), edge_count=distance[target], net_sign=net,
                    upstream_activity=activity, downstream_factor1_weight=weight,
                    input_sign_agreement=sign(activity)*net == sign(weight) if activity is not None else "not_assessed",
                    minimum_absolute_inferred_score=bottleneck[target],
                    node_score_levels_json=encoded([int(float(attributes[n]["level"])) for n in nodes]),
                    score_layer_aligned=all(edge_lookup[e]["score_layer_adjacent"] for e in ee),
                    edges_not_adjacent_in_score_layers=sum(not edge_lookup[e]["score_layer_adjacent"] for e in ee),
                    first_reaction_step=next((i for i, n in enumerate(nodes) if metadata[n]["node_type"] == "enzyme_reaction_instance"), ""),
                    first_species_rule_step=next((i for i, n in enumerate(nodes) if str(metadata[n].get("reaction_id", "")).startswith(("LIPID_RULE:", "EXPANDED_LIPID_RULE:"))), ""),
                    reaction_ids=joined(metadata[n].get("reaction_id") for n in nodes),
                    internal_currency_nodes=joined(n for n in nodes[1:-1] if n in currency),
                    contains_currency_bridge=any(n in bridges for n in nodes[1:-1]),
                    resources=joined(r.get("resource") for e in ee for r in support[e]),
                    node_ids="|".join(nodes), node_labels=" → ".join(label(n) for n in nodes), signed_edges_json=encoded(ee)))
        return paths
    no_currency = defaultdict(list)
    for s, i, t in keys - metabolic_currency_edges:
        no_currency[s].append(t)
    score_layers = defaultdict(list)
    for s, i, t in keys:
        if edge_lookup[s, i, t]["score_layer_adjacent"]:
            score_layers[s].append(t)
    catalyst = defaultdict(list)
    for s, i, t in sorted(keys):
        a, b = metadata[s]["node_type"], metadata[t]["node_type"]
        if (a in PROTEINS and b in PROTEINS or
            a in PROTEINS and b == "enzyme_reaction_instance" and s in metadata[t].get("catalysts", []) or
            a == "enzyme_reaction_instance" and t in downstream and
                any(e.get("mechanism") == "product_formation" for e in support[s, i, t])):
            catalyst[s].append(t)
    path_sets = {"all": routes(adjacency, "all"), "without_currency_bridges": routes(no_currency, "without_currency_bridges"),
                 "catalyst_only": routes(catalyst, "catalyst_only"), "score_layer_aligned": routes(score_layers, "score_layer_aligned")}
    steps = []
    for row in path_sets["all"]:
        for i, key in enumerate(json.loads(row["signed_edges_json"]), 1):
            steps.append(dict(path_id=row["path_id"], step=i, **edge_lookup[tuple(key)]))
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "network_labeled_ATT.csv", node_rows, list(dict.fromkeys(k for r in node_rows for k in r)) or ["node_id", "label"])
    write_csv(out / "network_labeled_SIF.csv", edge_rows, list(edge_rows[0]) if edge_rows else ["source", "interaction", "target"])
    graphml(out / "network.graphml", node_rows, edge_rows)
    write_jsonl(out / "network_edge_evidence.jsonl.gz", (r for k in sorted(keys) for r in support[k]))
    write_jsonl(out / "selected_reactions.jsonl.gz", (records[k] for k in sorted(records)))
    write_csv(out / "selected_reactions.csv", reaction_rows, list(reaction_rows[0]) if reaction_rows else ["reaction_id"])
    write_csv(out / "reaction_participant_completeness.csv", participants, list(participants[0]) if participants else ["reaction_id", "participant_node"])
    cofactor_fields = ["node_id", "label", "chemical_id", "incoming_edges", "outgoing_edges", "is_directed_bridge", "classification", "incoming_signed_edges_json", "outgoing_signed_edges_json"]
    write_csv(out / "generic_cofactor_bridge_audit.csv", cofactor_rows, cofactor_fields)
    path_fields = ["path_id", "route_category", "upstream_tf", "upstream_label", "downstream_node", "downstream_label", "edge_count", "net_sign", "upstream_activity", "downstream_factor1_weight", "input_sign_agreement", "minimum_absolute_inferred_score", "node_score_levels_json", "score_layer_aligned", "edges_not_adjacent_in_score_layers", "first_reaction_step", "first_species_rule_step", "reaction_ids", "internal_currency_nodes", "contains_currency_bridge", "resources", "node_ids", "node_labels", "signed_edges_json"]
    for name, rr in path_sets.items():
        write_csv(out / ("candidate_paths_" + name + ".csv"), rr, path_fields)
        write_jsonl(out / ("candidate_paths_" + name + ".jsonl"), rr)
    write_csv(out / "candidate_path_steps.csv", steps, ["path_id", "step"] + (list(edge_rows[0]) if edge_rows else ["source", "interaction", "target"]))
    hashes = [dict(path=str(p.resolve()), bytes=p.stat().st_size, sha256=sha256(p)) for p in dependencies]
    write_csv(out / "dependency_hashes.csv", hashes, ["path", "bytes", "sha256"])
    summary = dict(solution_edges=len(keys), solution_nodes=len(endpoints), selected_upstream_tfs=len(upstream),
        selected_observed_lipids=len(downstream), exact_evidence_records=sum(map(len, support.values())),
        complete_master_reactions=len(records), selected_reaction_instances=sum(map(len, selected_reactions.values())),
        participant_rows=len(participants), participant_connectors_missing_from_full=sum(not r["connector_in_full_pkn"] for r in participants),
        participant_connectors_missing_from_solution=sum(not r["connector_in_solution"] for r in participants),
        final_pkn_completeness_assessed=pruned is not None,
        currency_nodes=len(currency), currency_bridges=len(bridges), currency_role_counts=dict(Counter(r["classification"] for r in cofactor_rows)),
        route_counts={name: len(rr) for name, rr in path_sets.items()},
        observed_lipids_with_routes={name: len({r["downstream_node"] for r in rr}) for name, rr in path_sets.items()},
        all_representative_routes_with_currency_bridges=sum(r["contains_currency_bridge"] for r in path_sets["all"]),
        all_representative_routes_score_layer_aligned=sum(r["score_layer_aligned"] for r in path_sets["all"]),
        solution_edges_adjacent_in_score_layers=sum(r["score_layer_adjacent"] for r in edge_rows),
        independent_checks=dict(exact_pkn_membership=True, exact_eligible_signed_provenance=True,
            finite_unique_score_ids=True, edge_score_sign_consistency=True, all_nodes_on_input_to_observation_walk=True,
            downstream_observations_are_terminal=True),
        scoring_performed=False, inputs_or_fit_or_pkn_modified=False,
        route_policy="One shortest path per reachable TF/observation pair; ties maximize inferred-score bottleneck, then deterministic predecessor ID. Level-0 weights excluded from bottleneck.",
        currency_audit_policy="Explicit chemical-ID/label panel of water/proton, CoA, ATP/ADP/AMP/GTP/GDP, NAD(P)(H), phosphate/diphosphate, SAM/SAH, including H+/Pi/PPi aliases. A bridge requires incoming and outgoing selected edges. The currency-free route sensitivity removes only outgoing substrate_participation to reaction-instance edges and retains chemical-to-protein regulation. This is an audit panel, not exhaustive biochemical classification.",
        participant_policy="Each oriented master participant and selected catalyst has its original +1 connector checked separately in full PKN, optional final scoring PKN, and solution. Presence is bookkeeping; no joint availability/stoichiometric feasibility inferred.",
        score_layer_audit_policy="An edge can directly carry recursive scoring evidence only when source level equals target level plus one. The score-layer-aligned route audit keeps only these edges; sign-consistent extracted topology may additionally contain non-adjacent or backward-layer edges. This audit reuses fixed scores and does not quantify each child's contribution.",
        limitations=["No flux, kinetics, concentration control, joint co-substrate availability or atom balance is inferred.",
            "Generated species/catalyst alternatives are hypotheses, not independent observations.",
            "Shortest representative paths do not enumerate every alternative path/walk.",
            "The currency-free and catalyst-only route audits reuse fixed scores and do not refit MOON."])
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
