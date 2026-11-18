#!/usr/bin/env python3
"""Compare normalized curated and rule-expanded lipid-to-TF topology, not scores.

Run from any working directory. Writes only alongside this script. The original
PKN files and the 155 accepted input identifiers are never modified. Chemical
normalization uses exactly one lookup per endpoint in chemical_node_aliases.csv.
"""
import collections
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PKN = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
INPUTS = ROOT / "results/causal_modelling/cosmos_factor1_pilot/inputs"
DEPTHS = (6, 8, 12, None)


def read_csv(path, delimiter=","):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(name, rows, fields):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_edges(path):
    result = [(r["source"], int(r["interaction"]), r["target"])
              for r in read_csv(path)]
    assert all(s and t and sign in (-1, 1) for s, sign, t in result)
    return result


def chemical(node):
    return node[len("Metab__"):].rsplit("_", 1)[0] if node.startswith("Metab__") else ""


def protein(node):
    return node.startswith(("Gene__", "UniProt__"))


def reaction(node):
    return node.startswith("Enzyme__")


def node_set(edges):
    return {node for source, _, target in edges for node in (source, target)}


def adjacency(edges, terminal):
    result = collections.defaultdict(set)
    for source, _, target in edges:
        if source not in terminal:
            result[source].add(target)
    return result


def bfs(adj, starts):
    distance = {s: 0 for s in starts}
    queue = collections.deque(sorted(starts))
    while queue:
        source = queue.popleft()
        for target in sorted(adj.get(source, ())):
            if target not in distance:
                distance[target] = distance[source] + 1
                queue.append(target)
    return distance


def main():
    dependencies = [PKN / name for name in (
        "pkn_lipid_rules.csv", "pkn_lipid_curated.csv", "pkn_lipid_curated_ph73.csv",
        "chemical_node_aliases.csv", "rhea_lipid_identity_classification.tsv",
        "regulation_chebi_source_crosswalk.csv",
    )] + [INPUTS / name for name in (
        "factor1_lipid_node_mapping.csv", "factor1_hepatocyte_tf_activities_mapped.csv",
    )] + [Path(__file__).resolve()]
    alias_rows = read_csv(PKN / "chemical_node_aliases.csv")
    aliases = {r["original_node"]: r["normalized_node"] for r in alias_rows}
    assert len(aliases) == len(alias_rows), "Ambiguous/repeated alias source"
    assert all(chemical(s) and chemical(t) for s, t in aliases.items())
    alias_chain_count = sum(aliases.get(t, t) != t for t in aliases.values())

    accepted_rows = [r for r in read_csv(INPUTS / "factor1_lipid_node_mapping.csv")
                     if r["included"] == "TRUE"]
    accepted = {r["node_id"]: r["feature"] for r in accepted_rows}
    assert len(accepted) == len(accepted_rows) == 155
    changed_inputs = [(node, aliases[node]) for node in accepted
                      if aliases.get(node, node) != node]
    # Explicitly do not change measurement IDs, even if future alias files do.
    assert not changed_inputs, "An accepted input would change; do not silently remap it"
    tfs = {r["node_id"]: float(r["oriented_activity"])
           for r in read_csv(INPUTS / "factor1_hepatocyte_tf_activities_mapped.csv")
           if r["in_pkn"] == "TRUE"}
    terminal_sets = {
        "all_mapped_tfs": set(tfs),
        "abs_activity_gt_2": {node for node, value in tfs.items() if abs(value) > 2},
    }
    assert len(tfs) == 353 and len(terminal_sets["abs_activity_gt_2"]) == 49
    lipid_ids = {r["id"] for r in read_csv(PKN / "rhea_lipid_identity_classification.tsv", "\t")}
    lipid_ids |= {r["source"] for r in read_csv(PKN / "regulation_chebi_source_crosswalk.csv")
                  if r["lipid_is_a"] == "True"}

    def lipid(node):
        return node.startswith("Metab__MeasuredLipid:") or chemical(node) in lipid_ids

    def precursor(edges):
        return {e for e in edges if
                (protein(e[0]) and protein(e[2])) or
                (lipid(e[0]) and protein(e[2])) or
                (lipid(e[0]) and reaction(e[2])) or
                (reaction(e[0]) and lipid(e[2]))}

    clean_records, conflicts = [], []

    def cleanup(name, edges, project):
        projected = [(aliases.get(s, s), sign, aliases.get(t, t)) for s, sign, t in edges] if project else edges
        distinct = set(projected)
        selfloops = {e for e in distinct if e[0] == e[2]}
        pairs = collections.defaultdict(set)
        for source, sign, target in distinct - selfloops:
            pairs[source, target].add(sign)
        conflicting = {pair for pair, signs in pairs.items() if len(signs) > 1}
        for source, target in sorted(conflicting):
            conflicts.append({"graph": name, "source": source, "target": target, "signs": "-1|1"})
        clean_records.append({"graph": name, "input_signed_edges": len(edges),
            "one_step_projection": project,
            "changed_endpoint_occurrences": sum(s != ns or t != nt for (s, _, t), (ns, _, nt) in zip(edges, projected)),
            "duplicate_signed_edges_removed": len(projected) - len(distinct),
            "self_loop_signed_edges_removed": len(selfloops),
            "conflicting_sign_pairs": len(conflicting),
            "clean_signed_edges": len(distinct - selfloops)})
        return distinct - selfloops

    projected_curated = cleanup("projected_curated", read_edges(PKN / "pkn_lipid_curated.csv"), True)
    prebuilt_curated = cleanup("prebuilt_curated_ph73", read_edges(PKN / "pkn_lipid_curated_ph73.csv"), False)
    projected_rules = cleanup("projected_rules_ph73", read_edges(PKN / "pkn_lipid_rules.csv"), True)
    write_csv("cleanup_audit.csv", clean_records, list(clean_records[0]))
    write_csv("conflicting_signed_pairs.csv", conflicts, ["graph", "source", "target", "signs"])
    assert not conflicts, "Opposing signs require an explicit decision; no reachability results written"
    assert projected_curated == prebuilt_curated, "One-step projected curated graph differs from the prebuilt graph"
    assert prebuilt_curated <= projected_rules, "Rule projection lost a curated edge"
    graphs = {
        "prebuilt_curated_ph73": precursor(prebuilt_curated),
        "projected_rules_ph73": precursor(projected_rules),
    }
    added = graphs["projected_rules_ph73"] - graphs["prebuilt_curated_ph73"]
    assert not (graphs["prebuilt_curated_ph73"] - graphs["projected_rules_ph73"])
    edge_fields = ["source", "interaction", "target"]
    write_csv("added_rule_precursor_edges.csv", [dict(zip(edge_fields, e)) for e in sorted(added)], edge_fields)
    write_csv("pkn_lipid_rules_ph73_precursor.csv",
              [dict(zip(edge_fields, e)) for e in sorted(graphs["projected_rules_ph73"])], edge_fields)

    graph_counts, reaches, summaries, comparisons, input_status = [], [], [], [], []
    all_pairs, corridors = {}, {}
    for name, graph in graphs.items():
        nodes = node_set(graph)
        graph_counts.append({"graph": name, "signed_edges": len(graph), "nodes": len(nodes),
            "accepted_inputs_present": len(nodes & accepted.keys()),
            "accepted_inputs_absent": len(accepted.keys() - nodes),
            "protein_to_reaction_edges": sum(protein(s) and reaction(t) for s, _, t in graph),
            "nonlipid_chemical_nodes": sum(bool(chemical(n)) and not lipid(n) for n in nodes)})
        assert graph_counts[-1]["protein_to_reaction_edges"] == 0
        assert graph_counts[-1]["nonlipid_chemical_nodes"] == 0
        for tf_name, terminal in terminal_sets.items():
            adj = adjacency(graph, terminal)
            pairs = {}
            for source, feature in sorted(accepted.items()):
                distance = bfs(adj, [source]) if source in nodes else {}
                targets = terminal & distance.keys()
                proteins_reached = {n for n in distance if protein(n)}
                input_status.append({"graph": name, "tf_endpoint_set": tf_name,
                    "feature": feature, "source": source, "present": source in nodes,
                    "reachable_proteins": len(proteins_reached), "reachable_tf_endpoints": len(targets)})
                for target in sorted(targets):
                    pairs[source, target] = distance[target]
                    reaches.append({"graph": name, "tf_endpoint_set": tf_name,
                        "feature": feature, "source": source, "target_tf": target,
                        "shortest_path_edges": distance[target]})
            all_pairs[name, tf_name] = pairs
            for depth in DEPTHS:
                selected = {pair: distance for pair, distance in pairs.items()
                            if depth is None or distance <= depth}
                summaries.append({"graph": name, "tf_endpoint_set": tf_name,
                    "observed_tfs": len(terminal), "maximum_shortest_path_edges": depth or "unbounded",
                    "reachable_input_lipids": len({s for s, _ in selected}),
                    "reachable_tf_endpoints": len({t for _, t in selected}), "reachable_pairs": len(selected)})
            # This is a path corridor, not MOON's iterative pruning algorithm.
            forward = set(bfs(adj, accepted.keys() & nodes))
            reverse = collections.defaultdict(set)
            for source, targets in adj.items():
                for target in targets:
                    reverse[target].add(source)
            backward = set(bfs(reverse, terminal & nodes))
            corridor = {(s, sign, t) for s, sign, t in graph
                        if s not in terminal and s in forward and t in backward}
            corridors[name, tf_name] = corridor

    for tf_name in terminal_sets:
        base = all_pairs["prebuilt_curated_ph73", tf_name]
        full = all_pairs["projected_rules_ph73", tf_name]
        base_corridor = corridors["prebuilt_curated_ph73", tf_name]
        full_corridor = corridors["projected_rules_ph73", tf_name]
        comparisons.append({"tf_endpoint_set": tf_name,
            "added_endpoint_pairs": len(full.keys() - base.keys()),
            "lost_endpoint_pairs": len(base.keys() - full.keys()),
            "changed_shortest_distances": sum(base[p] != full[p] for p in base.keys() & full.keys()),
            "added_reachable_input_candidates": len({s for s, _ in full} - {s for s, _ in base}),
            "base_corridor_signed_edges": len(base_corridor), "rules_corridor_signed_edges": len(full_corridor),
            "added_corridor_signed_edges": len(full_corridor - base_corridor),
            "lost_corridor_signed_edges": len(base_corridor - full_corridor)})

    write_csv("graph_counts.csv", graph_counts, list(graph_counts[0]))
    write_csv("topology_summaries.csv", summaries, list(summaries[0]))
    write_csv("endpoint_reachability.csv", reaches, list(reaches[0]))
    write_csv("input_reachability.csv", input_status, list(input_status[0]))
    write_csv("comparison.csv", comparisons, list(comparisons[0]))
    manifest = [{"path": str(p.relative_to(ROOT)), "bytes": p.stat().st_size,
                 "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in dependencies]
    write_csv("dependency_hashes.csv", manifest, ["path", "bytes", "sha256"])
    report = {
        "analysis": "one_step_normalized_rules_lipid_to_tf_topology",
        "scientific_refit_performed": False, "source_pkn_modified": False,
        "measurement_mapping_modified": False, "accepted_input_nodes": len(accepted),
        "accepted_input_ids_changed_by_alias_map": len(changed_inputs),
        "nontrivial_alias_chains": alias_chain_count,
        "projected_curated_equals_prebuilt_curated": projected_curated == prebuilt_curated,
        "added_normalized_rule_edges_before_precursor_filter": len(projected_rules - prebuilt_curated),
        "added_normalized_rule_edges_after_precursor_filter": len(added),
        "graph_counts": graph_counts, "cleanup": clean_records, "comparisons": comparisons,
        "reachable_input_features": sorted({accepted[s] for pairs in all_pairs.values() for s, _ in pairs}),
        "semantics": "Directed endpoint shortest paths after making the specified TF observations terminal; unsigned reachability, not activity or significance. The additional corridor comparison is not iterative MOON pruning.",
        "policy": "Exactly the lipid precursor categories from audit_lipid_tf_coverage.py: protein->protein, lipid->protein, lipid->reaction, reaction->lipid. No protein->reaction or nonlipid chemicals. Explicit one-step node alias projection; no measured input ID changes or new biological edges.",
        "conflict_policy": "Deduplicate identical signed edges and remove selfloops; if opposing sign pairs exist, report them and refuse reachability analysis. No conflicts or selfloops occurred for these inputs.",
    }
    (HERE / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    (HERE / "README.md").write_text(
        "# Preserving measured-species rules during chemical normalization\n\n"
        "Run `python3 results/causal_modelling/factor1_lipid_to_tf/normalized_rules_topology/audit_normalized_rules_topology.py` "
        "from the project root. The script writes only in this directory and hashes its frozen dependencies. "
        "No original PKN, fit, or accepted measurement mapping is changed.\n\n"
        "The original rule-expanded graph is projected through exactly one frozen `chemical_node_aliases.csv` lookup "
        "per endpoint. Projecting the original curated graph reproduces the prebuilt pH7.3 curated graph exactly. "
        "No accepted input ID is changed. The lipid-only precursor policy is then applied.\n\n"
        "The normalized rules restore all 155 accepted input nodes, versus six in the curated graph, and add 898 "
        "precursor edges. They do not change any reachable lipid–TF pair, shortest distance, upstream candidate, "
        "or signed edge on an input-to-terminal-TF path corridor. Only LPC(16:0/0:0) and LPC(18:0/0:0) reach observed TFs. "
        "The other 149 rule-mapped inputs remain without a path to a protein under this policy. "
        "Thus omission of the measured-species rules affects input-node presence but does not explain the limited "
        "upstream signaling candidates in this particular topology. This audit is not a MOON refit and does not "
        "claim equality of scores from independently rerun pipelines.\n\n"
        "`cleanup_audit.csv` and `conflicting_signed_pairs.csv` report deduplication, selfloops and opposing signs. "
        "`comparison.csv` gives the exact pair/distance/corridor comparison for all 353 mapped TFs and the 49 TFs "
        "with absolute activity above two, made terminal separately. `input_reachability.csv` preserves all 155 "
        "inputs under both policies. `endpoint_reachability.csv` records each reachable pair; "
        "`topology_summaries.csv` also reports distance limits 6, 8, 12 and unbounded. These are topological "
        "distances, not iterative MOON pruning or scoring depth.\n\n"
        "`pkn_lipid_rules_ph73_precursor.csv` is a derived diagnostic copy of the normalized, filtered "
        "rule-expanded graph, not a replacement for the frozen PKN. `added_rule_precursor_edges.csv` isolates "
        "its additional edges relative to the normalized curated precursor graph.\n"
    )
    print(json.dumps({"graph_counts": graph_counts, "comparisons": comparisons,
                      "added_rule_precursor_edges": len(added)}, indent=2))


if __name__ == "__main__":
    main()
