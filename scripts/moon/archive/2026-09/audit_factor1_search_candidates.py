#!/usr/bin/env python3
"""Compare saved refined Factor 1 candidates; never score or edit search grids."""
import argparse
import collections
import csv
import hashlib
import io
import itertools
import json
import math
from pathlib import Path
import statistics
import subprocess

from export_factor1_network import edge_key, finite, read_csv, sign, write_csv


CANDIDATES = [
    ("U4", "radius6_depth4_uncompressed", "solution_p1_90_s1_60", 1.90, 1.60),
    ("U5", "radius6_depth5_uncompressed", "solution_p2_05_s1_60", 2.05, 1.60),
    ("C4", "radius6_depth4_compressed", "solution_p2_00_s0_85", 2.00, 0.85),
]


def quantile(values, probability):
    values = sorted(values)
    if not values:
        return None
    index = (len(values) - 1) * probability
    low, high = math.floor(index), math.ceil(index)
    return values[low] + (values[high] - values[low]) * (index - low)


def distribution(values):
    values = [value for value in values if value is not None and math.isfinite(value)]
    return {"n": len(values), "negative": sum(value < 0 for value in values), "positive": sum(value > 0 for value in values),
            "minimum": min(values) if values else None, "q25": quantile(values, 0.25), "median": quantile(values, 0.5),
            "q75": quantile(values, 0.75), "maximum": max(values) if values else None,
            "mean": statistics.mean(values) if values else None, "sd": statistics.stdev(values) if len(values) > 1 else None}


def pearson(x, y):
    if len(x) < 2:
        return None
    mx, my = statistics.mean(x), statistics.mean(y)
    numerator = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
    return numerator / denominator if denominator else None


def ranks(values):
    result = [None] * len(values)
    groups = collections.defaultdict(list)
    for index, value in enumerate(values):
        groups[value].append(index)
    start = 1
    for value in sorted(groups):
        indices = groups[value]
        for index in indices:
            result[index] = start + (len(indices) - 1) / 2
        start += len(indices)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", type=Path, default=Path("results/cosmos_factor1_parameter_search"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.search / "refined_candidate_audit"
    output.mkdir(parents=True, exist_ok=True)
    datasets, summary, endpoints, score_rows, score_distributions, support_rows, paths, source_hashes = {}, [], [], [], [], [], [], []
    downstream = {row["source"]: float(row["score"]) for row in read_csv(args.search / "topology_radius_6/downstream_input.csv")}
    catalog = {row["node_id"]: row for row in read_csv(Path("results/cosmos_mouse_pkn/nodes.csv"))}
    for alias, run, stem, primary, secondary in CANDIDATES:
        directory = args.search / "runs" / run
        att_path, sif_path = directory / (stem + "_ATT.csv"), directory / (stem + "_SIF.csv")
        attributes = {row["source"]: row for row in read_csv(att_path)}
        edges = {edge_key(row) for row in read_csv(sif_path)}
        all_scores = {row["source"]: row for row in read_csv(directory / "moon_scores.csv")}
        upstream = {node for node, row in attributes.items() if row["type"] == "upstream_input"}
        lipids = {node for node, row in attributes.items() if float(row["level"]) == 0}
        manifest = json.loads((directory / "run_manifest.json").read_text())
        audit = json.loads((directory / ("network_evidence_" + stem + "_summary.json")).read_text())
        if audit["inconsistent_scored_edges"] or audit["nodes_unreachable_from_intended_upstream"] or audit["nodes_unable_to_reach_observation"]:
            raise ValueError("Candidate did not pass independent solution audit: " + run)
        saved_graph = subprocess.run(["Rscript", str(Path(__file__).with_name("read_saved_candidate_pkn.R")),
                                      str(directory / "pkn_final_uncompressed.rds"), str(att_path)], capture_output=True, text=True, check=True)
        graph_rows = list(csv.DictReader(io.StringIO(saved_graph.stdout)))
        children = collections.defaultdict(list)
        for row in graph_rows:
            children[row["source"]].append(row)
        for filename in [att_path, sif_path, directory / "moon_scores.csv", directory / "run_manifest.json", directory / "pkn_final_uncompressed.rds"]:
            source_hashes.append({"candidate": alias, "file": str(filename), "sha256": hashlib.sha256(filename.read_bytes()).hexdigest()})
        path_rows = read_csv(directory / ("candidate_paths_" + stem + ".csv"))
        paths.extend(dict(row, candidate=alias) for row in path_rows)
        inferred_values = [float(row["score"]) for row in attributes.values() if float(row["level"]) > 0]
        minimum_margin = min(abs(value) - secondary for value in inferred_values)
        summary.append({"candidate": alias, "run_id": run, "solution_stem": stem, "compression": manifest["compression"],
                        "n_layers": manifest["n_layers"], "topology_radius": manifest["topology_radius"], "primary_threshold": primary, "secondary_threshold": secondary,
                        "nodes": len(attributes), "edges": len(edges), "upstream_tfs": len(upstream), "downstream_lipids": len(lipids),
                        "tf_identities": "|".join(sorted(upstream)), "lipid_identities": "|".join(sorted(lipids)),
                        "route_pairs": len(path_rows), "shortest_route_min_edges": min(int(row["edge_count"]) for row in path_rows),
                        "shortest_route_max_edges": max(int(row["edge_count"]) for row in path_rows),
                        "route_pairs_at_most_6": sum(int(row["edge_count"]) <= 6 for row in path_rows),
                        "minimum_inferred_secondary_margin": minimum_margin,
                        "inferred_nodes_within_0_05_of_secondary": sum(abs(value) - secondary < 0.05 for value in inferred_values),
                        "coherence_iterations": manifest["iterations"], "scored_nodes_in_full_run": len(all_scores),
                        "evidence_records": audit["edge_evidence_records"], "species_rule_edges": audit["edges_with_species_rule_hypotheses"],
                        "graphml_path": str(directory / (stem + ".graphml"))})
        for node in sorted(upstream | lipids):
            row = attributes[node]
            endpoints.append({"candidate": alias, "node_id": node, "label": row.get("label", node), "role": "upstream_tf" if node in upstream else "measured_lipid",
                              "inferred_moon_score": float(row["score"]) if node in upstream else None,
                              "tf_activity_input": finite(row.get("upstream_activity")), "factor1_weight": float(row["score"]) if node in lipids else None,
                              "secondary_margin": abs(float(row["score"])) - secondary if node in upstream else None})
        for node, row in attributes.items():
            score_rows.append({"candidate": alias, "node_id": node, "label": row.get("label", node), "level": float(row["level"]), "role": row["type"],
                               "score": float(row["score"]), "score_type": "Factor1_weight" if float(row["level"]) == 0 else "inferred_MOON",
                               "secondary_margin": abs(float(row["score"])) - secondary if float(row["level"]) > 0 else None})
        for scope, values in [
            ("candidate_inferred_signed", inferred_values),
            ("candidate_inferred_absolute", [abs(value) for value in inferred_values]),
            ("candidate_upstream_inferred", [float(attributes[node]["score"]) for node in upstream]),
            ("candidate_lipid_weights", [float(attributes[node]["score"]) for node in lipids]),
            ("full_run_inferred_signed", [float(row["score"]) for row in all_scores.values() if float(row["level"]) > 0])]:
            score_distributions.append(dict(distribution(values), candidate=alias, scope=scope))
        for level in sorted({float(row["level"]) for row in attributes.values()}):
            score_distributions.append(dict(distribution([float(row["score"]) for row in attributes.values() if float(row["level"]) == level]), candidate=alias, scope="candidate_level_" + str(int(level))))
        for node, row in attributes.items():
            if not node.startswith("Enzyme__"):
                continue
            measured_edges = [edge for edge in children[node] if edge["target"] in downstream]
            coefficients = {edge["target"]: int(edge["interaction"]) for edge in measured_edges}
            observed_ids = sorted(coefficients)
            score = float(row["score"])
            support_rows.append({"candidate": alias, "enzyme_node": node, "reaction_id": row.get("reaction_id", ""), "catalysts": row.get("catalysts", ""),
                                 "inferred_score": score, "level": float(row["level"]), "observed_target_count": len(observed_ids),
                                 "observed_targets": "|".join(observed_ids), "observation_background_n": len(downstream),
                                 "children_in_final_original_pkn": len({edge["target"] for edge in children[node]}),
                                 "single_observation_footprint": len(observed_ids) == 1,
                                 "interpretation": "Immediate observed support in saved final original PKN; target count is not independent replication or substrate-specific validation"})
        datasets[alias] = {"attributes": attributes, "edges": edges, "nodes": set(attributes), "upstream": upstream, "lipids": lipids, "all_scores": all_scores}

    pairwise, paired_scores = [], []
    for left, right in itertools.combinations(datasets, 2):
        a, b = datasets[left], datasets[right]
        result = {"left": left, "right": right}
        for name in ["nodes", "edges", "upstream", "lipids"]:
            shared, union = a[name] & b[name], a[name] | b[name]
            result[name + "_shared"] = len(shared)
            result[name + "_union"] = len(union)
            result[name + "_jaccard"] = len(shared) / len(union) if union else None
        for scope, scores_a, scores_b in [("candidate", a["attributes"], b["attributes"]), ("full_scored", a["all_scores"], b["all_scores"])]:
            common = sorted(set(scores_a) & set(scores_b))
            common = [node for node in common if float(scores_a[node]["level"]) > 0 and float(scores_b[node]["level"]) > 0]
            x = [float(scores_a[node]["score"]) for node in common]
            y = [float(scores_b[node]["score"]) for node in common]
            result[scope + "_shared_inferred_nodes"] = len(common)
            result[scope + "_score_pearson"] = pearson(x, y)
            result[scope + "_score_spearman"] = pearson(ranks(x), ranks(y))
            result[scope + "_sign_disagreements"] = sum(sign(v1) != sign(v2) for v1, v2 in zip(x, y))
            result[scope + "_median_absolute_score_change"] = statistics.median(abs(v1 - v2) for v1, v2 in zip(x, y)) if common else None
            result[scope + "_maximum_absolute_score_change"] = max((abs(v1 - v2) for v1, v2 in zip(x, y)), default=None)
            result[scope + "_level_disagreements"] = sum(scores_a[node]["level"] != scores_b[node]["level"] for node in common)
            if scope == "candidate":
                for node in common:
                    paired_scores.append({"left": left, "right": right, "node_id": node, "left_score": float(scores_a[node]["score"]),
                                          "right_score": float(scores_b[node]["score"]), "left_level": float(scores_a[node]["level"]), "right_level": float(scores_b[node]["level"])})
        pairwise.append(result)
    shared_nodes = set.intersection(*(data["nodes"] for data in datasets.values()))
    shared_edges = set.intersection(*(data["edges"] for data in datasets.values()))
    union_nodes = set.union(*(data["nodes"] for data in datasets.values()))
    membership = []
    for node in sorted(union_nodes):
        row = {"node_id": node, "label": catalog.get(node, {}).get("label", node)}
        for alias, data in datasets.items():
            row[alias + "_in_solution"] = node in data["nodes"]
            row[alias + "_is_upstream_endpoint"] = node in data["upstream"]
            score = data["all_scores"].get(node, {})
            row[alias + "_score"] = finite(score.get("score"))
            row[alias + "_level"] = finite(score.get("level"))
        membership.append(row)
    tables = {"candidate_summary.csv": summary, "endpoint_identities.csv": endpoints, "pairwise_overlap.csv": pairwise,
              "score_distributions.csv": score_distributions, "candidate_node_scores.csv": score_rows, "paired_candidate_scores.csv": paired_scores,
              "reaction_observation_support.csv": support_rows, "candidate_paths.csv": paths, "node_membership_and_scores.csv": membership,
              "source_artifact_hashes.csv": source_hashes}
    for name, records in tables.items():
        fields = list(dict.fromkeys(key for row in records for key in row))
        write_csv(output / name, records, fields)
    write_csv(output / "shared_signed_edges_all_three.csv", [{"source": a, "interaction": s, "target": b} for a, s, b in sorted(shared_edges)], ["source", "interaction", "target"])
    audit_manifest = {"candidates": summary, "shared_nodes_all_three": len(shared_nodes), "shared_signed_edges_all_three": len(shared_edges),
                      "shared_lipid_endpoints_all_three": sorted(set.intersection(*(data["lipids"] for data in datasets.values()))),
                      "shared_upstream_TF_endpoints_all_three": sorted(set.intersection(*(data["upstream"] for data in datasets.values()))),
                      "scope": "Read-only comparison of saved networks, scores, PKNs and provenance; no rescoring or threshold-grid changes",
                      "caution": "Candidates differ in propagation depth/compression and primary threshold, so comparisons do not isolate a single factor except where explicitly stated."}
    (output / "comparison_manifest.json").write_text(json.dumps(audit_manifest, indent=2) + "\n")
    lines = ["# Refined Factor 1 candidate comparison", "", "All three candidates pass the independent exact-PKN provenance, edge-sign and upstream/downstream connectivity audits. They retain the same two observed DG endpoints but different upstream TF endpoints. This is stability of a downstream reaction hypothesis, not stable identification of its upstream TF control.", "", "| Candidate | Scoring | Primary / secondary | Nodes / edges | Upstream TFs | Shortest routes |", "|---|---|---:|---:|---|---:|"]
    for row in summary:
        lines.append(f"| {row['candidate']} | {row['compression']}, depth {row['n_layers']} | {row['primary_threshold']:.2f} / {row['secondary_threshold']:.2f} | {row['nodes']} / {row['edges']} | {row['tf_identities'].replace('Gene__', '').replace('|', ', ')} | {row['shortest_route_min_edges']}-{row['shortest_route_max_edges']} edges |")
    lines += ["", "## Shared terminal hypothesis", "", "The retained observations are DG(24:0_18:2) and DG(24:1_18:2). Both terminate in Smpd1-associated, rule-instantiated PC-to-DG reactions supported by the generic RHEA:10605 template. The exact PC/DG acyl-chain specificity remains assumed. Each retained enzyme-instance score has only one directly linked measured observation in a background of 74 measured nodes; the two separate enzyme instances do not provide independent validation of the reaction rules.", "", "All three candidate graphs contain the Pdpk1 → Prkcd → Smpd1 → reaction-instance → DG tail. The complete selected shortest TF-to-DG routes are 7–11 edges long. None has a route of at most six edges, even though the score propagation depths are four or five; depth and route length are distinct quantities.", "", "## What changes between candidates", ""]
    for row in pairwise:
        lines.append(f"- {row['left']} versus {row['right']}: {row['nodes_shared']} shared nodes (Jaccard {row['nodes_jaccard']:.3f}), {row['edges_shared']} shared signed edges (Jaccard {row['edges_jaccard']:.3f}); TF endpoint Jaccard {row['upstream_jaccard']:.3f}.")
    depth_pair = next(row for row in pairwise if row["left"] == "U4" and row["right"] == "U5")
    depth_statement = (f"All {depth_pair['full_scored_shared_inferred_nodes']:,} inferred nodes shared by the saved U4 and U5 score tables have identical scores and levels. U5 adds scored nodes and uses a different primary cutoff; the TF endpoint change therefore cannot be described as changed scores on the shared nodes or attributed to depth alone." if depth_pair["full_scored_maximum_absolute_score_change"] == 0 and depth_pair["full_scored_level_disagreements"] == 0 else "The shared-score table records numerical and level differences explicitly; the selected candidates also use different primary cutoffs.")
    lines += ["", f"The three-way intersection contains {len(shared_nodes)} nodes and {len(shared_edges)} signed edges, but no common upstream TF endpoint. Score-distribution and shared-node score comparisons are exported separately; they do not imply numerical equivalence between differently compressed or propagated runs.", "", depth_statement, "", "U5 meets the higher primary-cutoff preference, whereas U4 is smaller and retains the same two lipid endpoints. Klf10 in U5 lies only 0.000186 above its secondary cutoff in absolute score. This threshold margin is a numerical sensitivity warning, not a confidence interval. U4's Nfe2l2 score also lies close to the 1.60 secondary cutoff. The compressed candidate changes the upstream TF endpoints completely and requires a much lower secondary cutoff.", "", "## Recommended interpretation", "", "Keep U4 and U5 as distinct uncompressed candidate hypotheses and C4 as a compression-sensitivity comparison. U4 is a useful compact network for inspection; U5 satisfies the stronger primary-cutoff preference. Neither is validated as the upstream mechanism. The common Smpd1-associated PC-to-DG tail is the reproducible part of these three extracted graphs, with measured-species reaction specificity still hypothetical. Treat the TF assignments and long regulatory routes as candidates requiring stronger context and perturbation evidence.", "", "## Auditable outputs", "", "- `candidate_summary.csv`, `endpoint_identities.csv`, `pairwise_overlap.csv`: exact identities, sizes, path bounds and signed overlaps.", "- `candidate_paths.csv`, `reaction_observation_support.csv`: auditable paths and immediate measured support.", "- `score_distributions.csv`, `candidate_node_scores.csv`, `paired_candidate_scores.csv`, `node_membership_and_scores.csv`: score scales, margins, memberships and comparisons.", "- `shared_signed_edges_all_three.csv`: three-way overlap; it is an intersection for comparison, not a newly fitted network.", "- `source_artifact_hashes.csv`: SHA-256 hashes of the saved source artifacts.", "", "GraphML and complete edge/reaction provenance are alongside each original SIF/ATT, using its exact solution stem. Original scores, grids and networks were not changed.", ""]
    (output / "refined_candidate_comparison.md").write_text("\n".join(lines))
    print(json.dumps(audit_manifest, indent=2))


if __name__ == "__main__":
    main()
