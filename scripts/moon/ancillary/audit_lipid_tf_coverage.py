#!/usr/bin/env python3
"""Audit measured lipid -> TF coverage without fitting MOON or changing a PKN.

All reachability counts are directed shortest-path endpoint counts, with the
specified observed TF set made terminal. They are not scores, significance tests,
or the output of MOON's alternating, iterative forward/backward pruning.
Only Python's standard library and frozen project resources are used.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
CURRENCY = {
    "CHEBI:15378": "H(+) / proton", "CHEBI:15377": "water",
    "CHEBI:18367": "inorganic phosphate", "CHEBI:33019": "diphosphate",
    "CHEBI:15379": "oxygen", "CHEBI:10545": "electron",
}
DEPTHS = (6, 8, 12, None)


def read_csv(path, delimiter=","):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_edges(path):
    rows = read_csv(path)
    edges = [(r["source"], int(r["interaction"]), r["target"]) for r in rows]
    assert all(s and t and sign in (-1, 1) for s, sign, t in edges)
    assert len(edges) == len(set(edges)), f"Duplicate signed edges: {path}"
    return edges


def chemical(node):
    return node[len("Metab__"):].rsplit("_", 1)[0] if node.startswith("Metab__") else ""


def protein(node):
    return node.startswith(("Gene__", "UniProt__"))


def reaction(node):
    return node.startswith("Enzyme__")


def nodes(edges):
    return {n for s, _, t in edges for n in (s, t)}


class ChemicalEquivalence:
    """Only explicit ChEBI alt_id and conjugate acid/base relationships."""
    def __init__(self):
        self.parent = {}

    def find(self, value):
        self.parent.setdefault(value, value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def join(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def adjacency(edges, terminal=()):
    terminal = set(terminal)
    result = collections.defaultdict(set)
    for source, _, target in edges:
        if source not in terminal:
            result[source].add(target)
    return {source: sorted(targets) for source, targets in result.items()}


def bfs(adj, source):
    distances, parents = {source: 0}, {}
    queue = collections.deque([source])
    while queue:
        source = queue.popleft()
        for target in adj.get(source, ()):
            if target not in distances:
                distances[target] = distances[source] + 1
                parents[target] = source
                queue.append(target)
    return distances, parents


def reconstruct(source, target, parents):
    result = [target]
    while result[-1] != source:
        result.append(parents[result[-1]])
    return result[::-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="New or empty directory; occupied directories are refused")
    parser.add_argument("--pkn", type=Path,
                        default=ROOT / "results/causal_modelling/cosmos_mouse_pkn",
                        help="Directory containing the frozen PKN and mapping resources")
    parser.add_argument("--inputs", type=Path,
                        default=ROOT / "miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs",
                        help="Directory containing the frozen Factor 1 input CSVs")
    args = parser.parse_args()
    output, pkn, inputs = args.output.resolve(), args.pkn.resolve(), args.inputs.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error(f"Output directory is occupied: {output}")

    pkn_dependencies = [
        "pkn_lipid_rules.csv", "pkn_lipid_curated.csv", "pkn_lipid_curated_ph73.csv",
        "regulation_edges.csv", "mapping_features.csv", "mapping_chebi_nodes.csv",
        "chemical_node_aliases.csv", "regulation_chebi_aliases.csv",
        "regulation_chebi_acid_base_relations.csv", "regulation_chebi_source_crosswalk.csv",
        "rhea_lipid_identity_classification.tsv", "nodes.jsonl",
    ]
    dependencies = [pkn / name for name in pkn_dependencies] + [
        inputs / "factor1_lipid_node_mapping.csv",
        inputs / "factor1_hepatocyte_tf_activities_mapped.csv",
        Path(__file__).resolve(),
    ]
    for path in dependencies:
        if not path.is_file():
            parser.error(f"Missing frozen dependency: {path}")

    original = read_edges(pkn / "pkn_lipid_rules.csv")
    curated = read_edges(pkn / "pkn_lipid_curated.csv")
    normalized = read_edges(pkn / "pkn_lipid_curated_ph73.csv")
    mapped = read_csv(inputs / "factor1_lipid_node_mapping.csv")
    accepted = {r["node_id"]: r["feature"] for r in mapped if r["included"] == "TRUE"}
    assert len(accepted) == sum(r["included"] == "TRUE" for r in mapped)
    assert all(node.startswith("Metab__") for node in accepted)
    tf_rows = read_csv(inputs / "factor1_hepatocyte_tf_activities_mapped.csv")
    tfs = {r["node_id"]: float(r["oriented_activity"]) for r in tf_rows if r["in_pkn"] == "TRUE"}
    tf_sets = {"all_mapped_tfs": set(tfs),
               "abs_activity_gt_2": {node for node, value in tfs.items() if abs(value) > 2}}
    features = read_csv(pkn / "mapping_features.csv")
    candidates = read_csv(pkn / "mapping_chebi_nodes.csv")
    regulations = read_csv(pkn / "regulation_edges.csv")
    crosswalk = {r["source"]: r for r in read_csv(pkn / "regulation_chebi_source_crosswalk.csv")}
    alias_rows = read_csv(pkn / "chemical_node_aliases.csv")
    node_aliases = {r["original_node"]: r["normalized_node"] for r in alias_rows}
    alias_basis = {r["original_node"]: r["basis"] for r in alias_rows}
    lipid_ids = {r["id"] for r in read_csv(pkn / "rhea_lipid_identity_classification.tsv", "\t")}
    lipid_ids |= {key for key, r in crosswalk.items() if r["lipid_is_a"] == "True"}

    def lipid(node):
        return node.startswith("Metab__MeasuredLipid:") or chemical(node) in lipid_ids

    def precursor_filter(edges):
        return [e for e in edges if
                (protein(e[0]) and protein(e[2])) or
                (lipid(e[0]) and protein(e[2])) or
                (lipid(e[0]) and reaction(e[2])) or
                (reaction(e[0]) and lipid(e[2]))]

    direct_regulation = [e for e in original if chemical(e[0]) and protein(e[2])]
    original_reg_keys = set(direct_regulation)
    current_reg_sources = {chemical(e[0]) for e in direct_regulation}
    strict = [e for e in original if (protein(e[0]) and protein(e[2])) or
              (chemical(e[0]) and protein(e[2]))]
    graphs = {
        "direct_regulatory": strict,
        "unrestricted_diagnostic": original,
        "currency_removed_diagnostic": [e for e in original if
            chemical(e[0]) not in CURRENCY and chemical(e[2]) not in CURRENCY],
        "lipid_precursor": precursor_filter(original),
        "ph73_lipid_precursor_sensitivity": precursor_filter(normalized),
    }
    assert not any(protein(s) and reaction(t) for s, _, t in graphs["lipid_precursor"])
    assert all(not chemical(n) or lipid(n) for e in graphs["lipid_precursor"] for n in (e[0], e[2]))
    equivalence = ChemicalEquivalence()
    for row in read_csv(pkn / "regulation_chebi_aliases.csv"):
        equivalence.join(row["alias"], row["primary_chebi"])
    for row in read_csv(pkn / "regulation_chebi_acid_base_relations.csv"):
        equivalence.join(row["acid"], row["base"])
    reg_components = collections.defaultdict(list)
    for row in regulations:
        reg_components[equivalence.find(row["source"])].append(row)
    candidate_opportunities = {}
    for row in candidates:
        # Read source identity once; pH7.3 duplicates cannot inflate opportunities.
        if row["mapping_relation"] != "SwissLipids identifier":
            continue
        mapped_id = row["original_chebi_id"]
        for reg in reg_components.get(equivalence.find(mapped_id), ()):
            source_id, tier = reg["source"], reg["evidence_tier"]
            key = (row["feature_id"], row["match_type"], mapped_id, source_id, tier)
            result = candidate_opportunities.setdefault(key, {
                "feature": row["feature_id"], "match_type": row["match_type"],
                "mapped_chemical": mapped_id, "regulatory_chemical": source_id,
                "identity_relation": "same_identifier" if mapped_id == source_id else
                    "explicit_ChEBI_alias_or_conjugate_acid_base",
                "evidence_tier": tier, "resource": reg["resource"],
                "regulatory_name": reg["source_name"], "targets": set(), "signed_edges": set(),
                "references": set(), "in_current_pkn": False,
            })
            result["targets"].add(reg["target_gene"] or reg["target"])
            result["signed_edges"].add((source_id, reg["interaction"], reg["target"]))
            result["references"].update(filter(None, reg["references"].split(";")))
            current_key = ("Metab__" + source_id + "_unknown", int(reg["interaction"]),
                           "Gene__" + reg["target_gene"])
            result["in_current_pkn"] |= current_key in original_reg_keys
    opportunities = []
    for key, row in sorted(candidate_opportunities.items()):
        row = dict(row)
        row["target_count"] = len(row["targets"])
        row["distinct_signed_regulatory_edges"] = len(row.pop("signed_edges"))
        row["targets"] = "|".join(sorted(row["targets"]))
        row["references"] = ";".join(sorted(row["references"]))
        opportunities.append(row)

    by_feature = collections.defaultdict(list)
    for row in opportunities:
        by_feature[row["feature"]].append(row)
    frozen_by_feature = {r["feature"]: r for r in mapped}
    outgoing = collections.Counter(s for s, _, _ in original)
    outgoing_reg = collections.Counter(s for s, _, _ in direct_regulation)
    mapping_coverage = []
    for feature in features:
        name = feature["feature_id"]
        frozen = frozen_by_feature.get(name, {})
        opps = by_feature[name]
        exact = [r for r in opps if r["match_type"] == "exact_reported_resolution"]
        mapping_coverage.append({
            "feature": name, "best_mapping_status": feature["best_mapping_status"],
            "exact_chebi_ids": feature["exact_chebi_ids"],
            "accepted_frozen_input": frozen.get("included", "FALSE"),
            "frozen_node": frozen.get("node_id", ""),
            "frozen_exclusion_reason": frozen.get("exclusion_reason", ""),
            "oriented_weight": frozen.get("oriented_weight", ""),
            "outgoing_pkn_edges": outgoing[frozen.get("node_id", "")],
            "outgoing_direct_regulatory_edges": outgoing_reg[frozen.get("node_id", "")],
            "exact_current_regulation_same_id": any(r["in_current_pkn"] and
                r["identity_relation"] == "same_identifier" for r in exact),
            "exact_current_regulation_including_acid_base": any(r["in_current_pkn"] for r in exact),
            "nonexact_current_regulation_opportunity": any(r["in_current_pkn"] and
                r["match_type"] != "exact_reported_resolution" for r in opps),
            "exact_cached_regulation_including_acid_base": any(
                r["evidence_tier"].startswith("cached_") for r in exact),
        })

    labels, metadata = {}, {}
    with (pkn / "nodes.jsonl").open() as handle:
        for line in handle:
            row = json.loads(line)
            labels[row["node_id"]] = row.get("label", row["node_id"])
            metadata[row["node_id"]] = row
    topology, reachability, examples, counts = [], [], [], []
    for policy, graph in graphs.items():
        graph_nodes = nodes(graph)
        starts = {node_aliases.get(s, s): f for s, f in accepted.items()} if policy.startswith("ph73") else accepted
        counts.append({"policy": policy, "edges": len(graph), "nodes": len(graph_nodes),
                       "accepted_inputs_present": len(set(starts) & graph_nodes),
                       "protein_to_protein": sum(protein(s) and protein(t) for s, _, t in graph),
                       "chemical_to_protein": sum(bool(chemical(s)) and protein(t) for s, _, t in graph),
                       "lipid_to_reaction": sum(lipid(s) and reaction(t) for s, _, t in graph),
                       "reaction_to_lipid": sum(reaction(s) and lipid(t) for s, _, t in graph),
                       "protein_to_reaction": sum(protein(s) and reaction(t) for s, _, t in graph)})
        for tf_set, terminal in tf_sets.items():
            adj = adjacency(graph, terminal)
            pairs = []
            for source, feature in sorted(starts.items()):
                if source not in graph_nodes:
                    continue
                distances, parents = bfs(adj, source)
                endpoints = sorted(terminal & distances.keys(), key=lambda n: (distances[n], n))
                for target in endpoints:
                    row = {"policy": policy, "tf_endpoint_set": tf_set, "feature": feature,
                           "source": source, "target_tf": target, "shortest_path_edges": distances[target]}
                    reachability.append(row)
                    pairs.append(row)
                selected_targets = []
                reason = ""
                if policy == "currency_removed_diagnostic" and feature.startswith("LPC("):
                    selected_targets, reason = endpoints[:1], "cofactor_shortcut_diagnostic_not_signalling_evidence"
                elif policy == "lipid_precursor":
                    selected_targets = sorted(n for n in distances if protein(n) and n not in terminal)[:2]
                    reason = "protein_dead_end_no_TF_path" if not endpoints else "precursor_path"
                elif policy == "ph73_lipid_precursor_sensitivity":
                    selected_targets, reason = endpoints[:5], "acid_base_normalized_precursor_sensitivity"
                elif policy == "unrestricted_diagnostic" and feature in ("PC(16:0_18:2)", "LPC(16:0/0:0)"):
                    selected_targets, reason = endpoints[:1], "unrestricted_shortcut_diagnostic"
                for target in selected_targets:
                    path = reconstruct(source, target, parents)
                    examples.append({"policy": policy, "tf_endpoint_set": tf_set,
                        "feature": feature, "source": source, "target": target,
                        "shortest_path_edges": len(path) - 1, "interpretation": reason,
                        "node_path": "|".join(path),
                        "label_path": " -> ".join(labels.get(n, n) for n in path),
                        "reaction_ids": "|".join(metadata.get(n, {}).get("reaction_id", "")
                            for n in path if reaction(n))})
            for depth in DEPTHS:
                selected = [r for r in pairs if depth is None or r["shortest_path_edges"] <= depth]
                topology.append({"policy": policy, "tf_endpoint_set": tf_set,
                    "n_tf_observations": len(terminal),
                    "maximum_shortest_path_edges": "unbounded" if depth is None else depth,
                    "reachable_measured_lipids": len({r["source"] for r in selected}),
                    "reachable_tf_endpoints": len({r["target_tf"] for r in selected}),
                    "reachable_lipid_tf_pairs": len(selected),
                    "distance_semantics": "directed shortest path; TF observations terminal; not iterative MOON pruning or scoring"})

    output.mkdir(parents=True, exist_ok=True)
    edge_fields = ["source", "interaction", "target"]
    for name, policy in [("pkn_lipid_precursor.csv", "lipid_precursor"),
                         ("pkn_ph73_lipid_precursor.csv", "ph73_lipid_precursor_sensitivity")]:
        write_csv(output / name, [dict(zip(edge_fields, edge)) for edge in graphs[policy]], edge_fields)
    ph_keys = set(graphs["ph73_lipid_precursor_sensitivity"])
    normalized_provenance = []
    for source, sign, target in curated:
        ns, nt = node_aliases.get(source, source), node_aliases.get(target, target)
        if (ns, sign, nt) in ph_keys:
            normalized_provenance.append({"original_source": source, "interaction": sign,
                "original_target": target, "normalized_source": ns, "normalized_target": nt,
                "source_identity_relation": alias_basis.get(source, "identity"),
                "target_identity_relation": alias_basis.get(target, "identity")})
    assert {(r["normalized_source"], r["interaction"], r["normalized_target"])
            for r in normalized_provenance} == ph_keys
    write_csv(output / "normalized_edge_provenance.csv", normalized_provenance,
              ["original_source", "interaction", "original_target", "normalized_source",
               "normalized_target", "source_identity_relation", "target_identity_relation"])
    exact_edges = [dict(feature=accepted[s], source=s, interaction=sign, target=t)
                   for s, sign, t in direct_regulation if s in accepted]
    write_csv(output / "exact_regulatory_edges.csv", exact_edges,
              ["feature", "source", "interaction", "target"])
    regulatory_rows = []
    for s, sign, t in direct_regulation:
        evidence = [r for r in regulations if r["source"] == chemical(s) and
                    "Gene__" + r["target_gene"] == t and int(r["interaction"]) == sign and
                    r["evidence_tier"] == "curated_mouse_current"]
        assert evidence, f"No curated evidence for {s}, {sign}, {t}"
        regulatory_rows.append({"source": s, "interaction": sign, "target": t,
            "source_name": crosswalk.get(chemical(s), {}).get("name", ""),
            "lipid_ontology_member": lipid(s),
            "evidence_tier": "curated_mouse_current", "resource": "GuideToPharmacology",
            "source_record_ids": "|".join(sorted({r["source_record_id"] for r in evidence})),
            "references": ";".join(sorted({v for r in evidence for v in r["references"].split(";") if v}))})
    write_csv(output / "curated_regulatory_edges.csv", regulatory_rows, list(regulatory_rows[0]))
    write_csv(output / "measured_mapping_coverage.csv", mapping_coverage, list(mapping_coverage[0]))
    write_csv(output / "any_mapping_opportunities.csv", opportunities,
              ["feature", "match_type", "mapped_chemical", "regulatory_chemical", "identity_relation",
               "evidence_tier", "resource", "regulatory_name", "target_count", "targets",
               "distinct_signed_regulatory_edges", "references", "in_current_pkn"])
    write_csv(output / "topology_summaries.csv", topology, list(topology[0]))
    write_csv(output / "endpoint_reachability.csv", reachability,
              ["policy", "tf_endpoint_set", "feature", "source", "target_tf", "shortest_path_edges"])
    write_csv(output / "example_paths.csv", examples,
              ["policy", "tf_endpoint_set", "feature", "source", "target", "shortest_path_edges",
               "interpretation", "node_path", "label_path", "reaction_ids"])
    write_csv(output / "policy_graph_counts.csv", counts, list(counts[0]))
    manifest = []
    for path in dependencies:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        manifest.append({"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()})
    write_csv(output / "dependency_hashes.csv", manifest, ["path", "bytes", "sha256"])
    opportunity_summary = []
    groups = collections.defaultdict(set)
    for row in opportunities:
        groups[(row["evidence_tier"], row["match_type"], row["identity_relation"])].add(row["feature"])
    for (tier, match_type, relation), values in sorted(groups.items()):
        opportunity_summary.append({"evidence_tier": tier, "match_type": match_type,
                                    "identity_relation": relation, "measured_features": len(values)})
    summary = {
        "analysis": "factor1_lipid_to_tf_coverage_audit", "scientific_refit_performed": False,
        "input_mapping_changed": False, "source_pkn_changed": False,
        "measured_features": len(features), "accepted_frozen_input_nodes": len(accepted),
        "accepted_exact_chebi_nodes": sum("CHEBI:" in s for s in accepted),
        "accepted_measured_species_rule_nodes": sum("MeasuredLipid:" in s for s in accepted),
        "all_mapped_tf_observations": len(tfs), "strong_tf_observations": len(tf_sets["abs_activity_gt_2"]),
        "current_regulatory_edges": len(direct_regulation),
        "current_regulatory_chemical_sources": len(current_reg_sources),
        "current_lipid_regulatory_edges": sum(lipid(s) for s, _, _ in direct_regulation),
        "current_lipid_regulatory_sources": len({s for s, _, _ in direct_regulation if lipid(s)}),
        "accepted_inputs_with_direct_regulatory_edges": len({r["source"] for r in exact_edges}),
        "all_measured_exact_current_regulation_same_id": sum(r["exact_current_regulation_same_id"] for r in mapping_coverage),
        "all_measured_exact_current_regulation_including_acid_base": sum(r["exact_current_regulation_including_acid_base"] for r in mapping_coverage),
        "graph_policies": counts, "topology_summaries": topology,
        "mapping_opportunity_counts": opportunity_summary,
        "currency_diagnostic_exclusions": CURRENCY,
        "distance_semantics": "Shortest directed endpoint distance after making the indicated TF observations terminal. Not MOON scoring, significance, or iterative pruning.",
        "precursor_policy": "Only lipid-to-reaction-to-lipid traversal, curated lipid-to-protein regulation and protein-to-protein edges. No protein-to-reaction edges, no nonlipid chemical nodes, no class mapping. Reactions are pairwise traversal hypotheses, not measured flux.",
        "ph73_policy": "Existing curated pH7.3 graph; explicit frozen acid/base alias map. Original six exact ChEBI inputs are used where present; measured-species rules are not appended. This is a sensitivity, not a replacement of the original PKN.",
        "interpretation_caveats": [
            "No direct signaling-ligand assignment follows from substrate participation.",
            "LPC16 -> PC(16:0/22:6) -> DHA is an acyl-remodeling cycle requiring unmeasured DHA-CoA and regenerating LPC. It does not establish net DHA synthesis or a causal effect of LPC abundance.",
            "The 49-strong-TF endpoint policy permits weak measured TFs such as Rxra to act as intermediates; the all-TF policy makes them observations and stops there.",
            "Broad-class and more-specific mapping candidates remain opportunities, never causal edges or automatic input assignments.",
            "Cached STITCH/BRENDA evidence is reported separately and is never incorporated into exported graphs.",
            "Unrestricted and currency-only graphs are diagnostics; cofactor shortcuts can create misleading connectivity.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Factor 1 lipid-to-TF coverage audit\n\n"
        "This is a read-only audit of frozen resource identities and topology. No MOON fit is run. "
        "The original resources and accepted measurement mapping are unchanged.\n\n"
        "`topology_summaries.csv` counts endpoint reachability by shortest directed path with the indicated "
        "TF observations made terminal. These distances differ from iterative forward/backward MOON pruning; "
        "they are not activity scores or significance results. `endpoint_reachability.csv` preserves every pair.\n\n"
        "`pkn_lipid_precursor.csv` and `pkn_ph73_lipid_precursor.csv` contain the filtered original and existing "
        "acid/base-normalized curated graph, respectively. The latter is a sensitivity, includes no measured-species "
        "reaction additions, and uses `normalized_edge_provenance.csv` to recover original edge IDs. "
        "Both omit nonlipid metabolites and protein-to-reaction catalysis so traversal cannot loop from signaling "
        "back into metabolism. Frozen classification admits explicit lipid descendants and their pH7.3 forms.\n\n"
        "`exact_regulatory_edges.csv` contains accepted measured nodes with direct curated regulatory edges; "
        "an empty table is a valid result. `curated_regulatory_edges.csv` gives the complete current regulatory layer. "
        "`any_mapping_opportunities.csv` reports exact, compatible, class and alias candidates separately, including "
        "cached evidence tiers that are not used in the exported graphs.\n\n"
        "The normalized LPC16-to-DHA route is acyl-remodeling participation: LPC plus unmeasured DHA-CoA makes PC, "
        "and phospholipase hydrolysis releases DHA while regenerating LPC. It does not prove net DHA production "
        "or a direct signaling effect of measured LPC. `example_paths.csv` also records precursor dead ends and "
        "unsafe cofactor shortcuts for inspection. Full reaction stoichiometry and catalyst provenance remain in "
        "the frozen PKN `rhea_reactions.jsonl` and `nodes.jsonl`.\n\n"
        "`dependency_hashes.csv` includes the script and every consumed dependency.\n"
    )
    print(json.dumps({"output": str(output), "accepted_input_nodes": len(accepted),
                      "direct_regulatory_input_nodes": len({r["source"] for r in exact_edges}),
                      "topology_configurations": len(topology), "endpoint_pairs": len(reachability),
                      "source_resources_modified": False}, indent=2))


if __name__ == "__main__":
    main()
