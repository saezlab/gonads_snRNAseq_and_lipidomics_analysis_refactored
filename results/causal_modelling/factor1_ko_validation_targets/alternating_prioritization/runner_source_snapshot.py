#!/usr/bin/env python3
"""Export reviewed KO priorities, source-corrected hypotheses and exact paths.

Never refit MOON or mutate frozen graphs. Ranking is reviewed scientific
judgement. Conditional local directions are not simulated interventions.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import copy
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_graph_hashes(review):
    observed = {path: sha(ROOT / path) for path in review["frozen_graph_hashes"]}
    assert observed == review["frozen_graph_hashes"], "Frozen graph hash mismatch"
    return observed


def graph(directory, review, mapping):
    nodefile = directory / "network_labeled_ATT.csv"
    edgefile = directory / "network_labeled_SIF.csv"
    attrs = {r["node_id"]: r for r in read_csv(nodefile)}
    edges = read_csv(edgefile)
    bad_edges = {(r["source"], r["target"]) for r in review["excluded_gene_edges"]}
    bad_genes = {"Gene__" + r["gene"] for r in review["excluded_catalyst_genes"]}
    bad_nodes = set(bad_genes)
    for n, r in attrs.items():
        if set(r["catalysts"].split("|")) & bad_genes:
            bad_nodes.add(n)
    adjacency = defaultdict(list)
    kept, removed = [], []
    for e in edges:
        s, t = e["source"], e["target"]
        sign = int(e["interaction"])
        assert sign in {-1, 1}
        # These checks apply only to the original scored, frozen graph.
        assert int(attrs[s]["level"]) == int(attrs[t]["level"]) + 1
        assert (float(attrs[s]["score"]) * float(attrs[t]["score"]) > 0) == (sign == 1)
        if (s, t) in bad_edges or s in bad_nodes or t in bad_nodes:
            removed.append(e)
        else:
            kept.append(e)
            adjacency[s].append((t, sign))

    def paths(node, nodes=(), signs=()):
        current = nodes + (node,)
        if attrs[node]["type"] == "level0":
            yield current, signs
        else:
            for child, sign in sorted(adjacency[node]):
                assert child not in current
                yield from paths(child, current, signs + (sign,))

    signatures = {}
    for gene in sorted(n for n in attrs if n.startswith("Gene__")):
        endpoints = defaultdict(list)
        for nodes, signs in paths(gene):
            net = 1
            for sign in signs:
                net *= sign
            feature = mapping[nodes[-1]]
            minimum = min(abs(float(attrs[n]["score"])) for n in nodes[:-1])
            endpoints[feature].append(dict(
                nodes=list(nodes), signs=list(signs), path_sign=net,
                conditional_KO_direction="down" if net == 1 else "up",
                levels=[int(attrs[n]["level"]) for n in nodes],
                path_kind="frozen_screened_path", source_hypothesis_id=None,
                source_corrected_prefix=None,
                minimum_inferred_absolute_score=minimum,
                minimum_full_path_inferred_absolute_score=minimum,
                minimum_frozen_portion_inferred_absolute_score=minimum,
                frozen_portion_nodes=list(nodes), frozen_portion_signs=list(signs),
                source_corrected_by_MOON=False))
        for rows in endpoints.values():
            assert len({r["path_sign"] for r in rows}) == 1
        signatures[gene.removeprefix("Gene__")] = dict(endpoints)
    return dict(attrs=attrs, retained_edges=kept, excluded_edges=removed,
                signatures=signatures, sources=[nodefile, edgefile])


def corrected_signatures(network, hypotheses):
    """Prefix proposed source-specific TF regulation; never insert scored nodes."""
    corrected = {}
    for h in hypotheses:
        assert h["path_kind"] == "source_corrected_extension"
        assert h["scored_by_MOON"] is False
        assert h["interaction"] in {-1, 1}
        source, mediator = h["source_gene"], h["mediator_gene"]
        assert source != mediator and mediator in network["signatures"]
        endpoints = {}
        for feature, frozen_paths in network["signatures"][mediator].items():
            rows = []
            for original in frozen_paths:
                assert "Gene__" + source not in original["nodes"]
                p = copy.deepcopy(original)
                net = h["interaction"] * original["path_sign"]
                p.update(nodes=["Gene__" + source] + original["nodes"],
                         signs=[h["interaction"]] + original["signs"],
                         levels=[None] + original["levels"],
                         path_sign=net,
                         conditional_KO_direction="down" if net == 1 else "up",
                         path_kind="source_corrected_extension",
                         source_hypothesis_id=h["hypothesis_id"],
                         source_corrected_prefix=dict(
                             source="Gene__" + source, target="Gene__" + mediator,
                             interaction=h["interaction"], source_ids=h["source_ids"],
                             context=h["context"], source_level=None,
                             source_inferred_score=None, scored_by_MOON=False),
                         minimum_inferred_absolute_score=None,
                         minimum_full_path_inferred_absolute_score=None,
                         minimum_frozen_portion_inferred_absolute_score=
                             original["minimum_frozen_portion_inferred_absolute_score"])
                assert p["levels"][0] is None
                assert p["minimum_full_path_inferred_absolute_score"] is None
                assert p["frozen_portion_nodes"] == original["nodes"]
                rows.append(p)
            endpoints[feature] = rows
        assert endpoints, (source, mediator, "No frozen downstream support")
        corrected[source] = endpoints
    return corrected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path,
                        default=ROOT / "support/cosmos_mouse/factor1_ko_validation_review.json")
    parser.add_argument("--expression", type=Path,
                        default=ROOT / "results/causal_modelling/factor1_ko_validation_targets/expression_audit/gene_expression.csv")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/causal_modelling/factor1_ko_validation_targets/alternating_prioritization")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    review_bytes = args.review.read_bytes()
    review = json.loads(review_bytes)
    assert [r["rank"] for r in review["sets"]] == list(range(1, 11))
    assert len({tuple(sorted(r["targets"])) for r in review["sets"]}) == 10
    assert review["strict_alternation"] is True
    assert review["combination_policy"].startswith("parallel_compensation_only")
    before_hashes = frozen_graph_hashes(review)
    mapping_file = ROOT / "results/causal_modelling/factor1_tf_to_lipid_expanded_species/inputs/input_mapping.csv"
    mapping_rows = read_csv(mapping_file)
    mapping = {r["node_id"]: r["feature_id"] for r in mapping_rows}
    measured_weights = {r["feature_id"]: float(r["weight"]) for r in mapping_rows}
    assert len(mapping) == len(mapping_rows)
    networks = {name: graph(ROOT / path, review, mapping)
                for name, path in review["networks"].items()}
    expression = {r["gene"]: r for r in read_csv(args.expression)}
    sources_file = ROOT / review["sources_file"]
    literature_bytes = sources_file.read_bytes()
    literature = json.loads(literature_bytes)
    sources = {r["id"]: r for r in literature["sources"]}
    tf_file = ROOT / review["frozen_tf_activity_file"]
    tf_rows = {r["source"]: r for r in read_csv(tf_file)}
    tf_targets = set(review["tf_targets"])
    for gene, evidence in review["tf_evidence_qualifiers"].items():
        assert gene in tf_targets and gene in tf_rows
        observed = float(tf_rows[gene]["oriented_activity"])
        assert abs(observed - evidence["frozen_factor1_activity"]) < 1e-10
        expected_input = (abs(observed) > review["original_upstream_absolute_activity_cutoff"]
                          and tf_rows[gene]["in_pkn"] == "TRUE")
        assert expected_input == evidence["supplied_upstream_input"]
        assert evidence["inherited_score_from_other_TF"] is False
    hypotheses = review["source_corrected_hypotheses"]
    assert len({h["source_gene"] for h in hypotheses}) == len(hypotheses)
    assert len({h["hypothesis_id"] for h in hypotheses}) == len(hypotheses)
    for h in hypotheses:
        assert h["source_gene"] in tf_targets
        assert all(s in sources for s in h["source_ids"])
    corrected = {layer: corrected_signatures(network, hypotheses)
                 for layer, network in networks.items()}

    def signatures(layer, gene):
        if gene in corrected[layer]:
            return corrected[layer][gene]
        return networks[layer]["signatures"].get(gene, {})

    def available(layer, gene):
        return (gene in corrected[layer]
                or "Gene__" + gene in networks[layer]["attrs"])

    ranked, endpoint_rows, path_rows = [], [], []
    for item in review["sets"]:
        targets = item["targets"]
        assert 1 <= len(targets) <= 3 and len(set(targets)) == len(targets)
        assert all(g in expression for g in targets)
        is_tf = bool(set(targets) & tf_targets)
        expected_category = "TF-containing" if item["rank"] % 2 == 0 else "Enzyme-only"
        assert is_tf == (item["rank"] % 2 == 0)
        assert item["category"] == expected_category
        assert item["type"] == {1: "Single", 2: "Pair", 3: "Triplet"}[len(targets)]
        name = "layers4" if all(available("layers4", g) for g in targets) else "layers6"
        assert all(available(name, g) for g in targets)
        features = [feature for group in item["primary_readout_groups"]
                    for feature in review["readout_groups"][group]]
        assert len(features) == len(set(features))
        assert item["conditional_gates"]
        compensation = item["compensation"]
        if len(targets) > 1:
            assert compensation and compensation["rationale"]
            assert compensation["compensation_established"] is False
            assert {compensation["primary_route_target"], compensation["backup_route_target"]} == set(targets)
            assert set(review["readout_groups"][compensation["shared_output_group"]]) == set(features)
        else:
            assert compensation is None
        for feature in features:
            supporting, signs, count, kinds, agreements = [], set(), 0, set(), set()
            for gene in targets:
                rows = signatures(name, gene).get(feature, [])
                assert rows, (item["set_id"], gene, feature, "Missing shared endpoint")
                supporting.append(gene)
                for p in rows:
                    serial = bool(set(p["nodes"][1:]) & {"Gene__" + t for t in targets})
                    assert not serial, (item["set_id"], gene, feature, "Serial co-target is not compensation")
                    signs.add(p["conditional_KO_direction"])
                    kinds.add(p["path_kind"])
                    count += 1
                    tf_evidence = review["tf_evidence_qualifiers"].get(gene)
                    agreement = None
                    if tf_evidence:
                        agreement = tf_evidence["frozen_factor1_activity"] * p["path_sign"] * measured_weights[feature] > 0
                        agreements.add(agreement)
                    path_rows.append(dict(
                        set_id=item["set_id"], rank=item["rank"], source_gene=gene,
                        readout=feature, network=name,
                        direction=p["conditional_KO_direction"], path_sign=p["path_sign"],
                        node_ids="|".join(p["nodes"]), levels=json.dumps(p["levels"]),
                        signs=json.dumps(p["signs"]), another_target_on_path=serial,
                        path_kind=p["path_kind"], source_hypothesis_id=p["source_hypothesis_id"],
                        source_prefix_scored_by_MOON=False,
                        minimum_full_path_inferred_absolute_score=p["minimum_full_path_inferred_absolute_score"],
                        minimum_frozen_portion_inferred_absolute_score=p["minimum_frozen_portion_inferred_absolute_score"],
                        frozen_portion_node_ids="|".join(p["frozen_portion_nodes"]),
                        factor1_TF_measurement_sign_agrees=agreement,
                        TF_evidence_qualifier=tf_evidence["evidence_status"] if tf_evidence else "",
                        source_ids="|".join(p["source_corrected_prefix"]["source_ids"])
                                   if p["source_corrected_prefix"] else "frozen_graph_edge_ledger"))
            assert set(supporting) == set(targets)
            assert signs == {item["local_direction"]}, (item["set_id"], feature, signs)
            endpoint_rows.append(dict(
                set_id=item["set_id"], rank=item["rank"], targets="|".join(targets),
                category=item["category"], readout=feature,
                local_direction=item["local_direction"], supporting_targets="|".join(supporting),
                supporting_path_count=count, network=name, path_kinds="|".join(sorted(kinds)),
                factor1_TF_measurement_sign_agrees=next(iter(agreements)) if len(agreements) == 1 else None,
                caveat=item["main_caveat"], prediction_status=item["prediction_status"],
                conditional_gates="|".join(item["conditional_gates"]),
                compensation_rationale=compensation["rationale"] if compensation else ""))
        assert all(s in sources for s in item["source_ids"])
        ranked.append(dict(
            rank=item["rank"], set_id=item["set_id"], targets=" + ".join(targets),
            type=item["type"], category=item["category"], module=item["module"],
            primary_readouts="|".join(features), local_direction=item["local_direction"],
            evidence_network=name,
            source_urls="|".join(sources[s]["url"] for s in item["source_ids"]),
            rationale=item["rationale"], caveat=item["main_caveat"],
            conditional_gates="|".join(item["conditional_gates"]),
            primary_route=compensation["primary_route"] if compensation else "",
            backup_route=compensation["backup_route"] if compensation else "",
            shared_output=compensation["shared_output_group"] if compensation else "",
            compensation_rationale=compensation["rationale"] if compensation else "",
            prediction_status=item["prediction_status"]))
    target_sets = {frozenset(r["targets"]) for r in review["sets"]}
    all_genes = set().union(*target_sets)
    missing_singles = {frozenset([g]) for g in all_genes} - target_sets
    after_hashes = frozen_graph_hashes(review)
    assert before_hashes == after_hashes
    summary = dict(
        passed=True, panel_id=review["panel_id"], sets=len(ranked), endpoint_rows=len(endpoint_rows),
        audited_path_rows=len(path_rows),
        path_kind_counts=dict(Counter(p["path_kind"] for p in path_rows)),
        target_set_type_counts=dict(Counter(r["type"] for r in ranked)),
        category_counts=dict(Counter(r["category"] for r in ranked)),
        strict_alternation_passed=True, all_pairs_share_endpoints=True,
        serial_co_target_paths_rejected=True, ranked_targets=sorted(all_genes),
        additional_constituent_singles=sorted(next(iter(s)) for s in missing_singles),
        distinct_target_conditions_with_single_controls=len(target_sets | missing_singles),
        corrected_hypotheses=len(hypotheses),
        corrected_source_levels_and_full_path_scores_null=True,
        sign_rule="Functional loss(-1) times directed path signs; not baseline MOON score sign or Factor1 alignment",
        biological_screen="Unsupported identities and ENPP2 omitted from frozen path copy; separately declared source-corrected prefixes are unscored hypotheses",
        scoring_performed=False, original_graphs_modified=False, frozen_graph_hashes_verified=True,
        ranking_policy=review["ranking_policy"],
        experimental_system_assumption=review["experimental_system_assumption"],
        finite_lipid_pool_or_flux_simulation=False,
        quantitative_additivity_or_compensation_predicted=False)
    if args.validate_only:
        print(json.dumps(summary, indent=2))
        return
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit("Nonempty output; choose a new --output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, rows in [("ranked_target_sets.csv", ranked),
                           ("expected_lipid_directions.csv", endpoint_rows),
                           ("audited_paths.csv", path_rows),
                           ("path_provenance.csv", path_rows)]:
        write_csv(args.output / filename, rows, list(rows[0]))
    gene_rows = []
    for layer, network in networks.items():
        combined = dict(network["signatures"])
        combined.update(corrected[layer])
        for gene, endpoints in combined.items():
            for feature, paths in endpoints.items():
                gene_rows.append(dict(
                    network=layer, gene=gene, readout=feature,
                    local_KO_direction=paths[0]["conditional_KO_direction"],
                    supporting_paths=len(paths), path_kind=paths[0]["path_kind"],
                    source_hypothesis_id=paths[0]["source_hypothesis_id"],
                    minimum_full_path_inferred_absolute_score=paths[0]["minimum_full_path_inferred_absolute_score"]))
    write_csv(args.output / "all_identity_screened_gene_predictions.csv", gene_rows, list(gene_rows[0]))
    exclusions = [dict(network=layer, source=e["source"], target=e["target"], interaction=e["interaction"])
                  for layer, network in networks.items() for e in network["excluded_edges"]]
    write_csv(args.output / "excluded_path_edges.csv", exclusions, list(exclusions[0]))
    hypothesis_rows = []
    for h in hypotheses:
        hypothesis_rows.append(dict(
            hypothesis_id=h["hypothesis_id"], source_gene=h["source_gene"],
            mediator_gene=h["mediator_gene"], interaction=h["interaction"],
            path_kind=h["path_kind"], scored_by_MOON=False, source_level=None,
            source_inferred_score=None, inherited_score_from_other_TF=False,
            minimum_full_path_inferred_absolute_score=None,
            original_misassigned_source_gene=h["original_misassigned_source_gene"],
            pmid=h["pmid"], source_ids="|".join(h["source_ids"]),
            source_urls="|".join(sources[s]["url"] for s in h["source_ids"]),
            context=h["context"], required_gate=h["required_gate"],
            frozen_factor1_activity=review["tf_evidence_qualifiers"][h["source_gene"]]["frozen_factor1_activity"],
            factor1_alignment=review["tf_evidence_qualifiers"][h["source_gene"]]["factor1_alignment"]))
    write_csv(args.output / "source_corrected_hypothesis_edges.csv", hypothesis_rows, list(hypothesis_rows[0]))
    source_paths = [args.review, args.expression, mapping_file, tf_file, sources_file, Path(__file__).resolve()]
    source_paths.extend(ROOT / path for path in review["frozen_graph_hashes"])
    assert frozen_graph_hashes(review) == before_hashes
    write_csv(args.output / "source_manifest.csv",
              [dict(path=str(p), sha256=sha(p)) for p in dict.fromkeys(source_paths)], ["path", "sha256"])
    (args.output / "review_snapshot.json").write_bytes(review_bytes)
    (args.output / "runner_source_snapshot.py").write_bytes(Path(__file__).read_bytes())
    (args.output / "literature_source_snapshot.json").write_bytes(literature_bytes)
    (args.output / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.output / "source_corrected_hypotheses.json").write_text(json.dumps(hypotheses, indent=2) + "\n")
    (args.output / "frozen_graph_hash_verification.json").write_text(json.dumps(
        dict(before=before_hashes, after=after_hashes, unchanged=True), indent=2) + "\n")
    (args.output / "path_audit.json").write_text(json.dumps(
        dict(frozen_screened={k: n["signatures"] for k, n in networks.items()},
             source_corrected_extensions=corrected), indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
