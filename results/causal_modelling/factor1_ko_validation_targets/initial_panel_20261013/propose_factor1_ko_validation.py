#!/usr/bin/env python3
"""Export reviewed KO priorities and exact signed paths; never refit MOON.

The ranking is explicit scientific judgement in a reviewed configuration.
This program verifies identities and derives conditional local intervention
directions; it is not a dynamical, flux, synergy or knockout simulation.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
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
        writer.writeheader(); writer.writerows(rows)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

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
        assert int(attrs[s]["level"]) == int(attrs[t]["level"]) + 1
        assert (float(attrs[s]["score"]) * float(attrs[t]["score"]) > 0) == (sign == 1)
        if (s, t) in bad_edges or s in bad_nodes or t in bad_nodes:
            removed.append(e)
        else:
            kept.append(e); adjacency[s].append((t, sign))
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
            for s in signs: net *= s
            feature = mapping[nodes[-1]]
            endpoints[feature].append(dict(nodes=list(nodes), signs=list(signs),
                path_sign=net, conditional_KO_direction="down" if net == 1 else "up",
                levels=[int(attrs[n]["level"]) for n in nodes],
                minimum_inferred_absolute_score=min(abs(float(attrs[n]["score"])) for n in nodes[:-1])))
        for rows in endpoints.values():
            assert len({r["path_sign"] for r in rows}) == 1
        signatures[gene.removeprefix("Gene__")] = dict(endpoints)
    return dict(attrs=attrs, retained_edges=kept, excluded_edges=removed,
                signatures=signatures, sources=[nodefile, edgefile])

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, default=ROOT / "support/cosmos_mouse/factor1_ko_validation_review.json")
    parser.add_argument("--expression", type=Path, default=ROOT / "results/causal_modelling/factor1_ko_validation_targets/expression_audit/gene_expression.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results/causal_modelling/factor1_ko_validation_targets/prioritization")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    review = json.loads(args.review.read_text())
    assert [r["rank"] for r in review["sets"]] == list(range(1, 11))
    assert len({tuple(sorted(r["targets"])) for r in review["sets"]}) == 10
    mapping_file = ROOT / "results/causal_modelling/factor1_tf_to_lipid_expanded_species/inputs/input_mapping.csv"
    mapping = {r["node_id"]: r["feature_id"] for r in read_csv(mapping_file)}
    networks = {name: graph(ROOT / path, review, mapping) for name, path in review["networks"].items()}
    expression_rows = read_csv(args.expression)
    expression = {r["gene"]: r for r in expression_rows}
    sources_file = ROOT / review["sources_file"]
    literature = json.loads(sources_file.read_text())
    sources = {r["id"]: r for r in literature["sources"]}
    ranked, endpoint_rows, path_rows = [], [], []
    for item in review["sets"]:
        targets = item["targets"]
        assert 1 <= len(targets) <= 3 and len(set(targets)) == len(targets)
        assert all(g in expression for g in targets)
        name = "layers4" if all("Gene__" + g in networks["layers4"]["attrs"] for g in targets) else "layers6"
        network = networks[name]
        assert all("Gene__" + g in network["attrs"] for g in targets)
        features = [f for group in item["primary_readout_groups"] for f in review["readout_groups"][group]]
        for feature in features:
            supporting, signs, count = [], set(), 0
            for g in targets:
                rows = network["signatures"][g].get(feature, [])
                if rows: supporting.append(g)
                for p in rows:
                    signs.add(p["conditional_KO_direction"]); count += 1
                    path_rows.append(dict(set_id=item["set_id"], rank=item["rank"], source_gene=g,
                        readout=feature, network=name, direction=p["conditional_KO_direction"],
                        path_sign=p["path_sign"], node_ids="|".join(p["nodes"]),
                        levels=json.dumps(p["levels"]), signs=json.dumps(p["signs"]),
                        another_target_on_path=bool(set(p["nodes"][1:]) & {"Gene__" + t for t in targets}),
                        minimum_inferred_absolute_score=p["minimum_inferred_absolute_score"]))
            assert supporting and signs == {item["local_direction"]}, (item["set_id"], feature, signs)
            endpoint_rows.append(dict(set_id=item["set_id"], rank=item["rank"], targets="|".join(targets),
                readout=feature, local_direction=item["local_direction"], supporting_targets="|".join(supporting),
                supporting_path_count=count, network=name,
                caveat=item["main_caveat"], prediction_status=item["prediction_status"]))
        assert all(s in sources for s in item["source_ids"])
        ranked.append(dict(rank=item["rank"], set_id=item["set_id"], targets=" + ".join(targets),
            type=item["type"], module=item["module"], primary_readouts="|".join(features),
            local_direction=item["local_direction"], evidence_network=name,
            source_urls="|".join(sources[s]["url"] for s in item["source_ids"]),
            rationale=item["rationale"], caveat=item["main_caveat"],
            prediction_status=item["prediction_status"]))
    summary = dict(passed=True, sets=len(ranked), endpoint_rows=len(endpoint_rows),
        ranked_targets=sorted({g for r in review["sets"] for g in r["targets"]}),
        sign_rule="A local loss(-1) times the product of directed path signs; independent of baseline MOON-score sign",
        biological_screen="Remove explicitly unsupported identities/directions and ENPP2 catalyst instances from the audited path copy; conditional/context-dependent edges retained with qualifications",
        scoring_performed=False, original_graphs_modified=False,
        ranking_policy=review["ranking_policy"], experimental_system_assumption=review["experimental_system_assumption"],
        finite_lipid_pool_or_flux_simulation=False, quantitative_additivity_or_epistasis_predicted=False)
    if args.validate_only:
        print(json.dumps(summary, indent=2)); return
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit("Nonempty output; choose a new --output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    for filename, rows in [("ranked_target_sets.csv", ranked), ("expected_lipid_directions.csv", endpoint_rows),
                           ("audited_paths.csv", path_rows)]:
        write_csv(args.output / filename, rows, list(rows[0]))
    gene_rows = []
    for layer, network in networks.items():
        for gene, endpoints in network["signatures"].items():
            for feature, paths in endpoints.items():
                gene_rows.append(dict(network=layer, gene=gene, readout=feature,
                    local_KO_direction=paths[0]["conditional_KO_direction"], supporting_paths=len(paths)))
    write_csv(args.output / "all_identity_screened_gene_predictions.csv", gene_rows, list(gene_rows[0]))
    exclusions = [dict(network=layer, source=e["source"], target=e["target"],
                       interaction=e["interaction"]) for layer, n in networks.items() for e in n["excluded_edges"]]
    write_csv(args.output / "excluded_path_edges.csv", exclusions, list(exclusions[0]))
    source_paths = [args.review, args.expression, mapping_file, sources_file, Path(__file__).resolve()]
    for n in networks.values(): source_paths.extend(n["sources"])
    write_csv(args.output / "source_manifest.csv", [dict(path=str(p), sha256=sha(p)) for p in source_paths], ["path", "sha256"])
    (args.output / "review_snapshot.json").write_text(json.dumps(review, indent=2) + "\n")
    (args.output / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.output / "path_audit.json").write_text(json.dumps({k: n["signatures"] for k, n in networks.items()}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__": main()
