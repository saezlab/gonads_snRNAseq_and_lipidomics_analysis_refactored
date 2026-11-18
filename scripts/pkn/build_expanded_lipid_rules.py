#!/usr/bin/env python3
"""Build separate chain-compatible and glycan-conditional species-rule PKNs.

Uses frozen directional mouse Rhea chemistry; never modifies the baseline PKN.
The graph represents positive reaction participation, not flux or concentration
regulation. Candidate chemical identities and unmeasured donors are explicit.
No measurements are assigned and no MOON scoring is performed.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import copy
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

from expanded_lipid_matching import (
    parse_feature, match_headgroup, match_lyso_diacyl, match_sphingolipid,
)

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
CATALOG = ROOT / "support/cosmos_mouse/expanded_species_templates.json"
DEFAULT_OUT = ROOT / "results/causal_modelling/cosmos_mouse_pkn_expanded_species"
INVENTORY = ROOT / "results/causal_modelling/lipid_reaction_target_mapping_audit/feature_reaction_target_inventory.csv"
SEEDS = ROOT / "results/causal_modelling/factor1_tf_to_lipid/upstream_input.csv"
UNIPROT = ROOT / "support/cosmos_mouse/raw/uniprot_mouse.tsv.gz"
TIERS = ("chain_compatible", "glycan_conditional")


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def jsonlines(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_rows(path, data, fields=None):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def measured_identity(feature, recipe):
    condition = None
    if feature["class"] in {"HexCer", "Hex2Cer"}:
        condition = ("lactosyl" if feature["class"] == "Hex2Cer" else
                     "glucose" if recipe["glycan_conditions"].startswith("glucose") else "galactose")
    prefix = "CandidateMeasuredLipid:" if condition else "MeasuredLipid:"
    return prefix + feature["feature_id"] + ("|" + condition if condition else ""), condition


def substitution(side, index, participant, identifier, label, role, **metadata):
    changed = copy.deepcopy(participant)
    changed.update(
        id=identifier, name=label, formula=None, rhea_compound=None,
        template_chemical_id=participant["id"], template_name=participant["name"],
        template_formula=participant["formula"], template_rhea_compound=participant["rhea_compound"],
        chemical_type=["conditional measured-species rule substitution"],
        measurement_resolution="reported chains; unresolved positions/isomers remain hypotheses",
        **metadata,
    )
    return (side, index, changed, dict(
        side=side, index=index, template_chemical_id=participant["id"],
        new_id=identifier, role=role, **metadata,
    ))


def measured_sub(template, side, index, feature, recipe):
    identifier, condition = measured_identity(feature, recipe)
    return substitution(
        side, index, template[side][index], identifier, feature["feature_id"], "measured_lipid",
        measurement_feature_id=feature["feature_id"], candidate_identity_condition=condition,
        reported_chains=feature["chains"],
    )


def acyl_sub(template, recipe, details):
    spec = recipe["acyl_transfer_participant"]
    side, index = spec["side"], spec["index"]
    chain = details["added_chain"]
    n, db = map(int, chain.split(":"))
    carrier = "ACoA" if side == "substrates" else "FA"
    value = substitution(
        side, index, template[side][index], f"HypotheticalAcyl:{carrier}({chain})",
        f"unmeasured {carrier}({chain})", "conditional_acyl_pool",
        conditional_acyl_chain=chain, unmeasured=True,
    )
    # Formula is only a composition annotation, not an exact positional identity.
    value[2]["formula"] = (f"C{n+21}H{2*n+30-2*db}N7O17P3S" if carrier == "ACoA"
                            else f"C{n}H{2*n-1-2*db}O2")
    value[2]["formula_basis"] = "generic Rhea carrier with reported unmodified acyl composition"
    return value


def instantiate(template, recipe, substitutions, details, coupled=None, context=None):
    signature = [template["reaction_id"],
                 [[s, i, p["id"]] for s, i, p, _ in substitutions]]
    record = copy.deepcopy(template)
    record.update(
        reaction_id="EXPANDED_LIPID_RULE:" + hashlib.sha256(json.dumps(signature).encode()).hexdigest()[:20],
        resource="Measured_species_rule_expanded",
        supporting_reaction_id=template["reaction_id"],
        supporting_master_reaction_id=template["master_reaction_id"],
        rule_family=recipe["family"],
        overlay_tier="glycan_conditional" if recipe.get("glycan_conditions") else "chain_compatible",
        template_recipe={k: v for k, v in recipe.items() if k != "participants"},
        match_details=details, coupled_matches=coupled or [],
        context_pathway_feature_ids=context or [],
        source_template_equation=template["equation"], equation=None,
        template_chemically_balanced=template["chemically_balanced"], chemically_balanced=None,
        chemical_balance_scope="Formal template participants/coefficients and matched chains preserved; no independent atom/flux balance test",
        reference_scope="Generic reaction/catalyst references do not establish activity on these measured species",
        substitutions=[s[3] for s in substitutions],
        matched_feature_ids=sorted({s[3]["measurement_feature_id"] for s in substitutions
                                   if s[3]["role"] == "measured_lipid"}),
        assumptions=list(dict.fromkeys(details.get("assumptions", []) + recipe["constraints"] + [
            "All measured-species enzyme specificities are hypotheses.",
            "Positive graph participation does not encode substrate depletion, flux, or joint substrate availability.",
            "Synthetic acyl pools/intermediates have no measured values or asserted exact ChEBI identities.",
        ])),
    )
    for side, index, participant, metadata in substitutions:
        record[side][index] = participant
    return record


def build_reactions(features, templates, recipes):
    by_class = defaultdict(list)
    for feature in features:
        by_class[feature["class"]].append(feature)
    records, exclusions, template_counts = {}, Counter(), Counter()
    for recipe in recipes:
        template = templates[recipe["reaction_id"]]
        assert template["direction"] in {"LR", "RL"} and template["reversible"] is False
        assert not template["is_transport"] and template["catalysts"]
        assert all(isinstance(p["stoichiometry"], (int, float))
                   for p in template["substrates"] + template["products"])
        if recipe["family"] == "Pemt serial methylation":
            continue
        for a in by_class[recipe["source_class"]]:
            for b in by_class[recipe["target_class"]]:
                si, pi = recipe["source_index"], recipe["target_index"]
                family = recipe["family"]
                if family in {"acyl-CoA acylation", "hydrolysis"}:
                    lyso, diacyl = (a, b) if family == "acyl-CoA acylation" else (b, a)
                    details = match_lyso_diacyl(lyso, diacyl)
                    if details:
                        sn = recipe.get("required_retained_sn_position", recipe["retained_chain_sn"])
                        details["positional_options"] = [p for p in details["positional_options"]
                                                         if sn is None or p["retained_sn"] == sn]
                        if not details["positional_options"]:
                            exclusions["template_sn_incompatible"] += 1
                            continue
                    if not details:
                        exclusions["retained_chain_incompatible"] += 1
                        continue
                    subs = [measured_sub(template, "substrates", si, a, recipe),
                            measured_sub(template, "products", pi, b, recipe),
                            acyl_sub(template, recipe, details)]
                    r = instantiate(template, recipe, subs, details)
                    records[r["reaction_id"]] = r
                    template_counts[recipe["reaction_id"]] += 1
                elif family == "diacyl headgroup":
                    details = match_headgroup(a, b)
                    if details:
                        r = instantiate(template, recipe, [
                            measured_sub(template, "substrates", si, a, recipe),
                            measured_sub(template, "products", pi, b, recipe)], details)
                        records[r["reaction_id"]] = r
                        template_counts[recipe["reaction_id"]] += 1
                else:
                    details = match_sphingolipid(a, b)
                    if not details:
                        continue
                    if details["n_acyl_oxygen"] or details["sphingobase_oxygen"] != 2:
                        exclusions["sphingolipid_oxygen_template_incompatible"] += 1
                        continue
                    if recipe["sphingobase_scope"] == "fixed18:1;O2" and details["sphingobase"] != "18:1;O2":
                        exclusions["fixed_sphingobase_template_incompatible"] += 1
                        continue
                    subs = [measured_sub(template, "substrates", si, a, recipe),
                            measured_sub(template, "products", pi, b, recipe)]
                    requirements = recipe.get("coupled_lipid_requirements")
                    if requirements:
                        donor = next(p for p in requirements["participants"] if p["side"] == "substrates")
                        coproduct = next(p for p in requirements["participants"] if p["side"] == "products")
                        for c in by_class[donor["class"]]:
                            for d in by_class[coproduct["class"]]:
                                coupled = match_headgroup(c, d)
                                if coupled:
                                    full = subs + [
                                        measured_sub(template, "substrates", donor["index"], c, recipe),
                                        measured_sub(template, "products", coproduct["index"], d, recipe)]
                                    r = instantiate(template, recipe, full, details, [coupled])
                                    records[r["reaction_id"]] = r
                                    template_counts[recipe["reaction_id"]] += 1
                    else:
                        r = instantiate(template, recipe, subs, details)
                        records[r["reaction_id"]] = r
                        template_counts[recipe["reaction_id"]] += 1
    # Preserve three Pemt reactions, sharing chain-specific unmeasured intermediates.
    serial = [recipe for recipe in recipes if recipe["family"] == "Pemt serial methylation"]
    assert [r["reaction_id"] for r in serial] == ["RHEA:11165", "RHEA:32736", "RHEA:32740"]
    for a in by_class["PE"]:
        for b in by_class["PC"]:
            details = match_headgroup(a, b)
            if not details:
                continue
            chain_key = ("/" if a["sn_order_reported"] or b["sn_order_reported"] else "_").join(
                a["chains"] if a["sn_order_reported"] else b["chains"] if b["sn_order_reported"] else a["chain_multiset"])
            context = [a["feature_id"], b["feature_id"]]
            for step, recipe in enumerate(serial):
                template = templates[recipe["reaction_id"]]
                subs = []
                for side, cls, feature in [("substrates", recipe["source_class"], a if step == 0 else None),
                                           ("products", recipe["target_class"], b if step == 2 else None)]:
                    index = recipe["source_index"] if side == "substrates" else recipe["target_index"]
                    if feature:
                        subs.append(measured_sub(template, side, index, feature, recipe))
                    else:
                        subs.append(substitution(
                            side, index, template[side][index], f"HypotheticalLipid:{cls}({chain_key})",
                            f"unmeasured {cls}({chain_key})", "conditional_intermediate",
                            reported_chains=a["chain_multiset"], unmeasured=True))
                r = instantiate(template, recipe, subs, details, context=context)
                # Shared serial steps carry multiple contextual endpoints; no duplicate reaction evidence.
                if r["reaction_id"] in records:
                    records[r["reaction_id"]]["context_pathway_feature_ids"] = sorted(set(
                        records[r["reaction_id"]]["context_pathway_feature_ids"] + context))
                else:
                    records[r["reaction_id"]] = r
                    template_counts[recipe["reaction_id"]] += 1
    return list(records.values()), exclusions, template_counts


def project(records, baseline_nodes):
    nodes = {n["node_id"]: n for n in baseline_nodes}
    evidence = []
    with gzip.open(UNIPROT, "rt") as handle:
        annotations = list(csv.DictReader(handle, delimiter="\t"))
    accessions = {r["Entry"] for r in annotations}
    gene_map = {r["Entry"]: r["Gene Names (primary)"] for r in annotations
                if re.fullmatch(r"[A-Za-z0-9_.-]+", r["Gene Names (primary)"])}
    def protein(accession):
        base_accession = re.sub(r"-\d+$", "", accession)
        assert accession in accessions or base_accession in accessions, accession
        symbol = gene_map.get(accession, gene_map.get(base_accession))
        if symbol:
            ident = "Gene__" + symbol
            nodes.setdefault(ident, dict(node_id=ident, node_type="gene_product", label=symbol,
                                          gene_symbol=symbol, organism=10090))
        else:
            ident = "UniProt__" + accession
            nodes.setdefault(ident, dict(node_id=ident, node_type="protein_unmapped", label=accession,
                                          uniprot=accession, organism=10090, mouse_accession_verified=True))
        return ident
    def metabolite(p):
        ident = "Metab__" + p["id"] + "_" + p.get("compartment", "unknown")
        nodes.setdefault(ident, dict(
            node_id=ident, node_type="metabolite", label=p["name"], chemical_id=p["id"],
            compartment=p.get("compartment", "unknown"),
            measurement_feature_id=p.get("measurement_feature_id"),
            candidate_identity_condition=p.get("candidate_identity_condition"),
            unmeasured=p.get("unmeasured", False)))
        return ident
    for record in records:
        catalysts = sorted({protein(p) for p in record["catalysts"]})
        for gene in catalysts:
            signature = [record["resource"], record["reaction_id"], False, [gene]]
            instance = "Enzyme__" + hashlib.sha256(json.dumps(signature).encode()).hexdigest()[:16]
            nodes[instance] = dict(
                node_id=instance, node_type="enzyme_reaction_instance", label=record["reaction_id"],
                reaction_id=record["reaction_id"], resource=record["resource"], reverse=False,
                catalysts=[gene], complex_required=False, uncertain_direction=False,
                overlay_tier=record["overlay_tier"])
            metadata = dict(reaction_id=record["reaction_id"], supporting_reaction_id=record["supporting_reaction_id"],
                            resource=record["resource"], overlay_tier=record["overlay_tier"],
                            master_record="expanded_species_reactions.jsonl",
                            source_version=record["source_version"], annotated_direction=record["direction"])
            evidence.append(dict(source=gene, target=instance, interaction=1, mechanism="catalyst_enables_reaction", **metadata))
            for side, mechanism in [("substrates", "substrate_participation"), ("products", "product_formation")]:
                for p in record[side]:
                    chemical = metabolite(p)
                    source, target = (chemical, instance) if side == "substrates" else (instance, chemical)
                    evidence.append(dict(source=source, target=target, interaction=1, mechanism=mechanism,
                                         stoichiometry=p["stoichiometry"], **metadata))
    return nodes, evidence


def distances(edges, seeds):
    adjacency = defaultdict(list)
    for source, sign, target in edges:
        adjacency[source].append(target)
    dist = {n: 0 for n in seeds}
    queue = deque(seeds)
    while queue:
        source = queue.popleft()
        for target in adjacency[source]:
            if target not in dist:
                dist[target] = dist[source] + 1
                queue.append(target)
    return dist


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Use a new/empty output directory: " + str(out))
    catalog = json.loads(CATALOG.read_text())
    templates = {r["reaction_id"]: r for r in jsonlines(BASE / "rhea_reactions.jsonl")}
    assert all(digest(BASE / "rhea_reactions.jsonl") == r["full_record_source_sha256"]
               for r in catalog["templates"])
    mapping = rows(BASE / "mapping_features.csv")
    features = [p for r in mapping if (p := parse_feature(r))]
    records, exclusions, counts = build_reactions(features, templates, catalog["templates"])
    baseline = {(r["source"], int(r["interaction"]), r["target"]) for r in rows(BASE / "pkn_lipid_rules.csv")}
    active_base_nodes = {n for source, sign, target in baseline for n in (source, target)}
    nodes, evidence = project(records, [n for n in jsonlines(BASE / "nodes.jsonl") if n["node_id"] in active_base_nodes])
    out.mkdir(parents=True, exist_ok=True)
    with (out / "expanded_species_reactions.jsonl").open("w") as handle:
        for record in sorted(records, key=lambda r: r["reaction_id"]):
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    with gzip.open(out / "extension_edge_evidence.jsonl.gz", "wt") as handle:
        for row in evidence:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    with (out / "nodes.jsonl").open("w") as handle:
        for node in sorted(nodes.values(), key=lambda n: n["node_id"]):
            handle.write(json.dumps(node, sort_keys=True) + "\n")
    graph_summary, endpoint_rows = {}, []
    inventory = rows(INVENTORY)
    seeds = [r["source"] for r in rows(SEEDS)]
    baseline_products = {r["feature_id"]: r["accepted_node"] for r in inventory
                         if r["accepted_node_is_active_product"] == "True"}
    baseline_distances = distances(baseline, seeds)
    baseline_reachable = {str(depth): sum(baseline_distances.get(node, float("inf")) <= depth
                                        for node in baseline_products.values()) for depth in [4, 6, 8, 12]}
    focus_rows = []
    for tier_index, (tier, filename) in enumerate(zip(TIERS, ["pkn_expanded_chain.csv", "pkn_expanded_glycan_conditional.csv"])):
        allowed = set(TIERS[:tier_index+1])
        edges = baseline | {(r["source"], r["interaction"], r["target"]) for r in evidence if r["overlay_tier"] in allowed}
        pairs = defaultdict(set)
        for source, sign, target in edges:
            pairs[source, target].add(sign)
        assert all(len(signs) == 1 for signs in pairs.values())
        assert all(source != target for source, sign, target in edges)
        write_rows(out / filename, [dict(source=s, interaction=i, target=t) for s, i, t in sorted(edges)],
                   ["source", "interaction", "target"])
        products = defaultdict(set)
        for row in inventory:
            if row["accepted_node_is_active_product"] == "True":
                products[row["feature_id"]].add(row["accepted_node"])
        for record in records:
            if record["overlay_tier"] in allowed:
                for p in record["products"]:
                    if p.get("measurement_feature_id"):
                        products[p["measurement_feature_id"]].add("Metab__" + p["id"] + "_" + p.get("compartment", "unknown"))
        dist = distances(edges, seeds)
        reachable = {feature: min([dist[n] for n in identifiers if n in dist], default=None)
                     for feature, identifiers in products.items()}
        for row in mapping:
            feature = row["feature_id"]
            endpoint_rows.append(dict(
                overlay=tier, feature_id=feature, lipid_class=row["class"],
                in_historical_factor1_64=row["in_factor1"], product_nodes="|".join(sorted(products[feature])),
                product_node_count=len(products[feature]), has_product_candidate=bool(products[feature]),
                minimum_unfiltered_distance_from_frozen_49_tfs=reachable.get(feature),
                glycan_identity_conditional=any("CandidateMeasuredLipid:" in n for n in products[feature])))
        graph_summary[tier] = dict(
            edges=len(edges), nodes=len({n for s, i, t in edges for n in (s, t)}),
            added_edges=len(edges-baseline), candidate_measured_product_features=sum(bool(v) for v in products.values()),
            new_candidate_measured_product_features=sum(bool(v) and
                not any(r["feature_id"] == f and r["accepted_node_is_active_product"] == "True" for r in inventory)
                for f, v in products.items()),
            unfiltered_reachable_product_features={str(depth): sum(d is not None and d <= depth for d in reachable.values())
                                                  for depth in [4, 6, 8, 12]},
            mapping_is_proposed_input_vector=False,
        )
        for tf in ["Gene__Esr1", "Gene__Nfic", "Gene__Dbp"]:
            focus_distances = distances(edges, [tf])
            for depth in [4, 6, 8, 12]:
                focus_rows.append(dict(overlay=tier, tf_node=tf, maximum_directed_distance=depth,
                    candidate_product_features_reachable=sum(
                        any(focus_distances.get(n, float("inf")) <= depth for n in values)
                        for values in products.values())))
    write_rows(out / "measured_product_candidate_inventory.csv", endpoint_rows)
    write_rows(out / "focus_tf_unfiltered_reachability.csv", focus_rows)
    write_rows(out / "template_instantiation_audit.csv", [
        dict(reaction_id=r["reaction_id"], family=r["family"], transformation=r["transformation"],
             instances=counts[r["reaction_id"]], sphingobase_scope=r.get("sphingobase_scope"),
             glycan_conditions=r.get("glycan_conditions")) for r in catalog["templates"]])
    summary = dict(
        baseline_edges=len(baseline), baseline_accepted_active_product_features=len(baseline_products),
        baseline_unfiltered_reachable_product_features=baseline_reachable, frozen_tf_seeds=len(seeds),
        eligible_explicit_chain_features=len(features), expanded_reaction_hypotheses=len(records),
        new_reaction_families=dict(Counter(r["rule_family"] for r in records)),
        tiers=dict(Counter(r["overlay_tier"] for r in records)),
        exclusions=dict(exclusions), graphs=graph_summary,
        topology_scope="Unfiltered directed reachability from frozen49TFseed IDs; ignores activity signs, RNA filtering, cofactors/AND availability and scores; not MOON preparation",
        no_omics_scoring_performed=True, no_measurement_assignment_finalized=True,
    )
    write_json(out / "build_summary.json", summary)
    dependencies = [Path(__file__), Path(__file__).with_name("expanded_lipid_matching.py"), CATALOG, UNIPROT,
                    INVENTORY, SEEDS, *[BASE / f for f in ["rhea_reactions.jsonl", "mapping_features.csv", "pkn_lipid_rules.csv", "nodes.jsonl"]]]
    write_json(out / "build_manifest.json", dict(
        inputs=[dict(path=str(p.relative_to(ROOT)), sha256=digest(p)) for p in dependencies],
        frozen_baseline_unchanged=True, hypothesis_scope="Conditional species chemistry with explicit mapping/isomer assumptions",
        measurement_assignment="No new MOON input vectors created; audit candidate ledger before scoring",
        node_projection="Positive catalyst/substrate participation and product formation; no inferred flux or regulatory sign",
    ))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
