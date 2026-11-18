#!/usr/bin/env python3
"""Export complete frozen reaction/catalyst evidence for an existing solution.

Read-only with respect to the PKN, selected network and fitted scores. No fitting,
new biological edges, input remapping or literature downloading is performed.
"""
import argparse
import collections
import csv
import gzip
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


def encoded(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Existing selected-network directory")
    parser.add_argument("--output", type=Path, required=True, help="New or empty output directory")
    parser.add_argument("--pkn", type=Path, default=ROOT / "results/causal_modelling/cosmos_mouse_pkn")
    args = parser.parse_args()
    results, out, pkn = args.results.resolve(), args.output.resolve(), args.pkn.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error(f"Output directory is occupied: {out}")
    raw_interactions = ROOT / "support/cosmos_mouse/raw/regulation_gtopdb_interactions.csv.gz"
    dependencies = [results / "network_SIF.csv", results / "network_regulatory_transitions.csv",
                    pkn / "nodes_ph73.jsonl", pkn / "rhea_reactions.jsonl",
                    pkn / "chemical_node_aliases.csv", pkn / "regulation_edges.csv",
                    raw_interactions, Path(__file__).resolve()]
    for path in dependencies:
        if not path.is_file():
            parser.error(f"Missing input: {path}")
    sif = read_csv(results / "network_SIF.csv")
    selected_nodes = {r[k] for r in sif for k in ("source", "target")}
    enzymes = {n for n in selected_nodes if n.startswith("Enzyme__")}
    node_catalog = {r["node_id"]: r for r in (json.loads(line) for line in (pkn / "nodes_ph73.jsonl").open())}
    assert enzymes <= node_catalog.keys()
    reaction_nodes = collections.defaultdict(list)
    for node in sorted(enzymes):
        reaction_nodes[node_catalog[node]["reaction_id"]].append(node)
    reaction_records = {r["reaction_id"]: r for r in (json.loads(line) for line in (pkn / "rhea_reactions.jsonl").open())
                        if r["reaction_id"] in reaction_nodes}
    assert reaction_records.keys() == reaction_nodes.keys(), "Selected reaction lacks a complete Rhea record"
    node_aliases = {r["original_node"]: r["normalized_node"] for r in read_csv(pkn / "chemical_node_aliases.csv")}

    def participant(row):
        original = "Metab__" + row["id"] + "_unknown"
        normalized = node_aliases.get(original, original)
        return {"id": row["id"], "name": row["name"], "stoichiometry": row["stoichiometry"],
                "unknown_compartment_normalized_node": normalized,
                "node_present_in_selected_network": normalized in selected_nodes}

    rows = []
    for reaction_id, ids in sorted(reaction_nodes.items()):
        record = reaction_records[reaction_id]
        selected_catalysts = sorted({n for node in ids for n in node_catalog[node].get("catalysts", [])})
        annotated_genes = sorted({r["gene_symbol"] for r in record.get("provenance", []) if r.get("gene_symbol")})
        substrates = [participant(row) for row in record["substrates"]]
        products = [participant(row) for row in record["products"]]
        rows.append({
            "reaction_id": reaction_id, "master_reaction_id": record.get("master_reaction_id", ""),
            "direction": record.get("direction", ""), "complete_equation": record["equation"],
            "rhea_url": "https://www.rhea-db.org/rhea/" + reaction_id.split(":")[1],
            "organism": record.get("organism", ""),
            "selected_reaction_nodes_json": encoded(ids),
            "selected_node_catalysts_json": encoded(selected_catalysts),
            "all_frozen_mouse_catalyst_genes_json": encoded(annotated_genes),
            "all_frozen_mouse_catalyst_uniprot_json": encoded(record.get("catalysts", [])),
            "substrates_json": encoded(substrates), "products_json": encoded(products),
            "substrates_absent_from_selected_network_json": encoded([r for r in substrates if not r["node_present_in_selected_network"]]),
            "products_absent_from_selected_network_json": encoded([r for r in products if not r["node_present_in_selected_network"]]),
            "selected_incoming_edges_json": encoded([r for r in sif if r["target"] in ids]),
            "selected_outgoing_edges_json": encoded([r for r in sif if r["source"] in ids]),
            "complete_frozen_annotation_provenance_json": encoded(record.get("provenance", [])),
            "interpretation": "Annotated reaction participation; no flux, catalyst-activity or full co-substrate constraint is inferred from this pairwise graph.",
        })
    regulation = {r["edge_id"]: r for r in read_csv(pkn / "regulation_edges.csv")}
    with gzip.open(raw_interactions, "rt", newline="") as handle:
        version = next(handle).strip()
        gtopdb = {i: r for i, r in enumerate(csv.DictReader(handle), start=1)}
    regulatory_rows = []
    for edge in read_csv(results / "network_regulatory_transitions.csv"):
        for edge_id in edge["original_edge_ids"].split("|"):
            r = regulation[edge_id]
            index = int(r["source_record_id"].rsplit(":", 1)[1])
            raw = gtopdb[index]
            assert raw["Target UniProt ID"] == r["target"]
            assert "GtoPdb:" + raw["Ligand ID"] == r["mapping_group"]
            regulatory_rows.append({
                "selected_source": edge["source"], "selected_target": edge["target"],
                "interaction": edge["interaction"], "source_label": edge["source_label"],
                "original_edge_id": edge_id, "original_chemical_id": r["source"],
                "target_uniprot": r["target"], "frozen_record_enumeration": index,
                "persistent_interaction_id": "not_supplied; enumeration is local to frozen export",
                "gtopdb_version": version, "ligand_id": raw["Ligand ID"],
                "target_id": raw["Target ID"], "target_species": raw["Target Species"],
                "curated_action": raw["Action"], "affinity_units": raw["Affinity Units"],
                "affinity_value": raw["Affinity Median"],
                "original_affinity_units": raw["Original Affinity Units"],
                "original_affinity_median_nM": raw["Original Affinity Median nm"],
                "pubmed_ids": raw["PubMed ID"], "assay_description_in_export": raw["Assay Description"],
                "complete_frozen_record_json": encoded(raw),
                "interpretation": "Positive action is inherited from GtoPdb; pKi/Ki are binding affinity, not functional efficacy or hepatocyte exposure.",
            })
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "selected_reaction_evidence.csv", rows, list(rows[0]))
    write_csv(out / "selected_regulatory_evidence.csv", regulatory_rows, list(regulatory_rows[0]))
    hashes = [{"path": str(path), "bytes": path.stat().st_size,
               "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in dependencies]
    write_csv(out / "dependency_hashes.csv", hashes, ["path", "bytes", "sha256"])
    summary = {"selected_network_edges": len(sif), "selected_network_nodes": len(selected_nodes),
               "selected_reaction_nodes": len(enzymes), "distinct_reaction_accessions": len(rows),
               "selected_regulatory_evidence_rows": len(regulatory_rows),
               "scoring_performed": False, "source_PKN_or_fit_modified": False,
               "reaction_equations_source": str(pkn / "rhea_reactions.jsonl"),
               "participant_node_presence_semantics": "One-step frozen alias lookup for the participant's unknown-compartment node; absence/presence is graph bookkeeping, not a concentration or co-substrate measurement."}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "README.md").write_text(
        "# Complete reaction and regulatory evidence for the selected solution\n\n"
        "Generated without fitting or changing the selected network. `selected_reaction_evidence.csv` includes "
        "every selected enzyme-reaction node, grouped by accession, complete frozen Rhea equations and all "
        "annotated mouse catalyst alternatives. JSON fields preserve the original participant stoichiometries, "
        "source annotations and selected incident edges. Graph presence does not establish measured availability "
        "or full multi-substrate reaction feasibility.\n\n"
        "`selected_regulatory_evidence.csv` recovers the original frozen GtoPdb records and their binding affinities. "
        "The record enumeration is not an official persistent interaction accession. Ki is not EC50 or agonist efficacy. "
        "The interpretive primary-source audit is in `docs/cosmos_factor1_lipid_tf_selected_evidence.md`.\n\n"
        "Reproduce with `python3 scripts/moon/ancillary/audit_lipid_tf_selected_reactions.py "
        "--results=results/causal_modelling/factor1_lipid_to_tf/refinement/selected_all_tfs --output=NEW_EMPTY_DIRECTORY`. "
        "Occupied output directories are refused.\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
