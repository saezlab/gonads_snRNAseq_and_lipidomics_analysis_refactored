#!/usr/bin/env python3
"""Audit the source-direction and graph invariants of a STITCH overlay."""

import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn_stitch_directional"
CID = re.compile(r"^CID[ms]\d+$")
MOUSE_PROTEIN = re.compile(r"^10090\.ENSMUSP\d+$")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonl(path):
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def graph_node(row):
    target = "Gene__" + row["target_gene"] if row["target_gene"] else "UniProt__" + row["target"]
    return "Metab__" + row["source"] + "_unknown", int(row["interaction"]), target


def validate(output):
    evidence_path = output / "stitch_directional_evidence.csv"
    ledger_path = output / "edge_evidence.jsonl.gz"
    graph_path = output / "pkn_exploratory_stitch_directional.csv"
    regulation_path = output / "regulation_edges.csv"
    manifest_path = output / "stitch_directional_manifest.json"
    for path in [evidence_path, ledger_path, graph_path, regulation_path, manifest_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    with evidence_path.open(newline="") as handle:
        evidence = list(csv.DictReader(handle))
    with regulation_path.open(newline="") as handle:
        regulation = list(csv.DictReader(handle))
    with graph_path.open(newline="") as handle:
        graph_rows = list(csv.DictReader(handle))
    ledger = list(jsonl(ledger_path))
    manifest = json.loads(manifest_path.read_text())

    signs_match = all(int(row["interaction"]) == (1 if row["original_action"] == "activation" else -1)
                       for row in evidence)
    checks = {
        "nonempty_directional_evidence": bool(evidence),
        "raw_actor_is_chemical": all(CID.fullmatch(row["raw_item_id_a"]) for row in evidence),
        "raw_receiver_is_mouse_ENSMUSP": all(MOUSE_PROTEIN.fullmatch(row["raw_item_id_b"]) for row in evidence),
        "raw_chemical_is_acting_partner": all(row["raw_a_is_acting"] == "t" for row in evidence),
        "raw_action_is_explicitly_signed": all(row["original_action"] in {"activation", "inhibition"} for row in evidence),
        "sign_derived_from_raw_action": signs_match,
        "action_score_meets_cutoff": all(int(row["raw_action_score"]) >= 700 for row in evidence),
        "detailed_textmining_below_cutoff": all(int(row["raw_link_textmining"]) < 700 for row in evidence),
        "no_cached_STITCH_in_regulation_ledger": not any(row["resource"] == "STITCH" for row in regulation),
        "no_cached_STITCH_in_full_evidence_ledger": not any(row.get("resource") == "STITCH" for row in ledger),
        "directional_evidence_count_matches_full_ledger": sum(
            row.get("resource") == "STITCH_directional_v5" for row in ledger) == len(evidence),
        "no_self_loops": all(row["source"] != row["target"] for row in graph_rows),
        "no_duplicate_source_target_pairs": len({(row["source"], row["target"]) for row in graph_rows}) == len(graph_rows),
    }

    pair_signs = collections.defaultdict(set)
    for row in ledger:
        pair_signs[(row["source"], row["target"])].add(int(row["interaction"]))
    expected_graph = {(source, next(iter(signs)), target)
                      for (source, target), signs in pair_signs.items() if len(signs) == 1}
    observed_graph = {(row["source"], int(row["interaction"]), row["target"]) for row in graph_rows}
    checks["graph_is_exact_projection_of_full_evidence_ledger"] = observed_graph == expected_graph
    checks["opposing_sign_pairs_excluded"] = not any(
        (source, next(iter(signs)), target) in observed_graph
        for (source, target), signs in pair_signs.items() if len(signs) > 1
    )
    base_root = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn"
    checks["base_lipid_rules_hash_unchanged"] = (
        sha256(base_root / "pkn_lipid_rules.csv") == manifest["base_lipid_rules_sha256"]
    )
    checks["base_evidence_hash_unchanged"] = (
        sha256(base_root / "edge_evidence.jsonl.gz") == manifest["base_edge_evidence_sha256"]
    )
    for filename, expected_hash in manifest["base_graph_csv_sha256"].items():
        checks["base_" + filename + "_hash_unchanged"] = sha256(base_root / filename) == expected_hash
    fresh_triples = {graph_node(row) for row in evidence}
    report = {
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "counts": {
            "directional_evidence_rows": len(evidence),
            "unique_directional_projected_signed_triples": len(fresh_triples),
            "directional_projected_signed_triples_surviving_graph_conflict_filter": len(fresh_triples & observed_graph),
            "directional_projected_signed_triples_quarantined": len(fresh_triples - observed_graph),
            "final_exploratory_edges": len(observed_graph),
            "final_conflicting_pairs": sum(len(signs) > 1 for signs in pair_signs.values()),
        },
        "no_omics_scoring_performed": True,
    }
    (output / "stitch_directional_validation.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if not report["all_checks_pass"]:
        failed = [name for name, passed in checks.items() if not passed]
        raise AssertionError("Failed overlay checks: " + ", ".join(failed))
    print(json.dumps(report, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    validate(args.output)


if __name__ == "__main__":
    main()
