#!/usr/bin/env python3
"""Export provenance and readable node metadata for one relaxed lipid-to-TF solution."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> List[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Iterable[dict], fields: List[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def edge_key(row: dict) -> Tuple[str, str, str]:
    return (row["source"], row["target"], str(row["interaction"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--pkn", required=True, type=Path)
    parser.add_argument("--mapping-evidence", required=True, type=Path)
    args = parser.parse_args()

    results = args.results.resolve()
    pkn = args.pkn.resolve()
    mapping_evidence = args.mapping_evidence.resolve()
    sif_path = results / "network_SIF.csv"
    att_path = results / "network_ATT.csv"
    nodes_path = pkn.parent / "nodes.csv"
    evidence_path = pkn.parent / "edge_evidence.jsonl.gz"
    for path in (sif_path, att_path, nodes_path, evidence_path, mapping_evidence):
        if not path.is_file():
            raise SystemExit(f"Required export input is absent: {path}")

    sif = read_csv(sif_path)
    att = read_csv(att_path)
    keys = {edge_key(row) for row in sif}
    if not keys:
        raise SystemExit("Cannot export an empty relaxed network")

    provenance: Dict[Tuple[str, str, str], List[dict]] = {}
    with gzip.open(evidence_path, "rt") as handle:
        for line in handle:
            row = json.loads(line)
            key = (row["source"], row["target"], str(row["interaction"]))
            if key in keys:
                provenance.setdefault(key, []).append(row)
    missing = sorted(keys.difference(provenance))
    if missing:
        raise SystemExit(f"Selected edges lack assembled PKN provenance: {missing[:5]}")

    evidence_fields = [
        "source", "target", "interaction", "consistency", "resource", "resources", "mechanism",
        "references", "source_accession", "target_accession", "species_evidence", "evidence_tier",
        "source_version", "original_edge_id",
    ]
    evidence_rows: List[dict] = []
    for edge in sif:
        for record in provenance[edge_key(edge)]:
            row = dict(record)
            row["consistency"] = edge.get("consistency", "")
            evidence_rows.append(row)
    write_csv(results / "network_edge_evidence.csv", evidence_rows, evidence_fields)

    mapping_rows = read_csv(mapping_evidence)
    direct_edges = [row for row in sif if row["source"].startswith("Metab__")]
    direct_keys = {edge_key(row) for row in direct_edges}
    direct_mapping_rows = [
        row for row in mapping_rows
        if (row.get("pkn_source_node"), row.get("target_node"), str(row.get("interaction"))) in direct_keys
    ]
    if len(direct_mapping_rows) != len(direct_keys):
        present = {
            (row.get("pkn_source_node"), row.get("target_node"), str(row.get("interaction")))
            for row in direct_mapping_rows
        }
        raise SystemExit("A selected direct chemical edge lacks relaxed mapping evidence: " + str(sorted(direct_keys - present)))
    direct_fields = list(mapping_rows[0]) if mapping_rows else []
    write_csv(results / "network_direct_chemical_evidence.csv", direct_mapping_rows, direct_fields)

    nodes = {row["node_id"]: row for row in read_csv(nodes_path)}
    node_rows = []
    for row in att:
        meta = nodes.get(row["source"], {})
        output = dict(row)
        output["label"] = meta.get("label", "")
        output["node_type"] = meta.get("node_type", "")
        output["resource"] = meta.get("resource", "")
        node_rows.append(output)
    node_fields = list(att[0]) + ["label", "node_type", "resource"]
    write_csv(results / "network_nodes.csv", node_rows, node_fields)

    summary = {
        "analysis": "factor1_lipid_to_tf_relaxed_selected_export",
        "network_edges": len(sif),
        "network_nodes": len(att),
        "provenance_rows": len(evidence_rows),
        "direct_chemical_edges": len(direct_edges),
        "direct_mapping_evidence_rows": len(direct_mapping_rows),
        "pkn_evidence_md5": md5(evidence_path),
        "mapping_evidence_md5": md5(mapping_evidence),
        "interpretation": "Edge records preserve assembled PKN and relaxed-mapping provenance. They do not establish lipid identity, liver exposure, receptor occupancy, or causal regulation in this experiment.",
    }
    (results / "network_evidence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Exported provenance for {len(sif)} relaxed network edges to {results}")


if __name__ == "__main__":
    main()
