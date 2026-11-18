#!/usr/bin/env python3
"""Attach edge-level PKN provenance to strict relaxed lipid-to-TF score-layer outputs."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ANALYSIS = ROOT / "results/causal_modelling/factor1_lipid_to_tf_relaxed_regulatory"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path: Path) -> List[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: List[dict], fields: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def edge_key(row: dict) -> Tuple[str, str, str]:
    return (row["source"], row["target"], str(row["interaction"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, default=DEFAULT_ANALYSIS)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    analysis = args.analysis.resolve()
    layers = analysis / "score_layer_constrained"
    output = (args.output or layers / "evidence_audit").resolve()
    if not layers.is_dir():
        raise SystemExit(f"Strict score-layer output is absent: {layers}")
    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"Output directory is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((analysis / "run_manifest.json").read_text())
    pkn_path = ROOT / manifest["config"]["pkn"]
    evidence_path = pkn_path.parent / "edge_evidence.jsonl.gz"
    if not evidence_path.is_file():
        raise SystemExit(f"PKN evidence ledger is absent: {evidence_path}")

    solution_files = sorted(layers.glob("solutions/*/network_SIF.csv"))
    if not solution_files:
        raise SystemExit("No materialized score-layer solution SIFs were found")
    by_solution: Dict[str, List[dict]] = {}
    keys = set()
    for sif in solution_files:
        solution = sif.parent.name
        rows = read_csv(sif)
        by_solution[solution] = rows
        keys.update(edge_key(row) for row in rows)

    # Several assembled evidence records can support the same retained signed
    # PKN triple. Preserve all of them rather than choosing one arbitrarily.
    provenance: Dict[Tuple[str, str, str], List[dict]] = {}
    with gzip.open(evidence_path, "rt") as handle:
        for line in handle:
            record = json.loads(line)
            key = (record["source"], record["target"], str(record["interaction"]))
            if key in keys:
                provenance.setdefault(key, []).append(record)
    missing = sorted(keys.difference(provenance))
    if missing:
        raise SystemExit(f"Selected score-layer edges absent from provenance ledger: {missing[:5]}")

    evidence_fields = [
        "solution", "source", "target", "interaction", "consistency",
        "first_variant", "resource", "resources", "mechanism", "references",
        "source_accession", "target_accession", "species_evidence", "evidence_tier",
        "source_version", "original_edge_id",
    ]
    all_rows: List[dict] = []
    for solution, edges in by_solution.items():
        records: List[dict] = []
        for edge in edges:
            for evidence in provenance[edge_key(edge)]:
                record = dict(evidence)
                record["solution"] = solution
                record["consistency"] = edge.get("consistency", "")
                records.append(record)
                all_rows.append(record)
        write_csv(layers / "solutions" / solution / "network_edge_evidence.csv", records, evidence_fields)
    write_csv(output / "all_score_layer_network_edge_evidence.csv", all_rows, evidence_fields)
    summary = {
        "analysis": str(analysis),
        "pkn_evidence": str(evidence_path),
        "solutions": len(by_solution),
        "unique_selected_edges": len(keys),
        "edge_rows_written": len(all_rows),
        "pkn_evidence_sha256": sha256(evidence_path),
        "script_sha256": sha256(Path(__file__).resolve()),
        "interpretation": "Each retained score-layer edge is linked to its assembled PKN evidence. Evidence provenance does not establish ligand exposure, tissue context, or causal regulation in this experiment.",
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Score-layer network edge evidence\n\n"
        "This audit attaches the assembled PKN evidence ledger to every edge of each "
        "materialized layer-constrained relaxed lipid-to-TF solution. The output is "
        "provenance, not a new fit or independent biological validation.\n"
    )
    print(
        f"Attached provenance for {len(keys)} unique edges across {len(by_solution)} "
        f"score-layer solutions: {output}"
    )


if __name__ == "__main__":
    main()
