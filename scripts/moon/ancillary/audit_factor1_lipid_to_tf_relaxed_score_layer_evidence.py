#!/usr/bin/env python3
"""Attach assembled PKN provenance to strict relaxed lipid-to-TF search solutions.

The parameter search deliberately stores only the signed network selected by the
strict score-layer extractor.  This ancillary audit resolves every retained
source/target/sign triple against the frozen PKN evidence ledger.  It does not
re-fit MOON or rank biological hypotheses.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, List, Mapping, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SEARCH = ROOT / "results/causal_modelling/factor1_lipid_to_tf_relaxed_parameter_search"


def hash_file(path: Path, algorithm: str) -> str:
    """Return a streaming checksum without loading large PKN ledgers into memory."""
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(f"Could not read JSON file {path}: {error}") from error
    if not isinstance(value, dict):
        raise SystemExit(f"Expected a JSON object in {path}")
    return value


def read_csv(path: Path) -> List[dict]:
    try:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))
    except OSError as error:
        raise SystemExit(f"Could not read CSV file {path}: {error}") from error


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    """Write a deterministic, rectangular CSV and preserve empty tables."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: csv_value(row.get(field, "")) for field in fields})


def csv_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def repo_path(value: str) -> Path:
    """Resolve frozen repository-relative paths recorded by the R workflows."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else ROOT / path


def canonical_interaction(value: Any) -> str:
    """Make signed interactions comparable across R CSV output and JSON ledger."""
    text = str(value).strip()
    try:
        numeric = float(text)
    except ValueError:
        return text
    if numeric.is_integer():
        return str(int(numeric))
    return text


def edge_key(row: Mapping[str, Any]) -> Tuple[str, str, str]:
    missing = [field for field in ("source", "target", "interaction") if not str(row.get(field, "")).strip()]
    if missing:
        raise SystemExit(f"An edge is missing required SIF columns {missing}: {dict(row)}")
    return (str(row["source"]), str(row["target"]), canonical_interaction(row["interaction"]))


def load_manifest(path: Path) -> Dict[str, str]:
    if not path.is_file():
        raise SystemExit(f"Parameter-search dependency manifest is absent: {path}")
    rows = read_csv(path)
    if not rows or not {"path", "md5"}.issubset(rows[0]):
        raise SystemExit(f"Dependency manifest has no path/md5 columns: {path}")
    return {str(row["path"]): str(row["md5"]).lower() for row in rows}


def identify_base_config(search: Path) -> Tuple[Path, dict, Path, dict]:
    """Return the copied search config plus the frozen base analysis config.

    `analysis_config.json` is copied into every search result by
    `relaxed_write_common_inputs()`.  It is the primary source because it makes
    an archived result self-describing.  The parameter-search config remains a
    second explicit link to the original source config and is checked when it is
    available.
    """
    analysis_config_path = search / "analysis_config.json"
    search_config_path = search / "parameter_search_config.json"
    if not analysis_config_path.is_file():
        raise SystemExit(f"Search analysis config is absent: {analysis_config_path}")
    if not search_config_path.is_file():
        raise SystemExit(f"Copied parameter-search config is absent: {search_config_path}")

    base_config = read_json(analysis_config_path)
    search_config = read_json(search_config_path)
    if not isinstance(base_config.get("pkn"), str) or not base_config["pkn"]:
        raise SystemExit(f"No PKN path is recorded in {analysis_config_path}")
    if not isinstance(search_config.get("base_config"), str) or not search_config["base_config"]:
        raise SystemExit(f"No base_config is recorded in {search_config_path}")

    original_base_path = repo_path(search_config["base_config"])
    if original_base_path.is_file():
        original_base = read_json(original_base_path)
        original_pkn = original_base.get("pkn")
        if original_pkn != base_config["pkn"]:
            raise SystemExit(
                "The copied search analysis config and its recorded base config disagree on the PKN: "
                f"{analysis_config_path} vs {original_base_path}"
            )
    return analysis_config_path, base_config, search_config_path, search_config


def locate_evidence_ledger(base_config: Mapping[str, Any], pkn_path: Path) -> Path:
    """Use the frozen input ledger path when present; otherwise derive it from the PKN."""
    expected_inputs = base_config.get("input_md5", {})
    if not isinstance(expected_inputs, dict):
        raise SystemExit("analysis_config.json input_md5 must be an object")
    candidates = [repo_path(path) for path in expected_inputs if Path(path).name == "edge_evidence.jsonl.gz"]
    if not candidates:
        candidates = [pkn_path.parent / "edge_evidence.jsonl.gz"]
    if len(candidates) != 1:
        raise SystemExit(f"Could not identify one PKN evidence ledger from configuration: {candidates}")
    evidence_path = candidates[0]
    if not evidence_path.is_file():
        raise SystemExit(f"PKN evidence ledger is absent: {evidence_path}")
    return evidence_path


def verify_frozen_inputs(
    base_config: Mapping[str, Any], manifest: Mapping[str, str], pkn_path: Path, evidence_path: Path
) -> Dict[str, str]:
    """Verify that the PKN and provenance ledger match both frozen hash records."""
    expected_inputs = base_config.get("input_md5", {})
    paths = (pkn_path, evidence_path)
    actual: Dict[str, str] = {}
    for path in paths:
        relative = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        expected_config = str(expected_inputs.get(relative, "")).lower()
        expected_manifest = str(manifest.get(relative, "")).lower()
        if not expected_config:
            raise SystemExit(f"Frozen analysis_config.json has no MD5 for required input: {relative}")
        if not expected_manifest:
            raise SystemExit(f"Search dependency manifest has no MD5 for required input: {relative}")
        if expected_config != expected_manifest:
            raise SystemExit(
                f"Frozen PKN hash disagreement for {relative}: analysis config={expected_config}, "
                f"search manifest={expected_manifest}"
            )
        observed = hash_file(path, "md5")
        if observed != expected_config:
            raise SystemExit(
                f"Frozen PKN input changed since the search for {relative}: expected {expected_config}, observed {observed}"
            )
        actual[relative] = observed
    return actual


def iter_provenance(path: Path, selected: set[Tuple[str, str, str]]) -> Iterable[Tuple[Tuple[str, str, str], dict]]:
    """Yield only ledger records whose signed triple occurs in a materialized SIF."""
    try:
        with gzip.open(path, "rt") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                try:
                    record = json.loads(line)
                except json.JSONDecodeError as error:
                    raise SystemExit(f"Invalid JSONL record in {path} at line {line_number}: {error}") from error
                if not isinstance(record, dict):
                    raise SystemExit(f"Non-object provenance record in {path} at line {line_number}")
                key = edge_key(record)
                if key in selected:
                    yield key, record
    except OSError as error:
        raise SystemExit(f"Could not read compressed PKN evidence ledger {path}: {error}") from error


def solution_metadata(solution_dir: Path) -> dict:
    audit_path = solution_dir / "extraction_audit.json"
    if not audit_path.is_file():
        return {}
    return read_json(audit_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--search",
        type=Path,
        default=DEFAULT_SEARCH,
        help="Parameter-search result directory (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="New audit directory (default: <search>/evidence_audit)",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate frozen inputs and provenance coverage without writing any files",
    )
    args = parser.parse_args()

    search = args.search.expanduser().resolve()
    if not search.is_dir():
        raise SystemExit(f"Parameter-search result directory is absent: {search}")
    output = (args.output.expanduser() if args.output else search / "evidence_audit").resolve()
    if not args.validate_only and output.exists() and any(output.iterdir()):
        raise SystemExit(f"Audit output directory is not empty: {output}")

    analysis_config_path, base_config, search_config_path, _ = identify_base_config(search)
    manifest = load_manifest(search / "dependency_manifest.csv")
    pkn_path = repo_path(base_config["pkn"])
    if not pkn_path.is_file():
        raise SystemExit(f"PKN recorded by the search is absent: {pkn_path}")
    evidence_path = locate_evidence_ledger(base_config, pkn_path)
    frozen_md5 = verify_frozen_inputs(base_config, manifest, pkn_path, evidence_path)

    solution_files = sorted(search.glob("runs/*/solutions/*/network_SIF.csv"))
    if not solution_files:
        raise SystemExit(
            "No materialized strict solution SIFs were found below "
            f"{search / 'runs'}. Run the parameter search before this provenance audit."
        )

    by_solution: Dict[str, List[dict]] = {}
    selected_keys: set[Tuple[str, str, str]] = set()
    for sif_path in solution_files:
        solution_ref = str(sif_path.parent.relative_to(search))
        sif_rows = read_csv(sif_path)
        if not sif_rows:
            raise SystemExit(f"Materialized strict solution has an empty SIF: {sif_path}")
        by_solution[solution_ref] = sif_rows
        selected_keys.update(edge_key(row) for row in sif_rows)

    provenance: DefaultDict[Tuple[str, str, str], List[dict]] = defaultdict(list)
    for key, record in iter_provenance(evidence_path, selected_keys):
        provenance[key].append(record)
    missing = sorted(selected_keys.difference(provenance))
    if missing:
        preview = "; ".join(" -> ".join(key) for key in missing[:5])
        raise SystemExit(
            "A materialized strict solution contains edge(s) absent from the frozen PKN evidence ledger: " + preview
        )

    if args.validate_only:
        print(
            f"Validated provenance for {len(selected_keys)} signed edges across {len(by_solution)} "
            f"materialized strict solutions. No files written."
        )
        return

    output.mkdir(parents=True, exist_ok=True)
    provenance_fields = sorted({field for records in provenance.values() for row in records for field in row})
    context_fields = [
        "solution_ref",
        "run_id",
        "topology_hash",
        "sif_edge_index",
        "sif_consistency",
        "sif_source",
        "sif_target",
        "sif_interaction",
    ]
    evidence_fields = context_fields + [field for field in provenance_fields if field not in context_fields]
    all_rows: List[dict] = []
    solution_rows: List[dict] = []
    for solution_ref, sif_rows in by_solution.items():
        solution_dir = search / solution_ref
        audit = solution_metadata(solution_dir)
        records: List[dict] = []
        for edge_index, edge in enumerate(sif_rows, start=1):
            key = edge_key(edge)
            context = {
                "solution_ref": solution_ref,
                "run_id": audit.get("run_id", solution_dir.parents[1].name),
                "topology_hash": audit.get("topology_hash", ""),
                "sif_edge_index": edge_index,
                "sif_consistency": edge.get("consistency", ""),
                "sif_source": key[0],
                "sif_target": key[1],
                "sif_interaction": key[2],
            }
            for evidence in provenance[key]:
                records.append({**context, **evidence})
        per_solution_path = solution_dir / "network_edge_evidence.csv"
        write_csv(per_solution_path, records, evidence_fields)
        all_rows.extend(records)
        solution_rows.append(
            {
                "solution_ref": solution_ref,
                "run_id": audit.get("run_id", solution_dir.parents[1].name),
                "topology_hash": audit.get("topology_hash", ""),
                "sif_edges": len(sif_rows),
                "unique_signed_edges": len({edge_key(row) for row in sif_rows}),
                "provenance_rows": len(records),
                "per_solution_evidence_csv": str(per_solution_path.relative_to(search)),
                "sif_md5": hash_file(solution_dir / "network_SIF.csv", "md5"),
            }
        )

    write_csv(output / "all_score_layer_network_edge_evidence.csv", all_rows, evidence_fields)
    write_csv(
        output / "per_solution_summary.csv",
        solution_rows,
        [
            "solution_ref", "run_id", "topology_hash", "sif_edges", "unique_signed_edges",
            "provenance_rows", "per_solution_evidence_csv", "sif_md5",
        ],
    )
    summary = {
        "search": str(search),
        "analysis_config": str(analysis_config_path),
        "parameter_search_config": str(search_config_path),
        "pkn": str(pkn_path),
        "pkn_evidence_ledger": str(evidence_path),
        "materialized_solutions": len(by_solution),
        "unique_selected_signed_edges": len(selected_keys),
        "provenance_rows_written": len(all_rows),
        "frozen_input_md5": frozen_md5,
        "pkn_evidence_sha256": hash_file(evidence_path, "sha256"),
        "script_sha256": hash_file(Path(__file__).resolve(), "sha256"),
        "interpretation": (
            "Each retained strict score-layer edge is linked to the assembled PKN evidence ledger. "
            "This provenance audit does not establish lipid exposure, tissue context, receptor engagement, "
            "or causal regulation in this experiment."
        ),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Strict relaxed lipid-to-TF network edge evidence\n\n"
        "This ancillary audit attaches the frozen assembled-PKN provenance ledger to every edge of every "
        "materialized strict score-layer solution. `all_score_layer_network_edge_evidence.csv` aggregates "
        "all records, while each solution directory contains a matching `network_edge_evidence.csv`. "
        "The audit is provenance only; it is not a MOON refit or independent biological validation.\n"
    )
    print(
        f"Attached provenance for {len(selected_keys)} signed edges across {len(by_solution)} "
        f"materialized strict solutions: {output}"
    )


if __name__ == "__main__":
    main()
