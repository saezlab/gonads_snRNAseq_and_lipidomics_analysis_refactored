#!/usr/bin/env python3
"""Audit a prospective GitHub commit without staging, deleting or committing.

Uses only Python's standard library and read-only Git commands. Outputs are
local audit artifacts, ignored by Git by default. Existing scientific files,
archives, manifests and Git's index are never modified.
"""
from __future__ import annotations

import argparse
from array import array
import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

from repository_paths import resolve_recorded_path

ROOT = Path(__file__).resolve().parents[2]
MIB = 1024 ** 2
CURRENT = "results/zonation_marker_update_20261116/"
CSV_MAX_BYTES = 50_000_000  # User policy is strict decimal MB, not MiB.
OBJECT_SUFFIXES = {".rds", ".rdata", ".rda", ".h5", ".hdf5", ".h5ad", ".h5seurat",
                   ".h5mu", ".loom", ".qs", ".qs2", ".pkl", ".pickle", ".mtx", ".zarr"}
ARCHIVE_SUFFIXES = {".gz", ".bz2", ".xz", ".zst", ".zip", ".tar", ".tgz"}
SINGLE_CELL_FORMATS = {".h5ad", ".h5seurat", ".h5mu", ".loom", ".zarr"}
SINGLE_CELL_NAME = re.compile(r"seurat|single.?cell|snrna|scrna|pseudobulk|anndata|cell.?ranger|(?:raw|filtered)_feature_bc_matrix", re.I)
TRANSPORT_ROLES = ("current_fcg_zonation", "current_gdx_zonation", "current_wdsd_zonation",
                   "wdsd_hepatocytes", "historical_wdsd_zonation", "wd_cell_types", "wd_zonation")
LOCAL = re.compile(
    r"(^|/)(?:\.DS_Store|\.Rhistory|\.Rapp\.history|\.RData|Rplots\.pdf)$"
    r"|(^|/)(?:__pycache__|\.pytest_cache|\.Rproj\.user|\.venv)(?:/|$)"
    r"|^\.local/R-library/|^\.claude/settings\.local\.json$|\.py[co]$"
    r"|\.tmp(?:\.|$)"
)
ABSOLUTE_LITERAL = re.compile(r"(?:/Users/|/Volumes/|/omics/|/home/)[^\s\"'<>)]*")
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]\n]+\]\((<[^>]+>|[^\s)]+)(?:\s+[\"'][^\n]*?[\"'])?\)")


def git(*args: str, binary: bool = False):
    result = subprocess.run(["git", *args], cwd=ROOT, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.stdout if binary else result.stdout.decode("utf-8", "surrogateescape")


def git_paths(*args: str) -> list[str]:
    return [x for x in git(*args, "-z").split("\0") if x]


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def relative(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def rna_input_policy(entries: list[dict], declared_models: set[str] | None = None) -> dict:
    """Classify candidate data only; never open or remove serialized objects."""
    forbidden_rna = []; oversized_csv = []; other_data_objects = []; single_cell_elsewhere = []
    obvious_anywhere = []; hdf5_review = []; archived_payloads = []
    declared_models = declared_models or set()
    filenames = {}
    for row in entries:
        path = Path(row["path"])
        filenames.setdefault(path.parent.as_posix(), set()).add(path.name.lower())
    for row in entries:
        name = row["path"]
        suffixes = {p.lower() for p in Path(name).suffixes}
        object_format = bool(suffixes & OBJECT_SUFFIXES) or any(p.lower().endswith(".zarr") for p in Path(name).parts)
        if object_format and name.startswith(("miscellaneous/older_code/", "miscellaneous/quick_tests/")) and name not in declared_models:
            archived_payloads.append(row)
        reason = ""
        if suffixes & SINGLE_CELL_FORMATS:
            reason = "Single-cell object/container format"
        elif object_format and SINGLE_CELL_NAME.search(Path(name).name):
            reason = "Serialized object with a single-cell/Seurat/pseudobulk input name"
        else:
            sibling_names = filenames[Path(name).parent.as_posix()]
            has_matrix = bool(sibling_names & {"matrix.mtx", "matrix.mtx.gz", "matrix.mtx.bz2"})
            has_barcodes = bool(sibling_names & {"barcodes.tsv", "barcodes.tsv.gz", "barcodes.tsv.bz2"})
            has_features = bool(sibling_names & {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"})
            raw_sparse_name = re.fullmatch(r"(?:matrix[.]mtx|barcodes[.]tsv|features[.]tsv|genes[.]tsv)(?:[.](?:gz|bz2))?", Path(name).name, re.I)
            if has_matrix and has_barcodes and has_features and raw_sparse_name:
                reason = "Raw barcode/feature sparse-matrix bundle"
        if reason:
            obvious_anywhere.append(dict(**row, reason=reason))
        if suffixes & {".h5", ".hdf5"} and not name.startswith("data/RNA/"):
            classification = "MOFA fitted model" if name.startswith("results/mofa/RNA/") else "Other HDF5 container; review provenance"
            hdf5_review.append(dict(**row, classification=classification))
        if not name.startswith("data/"):
            continue
        if name.startswith("data/RNA/"):
            if object_format or suffixes & ARCHIVE_SUFFIXES:
                forbidden_rna.append(row)
            if Path(name).suffix.lower() == ".csv" and row["bytes"] >= CSV_MAX_BYTES:
                oversized_csv.append(row)
        elif object_format:
            other_data_objects.append(row)
            if suffixes & SINGLE_CELL_FORMATS or SINGLE_CELL_NAME.search(name):
                single_cell_elsewhere.append(row)
    return dict(forbidden_rna_objects=forbidden_rna, rna_csv_at_or_above_50MB=oversized_csv,
                other_data_objects_for_review=other_data_objects,
                obvious_single_cell_objects_elsewhere=single_cell_elsewhere,
                obvious_single_cell_inputs_anywhere=obvious_anywhere, hdf5_containers_for_review=hdf5_review,
                undeclared_serialized_payloads_in_code_archives=archived_payloads)


def pseudobulk_transports(candidate: set[str]) -> dict:
    """Validate public CSV transports without reading their local source RDS."""
    registry = json.loads((ROOT/"scripts/config/analysis_release.json").read_text())
    paths = registry["paths"]
    transport_root = ROOT/paths.get("pseudobulk_csv_root", "data/RNA/pseudobulk_csv")
    datasets = []; files = []; layers = []; issues = []; declared_csv = set()
    seen_manifests = set()

    def demand(condition: bool, message: str) -> None:
        if not condition:
            raise ValueError(message)

    for role in TRANSPORT_ROLES:
        outcome = dict(dataset=role, manifest="", passed=False, features=None, profiles=None)
        datasets.append(outcome)
        try:
            key = role+"_transport"
            demand(key in paths, f"Missing release-registry key: {key}")
            manifest_path = ROOT/paths[key]
            if manifest_path.suffix.lower() != ".json":
                manifest_path /= "manifest.json"
            outcome["manifest"] = relative(manifest_path)
            demand(manifest_path.resolve() not in seen_manifests, "Distinct dataset roles must not share a transport manifest")
            seen_manifests.add(manifest_path.resolve())
            demand(manifest_path.resolve().is_relative_to(transport_root.resolve()), "Manifest is outside the registered transport root")
            demand(relative(manifest_path) in candidate, "Transport manifest is missing or excluded from the commit")
            manifest = json.loads(manifest_path.read_text())
            demand(manifest["format"] == "pseudobulk_csv" and manifest["format_version"] == 1, "Unsupported transport format/version")
            demand(manifest["individual_nuclei"] is False and manifest["profile_unit"]=="sample-level pseudobulk", "Transport must explicitly describe sample-level pseudobulk, not individual nuclei")
            demand(manifest["maximum_csv_bytes"] == CSV_MAX_BYTES, "Manifest CSV-byte limit differs from repository policy")
            profiles = manifest["profile_ids"]; features = manifest["features"]
            demand(profiles and len(set(profiles))==len(profiles) and all(isinstance(v,str) and v for v in profiles), "Profile IDs must be nonempty, unique strings")
            demand(features and len(set(features))==len(features) and all(isinstance(v,str) and v for v in features), "Feature IDs must be nonempty, unique strings")
            demand(re.fullmatch(r"[0-9a-f]{64}",manifest["source"]["sha256"]) is not None, "Missing source-object SHA256 provenance")
            outcome.update(features=len(features), profiles=len(profiles))
            descriptors = {}
            for descriptor in manifest["files"]:
                name = descriptor["file"]
                demand(name not in descriptors, f"Duplicate file descriptor: {name}")
                path = (manifest_path.parent/name).resolve()
                demand(path.is_relative_to(manifest_path.parent.resolve()) and path.suffix.lower()==".csv", f"Payload must be a CSV within its transport folder: {name}")
                demand(relative(path) in candidate, f"CSV is missing or excluded from the commit: {name}")
                actual_bytes = path.stat().st_size
                actual_sha = hashlib.sha256(path.read_bytes()).hexdigest()
                files.append(dict(dataset=role, role=descriptor["role"], file=relative(path),
                                  expected_bytes=descriptor["bytes"], actual_bytes=actual_bytes,
                                  sha256_matches=actual_sha==descriptor["sha256"], included=True))
                demand(actual_bytes < CSV_MAX_BYTES, f"CSV is not strictly below {CSV_MAX_BYTES} bytes: {name}")
                demand(actual_bytes == descriptor["bytes"] and actual_sha == descriptor["sha256"], f"CSV size/SHA256 mismatch: {name}")
                descriptors[name] = descriptor
                declared_csv.add(relative(path))
            demand("manifest.csv" in descriptors, "Companion manifest.csv index is not covered by file hashes")
            with (manifest_path.parent/"manifest.csv").open(newline="", encoding="utf-8-sig") as handle:
                index = list(csv.DictReader(handle))
            demand(len(index)==len(descriptors)-1 and {r["file"] for r in index}==set(descriptors)-{"manifest.csv"}, "Companion file index is incomplete or duplicated")
            for entry in index:
                descriptor = descriptors[entry["file"]]
                demand(entry["role"]==descriptor["role"] and int(entry["bytes"])==descriptor["bytes"] and entry["sha256"]==descriptor["sha256"], "Companion file-index descriptor mismatch")
                demand(entry["source_rds_sha256"]==manifest["source"]["sha256"] and int(entry["source_rds_bytes"])==manifest["source"]["bytes"], "Companion file-index source provenance mismatch")
            demand(set(manifest["layers"])=={"counts","data","scale.data"}, "All source layers counts/data/scale.data must be declared")
            for layer_name, layer in manifest["layers"].items():
                row_ids = layer["row_ids"]; column_ids = layer["column_ids"]
                demand(layer["dimensions"] == [len(row_ids),len(column_ids)], f"Layer dimensions disagree with IDs: {layer_name}")
                demand(len(set(row_ids))==len(row_ids) and set(row_ids).issubset(features), f"Layer row IDs are duplicated or outside source features: {layer_name}")
                demand(column_ids==profiles, f"Layer profile order differs from source profile IDs: {layer_name}")
                demand(re.fullmatch(r"[0-9a-f]{64}",layer["content_sha256"]) is not None, f"Missing layer numerical-content provenance: {layer_name}")
                demand(layer["content_encoding"]=="IEEE754-binary64-little-endian-column-major", "Unsupported matrix-content hash encoding")
                content = hashlib.sha256(); first_layer_record = len(layers)
                next_column = 1
                for chunk in layer["chunks"]:
                    name = chunk["file"]
                    demand(name in descriptors, f"Layer chunk is absent from the hashed file index: {name}")
                    for field in ("bytes","sha256"):
                        demand(chunk[field] == descriptors[name][field], f"Chunk/file-index {field} disagreement: {name}")
                    start, end = chunk["column_start"], chunk["column_end"]
                    demand(start==next_column and start<=end<=len(column_ids), f"Noncontiguous/overlapping layer column coverage: {name}")
                    demand(chunk["n_rows"]==len(row_ids) and chunk["n_columns"]==end-start+1, f"Incorrect declared chunk dimensions: {name}")
                    with (manifest_path.parent/name).open(newline="", encoding="utf-8-sig") as handle:
                        reader = csv.reader(handle); header = next(reader)
                        demand(header[1:]==column_ids[start-1:end], f"CSV header/profile order mismatch: {name}")
                        columns = [array("d") for _ in header[1:]]
                        count = 0
                        for count, row in enumerate(reader, start=1):
                            demand(count<=len(row_ids) and len(row)==len(header) and row[0]==row_ids[count-1], f"CSV row IDs/dimensions mismatch: {name},row{count+1}")
                            values = [float(value) for value in row[1:]]
                            demand(all(math.isfinite(value) for value in values), f"Nonfinite matrix value: {name},row{count+1}")
                            for column, value in zip(columns,values):
                                column.append(value)
                    demand(count==len(row_ids), f"CSV row count mismatch: {name}")
                    for column in columns:
                        demand(column.itemsize==8, "Runtime double is not binary64")
                        if sys.byteorder != "little":
                            column.byteswap()
                        content.update(column.tobytes())
                    layers.append(dict(dataset=role, layer=layer_name, file=relative(manifest_path.parent/name),
                                       rows=count, columns=end-start+1, profile_start=start, profile_end=end))
                    next_column = end+1
                demand(next_column==len(column_ids)+1, f"Layer chunks do not cover every profile: {layer_name}")
                demand(content.hexdigest()==layer["content_sha256"], f"Binary64 numerical-content SHA256 mismatch: {layer_name}")
                for record in layers[first_layer_record:]:
                    record.update(layer_content_sha256=content.hexdigest(), layer_content_sha256_matches=True)
            for attribute in ("metadata","identities"):
                descriptor = manifest[attribute]
                name = descriptor["file"]
                demand(name in descriptors, f"{attribute} is absent from the hashed file index")
                demand(descriptors[name]["role"]==attribute and all(descriptor[field]==descriptors[name][field] for field in ("bytes","sha256")), f"{attribute} descriptor disagrees with file index")
                with (manifest_path.parent/name).open(newline="", encoding="utf-8-sig") as handle:
                    reader = csv.reader(handle); header = next(reader); rows = list(reader)
                demand(all(len(row)==len(header) for row in rows), f"Malformed {attribute} CSV")
                demand([row[0] for row in rows]==profiles, f"{attribute} profile IDs/order do not match expression")
                if attribute=="metadata":
                    schema = descriptor["columns"]
                    demand(header[1:]==[c["name"] for c in schema] and descriptor["n_columns"]==len(schema) and descriptor["n_rows"]==len(rows), "Typed metadata schema/dimensions disagree with CSV")
                    demand(all(c["class"] and c["storage_type"] for c in schema), "Metadata type schema is incomplete")
                else:
                    demand(descriptor["names"]==profiles and header==["profile_id","identity"], "Identity names/schema disagree with profiles")
                    demand(all(row[1] in descriptor["levels"] or row[1]==manifest["na_token"] for row in rows), "Identity value is outside declared factor levels")
            outcome["passed"] = True
        except (KeyError,ValueError,TypeError,OSError,StopIteration) as error:
            issues.append(dict(dataset=role, issue=str(error)))
    extra_csv = sorted(p for p in candidate if p.startswith(relative(transport_root)+"/")
                       and Path(p).suffix.lower()==".csv" and p not in declared_csv)
    # A CSV outside a completed manifest must not silently enter the public tree.
    if all(row["passed"] for row in datasets) and extra_csv:
        issues.extend(dict(dataset="unregistered_csv", issue=p) for p in extra_csv)
    return dict(root=relative(transport_root), datasets=datasets, files=files, layers=layers, issues=issues,
                complete=all(row["passed"] for row in datasets) and not issues,
                numerical_content_note="File SHA256 and IEEE754 binary64 little-endian column-major matrix-content SHA256 are independently verified, together with numeric finiteness, dimensions and IDs.")


def markdown_links(files: list[str], candidate: set[str]) -> list[dict]:
    rows = []
    for name in files:
        if not name.lower().endswith(".md"):
            continue
        path = ROOT / name
        # Fenced examples are not rendered links.
        text = re.sub(r"```.*?```|~~~.*?~~~", "", path.read_text(errors="replace"), flags=re.S)
        for match in MARKDOWN_LINK.finditer(text):
            raw = match.group(1).strip("<>")
            target = urlsplit(raw)
            if target.scheme or not target.path:
                continue
            # Codex file links can include a one-based :line suffix. It does
            # not form part of the actual filename.
            destination = Path(re.sub(r":\d+$", "", unquote(target.path)))
            if not destination.is_absolute():
                destination = path.parent / destination
            destination = destination.resolve()
            resolved = relative(destination)
            state = "present_in_commit" if resolved in candidate else (
                "present_local_only" if destination.exists() else "missing")
            if state == "present_in_commit":
                continue
            # Frozen historical notes can retain old paths. They are reported,
            # not rewritten to imply that provenance was originally different.
            active = not (name.startswith(("miscellaneous/", "scripts/moon/archive/", "refactoring_docs/", "docs/research_plans/"))
                          or "/source_snapshots/" in name or re.search(r"/initial_panel_\d+/", name))
            rows.append(dict(document=name, link=raw, target=resolved, state=state,
                             active_document=active))
    return rows


def manifest_references(files: list[str], candidate: set[str]) -> list[dict]:
    rows = []
    path_column = re.compile(r"(?:^path$|^file$|^source$|^staged$|^generated_file$|^source_table$|_path$|_file$)")
    for name in files:
        p = Path(name)
        if p.suffix.lower() != ".csv" or not any(k in p.name.lower() for k in ("manifest", "provenance", "panel_trace", "sha256")):
            continue
        with (ROOT / name).open(encoding="utf-8-sig", errors="replace", newline="") as handle:
            for number, row in enumerate(csv.DictReader(handle), start=2):
                for column, raw in row.items():
                    if not column or not path_column.search(column) or not raw:
                        continue
                    for value in raw.split(";"):
                        value = value.strip()
                        if not value or "://" in value or not ("/" in value or Path(value).suffix):
                            continue
                        # The release registry maps old locations without
                        # rewriting immutable recorded paths in provenance.
                        path = resolve_recorded_path(value, root=ROOT).resolve()
                        resolved = relative(path)
                        state = "present_in_commit" if resolved in candidate else (
                            "present_local_only" if path.exists() else "missing_or_archived")
                        if state != "present_in_commit":
                            rows.append(dict(manifest=name, row=number, column=column,
                                             recorded_path=value, resolved_path=resolved, state=state))
    return rows


def history_sizes() -> dict:
    objects = git("rev-list", "--objects", "--all", binary=True)
    checked = subprocess.run(
        ["git", "cat-file", "--batch-check=%(objectname) %(objecttype) %(objectsize) %(rest)"],
        cwd=ROOT, input=objects, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    blobs = []
    for line in checked.stdout.decode("utf-8", "replace").splitlines():
        parts = line.split(" ", 3)
        if len(parts) >= 3 and parts[1] == "blob":
            blobs.append(dict(bytes=int(parts[2]), path=parts[3] if len(parts) > 3 else "", oid=parts[0]))
    return dict(blob_count=len(blobs), uncompressed_blob_bytes=sum(r["bytes"] for r in blobs),
                over_50MiB=sum(r["bytes"] > 50*MIB for r in blobs),
                over_100MiB=[r for r in blobs if r["bytes"] > 100*MIB])


def archive_recovery(head: str) -> dict:
    """Compare archived working bytes with pre-cleanup HEAD; write no Git objects."""
    destination = ROOT/"docs/provenance/archive_recovery.csv"
    if destination.exists():
        with destination.open(encoding="utf-8", newline="") as handle:
            recorded_heads = {r["pre_cleanup_HEAD"] for r in csv.DictReader(handle)}
        if len(recorded_heads) != 1:
            raise ValueError("Archive recovery ledger must retain one pre-cleanup commit")
        head = recorded_heads.pop()
        if not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", head):
            raise ValueError("Invalid pre-cleanup Git object ID")
        git("cat-file", "-e", head+"^{commit}")
    ledgers = ["docs/provenance/repository_migration_20261116.csv",
               "docs/provenance/causal_cleanup_20261116.csv"]
    moved = {}
    for ledger in ledgers:
        with (ROOT/ledger).open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                key = (row["old_path"], row["new_path"])
                if key in moved and moved[key]["recorded_sha256"] != row["sha256"]:
                    raise ValueError(f"Conflicting migration checksums: {key}")
                moved[key] = dict(old_path=row["old_path"], new_path=row["new_path"],
                                  action=row.get("action", "move"), recorded_sha256=row["sha256"],
                                  source_ledger=ledger)
    head_blobs = {}
    for entry in git("ls-tree", "-r", "-z", head).split("\0"):
        if entry:
            attributes, path = entry.split("\t", 1)
            _, kind, oid = attributes.split()
            if kind == "blob":
                head_blobs[path] = oid
    comparable = [r for r in moved.values() if r["old_path"] in head_blobs and (ROOT/r["new_path"]).is_file()]
    # Quoted stdin paths handle spaces and quotes without invoking a shell.
    input_paths = "".join(json.dumps(str(ROOT/r["new_path"]), ensure_ascii=False)+"\n" for r in comparable)
    hashed = subprocess.run(["git", "hash-object", "--no-filters", "--stdin-paths"], cwd=ROOT,
                            input=input_paths.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    working_oids = hashed.stdout.decode().splitlines()
    if len(working_oids) != len(comparable):
        raise ValueError("Git returned an unexpected number of read-only blob hashes")
    raw_oids = {r["new_path"]:oid for r,oid in zip(comparable, working_oids)}
    rows = []
    for row in moved.values():
        path = ROOT/row["new_path"]
        hasher = hashlib.sha256()
        exists = path.is_file()
        if exists:
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(8*1024**2), b""):
                    hasher.update(block)
        actual_sha = hasher.hexdigest() if exists else ""
        original_oid = head_blobs.get(row["old_path"], "")
        working_oid = raw_oids.get(row["new_path"], "")
        identical = bool(original_oid and working_oid == original_oid)
        rows.append(dict(**row, bytes=path.stat().st_size if exists else "", working_sha256=actual_sha,
                         working_sha256_matches_migration=exists and actual_sha==row["recorded_sha256"],
                         pre_cleanup_HEAD=head, head_blob_oid=original_oid, working_blob_oid=working_oid,
                         git_byte_identical=identical,
                         recovery_source="git_history_byte_identical" if identical else "local_archive_required",
                         git_object_spec=f"{head}:{row['old_path']}" if identical else ""))
    write_csv(destination, rows, ["old_path", "new_path", "action", "bytes", "working_sha256", "recorded_sha256",
              "working_sha256_matches_migration", "pre_cleanup_HEAD", "head_blob_oid", "working_blob_oid",
              "git_byte_identical", "recovery_source", "git_object_spec", "source_ledger"])
    return dict(path=relative(destination), pre_cleanup_HEAD=head, rows=len(rows),
                git_history_byte_identical=sum(r["git_byte_identical"] for r in rows),
                local_archive_required=sum(not r["git_byte_identical"] for r in rows),
                sha256_mismatches=sum(not r["working_sha256_matches_migration"] for r in rows))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "miscellaneous/local/github_readiness")
    parser.add_argument("--strict", action="store_true", help="Exit1 for size, local-metadata or canonical-input blockers")
    parser.add_argument("--skip-history", action="store_true", help="Skip the read-only reachable-blob size scan")
    parser.add_argument("--write-archive-recovery", action="store_true",
                        help="Write docs/provenance/archive_recovery.csv after exact working-byte/HEAD comparison")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT):
        parser.error("Audit output must be inside this checkout")
    # Snapshot Git's index and HEAD so the audit can assert it never stages.
    index_before = git("ls-files", "--stage", "-z", binary=True)
    head_before = git("rev-parse", "HEAD").strip()
    recovery = archive_recovery(head_before) if args.write_archive_recovery else None
    tracked = set(git_paths("ls-files"))
    untracked = set(git_paths("ls-files", "--others", "--exclude-standard"))
    present_tracked = {p for p in tracked if (ROOT / p).is_file()}
    candidate = present_tracked | untracked
    entries = [dict(path=p, bytes=(ROOT/p).stat().st_size,
                    git_state="tracked" if p in tracked else "untracked_nonignored")
               for p in sorted(candidate) if (ROOT/p).is_file()]
    over100 = [r for r in entries if r["bytes"] > 100*MIB]
    over50 = [r for r in entries if r["bytes"] > 50*MIB]
    local = [r for r in entries if LOCAL.search(r["path"])]
    frozen_manifest = ROOT / "docs/provenance/active_frozen_dependencies.json"
    declared_frozen = json.loads(frozen_manifest.read_text()) if frozen_manifest.exists() else {}
    rna_policy = rna_input_policy(entries, {row["path"] for row in declared_frozen.get("files", [])})
    transports = pseudobulk_transports(candidate)
    rna_csv_files = [r for r in entries if r["path"].startswith("data/RNA/") and Path(r["path"]).suffix.lower()==".csv"]
    largest_rna_csv = max(rna_csv_files, key=lambda row:row["bytes"], default=None)
    frozen = []
    for row in declared_frozen.get("files", []):
        p = ROOT / row["path"]
        exists = p.is_file()
        frozen.append(dict(path=row["path"], exists=exists, included=row["path"] in candidate,
                           expected_bytes=row["bytes"], actual_bytes=p.stat().st_size if exists else None,
                           sha256_matches=exists and hashlib.sha256(p.read_bytes()).hexdigest()==row["sha256"]))
    frozen_complete = (bool(frozen) and len(frozen)==declared_frozen.get("total_files")
                       and len({r["path"] for r in frozen})==len(frozen)
                       and relative(frozen_manifest) in candidate
                       and all(r["exists"] and r["included"] and r["sha256_matches"]
                               and r["actual_bytes"]==r["expected_bytes"] for r in frozen))
    links = markdown_links([r["path"] for r in entries], candidate)
    references = manifest_references([r["path"] for r in entries], candidate)
    absolute = []
    for row in entries:
        name = row["path"]
        if name == relative(Path(__file__).resolve()):
            continue  # The regex definitions above are not runtime paths.
        if not name.startswith(("scripts/", "tests/")) or name.startswith(("scripts/legacy/", "scripts/moon/archive/")):
            continue
        if Path(name).suffix not in (".R", ".r", ".py", ".json"):
            continue
        for line_no, line in enumerate((ROOT/name).read_text(errors="replace").splitlines(), 1):
            for match in ABSOLUTE_LITERAL.finditer(line):
                absolute.append(dict(file=name, line=line_no, literal=match.group()))
    history = {} if args.skip_history else history_sizes()
    checks = [dict(check="No prospective regular Git file exceeds100MiB", passed=not over100),
              dict(check="No local metadata or caches in prospective tree", passed=not local),
              dict(check="No RNA input objects or compressed containers in candidate tree", passed=not rna_policy["forbidden_rna_objects"]),
              dict(check="Every RNA input CSV is strictly below 50000000 bytes", passed=not rna_policy["rna_csv_at_or_above_50MB"]),
              dict(check="No obvious single-cell objects elsewhere under data", passed=not rna_policy["obvious_single_cell_objects_elsewhere"]),
              dict(check="No obvious single-cell containers or raw sparse bundles anywhere in candidate tree", passed=not rna_policy["obvious_single_cell_inputs_anywhere"]),
              dict(check="No undeclared serialized payloads hidden in archived code/quick-test areas", passed=not rna_policy["undeclared_serialized_payloads_in_code_archives"]),
              dict(check="All seven CSV transports have complete included/hash-verified layers and metadata", passed=transports["complete"]),
              dict(check="All declared active MOON frozen dependencies included with exact hashes", passed=frozen_complete),
              dict(check="All reachable historical blobs below100MiB", passed=not history.get("over_100MiB", [])),
              dict(check="Git index unchanged by audit", passed=index_before == git("ls-files", "--stage", "-z", binary=True)),
              dict(check="HEAD unchanged by audit", passed=head_before == git("rev-parse", "HEAD").strip())]
    if recovery is not None:
        checks.append(dict(check="All archive working SHA256 hashes match migration ledgers", passed=recovery["sha256_mismatches"]==0))
    # Links/manifests can legitimately reference regenerable, archived or private
    # files. Surface each one separately; never pretend they ship in a clone.
    summary = dict(
        generated_utc=dt.datetime.now(dt.timezone.utc).isoformat(), checkout=str(ROOT),
        head=head_before, ready_for_review=all(r["passed"] for r in checks),
        scope="Working-tree candidate for broad git add; no staging/commit/push performed",
        tracked_paths=len(tracked), tracked_present=len(present_tracked),
        tracked_deleted=len(tracked-present_tracked), untracked_nonignored=len(untracked),
        candidate_files=len(entries), candidate_bytes=sum(r["bytes"] for r in entries),
        candidate_over_50MiB=over50, candidate_over_100MiB=over100, local_metadata_or_caches=local,
        rna_input_policy=rna_policy, pseudobulk_transports=transports,
        rna_csv_file_count=len(rna_csv_files), largest_rna_input_csv=largest_rna_csv,
        active_frozen_dependencies=frozen, archive_recovery=recovery, history=history, checks=checks,
        nonshipping_markdown_links=len(links), active_missing_markdown_links=sum(r["active_document"] and r["state"]=="missing" for r in links),
        nonshipping_manifest_references=len(references), machine_specific_source_literals=len(absolute),
        dependency_files=[p for p in ("DESCRIPTION", "renv.lock", "requirements.txt", "pyproject.toml", "environment.yml") if (ROOT/p).exists()],
        interpretation="Size/cache/input checks are commit gates. Link, manifest and absolute-path ledgers require review; immutable historical provenance and optional private inputs can intentionally be absent from a public clone.")
    output.mkdir(parents=True, exist_ok=True)
    (output/"github_readiness.json").write_text(json.dumps(summary, indent=2)+"\n")
    write_csv(output/"candidate_files.csv", entries, ["path", "bytes", "git_state"])
    write_csv(output/"checks.csv", checks, ["check", "passed"])
    write_csv(output/"single_cell_input_candidates.csv", rna_policy["obvious_single_cell_inputs_anywhere"],
              ["path", "bytes", "git_state", "reason"])
    write_csv(output/"archived_serialized_payload_candidates.csv", rna_policy["undeclared_serialized_payloads_in_code_archives"],
              ["path", "bytes", "git_state"])
    write_csv(output/"hdf5_container_review.csv", rna_policy["hdf5_containers_for_review"],
              ["path", "bytes", "git_state", "classification"])
    write_csv(output/"pseudobulk_transport_files.csv", transports["files"],
              ["dataset", "role", "file", "expected_bytes", "actual_bytes", "sha256_matches", "included"])
    write_csv(output/"pseudobulk_transport_layers.csv", transports["layers"],
              ["dataset", "layer", "file", "rows", "columns", "profile_start", "profile_end", "layer_content_sha256", "layer_content_sha256_matches"])
    write_csv(output/"pseudobulk_transport_issues.csv", transports["issues"], ["dataset", "issue"])
    write_csv(output/"active_frozen_dependencies.csv", frozen,
              ["path", "exists", "included", "expected_bytes", "actual_bytes", "sha256_matches"])
    write_csv(output/"nonshipping_markdown_links.csv", links, ["document", "link", "target", "state", "active_document"])
    write_csv(output/"nonshipping_manifest_references.csv", references, ["manifest", "row", "column", "recorded_path", "resolved_path", "state"])
    write_csv(output/"machine_specific_source_literals.csv", absolute, ["file", "line", "literal"])
    lines = ["# GitHub readiness audit", "", f"Ready for commit review: **{summary['ready_for_review']}**.", "",
             f"Prospective tree: {len(entries):,} files; {summary['candidate_bytes']:,} bytes ({summary['candidate_bytes']/1024**3:.3f} GiB).",
             f"Tracked: {len(tracked):,} paths, {len(present_tracked):,} present; untracked and nonignored: {len(untracked):,}.", "",
             "| Check | Passed |", "| --- | --- |"]
    lines += [f"| {r['check']} | {r['passed']} |" for r in checks]
    lines += ["", f"Files above50MiB: {len(over50)}; above100MiB: {len(over100)}.",
              f"CSV-only RNA transports complete: {transports['complete']}; strict per-CSV limit: {CSV_MAX_BYTES:,} bytes.",
              f"Largest candidate RNA input CSV: {largest_rna_csv['bytes']:,} bytes ({largest_rna_csv['path']})." if largest_rna_csv else "No RNA input CSV found.",
              f"Declared archived dependencies retained for active MOON: {len(frozen)}; {sum(r['expected_bytes'] for r in frozen):,} bytes.",
              f"Nonshipping documentation links: {len(links)} ({summary['active_missing_markdown_links']} missing in active documents).",
              f"Nonshipping manifest references: {len(references)}; machine-specific source literals: {len(absolute)}.", "",
              "Review the companion CSVs. Files present only locally may be generated results, archived snapshots or optional private manuscript inputs. Missing historical paths are preserved as provenance, not silently rewritten.", "",
              "No analysis, staging, commit, push, deletion or history rewrite was performed."]
    (output/"README.md").write_text("\n".join(lines)+"\n")
    print(f"GitHub audit: {len(entries):,} candidate files, {summary['candidate_bytes']:,} bytes; "
          f"{len(over100)} above100MiB; {len(local)} local metadata/cache files; ready={summary['ready_for_review']}.")
    print(f"Reports: {relative(output)}")
    return int(args.strict and not summary["ready_for_review"])


if __name__ == "__main__":
    raise SystemExit(main())
