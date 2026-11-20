#!/usr/bin/env python3
"""Stage reproducible Workflow 17 tables and immutable manuscript references."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_path, resolve_recorded_path


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_ROOT = analysis_path("current_input_root")
TABLE_NOTE = (
    "Updated source table; preserve independent S11/S12. Submission formatting may differ."
)
TABLES = [
    ("S10", "fcg05/ttop_full_long.csv", "TableS10_zonation_differential_results.csv"),
    ("S13", "fcg05/sex_specific_zonation_genes_summary.csv", "TableS13_selected_zonation_genes.csv"),
    ("S13", "fcg05/t_table_zonation.csv", "TableS13_all_gene_t_statistics.csv"),
    ("S15", "fcg07/zonation_network_validation_scores.csv", "TableS15_zonation_network_validation_scores.csv"),
    ("S15", "fcg07/Mtmr11_male_PP_prespecified_source.csv", "TableS15_Mtmr11_male_PP_matched_lipid_source.csv"),
]
for _stem in (
    "cormat_Ovaries_portal", "cormat_ova_central", "cormat_testes_central", "cormat_testes_portal"
):
    for _kind in ("ATT", "SIF"):
        _name = f"{_stem}_full_{_kind}.csv"
        TABLES.append(("S14", f"fcg05/networks/{_name}", f"TableS14_{_name}"))

# Origins, roles and frozen hashes transcribed from the completed run's
# inputs/preprocessing_reference/provenance_manifest.csv. Existing per-run
# provenance is authoritative when refreshing a completed output directory.
REFERENCES = [
    ("preprocessing_02_fcg.R", "217fc9a01f0cab0b311f20f825c110e28172933286b63ffb385d76f4f9e6bc52", "supplied preprocessing source"),
    ("preprocessing_04_wdsd.R", "06232b157720856d009e46890ecf35ff5e968bd9e18ca512427f079a3b3f7769", "supplied preprocessing source"),
    ("preprocessing_03_gdx.R", "1cc9a73de59662aa540b9163124c094a5ae4076f196b34cdcb1f2cf7a8f32db4", "supplied preprocessing source"),
    ("Zonation_Genes.csv", "54f0c106a8bb9c026f239ae9546aa0da2fdaff253484246add5342fba918cbee", "candidate marker list; exact HPC runtime identity unverified"),
    ("Cleland and Dugourd et al.docx", "ead9393bd639361f5a999b29c715c0a4a3deb182c7080cc6b83aafc01ff1c244", "manuscript before zonation update"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_copy(source: Path, target: Path, expected: str) -> bool:
    """Copy exact bytes only when needed and verify before atomic replacement."""
    if target.exists() and sha256(target) == expected:
        return False
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    os.close(fd)
    try:
        shutil.copyfile(source, temporary)
        if sha256(Path(temporary)) != expected:
            raise RuntimeError(f"Source changed during copy: {source}")
        os.replace(temporary, target)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return True


def write_manifest(path: Path, fields: list[str], rows: list[dict]) -> None:
    import io

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    content = stream.getvalue().encode("utf-8")
    if path.exists() and path.read_bytes() == content:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as output:
            output.write(content)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def reference_rows(output: Path, source_root: Path | None) -> list[dict]:
    directory = output / "inputs/preprocessing_reference"
    manifest = directory / "provenance_manifest.csv"
    if manifest.exists():
        with manifest.open(encoding="utf-8-sig", newline="") as source:
            rows = list(csv.DictReader(source))
        expected_names = {name for name, _, _ in REFERENCES}
        names = [Path(row["snapshot"]).name for row in rows]
        if len(names) != len(expected_names) or set(names) != expected_names:
            raise RuntimeError(f"Unexpected frozen reference manifest: {manifest}")
        for row in rows:
            name = Path(row["snapshot"]).name
            row["snapshot"] = str(directory / name)
            if source_root is not None and name.startswith("preprocessing_"):
                row["source"] = str(source_root / "james_scripts_20261104" / name)
            elif not Path(row["source"]).exists():
                row["source"] = str(resolve_recorded_path(row["source"]))
        return rows
    origin = source_root or DEFAULT_SOURCE_ROOT
    rows = []
    for name, digest, role in REFERENCES:
        if name.startswith("preprocessing_"):
            source = origin / "james_scripts_20261104" / name
        elif name == "Zonation_Genes.csv":
            source = origin / name
        else:
            source = REPOSITORY / "manuscript" / name
        rows.append(dict(source=str(source), snapshot=str(directory / name), sha256=digest, role=role))
    return rows


def stage(output: Path, source_root: Path | None = None, tables_only: bool = False) -> dict:
    output = output.resolve()
    if source_root is not None:
        source_root = source_root.resolve()
    references = [] if tables_only else reference_rows(output, source_root)
    tables = []
    # Preflight everything before changing either manifests or derived copies.
    for table_id, relative, filename in TABLES:
        source = output / relative
        if not source.is_file():
            raise FileNotFoundError(f"Missing derived table source: {source}")
        tables.append(dict(table_id=table_id, source=str(source), file=str(output / "tables" / filename),
                           sha256=sha256(source), note=TABLE_NOTE))
    for row in references:
        source, target = Path(row["source"]), Path(row["snapshot"])
        if target.exists() and sha256(target) != row["sha256"]:
            raise RuntimeError(f"Refusing to overwrite differing frozen reference bytes: {target}")
        if not source.is_file():
            raise FileNotFoundError(f"Missing authoritative reference origin: {source}")
        if sha256(source) != row["sha256"]:
            raise RuntimeError(f"Authoritative reference origin differs from its frozen SHA256: {source}")
    table_copies = sum(atomic_copy(Path(row["source"]), Path(row["file"]), row["sha256"]) for row in tables)
    reference_copies = sum(atomic_copy(Path(row["source"]), Path(row["snapshot"]), row["sha256"]) for row in references)
    write_manifest(output / "tables/table_manifest.csv", ["table_id", "source", "file", "sha256", "note"], tables)
    if not tables_only:
        write_manifest(output / "inputs/preprocessing_reference/provenance_manifest.csv",
                       ["source", "snapshot", "sha256", "role"], references)
    return dict(tables_verified=len(tables), derived_tables_refreshed=table_copies,
                references_verified=len(references), references_copied=reference_copies)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, default=None,
                        help="Optional snRNAseq_seurat_objects folder containing james_scripts_20261104.")
    parser.add_argument("--tables-only", action="store_true",
                        help="Stage public scientific tables without accessing the private manuscript or its provenance.")
    args = parser.parse_args()
    result = stage(args.output, args.source_root, args.tables_only)
    print("; ".join(f"{key}={value}" for key, value in result.items()))


if __name__ == "__main__":
    main()
