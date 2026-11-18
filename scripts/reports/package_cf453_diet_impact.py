"""Package exactly four direct-diet figure panels, source data and legends."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_path, resolve_manifest_path
import zipfile


def package(output):
    output = Path(output).resolve()
    rows = list(csv.DictReader((output / "panel_manifest.csv").open()))
    assert len(rows) == 12 and len({r["panel"] for r in rows}) == 4
    payload = {}
    for row in rows:
        for key, folder in (("file", "panels"), ("source_table", "source_data")):
            path = resolve_manifest_path(row[key], output)
            relative = f"{folder}/{path.name}"
            payload[relative] = path.read_bytes()
            row[key] = relative
    for name in ("README.md", "METHODS.md", "PANEL_LEGENDS.md"):
        payload[name] = (output / name).read_bytes()
    for name in ("primary_contrasts", "conditional_contrasts", "example_selection", "heatmap_selection"):
        path = output / "tables" / f"{name}.csv"
        payload[f"selection_and_contrasts/{path.name}"] = path.read_bytes()
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    payload["panel_manifest.csv"] = buffer.getvalue().encode()
    if (output / "completion_manifest.json").exists():
        payload["verification_summary.json"] = (output / "completion_manifest.json").read_bytes()
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(("path", "sha256", "bytes"))
    for name, contents in sorted(payload.items()):
        writer.writerow((name, hashlib.sha256(contents).hexdigest(), len(contents)))
    payload["SHA256SUMS.csv"] = buffer.getvalue().encode()
    destination = output / "CF453_Figure6_direct_diet_panels.zip"
    temporary = destination.with_suffix(".tmp.zip")
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, contents in sorted(payload.items()):
            archive.writestr(name, contents)
    with zipfile.ZipFile(temporary) as archive:
        assert archive.testzip() is None
        for name, contents in payload.items():
            assert hashlib.sha256(archive.read(name)).digest() == hashlib.sha256(contents).digest()
    temporary.replace(destination)
    summary = dict(path=str(destination), panel_groups=4, graphic_files=12,
                   source_tables=4, bytes=destination.stat().st_size, verification="CRC and SHA256 passed")
    (output / "bundle_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Explicit full Workflow15 output; current Figure6 uses package_cf453_diet_rewiring.py")
    arguments = parser.parse_args()
    package(arguments.output)
