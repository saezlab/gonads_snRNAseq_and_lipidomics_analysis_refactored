"""Package exactly three standalone zonation/diet panel groups and source data."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "maintenance"))
from repository_paths import analysis_current, resolve_manifest_path
import zipfile


def package(output):
    output = Path(output).resolve()
    rows = list(csv.DictReader((output / "panel_manifest.csv").open()))
    assert len(rows) == 9 and len({r["panel"] for r in rows}) == 3
    payload = {}
    for row in rows:
        for key, folder in (("file", "panels"), ("source_table", "source_data")):
            path = resolve_manifest_path(row[key], output)
            relative = f"{folder}/{path.name}"
            payload[relative] = path.read_bytes()
            row[key] = relative
    for name in ("README.md", "METHODS.md", "PANEL_LEGENDS.md"):
        payload[name] = (output / name).read_bytes()
    for name in ("programmes", "network_genes", "sample_order", "threshold_sensitivity", "broad_sensitivity"):
        path = output / "tables" / f"{name}.csv"
        payload[f"selection_and_sensitivity/{path.name}"] = path.read_bytes()
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
    destination = output / "CF453_Figure6_zonation_diet_panels.zip"
    temporary = destination.with_suffix(".tmp.zip")
    with zipfile.ZipFile(temporary, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, contents in sorted(payload.items()):
            archive.writestr(name, contents)
    with zipfile.ZipFile(temporary) as archive:
        assert archive.testzip() is None
        for name, contents in payload.items():
            assert hashlib.sha256(archive.read(name)).digest() == hashlib.sha256(contents).digest()
    temporary.replace(destination)
    summary = dict(path=str(destination), panel_groups=3, graphic_files=9, source_tables=3,
                   bytes=destination.stat().st_size, verification="CRC and SHA256 passed")
    (output / "bundle_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(analysis_current("cf453_rewiring")))
    package(parser.parse_args().output)
