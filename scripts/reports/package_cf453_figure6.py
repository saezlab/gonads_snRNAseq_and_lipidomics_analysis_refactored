"""Package independent scientific panels and their source data, without assembly."""
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
    payload = {}
    for row in rows:
        source = resolve_manifest_path(row["file"], output)
        relative = "panels/" + source.name
        payload[relative] = source.read_bytes()
        row["file"] = relative
        source_table = resolve_manifest_path(row["source_table"], output)
        relative_table = "source_data/" + source_table.name
        payload[relative_table] = source_table.read_bytes()
        row["source_table"] = relative_table
    for name in ("METHODS.md", "PANEL_LEGENDS.md"):
        payload[name] = (output / name).read_bytes()
    for name in ("selection_ledger", "representatives", "representative_delta_summary", "class_pair_density"):
        path = output / "tables" / f"{name}.csv"
        if path.exists():
            payload[f"selection_and_summaries/{path.name}"] = path.read_bytes()
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    payload["panel_manifest.csv"] = buffer.getvalue().encode()
    text = """CF453 standalone Figure 6: independent panel files

No assembled manuscript figure is included. Each panel is supplied as editable
vector PDF/SVG and a 300-dpi PNG preview. Matching source data are in source_data/.
The portable panel_manifest.csv maps every graphic to its source table and size.
PANEL_LEGENDS.md and METHODS.md explain selection, statistics and interpretation.
selection_and_summaries/ includes representative selection and descriptive data.

Interpretation for manuscript assembly:
- RNA interactions use the prespecified multiple-testing correction.
- Individual gene-lipid correlations do not pass global FDR <=0.05; network,
  density and selected-example panels are exploratory at seven mice per group.
- Fig6C_SD/WD retain the saturated any-edge lipid membership result. The
  Fig6C_pair_density_SD/WD files provide a complementary descriptive view using
  the same candidate gene-lipid universe in both diets.
- Lipid measurements are whole-liver abundance, paired to zonated hepatocyte RNA.
- E/F replace the original perturbation panels with within-CF453 analyses and
  must not be presented as knockdown validation.

Full models, all coefficients, bootstrap indices and independent numerical and
visual QA remain in the explicitly selected Workflow14 output directory.
SHA256SUMS.csv verifies every other file in this bundle.
"""
    payload["README.txt"] = text.encode()
    checksums = io.StringIO()
    writer = csv.writer(checksums)
    writer.writerow(("path", "sha256", "bytes"))
    for name, contents in sorted(payload.items()):
        writer.writerow((name, hashlib.sha256(contents).hexdigest(), len(contents)))
    payload["SHA256SUMS.csv"] = checksums.getvalue().encode()
    destination = output / "CF453_Figure6_panels.zip"
    temporary = destination.with_suffix(".tmp.zip")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, contents in sorted(payload.items()):
            archive.writestr(name, contents)
    with zipfile.ZipFile(temporary) as archive:
        assert archive.testzip() is None
        assert len(archive.namelist()) == len(payload)
        for name, contents in payload.items():
            assert hashlib.sha256(archive.read(name)).digest() == hashlib.sha256(contents).digest()
    temporary.replace(destination)
    summary = dict(path=str(destination), panel_groups=len({r["panel"] for r in rows}),
                   graphic_files=len(rows), files=len(payload), bytes=destination.stat().st_size,
                   crc_and_sha256="passed")
    (output / "bundle_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Explicit full Workflow14 output; current Figure6 uses package_cf453_diet_rewiring.py")
    arguments = parser.parse_args()
    package(arguments.output)
