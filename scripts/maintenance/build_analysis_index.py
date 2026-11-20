"""Build the result navigation from the maintained current panel manifest."""
import csv
from pathlib import Path
from urllib.parse import quote

from repository_paths import ROOT


def main():
    manifest = ROOT / "docs/provenance/current_panel_manifest.csv"
    with manifest.open(newline="", encoding="utf-8") as stream:
        panels = list(csv.DictReader(stream))
    results = ROOT / "results"
    lines = [
        "# Current panels and source data", "",
        "The panel labels follow the Results layout. Saved Figure 6 filenames "
        "retain the workflow numbering. See the [Results sections](README.md) "
        "for additional panels and supplementary tables.", "",
        "| Panel | Workflow | Panel files | Source data |",
        "| --- | --- | --- | --- |",
    ]

    def links(names):
        links = []
        for name in names.split("; "):
            path = ROOT / name
            if not path.is_file() or not path.is_relative_to(results / "current"):
                raise FileNotFoundError(f"Missing published panel artifact: {name}")
            links.append(f"[{path.name}]({quote(path.relative_to(results).as_posix())})")
        return "<br>".join(links)

    for row in panels:
        lines.append(f"| {row['panel']} | {row['workflow']} | {links(row['files'])} | {links(row['source_tables'])} |")
    lines += ["", "Figure 6D shows selected association examples. The class "
              "and network panels are exploratory.", "",
              "[Panel manifest](../docs/provenance/current_panel_manifest.csv) · "
              "[Provenance](../docs/REPOSITORY_ORGANIZATION.md)"]
    (results / "current_panels.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    from build_workflow_documentation import build
    summary = build()
    print(f"Portable indices built: {len(panels)} panel rows and {summary['workflows']} workflows")


if __name__ == "__main__":
    main()
