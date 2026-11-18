"""Build portable public result navigation from the approved panel manifest."""
import csv
import os
from pathlib import Path
from repository_paths import ROOT, analysis_path, resolve_recorded_path


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    current = analysis_path("current_root")
    panels = list(csv.DictReader((current / "report/panel_trace.csv").open()))
    portable = []
    for row in panels:
        item = {k: row[k] for k in ("panel", "status", "workflow")}
        for before, after in (("generated_files", "files"), ("source_tables", "source_tables")):
            paths = [resolve_recorded_path(p) for p in row[before].split("; ") if p]
            assert all(p.is_file() for p in paths), item["panel"]
            item[after] = "; ".join(p.relative_to(ROOT).as_posix() for p in paths)
        item["notes"] = row["change_to_review"]
        portable.append(item)
    write_csv(ROOT / "docs/provenance/current_panel_manifest.csv", portable)
    lines = ["# Current panels and source data", "", "Default release: `expanded_markers_20261116`. Panels remain separate for manual assembly. This navigation uses repository-relative links; original production manifests are retained unchanged.", "",
             "| Panel | Status | Panel files | Source data |", "| --- | --- | --- | --- |"]
    def links(names):
        return "<br>".join(f"[{Path(n).name}]({os.path.relpath(ROOT / n, ROOT / 'results')})" for n in names.split("; ") if n)
    for row in portable:
        lines.append(f"| {row['panel']} | {row['status']} | {links(row['files'])} | {links(row['source_tables'])} |")
    lines += ["", "Figure4A/S4B–C and their nucleus-level statistics are deferred to James. Figure5C exports include both the faithful historical lipid-sex encoding and the explicit gene-context alternative; use the corresponding legend. The prespecified male-periportal Figure5F is the manuscript source. Figure6 class/edge associations remain exploratory.", "",
              "See [the pipeline](../PIPELINE.md), [portable panel manifest](../docs/provenance/current_panel_manifest.csv), [replacement ledger](zonation_marker_update_20261116/report/new_replacements.csv) and [provenance map](../docs/REPOSITORY_ORGANIZATION.md)."]
    (ROOT / "results/current_panels.md").write_text("\n".join(lines) + "\n")
    # Keep the full workflow catalog synchronized with the code-audited contracts.
    from build_workflow_documentation import build
    summary = build()
    print(f"Portable indices built: {len(portable)} panel rows; {summary['workflows']} workflows")


if __name__ == "__main__":
    main()
