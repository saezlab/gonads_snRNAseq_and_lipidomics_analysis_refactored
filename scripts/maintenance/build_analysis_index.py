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
    catalog = []
    def add(key, state, entry, inputs, outputs, figures, dependency):
        catalog.append(dict(id=key, status=state, entry=entry, inputs=inputs, outputs=outputs, figures=figures, dependency_note=dependency))
    add("whole_cell_rna", "main", "scripts/workflows/02_fcg_pseudobulk_analysis.R", "01 RNA CSVs; proteomics/GDX validation inputs", "results/pseudobulk", "main RNA / independent validation", "Whole-hepatocyte model, distinct from zonation-model sex statistics")
    add("lipidomics", "main", "scripts/workflows/03_fcg_lipidomic_analysis.R", "cleaned lipidomics; overview metadata", "results/lipidomic", "main lipid panels", "Feeds zonation and structural workflows")
    add("mofa", "main", "scripts/workflows/04_mofa_integration.R", "whole-cell RNA; lipids; trained current model set", "results/mofa", "global networks", "Training skipped by normal analysis")
    add("expanded_zonation", "current", "scripts/workflows/17_expanded_marker_zonation_update.R", "data/RNA/zonation/current; unchanged lipid/KD inputs", "results/zonation_marker_update_20261116/fcg05; fcg07; gdx_zonation", "4B;5A-D/F;S5A;S10/S13/S14/S15", "Expanded marker default; original scientific artifacts unchanged")
    add("standalone_diet", "current", "scripts/workflows/16_cf453_diet_rewiring_figure6.R", "current14models; current15direct tests; CF453 lipids", "results/zonation_marker_update_20261116/cf453_rewiring; cf453_layouts", "6A-C; two6Aalternatives", "Current14/15 model stages are dependencies; associations exploratory")
    add("cyp2c38_in_vivo", "main independent validation", "scripts/workflows/11_cyp2c38_in_vivo_validation.R", "local KO lipids; whole-hepatocyte RNA; frozen Figure3 pairs", "results/validation/cyp2c38_in_vivo", "independent in-vivo validation", "Does not consume expanded-marker zonation")
    for key, entry in (("historical_wf12", "12_wdsd_network_comparison.R"), ("historical_wf13", "13_wd_crosscohort_robustness.R")):
        add(key, "optional historical comparison", "scripts/workflows/" + entry, "original two-marker inputs and frozen networks", analysis_path(key).relative_to(ROOT).as_posix(), "historical comparisons", "Do not rebase one side to current zones")
    add("causal_models", "active optional extension", "scripts/moon/README.md", "frozen factor inputs; PKN; selection references", "results/causal_modelling", "causal hypotheses and sensitivities", "See CAUSAL_ANALYSIS_MAP; pinned full-fit environment")
    add("source_preprocessing", "external/reference", "scripts/workflows/00_snrnaseq_preprocessing.R", "external raw sci-RNA alignment environment", "original Seurat objects", "upstream provenance", "Not rerun during cleanup")
    add("human_xenium", "optional external validation", "scripts/workflows/09_human_analysis.py; scripts/workflows/10_xenium_analysis", "external human/Xenium raw data", "results/human; results/xenium", "independent validation", "Requires external input restoration")
    write_csv(ROOT / "docs/analysis_catalog.csv", catalog)
    print(f"Portable indices built: {len(portable)} panel rows; {len(catalog)} analysis branches")


if __name__ == "__main__":
    main()
