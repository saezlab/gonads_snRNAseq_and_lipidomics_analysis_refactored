"""Recorded, lossless repository organization; never deletes analysis artifacts.

This dated migration is intentionally explicit. --plan only inventories; --apply
also moves the listed units and stages byte-identical public current inputs.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
PROVENANCE = ROOT / "docs/provenance"
MOVES = [
    ("results/zonation", "miscellaneous/older_results/zonation_two_marker", "Historical two-marker FCG reference; WF12/13 and old/new replay"),
    *[(f"results/{n}", f"miscellaneous/older_results/{n}", "Superseded Figure6 analysis/presentation; retained comparison baseline") for n in
      ("cf453_standalone_figure6", "cf453_diet_impact_figure6", "cf453_diet_rewiring_figure6")],
    *[(f"results/{n}", f"miscellaneous/comparisons/{n}", "Frozen cross-cohort comparison; does not define current standalone Figure6") for n in
      ("wdsd_network_comparison", "wd_crosscohort_robustness")],
    *[(f"data/{n}", f"miscellaneous/older_inputs/two_marker_zonation/{n}", "Historical zone definition; required for exact old-reference replay") for n in
      ("RNA/fcg_zonation_pseudobulks.rds", "RNA/gdx_zonation_pseudobulks.rds", "RNA/wdsd/wdsd_zonation_pseudobulks.rds",
       "RNA/scaled_RNA_zonation.csv", "RNA/meta_data_df_zonation.csv")],
    ("data/RNA/legacy", "miscellaneous/older_inputs/legacy_rna", "Historical RNA/culture inputs; optional KO expression audit retains an explicit dependency"),
    *[(f"results/mofa/RNA/{n}", f"miscellaneous/older_results/mofa/{n}", "Superseded MOFA model generation; current unversioned model set retained") for n in
      ("old_20241023", "old_20260313")],
    ("results/multi-omic/legacy", "miscellaneous/older_results/multi_omic", "Superseded gene-lipid structural analysis"),
    ("refactoring_docs", "miscellaneous/refactoring", "Historical refactor notes, debug scripts and terminal logs"),
    ("scripts/legacy", "miscellaneous/older_code/legacy", "Original notebook-style analyses; not current executable entries"),
    ("scripts/old", "miscellaneous/older_code/old", "Original code retained for provenance"),
    ("unused", "miscellaneous/older_code/unused", "Previously unused scripts"),
    ("output", "miscellaneous/quick_tests/output", "Ad hoc generated outputs"),
    ("tmp", "miscellaneous/quick_tests/tmp", "Temporary exploratory exports"),
    ("omnipathr-log", "miscellaneous/logs/omnipathr", "Historical package logs"),
    ("results/comment_responses", "miscellaneous/comparisons/comment_responses", "Earlier reviewer-response explorations"),
    *[(f"results/validation/{n}", f"miscellaneous/older_results/validation_zonation_two_marker/{n}", "Superseded zonation validation; current outputs live in approved release/fcg07") for n in
      ("mtmr11_zonation_validation_scatter.pdf", "sult3a1_zonation_validation_scatter.pdf", "zonation_network_validation_scores.csv",
       "zonation_network_validation_scores_jpc.csv", "mtmr11_top_lipid_joint_kd_correlation.csv")],
    *[(f"results/{n}", f"miscellaneous/older_results/collaborator_exports/{n}", "Historical collaborator export; active WF02 results retained") for n in
      ("DGE_hepatocytes_results.xlsx", "Doublets_and_Sex_DEG_results_all_hepatocytes.xlsx")],
    *[(n, "miscellaneous/local/" + n.replace("/", "_"), "Local IDE/history state; excluded from Git") for n in
      (".claude/settings.local.json", ".Rhistory", "results/lipidomic/.Rapp.history", "results/pseudobulk/.Rapp.history")],
    ("Rplots.pdf", "miscellaneous/quick_tests/Rplots.pdf", "Unattributed default graphics-device export"),
    ("tmp_figure_legend.txt", "miscellaneous/quick_tests/tmp_figure_legend.txt", "Earlier scratch writing"),
]
SNAPSHOTS = ["README.md", "PIPELINE.md", "AGENTS.md", "CLAUDE.md", "docs/AGENT_GUIDE.md",
             "docs/workflow_plans/_meta_cross_workflow_overview.md"]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def files(path):
    return sorted(p for p in path.rglob("*") if p.is_file()) if path.is_dir() else [path]


def relative(path):
    return path.relative_to(ROOT).as_posix()


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def organize(apply=False):
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    ledger = PROVENANCE / "repository_migration_20261116.csv"
    if ledger.exists() and apply:
        raise RuntimeError("Dated migration already applied; use the manifest verifier, not a second move")
    moves = list(MOVES)
    moves += [(relative(p), "miscellaneous/quick_tests/code/" + p.name, "Ad hoc helper outside the numbered pipeline")
              for p in sorted((ROOT / "scripts").glob("quick_*.R")) if p.name != "quick_lookup.R"]
    moves += [(relative(p), "miscellaneous/quick_tests/logs/" + p.name, "Generated historical test output")
              for p in sorted((ROOT / "tests").glob("test_output*.txt"))]
    moves += [(relative(p), "miscellaneous/quick_tests/logs/" + p.name, "Historical test report; current tests remain under tests/")
              for p in sorted((ROOT / "tests").glob("TEST_RESULTS*.md"))]
    moves += [(relative(p), "miscellaneous/logs/workflows/" + p.name, "Historical workflow terminal log")
              for p in sorted((ROOT / "results").glob("wf*.log"))]
    moves = [(s, d, why) for s, d, why in moves if (ROOT / s).exists()]
    for s, d, _ in moves:
        if (ROOT / d).exists():
            raise FileExistsError(d)
    records = []
    for s, d, why in moves:
        source, target = ROOT / s, ROOT / d
        for p in files(source):
            dest = target / p.relative_to(source) if source.is_dir() else target
            records.append(dict(action="move", old_path=relative(p), new_path=relative(dest),
                                sha256=sha(p), bytes=p.stat().st_size, role=why))
    for name in SNAPSHOTS:
        p = ROOT / name
        if p.exists():
            dest = ROOT / "miscellaneous/documentation/pre_cleanup_20261116" / name
            if dest.exists(): raise FileExistsError(dest)
            records.append(dict(action="snapshot", old_path=name, new_path=relative(dest),
                                sha256=sha(p), bytes=p.stat().st_size, role="Documentation before default-release promotion"))
    if not apply:
        write_csv(PROVENANCE / "repository_migration_plan_20261116.csv", records)
        print(json.dumps(dict(planned_files=len(records), move_units=len(moves), bytes=sum(x["bytes"] for x in records)), indent=2))
        return
    # Durable inventory exists before the first mutation.
    write_csv(ledger, records)
    for s, d, _ in moves:
        dest = ROOT / d; dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(ROOT / s), str(dest))
    for row in records:
        if row["action"] == "snapshot":
            dest = ROOT / row["new_path"]; dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / row["old_path"], dest)
        assert sha(ROOT / row["new_path"]) == row["sha256"], row["new_path"]
    registry_path = ROOT / "scripts/config/analysis_release.json"
    registry = json.loads(registry_path.read_text())
    previous = {r["from"] for r in registry["relocations"]}
    for s, d, _ in moves:
        if s not in previous:
            registry["relocations"].append({"from": s, "to": d})
    registry_path.write_text(json.dumps(registry, indent=2) + "\n")
    print(json.dumps(dict(verified_files=len(records), moved_units=len(moves), preserved_bytes=sum(x["bytes"] for x in records)), indent=2))


def stage_current_inputs():
    source = ROOT / "results/zonation_marker_update_20261116/inputs"
    target = ROOT / "data/RNA/zonation/current"
    pairs = [(source / c / f"{c}_zonation_pseudobulks_score.rds", target / f"{c}_zonation_pseudobulks_score.rds") for c in ("fcg", "gdx", "wdsd")]
    pairs += [(p, target / c / p.name) for c in ("fcg", "gdx") for p in sorted((source / c).glob("*.csv"))]
    pairs += [(p, target / "james_scripts_20261104" / p.name) for p in sorted((source / "preprocessing_reference").glob("*.R"))]
    pairs += [(source / "preprocessing_reference/Zonation_Genes.csv", target / "Zonation_Genes.csv")]
    records = []
    for s, d in pairs:
        digest = sha(s)
        if d.exists() and sha(d) != digest: raise RuntimeError(f"Refusing differing current input: {d}")
        if not d.exists():
            d.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(s, d)
        assert sha(d) == digest
        records.append(dict(path=relative(d), sha256=digest, bytes=d.stat().st_size,
                            approved_source=relative(s), release_id="expanded_markers_20261116"))
    write_csv(PROVENANCE / "current_input_manifest_20261116.csv", records)
    print(f"Local source copies: {len(records)} byte-identical files staged; serialized objects remain excluded from Git")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group()
    modes.add_argument("--apply", action="store_true")
    modes.add_argument("--plan", action="store_true")
    modes.add_argument("--stage-current-inputs", action="store_true")
    a = p.parse_args()
    if a.stage_current_inputs: stage_current_inputs()
    else: organize(a.apply)
