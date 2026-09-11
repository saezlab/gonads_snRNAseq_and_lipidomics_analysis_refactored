#!/usr/bin/env python3
"""Summarize a Factor 1 MOON parameter grid without fitting or changing solutions.

Ranking is a declared threshold preference, not model selection or a claim that
different scoring depths/compression configurations have equivalent biology.
"""
import argparse
import collections
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess


REQUIRED = ["run_id", "run_dir", "compression", "n_layers", "topology_radius",
            "primary_threshold", "secondary_threshold", "status", "converged", "valid",
            "nodes", "edges", "upstream_tfs", "downstream_nodes", "solution_stem"]
OPTIONAL = ["phase", "realized_max_level", "candidate_edges", "consistent_edges",
            "short_route_pairs", "route_pairs", "run_manifest"]


def boolean(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def compact(value):
    return "" if value is None else format(value, ".12g")


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def stratum(row):
    return row["compression"], row["n_layers"], row["topology_radius"]


def rank_key(row):
    return (-row["secondary_threshold"], -row["primary_threshold"],
            -row["downstream_nodes"], -row["upstream_tfs"], row["edges"], row["candidate_id"])


def mark_frontier(rows):
    """Pareto frontier within one configuration; never pool scoring configurations."""
    frontier = []
    for row in sorted(rows, key=rank_key):
        objective = (row["secondary_threshold"], row["primary_threshold"],
                     row["downstream_nodes"], row["upstream_tfs"], -row["edges"])
        dominated = any(all(left >= right for left, right in zip(previous, objective)) and previous != objective for previous in frontier)
        row["pareto_within_stratum"] = not dominated
        if not dominated and objective not in frontier:
            frontier.append(objective)


def resolve_path(value, grid_path):
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    for candidate in [Path.cwd() / path, grid_path.parent / path]:
        if candidate.exists():
            return candidate.resolve()
    return ((Path.cwd() if path.parts and path.parts[0] == "results" else grid_path.parent) / path).resolve()


def markdown_table(rows):
    lines = ["| Candidate | Compression | Depth / radius | Secondary / primary | TFs / lipids | Edges |",
             "|---|---|---:|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['candidate_id']} | {row['compression']} | {compact(row['n_layers'])} / {compact(row['topology_radius'])} | {compact(row['secondary_threshold'])} / {compact(row['primary_threshold'])} | {compact(row['upstream_tfs'])} / {compact(row['downstream_nodes'])} | {compact(row['edges'])} |")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, default=Path("results/cosmos_factor1_parameter_search/parameter_grid.csv"))
    parser.add_argument("--output", type=Path, help="Defaults to GRID_PARENT/summary; never writes into solution directories")
    parser.add_argument("--pkn", type=Path, default=Path("results/cosmos_mouse_pkn"))
    parser.add_argument("--max-edges", type=int, default=500)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--write-template", type=Path, help="Write an empty grid header and schema instructions, then exit")
    args = parser.parse_args()
    if args.write_template:
        args.write_template.mkdir(parents=True, exist_ok=True)
        template = args.write_template / "parameter_grid_template.csv"
        if template.exists():
            raise FileExistsError(template)
        write_csv(template, [], REQUIRED + OPTIONAL)
        (args.write_template / "parameter_grid_schema.md").write_text(
            "# Factor 1 MOON search summary inputs\n\n" +
            "Required columns: `" + "`, `".join(REQUIRED) + "`.\n\n" +
            "Optional columns: `" + "`, `".join(OPTIONAL) + "`.\n\n" +
            "A row describes one score configuration and threshold pair. `valid` certifies the scoring/extraction pipeline's provenance, sign and connectivity checks; `converged` certifies its own coherence loop. Only `status=complete` passes. Missing certification does not pass. `run_dir` and `solution_stem` locate eventual `STEM_SIF.csv` and `STEM_ATT.csv` files. Solutions may be materialized after selection; the summary records that availability separately.\n\n" +
            "Rank secondary cutoff first, primary cutoff second, then lipid and TF coverage, then fewer edges. Require nonempty valid/converged solutions with at least one TF and one lipid and <=500 edges. A separate interpretable subset requires at least two TFs and two lipids. Pareto frontiers are computed within each compression/depth/topology-radius configuration. Cross-configuration counts describe eligibility reproducibility, not biological network stability.\n")
        print(template)
        return
    if args.max_edges < 1:
        raise ValueError("--max-edges must be positive")
    grid_path = args.grid.resolve()
    with grid_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = set(REQUIRED) - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Missing required grid columns: " + ", ".join(sorted(missing)))
        original_rows = list(reader)
    output = args.output.resolve() if args.output else grid_path.parent / "summary"
    if any(output == resolve_path(row["run_dir"], grid_path) for row in original_rows if row.get("run_dir")):
        raise ValueError("Summary output must not be a solution run directory")
    output.mkdir(parents=True, exist_ok=True)
    normalized, duplicate_count = {}, 0
    for index, raw in enumerate(original_rows, start=2):
        row = {column: raw.get(column, "") for column in REQUIRED + OPTIONAL}
        row["grid_line"] = index
        row["compression"] = str(row["compression"]).strip()
        for column in ["n_layers", "topology_radius", "primary_threshold", "secondary_threshold", "nodes", "edges", "upstream_tfs", "downstream_nodes", "realized_max_level", "candidate_edges", "consistent_edges", "short_route_pairs", "route_pairs"]:
            row[column] = number(row[column])
        row["converged"], row["valid"] = boolean(row["converged"]), boolean(row["valid"])
        identity = (row["run_id"], row["compression"], row["n_layers"], row["topology_radius"], row["primary_threshold"], row["secondary_threshold"], row["solution_stem"])
        row["candidate_id"] = str(row["run_id"]) + "__" + hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:10]
        if identity in normalized:
            duplicate_count += 1
        # An appended retry/refinement supersedes the earlier row for the same
        # configuration, thresholds and artifact identity; all raw rows remain intact.
        normalized[identity] = row
    rows = list(normalized.values())
    count_columns = ["nodes", "edges", "upstream_tfs", "downstream_nodes"]
    for row in rows:
        reasons = []
        if str(row["status"]).strip().lower() != "complete": reasons.append("status_not_complete")
        if not row["converged"]: reasons.append("not_converged")
        if not row["valid"]: reasons.append("validity_not_certified")
        if any(row[column] is None or row[column] < 0 or row[column] != int(row[column]) for column in count_columns): reasons.append("invalid_network_counts")
        if row["nodes"] is None or row["nodes"] <= 0 or row["edges"] is None or row["edges"] <= 0: reasons.append("empty_network")
        if row["edges"] is not None and row["edges"] > args.max_edges: reasons.append("edge_cap_exceeded")
        if row["upstream_tfs"] is None or row["upstream_tfs"] < 1: reasons.append("no_upstream_tf")
        if row["downstream_nodes"] is None or row["downstream_nodes"] < 1: reasons.append("no_downstream_lipid")
        if row["primary_threshold"] is None or row["secondary_threshold"] is None or row["primary_threshold"] < row["secondary_threshold"] or row["secondary_threshold"] < 0: reasons.append("invalid_threshold_order")
        if not row["compression"] or row["n_layers"] is None or row["n_layers"] < 1 or row["topology_radius"] is None or row["topology_radius"] < 1: reasons.append("invalid_configuration")
        row["eligible"] = not reasons
        row["eligibility_reasons"] = "|".join(reasons)
        row["interpretable_eligible"] = row["eligible"] and row["upstream_tfs"] >= 2 and row["downstream_nodes"] >= 2
        row.update(priority_rank="", within_stratum_rank="", interpretable_priority_rank="", interpretable_within_stratum_rank="", pareto_within_stratum=False)
        directory = resolve_path(row["run_dir"], grid_path)
        row["resolved_run_dir"] = str(directory)
        row["solution_materialized"] = bool(row["solution_stem"]) and (directory / (row["solution_stem"] + "_SIF.csv")).exists() and (directory / (row["solution_stem"] + "_ATT.csv")).exists()
    eligible = sorted([row for row in rows if row["eligible"]], key=rank_key)
    interpretable = sorted([row for row in rows if row["interpretable_eligible"]], key=rank_key)
    groups = collections.defaultdict(list)
    for rank, row in enumerate(eligible, start=1):
        row["priority_rank"] = rank
        groups[stratum(row)].append(row)
    for rank, row in enumerate(interpretable, start=1): row["interpretable_priority_rank"] = rank
    shortlist = {}
    for group in groups.values():
        for rank, row in enumerate(group, start=1): row["within_stratum_rank"] = rank
        readable = [row for row in group if row["interpretable_eligible"]]
        for rank, row in enumerate(readable, start=1): row["interpretable_within_stratum_rank"] = rank
        mark_frontier(group)
        for row, reason in [(group[0], "highest_threshold_priority_within_configuration")] + ([(readable[0], "highest_threshold_priority_with_at_least_two_TFs_and_two_lipids")] if readable else []):
            entry = shortlist.setdefault(row["candidate_id"], dict(row, shortlist_reasons=[]))
            entry["shortlist_reasons"].append(reason)
    short_rows = sorted(shortlist.values(), key=rank_key)
    for row in short_rows: row["shortlist_reasons"] = "|".join(row["shortlist_reasons"])
    frontier = sorted([row for row in eligible if row["pareto_within_stratum"]], key=rank_key)
    fields = ["candidate_id", "priority_rank", "within_stratum_rank", "interpretable_priority_rank", "interpretable_within_stratum_rank", "eligible", "interpretable_eligible", "pareto_within_stratum", "eligibility_reasons"] + REQUIRED + OPTIONAL + ["grid_line", "resolved_run_dir", "solution_materialized"]
    write_csv(output / "grid_assessment.csv", rows, fields)
    write_csv(output / "ranked_eligible.csv", eligible, fields)
    write_csv(output / "ranked_interpretable.csv", interpretable, fields)
    write_csv(output / "pareto_within_configuration.csv", frontier, fields)
    write_csv(output / "configuration_shortlist.csv", short_rows, fields + ["shortlist_reasons"])

    cutoff_groups = collections.defaultdict(list)
    for row in rows:
        if row["primary_threshold"] is not None and row["secondary_threshold"] is not None:
            cutoff_groups[row["primary_threshold"], row["secondary_threshold"]].append(row)
    stability = []
    for (primary, secondary), group in sorted(cutoff_groups.items(), key=lambda item: (-item[0][1], -item[0][0])):
        tested = {stratum(row) for row in group}
        passed = {stratum(row) for row in group if row["eligible"]}
        readable = {stratum(row) for row in group if row["interpretable_eligible"]}
        stability.append({"primary_threshold": primary, "secondary_threshold": secondary,
                          "tested_configurations": len(tested), "eligible_configurations": len(passed), "interpretable_configurations": len(readable),
                          "eligible_fraction_of_tested": len(passed) / len(tested), "interpretable_fraction_of_tested": len(readable) / len(tested),
                          "eligible_compressions": "|".join(sorted({row["compression"] for row in group if row["eligible"]})),
                          "eligible_depths": "|".join(compact(value) for value in sorted({row["n_layers"] for row in group if row["eligible"]})),
                          "interpretation": "Eligibility across tested configurations; not edge, endpoint or biological stability"})
    stability_fields = ["primary_threshold", "secondary_threshold", "tested_configurations", "eligible_configurations", "interpretable_configurations", "eligible_fraction_of_tested", "interpretable_fraction_of_tested", "eligible_compressions", "eligible_depths", "interpretation"]
    write_csv(output / "eligibility_stability_counts.csv", stability, stability_fields)
    exporter = Path(__file__).with_name("export_factor1_network.py").resolve()
    export_plan = []
    for row in short_rows:
        command = ["python3", str(exporter), "--results", row["resolved_run_dir"], "--stem", row["solution_stem"], "--pkn", str(args.pkn.resolve())]
        export_plan.append({"candidate_id": row["candidate_id"], "run_dir": row["resolved_run_dir"], "solution_stem": row["solution_stem"],
                            "ready_for_export": row["solution_materialized"], "command": shlex.join(command),
                            "graphml_path": str(Path(row["resolved_run_dir"]) / (row["solution_stem"] + ".graphml")),
                            "validation": "Before delivery, independently join exact PKN edge provenance and check signs, upstream reachability and downstream observability using the exporter; retain separate TF activities, Factor 1 weights and inferred scores."})
    write_csv(output / "candidate_export_plan.csv", export_plan, ["candidate_id", "run_dir", "solution_stem", "ready_for_export", "command", "graphml_path", "validation"])
    exclusions = collections.Counter(reason for row in rows for reason in row["eligibility_reasons"].split("|") if reason)
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "source_grid": str(grid_path),
                "source_sha256": hashlib.sha256(grid_path.read_bytes()).hexdigest(), "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "raw_rows": len(original_rows), "unique_rows": len(rows), "superseded_duplicate_rows": duplicate_count,
                "eligible_rows": len(eligible), "interpretable_rows": len(interpretable), "pareto_rows": len(frontier), "shortlist_rows": len(short_rows),
                "max_edges": args.max_edges, "exclusion_reason_counts": dict(exclusions),
                "ranking": "secondary descending, primary descending, downstream lipid endpoints descending, upstream TFs descending, edges ascending; candidate ID only breaks exact ties",
                "pareto_strata": ["compression", "n_layers", "topology_radius"], "pareto_objectives": ["maximize secondary", "maximize primary", "maximize lipid endpoints", "maximize TF endpoints", "minimize edges"],
                "grid_validity_scope": "The input valid flag certifies the producer's audits; this script does not recalculate networks or recertify unmaterialized artifacts.",
                "duplicate_policy": "Latest grid line supersedes an identical run/configuration/threshold/artifact identity.",
                "biological_stability": "Not established by eligibility counts. Materialize comparable candidate networks and compare exact endpoint/edge identities, signs and route support before claiming stability."}
    (output / "summary_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    report = (
        "# Factor 1 MOON parameter-search comparison\n\n"
        f"The grid contains **{len(rows):,} unique configurations/cutoff rows** ({len(original_rows):,} input rows). "
        f"**{len(eligible):,} rows** meet the declared validity, convergence, endpoint and ≤{args.max_edges}-edge criteria; "
        f"**{len(interpretable):,} rows** additionally retain at least two TFs and two lipid endpoints. "
        "These are candidate networks for comparison, not a single selected biological model.\n\n"
        "## Declared ranking and validity\n\n"
        "Threshold priority is lexicographic: higher secondary cutoff first, then higher primary cutoff, then more lipid endpoints, more upstream TFs and fewer edges. "
        "The interpretable ranking applies the same rule after requiring ≥2 TFs and ≥2 lipids. A complete, converged, valid, nonempty network with an upstream TF and downstream lipid is required. "
        "Missing or false validity certification fails eligibility. Materialization is reported separately because the scoring workflow may write candidate solutions after grid selection.\n\n"
        "Depth, compression and topology radius change the scoring problem. Priority ranks across those configurations reflect the requested threshold preference; they do not make the scores or mechanisms biologically equivalent. "
        "Pareto frontiers and first-ranked shortlists are therefore computed separately within each configuration.\n\n"
        "## Configuration-specific shortlist\n\n" +
        (markdown_table(short_rows[:30]) if short_rows else "No candidates meet the eligibility rules. Review the exclusion reasons and failed/convergence records before changing cutoffs.") +
        (f"\n\nShowing 30 of {len(short_rows)} rows; the full shortlist is in `configuration_shortlist.csv`." if len(short_rows) > 30 else "") +
        "\n\n## Pareto comparison and sensitivity\n\n"
        f"There are **{len(frontier):,} nondominated rows within their own configurations**. The objectives retain cutoff stringency, lipid/TF coverage and network size explicitly, without collapsing them into an arbitrary composite score. "
        "`eligibility_stability_counts.csv` counts how many tested depth/compression/radius configurations pass at each matched cutoff pair. Its denominators can differ when coarse and refined searches cover different cells. "
        "Those counts measure reproducibility of eligibility, not agreement in TF identities, lipid identities, edges, signs or mechanisms.\n\n"
        "## Candidate delivery and provenance\n\n"
        "`candidate_export_plan.csv` lists materialization status and the exact exporter command for each configuration-specific shortlist row. Commands are plans only; this summarizer never fits models or writes into solution directories. "
        "Once the corresponding SIF/ATT files exist, the exporter independently checks exact PKN membership and evidence, score-sign consistency, upstream controllability and downstream observability. "
        "It writes Cytoscape-ready GraphML with distinct Factor 1 weight, inferred MOON score and TF-activity fields, route tables and figures. Do not call an unmaterialized or unaudited candidate a validated network.\n\n"
        "Measured-species reaction rules remain substrate-specificity hypotheses. Positive reaction traversal is not a flux or stoichiometric production prediction. "
        "Scoring depth does not cap every TF-to-lipid route in a merged graph; compare bounded catalyst-only routes as well as unrestricted connectivity.\n\n"
        "## Files\n\n"
        "- `grid_assessment.csv`: every deduplicated row and explicit exclusion reasons.\n"
        "- `ranked_eligible.csv` and `ranked_interpretable.csv`: threshold-priority rankings, with configuration-specific ranks.\n"
        "- `pareto_within_configuration.csv` and `configuration_shortlist.csv`: within-configuration tradeoffs and shortlisted candidates.\n"
        "- `eligibility_stability_counts.csv`: matched-cutoff eligibility counts and tested denominators.\n"
        "- `candidate_export_plan.csv`: candidate artifact locations and auditable export commands.\n"
        "- `summary_manifest.json`: exact source/script hashes, criteria and row counts.\n"
    )
    (output / "parameter_search_report.md").write_text(report)
    print(json.dumps(manifest, indent=2))
    if not args.no_plot:
        subprocess.run(["Rscript", str(Path(__file__).with_name("plot_factor1_parameter_search.R")), str(output)], check=True)


if __name__ == "__main__":
    main()
