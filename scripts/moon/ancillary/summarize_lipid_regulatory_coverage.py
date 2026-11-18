#!/usr/bin/env python3
"""Consolidate per-measurement regulatory coverage annotations without remapping.

The frozen crosswalk and supplementary specificity audit remain separate evidence
layers. Ontology-descendant/compatible candidates are not accepted identities.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/causal_modelling"
GROUPS = {
    ("GuideToPharmacology", "curated_mouse_current"): "GtoPdb_primary_curated_signed",
    ("GuideToPharmacology", "curated_mouse_current_secondary_chemical_mapping"): "GtoPdb_secondary_mapping_signed",
    ("BRENDA", "cached_brenda_fallback"): "BRENDA_cached_signed",
    ("STITCH", "cached_stitch_fallback"): "STITCH_cached_signed_orientation_uncertain",
    ("MRCLinksDB", "excluded_unsigned_association"): "MRCLinksDB_unsigned_association",
}
UNSIGNED = "MRCLinksDB_unsigned_association"


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def truth(value):
    return str(value).lower() in {"true", "1", "yes"}


def joined(values):
    return "|".join(sorted(values))


def group(row, raw_mrc=False):
    key = (row["resource"], "excluded_unsigned_association" if raw_mrc else row["evidence_tier"])
    if key not in GROUPS:
        raise ValueError(f"Unrecognized resource/evidence tier: {key}")
    return GROUPS[key]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BASE / "lipid_regulatory_coverage_overview")
    parser.add_argument("--frozen-audit", type=Path, default=BASE / "lipid_regulatory_mapping_audit")
    parser.add_argument("--resolution-audit", type=Path, default=BASE / "lipid_regulatory_resolution_audit")
    args = parser.parse_args()
    out, frozen_dir, resolution_dir = (p.resolve() for p in (args.output, args.frozen_audit, args.resolution_audit))
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error(f"Output directory is occupied: {out}")
    inventory_path = frozen_dir / "measurement_inventory.csv"
    frozen_path = frozen_dir / "feature_edge_matches.csv"
    resolution_path = resolution_dir / "specific_regulatory_candidates.csv"
    dependencies = [inventory_path, frozen_path, resolution_path,
                    frozen_dir / "summary.json", resolution_dir / "summary.json", Path(__file__).resolve()]
    for path in dependencies:
        if not path.is_file():
            parser.error(f"Missing detailed audit dependency: {path}")
    inventory = read_csv(inventory_path)
    features = {r["feature_id"] for r in inventory}
    assert len(inventory) == len(features) == 1254
    frozen_rows = read_csv(frozen_path)
    resolution_rows = read_csv(resolution_path)
    # Do not count the duplicate raw-cache copy of an accepted BRENDA/STITCH row.
    frozen_rows = [r for r in frozen_rows if r["collection"] == "retained_regulation_export" or
                   (r["collection"] == "raw_mouse_service_cache" and r["resource"] == "MRCLinksDB")]
    state = {name: collections.defaultdict(set) for name in features}
    current_exact = set()
    incompatible = collections.Counter()
    evidence_groups = set(GROUPS.values())
    for row in frozen_rows:
        feature = row["feature"]
        assert feature in state
        eg = group(row, row["collection"] == "raw_mouse_service_cache")
        category = row["coverage_category"]
        data = state[feature]
        if category.startswith("exact_"):
            data["exact"].add(eg)
            data[category].add(eg)
            if truth(row["in_current_curated_pkn"]) and eg == "GtoPdb_primary_curated_signed":
                current_exact.add(feature)
        elif category == "compatible_resolution":
            data["compatible"].add(eg)
            data["frozen_compatible"].add(eg)
            if truth(row["in_current_curated_pkn"]):
                data["compatible_current_curated_edges"].add(eg)
        elif category == "class_or_ancestor":
            data["broad_class"].add(eg)
        else:
            raise ValueError(f"Unknown frozen coverage category: {category}")
        data["chemical_ids"].add(row["chemical_id"])
    for row in resolution_rows:
        feature = row["feature"]
        assert feature in state
        if row["formula_check"] == "formula_incompatible":
            incompatible[feature] += 1
            continue
        if truth(row["identity_confirmed"]):
            raise ValueError("Resolution rows must remain candidates, not confirmed identities")
        eg = group(row)
        data = state[feature]
        data["compatible"].add(eg)
        data["resolution_compatible"].add(eg)
        data["chemical_ids"].add(row["regulatory_chemical"])
        data["resolution_candidate_keys"].add((eg, row["regulatory_chemical"]))
        if truth(row["absent_from_frozen_candidate_table"]):
            data["new_resolution"].add(eg)
            data["new_candidate_keys"].add((eg, row["regulatory_chemical"]))
        if int(row["retained_curated_signed_edges"]) > 0:
            data["compatible_current_curated_edges"].add(eg)
        if int(row["retained_exploratory_signed_edges"]) > 0:
            data["compatible_exploratory_edges"].add(eg)
    rows = []
    for feature in inventory:
        name, data = feature["feature_id"], state[feature["feature_id"]]
        specific = data["exact"] | data["compatible"]
        signed_specific = specific - {UNSIGNED}
        rows.append({
            "feature": name, "lipid_class": feature["class"],
            "accepted_frozen_input": truth(feature["accepted_frozen_input"]),
            "in_historical_factor1_64": truth(feature["in_factor1"]),
            "two_lpc_candidate": name in {"LPC(16:0/0:0)", "LPC(18:0/0:0)"},
            "frozen_node_id": feature["frozen_node_id"],
            "frozen_mapping_exclusion": feature["frozen_exclusion_reason"],
            "oriented_factor1_weight": feature["oriented_weight"],
            "frozen_exact_chebi_ids": feature["exact_chebi_ids"],
            "frozen_exact_signed_evidence_groups": joined(data["exact"] - {UNSIGNED}),
            "frozen_exact_unsigned_evidence_groups": joined(data["exact"] & {UNSIGNED}),
            "frozen_exact_same_identifier_groups": joined(data["exact_same_identifier"]),
            "frozen_exact_explicit_alias_groups": joined(data["exact_explicit_alias"]),
            "frozen_exact_conjugate_acid_base_groups": joined(data["exact_conjugate_acid_base"]),
            "frozen_exact_current_curated_pkn_match": name in current_exact,
            "compatible_signed_candidate_groups": joined(data["compatible"] - {UNSIGNED}),
            "compatible_unsigned_candidate_groups": joined(data["compatible"] & {UNSIGNED}),
            "compatible_candidate_groups_union": joined(data["compatible"]),
            "compatible_frozen_candidate_groups": joined(data["frozen_compatible"]),
            "compatible_resolution_audit_groups": joined(data["resolution_compatible"]),
            "broad_class_signed_groups": joined(data["broad_class"] - {UNSIGNED}),
            "broad_class_unsigned_groups": joined(data["broad_class"] & {UNSIGNED}),
            "new_resolution_candidate": bool(data["new_resolution"]),
            "new_resolution_signed_groups": joined(data["new_resolution"] - {UNSIGNED}),
            "new_resolution_unsigned_groups": joined(data["new_resolution"] & {UNSIGNED}),
            "new_chemical_resource_candidate_pairs": len(data["new_candidate_keys"]),
            "compatible_candidates_with_curated_pkn_edges": joined(data["compatible_current_curated_edges"]),
            "compatible_candidates_with_exploratory_pkn_edges": joined(data["compatible_exploratory_edges"]),
            "has_any_exact_or_compatible_candidate": bool(specific),
            "has_any_signed_exact_or_compatible_candidate": bool(signed_specific),
            "has_any_broad_class_annotation": bool(data["broad_class"]),
            "annotated_chemical_ids": joined(data["chemical_ids"]),
            "formula_incompatible_resolution_rows_excluded": incompatible[name],
            "measurement_identity_changed": False,
            "coverage_scope": "union_of_frozen_crosswalk_and_limited_ontology_specificity_audit_not_exhaustive",
        })
    panels = {
        "all_measurements": features,
        "accepted_155_inputs": {r["feature"] for r in rows if r["accepted_frozen_input"]},
        "historical_64_factor1_lipids": {r["feature"] for r in rows if r["in_historical_factor1_64"]},
        "two_lpc_candidates": {r["feature"] for r in rows if r["two_lpc_candidate"]},
    }
    assert [len(v) for v in panels.values()] == [1254, 155, 64, 2]
    assert not (current_exact & panels["accepted_155_inputs"])
    summary_rows, union_rows = [], []
    for panel, selected in panels.items():
        for eg in sorted(evidence_groups):
            for category in ("exact", "compatible", "broad_class", "new_resolution", "exact_or_compatible_union"):
                hits = {feature for feature in selected if eg in
                        (state[feature]["exact"] | state[feature]["compatible"] if category == "exact_or_compatible_union"
                         else state[feature][category])}
                summary_rows.append({"panel": panel, "panel_features": len(selected), "evidence_group": eg,
                    "category": category, "unique_features": len(hits), "feature_ids": joined(hits)})
        exact = {f for f in selected if state[f]["exact"] - {UNSIGNED}}
        compatible = {f for f in selected if state[f]["compatible"] - {UNSIGNED}}
        broad = {f for f in selected if state[f]["broad_class"] - {UNSIGNED}}
        unsigned = {f for f in selected if UNSIGNED in state[f]["exact"] | state[f]["compatible"]}
        new = {f for f in selected if state[f]["new_resolution"]}
        union_rows.append({"panel": panel, "panel_features": len(selected),
            "signed_frozen_exact_features": len(exact), "signed_compatible_candidate_features": len(compatible),
            "signed_exact_or_compatible_unique_union": len(exact | compatible),
            "unsigned_exact_or_compatible_features": len(unsigned),
            "all_exact_or_compatible_unique_union": len(exact | compatible | unsigned),
            "signed_broad_class_features": len(broad),
            "signed_any_annotation_unique_union": len(exact | compatible | broad),
            "new_resolution_candidate_features": len(new),
            "current_curated_exact_features": len(current_exact & selected),
            "specific_candidate_feature_ids": joined(exact | compatible | unsigned)})
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "all_measured_lipid_regulatory_coverage.csv", rows, list(rows[0]))
    write_csv(out / "per_resource_coverage_counts.csv", summary_rows, list(summary_rows[0]))
    write_csv(out / "unique_feature_unions.csv", union_rows, list(union_rows[0]))
    dependency_rows = [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                       for path in dependencies]
    write_csv(out / "dependency_hashes.csv", dependency_rows, ["path", "sha256"])
    compact_resource = []
    for eg in sorted(evidence_groups):
        values = {r["category"]: r["unique_features"] for r in summary_rows
                  if r["panel"] == "all_measurements" and r["evidence_group"] == eg}
        compact_resource.append(dict(evidence_group=eg, **values))
    summary = {"analysis": "consolidated_lipid_regulatory_coverage_annotations",
        "features": len(rows), "annotation_only": True, "mapping_changed": False, "pkn_changed": False,
        "moon_refit": False, "panel_unions": union_rows, "per_evidence_group_all_measurements": compact_resource,
        "formula_incompatible_candidate_rows_excluded": sum(incompatible.values()),
        "cellphonedb": "Not assessed: no CellPhoneDB records were available in the frozen source collections.",
        "detailed_audits": {"frozen_crosswalk": str(frozen_dir), "specificity_candidates": str(resolution_dir)},
        "interpretation": [
            "Exact denotes the frozen mapper's reported-resolution classification, with original identity and acid/base labels retained. FA16 is a nested generic/specific identity ambiguity and was not admitted to the frozen155 inputs.",
            "Compatible combines frozen compatible matches and formula-compatible supplementary ontology descendants. These are candidate molecules, not confirmed measured identities or automatic mappings.",
            "GtoPdb primary and secondary chemical mapping tiers remain distinct; BRENDA/STITCH are cached evidence; MRCLinksDB is unsigned association data.",
            "A feature may have both exact and compatible records and may match multiple resources. Union counts deduplicate features and must replace sums across categories/resources.",
            "Pkn-edge availability columns preserve the detailed specificity audit's compartment-aware checks; candidate identity is still unconfirmed.",
            "No accepted155 feature has an exact direct current-curated-PKN regulatory match. Candidate-level coverage across all1254 measurements must not be presented as coverage of the accepted input panel.",
            "Coverage remains incomplete: the supplementary search only examines measured features already carrying exact reported-resolution ChEBI terms; neither audit establishes exhaustive biological coverage.",
        ]}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "README.md").write_text(
        "# One row per measured lipid: regulatory-coverage overview\n\n"
        "`all_measured_lipid_regulatory_coverage.csv` is the consolidated lookup for all1,254 measurements. "
        "It preserves the accepted155 and historical64 panel flags, original measurement identity and weight, "
        "frozen exact matches, conditional compatible candidates, broad-class annotations and newly detected "
        "specificity candidates. Formula-incompatible supplementary candidates are excluded. This is annotation "
        "only: no measurement, PKN or fit is changed.\n\n"
        "Detailed provenance remains in the [frozen-crosswalk audit](../lipid_regulatory_mapping_audit/README.md) "
        "and [specificity audit](../lipid_regulatory_resolution_audit/README.md).\n\n"
        "Evidence-group names keep primary-curated GtoPdb, secondary chemical mapping, cached BRENDA/STITCH and "
        "unsigned MRCLinksDB associations separate. CellPhoneDB coverage was not assessable from the frozen files. "
        "`per_resource_coverage_counts.csv` reports each group; `unique_feature_unions.csv` deduplicates overlapping "
        "features across categories and groups. Do not sum exact and compatible counts: FA16, for example, occurs "
        "in both while remaining excluded from the accepted155 inputs.\n\n"
        "A compatible candidate is not a confirmed measured molecular identity. The current accepted155 inputs "
        "have no exact direct regulatory matches in the current curated PKN. Whole-dataset candidate opportunities "
        "are a different quantity. Both detailed audits have limited search scope, so unmatched features do not "
        "prove absent biological knowledge.\n\n"
        "Reproduce with `python3 scripts/moon/ancillary/summarize_lipid_regulatory_coverage.py --output NEW_EMPTY_DIRECTORY`. "
        "Occupied output directories are refused, and `dependency_hashes.csv` records both input audits.\n")
    print(json.dumps({"output": str(out), "features": len(rows), "panel_unions": union_rows,
                      "per_evidence_group_all_measurements": compact_resource}, indent=2))


if __name__ == "__main__":
    main()
