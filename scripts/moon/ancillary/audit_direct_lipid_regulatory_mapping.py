#!/usr/bin/env python3
"""Audit frozen measured-feature matches to chemical/protein resource records.

This is coverage of the EXISTING crosswalk, not an exhaustive chemical mapping.
No input identifiers, resources, prior-network edges or MOON fits are changed.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RESOURCES = ("BRENDA", "STITCH", "MRCLinksDB", "CellPhoneDB", "GuideToPharmacology")
CATEGORIES = ("exact_same_identifier", "exact_explicit_alias", "exact_conjugate_acid_base",
              "compatible_resolution", "class_or_ancestor")
SCOPE = "frozen_crosswalk_only_not_exhaustive_chemical_mapping"


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalize_identifier(value):
    value = str(value)
    return "CHEBI:" + value if value.isdigit() else value


def joined(values):
    return "|".join(sorted({str(v) for v in values if v not in (None, "")}))


class Equivalence:
    def __init__(self):
        self.parent = {}

    def find(self, value):
        value = normalize_identifier(value)
        self.parent.setdefault(value, value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def join(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def ontology_terms(path, wanted):
    result, block, current = {}, [], None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.strip() in ("[Term]", "[Typedef]"):
                if current in wanted:
                    result[current] = block
                block, current = [], None
            elif line.startswith("id: "):
                current = line[4:].strip()
                block.append(line.rstrip())
            elif current in wanted:
                block.append(line.rstrip())
        if current in wanted:
            result[current] = block
    assert set(result) == set(wanted)
    parsed = []
    for chemical, lines in sorted(result.items()):
        get = lambda prefix: next((x[len(prefix):] for x in lines if x.startswith(prefix)), "")
        prop = lambda name: next((re.search(r'"(.*?)"', x).group(1) for x in lines
                                 if x.startswith("property_value: chemrof:" + name + " ")), "")
        parsed.append({"chemical_id": chemical, "name": get("name: "),
            "definition": get("def: "), "is_a": joined(x[6:] for x in lines if x.startswith("is_a: ")),
            "relationships": joined(x[14:] for x in lines if x.startswith("relationship: ")),
            "formula": prop("generalized_empirical_formula"), "charge": prop("charge"),
            "smiles": prop("smiles_string"), "inchi_key": prop("inchi_key_string")})
    return parsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results/causal_modelling/lipid_regulatory_mapping_audit")
    parser.add_argument("--pkn", type=Path, default=ROOT / "results/causal_modelling/cosmos_mouse_pkn")
    parser.add_argument("--inputs", type=Path, default=ROOT / "miscellaneous/older_results/causal_modelling/cosmos_factor1_pilot/inputs")
    parser.add_argument("--cache", type=Path, default=ROOT / "results/causal_modelling/cosmos_pkn_audit/raw/om-mouse-pkn.json.gz")
    args = parser.parse_args()
    out, pkn, inputs, cache_path = (p.resolve() for p in (args.output, args.pkn, args.inputs, args.cache))
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error(f"Output directory is occupied: {out}")
    ontology_path = ROOT / "support/cosmos_mouse/raw/regulation_chebi.obo.gz"
    dependencies = [pkn / name for name in (
        "mapping_features.csv", "mapping_chebi_nodes.csv", "mapping_candidates.csv",
        "regulation_edges.csv", "regulation_exclusions.csv", "regulation_chebi_aliases.csv",
        "regulation_chebi_acid_base_relations.csv", "regulation_chebi_source_crosswalk.csv",
        "pkn_lipid_rules.csv",
    )] + [inputs / "factor1_lipid_node_mapping.csv", cache_path, ontology_path,
          ROOT / "scripts/pkn/build_lipid_regulation.py", Path(__file__).resolve()]
    for path in dependencies:
        if not path.is_file():
            parser.error(f"Missing frozen dependency: {path}")
    features = read_csv(pkn / "mapping_features.csv")
    feature_ids = {r["feature_id"] for r in features}
    frozen = {r["feature"]: r for r in read_csv(inputs / "factor1_lipid_node_mapping.csv")}
    panels = {
        "all_1254_measurements": feature_ids,
        "accepted_155_inputs": {name for name, r in frozen.items() if r["included"] == "TRUE"},
        "historical_64_factor1_lipids": {r["feature_id"] for r in features if r["in_factor1"] == "True"},
        "two_lpc_candidates": {"LPC(16:0/0:0)", "LPC(18:0/0:0)"},
    }
    assert [len(v) for v in panels.values()] == [1254, 155, 64, 2]
    candidates = [r for r in read_csv(pkn / "mapping_chebi_nodes.csv")
                  if r["mapping_relation"] == "SwissLipids identifier"]
    alias, acidbase = Equivalence(), Equivalence()
    for row in read_csv(pkn / "regulation_chebi_aliases.csv"):
        alias.join(row["alias"], row["primary_chebi"])
        acidbase.join(row["alias"], row["primary_chebi"])
    for row in read_csv(pkn / "regulation_chebi_acid_base_relations.csv"):
        acidbase.join(row["acid"], row["base"])
    regulations = read_csv(pkn / "regulation_edges.csv")
    retained_ids = {r["source_record_id"] for r in regulations}
    exclusions = {r["source_record_id"]: r for r in read_csv(pkn / "regulation_exclusions.csv")}
    crosswalk = {r["source"]: r for r in read_csv(pkn / "regulation_chebi_source_crosswalk.csv")}
    pkn_edges = {(r["source"], int(r["interaction"]), r["target"])
                 for r in read_csv(pkn / "pkn_lipid_rules.csv")}
    with gzip.open(cache_path, "rt") as handle:
        cache = json.load(handle)
    assert cache["meta"]["organism"] == 10090
    gene_by_uniprot = {r["target"]: r["target_gene"] for r in regulations if r["target_gene"]}
    records = []
    for index, row in enumerate(cache["network"]):
        resource = row["resource"]
        if resource not in RESOURCES:
            continue
        record_id = "OmniPath-metabo:mouse-cache:record:" + str(index)
        source = normalize_identifier(row["source"])
        sign = row.get("mor")
        mode = row.get("attrs", {}).get("stitch_mode")
        gate = record_id in retained_ids
        sign_status = ("adapter_assumed_positive_not_signed_regulatory_evidence" if resource == "MRCLinksDB" else
                       "cached_explicit_action_direction_not_independently_auditable" if resource == "STITCH" else
                       "cached_signed_allosteric_record_original_evidence_unavailable")
        records.append({"collection": "raw_mouse_service_cache", "resource": resource,
            "evidence_tier": "raw_cached_" + resource.lower(), "record_id": record_id,
            "chemical_id": source, "raw_source_identifier": str(row["source"]),
            "chemical_name": crosswalk.get(source, {}).get("name", ""),
            "target_uniprot": row["target"], "target_gene": gene_by_uniprot.get(row["target"], ""),
            "interaction": sign, "original_action": mode or "",
            "interaction_type": row.get("interaction_type", ""),
            "source_type": row.get("source_type", ""), "target_type": row.get("target_type", ""),
            "organism": 10090, "passes_existing_regulation_export_gate": gate,
            "in_current_curated_pkn": False, "sign_status": sign_status,
            "references": "", "raw_attributes_json": json.dumps(row.get("attrs", {}), sort_keys=True),
            "existing_exclusion_reason": exclusions.get(record_id, {}).get("reason", ""),
            "source_file": str(cache_path)})
    for row in regulations:
        source = row["source"]
        key = ("Metab__" + source + "_unknown", int(row["interaction"]), "Gene__" + row["target_gene"])
        records.append({"collection": "retained_regulation_export", "resource": row["resource"],
            "evidence_tier": row["evidence_tier"], "record_id": row["source_record_id"],
            "chemical_id": source, "raw_source_identifier": source,
            "chemical_name": row["source_name"] or crosswalk.get(source, {}).get("name", ""),
            "target_uniprot": row["target"], "target_gene": row["target_gene"],
            "interaction": int(row["interaction"]), "original_action": row["original_action"],
            "interaction_type": row["original_type"], "source_type": row["source_type"],
            "target_type": row["target_type"], "organism": row["organism"],
            "passes_existing_regulation_export_gate": True,
            "in_current_curated_pkn": row["evidence_tier"] == "curated_mouse_current" and key in pkn_edges,
            "sign_status": "curated_mouse_activity_action" if row["resource"] == "GuideToPharmacology" else
                "cached_explicit_action_direction_not_independently_auditable" if row["resource"] == "STITCH" else
                "cached_signed_allosteric_record_original_evidence_unavailable",
            "references": row["references"], "raw_attributes_json": "",
            "existing_exclusion_reason": "", "source_file": str(pkn / "regulation_edges.csv")})

    component_records = collections.defaultdict(list)
    for row in records:
        if re.fullmatch(r"CHEBI:[0-9]+", row["chemical_id"]):
            component_records[acidbase.find(row["chemical_id"])].append(row)
    matches = {}
    for candidate in candidates:
        chemical = candidate["original_chebi_id"]
        for record in component_records.get(acidbase.find(chemical), ()):
            source = record["chemical_id"]
            relation = "same_identifier" if source == chemical else "explicit_alias" if alias.find(source) == alias.find(chemical) else "conjugate_acid_base"
            resolution = candidate["match_type"]
            category = ("exact_" + relation if resolution == "exact_reported_resolution" else
                        "class_or_ancestor" if resolution.startswith("hierarchy_") else "compatible_resolution")
            key = (record["collection"], record["resource"], record["record_id"], candidate["feature_id"], chemical, resolution)
            result = matches.setdefault(key, dict(record, feature=candidate["feature_id"],
                mapped_chemical=chemical, match_type=resolution, identity_relation=relation,
                coverage_category=category, swisslipids_ids=set(),
                accepted_frozen_input=candidate["feature_id"] in panels["accepted_155_inputs"],
                historical_factor1_lipid=candidate["feature_id"] in panels["historical_64_factor1_lipids"],
                lpc_candidate=candidate["feature_id"] in panels["two_lpc_candidates"],
                frozen_mapping_exclusion=frozen[candidate["feature_id"]]["exclusion_reason"],
                coverage_scope=SCOPE))
            result["swisslipids_ids"].add(candidate["swisslipids_id"])
    match_rows = []
    for key, row in sorted(matches.items()):
        row["swisslipids_ids"] = joined(row["swisslipids_ids"])
        match_rows.append(row)
    collections_available = ("raw_mouse_service_cache", "retained_regulation_export")
    counts, sources, funnels, per_feature = [], [], [], []
    for collection in collections_available:
        for resource in RESOURCES:
            rr = [r for r in records if r["collection"] == collection and r["resource"] == resource]
            mm = [r for r in match_rows if r["collection"] == collection and r["resource"] == resource]
            availability = "present_in_frozen_collection" if rr else "not_assessable_no_resource_records_in_this_collection"
            sources.append({"collection": collection, "resource": resource, "availability": availability,
                "records": len(rr), "chemical_sources": len({r["chemical_id"] for r in rr}),
                "target_proteins": len({r["target_uniprot"] for r in rr}),
                "positive_numeric_mor_records": sum(r["interaction"] == 1 for r in rr),
                "negative_numeric_mor_records": sum(r["interaction"] == -1 for r in rr),
                "passes_existing_regulation_export_gate_records": sum(r["passes_existing_regulation_export_gate"] for r in rr),
                "mapped_features_any_resolution": len({r["feature"] for r in mm}),
                "numeric_sign_caveat": "Numeric mor does not by itself establish signed causal regulation; inspect sign_status and evidence tier."})
            by_feature = collections.defaultdict(list)
            for row in mm:
                by_feature[row["feature"]].append(row)
            for feature in sorted(feature_ids):
                fm = by_feature[feature]
                categories = {r["coverage_category"] for r in fm}
                exact = [r for r in fm if r["coverage_category"].startswith("exact_")]
                if not rr:
                    reason = "resource_not_available_in_frozen_collection"
                elif exact and resource == "MRCLinksDB":
                    reason = "exact_candidate_but_assumed_sign_rejected"
                elif exact and not any(r["passes_existing_regulation_export_gate"] for r in exact):
                    reason = "exact_candidate_but_existing_evidence_gate_rejected"
                elif exact:
                    reason = "exact_candidate_opportunity_requires_existing_mapping_and_resource_policy_review"
                elif "compatible_resolution" in categories:
                    reason = "only_compatible_or_more_specific_candidates"
                elif "class_or_ancestor" in categories:
                    reason = "only_broad_class_or_ancestor_candidates"
                else:
                    reason = "no_match_in_frozen_crosswalk_not_proof_of_absent_chemical_knowledge"
                per_feature.append({"collection": collection, "resource": resource, "feature": feature,
                    "availability": availability, "reason": reason,
                    "accepted_frozen_input": feature in panels["accepted_155_inputs"],
                    "historical_factor1_lipid": feature in panels["historical_64_factor1_lipids"],
                    "lpc_candidate": feature in panels["two_lpc_candidates"],
                    "coverage_categories": joined(categories),
                    "matched_resource_chemicals": joined(r["chemical_id"] for r in fm),
                    "matched_records": len({r["record_id"] for r in fm}),
                    "exact_candidate_records": len({r["record_id"] for r in exact}),
                    "frozen_mapping_exclusion": frozen[feature]["exclusion_reason"], "coverage_scope": SCOPE})
            for panel, panel_features in panels.items():
                subset = [r for r in mm if r["feature"] in panel_features]
                for category in CATEGORIES + ("any_resolution", "any_exact_including_acid_base"):
                    selected = [r for r in subset if category == "any_resolution" or
                                category == r["coverage_category"] or
                                (category == "any_exact_including_acid_base" and r["coverage_category"].startswith("exact_"))]
                    counts.append({"collection": collection, "resource": resource, "panel": panel,
                        "panel_features": len(panel_features), "coverage_category": category,
                        "availability": availability, "unique_features": len({r["feature"] for r in selected}),
                        "resource_chemicals": len({r["chemical_id"] for r in selected}),
                        "unique_resource_records": len({r["record_id"] for r in selected}),
                        "exact_current_curated_graph_features": len({r["feature"] for r in selected
                            if r["coverage_category"].startswith("exact_") and r["in_current_curated_pkn"]}),
                        "feature_ids": joined(r["feature"] for r in selected), "coverage_scope": SCOPE})
                feature_reasons = collections.Counter(r["reason"] for r in per_feature
                    if r["collection"] == collection and r["resource"] == resource and r["feature"] in panel_features)
                for reason, number in sorted(feature_reasons.items()):
                    funnels.append({"collection": collection, "resource": resource, "panel": panel,
                                    "reason": reason, "unique_features": number, "coverage_scope": SCOPE})

    palmitate = ontology_terms(ontology_path, {"CHEBI:78123", "CHEBI:7896", "CHEBI:15756"})
    assert any(r["chemical_id"] == "CHEBI:7896" and "CHEBI:78123" in r["is_a"] for r in palmitate)
    palmitate_candidates = [r for r in read_csv(pkn / "mapping_candidates.csv")
                            if r["feature_id"] == "FA(16:0)" and r["match_type"] == "exact_reported_resolution"]
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "feature_edge_matches.csv", match_rows, list(match_rows[0]))
    write_csv(out / "coverage_counts.csv", counts, list(counts[0]))
    write_csv(out / "resource_source_coverage.csv", sources, list(sources[0]))
    write_csv(out / "feature_reason_funnel.csv", per_feature, list(per_feature[0]))
    write_csv(out / "reason_funnel_counts.csv", funnels, list(funnels[0]))
    write_csv(out / "palmitate_ontology.csv", palmitate, list(palmitate[0]))
    write_csv(out / "palmitate_swisslipids_candidates.csv", palmitate_candidates, list(palmitate_candidates[0]))
    matched_record_keys = {(r["collection"], r["record_id"]) for r in match_rows}
    matched_records = [r for r in records if (r["collection"], r["record_id"]) in matched_record_keys]
    write_csv(out / "matched_source_records.csv", matched_records, list(records[0]))
    inventory = []
    for row in features:
        inventory.append(dict(row, accepted_frozen_input=row["feature_id"] in panels["accepted_155_inputs"],
            frozen_node_id=frozen[row["feature_id"]]["node_id"],
            frozen_exclusion_reason=frozen[row["feature_id"]]["exclusion_reason"],
            oriented_weight=frozen[row["feature_id"]]["oriented_weight"], coverage_scope=SCOPE))
    write_csv(out / "measurement_inventory.csv", inventory, list(inventory[0]))
    dependency_rows = [{"path": str(path), "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in dependencies]
    write_csv(out / "dependency_hashes.csv", dependency_rows, ["path", "bytes", "sha256"])
    summary = {"analysis": "direct_lipid_regulatory_frozen_mapping_audit", "coverage_scope": SCOPE,
        "panel_sizes": {name: len(values) for name, values in panels.items()},
        "resource_collections": sources, "feature_edge_match_rows": len(match_rows),
        "pk_network_or_measurement_mapping_changed": False, "refitting_performed": False,
        "key_limitations": [
            "Coverage uses existing mapping_chebi_nodes.csv candidates only. The original abbreviation-based SwissLipids search omits some structured entries with empty abbreviations; this audit is not an exhaustive ontology-descendant or synonym search.",
            "The separate lipid_regulatory_resolution_audit examines chemical specificity beyond the frozen candidate list; its findings must be considered before concluding that measured species lack prior knowledge.",
            "CellPhoneDB has no records in the frozen mouse cache or retained regulation export; a zero count here means not assessed from this collection, not absent biology or absence from CellPhoneDB.",
            "MRCLinksDB numeric +1 values are rejected by the existing builder as assumed signs, not signed activity evidence.",
            "Cached BRENDA/STITCH exports lack original chemical mappings and per-edge references. STITCH causal action orientation cannot be independently recovered from this cache.",
            "Class/ancestor and compatible more-specific matches never become measured identities or causal edges in this audit.",
            "FA(16:0)'s two exact SwissLipids hits are nested specificity levels: explicit hexadecanoate CHEBI:7896 is_a generic fatty-acid-anion16:0 CHEBI:78123. The frozen mapper treated their distinct IDs as ambiguity; that can be reviewed, but no remapping is performed here.",
        ]}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "README.md").write_text(
        "# Direct lipid-to-protein matches in the frozen mapping\n\n"
        "Run `python3 scripts/moon/ancillary/audit_direct_lipid_regulatory_mapping.py --output NEW_EMPTY_DIRECTORY`. "
        "This script compares the existing 1,254-feature crosswalk to raw cached BRENDA/STITCH/MRCLinksDB records and "
        "the retained regulation export, keeping those collections separate. No mapping, PKN or fit is changed.\n\n"
        "**Scope: frozen-crosswalk coverage, not an exhaustive chemical search.** The original abbreviation-based "
        "mapper misses some structured SwissLipids entries with empty abbreviations. The companion "
        "`lipid_regulatory_resolution_audit` examines specificity/ontology alternatives beyond these frozen candidates. "
        "These tables must not be used alone to claim that no prior chemical knowledge exists.\n\n"
        "`coverage_counts.csv` reports unique features for all measurements, the 155 accepted inputs, the historical "
        "64 Factor1 features and the two LPC candidates. Exact same-ID, explicit alias, explicit acid/base, compatible "
        "resolution and class/ancestor matches remain distinct. Categories can overlap; do not add counts to obtain "
        "a union. `feature_edge_matches.csv` preserves each matched source record with actions, signs, evidence gates "
        "and source IDs. `feature_reason_funnel.csv` includes unmatched features; its grouped version is "
        "`reason_funnel_counts.csv`.\n\n"
        "`resource_source_coverage.csv` separates raw and retained records. MRCLinksDB's positive numerical values "
        "are assumed by its adapter and were rejected as signed regulatory evidence. CellPhoneDB was not present "
        "in these frozen collections; its zero counts mean not assessable here. Raw cache provenance only establishes "
        "a mouse service export, not native-mouse evidence for every original source interaction.\n\n"
        "`palmitate_ontology.csv` and `palmitate_swisslipids_candidates.csv` show why the frozen FA16 ambiguity is a "
        "specificity issue: explicit straight-chain hexadecanoate is a subclass of the generic saturated C16 anion, "
        "rather than an unrelated alternative identity. Its acid/base overlap with palmitic acid is explicit; the "
        "existing input exclusion is preserved for this audit.\n\n"
        "All numerical counts are frozen-data audit results, not causal confidence or experimental evidence counts. "
        "`dependency_hashes.csv` identifies the consumed resources and implementation.\n")
    print(json.dumps({"output": str(out), "panel_sizes": summary["panel_sizes"],
                      "feature_edge_match_rows": len(match_rows), "coverage_scope": SCOPE}, indent=2))


if __name__ == "__main__":
    main()
