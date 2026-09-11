#!/usr/bin/env python3
"""Instantiate conservative measured-lipid headgroup reaction hypotheses.

Run from the project root: python3 scripts/pkn/build_lipid_rules.py
Standard library only. Requires frozen LINEX2 source, the current Rhea mouse
reaction export and mapping_features.csv. Use --self-test for matching checks.

The matching rule follows LINEX2's (1,1) headgroup-transformation criterion of
equal fatty-acid multisets, but preserves full current Rhea chemistry and uses
only measured lipids with two explicitly reported, unmodified ester acyl chains.
It does not infer acyl composition from sum species, use the saved adjacency,
promote class reactions to curated species facts, or infer missing reverses.
"""

import argparse
import collections
import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/cosmos_mouse_pkn"
RAW = ROOT / "support/cosmos_mouse/raw"
COMMIT = "04e0099655a6b7a8cb57aaeca7b1301b3524079d"

# Explicitly audited generic 1,2-diacyl classes, rather than blindly trusting
# LINEX2's FixedChain field (which is blank for some chemically fixed species).
# Rhea formula must retain exactly two R substituents before a class is used.
GENERIC = {
    "CHEBI:58608": "PA", "CHEBI:17815": "DG", "CHEBI:49172": "DG",
    "CHEBI:57643": "PC", "CHEBI:64482": "PC", "CHEBI:64612": "PE",
    "CHEBI:57262": "PS", "CHEBI:57880": "PI", "CHEBI:64716": "PG",
    "CHEBI:60523": "PG",
}
CLASS_ALIASES = {"DAG": "DG"}
FEATURE_PATTERN = re.compile(r"^(PC|PE|PS|PA|DG|PI|PG) (\d+:\d+)([_/])(\d+:\d+)$")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archive_rows(archive, suffix):
    members = [m for m in archive.getmembers() if m.name.endswith("/" + suffix)]
    if len(members) != 1:
        raise ValueError("Expected one pinned LINEX2 file: " + suffix)
    with io.TextIOWrapper(archive.extractfile(members[0]), encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_feature(row):
    match = FEATURE_PATTERN.fullmatch(row["canonical_name"])
    if not match or row.get("parse_error"):
        return None
    chains = [match[2], match[4]]
    return {"feature_id": row["feature_id"], "class": match[1],
            "chains": chains, "chain_multiset": sorted(chains),
            "sn_order_reported": match[3] == "/", "canonical_name": row["canonical_name"],
            "in_factor1": row["in_factor1"] == "True", "in_lipid_de": row["in_lipid_de"] == "True"}


def compatible(a, b):
    if a["chain_multiset"] != b["chain_multiset"]:
        return False
    # If both features report sn order, a pure headgroup reaction cannot swap it.
    return not (a["sn_order_reported"] and b["sn_order_reported"] and a["chains"] != b["chains"])


def replace_participant(participant, feature):
    result = copy.deepcopy(participant)
    result.update({"id": "MeasuredLipid:" + feature["feature_id"], "name": feature["feature_id"],
                   "template_chemical_id": participant["id"], "template_name": participant["name"],
                   "template_formula": participant["formula"], "formula": None,
                   "template_rhea_compound": participant["rhea_compound"], "rhea_compound": None,
                   "chemical_type": ["measured molecular species; rule-instantiated"],
                   "measurement_resolution": "reported acyl composition; unreported isomers unresolved",
                   "chain_multiset": feature["chain_multiset"], "sn_order_reported": feature["sn_order_reported"],
                   "lipid_classification_evidence": [{"basis": "measured lipid class plus exact-chain rule instantiation",
                                                       "template_chemical_id": participant["id"]}],
                   "is_lipid": True})
    return result


def instantiate(template, sub_index, prod_index, source, target, definitions):
    record = copy.deepcopy(template)
    signature = json.dumps([template["reaction_id"], source["feature_id"], target["feature_id"]])
    record.update({
        "reaction_id": "LIPID_RULE:" + hashlib.sha256(signature.encode()).hexdigest()[:20],
        "resource": "LINEX2_rule", "source_version": "LINEX2 " + COMMIT + "; Rhea " + template["source_version"],
        "supporting_reaction_id": template["reaction_id"], "supporting_master_reaction_id": template["master_reaction_id"],
        "template_resource": "Rhea", "catalyst_namespace": "UniProtKB",
        "hypothesis_type": "rule-instantiated species hypothesis; not a curated species reaction",
        "reference_scope": "Source generic reaction/catalyst annotations; references do not establish activity on these measured species",
        "chemical_balance_scope": "Source generic Rhea template; identical two-chain substitutions preserve the formal class stoichiometry",
        "headgroup_transform": source["class"] + " -> " + target["class"],
        "chain_transform": "identical unordered acyl-chain multiset; headgroup change only",
        "matched_feature_ids": [source["feature_id"], target["feature_id"]],
        "matched_factor1_features": [f["feature_id"] for f in (source, target) if f["in_factor1"]],
        "chain_multiset": source["chain_multiset"],
        "source_template_equation": template["equation"], "equation": None,
        "rule_definitions": definitions,
        "rule_implementation_reference": "https://github.com/LipiTUM/LINEX2_package/blob/" + COMMIT + "/linex2/reaction_evaluation.py",
        "assumptions": [
            "The generic headgroup reaction can act on the measured acyl-chain combination; this substrate specificity is hypothetical.",
            "Unreported sn positions, double-bond positions and stereochemistry remain unresolved.",
            "A graph projection represents permissive traversal; it does not enforce joint substrate availability or flux balance.",
        ],
    })
    record["substrates"][sub_index] = replace_participant(template["substrates"][sub_index], source)
    record["products"][prod_index] = replace_participant(template["products"][prod_index], target)
    record["substitutions"] = [
        {"side": "substrates", "index": sub_index, "template_chemical_id": template["substrates"][sub_index]["id"],
         "measured_id": "MeasuredLipid:" + source["feature_id"]},
        {"side": "products", "index": prod_index, "template_chemical_id": template["products"][prod_index]["id"],
         "measured_id": "MeasuredLipid:" + target["feature_id"]},
    ]
    # Every non-substituted participant, side count, coefficient and catalyst is
    # unchanged. The source Rhea record remains the authoritative hyperedge.
    for side, changed in [("substrates", sub_index), ("products", prod_index)]:
        assert len(record[side]) == len(template[side])
        for i, original in enumerate(template[side]):
            assert record[side][i]["stoichiometry"] == original["stoichiometry"]
            if i != changed:
                assert record[side][i] == original
    assert record["catalysts"] == template["catalysts"]
    return record


def main():
    commit_metadata = json.loads((RAW / "linex-github-commit.json").read_text())
    if commit_metadata["sha"] != COMMIT:
        raise ValueError("LINEX2 archive commit differs from the audited version")
    with tarfile.open(RAW / "linex-source.tar.gz") as archive:
        definitions = archive_rows(archive, "linex2/data/Rhea_reactions_for_curation.csv")
        molecules = archive_rows(archive, "linex2/data/Rhea_Others_all_unique.csv")
    generic = {}
    for row in molecules:
        ident = row["ChEBI"]
        if ident in GENERIC:
            if row["isLipid"] != "1" or row["FixedChain"] or CLASS_ALIASES.get(row["Abbreviation"], row["Abbreviation"]) != GENERIC[ident]:
                raise ValueError("Unexpected generic LINEX2 annotation: " + str(row))
            if "diacyl" not in row["name1"].lower():
                raise ValueError("Generic molecule not explicitly diacyl: " + str(row))
            generic[ident] = row
    by_master = collections.defaultdict(list)
    for i, row in enumerate(definitions, start=2):
        if row["reaction_type"] not in ("L_HGmodify", "L_HGdelete"):
            continue
        evidence = {"source_file": "linex2/data/Rhea_reactions_for_curation.csv", "source_row": i,
                    "source_commit": COMMIT, "reaction_type": row["reaction_type"],
                    "reaction_ids": row["reaction_IDs"], "class_equation": row["merged_reaction"],
                    "chemical_ids": row["merged_ids"], "notes": row["reaction_notes"]}
        for rid in row["reaction_IDs"].split(";"):
            by_master[rid.strip()].append(evidence)
    with (OUT / "mapping_features.csv").open() as handle:
        all_features = list(csv.DictReader(handle))
    features = [f for row in all_features if (f := parse_feature(row))]
    by_class = collections.defaultdict(list)
    for feature in features:
        by_class[feature["class"]].append(feature)
    reactions, templates, excluded = [], [], collections.Counter()
    for template in map(json.loads, (OUT / "rhea_reactions.jsonl").read_text().splitlines()):
        if template["master_reaction_id"] not in by_master:
            continue
        if template["direction"] not in ("LR", "RL") or template["is_transport"]:
            excluded["unknown_direction_or_transport"] += 1; continue
        candidates = [[i for i, p in enumerate(template[side]) if p["id"] in generic] for side in ("substrates", "products")]
        if len(candidates[0]) != 1 or len(candidates[1]) != 1:
            excluded["not_one_generic_diacyl_class_per_side"] += 1; continue
        sub_i, prod_i = candidates[0][0], candidates[1][0]
        sub, prod = template["substrates"][sub_i], template["products"][prod_i]
        source_class, target_class = GENERIC[sub["id"]], GENERIC[prod["id"]]
        if source_class == target_class:
            excluded["unchanged_headgroup"] += 1; continue
        if any(p["is_lipid"] for side, changed in [("substrates", sub_i), ("products", prod_i)]
               for i, p in enumerate(template[side]) if i != changed):
            excluded["additional_lipid_participant_requires_joint_instantiation"] += 1; continue
        if any(p["stoichiometry"] != 1 or not re.search(r"R2$", p["formula"] or "") for p in (sub, prod)):
            excluded["nonunit_or_non_generic_two_acyl_substituents"] += 1; continue
        matched = 0
        for source in by_class[source_class]:
            for target in by_class[target_class]:
                if compatible(source, target):
                    reactions.append(instantiate(template, sub_i, prod_i, source, target, by_master[template["master_reaction_id"]]))
                    matched += 1
        templates.append({"reaction_id": template["reaction_id"], "master_reaction_id": template["master_reaction_id"],
                          "transform": source_class + " -> " + target_class, "catalysts": template["catalysts"],
                          "instantiated_pairs": matched, "generic_substrate_id": sub["id"], "generic_product_id": prod["id"]})
    with (OUT / "lipid_rules_reactions.jsonl").open("w") as handle:
        for record in reactions:
            handle.write(json.dumps(record, sort_keys=True) + "\n")
    with (OUT / "lipid_rules_pairs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["reaction_id", "supporting_reaction_id", "headgroup_transform", "source_feature", "target_feature", "chain_multiset", "catalysts", "factor1_features"])
        writer.writeheader()
        for r in reactions:
            writer.writerow({"reaction_id": r["reaction_id"], "supporting_reaction_id": r["supporting_reaction_id"],
                             "headgroup_transform": r["headgroup_transform"], "source_feature": r["matched_feature_ids"][0],
                             "target_feature": r["matched_feature_ids"][1], "chain_multiset": "_".join(r["chain_multiset"]),
                             "catalysts": "|".join(r["catalysts"]), "factor1_features": "|".join(r["matched_factor1_features"])})
    summary = {"rule_instances": len(reactions), "supporting_directional_reactions": len({r["supporting_reaction_id"] for r in reactions}),
               "supporting_reaction_masters": len({r["supporting_master_reaction_id"] for r in reactions}),
               "unique_measured_features": len({f for r in reactions for f in r["matched_feature_ids"]}),
               "factor1_features": sorted({f for r in reactions for f in r["matched_factor1_features"]}),
               "mouse_catalysts": len({c for r in reactions for c in r["catalysts"]}),
               "headgroup_transforms": dict(collections.Counter(r["headgroup_transform"] for r in reactions)),
               "eligible_measured_diacyl_features": len(features), "excluded_candidate_reasons": dict(excluded),
               "templates": templates,
               "interpretation": "Separate exploratory species rules. Exact reported acyl chains matched; enzyme species specificity, isomers and unreported sn positions unresolved.",
               "not_implemented": "Acyl removal/addition, elongation/desaturation, ethers, sphingolipids, sum-composition expansion, multi-lipid coupled transformations or unmeasured intermediates."}
    (OUT / "lipid_rules_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    manifest = {"linex_commit": COMMIT, "linex_commit_date": commit_metadata["commit"]["committer"]["date"],
                "linex_repository": "https://github.com/LipiTUM/LINEX2_package",
                "archive_url": "https://api.github.com/repos/LipiTUM/LINEX2_package/tarball/" + COMMIT,
                "inputs": [{"file": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in [RAW / "linex-source.tar.gz", RAW / "linex-github-commit.json", OUT / "rhea_reactions.jsonl", OUT / "mapping_features.csv"]],
                "generic_chemical_allowlist": GENERIC,
                "data_provenance": "LINEX2 includes Rhea class annotations; current catalysts, directional chemistry and full participants supplied by Rhea 142 mouse export.",
                "license": "LINEX2 software AGPLv3; underlying Rhea data CC BY 4.0 (see pinned source license)."}
    (OUT / "lipid_rules_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k:v for k,v in summary.items() if k != "templates"}, indent=2))


def self_test():
    def f(name):
        return parse_feature({"feature_id": name, "canonical_name": name, "in_factor1": "False", "in_lipid_de": "True"})
    assert f("PC 36:2") is None and f("PC O-16:0_18:1") is None and f("PE 16:0_18:1;O") is None
    assert compatible(f("PC 16:0_18:1"), f("PE 18:1_16:0"))
    assert not compatible(f("PC 16:0_18:2"), f("PE 18:1_16:1"))  # same totals, different chains
    assert not compatible(f("PC 16:0/18:1"), f("PE 18:1/16:0"))  # reported sn swap
    assert compatible(f("PC 16:0/18:1"), f("PE 18:1_16:0"))  # compatible unknown sn
    print("Lipid rule matching self-tests passed: sums, ether/oxidation exclusions, exact chain multiset and reported sn consistency.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    if parser.parse_args().self_test:
        self_test()
    else:
        main()
