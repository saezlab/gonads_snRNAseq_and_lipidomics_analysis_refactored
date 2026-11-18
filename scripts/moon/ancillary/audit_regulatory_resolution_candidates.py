#!/usr/bin/env python3
"""Find specific regulatory chemicals beneath measured composition-level ChEBI terms.

Read-only annotation audit. Explicit ChEBI is_a relations, alt IDs and conjugate
acid/base records provide candidate membership, never identity or causal edges.
No input mapping, PKN or MOON fit is changed. Uses frozen project resources.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class Equivalence:
    def __init__(self):
        self.parent = {}

    def find(self, node):
        self.parent.setdefault(node, node)
        if self.parent[node] != node:
            self.parent[node] = self.find(self.parent[node])
        return self.parent[node]

    def join(self, left, right):
        a, b = self.find(left), self.find(right)
        self.parent[max(a, b)] = min(a, b)


def read_ontology(path):
    terms, current = {}, None
    with gzip.open(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("["):
                current = None
            elif line.startswith("id: CHEBI:"):
                current = {"id": line[4:], "name": "", "parents": [], "formula": "", "charge": ""}
                terms[current["id"]] = current
            elif current is not None:
                if line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("is_a: CHEBI:"):
                    current["parents"].append(line.split()[1])
                elif line.startswith('property_value: chemrof:generalized_empirical_formula '):
                    current["formula"] = line.split('"')[1]
                elif line.startswith('property_value: chemrof:charge '):
                    current["charge"] = line.split('"')[1]
    return terms


def fatty_acid_formula_status(feature, term):
    match = re.fullmatch(r"FA\((\d+):(\d+)\)", feature)
    if not match:
        return "not_a_simple_FA_formula_check"
    carbon, db = map(int, match.groups())
    try:
        charge = int(term["charge"])
    except ValueError:
        return "missing_formula_or_charge"
    if charge not in {-1, 0} or not term["formula"]:
        return "missing_formula_or_charge"
    expected = f"C{carbon}H{2 * carbon - 2 * db + charge}O2"
    return "matches_reported_FA_composition" if term["formula"] == expected else "formula_incompatible"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New or empty output directory")
    args = parser.parse_args()
    out = args.output.resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error(f"Occupied output directory: {out}")
    pkn = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
    raw_cache = ROOT / "results/causal_modelling/cosmos_pkn_audit/raw/om-mouse-pkn.json.gz"
    obo = ROOT / "support/cosmos_mouse/raw/regulation_chebi.obo.gz"
    uniprot = ROOT / "support/cosmos_mouse/raw/uniprot_mouse.tsv.gz"
    names = ["mapping_features.csv", "mapping_chebi_nodes.csv", "mapping_swisslipids_catalog.csv",
             "regulation_edges.csv", "regulation_chebi_aliases.csv", "regulation_chebi_acid_base_relations.csv",
             "pkn_lipid_rules.csv", "pkn_exploratory.csv"]
    dependencies = [pkn / name for name in names] + [obo, raw_cache, uniprot, Path(__file__).resolve()]
    assert all(path.is_file() for path in dependencies)
    features = read_csv(pkn / "mapping_features.csv")
    candidates = read_csv(pkn / "mapping_chebi_nodes.csv")
    catalog = read_csv(pkn / "mapping_swisslipids_catalog.csv")
    regulations = read_csv(pkn / "regulation_edges.csv")
    # Reproduce the assembler's accession/isoform handling without importing its
    # build entry point or mutating any graph.
    with gzip.open(uniprot, "rt") as handle:
        gene_map = {r["Entry"]: r["Gene Names (primary)"] for r in csv.DictReader(handle, delimiter="\t")
                    if re.fullmatch(r"[A-Za-z0-9_.-]+", r["Gene Names (primary)"])}
    def protein_node(accession):
        base = re.sub(r"-\d+$", "", accession)
        symbol = gene_map.get(accession, gene_map.get(base))
        return "Gene__" + symbol if symbol else "UniProt__" + accession
    graph_keys = {name: {(r["source"], int(r["interaction"]), r["target"]) for r in read_csv(pkn / filename)}
                  for name, filename in [("curated", "pkn_lipid_rules.csv"), ("exploratory", "pkn_exploratory.csv")]}
    # Preserve MRC associations as an unsigned source tier, never signed evidence.
    cache = json.load(gzip.open(raw_cache, "rt"))["network"]
    for index, row in enumerate(cache, 1):
        if row["resource"] == "MRCLinksDB":
            regulations.append(dict(source="CHEBI:" + row["source"].removeprefix("CHEBI:"),
                target=row["target"], target_gene="", interaction="", resource="MRCLinksDB",
                evidence_tier="excluded_unsigned_association", edge_id=f"cached_row_{index}",
                mechanism=row["interaction_type"], references=""))
    terms = read_ontology(obo)
    eq = Equivalence()
    for row in read_csv(pkn / "regulation_chebi_aliases.csv"):
        eq.join(row["alias"], row["primary_chebi"])
    for row in read_csv(pkn / "regulation_chebi_acid_base_relations.csv"):
        eq.join(row["acid"], row["base"])
    # Each collapsed graph edge retains the exact underlying is_a statement.
    parents = collections.defaultdict(dict)
    for node, term in terms.items():
        for parent in term["parents"]:
            source_component, target_component = eq.find(node), eq.find(parent)
            if source_component != target_component:
                parents[source_component].setdefault(target_component, (node, parent))
    exact_features = collections.defaultdict(list)
    for feature in features:
        for cid in filter(None, feature["exact_chebi_ids"].split("|")):
            exact_features[eq.find(cid)].append((feature, cid))
    baseline = collections.defaultdict(list)
    for row in candidates:
        if row["mapping_relation"] == "SwissLipids identifier":
            baseline[(row["feature_id"], eq.find(row["original_chebi_id"]))].append(row["match_type"])
    slm_by_component = collections.defaultdict(list)
    for row in catalog:
        for cid in filter(None, row["chebi_ids"].split("|")):
            slm_by_component[eq.find(cid)].append(row)
    reg_by_source = collections.defaultdict(list)
    for row in regulations:
        reg_by_source[row["source"]].append(row)
    output = []
    for source, records in sorted(reg_by_source.items()):
        start = eq.find(source)
        queue, paths = collections.deque([start]), {start: []}
        while queue:
            current = queue.popleft()
            for target, statement in sorted(parents.get(current, {}).items()):
                if target not in paths:
                    paths[target] = paths[current] + [statement]
                    queue.append(target)
        for ancestor in sorted(set(paths) & exact_features.keys() - {start}):
            for feature, cid in exact_features[ancestor]:
                existing = sorted(set(baseline.get((feature["feature_id"], start), [])))
                term = terms.get(source, dict(name="", formula="", charge=""))
                formula_status = fatty_acid_formula_status(feature["feature_id"], term)
                for (resource, tier), grouped in _group(records).items():
                    slm = slm_by_component.get(start, [])
                    signed = {("Metab__" + source + "_" + (r.get("compartment") or "unknown"),
                               int(r["interaction"]), protein_node(r["target"]))
                              for r in grouped if r["interaction"]}
                    # A graph intersection is an opportunity only. The measurement
                    # is not assigned to this more-specific source chemical.
                    curated_count = len(signed & graph_keys["curated"]) if tier == "curated_mouse_current" else 0
                    exploratory_count = len(signed & graph_keys["exploratory"])
                    output.append(dict(feature=feature["feature_id"], feature_chebi=cid,
                        feature_chebi_name=terms.get(cid, {}).get("name", ""),
                        in_historical_factor1_64=feature["in_factor1"],
                        regulatory_chemical=source, regulatory_name=term["name"],
                        formula=term["formula"], formula_check=formula_status,
                        candidate_relation="specific_regulatory_chemical_descendant_of_measured_composition",
                        identity_confirmed=False, resource=resource, evidence_tier=tier,
                        evidence_rows=len(grouped), unique_signed_records=len(signed),
                        retained_curated_signed_edges=curated_count,
                        retained_exploratory_signed_edges=exploratory_count,
                        target_proteins="|".join(sorted({r["target"] for r in grouped})),
                        target_genes="|".join(sorted({r.get("target_gene", "") for r in grouped} - {""})),
                        signs="|".join(sorted({r["interaction"] for r in grouped} - {""})),
                        evidence_ids="|".join(sorted(r["edge_id"] for r in grouped)),
                        baseline_match_types="|".join(existing),
                        absent_from_frozen_candidate_table=not bool(existing),
                        swisslipids_ids="|".join(sorted({r["swisslipids_id"] for r in slm})),
                        all_related_swisslipids_abbreviations_empty=bool(slm) and all(not r["abbreviation"] for r in slm),
                        is_a_steps=len(paths[ancestor]), is_a_statements_json=json.dumps(paths[ancestor]),
                        equivalence_policy="Only explicit alt IDs and conjugate acid/base groups; is_a is compatibility, never identity."))
    output.sort(key=lambda r: (r["feature"], r["resource"], r["regulatory_chemical"], r["feature_chebi"]))
    assert any(r["feature"] == "FA(22:6)" and r["regulatory_chemical"] == "CHEBI:28125" and
               r["formula_check"] == "matches_reported_FA_composition" and
               r["absent_from_frozen_candidate_table"] for r in output)
    assert any(r["feature"] == "FA(20:4)" and r["regulatory_chemical"] == "CHEBI:15843" for r in output)
    out.mkdir(parents=True, exist_ok=True)
    fields = list(output[0]) if output else []
    write_csv(out / "specific_regulatory_candidates.csv", output, fields)
    summary_rows = []
    for (resource, tier), rows in _group(output).items():
        passing = [r for r in rows if r["formula_check"] != "formula_incompatible"]
        summary_rows.append(dict(resource=resource, evidence_tier=tier,
            candidate_features=len({r["feature"] for r in passing}),
            newly_detected_features=len({r["feature"] for r in passing if r["absent_from_frozen_candidate_table"]}),
            regulatory_chemicals=len({r["regulatory_chemical"] for r in passing}),
            features_with_curated_edge_candidates=len({r["feature"] for r in passing if r["retained_curated_signed_edges"]}),
            features_with_exploratory_edge_candidates=len({r["feature"] for r in passing if r["retained_exploratory_signed_edges"]}),
            historical_factor1_features=len({r["feature"] for r in passing if r["in_historical_factor1_64"] == "True"}),
            candidate_features_list="|".join(sorted({r["feature"] for r in passing})),
            formula_incompatible_rows=sum(r["formula_check"] == "formula_incompatible" for r in rows)))
    write_csv(out / "resource_summary.csv", summary_rows, list(summary_rows[0]))
    passing = [r for r in output if r["formula_check"] != "formula_incompatible"]
    summary = dict(analysis="regulatory_composition_resolution_audit", measured_features=len(features),
        features_with_exact_chebi_terms=sum(bool(r["exact_chebi_ids"]) for r in features),
        candidate_features=len({r["feature"] for r in passing}),
        new_candidate_features=len({r["feature"] for r in passing if r["absent_from_frozen_candidate_table"]}),
        mapping_changed=False, pkn_changed=False, moon_refit=False,
        limitations=["These are conditional identities, not accepted measurements or signed species-specific edges.",
            "is_a may preserve unresolved positional, stereochemical or chain-branching alternatives.",
            "A missing source match is not evidence of biological absence.",
            "Only features already carrying exact reported-resolution ChEBI terms are screened.",
            "MRCLinksDB records are unsigned associations; their adapter's +1 is not used here.",
            "Simple FA candidates additionally require compatible elemental composition; other lipid classes need manual structure review."],
        resources=summary_rows)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    write_csv(out / "dependency_hashes.csv", [dict(path=str(p.relative_to(ROOT)),
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in dependencies], ["path", "sha256"])
    (out / "README.md").write_text("# Regulatory chemical resolution audit\n\n"
        "Run `python3 scripts/moon/ancillary/audit_regulatory_resolution_candidates.py --output NEW_EMPTY_DIRECTORY`. "
        "This audit uses explicit ChEBI subclass relations to detect specific regulatory molecules beneath measured composition-level terms. "
        "It does not change chemical identities, measurements, networks or scores.\n\n"
        "`specific_regulatory_candidates.csv` retains the exact subclass statements, original mapping status, formula checks, "
        "resource tier and evidence IDs. `resource_summary.csv` counts unique feature opportunities, excluding formula-incompatible rows. "
        "Empty SwissLipids abbreviations explain some candidates missed by the earlier abbreviation-only search. "
        "These candidates require explicit identity assumptions or analytical confirmation before use.\n")
    print(json.dumps(summary, indent=2))


def _group(rows):
    result = collections.defaultdict(list)
    for row in rows:
        result[(row["resource"], row["evidence_tier"])].append(row)
    return result


if __name__ == "__main__":
    main()
