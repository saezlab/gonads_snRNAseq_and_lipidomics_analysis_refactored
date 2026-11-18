#!/usr/bin/env python3
"""Build mouse lipid reactions from current Rhea RDF and curated annotations.

Run from the project root with Python >=3.9 (standard library only):
    python3 scripts/pkn/build_rhea_lipid_reactions.py
    python3 scripts/pkn/build_rhea_lipid_reactions.py --self-test

The frozen inputs are reused by default. --download retrieves Rhea inputs using
curl, while the shared mouse UniProt and SwissLipids enzyme files are supplied by
the overall PKN build. The schema retains every participant, coefficient,
annotation and catalyst. It deliberately produces no signed projection.

An unspecified Rhea reaction is NOT reversible. For UN records the compatibility
fields substrates/products mean canonical left/right sides only; direction is
unknown and reversible is null. Directional variants use explicit RDF substrate
and product properties, including reverse variants. Catalysts are alternatives,
never a protein complex inferred by concatenating accessions.
"""

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "support/cosmos_mouse/raw"
OUT = ROOT / "results/causal_modelling/cosmos_mouse_pkn"
RHEA = "http://rdf.rhea-db.org/"
RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
RDFS = "http://www.w3.org/2000/01/rdf-schema#"
CHEBI = "http://purl.obolibrary.org/obo/CHEBI_"
LIPID_QUERY = (
    "PREFIX rh: <http://rdf.rhea-db.org/> "
    "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
    "SELECT DISTINCT ?reaction ?chebi WHERE { "
    "?reaction rdfs:subClassOf rh:Reaction; "
    "rh:side/rh:contains/rh:compound/rh:chebi ?chebi . "
    "?chebi rdfs:subClassOf* <http://purl.obolibrary.org/obo/CHEBI_18059> . }"
)
ALL_LIPID_QUERY = (
    "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> "
    "SELECT DISTINCT ?chebi WHERE { "
    "?chebi rdfs:subClassOf* <http://purl.obolibrary.org/obo/CHEBI_18059> . }"
)
SOURCES = {
    "rhea-release.properties": "https://ftp.expasy.org/databases/rhea/rhea-release.properties",
    "rhea.rdf.gz": "https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz",
    "rhea-directions.tsv": "https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv",
    "rhea2uniprot_sprot.tsv": "https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv",
    "rhea-lipid-reactions.tsv": "https://sparql.rhea-db.org/sparql",
    "rhea-all-lipid-chebi.tsv": "https://sparql.rhea-db.org/sparql",
    "rhea-chebi_pH7_3_mapping.tsv": "https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv",
}


def open_text(path, encoding="utf-8"):
    return gzip.open(path, "rt", encoding=encoding, newline="") if str(path).endswith(".gz") else open(path, encoding=encoding, newline="")


def rows(path, encoding="utf-8"):
    with open_text(path, encoding) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            yield {k.strip().lstrip("?"): v.strip() for k, v in row.items()}


def write_table(path, records, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch_inputs():
    RAW.mkdir(parents=True, exist_ok=True)
    for filename, url in SOURCES.items():
        target = RAW / filename
        command = ["curl", "--fail", "--location", "--silent", "--show-error", "--max-time", "180"]
        if "lipid" in filename:
            query = ALL_LIPID_QUERY if "all-lipid" in filename else LIPID_QUERY
            command += ["--get", "--data-urlencode", "query=" + query,
                        "--data-urlencode", "format=text/tab-separated-values"]
        command += ["--dump-header", str(target) + ".headers.txt", "--output", str(target) + ".partial", url]
        subprocess.run(command, check=True)
        Path(str(target) + ".partial").replace(target)


def load_direction_map(path):
    result = {}
    for row in rows(path):
        master = row["RHEA_ID_MASTER"]
        for suffix, direction in [("MASTER", "UN"), ("LR", "LR"), ("RL", "RL"), ("BI", "BI")]:
            result[row["RHEA_ID_" + suffix]] = (master, direction)
    return result


def merge_description(graph, element):
    subject = element.get("{" + RDF + "}about")
    if subject is None:
        return
    props = graph.setdefault(subject, collections.defaultdict(list))
    for child in element:
        predicate = child.tag[1:].replace("}", "", 1)
        value = child.get("{" + RDF + "}resource", child.text or "")
        if value not in props[predicate]:
            props[predicate].append(value)


def parse_rdf(path, selected_ids, masters):
    """Stream RDF, merging repeated Description records instead of overwriting.

    Rhea encodes side coefficients on typed rh:containsN predicates. Reading only
    rh:contains loses stoichiometry; reading Description tags as unique entities
    loses participants because Rhea splits side metadata across many tags.
    """
    graph = {}
    with gzip.open(path, "rb") as handle:
        events = ET.iterparse(handle, events=("start", "end"))
        _, root = next(events)
        for event, element in events:
            if event != "end" or element.tag != "{" + RDF + "}Description":
                continue
            subject = element.get("{" + RDF + "}about", "")
            local = subject.removeprefix(RHEA)
            keep = local in selected_ids or local in masters or local.startswith("Compound_") or local.startswith("contains")
            if local.endswith(("_L", "_R")):
                keep = local.rsplit("_", 1)[0] in masters
            if local.startswith("Participant_"):
                keep = local.split("_", 2)[1] in masters
            if keep:
                merge_description(graph, element)
            element.clear()
            root.clear()
    return graph


def values(graph, node, property_name, namespace=RHEA):
    return graph.get(node, {}).get(namespace + property_name, [])


def one(graph, node, property_name, default=None, namespace=RHEA):
    result = values(graph, node, property_name, namespace)
    if len(result) > 1:
        raise ValueError(f"Expected one {property_name} for {node}, got {result}")
    return result[0] if result else default


def coefficient(value):
    if re.fullmatch(r"[0-9]+", value):
        return int(value)
    if re.fullmatch(r"[0-9]+\.[0-9]+", value):
        return float(value)
    return value  # Preserve symbolic polymer coefficients, never invent a value.


def molecule(graph, node, lipid_chebi):
    chebi = one(graph, node, "chebi")
    ident = chebi.replace(CHEBI, "CHEBI:") if chebi else "RHEA_COMPOUND:" + node.rsplit("_", 1)[-1]
    classes = [c.removeprefix(RHEA) for c in values(graph, node, "subClassOf", RDFS) if c.startswith(RHEA)]
    reactive = []
    for part in values(graph, node, "reactivePart"):
        reactive.append({"rhea_compound": part, "chebi": one(graph, part, "chebi"),
                         "name": one(graph, part, "name"), "position": one(graph, part, "position")})
    return {
        "id": ident, "name": one(graph, node, "name"), "rhea_compound": node,
        "chemical_type": classes, "formula": one(graph, node, "formula"),
        "charge": one(graph, node, "charge"), "is_lipid": ident in lipid_chebi,
        "reactive_parts": reactive, "underlying_chebi": one(graph, node, "underlyingChebi"),
        "polymerization_index": one(graph, node, "polymerizationIndex"),
    }


def side_participants(graph, side, lipid_chebi):
    participants = values(graph, side, "contains")
    if not participants:
        raise ValueError("Reaction side has no participants: " + side)
    typed = collections.defaultdict(set)
    for predicate, targets in graph[side].items():
        if predicate.startswith(RHEA + "contains") and predicate != RHEA + "contains":
            coeff = one(graph, predicate, "coefficient")
            if coeff is None:
                raise ValueError("No coefficient definition for " + predicate)
            for target in targets:
                typed[target].add(coeff)
    if set(typed) != set(participants):
        raise ValueError("Typed/untyped participant sets differ on " + side)
    result = []
    for participant in participants:
        if len(typed[participant]) != 1:
            raise ValueError("Ambiguous coefficient: " + participant)
        compound = one(graph, participant, "compound")
        if not compound or compound not in graph:
            raise ValueError("Unresolved compound: " + participant)
        info = molecule(graph, compound, lipid_chebi)
        info.update({"compartment": "unknown", "relative_location": one(graph, participant, "location"),
                     "stoichiometry": coefficient(next(iter(typed[participant]))),
                     "participant_id": participant, "side_id": side})
        result.append(info)
    return result


def reaction_record(graph, rid, master, direction, annotations, lipid_chebi, source_version="unspecified"):
    uri = RHEA + rid
    if one(graph, uri, "status") != RHEA + "Approved":
        raise ValueError("Reaction not approved or missing: " + rid)
    if direction in ("LR", "RL"):
        left = one(graph, uri, "substrates")
        right = one(graph, uri, "products")
        basis = "explicit Rhea substrates/products; annotated directional accession"
    else:
        sides = values(graph, uri, "side" if direction == "UN" else "substratesOrProducts")
        sides = sorted(sides, key=lambda s: int(one(graph, s, "curatedOrder")))
        if len(sides) != 2:
            raise ValueError("Expected exactly two reaction sides: " + rid)
        left, right = sides
        basis = "canonical left/right only; physiological direction unspecified" if direction == "UN" else "annotated bidirectional accession"
    if not left or not right:
        raise ValueError("Unresolved substrate/product sides: " + rid)
    substrates, products = side_participants(graph, left, lipid_chebi), side_participants(graph, right, lipid_chebi)
    return {
        "reaction_id": "RHEA:" + rid, "master_reaction_id": "RHEA:" + master,
        "resource": "Rhea", "source_version": source_version, "organism": 10090, "direction": direction,
        "direction_basis": basis, "reversible": None if direction == "UN" else direction == "BI",
        "reversibility_scope": "annotation-specific traversal permission; LR/RL false does not establish thermodynamic irreversibility",
        "substrates": substrates, "products": products,
        "catalysts": sorted({a["uniprot"] for a in annotations}),
        "catalyst_semantics": "alternative annotated catalysts; no complex inferred",
        "equation": one(graph, uri, "equation"),
        "is_transport": one(graph, uri, "isTransport") == "true",
        "chemically_balanced": one(graph, uri, "isChemicallyBalanced") == "true",
        "references": values(graph, uri, "citation"),
        "reference_scope": "Rhea reaction-level references; catalyst evidence retained separately",
        "provenance": annotations,
        "selection_basis": sorted({a["selection_basis"] for a in annotations}),
    }


def main(args):
    if args.download:
        fetch_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    release = dict(line.strip().split("=", 1) for line in (RAW / "rhea-release.properties").read_text().splitlines() if "=" in line)
    source_version = f"{release['rhea.release.number']} ({release['rhea.release.date']})"
    direction_map = load_direction_map(RAW / "rhea-directions.tsv")
    lipid_rows = list(rows(RAW / "rhea-lipid-reactions.tsv"))
    lipid_participants = {r["chebi"].replace(CHEBI, "CHEBI:") for r in lipid_rows}
    all_lipids = {r["chebi"].replace(CHEBI, "CHEBI:") for r in rows(RAW / "rhea-all-lipid-chebi.tsv") if r["chebi"].startswith(CHEBI)}
    for filename, ids in [("rhea_lipid_chebi.tsv", lipid_participants), ("rhea_all_lipid_chebi.tsv", all_lipids)]:
        write_table(OUT / filename, ({"id": x, "classification": "CHEBI:18059 descendant",
                                     "resource": "Rhea synchronized ChEBI hierarchy"} for x in sorted(ids)),
                    ["id", "classification", "resource"])
    ph_mapping = list(rows(RAW / "rhea-chebi_pH7_3_mapping.tsv"))
    write_table(OUT / "rhea_chebi_ph73.tsv", ({"source": "CHEBI:" + row["CHEBI"],
                "target": "CHEBI:" + row["CHEBI_PH7_3"], "origin": row["ORIGIN"],
                "relation": "major microspecies at pH 7.3"} for row in ph_mapping),
                ["source", "target", "origin", "relation"])
    lipid_identity_evidence = collections.defaultdict(list)
    for ident in all_lipids:
        lipid_identity_evidence[ident].append({"basis": "direct ChEBI subclass descendant", "origin": "ontology", "lipid_id": ident})
    for row in ph_mapping:
        source, target = "CHEBI:" + row["CHEBI"], "CHEBI:" + row["CHEBI_PH7_3"]
        if source in all_lipids and source != target:
            lipid_identity_evidence[target].append({"basis": "pH 7.3 microspecies of lipid descendant", "origin": row["ORIGIN"], "lipid_id": source})
    write_table(OUT / "rhea_lipid_identity_classification.tsv", ({"id": ident, **evidence}
                for ident in sorted(lipid_identity_evidence) for evidence in lipid_identity_evidence[ident]),
                ["id", "basis", "origin", "lipid_id"])
    lipid_identities = set(lipid_identity_evidence)
    mouse = {r["Entry"]: r for r in rows(args.mouse_uniprot)}
    annotations = collections.defaultdict(list)
    rejected = []
    mouse_swiss_rows = 0
    for index, row in enumerate(rows(args.swisslipids_enzymes, "latin-1"), start=2):
        if row["Protein taxon"] != "10090":
            continue
        mouse_swiss_rows += 1
        rid, accession = row["Rhea ID"], row["UniProtKB AC(s)"]
        # The SwissLipids enzyme table currently has one accession per row. Fail
        # closed on an unfamiliar format rather than treat a list as a complex.
        if not re.fullmatch(r"[A-Z0-9]+(?:-[0-9]+)?", accession):
            rejected.append({"resource": "SwissLipids", "row": index, "reaction": rid,
                             "reason": "ambiguous multi-accession field", "value": accession})
            continue
        if rid not in direction_map:
            rejected.append({"resource": "SwissLipids", "row": index, "reaction": rid,
                             "reason": "reaction missing from current Rhea quartet map", "value": accession})
            continue
        annotations[rid].append({
            "resource": "SwissLipids enzymes", "uniprot": accession, "organism": 10090,
            "gene_symbol": row["Gene name"], "swisslipids_enzyme_id": row["SwissLipids ID"],
            "evidence_tag_ids": [x.strip() for x in row["Evidence tag ID"].split("|")],
            "annotated_rhea_id": "RHEA:" + rid, "annotation_equation": html.unescape(row["Reaction text"]),
            "source_row": index, "source_file": str(args.swisslipids_enzymes.relative_to(ROOT)),
            "selection_basis": "SwissLipids mouse enzyme annotation",
            "source_url": "https://www.swisslipids.org/api/file.php?cas=download_files&file=enzymes.tsv",
        })
    for index, row in enumerate(rows(RAW / "rhea2uniprot_sprot.tsv"), start=2):
        if row["ID"] not in mouse:
            continue
        rid = row["RHEA_ID"]
        if direction_map.get(rid) != (row["MASTER_ID"], row["DIRECTION"]):
            raise ValueError("Rhea crosswalk direction mismatch: " + str(row))
        annotations[rid].append({
            "resource": "Rhea-UniProtKB/Swiss-Prot", "uniprot": row["ID"], "organism": 10090,
            "gene_symbol": mouse[row["ID"]]["Gene Names (primary)"],
            "annotated_rhea_id": "RHEA:" + rid, "annotated_direction": row["DIRECTION"],
            "selection_basis": "mouse curated UniProt annotation and ChEBI lipid participant (direct or pH7.3-related)",
            "source_file": "support/cosmos_mouse/raw/rhea2uniprot_sprot.tsv", "source_row": index,
            "source_url": SOURCES["rhea2uniprot_sprot.tsv"],
        })
    masters = {direction_map[rid][0] for rid in annotations}
    print(f"Selected {len(annotations)} annotated reaction variants / {len(masters)} masters; parsing RDF", flush=True)
    graph = parse_rdf(RAW / "rhea.rdf.gz", set(annotations), masters)
    reactions = []
    for rid, evidence in sorted(annotations.items(), key=lambda pair: int(pair[0])):
        if one(graph, RHEA + rid, "status") != RHEA + "Approved":
            rejected.append({"resource": "Rhea", "reaction": rid, "reason": "not currently approved",
                             "status": one(graph, RHEA + rid, "status"), "annotations": evidence})
            continue
        master, direction = direction_map[rid]
        record = reaction_record(graph, rid, master, direction, evidence, lipid_identities, source_version)
        for participant in record["substrates"] + record["products"]:
            participant["lipid_classification_evidence"] = lipid_identity_evidence.get(participant["id"], [])
        if not any(p["is_lipid"] for p in record["substrates"] + record["products"]) and not any(a["resource"] == "SwissLipids enzymes" for a in evidence):
            continue
        reactions.append(record)
    with (OUT / "rhea_reactions.jsonl").open("w") as handle:
        for reaction in reactions:
            handle.write(json.dumps(reaction, ensure_ascii=False, sort_keys=True) + "\n")
    molecules = {}
    for reaction in reactions:
        for participant in reaction["substrates"] + reaction["products"]:
            molecules[participant["id"]] = {key: participant[key] for key in ["id", "name", "formula", "charge", "is_lipid"]}
            molecules[participant["id"]]["chemical_type"] = ";".join(participant["chemical_type"])
    write_table(OUT / "rhea_molecules.tsv", [molecules[k] for k in sorted(molecules)],
                ["id", "name", "formula", "charge", "is_lipid", "chemical_type"])
    (OUT / "rhea_rejected_annotations.json").write_text(json.dumps(rejected, indent=2, ensure_ascii=False) + "\n")
    summary = {
        "reaction_variants": len(reactions), "reaction_masters": len({r["master_reaction_id"] for r in reactions}),
        "direction_counts": dict(collections.Counter(r["direction"] for r in reactions)),
        "catalysts": len({c for r in reactions for c in r["catalysts"]}),
        "unique_participants": len(molecules), "lipid_participants": sum(m["is_lipid"] for m in molecules.values()),
        "participant_records": sum(len(r["substrates"]) + len(r["products"]) for r in reactions),
        "unequal_side_participant_counts": sum(len(r["substrates"]) != len(r["products"]) for r in reactions),
        "symbolic_stoichiometry_records": sum(isinstance(p["stoichiometry"], str) for r in reactions for p in r["substrates"] + r["products"]),
        "non_chebi_participant_ids": sorted(i for i in molecules if not i.startswith("CHEBI:")),
        "swisslipids_mouse_rows": mouse_swiss_rows, "rejected_annotations": len(rejected),
        "mouse_annotated_variants_examined_before_lipid_filter": len(annotations),
        "swisslipids_enzyme_context_variants_without_classified_lipid": sum(not any(p["is_lipid"] for p in r["substrates"] + r["products"]) for r in reactions),
        "all_chebi_lipid_descendants": len(all_lipids), "rhea_lipid_participant_descendants": len(lipid_participants),
        "lipid_identities_including_ph73": len(lipid_identities),
        "ph73_mapping_origins": dict(collections.Counter(row["ORIGIN"] for row in ph_mapping)),
        "direction_policy": "UN has reversible=null; no signed projection. LR/RL use explicit RDF sides, and reversible=false means one annotated traversal direction, not thermodynamic irreversibility. BI requires annotated BI accession.",
        "species_policy": "Mouse SwissLipids taxonomy or exact mouse UniProt accession; Rhea chemistry is species independent.",
        "resolution_policy": "Chemical class participants remain classes; no unannotated species-level substitution. Generic compounds retain RHEA_COMPOUND IDs and reactive-part annotations.",
    }
    (OUT / "rhea_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    manifest = {
        "built_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "source_version": "Rhea release " + source_version,
        "sources": [{"file": str((RAW / filename).relative_to(ROOT)), "url": url,
                     "sha256_file": digest(RAW / filename),
                     **({"sparql_query": ALL_LIPID_QUERY if "all-lipid" in filename else LIPID_QUERY} if "lipid" in filename else {})}
                    for filename, url in SOURCES.items()],
        "shared_inputs": [{"file": str(p.relative_to(ROOT)), "sha256_file": digest(p)}
                          for p in [args.swisslipids_enzymes, args.mouse_uniprot]],
        "outputs": ["rhea_reactions.jsonl", "rhea_molecules.tsv", "rhea_lipid_chebi.tsv", "rhea_all_lipid_chebi.tsv", "rhea_lipid_identity_classification.tsv", "rhea_chebi_ph73.tsv", "rhea_rejected_annotations.json", "rhea_summary.json"],
        "rhea_direction_documentation": "https://www.rhea-db.org/help/reaction-side-direction",
        "rhea_download_documentation": "https://www.rhea-db.org/help/download",
    }
    (OUT / "rhea_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def self_test():
    """Exercise unequal sides, split RDF records, coefficients and UN semantics."""
    graph = {}
    def put(node, prop, val):
        graph.setdefault(RHEA + node, collections.defaultdict(list))[RHEA + prop].append(val)
    for n in range(1, 6):
        put(f"Compound_{n}", "chebi", CHEBI + str(n))
        put(f"Compound_{n}", "name", f"molecule {n}")
        put(f"Participant_10000_compound_{n}", "compound", RHEA + f"Compound_{n}")
    for side, indexes in [("L", (1, 2)), ("R", (3, 4, 5))]:
        put("10000_" + side, "curatedOrder", "1" if side == "L" else "2")
        for n in indexes:
            put("10000_" + side, "contains", RHEA + f"Participant_10000_compound_{n}")
            put("10000_" + side, "contains2" if n == 1 else "contains1", RHEA + f"Participant_10000_compound_{n}")
    put("contains1", "coefficient", "1")
    put("contains2", "coefficient", "2")
    for rid in ("10000", "10001", "10002", "10003"):
        put(rid, "status", RHEA + "Approved")
    put("10000", "side", RHEA + "10000_L"); put("10000", "side", RHEA + "10000_R")
    put("10001", "substrates", RHEA + "10000_L"); put("10001", "products", RHEA + "10000_R")
    put("10002", "products", RHEA + "10000_L"); put("10002", "substrates", RHEA + "10000_R")
    put("10003", "substratesOrProducts", RHEA + "10000_L"); put("10003", "substratesOrProducts", RHEA + "10000_R")
    evidence = [{"uniprot": "QTEST", "selection_basis": "test"}]
    un = reaction_record(graph, "10000", "10000", "UN", evidence, set())
    lr = reaction_record(graph, "10001", "10000", "LR", evidence, set())
    rl = reaction_record(graph, "10002", "10000", "RL", evidence, set())
    bi = reaction_record(graph, "10003", "10000", "BI", evidence, set())
    assert un["reversible"] is None and bi["reversible"] is True
    assert lr["reversible"] is False and rl["reversible"] is False
    assert len(lr["substrates"]) == 2 and len(lr["products"]) == 3
    assert len(rl["substrates"]) == 3 and len(rl["products"]) == 2
    assert lr["substrates"][0]["stoichiometry"] == 2
    # Rhea repeats Description tags for the same subject. A parser must merge.
    a = ET.fromstring(f'<rdf:Description xmlns:rdf="{RDF}" xmlns:rh="{RHEA}" rdf:about="{RHEA}side"><rh:contains rdf:resource="A"/></rdf:Description>')
    b = ET.fromstring(f'<rdf:Description xmlns:rdf="{RDF}" xmlns:rh="{RHEA}" rdf:about="{RHEA}side"><rh:contains rdf:resource="B"/></rdf:Description>')
    merged = {}; merge_description(merged, a); merge_description(merged, b)
    assert values(merged, RHEA + "side", "contains") == ["A", "B"]
    assert coefficient("2n") == "2n"
    print("Rhea parser self-tests passed: lossless unequal sides, repeated RDF subjects, coefficients, LR/RL/BI/UN.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Refresh Rhea sources (requires public network access)")
    parser.add_argument("--self-test", action="store_true", help="Run parser semantic regression checks only")
    parser.add_argument("--mouse-uniprot", type=Path, default=RAW / "uniprot_mouse.tsv.gz")
    parser.add_argument("--swisslipids-enzymes", type=Path, default=ROOT / "results/causal_modelling/cosmos_pkn_audit/raw/swisslipids-enzymes.tsv.gz")
    options = parser.parse_args()
    if options.self_test:
        self_test()
    else:
        main(options)
