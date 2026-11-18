#!/usr/bin/env python3
"""Replace orientation-unknown cached STITCH rows with raw directional evidence.

The frozen OmniPath-metabo mouse cache normalised STITCH endpoints to
chemical -> protein after parsing.  It did not retain the original A/B order
or ``a_is_acting`` flag, so it cannot distinguish a true chemical-to-protein
action from a reversed protein-to-chemical action.  This builder leaves the
published PKN untouched and creates a separate exploratory overlay from the
raw STITCH mouse action archive.

Admission is deliberately source-level and conservative:

* raw A must be a STITCH chemical (CIDm/CIDs) and raw B a mouse ENSMUSP ID;
* ``a_is_acting`` must be ``t``;
* ``action`` must be exactly activation or inhibition, from which the sign is
  derived (``mode`` is provenance only);
* action score must be at least 700 and detailed-link text mining below 700;
* chemical and protein mappings must each resolve uniquely in frozen maps.

It is a directionality correction, not a claim that STITCH v5.0 is current or
that every retained relationship is direct mouse experimental evidence.  The
raw source was last modified in 2016 and remains an exploratory sensitivity
layer only.
"""

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "support" / "cosmos_mouse" / "raw"
BASE = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn"
DEFAULT_OUT = ROOT / "results" / "causal_modelling" / "cosmos_mouse_pkn_stitch_directional"

ACTIONS = RAW / "stitch_actions_10090_v5.0.tsv.gz"
LINKS = RAW / "stitch_links_detailed_10090_v5.0.tsv.gz"
STRING_ALIASES = RAW / "string_aliases_10090_v12.0.txt.gz"
UNICHEM = RAW / "unichem_chebi_to_pubchem_20260912.tsv.gz"
UNIPROT = RAW / "uniprot_mouse.tsv.gz"
CHEBI = RAW / "regulation_chebi.obo.gz"
CHEBI_ALIASES = BASE / "regulation_chebi_aliases.csv"

ACTION_URL = "http://stitch-db.org/download/actions.v5.0/10090.actions.v5.0.tsv.gz"
LINK_URL = "http://stitch-db.org/download/protein_chemical.links.detailed.v5.0/10090.protein_chemical.links.detailed.v5.0.tsv.gz"
STRING_ALIAS_URL = "https://stringdb-static.org/download/protein.aliases.v12.0/10090.protein.aliases.v12.0.txt.gz"
UNICHEM_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id7/src7src22.txt.gz"

REGULATION_FIELDS = [
    "edge_id", "source", "target", "interaction", "resource", "evidence_tier",
    "mechanism", "organism", "source_organism", "references", "compartment",
    "source_version", "source_type", "target_type", "source_name", "target_gene",
    "target_class", "ligand_type", "endogenous", "source_record_id", "source_file",
    "source_sha256_uncompressed", "mapping_method", "mapping_count", "mapping_group",
    "source_pubchem", "original_action", "original_type", "sign_basis", "affinity_units",
    "affinity_value", "assay_description", "sign_conflict", "quality_notes",
]
RAW_PROVENANCE_FIELDS = [
    "raw_item_id_a", "raw_item_id_b", "raw_a_is_acting", "raw_action_score",
    "raw_action_line", "raw_link_experimental", "raw_link_prediction",
    "raw_link_database", "raw_link_textmining", "raw_link_combined_score",
    "source_sha256_compressed", "detailed_link_sha256_compressed",
]
STITCH_EVIDENCE_FIELDS = REGULATION_FIELDS + RAW_PROVENANCE_FIELDS

CHEMICAL = re.compile(r"^CID[ms](\d+)$")
MOUSE_ENSMUSP = re.compile(r"^10090\.ENSMUSP\d+$")
GENE_SYMBOL = re.compile(r"^[A-Za-z0-9_.-]+$")
VARIANT_ORDER = {"core": 0, "lipid_curated": 1, "lipid_rules": 2, "exploratory": 3}


def compressed_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def uncompressed_sha256(path):
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonl(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def raw_pubchem_cid(identifier):
    match = CHEMICAL.fullmatch(str(identifier))
    if not match:
        return None
    return str(int(match.group(1)))


def preliminary_action_reason(row, minimum_action_score):
    """Return a source-level rejection reason before identifier mapping."""
    if raw_pubchem_cid(row.get("item_id_a", "")) is None:
        return "raw_actor_is_not_STITCH_chemical"
    if not MOUSE_ENSMUSP.fullmatch(str(row.get("item_id_b", ""))):
        return "raw_receiver_is_not_mouse_ENSMUSP_protein"
    if str(row.get("a_is_acting", "")).lower() != "t":
        return "raw_chemical_is_not_acting_partner"
    if str(row.get("action", "")).lower() not in {"activation", "inhibition"}:
        return "raw_action_is_not_explicitly_signed"
    try:
        score = int(row["score"])
    except (KeyError, TypeError, ValueError):
        return "invalid_action_score"
    if score < minimum_action_score:
        return "action_score_below_threshold"
    return None


def detailed_link_reason(row, maximum_textmining):
    """Apply the pre-specified detailed-link text-mining filter without inferring sign."""
    try:
        textmining = int(row["textmining"])
    except (KeyError, TypeError, ValueError):
        return "invalid_detailed_link_textmining_score"
    if textmining >= maximum_textmining:
        return "detailed_link_textmining_at_or_above_threshold"
    return None


def current_chebi_terms(path):
    """Return current ChEBI names while excluding obsolete terms."""
    current, names = set(), {}
    identifier = name = None
    obsolete = False

    def flush():
        if identifier and not obsolete:
            current.add(identifier)
            names[identifier] = name or identifier

    with gzip.open(path, "rt") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("["):
                flush()
                identifier = name = None
                obsolete = False
            elif line.startswith("id: CHEBI:"):
                identifier = line.split(" ", 1)[1]
            elif line.startswith("name: "):
                name = line[6:]
            elif line == "is_obsolete: true":
                obsolete = True
        flush()
    return current, names


def chebi_primary_map(path):
    aliases = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            aliases[row["alias"]] = row["primary_chebi"]
    return aliases


def canonical_chebi(identifier, aliases):
    seen = set()
    while identifier in aliases and identifier not in seen:
        seen.add(identifier)
        identifier = aliases[identifier]
    return identifier


def pubchem_to_current_chebi(path, aliases, current):
    mapped = collections.defaultdict(set)
    with gzip.open(path, "rt") as handle:
        header = next(handle).rstrip("\n")
        if header != "From src:'7'\tTo src:'22'":
            raise ValueError("Unexpected UniChem ChEBI/PubChem mapping header: " + header)
        for line in handle:
            chebi, pubchem = line.rstrip("\n").split("\t")
            chebi = canonical_chebi(chebi, aliases)
            if chebi in current:
                mapped[pubchem].add(chebi)
    return mapped


def frozen_mouse_proteins(uniprot_path, alias_path):
    genes, accessions = {}, set()
    with gzip.open(uniprot_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["Entry"]
            accessions.add(accession)
            gene = row["Gene Names (primary)"].strip()
            if GENE_SYMBOL.fullmatch(gene):
                genes[accession] = gene
    aliases = collections.defaultdict(set)
    with gzip.open(alias_path, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source"] == "UniProt_AC" and row["alias"] in accessions:
                aliases[row["#string_protein_id"]].add(row["alias"])
    return genes, aliases


def candidate_from_action(row, chebi_by_cid, proteins, genes, source_hashes, line_number):
    """Map a source-level accepted action or return a mapping rejection reason."""
    pubchem = raw_pubchem_cid(row["item_id_a"])
    chebis = sorted(chebi_by_cid.get(pubchem, set()))
    if not chebis:
        return None, "chemical_has_no_current_unique_CHEBI_mapping"
    if len(chebis) != 1:
        return None, "chemical_has_ambiguous_current_CHEBI_mapping"
    uniprots = sorted(proteins.get(row["item_id_b"], set()))
    if not uniprots:
        return None, "mouse_protein_has_no_frozen_UniProt_mapping"
    if len(uniprots) != 1:
        return None, "mouse_protein_has_ambiguous_frozen_UniProt_mapping"
    accession = uniprots[0]
    gene = genes.get(accession, "")
    target_node = "Gene__" + gene if gene else "UniProt__" + accession
    sign = 1 if row["action"].lower() == "activation" else -1
    source = chebis[0]
    key = "|".join([row["item_id_a"], row["item_id_b"], row["action"], str(line_number), source, accession])
    edge_id = "reg_stitch_v5_" + hashlib.sha256(key.encode()).hexdigest()[:20]
    record = {
        "edge_id": edge_id,
        "source": source,
        "target": accession,
        "interaction": sign,
        "resource": "STITCH_directional_v5",
        "evidence_tier": "fresh_stitch_v5_chemical_actor_directional_exploratory",
        "mechanism": "chemical_protein_activity_modulation",
        "organism": "10090",
        "source_organism": "",
        "references": "",
        "compartment": "unknown",
        "source_version": "STITCH 5.0 mouse actions/detailed links; source files last modified 2016-08-16; retrieved 2026-09-15",
        "source_type": "small_molecule",
        "target_type": "protein",
        "source_name": "",
        "target_gene": gene,
        "target_class": "STITCH_chemical_protein_action",
        "ligand_type": "unclassified_chemical",
        "endogenous": "",
        "source_record_id": "STITCH:5.0:actions:line:" + str(line_number),
        "source_file": "support/cosmos_mouse/raw/stitch_actions_10090_v5.0.tsv.gz;support/cosmos_mouse/raw/stitch_links_detailed_10090_v5.0.tsv.gz",
        "source_sha256_uncompressed": source_hashes["actions_uncompressed"],
        "mapping_method": "raw_STITCH_CID_to_UniChem_current_ChEBI_unique;STRING_v12_UniProt_AC_exact;frozen_mouse_UniProt_validation",
        "mapping_count": "1",
        "mapping_group": "raw_STITCH:" + row["item_id_a"] + ";UniChem:7_to_22",
        "source_pubchem": pubchem,
        "original_action": row["action"].lower(),
        "original_type": row.get("mode", ""),
        "sign_basis": "raw_STITCH_action_with_chemical_actor",
        "affinity_units": "",
        "affinity_value": "",
        "assay_description": "",
        "sign_conflict": "",
        "quality_notes": "Raw STITCH partner A is a CIDm/CIDs chemical and partner B a mouse ENSMUSP protein; a_is_acting=t. Sign is derived from raw action, not mode. Detailed-link text mining is below 700. The current STRING v12 alias snapshot was used for exact ENSMUSP-to-UniProt identity because a public v5 alias snapshot was unavailable. This directionally auditable STITCH 5.0 evidence remains exploratory because the source release is historical and may include predictions or transferred evidence.",
        "raw_item_id_a": row["item_id_a"],
        "raw_item_id_b": row["item_id_b"],
        "raw_a_is_acting": row["a_is_acting"],
        "raw_action_score": row["score"],
        "raw_action_line": line_number,
        "raw_link_experimental": "",
        "raw_link_prediction": "",
        "raw_link_database": "",
        "raw_link_textmining": "",
        "raw_link_combined_score": "",
        "source_sha256_compressed": source_hashes["actions_compressed"],
        "detailed_link_sha256_compressed": source_hashes["links_compressed"],
        "_graph_source": "Metab__" + source + "_unknown",
        "_graph_target": target_node,
        "_target_accession": accession,
    }
    return record, None


def add_pair(pair_signs, source, target, sign):
    pair_signs[(source, target)].add(int(sign))


def project_signed_pairs(evidence):
    """Return unique signed triples and opposing-sign pairs; never majority vote."""
    pair_signs = collections.defaultdict(set)
    for edge in evidence:
        add_pair(pair_signs, edge["source"], edge["target"], edge["interaction"])
    graph = [(source, next(iter(signs)), target)
             for (source, target), signs in pair_signs.items() if len(signs) == 1]
    conflicts = [(source, target) for (source, target), signs in pair_signs.items() if len(signs) > 1]
    return sorted(graph), sorted(conflicts), pair_signs


def validate_inputs(paths):
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Required source(s) missing:\n" + "\n".join(missing))


def build(output, minimum_action_score=700, maximum_textmining=700):
    validate_inputs([ACTIONS, LINKS, STRING_ALIASES, UNICHEM, UNIPROT, CHEBI, CHEBI_ALIASES,
                     BASE / "edge_evidence.jsonl.gz", BASE / "nodes.jsonl", BASE / "regulation_edges.csv",
                     BASE / "pkn_lipid_rules.csv"])
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Refusing to overwrite an existing output directory: " + str(output))
    output.mkdir(parents=True, exist_ok=True)

    print("[STITCH] Hashing frozen source files")
    source_hashes = {
        "actions_compressed": compressed_sha256(ACTIONS),
        "links_compressed": compressed_sha256(LINKS),
        "string_aliases_compressed": compressed_sha256(STRING_ALIASES),
        "unichem_compressed": compressed_sha256(UNICHEM),
        "uniprot_compressed": compressed_sha256(UNIPROT),
        "chebi_compressed": compressed_sha256(CHEBI),
        "actions_uncompressed": uncompressed_sha256(ACTIONS),
        "links_uncompressed": uncompressed_sha256(LINKS),
    }
    print("[STITCH] Loading frozen identifier crosswalks")
    current, chebi_names = current_chebi_terms(CHEBI)
    aliases = chebi_primary_map(CHEBI_ALIASES)
    chebi_by_cid = pubchem_to_current_chebi(UNICHEM, aliases, current)
    genes, proteins = frozen_mouse_proteins(UNIPROT, STRING_ALIASES)
    print("[STITCH] Current ChEBI terms:", len(current), "; mapped PubChem CIDs:", len(chebi_by_cid),
          "; uniquely mapped mouse STRING proteins:", sum(len(x) == 1 for x in proteins.values()))

    action_counts = collections.Counter()
    mapping_inputs = {"chemical": set(), "protein": set()}
    pending = collections.defaultdict(list)
    with gzip.open(ACTIONS, "rt", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
            action_counts["raw_action_records"] += 1
            reason = preliminary_action_reason(row, minimum_action_score)
            if reason:
                action_counts[reason] += 1
                continue
            action_counts["source_level_directional_signed_actions"] += 1
            mapping_inputs["chemical"].add(raw_pubchem_cid(row["item_id_a"]))
            mapping_inputs["protein"].add(row["item_id_b"])
            record, reason = candidate_from_action(row, chebi_by_cid, proteins, genes, source_hashes, line_number)
            if reason:
                action_counts[reason] += 1
                continue
            pending[(row["item_id_a"], row["item_id_b"])].append(record)
            action_counts["mapped_action_records_before_detailed_link_filter"] += 1

    mapping_review = []
    for cid in sorted(mapping_inputs["chemical"], key=int):
        candidates = sorted(chebi_by_cid.get(cid, set()))
        mapping_review.append({
            "entity_type": "raw_STITCH_chemical", "input_identifier": "PubChem:" + cid,
            "candidate_mappings": ";".join(candidates), "mapping_count": len(candidates),
            "mapping_method": "UniChem whole-source mapping 7 (ChEBI) -> 22 (PubChem), current ChEBI only",
            "admitted": len(candidates) == 1,
            "reason": "unique_current_CHEBI" if len(candidates) == 1 else
                      "no_current_CHEBI" if not candidates else "ambiguous_current_CHEBI",
        })
    for string_id in sorted(mapping_inputs["protein"]):
        candidates = sorted(proteins.get(string_id, set()))
        mapping_review.append({
            "entity_type": "mouse_STRING_protein", "input_identifier": string_id,
            "candidate_mappings": ";".join(candidates), "mapping_count": len(candidates),
            "mapping_method": "STRING v12 UniProt_AC alias restricted to frozen mouse UniProt export",
            "admitted": len(candidates) == 1,
            "reason": "unique_mouse_UniProt" if len(candidates) == 1 else
                      "no_mouse_UniProt" if not candidates else "ambiguous_mouse_UniProt",
        })

    print("[STITCH] Joining detailed link evidence for", len(pending), "raw chemical/protein pairs")
    retained = []
    with gzip.open(LINKS, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = (row["chemical"], row["protein"])
            records = pending.pop(key, None)
            if not records:
                continue
            reason = detailed_link_reason(row, maximum_textmining)
            if reason:
                action_counts[reason] += len(records)
                continue
            for record in records:
                record["raw_link_experimental"] = row["experimental"]
                record["raw_link_prediction"] = row["prediction"]
                record["raw_link_database"] = row["database"]
                record["raw_link_textmining"] = row["textmining"]
                record["raw_link_combined_score"] = row["combined_score"]
                retained.append(record)
    missing_link_records = sum(len(records) for records in pending.values())
    action_counts["detailed_link_missing"] += missing_link_records
    action_counts["retained_directional_STITCH_evidence_rows"] = len(retained)
    if not retained:
        raise ValueError("No directional STITCH rows survived; refusing to write an empty overlay")
    for record in retained:
        record["source_name"] = chebi_names.get(record["source"], record["source"])

    # Build graph conflicts from all evidence, after replacing every cache-STITCH row.
    pair_signs = collections.defaultdict(set)
    cache_stitch_removed = 0
    base_evidence_rows = 0
    for edge in jsonl(BASE / "edge_evidence.jsonl.gz"):
        if edge.get("resource") == "STITCH":
            cache_stitch_removed += 1
            continue
        base_evidence_rows += 1
        if VARIANT_ORDER[edge["first_variant"]] <= VARIANT_ORDER["exploratory"]:
            add_pair(pair_signs, edge["source"], edge["target"], edge["interaction"])
    for record in retained:
        add_pair(pair_signs, record["_graph_source"], record["_graph_target"], record["interaction"])
    graph = [(source, next(iter(signs)), target)
             for (source, target), signs in pair_signs.items() if len(signs) == 1]
    conflicts = [(source, target) for (source, target), signs in pair_signs.items() if len(signs) > 1]
    graph.sort()
    conflicts.sort()
    conflict_pairs = set(conflicts)
    for record in retained:
        record["sign_conflict"] = (record["_graph_source"], record["_graph_target"]) in conflict_pairs
    fresh_projected_triples = {
        (record["_graph_source"], int(record["interaction"]), record["_graph_target"])
        for record in retained
    }
    final_graph_triples = set(graph)

    # Full evidence ledger: copy all non-STITCH cache evidence verbatim and append raw-derived rows.
    with gzip.open(output / "edge_evidence.jsonl.gz", "wt") as handle:
        for edge in jsonl(BASE / "edge_evidence.jsonl.gz"):
            if edge.get("resource") != "STITCH":
                handle.write(json.dumps(edge, sort_keys=True) + "\n")
        for record in sorted(retained, key=lambda x: x["edge_id"]):
            ledger = {key: value for key, value in record.items() if not key.startswith("_")}
            ledger.update({
                "source": record["_graph_source"], "target": record["_graph_target"],
                "interaction": int(record["interaction"]), "first_variant": "exploratory",
                "resource": "STITCH_directional_v5", "mechanism": record["mechanism"],
                "evidence_tier": record["evidence_tier"], "original_edge_id": record["edge_id"],
                "source_accession": record["source"], "target_accession": record["target"],
                "source_pubchem": record["source_pubchem"], "raw_item_id_a": record["raw_item_id_a"],
                "raw_item_id_b": record["raw_item_id_b"], "raw_a_is_acting": record["raw_a_is_acting"],
                "raw_action_score": record["raw_action_score"],
                "raw_link_textmining": record["raw_link_textmining"],
                "raw_link_combined_score": record["raw_link_combined_score"],
            })
            handle.write(json.dumps(ledger, sort_keys=True) + "\n")

    # Preserve the base regulation ledger apart from unsafe cache STITCH rows.
    base_regulation = []
    with (BASE / "regulation_edges.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["resource"] != "STITCH":
                base_regulation.append(row)
    write_csv(output / "regulation_edges.csv", base_regulation + retained, REGULATION_FIELDS)
    write_csv(output / "stitch_directional_evidence.csv", retained, STITCH_EVIDENCE_FIELDS)
    write_csv(output / "stitch_directional_mapping_review.csv", mapping_review,
              ["entity_type", "input_identifier", "candidate_mappings", "mapping_count", "mapping_method", "admitted", "reason"])

    exclusions = [{"stage": "raw_actions", "reason": reason, "records": count}
                  for reason, count in sorted(action_counts.items())
                  if reason not in {"raw_action_records", "source_level_directional_signed_actions",
                                    "mapped_action_records_before_detailed_link_filter",
                                    "retained_directional_STITCH_evidence_rows"}]
    write_csv(output / "stitch_directional_exclusions.csv", exclusions, ["stage", "reason", "records"])

    with (output / "pkn_exploratory_stitch_directional.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source", "interaction", "target"])
        writer.writerows(graph)
    with (output / "conflicts_exploratory_stitch_directional.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source", "target"])
        writer.writerows(conflicts)

    # Preserve the explicit node contract and add only nodes introduced by the directional layer.
    base_nodes = list(jsonl(BASE / "nodes.jsonl"))
    node_by_id = {row["node_id"]: row for row in base_nodes}
    for record in retained:
        source_node, target_node = record["_graph_source"], record["_graph_target"]
        if source_node not in node_by_id:
            chemical = record["source"]
            node_by_id[source_node] = {
                "node_id": source_node, "node_type": "metabolite",
                "label": chebi_names.get(chemical, chemical), "chemical_id": chemical,
                "compartment": "unknown",
            }
        if target_node not in node_by_id:
            accession, gene = record["target"], record["target_gene"]
            node_by_id[target_node] = ({
                "node_id": target_node, "node_type": "gene_product", "label": gene,
                "gene_symbol": gene, "organism": 10090,
            } if gene else {
                "node_id": target_node, "node_type": "protein_unmapped", "label": accession,
                "uniprot": accession, "organism": 10090, "mouse_accession_verified": True,
            })
    nodes = sorted(node_by_id.values(), key=lambda row: row["node_id"])
    with (output / "nodes.jsonl").open("w") as handle:
        for row in nodes:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    node_fields = ["node_id", "label", "node_type"] + sorted({key for row in nodes for key in row} - {"node_id", "label", "node_type"})
    write_csv(output / "nodes.csv", nodes, node_fields)

    base_pkn_rows = sum(1 for _ in csv.DictReader((BASE / "pkn_exploratory.csv").open()))
    base_node_ids = {row["node_id"] for row in base_nodes}
    used_new_nodes = {node for source, _, target in graph for node in (source, target)
                      if node in node_by_id and node not in base_node_ids}
    summary = {
        "purpose": "Separate exploratory overlay replacing all orientation-unknown cached STITCH evidence with raw directionally auditable STITCH records.",
        "directionality_policy": {
            "raw_actor": "item_id_a must be CIDm/CIDs chemical",
            "raw_receiver": "item_id_b must be 10090.ENSMUSP protein",
            "acting_flag": "a_is_acting == t",
            "sign": "raw action exactly activation/inhibition; mode is provenance only",
            "minimum_action_score": minimum_action_score,
            "maximum_detailed_link_textmining_exclusive": maximum_textmining,
            "mapping": "unique current ChEBI via frozen UniChem and unique mouse UniProt via frozen STRING aliases",
            "conflicts": "opposing signs for an assembled pair are excluded; no majority vote",
        },
        "counts": dict(action_counts,
                       cached_STITCH_rows_removed=cache_stitch_removed,
                       base_non_STITCH_evidence_rows=base_evidence_rows,
                       baseline_exploratory_edges=base_pkn_rows,
                       directional_overlay_exploratory_edges=len(graph),
                       directional_overlay_conflicting_pairs=len(conflicts),
                       directional_unique_projected_signed_triples=len(fresh_projected_triples),
                       directional_projected_signed_triples_surviving_graph_conflict_filter=len(
                           fresh_projected_triples & final_graph_triples),
                       directional_projected_signed_triples_quarantined=len(
                           fresh_projected_triples - final_graph_triples),
                       directional_layer_nodes_added=len(used_new_nodes),
                       directional_raw_pairs_pending_after_mapping=len(pending)),
        "caveat": "STITCH v5.0 source files were last modified 2016-08-16. The retained edges have auditable directionality but remain exploratory; their direction does not establish currentness, direct mouse experimental support, hepatocyte context, or lipid measurement identity. A current STRING v12 alias snapshot supplies exact ENSMUSP-to-UniProt identity because the publisher no longer exposes a matched v5 alias archive; this affects identifier provenance, not raw action orientation.",
        "no_omics_scoring_performed": True,
    }
    (output / "stitch_directional_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    sources = [
        {"file": str(ACTIONS.relative_to(ROOT)), "url": ACTION_URL, "sha256_compressed": source_hashes["actions_compressed"], "sha256_uncompressed": source_hashes["actions_uncompressed"], "version": "STITCH 5.0; HTTP last-modified 2016-08-16"},
        {"file": str(LINKS.relative_to(ROOT)), "url": LINK_URL, "sha256_compressed": source_hashes["links_compressed"], "sha256_uncompressed": source_hashes["links_uncompressed"], "version": "STITCH 5.0; HTTP last-modified 2016-08-16"},
        {"file": str(STRING_ALIASES.relative_to(ROOT)), "url": STRING_ALIAS_URL, "sha256_compressed": source_hashes["string_aliases_compressed"], "version": "STRING v12.0 exact UniProt_AC aliases"},
        {"file": str(UNICHEM.relative_to(ROOT)), "url": UNICHEM_URL, "sha256_compressed": source_hashes["unichem_compressed"], "version": "UniChem whole-source mapping, listed 2026-09-12"},
        {"file": str(UNIPROT.relative_to(ROOT)), "sha256_compressed": source_hashes["uniprot_compressed"], "version": "frozen mouse UniProt export used by base PKN"},
        {"file": str(CHEBI.relative_to(ROOT)), "sha256_compressed": source_hashes["chebi_compressed"], "version": "frozen current ChEBI snapshot used by base PKN"},
        {"file": str(CHEBI_ALIASES.relative_to(ROOT)), "sha256": compressed_sha256(CHEBI_ALIASES), "version": "base PKN ChEBI alternate-ID table"},
    ]
    outputs = {}
    for path in sorted(output.iterdir()):
        if path.is_file():
            outputs[path.name] = {"sha256": compressed_sha256(path), "bytes": path.stat().st_size}
    base_graph_hashes = {
        name: compressed_sha256(BASE / name)
        for name in ["pkn_core.csv", "pkn_lipid_curated.csv", "pkn_lipid_rules.csv", "pkn_exploratory.csv"]
    }
    manifest = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "builder": str(Path(__file__).relative_to(ROOT)),
        "builder_sha256": compressed_sha256(Path(__file__)),
        "base_pkn": str(BASE.relative_to(ROOT)),
        "base_lipid_rules_sha256": compressed_sha256(BASE / "pkn_lipid_rules.csv"),
        "base_edge_evidence_sha256": compressed_sha256(BASE / "edge_evidence.jsonl.gz"),
        "base_graph_csv_sha256": base_graph_hashes,
        "sources": sources,
        "outputs": outputs,
        "policy": summary["directionality_policy"],
        "historical_source_caveat": summary["caveat"],
    }
    (output / "stitch_directional_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print("[STITCH] Retained", len(retained), "evidence rows; exported", len(graph),
          "signed exploratory edges; quarantined", len(conflicts), "conflicting pairs")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT,
                        help="New empty overlay directory (default: %(default)s)")
    parser.add_argument("--minimum-action-score", type=int, default=700)
    parser.add_argument("--maximum-textmining", type=int, default=700,
                        help="Exclude detailed links with textmining >= this value")
    args = parser.parse_args()
    if args.minimum_action_score < 0 or args.maximum_textmining < 0:
        raise ValueError("Score thresholds must be non-negative")
    build(args.output, args.minimum_action_score, args.maximum_textmining)


if __name__ == "__main__":
    main()
