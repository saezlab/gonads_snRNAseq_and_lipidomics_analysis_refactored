#!/usr/bin/env python3
"""Build traceable mouse chemical -> protein activity regulation for COSMOS.

The current, mouse-specific Guide to Pharmacology (GtoPdb) layer is a fresh
addition to the OmniPath-metabo resource family. Cached OmniPath-metabo BRENDA
and STITCH edges remain separate fallback evidence. These are activity edges,
not direct transcriptional regulation. MRCLinksDB's unconditional +1 assignment
is NOT retained as evidence of activation.

Default execution is offline and deterministic. --refresh downloads the four
current public GtoPdb CSVs and optional PubChem synonym mapping. --refresh-chebi
downloads the ChEBI ontology used for identity crosswalks. No assay data,
correlations, or measured lipid abundance is used to invent regulatory signs.
"""

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import html
import io
import json
import re
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "support/cosmos_mouse/raw"
OUT = ROOT / "results/cosmos_mouse_pkn"
CACHED = ROOT / "results/cosmos_pkn_audit/raw/om-mouse-pkn.json.gz"
BASE = "https://www.guidetopharmacology.org/DATA/"
GTOP_FILES = ("interactions.csv", "ligand_id_mapping.csv", "ligands.csv",
              "targets_and_families.csv")
PUBCHEM_NAME = "regulation_pubchem_missing_metabolite_synonyms.json.gz"
UNIPROT = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$")

FIELDS = ["edge_id", "source", "target", "interaction", "resource",
          "evidence_tier", "mechanism", "organism", "source_organism",
          "references", "compartment", "source_version", "source_type",
          "target_type", "source_name", "target_gene", "target_class",
          "ligand_type", "endogenous", "source_record_id", "source_file",
          "source_sha256_uncompressed", "mapping_method", "mapping_count",
          "mapping_group", "source_pubchem", "original_action",
          "original_type", "sign_basis", "affinity_units", "affinity_value",
          "assay_description", "sign_conflict", "quality_notes"]


def get_bytes(path):
    with gzip.open(path, "rb") as handle:
        return handle.read()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_gtop(name):
    data = get_bytes(RAW / ("regulation_gtopdb_" + name + ".gz"))
    lines = data.decode("utf-8-sig").splitlines()
    version = lines[0].strip('"')
    assert version.startswith("# GtoPdb Version:"), name + " is not a GtoPdb CSV"
    return list(csv.DictReader(io.StringIO("\n".join(lines[1:])))), version, digest(data)


def fetch(url, destination):
    request = urllib.request.Request(url, headers={"User-Agent": "FCGL-COSMOS-PKN-build/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = response.read()
    if destination.name.startswith("regulation_gtopdb_") and b"GtoPdb Version:" not in data[:200]:
        raise ValueError("Expected public GtoPdb CSV, got another document: " + url)
    with gzip.GzipFile(filename=str(destination), mode="wb", mtime=0) as handle:
        handle.write(data)
    return data


def refresh():
    RAW.mkdir(parents=True, exist_ok=True)
    for name in GTOP_FILES:
        fetch(BASE + name, RAW / ("regulation_gtopdb_" + name + ".gz"))
    interactions, _, _ = read_gtop("interactions.csv")
    mappings, _, _ = read_gtop("ligand_id_mapping.csv")
    by_id = {row["Ligand id"]: row for row in mappings}
    cids = sorted({by_id[row["Ligand ID"]]["PubChem CID"] for row in interactions
                   if row["Target Species"] == "Mouse" and row["Ligand Type"] == "Metabolite"
                   and not by_id[row["Ligand ID"]]["Chebi ID"]
                   and by_id[row["Ligand ID"]]["PubChem CID"].isdigit()})
    if cids:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + ",".join(cids) + "/synonyms/JSON"
        fetch(url, RAW / PUBCHEM_NAME)


def refresh_chebi():
    RAW.mkdir(parents=True, exist_ok=True)
    url = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz"
    path = RAW / "regulation_chebi.obo.gz"
    with urllib.request.urlopen(url, timeout=120) as response:
        headers = dict(response.headers)
        with path.open("wb") as handle:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                handle.write(block)
    Path(str(path) + ".metadata.json").write_text(json.dumps(dict(url=url, headers=headers,
        retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()), indent=2) + "\n")


def chebi_identity_tables(regulatory_sources):
    """Export curated protonation relations and aliases; do not change graph edges.

    Current ChEBI uses RO:0018033/4 for conjugate base/acid. Enantiomer,
    tautomer, functional-parent, and is_a relations are not identity mappings.
    is_a is used only to annotate membership in the lipid hierarchy.
    """
    path = RAW / "regulation_chebi.obo.gz"
    if not path.exists():
        return None
    terms, aliases, raw_relations, header = {}, {}, [], []
    current = None
    in_term = False
    def finish(term):
        if term and term.get("id"):
            terms[term["id"]] = term
            for alias in term["aliases"]:
                assert alias not in aliases or aliases[alias] == term["id"]
                aliases[alias] = term["id"]
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("["):
                finish(current)
                current = None
                in_term = line == "[Term]"
                if in_term:
                    current = {"aliases": [], "parents": [], "charge": "", "formula": "", "obsolete": False}
                continue
            if not in_term:
                if not terms:
                    header.append(line)
                continue
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("alt_id: "):
                current["aliases"].append(line[8:])
            elif line.startswith("is_a: "):
                current["parents"].append(line.split()[1])
            elif line == "is_obsolete: true":
                current["obsolete"] = True
            elif line.startswith("relationship: "):
                _, relation, target = line.split()[:3]
                if relation in {"RO:0018033", "RO:0018034", "is_conjugate_acid_of", "is_conjugate_base_of"}:
                    raw_relations.append((current["id"], relation, target))
            elif line.startswith("property_value: chemrof:charge "):
                current["charge"] = line.split('"')[1]
            elif line.startswith("property_value: chemrof:generalized_empirical_formula "):
                current["formula"] = line.split('"')[1]
    finish(current)
    version = next((s.split(": ", 1)[1] for s in header if s.startswith("data-version: ")), "unknown")
    # Classify all terms by is_a only; no inheritance of regulatory edges follows.
    children = collections.defaultdict(set)
    for ident, term in terms.items():
        for parent in term["parents"]:
            children[parent].add(ident)
    lipids, stack = set(), ["CHEBI:18059"]
    while stack:
        ident = stack.pop()
        if ident not in lipids:
            lipids.add(ident)
            stack.extend(children.get(ident, ()))
    direct = collections.defaultdict(set)
    for source, relation, target in raw_relations:
        source, target = aliases.get(source, source), aliases.get(target, target)
        acid, base = (source, target) if relation in {"RO:0018034", "is_conjugate_acid_of"} else (target, source)
        direct[(acid, base)].add(relation)
    adjacency = collections.defaultdict(set)
    relation_rows = []
    for (acid, base), relations in sorted(direct.items()):
        ta, tb = terms.get(acid, {}), terms.get(base, {})
        # Explicit source relation is retained without substituting structural similarity.
        relation_rows.append(dict(acid=acid, base=base, acid_name=ta.get("name", ""), base_name=tb.get("name", ""),
                                  acid_charge=ta.get("charge", ""), base_charge=tb.get("charge", ""),
                                  acid_formula=ta.get("formula", ""), base_formula=tb.get("formula", ""),
                                  ontology_relations=";".join(sorted(relations)), relation="is_conjugate_acid_of",
                                  resource="ChEBI", source_version=version,
                                  changes_to_model="none; identity crosswalk only, not a causal interaction"))
        adjacency[acid].add(base)
        adjacency[base].add(acid)
    ph73_path = OUT / "rhea_chebi_ph73.tsv"
    ph73 = {}
    if ph73_path.exists():
        with ph73_path.open() as handle:
            ph73 = {r["source"]: r for r in csv.DictReader(handle, delimiter="\t")}
    crosswalk = []
    for original in sorted(regulatory_sources):
        primary = aliases.get(original, original)
        t = terms.get(primary, {})
        component, stack = set(), [primary]
        while stack:
            ident = stack.pop()
            if ident not in component:
                component.add(ident)
                stack.extend(adjacency.get(ident, ()))
        ph = ph73.get(original, ph73.get(primary, {}))
        crosswalk.append(dict(source=original, primary_chebi=primary, name=t.get("name", ""),
                              present_in_current_chebi=bool(t), obsolete=t.get("obsolete", "unknown"),
                              alias_mapping="explicit_ChEBI_alt_id" if original != primary else "identity",
                              charge=t.get("charge", ""), formula=t.get("formula", ""),
                              lipid_is_a=primary in lipids, acid_base_component=";".join(sorted(component)),
                              component_size=len(component), source_version=version,
                              rhea_ph73_target=ph.get("target", ""), rhea_ph73_origin=ph.get("origin", ""),
                              rhea_ph73_target_in_explicit_acid_base_component=ph.get("target", "") in component if ph else "unknown",
                              interpretation="Optional protonation-state equivalence; no compartment, enantiomer, tautomer, or structural-class merge."))
    alias_rows = [dict(alias=k, primary_chebi=v, resource="ChEBI", source_version=version, relation="alt_id")
                  for k, v in sorted(aliases.items())]
    write_csv("regulation_chebi_acid_base_relations.csv", relation_rows, list(relation_rows[0]))
    write_csv("regulation_chebi_aliases.csv", alias_rows, ["alias", "primary_chebi", "resource", "source_version", "relation"])
    write_csv("regulation_chebi_source_crosswalk.csv", crosswalk, list(crosswalk[0]))
    metadata_path = Path(str(path) + ".metadata.json")
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    summary = dict(source_version=version, terms=len(terms), aliases=len(aliases), lipid_descendants=len(lipids),
                   direct_acid_base_pairs=len(direct), regulatory_sources=len(crosswalk),
                   regulatory_sources_mapped_by_alt_id=sum(r["source"] != r["primary_chebi"] for r in crosswalk),
                   regulatory_sources_with_conjugate_forms=sum(r["component_size"] > 1 for r in crosswalk),
                   regulatory_sources_lipid_is_a=sum(r["lipid_is_a"] for r in crosswalk),
                   source_file=str(path.relative_to(ROOT)), sha256_compressed=digest(path.read_bytes()),
                   retrieval=metadata,
                   optional_rhea_ph73_input=(dict(file=str(ph73_path.relative_to(ROOT)), sha256=digest(ph73_path.read_bytes()))
                                            if ph73_path.exists() else None),
                   limitations=["These are explicit ontology protonation relations, not experimentally measured reaction direction or activity.",
                                "A connected component can span multiple protonation states; optionally collapsing it assumes the measured chemical represents that ensemble.",
                                "Only RO:0018033/RO:0018034 and alt_id mappings are identity candidates. No enantiomer/tautomer/functional-parent merge is made.",
                                "No PKN edges are modified by this export."])
    (OUT / "regulation_chebi_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
    assert ("CHEBI:37550", "CHEBI:60119") in direct, "Expected curated S1P conjugate pair missing"
    assert ("CHEBI:15843", "CHEBI:32395") in direct, "Expected arachidonate conjugate pair missing"
    assert ("CHEBI:16393", "CHEBI:46967") not in direct, "Enantiomer relation must never be an acid/base mapping"
    return summary


def clean_name(value):
    return html.unescape(re.sub(r"<[^>]*>", "", value))


def chebi_structures(wanted):
    """Audit foreign cross-references against primary ChEBI structures.

    Only compare exact InChIKeys; a stereo or protonation discrepancy is marked
    for review rather than repaired automatically. Retain exact-key candidates
    for a transparent identifier correction in a future build.
    """
    path = RAW / "regulation_chebi.obo.gz"
    if not path.exists():
        raise FileNotFoundError("ChEBI ontology is required to validate GtoPdb chemical mappings; run with --refresh-chebi.")
    lookup, by_key = {}, collections.defaultdict(list)
    current = None
    def finish(term):
        if not term or not term.get("id") or term.get("obsolete"):
            return
        if term.get("key"):
            by_key[term["key"]].append(term["id"])
        if term["id"] in wanted or wanted.intersection(term["aliases"]):
            for ident in [term["id"]] + term["aliases"]:
                lookup[ident] = dict(term)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("["):
                finish(current)
                current = {"aliases": []} if line == "[Term]" else None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("alt_id: "):
                    current["aliases"].append(line[8:])
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
                elif line.startswith("property_value: chemrof:inchi_key_string "):
                    current["key"] = line.split('"')[1]
        finish(current)
    return lookup, by_key


def chebi_ids(value):
    return sorted(set(re.findall(r"CHEBI:[0-9]+", str(value))))


def sign_from_action(action, kind):
    """Map explicit pharmacological effects; never equate binding with activation."""
    positive = {"Full agonist", "Agonist", "Partial agonist", "Activation", "Positive",
                "Potentiation", "Biased agonist", "Irreversible agonist", "Slows inactivation"}
    negative = {"Inhibition", "Antagonist", "Inverse agonist", "Negative", "Pore blocker",
                "Irreversible inhibition", "Feedback inhibition", "Voltage-dependent inhibition"}
    if action in positive:
        return 1, "explicit_action"
    if action in negative:
        return -1, "explicit_action"
    if action in {"Binding", "Neutral", "Mixed", "Biphasic", "Unknown"}:
        return None, "unsigned_or_context_dependent_action"
    # Some records explicitly curate Activator/Inhibitor in Type but leave Action empty.
    if action in {"", "None"}:
        if kind in {"Activator", "Agonist"}:
            return 1, "explicit_type_action_unspecified"
        if kind in {"Inhibitor", "Antagonist", "Channel blocker", "Gating inhibitor"}:
            return -1, "explicit_type_action_unspecified"
    return None, "no_unambiguous_effect"


def write_csv(name, rows, fields=None):
    fields = fields or FIELDS
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def blank_edge(**kwargs):
    result = dict.fromkeys(FIELDS, "")
    result.update(source_type="small_molecule", target_type="protein", organism=10090,
                  source_organism=10090, compartment="unknown", mapping_count=1)
    result.update(kwargs)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--refresh-chebi", action="store_true")
    args = parser.parse_args()
    if args.refresh:
        refresh()
    if args.refresh_chebi:
        refresh_chebi()
    OUT.mkdir(parents=True, exist_ok=True)
    interactions, version, interaction_sha = read_gtop("interactions.csv")
    mappings, mapping_version, _ = read_gtop("ligand_id_mapping.csv")
    ligand_structures, ligand_version, _ = read_gtop("ligands.csv")
    targets, target_version, _ = read_gtop("targets_and_families.csv")
    assert version == mapping_version == target_version == ligand_version, "Mixed GtoPdb versions"
    ligand_map = {r["Ligand id"]: r for r in mappings}
    ligand_structures = {r["Ligand ID"]: r for r in ligand_structures}
    target_map = {r["Target id"]: r for r in targets if not r["Subunit id"]}
    target_types = {r["Target id"]: r["Type"] for r in targets}
    with gzip.open(RAW / "uniprot_mouse.tsv.gz", "rt") as handle:
        mouse_proteins = {r["Entry"]: r["Gene Names (primary)"] for r in csv.DictReader(handle, delimiter="\t")}
    pubchem = {}
    pubchem_path = RAW / PUBCHEM_NAME
    if pubchem_path.exists():
        for rec in json.loads(get_bytes(pubchem_path))["InformationList"]["Information"]:
            pubchem[str(rec["CID"])] = sorted({s for s in rec.get("Synonym", [])
                                             if re.fullmatch(r"CHEBI:[0-9]+", s)})
    wanted = {s for r in mappings for s in chebi_ids(r["Chebi ID"])}
    wanted.update(s for values in pubchem.values() for s in values)
    chemical_structures, by_inchikey = chebi_structures(wanted)
    evidence, excluded, mapping_review = [], [], []
    def reject(resource, row_id, reason, source, target, detail=""):
        excluded.append(dict(resource=resource, source_record_id=row_id, reason=reason,
                             original_source=source, original_target=target, detail=detail))
    fresh_pharmacological = []
    for number, row in enumerate(interactions, start=1):
        if row["Target Species"] != "Mouse":
            continue
        row_id = "GtoPdb:interaction_record:" + str(number)
        if row["Ligand Type"] in {"Peptide", "Antibody", "Endogenous peptide"}:
            reject("GuideToPharmacology", row_id, "not_small_molecule", row["Ligand"], row["Target"])
            continue
        sign, basis = sign_from_action(row["Action"], row["Type"])
        if sign is None:
            reject("GuideToPharmacology", row_id, basis, row["Ligand"], row["Target"], row["Action"])
            continue
        ligand = ligand_map[row["Ligand ID"]]
        sources = chebi_ids(ligand["Chebi ID"])
        mapping = "GtoPdb_curated_ChEBI_cross_reference"
        if not sources and row["Ligand Type"] == "Metabolite":
            sources = pubchem.get(ligand["PubChem CID"], [])
            mapping = "GtoPdb_PubChem_CID_to_PubChem_ChEBI_synonym"
        if not sources:
            reject("GuideToPharmacology", row_id, "no_ChEBI_mapping", row["Ligand"], row["Target"],
                   "GtoPdb:" + row["Ligand ID"] + ";PubChem:" + ligand["PubChem CID"])
            continue
        target = row["Target UniProt ID"].strip()
        target_method = "interaction_target_UniProt"
        if not target:
            target = target_map.get(row["Target ID"], {}).get("Mouse SwissProt", "")
            target_method = "GtoPdb_target_table_Mouse_SwissProt"
        if not UNIPROT.fullmatch(target) or target not in mouse_proteins:
            reject("GuideToPharmacology", row_id, "unresolved_or_multisubunit_mouse_target", row["Ligand"], target,
                   row["Target"])
            continue
        target_class = target_types.get(row["Target ID"], "unknown")
        mechanism = ("nuclear_receptor_activity_modulation" if target_class == "nhr" else
                     "enzyme_activity_modulation" if target_class == "enzyme" else
                     "transporter_activity_modulation" if target_class == "transporter" else
                     "ligand_receptor_or_channel_activity_modulation")
        notes = ["Activity sign does not imply direct gene-expression regulation."]
        if row["Action"] in {"Partial agonist", "Biased agonist", "Antagonist"}:
            notes.append("Effect depends on receptor state, competing ligands, and assay context.")
        if mapping.startswith("GtoPdb_PubChem"):
            notes.append("Secondary identifier mapping; PubChem synonyms can be ambiguous or obsolete.")
        if len(sources) > 1:
            notes.append("Multiple ChEBI alternatives retained; not independent chemical evidence.")
        for source in sources:
            chemical = chemical_structures.get(source, {})
            chebi_key = chemical.get("key", "")
            gtop_key = ligand_structures[row["Ligand ID"]]["InChIKey"]
            if chebi_key and gtop_key and chebi_key != gtop_key:
                reason = ("chemical_connectivity_mismatch" if chebi_key.split("-")[0] != gtop_key.split("-")[0]
                          else "chemical_stereochemistry_or_protonation_discrepancy")
                reject("GuideToPharmacology", row_id, reason, row["Ligand"], target,
                       source + ";ChEBI_name=" + chemical.get("name", ""))
                mapping_review.append(dict(ligand_id=row["Ligand ID"], ligand_name=clean_name(row["Ligand"]),
                    target=target, proposed_chebi=source, chebi_name=chemical.get("name", ""),
                    gtop_inchikey=gtop_key, chebi_inchikey=chebi_key, reason=reason,
                    exact_inchikey_chebi_candidates=";".join(sorted(by_inchikey.get(gtop_key, []))),
                    action="quarantined; no automatic identifier correction", source_record_id=row_id))
                continue
            chemical_note = ("Exact GtoPdb and ChEBI InChIKeys agree." if chebi_key and gtop_key else
                             "Chemical cross-reference lacks a structure pair for exact verification; may denote a generic lipid class.")
            edge = blank_edge(source=source, target=target, interaction=sign,
                resource="GuideToPharmacology", evidence_tier=("curated_mouse_current" if mapping.startswith("GtoPdb_curated")
                    else "curated_mouse_current_secondary_chemical_mapping"), mechanism=mechanism,
                references=";".join("PMID:" + x for x in row["PubMed ID"].split("|") if x),
                source_version=version, source_name=clean_name(row["Ligand"]),
                target_gene=mouse_proteins[target], target_class=target_class,
                ligand_type=row["Ligand Type"], endogenous=row["Endogenous"], source_record_id=row_id,
                source_file="support/cosmos_mouse/raw/regulation_gtopdb_interactions.csv.gz",
                source_sha256_uncompressed=interaction_sha, mapping_method=mapping + ";" + target_method,
                mapping_count=len(sources), mapping_group="GtoPdb:" + row["Ligand ID"],
                source_pubchem=ligand["PubChem CID"], original_action=row["Action"], original_type=row["Type"],
                sign_basis=basis, affinity_units=row["Affinity Units"],
                affinity_value=row["Affinity Median"] or "|".join([row["Affinity Low"],row["Affinity High"]]),
                assay_description="; ".join(x for x in [row["Assay Description"], row["Action comment"],
                                                       row["Ligand Context"], row["Webpage URLs"], row["Patent Numbers"]] if x),
                quality_notes=" ".join(notes + [chemical_note, "Current ChEBI label: " + chemical.get("name", "unknown")]))
            if row["Ligand Type"] == "Metabolite":
                evidence.append(edge)
            else:
                edge["evidence_tier"] = "pharmacological_nonmetabolite_supplement"
                fresh_pharmacological.append(edge)

    cache_bytes = get_bytes(CACHED)
    cache_sha = digest(cache_bytes)
    cached_records = json.loads(cache_bytes)
    assert cached_records["meta"]["organism"] == 10090
    for number, row in enumerate(cached_records["network"]):
        resource = row["resource"]
        if resource not in {"BRENDA", "STITCH", "MRCLinksDB"}:
            continue
        row_id = "OmniPath-metabo:mouse-cache:record:" + str(number)
        if resource == "MRCLinksDB":
            reject(resource, row_id, "upstream_assumed_sign_not_regulatory_evidence", row["source"], row["target"],
                   "Upstream adapter assigns +1 to all ligand and transporter associations.")
            continue
        if row["source_type"] != "small_molecule" or row["target_type"] != "protein" or row["mor"] not in {-1, 1}:
            reject(resource, row_id, "not_signed_chemical_to_protein", row["source"], row["target"])
            continue
        mode = row.get("attrs", {}).get("stitch_mode")
        if resource == "STITCH" and (mode not in {"activation", "inhibition"} or row["mor"] != (1 if mode == "activation" else -1)):
            reject(resource, row_id, "STITCH_action_not_explicit_or_sign_mismatch", row["source"], row["target"], str(mode))
            continue
        source = str(row["source"])
        if source.isdigit():
            source = "CHEBI:" + source
        if not re.fullmatch(r"CHEBI:[0-9]+", source):
            reject(resource, row_id, "invalid_CHEBI_identifier", row["source"], row["target"])
            continue
        target = row["target"]
        if target not in mouse_proteins:
            reject(resource, row_id, "not_in_current_mouse_UniProt_export", row["source"], target)
            continue
        evidence.append(blank_edge(source=source, target=target, interaction=row["mor"], resource=resource,
            evidence_tier="cached_" + resource.lower() + "_fallback",
            mechanism="enzyme_activity_modulation" if resource == "BRENDA" else "chemical_protein_activity_modulation",
            references="", compartment=";".join(sorted(row.get("locations") or [])) or "unknown",
            source_version="OmniPath-metabo service cache component modified 2026-04-21; downloaded 2026-09-08; upstream release unspecified",
            target_gene=mouse_proteins[target], target_class=row["interaction_type"], ligand_type="unclassified_chemical",
            source_record_id=row_id, source_file=str(CACHED.relative_to(ROOT)),
            source_sha256_uncompressed=cache_sha, mapping_method="upstream_cached_mapping_not_independently_reconstructed",
            mapping_count="unknown", mapping_group="unknown", original_action=mode or ("activator" if row["mor"] == 1 else "inhibitor"),
            original_type=row["interaction_type"], sign_basis="upstream_explicit_action",
            quality_notes=("Frozen fallback; per-edge references, original chemical IDs and mapping ambiguity absent from service export. "
                + ("STITCH protein-to-chemical action orientation cannot be independently audited from this cache; exploratory only. " if resource == "STITCH" else "")
                + "Resource membership does not establish that the chemical is a lipid. Protein activity, not direct transcription.")))

    # Preserve contradictory evidence; the parent integrator chooses graph policy.
    pair_signs = collections.defaultdict(set)
    for row in evidence:
        pair_signs[(row["source"], row["target"])].add(row["interaction"])
    conflicts = []
    for row in evidence:
        row["sign_conflict"] = len(pair_signs[(row["source"], row["target"])]) > 1
        key = "|".join(str(row[k]) for k in ["resource", "source_record_id", "source", "target", "interaction"])
        row["edge_id"] = "reg_" + digest(key.encode())[:20]
        if row["sign_conflict"]:
            conflicts.append(row)
    for row in fresh_pharmacological:
        key = "|".join(str(row[k]) for k in ["resource", "source_record_id", "source", "target", "interaction"])
        row["edge_id"] = "reg_" + digest(key.encode())[:20]
    evidence.sort(key=lambda r: (r["resource"], r["source"], r["target"], r["interaction"], r["edge_id"]))
    write_csv("regulation_edges.csv", evidence)
    write_csv("regulation_sign_conflicts.csv", conflicts)
    write_csv("regulation_pharmacological_supplement.csv", fresh_pharmacological)
    write_csv("regulation_exclusions.csv", excluded, ["resource", "source_record_id", "reason", "original_source", "original_target", "detail"])
    write_csv("regulation_mapping_review.csv", mapping_review,
              ["ligand_id", "ligand_name", "target", "proposed_chebi", "chebi_name", "gtop_inchikey", "chebi_inchikey", "reason",
               "exact_inchikey_chebi_candidates", "action", "source_record_id"])
    sources = []
    for name in GTOP_FILES:
        path = RAW / ("regulation_gtopdb_" + name + ".gz")
        sources.append(dict(file=str(path.relative_to(ROOT)), url=BASE+name, sha256_uncompressed=digest(get_bytes(path)),
                            source_version=version, local_snapshot_mtime_utc=datetime.datetime.fromtimestamp(path.stat().st_mtime, datetime.timezone.utc).isoformat(),
                            license="GtoPdb database ODbL; contents CC BY-SA 4.0"))
    if pubchem_path.exists():
        cids = sorted(pubchem)
        sources.append(dict(file=str(pubchem_path.relative_to(ROOT)),
                            url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + ",".join(cids) + "/synonyms/JSON",
                            sha256_uncompressed=digest(get_bytes(pubchem_path)), source_version="live PUG REST",
                            local_snapshot_mtime_utc=datetime.datetime.fromtimestamp(pubchem_path.stat().st_mtime, datetime.timezone.utc).isoformat()))
    sources += [dict(file=str(CACHED.relative_to(ROOT)), url="https://metabo.omnipathdb.org/cosmos/pkn?organism=10090&categories=all",
                     sha256_uncompressed=cache_sha, source_version="cached component build 2026-04-21, retrieved 2026-09-08"),
                dict(file="support/cosmos_mouse/raw/uniprot_mouse.tsv.gz", url="https://rest.uniprot.org/uniprotkb/search",
                     query="organism_id:10090", sha256_uncompressed=digest(get_bytes(RAW / "uniprot_mouse.tsv.gz")),
                     source_version="mouse accession validation; parent builder manifest records retrieval"),
                dict(file="support/cosmos_mouse/raw/regulation_chebi.obo.gz", url="https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz",
                     sha256_compressed=digest((RAW / "regulation_chebi.obo.gz").read_bytes()),
                     source_version="Frozen ChEBI ontology; regulation_chebi_manifest.json records release and retrieval",
                     purpose="Audit GtoPdb cross-reference structures; export explicit protonation and alt_id crosswalks")]
    summary = dict(evidence_rows=len(evidence), unique_signed_edges=len({(r["source"],r["target"],r["interaction"]) for r in evidence}),
                   source_chemicals=len({r["source"] for r in evidence}), target_proteins=len({r["target"] for r in evidence}),
                   evidence_by_tier=dict(collections.Counter(r["evidence_tier"] for r in evidence)),
                   exclusion_reasons=dict(collections.Counter(r["reason"] for r in excluded)),
                   conflicting_pairs=sum(len(s) > 1 for s in pair_signs.values()),
                   pharmacological_supplement_rows=len(fresh_pharmacological),
                   direct_gene_expression_edges=0)
    limitations = ["Regulatory signs describe target protein activity; downstream transcription requires a separate signed GRN.",
                   "Only GtoPdb interactions explicitly measured/curated for Mouse are in the current tier; human data are not projected.",
                   "GtoPdb Metabolite classification is used for the metabolic tier and does not prove that the compound is a measured lipid.",
                   "Partial/biased agonist and antagonist signs are assay-context dependent; affinity and references are preserved.",
                   "PubChem synonym-derived chemical mappings are an explicit secondary tier; multi-ID mappings are alternatives.",
                   "Fresh chemical cross-references with different GtoPdb/ChEBI InChIKeys are quarantined, including stereo/protonation discrepancies; generic entities may lack InChIKeys.",
                   "Cached BRENDA/STITCH release, original chemical mappings, and per-edge references are unavailable in the service export.",
                   "Cached STITCH causal orientation cannot be independently reconstructed; use only as exploratory fallback.",
                   "MRCLinksDB associations are excluded because the upstream +1 is unconditional, including transport.",
                   "CellPhoneDB nonpeptidic adapter is human-only and assigns +1 without checking effect; not silently projected to mouse.",
                   "Sign conflicts are retained as evidence and reported, not resolved by majority vote."]
    (OUT / "regulation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (OUT / "regulation_manifest.json").write_text(json.dumps(dict(builder="scripts/pkn/build_lipid_regulation.py", sources=sources,
                                                                 limitations=limitations), indent=2) + "\n")
    # Meaningful invariants: signs, taxonomy validation, unique evidence IDs, and source integrity.
    assert all(r["interaction"] in {-1, 1} and r["target"] in mouse_proteins for r in evidence)
    assert all(re.fullmatch(r"CHEBI:[0-9]+", r["source"]) for r in evidence)
    assert all(r["resource"] != "MRCLinksDB" for r in evidence)
    assert len({r["edge_id"] for r in evidence}) == len(evidence)
    assert any(r["resource"] == "GuideToPharmacology" and r["source_name"] == "sphingosine 1-phosphate" for r in evidence)
    chebi_summary = chebi_identity_tables({r["source"] for r in evidence})
    if chebi_summary:
        summary["chebi_identity_mapping"] = {k: v for k, v in chebi_summary.items()
                                             if isinstance(v, (str, int)) and k != "sha256_compressed"}
        (OUT / "regulation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
