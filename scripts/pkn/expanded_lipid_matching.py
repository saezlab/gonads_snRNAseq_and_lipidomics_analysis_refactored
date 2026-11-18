#!/usr/bin/env python3
"""Conservative structural matching for exploratory measured-lipid rules.

This module checks compatibility of *reported* molecular features. It does not
invent a reaction, direction, catalyst, substrate specificity, reaction sign,
or omitted co-substrate. Callers must separately validate a complete, directional
Rhea template and preserve its participants and stoichiometry.

API (all objects are ordinary dictionaries; no third-party dependencies):

``parse_feature(mapping_features_row)`` returns ``None`` for unsupported or
unparsed chemistry, or a feature with these fields:
  feature_id, canonical_name, class, family, chains, chain_multiset,
  chain_objects [{notation,carbon,double_bonds,oxygen}],
  sn_order_reported, sn_positions, retained_sn_position,
  sphingobase, n_acyl_chain, oxygenated_acyl_requires_template_review,
  headgroup_isomer_unresolved, glycerol_occupancy_unresolved,
  in_factor1, in_lipid_de.
``chains`` omits an explicitly empty ``0:0`` lyso position; ``sn_positions``
retains its nonempty position. For sphingolipids, the slash distinguishes the
sphingobase from N-acyl chain and NEVER denotes glycerol sn order. Oxygen counts
are retained, including reported N-acyl oxygenation; the caller must reject it
when unsupported by its Rhea template. No ether or oxidized glycerol features
are accepted. Sum-composition diacyl and sphingolipid features are excluded.

``match_headgroup(a,b)`` checks two diacyl glycerol features, preserving both
chains and any reported sn order. ``match_sphingolipid(a,b)`` preserves the
role-specific backbone and oxygen counts; neither matcher proves that the
headgroup conversion exists. Results contain ``rule_family``, feature IDs,
``chain_transform``, positional/sugar assumptions, and retained identities.

``match_lyso_diacyl(lyso,diacyl)`` checks the corresponding class pair (e.g.
LPC/PC), returning retained_chain, added_chain, removed_chain and a list of
``positional_options``. Each option exposes retained_sn, added_or_removed_sn,
and whether its position is conditional. The same structural pair can be used
in either direction ONLY when that direction is in the source Rhea template.
Known lyso sn positions inconsistent with known diacyl positions are rejected.
Unknown lyso sn positions remain conditional; they are never promoted to an
observed sn1 or sn2 assignment. A caller can filter these options against the
template's retained sn position, which distinguishes 1-acyl and 2-acyl lyso.

Run ``python3 scripts/pkn/expanded_lipid_matching.py --self-test`` to exercise
the compatibility checks and audit eligibility against the frozen mapping.
The CLI is read-only and emits JSON; it does not alter a PKN or measurements.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MAPPING = ROOT / "results/causal_modelling/cosmos_mouse_pkn/mapping_features.csv"
DIACYL_CLASSES = frozenset({"PC", "PE", "PS", "PA", "PI", "PG", "DG"})
LYSO_TO_DIACYL = {"LPC": "PC", "LPE": "PE", "LPS": "PS", "LPA": "PA", "LPI": "PI", "LPG": "PG"}
SPHINGOLIPID_CLASSES = frozenset({"Cer", "SM", "HexCer", "Hex2Cer"})
CHAIN = re.compile(r"(\d+):(\d+)(?:;O(\d*))?\Z")


def _truth(value):
    return str(value).lower() == "true"


def _chain(value, allow_oxygen=False, allow_empty=False):
    """Parse explicit chain counts without removing modification information."""
    match = CHAIN.fullmatch(value)
    if not match:
        return None
    carbon, double_bonds = int(match[1]), int(match[2])
    oxygen = int(match[3] or 1) if ";O" in value else 0
    if carbon == 0:
        if not allow_empty or double_bonds or oxygen:
            return None
    elif carbon < 2 or double_bonds >= carbon:
        return None
    if oxygen and not allow_oxygen:
        return None
    # Normalize only the equivalent explicitly counted ;O and ;O1 spelling.
    notation = f"{carbon}:{double_bonds}" + (f";O{oxygen}" if oxygen > 1 else ";O" if oxygen else "")
    return {"notation": notation, "carbon": carbon, "double_bonds": double_bonds, "oxygen": oxygen}


def parse_feature(row):
    """Return a role-aware explicit-chain feature, or None when unsupported.

    ``canonical_name`` is authoritative: do not retry a failed frozen-parser
    result using the raw label, and do not expand a sum into possible chains.
    """
    if row.get("parse_error") or not row.get("canonical_name"):
        return None
    value = row["canonical_name"].strip()
    match = re.fullmatch(r"([^ ()]+) (.+)", value)
    if not match:
        return None
    cls, body = match.groups()
    if cls != row.get("class", cls):
        return None
    result = {
        "feature_id": row.get("feature_id", value), "canonical_name": value,
        "class": cls, "in_factor1": _truth(row.get("in_factor1", False)),
        "in_lipid_de": _truth(row.get("in_lipid_de", False)),
        "sn_order_reported": False, "sn_positions": [], "retained_sn_position": None,
        "sphingobase": None, "n_acyl_chain": None,
        "oxygenated_acyl_requires_template_review": False,
        "headgroup_isomer_unresolved": False, "glycerol_occupancy_unresolved": False,
    }
    if cls in DIACYL_CLASSES:
        parts = re.split(r"([/_])", body)
        if len(parts) != 3:
            return None
        parsed = [_chain(parts[0]), _chain(parts[2])]
        if any(chain is None for chain in parsed):
            return None
        result.update(family="glycerol_diacyl", sn_order_reported=parts[1] == "/",
                      sn_positions=[1, 2] if parts[1] == "/" else [None, None],
                      glycerol_occupancy_unresolved=cls == "DG")
    elif cls in LYSO_TO_DIACYL:
        if "/" in body:
            parts = body.split("/")
            if len(parts) != 2:
                return None
            parsed_slots = [_chain(part, allow_empty=True) for part in parts]
            if any(chain is None for chain in parsed_slots):
                return None
            nonempty = [(i + 1, chain) for i, chain in enumerate(parsed_slots) if chain["carbon"]]
            if len(nonempty) != 1:
                return None
            position, chain = nonempty[0]
            parsed = [chain]
            result.update(sn_order_reported=True, sn_positions=[position], retained_sn_position=position)
        else:
            chain = _chain(body)
            if chain is None:
                return None
            parsed = [chain]
            result["sn_positions"] = [None]
        result["family"] = "glycerol_monoacyl"
    elif cls in SPHINGOLIPID_CLASSES:
        # Explicit sphingobase/N-acyl roles require a slash; unordered '_' or
        # summed oxygen does not establish which chain is the sphingobase.
        parts = body.split("/")
        if len(parts) != 2:
            return None
        base, acyl = (_chain(part, allow_oxygen=True) for part in parts)
        if base is None or acyl is None or base["oxygen"] not in (2, 3):
            return None
        # One additional oxygen in N-acyl is explicitly reportable but must be
        # supported by the chosen template. Unspecified/excess modification is
        # outside this initial exploratory matcher.
        if acyl["oxygen"] not in (0, 1):
            return None
        parsed = [base, acyl]
        result.update(family="sphingolipid", sphingobase=base["notation"],
                      n_acyl_chain=acyl["notation"],
                      oxygenated_acyl_requires_template_review=bool(acyl["oxygen"]),
                      headgroup_isomer_unresolved=cls in {"HexCer", "Hex2Cer"})
    else:
        return None
    result["chain_objects"] = parsed
    result["chains"] = [chain["notation"] for chain in parsed]
    result["chain_multiset"] = sorted(result["chains"])
    return result


def _glycerol_assumptions(a, b):
    assumptions = ["Double-bond positions, stereochemistry and enzyme species specificity remain unresolved."]
    if not a["sn_order_reported"] or not b["sn_order_reported"]:
        assumptions.append("Unreported glycerol sn positions remain conditional and are not assigned as measured identities.")
    if a["glycerol_occupancy_unresolved"] or b["glycerol_occupancy_unresolved"]:
        assumptions.append("A generic 1,2-diacyl Rhea template assumes that glycerol occupancy; chain-only DG measurements do not establish it.")
    return assumptions


def match_headgroup(a, b):
    """Check chain-compatible diacyl headgroup change; no chemistry is implied."""
    if not a or not b or a["family"] != "glycerol_diacyl" or b["family"] != "glycerol_diacyl":
        return None
    if a["class"] == b["class"] or a["chain_multiset"] != b["chain_multiset"]:
        return None
    if a["sn_order_reported"] and b["sn_order_reported"] and a["chains"] != b["chains"]:
        return None
    return {
        "rule_family": "diacyl_headgroup", "source_feature_id": a["feature_id"],
        "target_feature_id": b["feature_id"], "source_class": a["class"], "target_class": b["class"],
        "chain_transform": "preserve both reported ester acyl chains",
        "chain_multiset": a["chain_multiset"],
        "sn_assignment_conditional": not (a["sn_order_reported"] and b["sn_order_reported"]),
        "assumptions": _glycerol_assumptions(a, b),
    }


def match_lyso_diacyl(lyso, diacyl):
    """Match a retained lyso acyl chain to its corresponding diacyl class.

    The match itself has no reaction direction. In acylation ``added_chain``
    must match the template acyl donor; in hydrolysis ``removed_chain`` must
    match the released fatty acid. Both aliases contain the same chain.
    """
    if not lyso or not diacyl or lyso["family"] != "glycerol_monoacyl" or diacyl["family"] != "glycerol_diacyl":
        return None
    if LYSO_TO_DIACYL.get(lyso["class"]) != diacyl["class"]:
        return None
    retained = lyso["chains"][0]
    if retained not in diacyl["chains"]:
        return None
    alternatives = []
    for index, chain in enumerate(diacyl["chains"]):
        if chain != retained:
            continue
        possible_positions = [index + 1] if diacyl["sn_order_reported"] else [1, 2]
        for position in possible_positions:
            if lyso["retained_sn_position"] is not None and position != lyso["retained_sn_position"]:
                continue
            added = diacyl["chains"][1 - index]
            option = {
                "retained_sn": position, "added_or_removed_sn": 3 - position,
                "retained_chain": retained, "added_chain": added, "removed_chain": added,
                "retained_chain_diacyl_index": index,
                "lyso_sn_conditional": lyso["retained_sn_position"] is None,
                "diacyl_sn_conditional": not diacyl["sn_order_reported"],
                "position_conditional": lyso["retained_sn_position"] is None or not diacyl["sn_order_reported"],
            }
            # A repeated chain may yield an equivalent assignment; preserve
            # distinct sn options but do not duplicate chemically identical ones.
            if not any(all(other[key] == option[key] for key in option if key != "retained_chain_diacyl_index")
                       for other in alternatives):
                alternatives.append(option)
    if not alternatives:
        return None
    added_chains = {option["added_chain"] for option in alternatives}
    assert len(added_chains) == 1
    added = next(iter(added_chains))
    assumptions = _glycerol_assumptions(lyso, diacyl)
    assumptions.append("Rhea template direction, retained sn specificity, acyl donor/product chemistry and mouse catalyst evidence must be checked separately.")
    return {
        "rule_family": "lyso_diacyl_acyl_exchange", "lyso_feature_id": lyso["feature_id"],
        "diacyl_feature_id": diacyl["feature_id"], "lyso_class": lyso["class"], "diacyl_class": diacyl["class"],
        "chain_transform": "add/remove exactly one reported ester acyl chain; preserve the other",
        "retained_chain": retained, "added_chain": added, "removed_chain": added,
        "lyso_retained_sn_position": lyso["retained_sn_position"],
        "positional_options": alternatives, "assumptions": assumptions,
    }


def match_sphingolipid(a, b):
    """Check explicit backbone conservation, keeping chain roles and oxygen."""
    if not a or not b or a["family"] != "sphingolipid" or b["family"] != "sphingolipid":
        return None
    if a["class"] == b["class"] or a["sphingobase"] != b["sphingobase"] or a["n_acyl_chain"] != b["n_acyl_chain"]:
        return None
    assumptions = [
        "Reported sphingobase and N-acyl identities, including oxygen counts, are conserved separately.",
        "Double-bond/hydroxyl positions, stereochemistry and enzyme species specificity remain unresolved.",
        "The Rhea template's fixed sphingobase, oxygenation and complete paired-lipid/co-substrate chemistry must be checked separately.",
    ]
    if a["headgroup_isomer_unresolved"] or b["headgroup_isomer_unresolved"]:
        assumptions.append("HexCer/Hex2Cer do not resolve sugar identity, stereochemistry or linkage; a named glycan template is an explicitly conditional hypothesis.")
    if a["oxygenated_acyl_requires_template_review"] or b["oxygenated_acyl_requires_template_review"]:
        assumptions.append("Reported N-acyl oxygenation is retained; an unmodified-acyl template is not sufficient evidence for this substitution.")
    return {
        "rule_family": "sphingolipid_headgroup", "source_feature_id": a["feature_id"],
        "target_feature_id": b["feature_id"], "source_class": a["class"], "target_class": b["class"],
        "chain_transform": "preserve role-specific sphingobase and N-acyl chain including oxygen",
        "sphingobase": a["sphingobase"], "n_acyl_chain": a["n_acyl_chain"],
        "sphingobase_oxygen": a["chain_objects"][0]["oxygen"],
        "n_acyl_oxygen": a["chain_objects"][1]["oxygen"],
        "oxygenated_acyl_requires_template_review": a["oxygenated_acyl_requires_template_review"],
        "headgroup_isomer_conditional": a["headgroup_isomer_unresolved"] or b["headgroup_isomer_unresolved"],
        "assumptions": assumptions,
    }


def self_test():
    def feature(value):
        cls = value.split(" ")[0]
        return parse_feature({"feature_id": value, "canonical_name": value, "class": cls})
    # Reject sums, ethers, oxidized glycerol, malformed/unknown chains, parsing
    # failures, and role-unresolved sphingolipid chain annotations.
    for value in ("PC 36:1", "SM 36:1;O2", "Cer 36:1;O2", "PC O-16:0/18:1",
                  "PE P-16:0/18:1", "PC 16:0_18:2;O", "LPC O-16:0", "PC 0:0/18:1",
                  "LPC 16:0/18:1", "LPC 0:0/0:0", "Cer 18:1_18:0", "Cer 18:1/18:0",
                  "Cer 18:1;O4/18:0", "Cer 18:1;O2/18:0;O2", "PC unknown", "TG 16:0_18:1_18:1"):
        assert feature(value) is None, value
    assert parse_feature({"canonical_name": "PC 16:0_18:1", "parse_error": "failed"}) is None
    assert match_headgroup(feature("PC 16:0/18:1"), feature("PE 18:1/16:0")) is None
    assert match_headgroup(feature("PC 16:0/18:1"), feature("PE 16:0/18:1"))
    assert match_headgroup(feature("PC 16:0_18:1"), feature("PE 18:1_16:0"))
    assert match_headgroup(feature("PC 16:0_18:1"), feature("PE 16:0_18:2")) is None
    assert match_headgroup(feature("PC 16:0_18:1"), feature("PC 16:0_18:1")) is None
    # Both reported positions must agree, while unknown position remains an
    # explicit conditional alternative suitable for a template-specific gate.
    lyso1, lyso2 = feature("LPC 16:0/0:0"), feature("LPC 0:0/16:0")
    ordered = feature("PC 16:0/18:1")
    assert lyso1["retained_sn_position"] == 1 and lyso2["retained_sn_position"] == 2
    assert match_lyso_diacyl(lyso2, ordered) is None
    hit = match_lyso_diacyl(lyso1, ordered)
    assert hit["added_chain"] == "18:1" and len(hit["positional_options"]) == 1
    assert not hit["positional_options"][0]["position_conditional"]
    hit = match_lyso_diacyl(feature("LPC 16:0"), ordered)
    assert hit["positional_options"][0]["lyso_sn_conditional"]
    assert hit["lyso_retained_sn_position"] is None
    hit = match_lyso_diacyl(feature("LPC 16:0"), feature("PC 18:1_16:0"))
    assert {option["retained_sn"] for option in hit["positional_options"]} == {1, 2}
    assert all(option["position_conditional"] for option in hit["positional_options"])
    assert match_lyso_diacyl(feature("LPC 16:0"), feature("PC 18:0_18:1")) is None
    assert match_lyso_diacyl(feature("LPC 16:0"), feature("PE 16:0_18:1")) is None
    repeat = match_lyso_diacyl(feature("LPC 16:0"), feature("PC 16:0_16:0"))
    assert len(repeat["positional_options"]) == 2 and repeat["added_chain"] == "16:0"
    # A slash in a sphingolipid encodes chain roles, not glycerol sn positions.
    cer = feature("Cer 18:1;O2/18:0")
    assert not cer["sn_order_reported"] and cer["sn_positions"] == []
    assert match_sphingolipid(cer, feature("SM 18:1;O2/18:0"))
    assert match_sphingolipid(cer, feature("SM 18:0;O2/18:0")) is None
    assert match_sphingolipid(cer, feature("SM 18:1;O3/18:0")) is None
    assert match_sphingolipid(cer, feature("SM 18:1;O2/18:0;O")) is None
    hit = match_sphingolipid(feature("Cer 18:1;O2/18:0;O1"), feature("HexCer 18:1;O2/18:0;O"))
    assert hit["n_acyl_chain"] == "18:0;O" and hit["n_acyl_oxygen"] == 1
    assert hit["oxygenated_acyl_requires_template_review"] and hit["headgroup_isomer_conditional"]
    assert match_sphingolipid(cer, feature("Cer 18:1;O2/18:0")) is None


def audit_mapping(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    parsed = [value for row in rows if (value := parse_feature(row))]
    families = Counter(value["family"] for value in parsed)
    return {
        "mapping_file": str(path.resolve()),
        "mapping_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "total_mapping_features": len(rows), "eligible_explicit_chain_features": len(parsed),
        "eligible_by_family": dict(sorted(families.items())),
        "eligible_by_class": dict(sorted(Counter(value["class"] for value in parsed).items())),
        "eligible_historical_factor1_64_features": sum(value["in_factor1"] for value in parsed),
        "glycerol_diacyl_reported_sn_order": sum(value["family"] == "glycerol_diacyl" and value["sn_order_reported"] for value in parsed),
        "monoacyl_reported_sn_position": sum(value["family"] == "glycerol_monoacyl" and value["retained_sn_position"] is not None for value in parsed),
        "monoacyl_unreported_sn_position": sum(value["family"] == "glycerol_monoacyl" and value["retained_sn_position"] is None for value in parsed),
        "sphingolipids_requiring_oxygenated_acyl_template_review": sum(value["oxygenated_acyl_requires_template_review"] for value in parsed),
        "interpretation": "Structural eligibility only: template chemistry, direction, catalyst specificity and full reaction participants still require validation. in_factor1 is the historical 64-feature subset, not availability of Factor 1 weights.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    result = audit_mapping(args.mapping)
    result["self_test"] = "passed" if args.self_test else "not_requested"
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
