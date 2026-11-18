# Direct lipid-to-protein matches in the frozen mapping

Run `python3 scripts/moon/ancillary/audit_direct_lipid_regulatory_mapping.py --output NEW_EMPTY_DIRECTORY`. This script compares the existing 1,254-feature crosswalk to raw cached BRENDA/STITCH/MRCLinksDB records and the retained regulation export, keeping those collections separate. No mapping, PKN or fit is changed.

**Scope: frozen-crosswalk coverage, not an exhaustive chemical search.** The original abbreviation-based mapper misses some structured SwissLipids entries with empty abbreviations. The companion `lipid_regulatory_resolution_audit` examines specificity/ontology alternatives beyond these frozen candidates. These tables must not be used alone to claim that no prior chemical knowledge exists.

`coverage_counts.csv` reports unique features for all measurements, the 155 accepted inputs, the historical 64 Factor1 features and the two LPC candidates. Exact same-ID, explicit alias, explicit acid/base, compatible resolution and class/ancestor matches remain distinct. Categories can overlap; do not add counts to obtain a union. `feature_edge_matches.csv` preserves each matched source record with actions, signs, evidence gates and source IDs. `feature_reason_funnel.csv` includes unmatched features; its grouped version is `reason_funnel_counts.csv`.

`resource_source_coverage.csv` separates raw and retained records. MRCLinksDB's positive numerical values are assumed by its adapter and were rejected as signed regulatory evidence. CellPhoneDB was not present in these frozen collections; its zero counts mean not assessable here. Raw cache provenance only establishes a mouse service export, not native-mouse evidence for every original source interaction.

`palmitate_ontology.csv` and `palmitate_swisslipids_candidates.csv` show why the frozen FA16 ambiguity is a specificity issue: explicit straight-chain hexadecanoate is a subclass of the generic saturated C16 anion, rather than an unrelated alternative identity. Its acid/base overlap with palmitic acid is explicit; the existing input exclusion is preserved for this audit.

All numerical counts are frozen-data audit results, not causal confidence or experimental evidence counts. `dependency_hashes.csv` identifies the consumed resources and implementation.
