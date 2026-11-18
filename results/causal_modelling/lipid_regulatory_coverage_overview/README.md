# One row per measured lipid: regulatory-coverage overview

`all_measured_lipid_regulatory_coverage.csv` is the consolidated lookup for all1,254 measurements. It preserves the accepted155 and historical64 panel flags, original measurement identity and weight, frozen exact matches, conditional compatible candidates, broad-class annotations and newly detected specificity candidates. Formula-incompatible supplementary candidates are excluded. This is annotation only: no measurement, PKN or fit is changed.

Detailed provenance remains in the [frozen-crosswalk audit](../lipid_regulatory_mapping_audit/README.md) and [specificity audit](../lipid_regulatory_resolution_audit/README.md).

Evidence-group names keep primary-curated GtoPdb, secondary chemical mapping, cached BRENDA/STITCH and unsigned MRCLinksDB associations separate. CellPhoneDB coverage was not assessable from the frozen files. `per_resource_coverage_counts.csv` reports each group; `unique_feature_unions.csv` deduplicates overlapping features across categories and groups. Do not sum exact and compatible counts: FA16, for example, occurs in both while remaining excluded from the accepted155 inputs.

A compatible candidate is not a confirmed measured molecular identity. The current accepted155 inputs have no exact direct regulatory matches in the current curated PKN. Whole-dataset candidate opportunities are a different quantity. Both detailed audits have limited search scope, so unmatched features do not prove absent biological knowledge.

Reproduce with `python3 scripts/moon/ancillary/summarize_lipid_regulatory_coverage.py --output NEW_EMPTY_DIRECTORY`. Occupied output directories are refused, and `dependency_hashes.csv` records both input audits.
