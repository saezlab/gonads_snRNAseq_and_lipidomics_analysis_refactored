# Frozen sources for the experimental mouse COSMOS prior

Raw downloads are under `raw/`. Preserve them with the output manifests for
offline reproduction. The complete build is documented in
`docs/cosmos_mouse_pkn_build.md`; source checksums are recorded in
`results/cosmos_mouse_pkn/build_manifest.json` and resource-specific manifests.

Most source adapters have explicit download/refresh options in their `--help`.
The parent build also fetched these sources directly:

- Mouse-GEM SBML and annotation tables from
  `https://raw.githubusercontent.com/SysBioChalmers/Mouse-GEM/fc769ac1da421c4e4e9e1cf2a9aec3fe4f02ef6a/model/`
  (`Mouse-GEM.xml`, `metabolites.tsv`, `reactions.tsv`). Locally the TSV files
  are `mouse_gem_metabolites.tsv` and `mouse_gem_reactions.tsv`.
- Mouse UniProt mapping from
  `https://rest.uniprot.org/uniprotkb/stream?query=organism_id%3A10090&format=tsv&fields=accession%2Cgene_primary%2Creviewed&compressed=true`,
  saved as `uniprot_mouse.tsv.gz` on 2026-09-08. Refreshing this endpoint changes
  the version of the mapping; the frozen hash is authoritative for this build.
- Direct protein-layer and cached fallback source URLs are in
  `results/cosmos_pkn_audit/manifest.json`. These source responses remain in
  `results/cosmos_pkn_audit/raw/` and are dependencies of this build.

New source versions should be archived as a new build before interpreting
changes in coverage. Fetch date is not a substitute for the underlying resource
release date; the individual manifests distinguish them when available.
