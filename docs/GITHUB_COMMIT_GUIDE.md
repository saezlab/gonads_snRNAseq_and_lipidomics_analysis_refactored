# Preparing the repository for GitHub

The seven public pseudobulk transports under `data/RNA/pseudobulk_csv/` contain
split CSV matrices and typed metadata/identity manifests. They include the
current expanded-marker FCG/GDX/diet inputs and the separate frozen inputs
needed by the documented comparison branches. The approved dated release remains in
`results/zonation_marker_update_20261116/`. Historical files moved to
`miscellaneous/` are retained locally with migration/provenance records;
moving an artifact does not authorize changing its scientific content.

This guide prepares a reviewable working tree. It does not authorize or perform
a Git commit, push, release upload, history rewrite, or deletion of local data.

## What belongs in the source commit

- Active workflows, shared functions, meaningful tests, configuration and guides.
- The registered CSV transports under `data/RNA/pseudobulk_csv/`, corresponding
  metadata/identity CSVs and JSON/CSV manifests, plus retained small preparation
  and compatibility files. Every RNA input CSV must be **strictly smaller than
  50,000,000 bytes** (decimal MB, not MiB). The transport preserves source
  counts/data/scale.data matrices, gene/profile order and metadata types without
  committing a Seurat object. Current and historical zonation definitions remain
  separate. Whole-hepatocyte data, lipid data and independent experiments retain
  their own provenance.
- Curated current-release vector panels, compact panel source tables,
  supplementary-table exports, numerical QA summaries, methods and manifests.
- Archived source code, documentation and migration/checksum manifests. These
  explain historical analyses even when their bulky generated files are not
  distributed in the Git tree.
- Existing active COSMOS/MOON material. Do not remove frozen causal-model states
  or change their line endings as a side effect of zonation cleanup.

The ignore rules keep local libraries, tool settings, histories, caches, render
images, replay/model/bootstrap states, runtime logs, delivery ZIPs, duplicate
input copies and selected bulk comparison exports out of an ordinary commit.
Fresh reproduction directories under `results/runs/` are generated and ignored until intentionally promoted and curated as a new release. The approved release remains complete on disk. `.gitignore` changes do not
remove files already tracked by Git; any tracked-file relocation must be
reviewed explicitly.

Single-cell/Seurat/pseudobulk RDS, RData, HDF5/H5AD and related RNA object formats
are excluded from ordinary commits, including all seven former RNA RDS
allowlist entries. The original objects remain on disk as local source evidence.
This does not exclude the separately declared causal/model RDS reference states
outside `data/RNA/`. The audit also flags obvious single-cell containers
elsewhere under `data/`; a lipid-only RData table is classified separately.

`miscellaneous/` is intentionally not ignored as a whole. Its source code,
documentation and migration manifests remain visible. `older_code/`,
`documentation/`, and source files in `refactoring/` and `quick_tests/` remain
available. `older_results/`, `older_inputs/`, and `comparisons/` exclude their
payloads while allowing Markdown, compact README/session records, manifest
CSV/JSON, provenance CSV and SHA256 CSV records.
Scoped `logs/`, `local/` and generated-output directories retain
machine-specific material locally. Generated quick-test figures and logs are
also excluded.
There is no repository-wide `**/legacy/` rule hiding archived code.

An explicit exception preserves 13 small archived files (4,074,361 bytes total)
still required by active MOON code: full-precision Factor 1 inputs, mapping and
observation companions, and selected reference score/compression states. Their
exact paths, roles and SHA256 hashes are in
[active_frozen_dependencies.json](provenance/active_frozen_dependencies.json).
Only those payloads are unignored; the remaining archived discovery/search
states stay local. They are independent of the expanded-marker zonation data.

## Frozen bytes and Git line endings

The scoped `.gitattributes` exceptions preserve byte-level checksums for frozen
CSV inputs, archival engine/source copies, scientific result snapshots and
provenance ledgers. This includes newly imported lipid data and the source
manifests used by validation/replay. Ordinary active source files retain the
repository's LF text convention. Verify staged blobs, not only working files:
normalizing a frozen CSV during `git add` would change its SHA256 even when
its parsed numerical values remain identical.

## Run the readiness audit

From the repository root, using any supported system Python3:

```sh
python3 scripts/maintenance/audit_github_readiness.py --strict
```

The audit uses the Python standard library and read-only Git commands. It does
not fit models, run pipelines, stage changes or contact a remote service. Its
default ignored output directory is
`miscellaneous/local/github_readiness/`, containing:

- `github_readiness.json` and `README.md`: totals, gates and interpretation.
- `candidate_files.csv`: every present tracked or nonignored untracked file,
  its byte size and Git state. This approximates the tree a broad `git add`
  would include, before staging anything.
- `checks.csv`: file-size, CSV-only RNA inputs, transport completeness and unchanged-index
  checks. `--strict` returns a nonzero exit status if a gate fails.
- `active_frozen_dependencies.csv`: inclusion, byte-size and SHA256 checks for
  every declared archived dependency of active MOON workflows.
- `pseudobulk_transport_files.csv` and `pseudobulk_transport_layers.csv`:
  descriptor/file hashes, strict size bounds, layer/chunk dimensions and IDs.
  `pseudobulk_transport_issues.csv` identifies missing or inconsistent files.
  Every chunk is parsed to check numeric finiteness, exact gene/profile order,
  contiguous column coverage, metadata joins and counts/data/scale.data
  completeness. Numerical matrix hashes are independently recomputed from
  IEEE754 binary64 values in little-endian column-major order, so a changed
  value cannot pass merely because a CSV descriptor was updated. Source-object
  provenance hashes remain recorded without requiring the local RDS in a
  public clone.
- `single_cell_input_candidates.csv` scans the entire candidate tree for
  obvious single-cell containers and barcode/feature sparse-matrix bundles.
  `hdf5_container_review.csv` distinguishes retained fitted MOFA models from
  single-cell input formats; the audit does not load large fitted models.
- `nonshipping_markdown_links.csv`: local document links whose targets would
  not ship, including missing targets and intentionally local-only outputs.
- `nonshipping_manifest_references.csv`: provenance paths not included in the
  candidate tree. These require interpretation; historical absolute paths are
  evidence of the original run, not portable runtime configuration.
- `machine_specific_source_literals.csv`: source lines containing common
  local/HPC path prefixes, for review rather than automatic replacement.

The audit also scans reachable Git history for oversized blobs. Existing
tracked scientific data and history can make the repository large even after
new generated files are excluded. Its byte totals are logical file sizes, not
an estimate of the compressed Git transfer size. The audit records missing
documentation links and nonshipping provenance separately from the size/cache
gates; inspect these ledgers before considering the public clone complete.

To regenerate the public archive-recovery ledger without writing Git objects:

```sh
python3 scripts/maintenance/audit_github_readiness.py --strict --write-archive-recovery
```

[archive_recovery.csv](provenance/archive_recovery.csv) records old/new paths,
fresh working-file SHA256 checks and the pre-cleanup HEAD/blob. A row is marked
`git_history_byte_identical` only when `git hash-object --no-filters` for the
archived working file exactly matches its original HEAD blob. Retrieve a
verified historical file using the row's `git_object_spec`, for example:

```sh
git show '<pre_cleanup_HEAD>:<old_path>' > /a/new/recovery/location
```

Create the destination directory first and choose a new destination to avoid
overwriting local files. Recovery requires the recorded commit to be present;
a shallow clone may need additional history. Regeneration retains the existing
ledger's pre-cleanup commit rather than substituting a later HEAD. Previously
untracked or byte-different files are
marked `local_archive_required`; Git history is not claimed to contain those
exact working bytes. Keep/deposit their local archive with its recorded hashes.
The command writes only the recovery CSV and local audit outputs; it neither
stages files nor rewrites history.

The repository's RNA input limit of less than 50,000,000 bytes is stricter than
GitHub's regular file limits. It does not apply to unrelated existing result
tables or frozen model resources. GitHub warns about regular Git files above
50 MiB and blocks files above 100 MiB.
Large downloadable artifacts can instead be distributed as release assets or
through a scientific data repository, with their URLs, sizes and SHA256 hashes
recorded in the source tree. Existing ordinary Git history need not be rewritten
merely because some valid files exceed the warning threshold. See
[GitHub's large-file guidance](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).

For this repository, the full historical CF453 edge CSVs and downloaded STITCH
databases are excluded explicitly because they exceed the ordinary Git limit.
Their import/generation code and source/checksum provenance remain relevant.
The audit will catch any future oversized candidate rather than assuming an
extension such as CSV or GZ is safe.

## Public-clone dependencies and private manuscript material

Run the scientific workflows in their documented dependency order. A model
cache excluded from Git can be regenerated from the canonical input; a plotting
or report-only entry that consumes that cache requires the corresponding
analysis stage first. Read the matching workflow plan before changing or
reexecuting it. A curated Git checkout is not the complete historical output
bundle.

Selected MOON replay additionally enforces the package versions and cosmosR
commit pinned in its configuration. The readiness file audit does not load
those packages or promise that the current runtime satisfies that contract.
Large upstream source databases and optional external KO culture counts remain
separately documented restoration requirements; the public tree does not claim
to contain every raw input for every optional extension.

The working Word manuscript and its embedded reference figures are private,
optional authoring inputs. They are not prerequisites for the scientific
pseudobulk fits. Regenerating the manuscript paragraph inventory or the exact
author-facing replacement report requires the appropriate local DOCX. The full
source extraction CSVs (`manuscript_paragraph_locations.csv`,
`manuscript_zonation_context_scan.csv`, `manuscript_numeric_inventory.csv`) and
`zonation_update_report.html` embed original manuscript wording and remain
local-only; the ignore rules do not delete or modify them. The reviewed public
PDF, `new_replacements.csv` numeric ledger and portable `current_panels.md`
index remain available for the commit. Use those public artifacts when the
private input is unavailable. The scientific report-source preparation supports
`--tables-only` without the DOCX. Do not claim that a public checkout can
reproduce an optional manuscript-audit command requiring an undistributed
manuscript. Scientific defaults and private report steps remain separate.

Raw HPC preprocessing, human-cohort and collaborator imaging workflows can
still require external datasets. Document those prerequisites rather than
inventing substitutes or silently treating an absolute-path provenance field
as an input present in the public repository. The provided expanded-marker
preprocessing scripts establish an algorithm; exact retained-marker/runtime
confirmation remains a separate methods-provenance item.

The repository currently documents package requirements and records per-run
session information; a complete cross-language dependency lock is not implied.
Use those version records when reproducing frozen results. Add an intentional
tested environment specification if environment reconstruction is in scope;
do not commit the local R library or a machine-specific Python installation.

## Final review before the user commits

```sh
git status --short
git diff --stat
git diff --check
python3 scripts/maintenance/audit_github_readiness.py --strict
git add --dry-run .
```

Review the actual intended file list, especially relocations of tracked files.
Confirm that canonical inputs and archived source/provenance are present,
that current release artifacts retain their recorded hashes, and that a broad
add would not include local settings or oversized generated data. Keep the
full local archive and distribution bundles until the independently reviewed
source commit and any external artifact deposit have been verified. Commit,
push and release publication remain separate user actions.
