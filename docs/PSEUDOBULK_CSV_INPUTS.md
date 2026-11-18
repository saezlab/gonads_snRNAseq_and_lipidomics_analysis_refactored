# Pseudobulk CSV inputs

## Commit policy

RNA input data are distributed as **CSV parts strictly smaller than 50,000,000 bytes**. Single-cell/Seurat input objects are not included in the commit. The original RDS objects stay local as source/provenance records; fitted model and causal-network states are a separate artifact category.

These exports contain aggregated profiles, not individual nuclei. A profile is a mouse-by-zone or mouse-by-cell-type pseudobulk. Exporting the pseudobulks does not provide per-nucleus assignments or recover the missing composition QC for Figure4A/S4B–C.

The default readers prefer the verified CSV manifests even when local original objects exist. An original RDS can still be supplied explicitly for local source verification or historical replay. A fresh scientific checkout does not need those RNA RDS objects.

## Where to find them

All seven exports live under `data/RNA/pseudobulk_csv/`. The release registry maps each `<role>_transport` key to that export's `manifest.json`.

| Role / folder | Aggregated profiles | Scientific use |
| --- | --- | --- |
| `current_fcg_zonation` | 86 | Expanded-marker FCG, 43 paired mice; current Figure5/zonation branch |
| `current_gdx_zonation` | 44 | Expanded-marker sham/GDX, 22 paired mice; supporting analysis |
| `current_wdsd_zonation` | 56 | Expanded-marker standard/Western-diet cohort; current standalone Figure6 |
| `wdsd_hepatocytes` | 28 | Whole-hepatocyte profiles in the standard/Western-diet cohort |
| `historical_wdsd_zonation` | 56 | Original zone definition, retained for frozen comparison workflows |
| `wd_cell_types` | 40 | Cell-type pseudobulks from the WD comparison cohort; this object is not hepatocyte-only |
| `wd_zonation` | 16 | Original zone profiles from the WD comparison cohort |

The current and historical diet-cohort zonation exports are different inputs. Never substitute one for the other to satisfy a missing-file check. All seven source objects have 20,123 genes; retained genes in a downstream model are determined by that model's documented eligibility rules.

## Contents and matrix orientation

Each export has a JSON manifest, a CSV file index, ordered CSV parts for the supplied `counts`, normalized `data`, and `scale.data` layers, and aggregate-profile metadata/identity CSVs. The manifests specify filenames, dimensions, ordering, types and hashes; use them rather than globbing and guessing part order.

Numeric matrices are genes × profiles. The first column identifies the gene and the remaining columns identify profiles. Parts partition profiles/columns, retaining the same ordered gene list. Small layers can have a single numbered part; larger layers are split as needed to meet the strict byte limit. Full metadata is preserved at the aggregate-profile level.

Export is a format conversion, not preprocessing. Counts, supplied normalization, scaling, gene/profile identifiers and pairing are retained. The pipeline does not recalculate zonation assignments or normalize the RNA again merely because CSV is used.

## Provenance and verification

The original source RDS path, size and SHA256 remain distinct from the CSV manifest/payload hashes. Numerical round-trip checks establish the export's relation to that source. The source RDS's hash is never relabelled as a CSV file hash.

New analyses record the actual CSV manifest and payload files they read. Approved older caches retain their original input identities. Where a frozen cache names an excluded RDS, the reader can use a verified export bridge: the source SHA256 must match the recorded original identity, and every CSV payload must pass verification. This does not rewrite the approved cache or assert that the old analysis was originally run from CSV.

Use [the public input index](provenance/current_input_manifest.csv) for current distributed files and [the original local input ledger](provenance/current_input_manifest_20261116.csv) for the previous local-object staging record. Detailed source identity and content checks are in each export manifest. Source objects and approved scientific results remain unchanged.

Run these checks from the repository root:

```sh
python3 scripts/maintenance/verify_repository_organization.py
python3 scripts/maintenance/audit_github_readiness.py --strict
Rscript scripts/workflows/17_expanded_marker_zonation_update.R --validate-only
```

The commit audit rejects serialized RNA input objects, incomplete/hash-inconsistent CSV exports, and RNA CSV files at or above 50,000,000 bytes. Maintenance checks do not stage, commit or push.

[Final distribution and Git-only checkout verification](provenance/pseudobulk_csv_distribution_20261118.json) records the policy gates and source-RDS-free smoke test.

## Local re-export and direct reading

`Rscript scripts/maintenance/export_pseudobulk_csv.R --validate-only` compares the existing exports against all seven original local objects without modifying them or running scientific models. Re-exporting or this source comparison requires the local originals; ordinary CSV reading and scientific default validation do not.

```r
source("scripts/config/analysis_release.R")
source("scripts/core/pseudobulk_csv.R")
x <- pb_components(analysis_path("current_wdsd_zonation_transport"))
# x$layers contains the three ordered gene-by-profile matrices.
# x$meta.data contains aggregate-profile metadata.
```

The exporter tests forced multi-part splits and difficult floating-point/metadata cases in addition to the actual datasets. CSV values were checked against source objects, and an independent Python audit checks numerical content hashes, file hashes, order, dimensions and metadata joins. See `docs/provenance/pseudobulk_csv_roundtrip.json` and the GitHub-readiness audit outputs.
