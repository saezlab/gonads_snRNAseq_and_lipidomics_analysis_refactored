# WP1 and WP2 frozen input files

These repository-local copies make the October 2026 WP1 and WP2 workflows
independent of the sibling `original_datasets` and `Revised_Figures` directories.
Original shared archive files were preserved. The copies retain their original
bytes; no normalization, imputation, Seurat update or annotation rewriting was
performed during import.

| Repository file | Role |
|---|---|
| `metabolomics/cyp2c38ko/lipidomics_cyp2c38ko.csv` | WP1 integrated Cyp2c38 KO lipidomics export, including biological and analytical-QC columns |
| `metabolomics/wdsd/lipidomics_cf453.csv` | WP2 paired CF453 SD/WD lipidomics and annotation export |
| `RNA/wdsd/wdsd_hep_pseudobulks.rds` | WP2 updated whole-hepatocyte pseudobulk Seurat object; 28 mice |
| `RNA/wdsd/wdsd_zonation_pseudobulks.rds` | WP2 updated zonation pseudobulk Seurat object; 56 paired zone profiles |
| `metabolomics/fcg/lipidomics_fcg.csv` | Original FCG raw lipid annotation export used to audit strict and broader molecular matching |
| `../docs/reference_figures/FCGL_Fig3.ai` | Frozen revised Figure 3 selection/provenance; optional for numerical network replay |

`wp1_wp2_input_migration_manifest.csv` records source and destination paths,
byte counts and independently verified SHA256 hashes. All six copies total
30,965,755 bytes. Source paths identify the original archive; workflows read the
repository destinations. Run `Rscript tests/test_wp1_wp2_local_inputs.R` from the
repository root to verify these bytes without access to the original archive.

The existing FCG matrices under `RNA/` and
`metabolomics/Lipid_data_file_JPC_cleaned.csv`, sample metadata and original
network exports under `results/mofa/` and `results/zonation/` remain the frozen
reference inputs. The raw FCG annotation CSV supplements the cleaned feature
matrix and does not replace it. Only the two WDSD input RDS files receive a
narrow `.gitignore` exception; the general generated-RDS rule remains intact.

CF459 WD-only inputs are now separately imported for Workflow 13/WP3; see
[WP3_INPUTS.md](WP3_INPUTS.md). Other datasets outside these work packages
remain in their existing archive locations. Collaborator statistical result exports
are not runtime inputs: WP1 models are refitted from the integrated KO export.
