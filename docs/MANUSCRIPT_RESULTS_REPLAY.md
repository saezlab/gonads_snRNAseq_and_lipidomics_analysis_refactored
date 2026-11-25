# Parallel manuscript results

The standalone replay writes a separate manuscript-shaped tree at `results/manuscript/`. The [release registry](../scripts/config/analysis_release.json) continues to select the approved scientific release. The replay writes neither in that release nor in the existing workflow result roots.

| Results subsection | Folder under `results/manuscript/` | Workflow results |
| --- | --- | --- |
| Sexual dimorphism in hepatic lipid homeostasis is driven by gonadal sex | `01_gonadal_sex_and_hepatic_lipid_homeostasis/` | 03 lipid DE and gonadectomy, 06 structural and 08 variability |
| Sexual dimorphism in hepatocyte gene expression is driven by gonadal sex | `02_gonadal_sex_and_hepatocyte_gene_expression/` | 01 QC and 02 whole-cell-type pseudobulk analyses |
| Gonadal sex-biased regulatory circuitry controlling lipid homeostasis | `03_gonadal_sex_biased_regulatory_circuitry/` | 04 MOFA, global 07 knockdown validation and 11 in-vivo KO |
| Distinct gonadal sex biases in pericentral and periportal hepatocytes | `04_pericentral_and_periportal_hepatocyte_bias/` | 17 paired gonadectomy and the gonad-versus-zone panel/source lookup |
| Gene-lipid regulatory interactions with gonadal sex-polarised zonation bias | `05_sex_polarised_zonation_gene_lipid_interactions/` | 17 current FCG zonation, zonal 07 validation and shared zonation table exports |
| Sex-specific impact of chronic Western diet on gene-lipid regulatory circuitry | `06_western_diet_and_gene_lipid_circuitry/` | 17 invokes 14, 15, 16, both anchor layouts and the Figure 6B/6E manuscript renderer |

## Running the replay

Run from the repository root in the analysis environment, with Python 3.11 or later and the workflow R dependencies available.

```sh
python3 scripts/maintenance/replay_manuscript_results.py --scratch results/runs/NAME
```

Use `--output results/runs/NAME_results` for another complete run. Use `--workflow 02 03 04 06 07 17` to select a subset in the specified order. Each workflow checks its immediate result inputs and requires a fresh destination for the results it produces. A workflow that fails before writing results can be retried in the same scratch directory, with the next attempt recorded separately.

The runner builds a private working directory with links to the current inputs and code. Within it, the normal published workflow paths point at the fresh subsection outputs. For example, the internal `results/current/01_gonadal_sex_and_hepatic_lipid_homeostasis/lipidomic/ttop_full_list.RData` resolves to the new lipid subsection, and `results/current/03_gonadal_sex_biased_regulatory_circuitry/in_vitro_validation/knockdown_t_table.csv` resolves to the new global KD output. The numbered diet handoffs retain their own adjacent provenance.

Workflow01 writes isolated RNA exports. The shipped split CSVs and metadata remain the inputs to the downstream whole-cell workflows, so a different local RNA object cannot silently replace the release inputs. The replay records the re-export checksums for comparison. Workflow04 uses the exact retained 5–15 factor HDF5 set in read-only reuse mode and verifies the files' checksums after execution.

Workflow00 documents the original external preprocessing. Human09 requires its GEO inputs, mouse tables and orthology resources. Xenium10 requires its collaborator imaging inputs and runtime. Check their plans and retain sourced independent assay exports where the raw analysis cannot be replayed. Figures assembled outside the repository should be checked against the newly generated component panels and tables.

## Local comparison

Execution paths, commands and logs are under the specified scratch directory. Keep the scratch path available while using model caches that record its input routes. Store private manuscript evidence, figure extracts and editing reports in ignored paths. The ignored `results/legacy/` tree holds the previous result directories and their original provenance. The existing repository history and relocation records remain the recovery source.

The publication tree is `results/current/`. The promotion ledger records which verified artifacts were selected from the release and the replay. For a future promotion, compare the new outputs and update the registry and manifests together. A file already tracked by Git stays tracked regardless of an ignore rule. Review existing tracked result paths during promotion and retain any frozen inputs still used by active workflows. Run the organization verifier and strict GitHub readiness audit before staging.
