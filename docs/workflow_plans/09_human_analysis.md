# Workflow 09: Human Liver Comparison

**Script:** `scripts/workflows/09_human_analysis.py`
**Language:** Python
**Purpose:** Compare the mouse FCG hepatocyte gonadal t-value signature with sex-associated gene expression in human bulk liver RNA-seq from `GSE179250`.

## Source Dependencies

This Python script defines its own helpers and does not source the R helper
stack used by workflows 00-08. It is not self-contained in a public checkout.

It is an optional external-data analysis, not a stage invoked by Workflow 17.
The mouse input is a whole-hepatocyte signature; expanded-marker zonation
inputs and the standalone diet cohort are not upstream inputs.

## Execution and side effects

The collaborator's HPC paths are the defaults. Set `--data-dir` for a local
input tree and `--output` for a fresh destination. `--de-only` runs the three
human differential expression analyses and their applicable plots while
skipping mouse orthology, signature scoring and cross-species ranking. It
still requires the human count, TPM and SOFT inputs. `--help` shows the CLI
after importing the base plotting and statistics libraries. There is no
`--validate-only` mode.

The `pypath` homology lookup in `main()` obtains mouse-to-human mappings from
the installed package's resources/cache; the script does not use the
repository's frozen `support/ortholog_resource.RData`. A suitable cache or
network access and a compatible Python environment are full-run prerequisites.
An empty mapping stops the full run. The package list below is not a version lock.

Python packages imported by the script:

- `anndata`
- `scanpy`
- `decoupler`
- `pydeseq2`
- `pypath`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `adjustText`

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Header and Configuration | 1-69 | Documents external inputs, sets `DATA_DIR`, `RESULTS_DIR`, plotting palette, and volcano thresholds. |
| 2 | I/O Helpers | 73-88 | Save matplotlib/seaborn figures to `RESULTS_DIR`. |
| 3 | Metadata and AnnData Helpers | 91-125 | Parse GEO SOFT metadata and build normalized `AnnData` with raw count and TPM layers. |
| 4 | Signature and Plot Helpers | 129-265 | Translate mouse signature genes to human orthologs, plot ULM boxplots, run ULM scoring, run PyDESeq2 DE, and prepare volcano tables. |
| 5 | Volcano and Dot Plot Helpers | 269-545 | Python plotting helpers for DE volcanoes and concordant mouse-human ranked-gene dot plots. |
| 6 | Main Load and QC | 549-592 | Load human count/TPM matrices, load mouse signature, parse sample metadata, build `AnnData`, and plot PCA by sex. |
| 7 | Signature Scoring | 595-616 | Translate the mouse signature, run decoupler ULM, save estimates, and plot ages 12-49 signature enrichment by sex. |
| 8 | Human Sex DE | 619-670 | Run male-vs-female PyDESeq2 across all ages and ages 12-49; write statistic matrices/full DE tables and volcano plots. |
| 9 | Mouse-Human Concordance | 673-698 | Merge mouse DE/signature results with human DE and plot top concordant female- and male-biased genes. |
| 10 | Female Age Contrast | 700-718 | Run old-vs-young DE among female samples and plot signature enrichment by age category. |

## Input Files

External files expected under `DATA_DIR`:

- `bulk/zhou/GSE179250_gene_reads.txt` - Human liver gene-level count matrix.
- `bulk/zhou/GSE179250_gene_tpm_normalized.txt` - Human liver TPM-normalized matrix.
- `bulk/zhou/GSE179250_family.soft.gz` - GEO SOFT metadata for sample sex, age, and titles.
- `hepa_gonad_signature.csv` - Copy of the mouse hepatocyte gonadal signature.
- `male_vs_female_DA_res.csv` - Separate mouse DE table used for the concordant
  ranked-gene dot plots, with `ID`, `t_Hepatocytes`, `adj.P.Val` and `logFC`
  columns. Workflow 09 does not create this input; it must be supplied explicitly.

Repository source for the signature:

- `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/hepa_gonad_signature.csv` - Retained two-column table
  (`ID`, `t_Hepatocytes`) containing the mouse hepatocyte Testes vs Ovaries
  t-value vector. The script reads an external copy, so merely running
  Workflow 02 does not provision `DATA_DIR/hepa_gonad_signature.csv`.

The human liver data come from GEO accession [`GSE179250`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179250), associated with Zhou et al. 2022, "A resource for integrated genomic analysis of the human liver."

## Output Files

Outputs are written to the external `RESULTS_DIR` by default or to the supplied
`--output`. Retained tables and figures under `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/human/` are saved study
artifacts. The explicit output option allows a parallel replay and comparison.

- `adata_processed.h5ad`
- `pca_by_sex.pdf`
- `ulm_sex_signature_estimates.csv`
- `ulm_enrichment_12_49.pdf`
- `deg_zhou.csv`
- `full_results_deg_zhou_malevsfemale.csv`
- `volcano_male_vs_female_highlight.pdf`
- `12_49_malevsfemale_deg_zhou.csv`
- `12_49_malevsfemale_full_results_deg_zhou.csv`
- `volcano_male_vs_female_highlight_1249.pdf`
- `top_ranked_genes_female_dot.pdf`
- `top_ranked_genes_male_dot.pdf`
- `olderthan11female_oldvsyoung_full_results_deg_zhou.csv`
- `olderthan11female_oldvsyoung_deg_statmatrix_zhou.csv`
- `ulm_enrichment_female__old_vs_young.pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `DATA_DIR` | `pathlib.Path` | External directory holding the large human inputs and copied mouse signature files. |
| `RESULTS_DIR` | `pathlib.Path` | External output directory for tables, `.h5ad`, and PDF figures. |
| `SOFT_PATH` | `pathlib.Path` | `GSE179250_family.soft.gz` metadata path. |
| `GENDER_PALETTE` | `dict` | Female/male plotting colors. |
| `bulk_data`, `tpm_bulk_data` | `pandas.DataFrame` | Human liver count and TPM matrices. |
| `meta_df` | `pandas.DataFrame` | Parsed GEO metadata aligned to count-matrix columns. |
| `adata` | `anndata.AnnData` | Human liver expression object with raw counts and normalized layers. |
| `signature` | `pandas.DataFrame` | Mouse hepatocyte gonad signature translated to human ortholog symbols for ULM scoring. |
| `estimate` | `pandas.DataFrame` | decoupler ULM signature estimates per human sample. |
| `results_df`, `results_df_sub` | `pandas.DataFrame` | PyDESeq2 male-vs-female results across all ages and ages 12-49. |

## Coherence Assessment

**Status: documented optional external analysis; not a public-checkout preflight.**
The workflow is isolated from the R rerun path. The human data, both copied mouse
tables, external paths, homology resources and compatible packages must be
provisioned before a full run. No analysis was rerun during the documentation
audit on 18 November 2026.

## Notes for Context Restoration

1. This workflow should be treated as optional validation/comparison, not a dependency of workflows 00-08.
2. The script defaults to collaborator HPC paths. Supply both path overrides for a local replay.
3. `results/current/02_gonadal_sex_and_hepatocyte_gene_expression/pseudobulk/hepa_gonad_signature.csv` is checked in as the stable seed for the mouse hepatocyte signature; the script reads a copy from external `DATA_DIR`.
4. The human DE contrast is `male` versus `female` with `female` as the reference in `calc_deg()`.
5. The ages 12-49 subset is used for one ULM enrichment plot and one human male-vs-female DE rerun.
