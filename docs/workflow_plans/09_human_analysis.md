# Workflow 09: Human Liver Comparison

**Script:** `scripts/workflows/09_human_analysis.py`
**Language:** Python
**Purpose:** Compare the mouse FCG hepatocyte gonadal t-value signature with sex-associated gene expression in human bulk liver RNA-seq from `GSE179250`.

## Source Dependencies

This workflow is standalone Python and does not source the R helper stack used by workflows 00-08.

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
- `male_vs_female_DA_res.csv` - Mouse DE/signature table used for the concordant ranked-gene dot plots.

Repository source for the signature:

- `results/pseudobulk/hepa_gonad_signature.csv` - Two-column table (`ID`, `t_Hepatocytes`) containing the mouse hepatocyte Testes vs Ovaries t-value vector. This is a near-identical earlier export of the current workflow 02 hepatocyte gonad t-values.

The human liver data come from GEO accession [`GSE179250`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179250), associated with Zhou et al. 2022, "A resource for integrated genomic analysis of the human liver."

## Output Files

Outputs are written to the external `RESULTS_DIR` configured in the script, not to this repository by default.

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

**Status: GOOD, external-data dependent.** The workflow is documented and isolated from the R rerun path. It cannot be rerun from a standard checkout unless the large `GSE179250` human liver files and copied mouse input tables are present under the script's external `DATA_DIR`.

## Notes for Context Restoration

1. This workflow should be treated as optional validation/comparison, not a dependency of workflows 00-08.
2. The script currently uses absolute collaborator/HPC paths for `DATA_DIR` and `RESULTS_DIR`.
3. `results/pseudobulk/hepa_gonad_signature.csv` is checked in as the stable seed for the mouse hepatocyte signature; the script reads a copy from external `DATA_DIR`.
4. The human DE contrast is `male` versus `female` with `female` as the reference in `calc_deg()`.
5. The ages 12-49 subset is used for one ULM enrichment plot and one human male-vs-female DE rerun.
