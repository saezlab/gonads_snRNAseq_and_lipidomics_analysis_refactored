# Workflow Testing Guide

This guide explains how to systematically test all refactored workflows to ensure they work correctly and produce results consistent with the manuscript.

---

## Quick Start

### Option 1: Interactive Testing (Recommended)

Run the interactive test script that checks each workflow one by one:

```r
source("tests/run_workflow_tests.R")
```

This will:
- Run each workflow sequentially
- Check that expected output files are created
- Validate key results against manuscript findings
- Provide a summary report

### Option 2: Comprehensive Testing (Advanced)

Run the full test suite using `testthat`:

```r
source("tests/test_all_workflows.R")
```

Set configuration options at the top of the script:
- `RUN_FULL_WORKFLOWS = TRUE` to run all workflows (slow, ~30-60 min)
- `RUN_FULL_WORKFLOWS = FALSE` to run quick validation tests only

---

## Testing Individual Workflows

### Workflow 1: Data Preparation & QC

**Script:** `workflows/01_data_preparation_qc.R`

**What it does:**
- Generates split RNA data files from Seurat object
- Creates metadata files
- Produces QC reports

**Expected outputs:**
```
data/RNA/scaled_RNA_1.csv          # RNA data part 1 (gene + 79 samples)
data/RNA/scaled_RNA_2.csv          # RNA data part 2 (gene + 80 samples)
data/RNA/scaled_RNA_3.csv          # RNA data part 3 (gene + 56 samples)
data/RNA/meta_data_df.csv          # Metadata (215 samples × variables)
results/qc/data_summary.txt        # QC summary report
```

**Key validations:**
- ✓ RNA data has 20,123 genes
- ✓ RNA data has 215 pseudobulk samples
- ✓ Metadata has required columns: `gonadal_type`, `allosome_type`, `group`, `cell_type`
- ✓ Sample names match between RNA data and metadata

**Test manually:**
```r
source("workflows/01_data_preparation_qc.R")

# Check outputs
file.exists("data/RNA/scaled_RNA_1.csv")
file.exists("data/RNA/meta_data_df.csv")

# Load and validate
source("scripts/core/data_loading.R")
rna <- load_rna_data(verbose = TRUE)
meta <- load_metadata(rna_metadata_file = "data/RNA/meta_data_df.csv", verbose = TRUE)

dim(rna)  # Should be: 20123 x 215
dim(meta) # Should be: 215 x n_columns
```

---

### Workflow 2: Pseudobulk RNA-seq Analysis

**Script:** `workflows/02_fcg_pseudobulk_analysis.R`

**What it does:**
- Differential expression analysis across 6 contrasts and 5 cell types
- PCA visualization per cell type
- Volcano plots for each contrast × cell type
- Pathway enrichment analysis

**Expected outputs:**
```
results/pseudobulk/ttop_full_list.RData             # All DE results (6 contrasts × 5 cell types)
results/pseudobulk/tval_full_df.RData               # t-statistic matrix
results/pseudobulk/pca/pca_*.pdf                    # PCA plots (5 cell types)
results/pseudobulk/volcanos/volcano_*.pdf           # Volcano plots (30 total: 6 × 5)
results/pseudobulk/pathways/pathway_heatmap_*.pdf   # Pathway heatmaps
```

**Key validations:**

✓ **All 6 contrasts analyzed:**
1. `gonadal_type` (Testes vs Ovaries) - PRIMARY
2. `allosome_type` (XY vs XX)
3. `XXT_vs_XXO` (gonadal effect in XX background)
4. `XYT_vs_XYO` (gonadal effect in XY background)
5. `XYO_vs_XXO` (chromosome effect in ovaries background)
6. `XYT_vs_XXT` (chromosome effect in testes background)

✓ **All 5 cell types analyzed:**
- Hepatocytes
- Endothelial cells
- Immune cells
- Stellate cells
- Cholangiocytes

✓ **Expected gene counts (Hepatocytes gonadal_type):**
- Total genes: 20,123
- Significant DEGs (FDR ≤ 0.05): **~4,000** (manuscript reports >4,000)
- Range: 3,000-5,000 is acceptable

✓ **Minimal chromosome effect:**
- Allosome_type should have few DEGs (<100 in hepatocytes)
- Mainly Y-linked genes and X-escapees

**Test manually:**
```r
source("workflows/02_fcg_pseudobulk_analysis.R")

# Load results
load("results/pseudobulk/ttop_full_list.RData")

# Check structure
names(ttop_full_list)  # Should have all 6 contrasts
names(ttop_full_list$gonadal_type)  # Should have all 5 cell types

# Check hepatocyte gonadal_type (PRIMARY FINDING)
hepa_gonad <- ttop_full_list$gonadal_type$Hepatocytes
nrow(hepa_gonad)  # Should be 20,123
sum(hepa_gonad$adj.P.Val <= 0.05)  # Should be ~4,000

# Check top genes
head(hepa_gonad[order(hepa_gonad$adj.P.Val), ], 20)
# Should include: Cyp2d9, Cyp2b9, Cyp2b13, Cyp2c38, Mup genes, etc.
```

---

### Workflow 3: Lipidomic Analysis

**Script:** `workflows/03_fcg_lipidomic_analysis.R`

**What it does:**
- Differential abundance analysis for lipids
- Lipid class enrichment analysis
- LINEX2 network analysis
- Bubble plots showing structural heterogeneity

**Expected outputs:**
```
results/lipidomic/ttop_full_list.RData             # All lipid DE results
results/lipidomic/gonade_type_ttop.csv             # Gonad-dependent lipids table
results/lipidomic/pca_plot.pdf                     # PCA plot (samples separate by gonad)
results/lipidomic/lipid_class_enrichment.pdf       # Class enrichment lollipop plot
results/lipidomic/significant_lipids_network_SIF.csv  # LINEX2 network edges
results/lipidomic/significant_lipids_network_ATT.csv  # LINEX2 network attributes
```

**Key validations:**

✓ **Gonadal effect (PRIMARY FINDING):**
- Total lipids: 1,386
- Significant (FDR ≤ 0.05): **733** (53.9%)
- Range: 650-850 is acceptable
- Male-biased: **448 lipids** (TG, DG, glycerophospholipids)
- Female-biased: **285 lipids** (SM, Cer, HexCer)

✓ **Chromosome effect (MINIMAL):**
- Significant (FDR ≤ 0.05): **0** (manuscript reports 0)
- Range: 0-50 is acceptable

✓ **Lipid class enrichment:**
- Male-biased classes: TG, DG, PI, PE, LPE, PA, LPC, PG, PC
- Female-biased classes: SM, Cer, HexCer

**Test manually:**
```r
source("workflows/03_fcg_lipidomic_analysis.R")

# Load results
load("results/lipidomic/ttop_full_list.RData")

# Check gonadal_type
gonad_results <- ttop_full_list$gonadal_type
nrow(gonad_results)  # Should be 1,386
sum(gonad_results$adj.P.Val <= 0.05)  # Should be ~733

# Check allosome_type (should be minimal)
allosome_results <- ttop_full_list$allosome_type
sum(allosome_results$adj.P.Val <= 0.05)  # Should be ~0

# Check top lipids
head(gonad_results[order(gonad_results$adj.P.Val), ], 20)
```

---

### Workflow 4: MOFA Integration

**Script:** `workflows/04_mofa_integration.R`

**What it does:**
- Multi-omic factor analysis (RNA + lipidomics)
- Identifies Factor 1 (strongest gonad association)
- Constructs gene-lipid-TF correlation networks
- Identifies key regulatory genes: Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9

**Expected outputs:**
```
results/mofa/main_model_RNA_lipid_weight.RData                          # MOFA model weights
results/mofa/factor_scores.pdf                                          # Factor scores plot
results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata  # Network
```

**Key validations:**

✓ **Factor 1 associations:**
- Strongest association with gonadal type
- Explains most variance in hepatocyte RNA + lipidomics

✓ **Top genes identified:**
- **Cyp2b9** (3rd top female-biased lipid metabolism gene)
- **Cyp2d9** (top male-biased lipid metabolism gene)
- **Cyp2b13** (top female-biased lipid metabolism gene)
- **Cyp2c38** (top 5 female-biased lipid metabolism gene)

✓ **Top TFs identified:**
- **Nfic, Dbp, Esr1** (correlated with Cyp genes)

✓ **Network structure:**
- 94 nodes (lipids, genes, TFs)
- 246 edges
- Cyp2b9, Cyp2d9, Cyp2b13 correlated with multiple lipid classes

**Test manually:**
```r
source("workflows/04_mofa_integration.R")

# Load MOFA results
load("results/mofa/main_model_RNA_lipid_weight.RData")

# Check Factor 1 enrichment for gonad type
# (Implementation depends on MOFA object structure)

# Load network
load("results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata")
```

---

### Workflow 5: Zonation Analysis

**Script:** `workflows/05_zonation_analysis.R`

**What it does:**
- Assigns hepatocytes to pericentral vs. periportal zones
- Differential expression within zonation subtypes
- Identifies zonation-specific gonad-dependent genes

**Expected outputs:**
```
results/pseudobulk/zonation/ttop_pericentral_vs_periportal.RData
results/pseudobulk/zonation/zonation_gonadal_pericentral.RData
results/pseudobulk/zonation/zonation_gonadal_periportal.RData
results/pseudobulk/zonation/zonation_specific_networks/
```

**Key validations:**

✓ **Zonation markers:**
- Pericentral: Cyp2e1+, Glul+
- Periportal: Cyp2f2+, Gls2+, Hal+, Sds+

✓ **Gonad-dependent genes by zonation:**
- **43%** pericentral-biased: Cyp2c38, Cyp2b13, Cyp2d9
- **18%** periportal-biased: Cyp2b9

✓ **Zonation-specific regulators:**
- Female pericentral: Sult3a1, Abcb1a, Cyp2a4
- Male periportal: Mtmr11, Mindy3, Esco2, Kitl

---

### Workflow 6: Structural Analysis

**Script:** `workflows/06_structural_analysis.R`

**What it does:**
- Correlates gene expression with lipid abundance
- Regresses gene-lipid correlations on lipid structural features
- Identifies genes with substrate specificity for chain length/saturation

**Expected outputs:**
```
results/multi-omic/structural_model/hepatocyte_cor_mat_RNA_lip.csv
results/multi-omic/structural_model/lipids_with_structure.csv
results/multi-omic/structural_model/significant_gene_lipid_class_associations.csv
results/multi-omic/structural_model/genes_associated_with_chain_length.csv
results/multi-omic/structural_model/genes_associated_with_db_position.csv
results/multi-omic/structural_model/gene_plots/*.pdf
```

**Key validations:**

✓ **Structural features analyzed:**
- Chain length (C)
- Double bond position (DB)
- Interaction effects (abundance × structure)

✓ **Regression models:**
- Main model: `correlation ~ class + abundance + C + DB + abundance:C + abundance:DB`
- Class-specific models: Same, but within lipid classes (≥10 lipids/class)

✓ **Expected genes with structural associations:**
- Cyp2c38, Cyp2b9, Cyp2d9, etc.
- Genes showing preference for specific chain lengths or saturation levels

---

## Common Issues and Troubleshooting

### Issue 1: "Contrast not found" Error

**Error:**
```
Error in get_contrast(contrast_type) :
  Contrast 'group' not found
```

**Cause:** Trying to use `"group"` as a contrast name instead of a specific FCG comparison.

**Fix:** Use specific contrast names:
- ✓ `"XXT_vs_XXO"` (correct)
- ✗ `"group"` (incorrect)

The workflow should extract `group_variable` from the contrast configuration:
```r
contrast_config <- get_contrast("XXT_vs_XXO")
group_variable <- contrast_config$group_variable  # Returns "group"
```

### Issue 2: "undefined columns selected" Error

**Error:**
```
Error in `[.data.frame`(targets, , c("sample", group_variable)) :
  undefined columns selected
```

**Cause:** Using contrast name as `group_variable` instead of the actual metadata column.

**Fix:** Same as Issue 1 - extract `group_variable` from contrast config.

### Issue 3: Zero Genes in RNA Data

**Error:**
```
RNA data loaded successfully:
  Genes: 0
  Samples: 645
```

**Cause:** Split RNA files have same columns (samples) but different rows (genes).

**Fix:** Regenerate split files:
```r
source("scripts/seurat_to_df.R")
```

This should create files with:
- File 1: gene + samples 1-79
- File 2: gene + samples 80-159
- File 3: gene + samples 160-215

### Issue 4: Missing Expected DEGs

**Symptom:** Hepatocyte gonadal_type has <3,000 or >5,000 DEGs

**Possible causes:**
- Incorrect contrast vector
- Wrong group_variable used
- Metadata mismatch
- Missing samples

**Debug:**
```r
source("scripts/core/data_loading.R")
rna <- load_rna_data()
meta <- load_metadata(rna_metadata_file = "data/RNA/meta_data_df.csv")

# Check samples match
colnames(rna) %in% meta$sample
sum(colnames(rna) %in% meta$sample)  # Should be 215

# Check grouping variable
table(meta$gonadal_type)  # Should show Ovaries and Testes
table(meta$allosome_type) # Should show XX and XY
table(meta$group)         # Should show XXO, XXT, XYO, XYT
```

---

## Expected Runtime

Typical runtimes on a standard workstation:

| Workflow | Runtime | Memory |
|----------|---------|--------|
| 01 - Data Preparation | 5-10 min | 8 GB |
| 02 - Pseudobulk Analysis | 10-20 min | 16 GB |
| 03 - Lipidomic Analysis | 5-10 min | 8 GB |
| 04 - MOFA Integration | 15-30 min | 16 GB |
| 05 - Zonation Analysis | 15-25 min | 16 GB |
| 06 - Structural Analysis | 20-40 min | 16 GB |
| **Total** | **70-135 min** | **16 GB** |

---

## Validation Checklist

Use this checklist to verify that all workflows are working correctly:

### Data Preparation
- [ ] RNA data has 20,123 genes × 215 samples
- [ ] Metadata has 215 samples with required columns
- [ ] Sample names match between RNA and metadata

### Pseudobulk Analysis
- [ ] All 6 contrasts analyzed
- [ ] All 5 cell types analyzed
- [ ] Hepatocyte gonadal_type: 3,000-5,000 DEGs (FDR ≤ 0.05)
- [ ] Allosome_type has minimal DEGs (<100)
- [ ] 30 volcano plots generated (6 contrasts × 5 cell types)
- [ ] PCA plots show gonad separation (not chromosome)

### Lipidomic Analysis
- [ ] Gonadal_type: 650-850 significant lipids (expected ~733)
- [ ] Allosome_type: 0-50 significant lipids (expected 0)
- [ ] PCA shows gonad clustering
- [ ] Male-biased: TG, DG, glycerophospholipids
- [ ] Female-biased: SM, Cer, HexCer

### MOFA Integration
- [ ] Factor 1 associates with gonadal type
- [ ] Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9 emerge as top genes
- [ ] Nfic, Dbp, Esr1 identified as top TFs
- [ ] Gene-lipid network has ~94 nodes, ~246 edges

### Zonation Analysis
- [ ] 43% of gonad-dependent genes are pericentral-biased
- [ ] 18% of gonad-dependent genes are periportal-biased
- [ ] Cyp2c38, Cyp2b13, Cyp2d9 in pericentral
- [ ] Cyp2b9 in periportal
- [ ] Mtmr11 identified as zonation-specific regulator

### Structural Analysis
- [ ] Gene-lipid correlation matrix generated
- [ ] Lipid structural features (C, DB) extracted
- [ ] Regression models identify genes with structural specificity
- [ ] Visualization plots generated for key genes

---

## Contact & Support

If you encounter issues not covered in this guide:

1. Check `ANALYSIS_SCOPE_SUMMARY.md` for expected analysis details
2. Review error messages carefully
3. Ensure all required input files exist
4. Verify R package versions match those used in original analysis

---

**Last Updated:** 2025-01-28
