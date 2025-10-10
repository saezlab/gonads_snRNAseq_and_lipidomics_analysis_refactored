# Refactoring Plan: FCG Liver Multi-Omic Analysis

**Goal:** Reduce script redundancy, establish clear structure, and document complex analysis steps

**Date:** 2025-10-10

---

## 1. Current State Analysis

### Script Inventory (48 scripts total)

**By Function:**
- **Data Conversion:** 7 scripts (seurat_to_df variants)
- **Differential Analysis:** 10 scripts (6 pseudobulk, 4 lipidomic)
- **Pathway/TF/Enrichment:** 5 scripts
- **Multi-omic Integration:** 3 scripts (MOFA)
- **Structural Analysis:** 3 scripts (lipid structure)
- **Validation:** 8 scripts
- **Visualization:** 7 scripts (figures, bubbles, volcanos)
- **Comparison:** 4 scripts
- **Support:** 1 script (support_functions.R)

### Identified Redundancy Patterns

#### Pattern 1: Dataset Variants
**Problem:** Same analysis repeated for different datasets with minimal changes

- **Pseudobulk DA:** 6 variants
  - `pseudobulk_DA.R` (BASE - FCG main)
  - `pseudobulk_DA_gdhxr.R` (gonadectomy + hormone replacement)
  - `pseudobulk_DA_JPC.R` (specific contrast version)
  - `pseudobulk_DA_WT.R` (wild-type comparison)
  - `pseudobulk_DA_zonation.R` (spatial zonation)
  - `gdxhr_pseudobulk_DA_zonation.R` (zonation + gdxhr)

- **Lipidomic DA:** 4 variants
  - `lipidomic_DA.R` (BASE - FCG main)
  - `lipidomic_DA_JPC.R`
  - `gdxhr_lipidomic_DA.R`
  - `lipidomic_validation_DA.R`

- **Seurat to DataFrame:** 5 variants
  - `seurat_to_df.R` (BASE)
  - `seurat_to_df_withWT.R`
  - `seurat_to_df_zonation.R`
  - `seurat_to_df_pnpla3ko_zonation.R`
  - `gdxhr_seurat_to_df.R` + `gdxhr_seurat_to_df_zonation.R`

**Impact:** ~40% of scripts are variants (19/48)

#### Pattern 2: Iterative Versions
**Problem:** Multiple numbered versions suggesting development iterations

- Validation scripts:
  - `lipidomic_gonad_validation.R` / `lipidomic_gonad_validation_2.R`
  - `lipidomic_gonad_sszb_validation.R` / `lipidomic_gonad_sszb_validation_2.R`
  - Each with additional `_JPC.R` variant

- MOFA downstream:
  - `MOFA_lipidcell_downstream.R`
  - `MOFA_lipidcell_downstream_largenet.R`

- Bubble plots:
  - `lipidomic_bubbles.R`
  - `lipidomic_bubbles_advanced.R`

**Impact:** 8 scripts are likely redundant iterations

#### Pattern 3: Hardcoded Parameters
**Problem:** Scripts have dataset-specific hardcoded values instead of parameters

Common hardcoded elements:
- File paths
- Contrast definitions
- Sample filtering criteria
- Color schemes
- Plot dimensions

**Impact:** Difficult to reuse analysis logic across datasets

---

## 2. Proposed Refactored Structure

### Module Organization

```
scripts/
├── core/                           # Core analysis modules
│   ├── data_loading.R              # Unified data loading functions
│   ├── differential_analysis.R     # Generic DA functions for RNA/lipid
│   ├── enrichment_analysis.R       # Pathway/TF/LR enrichment
│   ├── multiomics.R                # MOFA integration functions
│   └── structural_analysis.R       # Lipid structural analysis
│
├── visualization/                  # Visualization modules
│   ├── volcano_plots.R             # Volcano plot functions
│   ├── bubble_plots.R              # Bubble plot functions
│   ├── dimred_plots.R              # PCA/UMAP plotting
│   └── heatmaps.R                  # Heatmap utilities
│
├── workflows/                      # High-level analysis workflows
│   ├── 01_data_preparation.R       # Load and prepare all datasets
│   ├── 02_fcg_pseudobulk_analysis.R
│   ├── 03_fcg_lipidomic_analysis.R
│   ├── 04_gdxhr_validation.R       # GDX+HR validation
│   ├── 05_zonation_analysis.R
│   ├── 06_mofa_integration.R
│   ├── 07_structural_analysis.R
│   ├── 08_validation_experiments.R
│   └── 09_figures_generation.R
│
├── config/                         # Configuration files
│   ├── datasets.R                  # Dataset definitions
│   ├── contrasts.R                 # Contrast definitions
│   ├── plotting_themes.R           # Color schemes, themes
│   └── paths.R                     # File path configurations
│
└── support_functions.R             # Keep existing support functions
```

### Key Improvements

1. **Separation of Concerns:**
   - **Core modules:** Reusable analysis logic
   - **Workflows:** Dataset-specific pipelines that call core modules
   - **Config:** Centralized parameters and settings
   - **Visualization:** Plotting separated from analysis

2. **Parameterization:**
   - Dataset configs define: file paths, samples, metadata
   - Contrast configs define: comparison groups
   - Core functions accept parameters instead of hardcoding

3. **Clarity:**
   - Numbered workflows show analysis order
   - Each workflow script corresponds to a section in the paper
   - Clear naming conventions

---

## 3. Detailed Refactoring Steps

### Phase 1: Core Module Creation (Foundation)

#### Step 1.1: Create `core/data_loading.R`
**Extract from:** `seurat_to_df*.R`, `pseudobulk_DA*.R`, `lipidomic_DA*.R`

**Functions to create:**
```r
load_rna_data(dataset_name, split_files = TRUE)
load_lipidomic_data(dataset_name)
load_metadata(dataset_name)
filter_samples(data, criteria)
prepare_pseudobulk(seurat_obj, group_by, cell_type)
```

**Consolidates:** 5 seurat_to_df variants → 1 function with parameters

#### Step 1.2: Create `core/differential_analysis.R`
**Extract from:** `pseudobulk_DA*.R`, `lipidomic_DA*.R`

**Functions to create:**
```r
run_differential_analysis(data, metadata, contrasts, method = "limma")
  # Unified DA function for both RNA and lipid data
  # Parameters:
  #   - data: expression/abundance matrix
  #   - metadata: sample information
  #   - contrasts: list of contrasts to test
  #   - method: "limma" (default), other methods possible

define_contrasts(design, comparison_type)
  # Helpers to define contrast vectors
  # comparison_type: "gonadal", "chromosomal", "group", "custom"

filter_significant_features(results, fdr_cutoff = 0.05, fc_cutoff = 0)
  # Consistent significance filtering

save_da_results(results, output_dir, prefix)
  # Standardized output format
```

**Consolidates:** 10 DA scripts → parameterized functions

#### Step 1.3: Create `core/enrichment_analysis.R`
**Extract from:** `pseudobulk_pathways*.R`, `pseudobulk_TFs.R`, `pseudobulk_LRs.R`

**Functions to create:**
```r
run_pathway_enrichment(gene_stats, pathway_db, method = "ulm")
run_tf_activity(gene_stats, regulons)
run_lr_analysis(expression_data, lr_database)
run_lipid_class_enrichment(lipid_stats)
```

**Consolidates:** 5 enrichment scripts → unified enrichment module

#### Step 1.4: Create `core/multiomics.R`
**Extract from:** `MOFA_lipidcell_RNA.R`, `MOFA_lipidcell_downstream*.R`

**Functions to create:**
```r
prepare_mofa_input(rna_data, lipid_data, sample_mapping)
train_mofa_models(data, factor_range = 5:15)
select_optimal_model(model_list, criteria = "r2")
extract_mofa_features(model, top_percent = 0.05, factor_n = 1)
compute_feature_correlations(features, cor_threshold = 0.65)
build_correlation_network(mofa_weights, correlations, filters)
```

**Consolidates:** 3 MOFA scripts → unified pipeline

#### Step 1.5: Create `core/structural_analysis.R`
**Extract from:** `lipid_rna_correaltion_for_structural_analysis.R`, `lipid_structure_analysis.R`

**Functions to create:**
```r
compute_gene_lipid_correlations(rna_data, lipid_data, gene_filter)
load_lipid_structural_features(lipid_list)
run_structural_regression(correlations, structural_features)
identify_structural_associations(regression_results, p_threshold = 0.05)
```

**Consolidates:** 2 of 3 structural scripts (keep visualization separate)

### Phase 2: Visualization Module

#### Step 2.1: Create `visualization/volcano_plots.R`
**Extract from:** `pseudobulk_DA.R`, `lipidomic_DA.R`, `misc_volcanos.R`

**Functions to create:**
```r
volcano_plot(da_results, title, colors, label_top_n = 10)
  # Unified volcano plot function
  # Uses existing volcano_nice() from support_functions.R as base

save_volcano_plot(plot, output_path, width = 8, height = 6)
```

#### Step 2.2: Create `visualization/bubble_plots.R`
**Extract from:** `lipidomic_bubbles*.R`

**Functions to create:**
```r
lipid_bubble_plot(lipid_data, x = "C", y = "DB", color = "LogFC")
plot_all_lipid_classes(data_dir, output_dir, plot_function)
```

**Consolidates:** 2 bubble scripts → parameterized functions

#### Step 2.3: Create `visualization/dimred_plots.R`
**Extract from:** `figure_1_UMAP.R`, `pseudobulk_DA.R` (PCA sections)

**Functions to create:**
```r
plot_pca(data, metadata, color_by, shape_by = NULL)
plot_umap(seurat_obj, group_by, colors = NULL)
plot_variance_explained(pca_obj, n_pcs = 20)
```

#### Step 2.4: Create `visualization/heatmaps.R`
**Extract from:** Various scripts with pheatmap/ComplexHeatmap calls

**Functions to create:**
```r
expression_heatmap(data, annotation, scale = "row", colors = NULL)
pathway_heatmap(enrichment_scores, cluster = TRUE)
```

### Phase 3: Configuration Files

#### Step 3.1: Create `config/datasets.R`
**Define all dataset configurations:**

```r
DATASETS <- list(
  fcg_main = list(
    rna = list(
      files = c("data/RNA/scaled_RNA_1.csv",
                "data/RNA/scaled_RNA_2.csv",
                "data/RNA/scaled_RNA_3.csv"),
      metadata = "data/RNA/meta_data_df.csv",
      split = TRUE
    ),
    lipid = list(
      file = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
      metadata_col = "sample_id"
    ),
    description = "Four Core Genotypes main dataset"
  ),

  gdxhr = list(
    rna = list(
      files = c("data/RNA/scaled_RNA_gdxhr_1.csv",
                "data/RNA/scaled_RNA_gdxhr_2.csv",
                "data/RNA/scaled_RNA_gdxhr_3.csv"),
      metadata = "data/RNA/meta_data_df_gdxhr.csv",
      split = TRUE
    ),
    lipid = list(
      file = "data/metabolomics/lipidomics_gdxhr.csv"
    ),
    description = "Gonadectomy + hormone replacement validation"
  ),

  # ... other datasets ...
)
```

#### Step 3.2: Create `config/contrasts.R`
**Define all contrast configurations:**

```r
CONTRASTS <- list(
  gonadal_type = list(
    name = "Testes vs Ovaries",
    vector = c(-1, 2),
    description = "Primary gonadal hormone effect",
    colors = c(Ovaries = "red", Testes = "royalblue")
  ),

  allosome_type = list(
    name = "XY vs XX",
    vector = c(-1, 2),
    description = "Sex chromosome effect",
    colors = c(XX = "red", XY = "royalblue")
  ),

  # ... other contrasts ...
)
```

#### Step 3.3: Create `config/plotting_themes.R`
**Centralize color schemes and plot themes:**

```r
# Cell type colors
CELLTYPE_COLORS <- c(
  "Hepatocytes" = "#ca3659",
  "Endothelial cells" = "#8f7b2f",
  "Immune cells" = "#0d934d",
  "Stellate cells" = "#1892b1",
  "Cholangiocytes" = "#834f9f"
)

# Standard plot dimensions
PLOT_DIMS <- list(
  volcano = c(width = 8, height = 6),
  heatmap = c(width = 10, height = 8),
  pca = c(width = 7, height = 6)
)
```

### Phase 4: Workflow Scripts

#### Step 4.1: Create numbered workflow scripts

Each workflow script:
1. Sources necessary core modules
2. Loads dataset configs
3. Runs analysis pipeline
4. Saves results to standardized locations
5. Includes inline documentation explaining complex steps

**Example: `workflows/02_fcg_pseudobulk_analysis.R`**
```r
# FCG Pseudobulk Differential Expression Analysis
#
# This workflow performs differential expression analysis on the Four Core
# Genotypes dataset to identify genes regulated by gonadal hormones vs.
# sex chromosomes.

source("scripts/core/data_loading.R")
source("scripts/core/differential_analysis.R")
source("scripts/visualization/volcano_plots.R")
source("scripts/config/datasets.R")
source("scripts/config/contrasts.R")

# Load data
cat("Loading FCG RNA data...\n")
rna_data <- load_rna_data("fcg_main")
metadata <- load_metadata("fcg_main", type = "rna")

# Define contrasts
contrasts_to_test <- list(
  CONTRASTS$gonadal_type,
  CONTRASTS$allosome_type
)

# Run analysis for each cell type
cell_types <- unique(metadata$cell_type)
results <- list()

for (cell_type in cell_types) {
  cat("Analyzing", cell_type, "...\n")

  # Prepare pseudobulk
  pseudobulk <- prepare_pseudobulk(rna_data, metadata, cell_type)

  # Run differential analysis
  da_results <- run_differential_analysis(
    data = pseudobulk,
    metadata = metadata,
    contrasts = contrasts_to_test
  )

  results[[cell_type]] <- da_results
}

# Save results
save_da_results(results, "results/pseudobulk/", prefix = "fcg_main")

# Generate visualizations
for (contrast in names(results$Hepatocytes)) {
  plot <- volcano_plot(
    results$Hepatocytes[[contrast]],
    title = contrast,
    colors = CONTRASTS[[contrast]]$colors
  )
  save_volcano_plot(plot, paste0("results/pseudobulk/volcanos/", contrast, ".pdf"))
}
```

---

## 4. Documentation Priorities

### Complex Analysis Steps Requiring Documentation

#### Priority 1: MOFA Integration (HIGH complexity)
**Location:** `core/multiomics.R` + dedicated doc

**Document:**
1. **Factor selection rationale**
   - Why test 5-15 factors?
   - How to interpret model comparison
   - What does R² per factor mean?

2. **Weight matrix interpretation**
   - What are MOFA weights?
   - Why filter top 50% by weight?
   - How to identify factor-associated features

3. **Network construction**
   - Correlation threshold selection (why 0.65?)
   - Filtering strategy for lipid metabolism genes
   - TF score calculation

**Documentation file:** `docs/MOFA_INTEGRATION_EXPLAINED.md`

#### Priority 2: Structural Analysis (HIGH complexity)
**Location:** `core/structural_analysis.R` + dedicated doc

**Document:**
1. **Regression model interpretation**
   - What does gene-lipid correlation mean?
   - How structural features predict correlations
   - Interaction term interpretation (abundance × structure)

2. **Biological significance**
   - Chain length associations → substrate specificity
   - Double bond position → saturation preferences
   - Why this reveals enzyme function

**Documentation file:** `docs/STRUCTURAL_ANALYSIS_EXPLAINED.md`

#### Priority 3: Contrast Definitions (MEDIUM complexity)
**Location:** `config/contrasts.R` + inline documentation

**Document:**
1. **Contrast vector mechanics**
   - How numeric vectors define contrasts in limma
   - Order of conditions in design matrix
   - Complex contrasts (averaging groups)

2. **Biological interpretation**
   - Gonadal vs chromosomal effects
   - FCG model design and interpretation

**Documentation file:** `docs/CONTRAST_DEFINITIONS_EXPLAINED.md`

#### Priority 4: decoupleR Usage (MEDIUM complexity)
**Location:** `core/enrichment_analysis.R` + inline comments

**Document:**
1. **run_ulm() method**
   - What is Univariate Linear Model?
   - Input format requirements
   - Output interpretation (scores, p-values)

2. **Network format**
   - Source-target-mor structure
   - Mode of regulation (+1 / -1)
   - Why this differs from traditional ORA/GSEA

**Inline documentation** in functions

#### Priority 5: Zonation Analysis (LOW-MEDIUM complexity)
**Location:** Workflow scripts + inline comments

**Document:**
1. **Marker-based zonation**
   - Cyp2e1 (pericentral) vs Cyp2f2 (periportal)
   - How to define zonation-biased genes

2. **Integration with gonad effects**
   - Why zonation modulates hormone responses

**Inline documentation** in workflow scripts

---

## 5. Implementation Order

### Week 1: Foundation
- [ ] Create module directory structure
- [ ] Implement `core/data_loading.R`
- [ ] Implement `config/datasets.R`
- [ ] Test data loading with existing scripts

### Week 2: Core Analysis
- [ ] Implement `core/differential_analysis.R`
- [ ] Implement `config/contrasts.R`
- [ ] Create `workflows/02_fcg_pseudobulk_analysis.R`
- [ ] Validate against original `pseudobulk_DA.R` results

### Week 3: Enrichment & Visualization
- [ ] Implement `core/enrichment_analysis.R`
- [ ] Implement `visualization/` modules
- [ ] Create `workflows/03_fcg_lipidomic_analysis.R`
- [ ] Validate lipidomic workflow

### Week 4: Multi-omics & Advanced
- [ ] Implement `core/multiomics.R`
- [ ] Implement `core/structural_analysis.R`
- [ ] Create `workflows/06_mofa_integration.R`
- [ ] Create `workflows/07_structural_analysis.R`

### Week 5: Documentation & Validation Workflows
- [ ] Create `docs/MOFA_INTEGRATION_EXPLAINED.md`
- [ ] Create `docs/STRUCTURAL_ANALYSIS_EXPLAINED.md`
- [ ] Create `docs/CONTRAST_DEFINITIONS_EXPLAINED.md`
- [ ] Implement remaining workflow scripts
- [ ] Full pipeline validation

### Week 6: Testing & Finalization
- [ ] Update test framework for new structure
- [ ] Add tests for core modules
- [ ] Final validation: reproduce all main results
- [ ] Update CLAUDE.md with new structure
- [ ] Archive old scripts (don't delete)

---

## 6. Validation Strategy

**Critical:** Must reproduce existing results exactly

### Validation Checklist
- [ ] FCG pseudobulk: Same number of DEGs at FDR < 0.05
- [ ] FCG lipidomics: Same number of significant lipids
- [ ] MOFA model: Same factor scores and weights
- [ ] Structural analysis: Same regression coefficients
- [ ] All figures: Visual comparison to originals

### Validation Approach
1. Run old script → save results as reference
2. Run new workflow → save results as test
3. Compare numerically (ttop tables, scores, etc.)
4. Document any differences and explain

---

## 7. Migration Plan for Old Scripts

**DO NOT DELETE old scripts immediately**

### Archive Strategy
```
scripts/
├── archive/                    # Old scripts moved here
│   ├── README.md               # Explains mapping to new structure
│   ├── original_scripts/       # Original unmodified scripts
│   └── validation_results/     # Outputs from old scripts for comparison
```

### Mapping Document
Create `scripts/archive/MIGRATION_MAP.md`:

```markdown
# Script Migration Map

## Differential Analysis
- `pseudobulk_DA.R` → `workflows/02_fcg_pseudobulk_analysis.R`
- `pseudobulk_DA_gdhxr.R` → `workflows/04_gdxhr_validation.R`
- `pseudobulk_DA_zonation.R` → `workflows/05_zonation_analysis.R`
...
```

---

## 8. Expected Benefits

### Quantitative Improvements
- **Script count reduction:** 48 → ~20 (58% reduction)
  - 10 core modules
  - 9 workflow scripts
  - 1 support functions (retained)

- **Code reuse:** Eliminate ~12,000 lines of redundant code

- **Maintenance:** Single place to fix bugs/update methods

### Qualitative Improvements
- **Clarity:** Numbered workflows show analysis order
- **Reproducibility:** Parameterized configs make it easy to rerun with different settings
- **Documentation:** Complex steps explained in dedicated docs
- **Extensibility:** Easy to add new datasets or contrasts
- **Testing:** Modular functions easier to test

---

## 9. Risks and Mitigation

### Risk 1: Results Don't Match Original
**Mitigation:**
- Implement validation checklist at each phase
- Keep old scripts until full validation complete
- Document any necessary differences

### Risk 2: Breaking Changes During Refactor
**Mitigation:**
- Use git branches for refactoring work
- Frequent commits with clear messages
- Maintain both old and new structures temporarily

### Risk 3: Over-abstraction
**Mitigation:**
- Prioritize readability over DRY principle
- Keep workflow scripts explicit rather than overly generic
- Inline documentation for complex parameters

---

## 10. Success Criteria

Refactoring is complete when:

1. ✅ All main analyses reproducible via numbered workflows
2. ✅ Core modules tested and validated
3. ✅ Complex analyses documented in detail
4. ✅ Test framework updated and passing
5. ✅ CLAUDE.md updated with new structure
6. ✅ Old scripts archived with migration map
7. ✅ All main results reproduced and validated
8. ✅ Code review by collaborators passed

---

## Notes

- This is an ambitious refactoring plan
- Can be implemented incrementally
- Priority should be given to most redundant scripts first
- Validation is critical - scientific reproducibility over code elegance
