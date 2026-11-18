# Workflow 06: Structural Analysis

**Script:** `scripts/workflows/06_structural_analysis.R`  
**Lines:** ~534  
**Purpose:** Current active workflow performs lipid-only structural analysis. It parses lipid names from workflow 03 limma results, tests whether chain length and double-bond count explain gonadal logFC within lipid classes, and writes the structural feature table consumed by workflow 07.

**Current-status note:** The script header and `scripts/core/structural_analysis.R` still document the older gene-lipid structural-regression concept. In the active workflow, those older gene-lipid correlation sections are not run; the script ends with a note that sections 6-14 are commented out.

## Current execution contract (audited 18 November 2026)

This direct script executes under both `source()` and `Rscript`; there is no CLI
parser, source guard, validate-only or fresh-output override. It overwrites named
outputs under `results/lipidomic/structural_model/`, so reproduce in an isolated
checkout/copy when retaining curated outputs. The current public WF03 result
RData is sufficient with the R environment; no RNA/Seurat object, MOFA result,
zonation result or bubble-annotation CSV is loaded by this active workflow.

Run after03 when recomputing its input. Both maintained07 modes read this
workflow's `lipids_with_structure.csv` for structural reference axes: default zonal
validation and explicit `--global-only`. Config/dataset files are sourced, but
input/output paths and the four plotted focus classes are hardcoded locally.
Conditional significance CSVs are not removed when a later run has no qualifying
rows, so an old file in a reused output root is not proof of a current result.

## Source Dependencies

```r
source("scripts/config/analysis_config.R")
source("scripts/config/datasets.R")
source("scripts/support_functions.R")
source("scripts/core/data_loading.R")
source("scripts/core/structural_analysis.R")
```

## Section Breakdown

| # | Section | Lines (approx) | Description |
|---|---------|----------------|-------------|
| 1 | Configuration | 48-68 | Load config, dataset definitions, support/core helpers, and create `results/lipidomic/structural_model/`. |
| 2 | Load Lipid Structural Features | 71-126 | Load `results/lipidomic/ttop_full_list.RData`, extract `Testes vs Ovaries`, parse lipid names into class, chain length, DB count, and logFC. |
| 3 | Lipid-Only Structural Analysis | 130-198 | Run class-specific regressions with `run_lipid_only_class_regression()` and extract p < 0.05 effects. |
| 3b | Lipid-Only Plots | 201-442 | Plot focus classes with shared axes and fixed color scale. |
| 3c | Summary Heatmap | 445-503 | Plot signed `-log10(p)` heatmap for chain length, DB, and interaction effects across lipid classes. |
| - | Complete / Historical Note | 510-534 | Print outputs and record that older gene-lipid structural sections are commented out. |

## Input Files

- `results/lipidomic/ttop_full_list.RData` - Lipid differential abundance results from workflow 03. The active script expects nested `ttop_full_list[["Testes vs Ovaries"]][["Testes vs Ovaries"]]`.

## Output Files

- `results/lipidomic/structural_model/lipids_with_structure.csv` - Parsed lipid structural table; consumed by workflow 07 section 11.
- `results/lipidomic/structural_model/lipid_only_class_coefficients.csv`
- `results/lipidomic/structural_model/lipid_only_class_summaries.csv`
- `results/lipidomic/structural_model/lipid_only_chain_length_effects.csv` - Written only if significant chain-length effects exist.
- `results/lipidomic/structural_model/lipid_only_db_position_effects.csv` - Written only if significant DB effects exist.
- `results/lipidomic/structural_model/lipid_only_interaction_effects.csv` - Written only if significant `chain_length:db_position` effects exist.
- `results/lipidomic/structural_model/lipid_only_summary_heatmap.pdf`
- `results/lipidomic/structural_model/lipid_only_plots/[CLASS]_structural_effects.pdf`

## Key Variables

| Variable | Type | Description |
|----------|------|-------------|
| `structural_dir` | character | `results/lipidomic/structural_model` |
| `ttop_gonadal` | data.frame | Workflow 03 `Testes vs Ovaries` lipid limma table |
| `lipid_structure` | data.frame | `lipid_name`, `lipid_class`, `chain_length`, `db_position`, `logfc` |
| `lipid_only_results` | list | Regression coefficients and class summaries from `run_lipid_only_class_regression()` |
| `lipid_only_sig` | list | Significant chain, DB, and interaction effect tables |
| `focus_lipid_classes` | character | `c("SM", "BMP", "PG", "TG")` |

## Coherence Assessment

**Status: GOOD, with intentional historical comments.** The active script is linear and self-contained after workflow 03. The main documentation hazard is conceptual: older gene-lipid structural analysis is described in comments/modules but is not currently executed.

## Vulnerabilities to Keep in Mind When Modifying

1. **Nested lipid result assumption:** Section 2 assumes the workflow 03 RData has `ttop_full_list[["Testes vs Ovaries"]][["Testes vs Ovaries"]]`. If workflow 03's list structure is flattened, section 2 must be made tolerant like workflow 05 section 14.
2. **Conditional output files:** The three `lipid_only_*effects.csv` files are written only when significant effects exist. Downstream docs/agents should not assume all three exist after every run.
3. **WF06 -> WF07 dependency:** Workflow 07 section 11 expects `results/lipidomic/structural_model/lipids_with_structure.csv` and focus classes `SM`, `BMP`, `PG`, `TG`. Keep these synchronized.
4. **Conceptual drift:** If the older gene-lipid structural regression is re-enabled, update this plan, README, AGENT_GUIDE, and the meta dependency map because WF06 would again depend on RNA/MOFA files.

## Active Regression Model

```r
# Per lipid class with at least 10 lipids:
logfc ~ chain_length + db_position + chain_length:db_position
```

## Notes for Context Restoration

1. Current outputs live under `results/lipidomic/structural_model/`, not `results/multi-omic/structural_model/`.
2. The workflow tests structural effects on lipid gonadal logFC, not gene-lipid correlations.
3. Plot focus classes are hardcoded to `SM`, `BMP`, `PG`, and `TG`.
4. The structural parser comes from `scripts/core/structural_analysis.R` and extracts total carbons and total double bonds from lipid names.
