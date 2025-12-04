# Workflow 05 Update: Lipid Class Association Visualization

## Summary

Added a new section 11 to `workflows/05_zonation_analysis.R` that creates a visualization showing which lipid classes are associated with genes from each of the four zonation/gonad categories.

## What was added

### New Section 11: Lipid Class Association by Zonation/Gonad Category

This section:

1. **Loads lipidomic differential analysis** results to determine gonadal bias (Testes vs Ovaries) for each lipid class
2. **Extracts lipid class associations** from the four correlation networks created in section 10:
   - **TPC**: Testes Pericentral genes
   - **TPP**: Testes Periportal genes
   - **OPC**: Ovaries Pericentral genes
   - **OPP**: Ovaries Periportal genes
3. **Counts associations**: For each category, counts how many genes correlate with each lipid class
4. **Determines gonadal bias**: Assigns each lipid class as "Testes-biased" (blue) or "Ovaries-biased" (red) based on average logFC across all lipids in that class
5. **Creates visualization**: Generates a faceted bar plot showing:
   - **X-axis**: Lipid classes (top 20 by total associations)
   - **Y-axis**: Number of associated genes
   - **Facets**: Four panels for TPC, TPP, OPC, OPP
   - **Colors**: Blue for testes-biased lipid classes, red for ovaries-biased

## Output Files

### 1. `results/zonation/lipid_class_by_category.pdf`
A publication-quality plot showing lipid class associations across all four categories, matching the sketch provided.

### 2. `results/zonation/lipid_class_category_summary.csv`
A summary table containing:
- `lipid_class`: Lipid class name (PC, PE, TG, etc.)
- `n_genes`: Number of genes associated with this lipid class
- `category`: Gene category (TPC, TPP, OPC, OPP)
- `mean_logFC`: Average log fold change for this lipid class (Testes vs Ovaries)
- `gonad_bias`: "Testes" or "Ovaries" based on mean_logFC

## Results Summary

From the latest run:
- **39 lipid class × category combinations**
- **18 unique lipid classes** associated with genes
- **Distribution by category:**
  - TPC: 11 lipid classes (7 testes-biased, 4 ovaries-biased)
  - TPP: 12 lipid classes (7 testes-biased, 2 ovaries-biased)
  - OPC: 7 lipid classes (4 testes-biased, 1 ovaries-biased)
  - OPP: 9 lipid classes (5 testes-biased, 2 ovaries-biased)

## Key Biological Insights

This visualization allows you to identify:

1. **Shared lipid classes**: Which lipid classes are associated with genes across multiple categories
2. **Category-specific lipid classes**: Which lipid classes are unique to specific zonation/gonad combinations
3. **Gonadal bias patterns**: Whether certain categories preferentially associate with testes-biased or ovaries-biased lipid classes

For example:
- If PC (phosphatidylcholine) appears in all four categories, it suggests this lipid class is broadly regulated by zonation-specific genes regardless of gonad type
- If certain lipid classes only appear in TPC or TPP, they may be specifically regulated by testes-specific zonation patterns

## Implementation Details

### Key Functions
- `summarize_lipid_classes()`: Counts gene associations per lipid class for each network
- Filters out non-lipid entries (gene names) by matching against known lipid classes from the lipidomics data

### Data Flow
1. Load lipidomic differential analysis → extract lipid class and gonadal bias
2. Load compressed correlation networks (already aggregated by lipid class in section 10)
3. For each network, count genes per lipid class
4. Combine all four summaries and filter to actual lipid classes
5. Merge with gonadal bias information
6. Create faceted bar plot

## Testing

The workflow was tested end-to-end and successfully:
- Loaded 1,254 lipids across 31 lipid classes
- Processed 4 correlation networks
- Generated both output files
- Filtered out 25 non-lipid entries (gene names) to keep only 39 valid lipid class associations

## Next Steps

You can now use this visualization to:
1. Identify key lipid classes for follow-up experimental validation
2. Compare with lipid structural analysis results (workflow scripts for structural analysis)
3. Include in your manuscript as a figure panel showing lipid-gene coordination across zonation and gonad axes
