# Analysis Scope Summary - FCG Liver Study

**Based on manuscript:** "Gonadal hormone coordination of sex-specific adult liver physiology"

---

## 1. EXPERIMENTAL DATASETS

### Primary Datasets (FCG Model)
- **sci-RNA-seq:** 466,986 high-quality nuclei from 48 livers
  - 43 FCG livers: ~10-11 per group (XXO, XXT, XYO, XYT)
  - 5 wild-type male (XYT-WT) controls
  - 5 major cell types: Hepatocytes, Immune cells, Endothelial cells, Stellate cells, Cholangiocytes

- **Lipidomics:** Untargeted lipidomics on same 48 livers
  - 1,386 lipid species across 42 lipid classes
  - Analytical precision ~10%

### Validation Datasets

**Gonadectomy (GDX) Experiments:**
- Wild-type C57BL/6J male and female mice
- 6 weeks old at GDX, analyzed 5 weeks later
- sci-RNA-seq + lipidomics on same samples
- Used to validate activational hormone effects

**Xenium Spatial Transcriptomics:**
- Wild-type male and female livers (n=5 per sex)
- Custom 100-gene panel:
  - 16 liver cell type markers
  - 81 gonad-dependent genes (from sci-RNA-seq)
  - 1 hormone receptor (Ar)
  - 2 X chromosome genes (Xist, Kdm6a)
- In situ validation of gene expression patterns

**Bulk Proteomics:**
- Wild-type male and female liver tissue
- TMT16plex labeling
- Validates gene expression at protein level

**RNA FISH:**
- Wild-type male and female liver cryosections (n=5 per sex)
- Validates specific genes: Cyp2d9 (male-biased), Cyp2b9 (female-biased)
- HCR v3.0 approach

**siRNA Knockdown Experiments:**
- Primary hepatocytes from wild-type female mice
- Targets validated:
  - **MOFA-identified:** Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9
  - **Zonation-specific:** Sult3a1, Mtmr11, Abcb1a, Mindy3, Esco2, Kitl, Cyp2a4
- 48 hours post-transfection
- Bulk RNA-seq + lipidomics readout

---

## 2. KEY CONTRASTS TESTED

### Main Contrasts (FCG Model)

**A. Gonadal Type (PRIMARY):**
- **Comparison:** Testes vs. Ovaries (pooling XXT+XYT vs. XXO+XYO)
- **Expected result:** ~4,000 DEGs in hepatocytes (FDR ≤ 0.05)
- **Biological interpretation:** Effect of gonadal hormones
- **Lipids:** 733/1,386 lipids (53.9%) gonad-dependent

**B. Sex Chromosome Complement (allosome):**
- **Comparison:** XY vs. XX (pooling XYT+XYO vs. XXT+XXO)
- **Expected result:** Few DEGs (mainly Y-linked, X-escapees)
- **Biological interpretation:** Effect of sex chromosomes (minimal)
- **Lipids:** 0/1,386 lipids significantly different

**C. FCG Group-Specific Comparisons:**
These isolate hormone or chromosome effects with the other held constant:

1. **XXT vs. XXO** - Gonadal effect in XX background
2. **XYT vs. XYO** - Gonadal effect in XY background
3. **XYO vs. XXO** - Chromosome effect in ovaries background
4. **XYT vs. XXT** - Chromosome effect in testes background

**D. XYT vs. XYT-WT:**
- Controls for FCG model artifacts
- Expected: ~7% of lipids differ, ~93% concordant

### Validation Contrasts

**E. Gonadectomy (GDX):**
- **Male:** GDX vs. Sham (activational testosterone effect)
- **Female:** GDX vs. Sham (activational estrogen effect)
- **Expected correlation with FCG:** r = -0.48 (inverse, as expected)

**F. Zonation-Specific:**
- **Pericentral vs. Periportal** hepatocytes
- **Gonadal effect within zonation:**
  - Testes vs. Ovaries in pericentral hepatocytes
  - Testes vs. Ovaries in periportal hepatocytes
- **Expected:** 43% of gonad-dependent genes are pericentral-biased, 18% periportal-biased

**G. siRNA Knockdown:**
- **KD vs. Control** for each target gene
- Validates causal gene-lipid relationships

---

## 3. ANALYSES PERFORMED

### 3.1 RNA-seq Analyses

**Preprocessing & QC:**
- Seurat v5 pipeline
- Quality filters: >500 UMIs, <3 SD above mean
- Doublet removal: scDblFinder (2% doublet rate)
- Cell type annotation using canonical markers
- Zonation assignment: Cyp2e1+ (pericentral), Cyp2f2+ (periportal)

**Pseudobulk Aggregation:**
- Per cell type × biological replicate (mouse)
- ~215 pseudobulk samples (43 FCG × 5 cell types)
- Log-normalized, scaled for differential analysis

**Differential Expression:**
- **Method:** limma (moderated t-tests)
- **Contrasts:** All contrasts listed in Section 2
- **Cell types analyzed:** All 5 major types (focus on hepatocytes)
- **Thresholds:** FDR ≤ 0.05 for significance
- **Output:** Volcano plots, heatmaps, gene lists

**Pathway Enrichment:**
- **Method:** decoupleR + ULM (Univariate Linear Model)
- **Gene sets:**
  - MSigDB GO biological processes (m5.all.v2023.2.Mm)
  - CollecTRI TF regulons
  - Custom ESR1/AR targets (Zheng et al. 2018)
- **Key pathways identified:**
  - Testes-biased: Ribosomal complexes, translation
  - Ovaries-biased: Steroid/sex hormone metabolism

**TF Activity Inference:**
- **Method:** decoupleR + ULM
- **TF regulons:** CollecTRI + custom (ESR1, AR, STAT5a, STAT5b)
- **Key findings:**
  - ESR1 activity: -5.5 (p = 4.17 × 10⁻⁸)
  - AR activity: -2.9 (p = 3.79 × 10⁻³)
  - STAT5a activity: -3.12 (p = 1.7 × 10⁻³)

**PCA:**
- Per cell type, colored by experimental factors
- Ridgeline plots showing PC1 distribution by group
- Assesses variance explained by gonads vs. chromosomes

### 3.2 Lipidomic Analyses

**Preprocessing:**
- MS-DIAL 4.9 processing
- Total ion count normalization
- Feature identification: accurate mass, MS/MS, retention time
- Deduplication of multi-mode detections

**Statistical Testing:**
- **Method:** Wilcoxon Rank-Sum tests (non-parametric, many lipids non-normal)
- **Multiple testing correction:** Benjamini-Hochberg FDR
- **Contrasts:** Same as RNA (gonadal, allosome, group-specific)
- **Threshold:** FDR ≤ 0.05

**Lipid Class Enrichment:**
- **Method:** decoupleR + ULM
- Creates lipid class → lipid species network
- Tests which classes are enriched in testes vs. ovaries
- **Key findings:**
  - Male-biased: TG, DG, glycerophospholipids (PI, PE, LPE, PA, LPC, PG, PC)
  - Female-biased: Sphingolipids (SM, Cer, HexCer)

**Network Analysis:**
- **Tool:** LINEX2 (Rose et al. 2023)
- Lipid-lipid conversion reactions
- ~300 lipids, 194 edges
- Top 5 testes + top 5 ovaries lipids form interconnected subnetwork

**Structural Analysis:**
- **Features:** Chain length (C), double bond position (DB)
- **Method:** Linear regression modeling gene-lipid correlations as function of lipid structure
- **Models:**
  - Main: `correlation ~ class + abundance + C + DB + abundance:C + abundance:DB`
  - Class-specific: Same model within each lipid class (≥10 lipids/class)
- **Output:** Identifies genes with substrate specificity for structural features

**PCA:**
- All 1,386 lipids, all 48 samples
- Samples cluster by gonad identity (not chromosome)

### 3.3 Multi-Omic Integration

**MOFA (Multi-Omics Factor Analysis):**
- **Input views:**
  - Hepatocyte RNA expression (top variable genes)
  - Lipidomics (all 1,386 lipids)
  - Additional cell types (separate views)
- **Model selection:** Tested 5-15 factors, selected 12 factors (model12)
- **Factor 1 analysis:**
  - Strongest association with gonadal type
  - Explains most variance in hepatocyte RNA + lipidomics
  - Used to identify top contributing features

**Gene-Lipid Correlation Networks:**
- **Approach 1 (MOFA-based):**
  - Top 2.5% MOFA Factor 1 features:
    - 503 genes (hepatocyte RNA)
    - 31 lipids
    - 9 TFs (from TF activity scores)
  - Filter genes to "Cellular lipid metabolic process" (GO term)
  - For each lipid: find top 3 correlated genes (Pearson r > 0.65)
  - For each gene: find top 3 correlated TFs
  - **Result:** 94 nodes, 246 edges
  - **Key hubs:** Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9 (correlated with Nfic, Dbp, Esr1)

- **Approach 2 (Zonation-specific):**
  - Identify pericentral vs. periportal gonad-dependent genes
  - 4 categories:
    - Male-biased pericentral
    - Male-biased periportal
    - Female-biased pericentral
    - Female-biased periportal
  - For each gene: find correlated TFs + lipids (r > 0.65)
  - **Key genes:** Sult3a1, Abcb1a, Cyp2a4 (female pericentral), Mtmr11, Mindy3, Esco2, Kitl (male periportal)

### 3.4 Zonation Analyses

**Hepatocyte Zonation Assignment:**
- **Method:** Expression-based inference (Halpern et al. 2017 markers)
- **Pericentral:** Cyp2e1+ (Glul+ co-marker)
- **Periportal:** Cyp2f2+ (Gls2+, Hal+, Sds+ co-markers)

**Sex Differences in Zonation:**
- Higher periportal:pericentral ratio in males
- More portal veins + cholangiocytes in females (despite fewer periportal hepatocytes)

**Zonation-Specific Differential Expression:**
- **Contrasts:**
  - Testes vs. Ovaries (within pericentral)
  - Testes vs. Ovaries (within periportal)
  - Pericentral vs. Periportal (within sex)
- **Key findings:**
  - 35/82 gonad-dependent genes (43%) pericentral-biased: Cyp2c38, Cyp2b13, Cyp2d9
  - 18/82 genes (18%) periportal-biased: Cyp2b9

**TF Activity in Zonation:**
- Portal-biased: GATA6, GATA4, ONECUT1/HNF6, FOXA1, FOS, JUN, FXR, C/EBPα, STAT3
- Central-biased: PPARγ, PPARα, ATF4, PXR, CAR, LXR, NRF2

### 3.5 Validation Analyses

**GDX Correlation:**
- FCG (testes vs. ovaries) vs. GDX male (GDX vs. sham)
- Pearson r = -0.48 (expected inverse correlation)
- Scatter plots showing concordance

**Xenium Validation:**
- 82/82 gonad-dependent genes tested
- XX/82 significantly sex-biased in wild-type tissue
- Spatial visualization of zonation patterns

**Proteomics Validation:**
- ~33% of hepatocyte gonad-dependent genes also protein-level dimorphic
- FDR ≤ 0.05
- Validates Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9 at protein level

**siRNA Validation:**
- KD efficiency: Bulk RNA-seq (target expression reduced)
- Lipid consequences: Untargeted lipidomics
- "Validation score" metric: consistency across predicted lipid classes
- Top validated: Cyp2c38, Mtmr11

---

## 4. KEY OUTPUTS & FIGURES

### Figure 1: FCG Dataset Overview
- **A.** FCG model schematic
- **B.** UMAP showing sample contributions
- **C.** UMAP with cell type annotations
- **D.** Top cell type markers
- **E.** PCA (cell types separate, variance driven by cell type)
- **F.** Cell type proportions

**Output files:**
- `results/pseudobulk/ttop_full_list.RData` - All DE results
- `results/pseudobulk/tval_full_df.RData` - t-statistic matrix
- `results/pseudobulk/PCA_sample_coord_RNA_HEPA.csv`

### Figure 2: Gonadal Hormones Shape Gene Expression
- **A.** Ridgeline plots: PC1 by allosome type (per cell type)
- **B.** Ridgeline plots: PC1 by gonadal type (per cell type)
- **C.** Volcano plot: Hepatocyte allosome effect (few DEGs)
- **D.** Volcano plot: Hepatocyte gonadal effect (~4,000 DEGs)
- **E.** Pathway enrichment heatmap (hepatocytes)
- **F.** FCG vs. GDX scatter plot (r = -0.48)
- **G-I.** Xenium spatial validation

**Output files:**
- `results/pseudobulk/volcanos/` - Volcano plots per contrast × cell type
- `results/pseudobulk/pathways/pathway_scores_df_all.RData`

### Figure 3: Hepatic Lipid Abundance
- **A.** Lipidomics experimental cartoon
- **B.** PCA (samples separate by gonad, not chromosome)
- **C.** LINEX2 network: Top 5 testes + top 5 ovaries lipids
- **D.** LINEX2 subnetwork: TG/DG connections
- **E.** Lipid class enrichment lollipop plot
- **F.** Bubble plot: SM structural heterogeneity
- **G.** Bubble plot: TG structural heterogeneity
- **H.** FCG vs. GDX lipid correlation

**Output files:**
- `results/lipidomic/ttop_full_list.RData`
- `results/lipidomic/gonade_type_ttop.csv`
- `results/lipidomic/significant_lipids_network_SIF.csv`
- `results/lipidomic/significant_lipids_network_ATT.csv`

### Figure 4: Gene-Lipid Interactions (MOFA)
- **A.** MOFA → network approach cartoon
- **B.** MOFA results (factor scores, variance explained)
- **C.** Gene-lipid-TF network (focus on 4 Cyps)
- **D.** siRNA validation cartoon
- **E.** KD efficiency (RNA-seq)
- **F-G.** Lipidomics after KD (consistency scores)

**Output files:**
- `results/mofa/main_model_RNA_lipid_weight.RData`
- `results/multi-omic/hepatocyte_TF_spec_gene_lipid_network_mofa_full50percent.Rdata`

### Figure 5: Spatial Compartmentalization
- **A.** Zonation analysis cartoon
- **B.** Xenium: zonation-specific gene expression
- **C.** Xenium volcano plot (pericentral vs. periportal)
- **D.** Sex-specific zonation-dependent genes
- **E-F.** Zonation-specific subnetworks
- **G.** Validation cartoon
- **H.** Mtmr11 KD efficiency
- **I.** Mtmr11 KD lipidomics
- **J.** Model figure

**Output files:**
- `results/pseudobulk/zonation/` - Zonation-specific DE results
- `results/multi-omic/structural_model/` - Structural analysis results

### Supplementary Figures
- Genotyping validation (Xist, Y-linked genes)
- QC metrics (UMIs, genes, mitochondrial %)
- Cell type-specific DE results (all contrasts)
- GDX validation (all cell types)
- Proteomics validation
- siRNA full results (all targets)
- Lipid structural analysis (all genes)

---

## 5. ANALYSIS PRIORITIES

### Tier 1 (Core Analyses - Must Include)
1. **Pseudobulk differential analysis** (6 contrasts × 5 cell types)
   - gonadal_type, allosome_type, XXT_vs_XXO, XYT_vs_XYO, XYO_vs_XXO, XYT_vs_XXT
2. **Lipidomic differential analysis** (same 6 contrasts)
3. **Pathway enrichment** (GO + TF regulons)
4. **MOFA integration** (hepatocyte RNA + lipidomics)
5. **Gene-lipid correlation networks** (MOFA-based)

### Tier 2 (Important Validations)
6. **GDX analysis** (RNA + lipids)
7. **Zonation analysis** (pericentral vs. periportal)
8. **Zonation-specific gene-lipid networks**
9. **Proteomics comparison**

### Tier 3 (Advanced/Optional)
10. **Lipid structural analysis** (regression models)
11. **Xenium spatial analysis** (spatial patterns)
12. **siRNA validation analysis** (KD + lipidomics)

---

## 6. SOFTWARE & METHODS

### R Packages
- **Core:** Seurat v5, limma, dplyr, tidyr, ggplot2
- **Integration:** MOFA2
- **Enrichment:** decoupleR (ULM method)
- **Visualization:** pheatmap, ComplexHeatmap, ggrepel, cowplot, patchwork
- **Networks:** igraph (for LINEX2 network handling)

### Key Analysis Parameters
- **DE significance:** FDR ≤ 0.05
- **Correlation threshold:** Pearson r > 0.65
- **MOFA:** 12 factors, top 2.5% features
- **Pathway enrichment:** min set size 10-20 genes
- **Lipid structural models:** ≥10 lipids per class for class-specific models

### Data Formats
- **RNA:** Seurat objects → pseudobulk matrices → limma
- **Lipids:** MS-DIAL output → normalized intensity matrix
- **MOFA:** Multi-view list format
- **Networks:** SIF (edge list) + ATT (node attributes)

---

## 7. CRITICAL BIOLOGICAL FINDINGS TO REPLICATE

### Gene Expression
- **Gonadal effect:** >4,000 hepatocyte DEGs (FDR ≤ 0.05)
- **Chromosome effect:** Minimal (Y-linked + X-escapees only)
- **Cell type specificity:** Hepatocytes most responsive, cholangiocytes least
- **Top male-biased:** Cyp2d9, Mup genes, ribosomal genes
- **Top female-biased:** Cyp2b9, Cyp2b13, Cyp2c38, steroid metabolism genes

### Lipid Abundance
- **Gonadal effect:** 733/1,386 lipids (53.9%)
- **Chromosome effect:** 0/1,386 lipids
- **Male-biased classes:** TG, DG, glycerophospholipids
- **Female-biased classes:** SM, Cer, HexCer

### Integration
- **Key regulatory genes:** Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9
- **Key TFs:** Nfic, Dbp, Esr1 (correlated with Cyps)
- **Zonation patterns:**
  - Pericentral: Cyp2c38, Cyp2b13, Cyp2d9 (43% of gonad-dep genes)
  - Periportal: Cyp2b9 (18% of gonad-dep genes)
  - Zonation-specific regulators: Mtmr11 (male periportal), Sult3a1 (female pericentral)

---

## 8. CHECKLIST FOR WORKFLOW VALIDATION

When building/testing workflows, verify:

- [ ] All 6 main contrasts are run (not just 3)
- [ ] All 5 cell types are analyzed (not just hepatocytes)
- [ ] Gonadal effect shows ~4,000 hepatocyte DEGs
- [ ] Allosome effect shows minimal DEGs
- [ ] 733 gonad-dependent lipids identified
- [ ] 0 chromosome-dependent lipids
- [ ] MOFA Factor 1 associates with gonadal type
- [ ] Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9 emerge as top genes
- [ ] Nfic, Dbp, Esr1 identified as top correlated TFs
- [ ] Zonation analysis identifies pericentral/periportal biases
- [ ] GDX correlation is negative (r ≈ -0.48)
- [ ] All expected output files are generated

---

**Last updated:** 2025-01-28 (based on manuscript draft)
