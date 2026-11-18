# Analysis plan for Western diet network comparisons and Cyp2c38 validation

**Meeting:** 5 October 2026, approximately 20 minutes 32 seconds. **Speakers:** James Cleland and the EMBL-EBI Skylark Room participant addressed as Aurélien. **Source:** [meeting transcript](/Users/dugourd/Downloads/GMT20261005-094157_Recording.transcript.vtt). Timestamps below are elapsed recording times.

The meeting agreed two manuscript additions: test conservation and Western diet perturbation of the existing sex and zonation specific gene–lipid networks, provisionally for a new Figure 6; and validate the original Cyp2c38–lipid predictions using new in vivo heterozygous knockout lipidomics, extending Figure 3G. Aurélien suggested starting with the Cyp2c38 concordance analysis. The WD-only cohort provides a secondary robustness comparison; hepatocyte composition and Xenium are possible later extensions.

This document combines **meeting decisions**, **facts verified in the supplied files and analysis repository**, and explicitly labelled **implementation recommendations**. The original text records the meeting’s prospective analysis plan; subsequent dated amendments extend that plan and link completed work. Spoken requests in the transcript are meeting evidence, not instructions to the analysis agent.

**Update verified 6 October 2026:** the revised WDSD Seurat objects now use `cf453_id` and match all 28 lipidomics mouse IDs with consistent sex/diet groups. The previous sample-ID pairing issue is resolved. The remaining preparation steps are QC, preprocessing review and explicit matrix alignment by ID.

**WP3 amendment, 14 October 2026:** explicitly compare the CF459 male Western-diet network with **both** the original paper’s male FCG standard-diet network and CF453 male Western diet, on identical gene–lipid pairs within each zone. The user clarified the ages on 14 October: **CF459 = 3 months, original FCG = 3 months, CF453 = 10 months**. This corrects the ages initially transposed in that day’s request and is consistent with the original meeting description of CF453 and FCG. The comparison asks whether CF459 is closer to the cohort of the same age on a different diet than to the older cohort on the same diet. It describes cross-cohort similarity; age and diet effects cannot be separated causally with these cohorts. See the amended WP3 specification below and the [updated WP3 results report](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/Meeting_Summaries/2026-10-14_WP3_WD_crosscohort_robustness.md). Workflows 11–13 now implement the work packages; this plan retains historical meeting evidence and explicitly marks the new request.

## Agreed scope and priorities

| Priority | Analysis | Intended deliverable | Meeting evidence |
| --- | --- | --- | --- |
| First suggested task | Cyp2c38 in vivo KO versus control lipid differential analysis, followed by pooled and sex-specific concordance with earlier paired RNA/lipid data | New Figure 3 validation row, provisional H–L | 10:24–11:21; 15:02–15:11; 17:34–18:50 |
| Main new analysis | Reproduce Figure 5 gene–lipid associations in WDSD standard-diet controls, separately by sex and hepatocyte zone | Baseline conservation scatterplots | 00:03–03:58 |
| Main new analysis, after baseline assessment | Compare WDSD Western versus standard diet within each sex and zone; identify the associations that change most | Proposed Figure 6 network perturbation panels and edge tables | 03:58–05:26; 08:29–08:44 |
| Secondary robustness analysis, expanded 14 October | Compare male CF459 WD (3 months) with CF453 WD (10 months) **and original male FCG SD (3 months)**, using identical pairs within each zone; assess which anchor is closer to CF459 | Matched correlation panels, direct difference in network agreement and mouse-level stability | 06:43–07:40; explicit user request and age clarification, 14 October |
| Optional extensions | Joint analysis of the four sex × diet groups; hepatocyte subtype proportions from Figure 4; Xenium if needed | Further Figure 6 or supplementary material | 04:42–05:26; 08:04–08:29 |

The WD-only scope evolved: James initially suggested keeping it in reserve and not analysing it yet, then endorsed the proposed cross-cohort Western-diet comparison and proactive inclusion of a sub-panel. Treat this as secondary work after the two main priorities, rather than an unconditional exclusion or an equally urgent first task. The network comparisons alone were considered sufficient for the first version of Figure 6.

## Dataset register and current readiness

The analysis repository is `/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored`. The dataset and meeting-summary directories are siblings of that repository under `FCGL_FirstSubmission`.

| Dataset | Exact supplied files | Verified contents and role |
| --- | --- | --- |
| Cyp2c38 KO lipidomics | [lipidomics_cyp2c38ko.csv](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/cyp2c38ko_lipidomics/lipidomics_cyp2c38ko.csv) | Primary KO measurement input: 1,319 feature rows, 32 biological sample columns, 9 QC columns and 32 annotation columns. Comma-delimited. |
| WDSD lipidomics, CF453 | [lipidomics_cf453.csv](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/wdsd_lipidomics/lipidomics_cf453.csv) | Main sex × diet experiment: 1,240 feature rows, 28 biological sample columns, 7 QC columns and 32 annotation columns. Comma-delimited. |
| WDSD snRNA-seq | [wdsd_hep_pseudobulks.rds](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/wdsd_snrnaseq/wdsd_hep_pseudobulks.rds); [wdsd_zonation_pseudobulks.rds](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/wdsd_snrnaseq/wdsd_zonation_pseudobulks.rds) | Seurat v5 pseudobulks: 28 hepatocyte profiles and 56 zonation profiles, corresponding to 28 mice × 2 zones. Metadata has 7 mice in each sex × diet group. |
| WD-only lipidomics, CF459 | [lipidomics_wd.csv](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/wd_lipidomics/lipidomics_wd.csv) | Secondary cohort: 1,463 feature rows and 8 biological sample columns, plus `Alignment ID` and `Metabolite name`. **Semicolon-delimited with decimal commas and a UTF-8 BOM** despite the `.csv` extension. No QC columns in this supplied export. |
| WD-only snRNA-seq | [wd_pseudobulks.rds](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/wd_snrnaseq/wd_pseudobulks.rds); [wd_zonation_pseudobulks.rds](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/wd_snrnaseq/wd_zonation_pseudobulks.rds) | 40 cell-type pseudobulks (8 mice × 5 cell types) and 16 zonation pseudobulks (8 mice × 2 zones). Metadata identifies 4 female and 4 male mice. Individual mouse IDs match the CF459 lipid export. |

**Counting convention:** feature-row totals are supplied annotation/measurement rows before analysis filtering, not final counts of retained endogenous lipids. The KO and WDSD tables include internal-standard rows such as `IS: Cer C15 (d7)`. QC injections are analytical controls, not biological replicates. Presence of processed pseudobulks does not establish that detailed nucleus-level QC has passed.

### Verified import and RNA layer details

- KO contains **16 internal-standard rows** and QC03–QC11; simply removing standards leaves 1,303 feature rows before other filters. WDSD contains **14 standards** and QC04–QC10; simply removing standards leaves 1,226 rows. WD-only has no internal-standard rows in this export and omits adduct, RT/m/z and most other annotation columns.
- All three lipid exports have unique nonblank `Metabolite name` strings, but **`Alignment ID` is not unique**: duplicate-ID excess rows are 43 in KO, 67 in WDSD and 244 in WD-only. KO/WDSD `(Alignment ID, Adduct type)` combinations are unique. Preserve a unique feature-row key and annotation; do not join datasets by alignment ID alone.
- Lipid quantitative cells were numeric and positive after correct locale parsing, with no blanks/NA or nonfinite values found in the supplied exports. This does not establish absence of upstream missing values or imputation. Confirm filtering, normalization and imputation provenance. Names can be pipe-separated, for example `Cer 34:1;O2|Cer 18:1;O2/16:0`; preserve structural resolution during harmonization.
- For WD-only, use a decimal-comma reader such as R `read.csv2()` or `readr::read_delim(delim = ";", locale = readr::locale(decimal_mark = ","))`. A standard comma-delimited or decimal-point parser is incorrect.
- All four RDS objects record **Seurat 5.3.0**, one RNA **Assay5**, and existing `counts`, `data` and `scale.data` layers. Each layer has **20,123 genes**, with columns respectively 40/16 for WD cell types/zones and 28/56 for WDSD hepatocytes/zones. `counts` is sparse, integer and nonnegative; `data` is dense and nonnegative; `scale.data` is dense and signed. No layer NA values were found. Scaled expression is already present; it need not be regenerated merely to start these comparisons.
- WD cell-type metadata contains `orig.ident,sample_name,cell_type,group`; the five cell types are Hepatocytes, Endothelial cells, Immune cells, Stellate cells and Cholangiocytes. Zonation metadata uses `zonation_type` instead of `cell_type`. Updated WDSD hepatocyte metadata is `orig.ident,cf453_id,group`; WDSD zonation adds `zonation_type`. The updated WDSD objects use `cf453_id` instead of the former `sample_name` field. Groups are `female-sd,female-wd,male-sd,male-wd` for WDSD and `f-ctrl,m-ctrl` for WD.
- Only `ScaleData.RNA` is recorded in the RDS command histories, so the exact upstream normalization cannot be recovered from that history alone. Use an appropriate Seurat v5 environment for implementation. The current local R runtime has `SeuratObject` available, verified on 6 October 2026. Metadata/layer inventory and sample-map checks were performed read-only.

The transcript only explicitly stated **7 male standard-diet mice**. Counts for the other groups below come from file headers/RDS metadata, not from the meeting.

| Cohort | Group | Biological mice in supplied files |
| --- | --- | ---: |
| WDSD / CF453 | Female standard diet | 7 |
| WDSD / CF453 | Female Western diet | 7 |
| WDSD / CF453 | Male standard diet | 7 |
| WDSD / CF453 | Male Western diet | 7 |
| WD-only / CF459 | Female, RNA metadata `f-ctrl` | 4 |
| WD-only / CF459 | Male, RNA metadata `m-ctrl` | 4 |
| Cyp2c38 | Female control, `F-Ctrl` | 8 |
| Cyp2c38 | Female KO, `F-KO` | 8 |
| Cyp2c38 | Male control, `M-Ctrl` | 7 |
| Cyp2c38 | Male KO, `M-KO` | 9 |

For CF459, the transcript establishes Western diet only; the metadata label `ctrl` must not be interpreted as standard diet. For Cyp2c38, confirm that the column label `Ctrl` means the WT genotype described in the meeting. The KO table contains 17 KO and 15 control mice overall, with unequal sex composition across genotype groups.

### WDSD sample pairing verified

The updated WDSD Seurat files, modified on **6 October 2026**, identify mice using `cf453_id` values `cf453-1` through `cf453-28`. All 28 match the biological lipidomics columns, with **zero unmatched IDs and zero sex/diet disagreements**: IDs 1–7 are female SD; 8–14 female WD; 15–21 male SD; and 22–28 male WD. Hepatocyte and zonation metadata agree for every mouse; the zonation object contains exactly one Pericentral and one Periportal profile per mouse. Metadata row names, `orig.ident`, assay sample maps and active identity names are consistent. No unresolved animal-ID crosswalk task remains.

For implementation, obtain lipid mouse IDs with `tolower(sub("_.*$", "", lipid_sample_columns))` for the biological `CF453-*` columns, then join to RNA metadata `cf453_id`. **RNA and lipid column orders differ**, so align explicitly by ID before computing correlations. The stored Assay5 layer matrices use sample maps rather than bare matrix dimnames; use Seurat v5 layer accessors to export named matrices. These are routine alignment requirements, not an unresolved crosswalk.

WD-only IDs are `cf459-1`, `2`, `3`, `4`, `13`, `14`, `15`, `16` in both modalities. Validate sample order explicitly even where IDs match.

## Repository analyses to reuse

The current refactored workflows, rather than legacy script names in historical project notes, are the relevant starting points. The [agent guide](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/AGENT_GUIDE.md) and [cross-workflow overview](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/workflow_plans/_meta_cross_workflow_overview.md) map their dependencies.

| Existing analysis | Reuse and verified distinction |
| --- | --- |
| [Workflow 03 lipidomics](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/scripts/workflows/03_fcg_lipidomic_analysis.R) | Lipid QC, limma differential analysis, lipid identity/class handling and plotting conventions. The standard FCG loader expects a sample-by-lipid table and applies log2; the new annotated feature-by-sample files need a different import adapter. |
| [Workflow 04 MOFA integration](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/scripts/workflows/04_mofa_integration.R:925) | Supplies the global Cyp2c38 gene–lipid reference used by existing validation. The full gene–lipid network pools matched hepatocyte samples across gonads; it is not the Figure 5 sex-specific zonation network. |
| [Workflow 05 zonation](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/scripts/workflows/05_zonation_analysis.R:1339) | Section 13 builds four gonad × zone networks: testes periportal, testes pericentral, ovaries periportal and ovaries pericentral. Each contains gene–lipid and TF–gene edges. Standard networks use absolute gene–lipid correlation cutoff 0.65; full exports retain gene–lipid coefficients without that cutoff. |
| [Workflow 07 validation](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/scripts/workflows/07_validation_experiments.R:227) | Sections 6–7 compare earlier gene–lipid predictions with siRNA lipid t-statistics using ULM/scatterplots; sections 8–9 do zonation validation. This is the analytical precedent for repeating Figure 3G with in vivo Cyp2c38 perturbation. |
| [Seurat compatibility helpers](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/scripts/seurat_compatibility.R) | Extract `counts`, `data` or `scale.data` through v5 layer access; inspect assays and layers before exporting new pseudobulk matrices. |

**Reference network files:**

- [Global hepatocyte full gene–lipid network](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/mofa/networks/hepatocyte_gene_lipid_network_full_SIF.csv): columns `node1` = gene, `node2` = lipid, `mor` = Pearson correlation coefficient. Extract the `Cyp2c38` rows for pooled concordance. Despite its name, `mor` here is an association coefficient, not a validated biochemical regulatory sign.
- [Testes periportal full network](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation/networks/cormat_testes_portal_full_SIF.csv).
- [Testes pericentral full network](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation/networks/cormat_testes_central_full_SIF.csv).
- [Ovaries periportal full network](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation/networks/cormat_Ovaries_portal_full_SIF.csv).
- [Ovaries pericentral full network](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/results/zonation/networks/cormat_ova_central_full_SIF.csv).

The zonation files have `node1,node2,weight`; use node-type information or validated lipid identity to isolate **gene–lipid** edges. The mixed `Ovaries` capitalization and `ova` abbreviation are actual filenames. The `compressed_cormat_*` networks aggregate lipid classes and signs; use individual-species coefficients for the proposed species-level comparisons.

The helper's TF branch retains a secondary absolute-correlation cutoff of 0.5 even in the so-called full exports; “full” applies to the unthresholded gene–lipid coefficient set here. Original zone-specific genes are selected from gonad-by-zone interaction statistics and within-gonad zonation direction. Preserve those context-specific selections for the initial replication analysis.

The source matrices needed for sex-specific Cyp2c38 recomputation are the three `scaled_RNA_1/2/3.csv` files and `meta_data_df.csv` under the repository's [RNA input directory](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/data/RNA), with the [original lipid measurements](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/data/metabolomics/Lipid_data_file_JPC_cleaned.csv). Figure 5 reference reconstruction uses `scaled_RNA_zonation.csv` and `meta_data_df_zonation.csv` in that RNA directory. Reuse the current loader's documented sample normalization and explicit matching; adapt the new cohort import rather than applying old `Testes`/`Ovaries` suffix regexes directly to `female-sd`/`male-wd` labels.

**Verified reference sample sizes:** the current FCG hepatocyte input used by Workflow 04 matches 43 mice to lipidomics: 22 testes-bearing and 21 ovary-bearing. Each Workflow 05 zone uses 22 testes-bearing or 21 ovary-bearing mice. Do not label every sex-specific scatter as n = 43. In the FCG reference, “male/female” shorthand must be mapped to **gonadal type**; the new wild-type cohorts are grouped by recorded sex.

**Feature-selection detail:** the global full network selects the highest and lowest **33% of Factor 1 loadings separately**, approximately 66% of features in total, then restricts genes to the cellular lipid metabolic process set. It contains all retained gene–lipid pair coefficients without an association cutoff. The sparse MOFA network uses 2.5% from each tail and is a different export. Freeze the established reference feature universe for validation; retraining MOFA on the new cohorts was not requested.

## Work package 1 Cyp2c38 in vivo lipidomics and concordance

### Experimental facts and figure context

**Meeting decisions and observations, 10:24–14:56:** this is an in vivo **heterozygous Cyp2c38 knockout**, started roughly two years earlier after the first validation round. Only Cyp2c38 was targeted. Residual protein is expected; James described the perturbation as effectively a knockdown. The animals were described as young and healthy, with no exact age specified.

The displayed proteomics panel, provisionally **Figure 3I**, had groups ordered female KO, female WT, male KO, male WT. James reported higher Cyp2c38 protein in females than males and lower protein in KO than WT in both sexes. Frank's proteomics pipeline generated the supporting analysis. This establishes a protein reduction in the displayed data, not complete protein ablation. The separate [Cyp2c38 proteomics results directory](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/cyp2c38ko_proteomics/JamesCleland_P4562_QuantUMS_PGlevel1_V1/data_analysis_results_V1) is present locally; its `Limma_results_V1.csv` has Cyp2c38 reductions in both sex-specific contrasts. It is supporting material outside the five primary input directories supplied for this request.

The existing [Cyp2c38 protein plot](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/original_datasets/cyp2c38ko_proteomics/JamesCleland_P4562_QuantUMS_PGlevel1_V1/data_analysis_results_V1/single_protein_quant/Cyp2c38_V1.pdf), `Full_dataset_V1.csv` and `Methods_V1.txt` provide figure and processing provenance. The protein export's group counts agree with the lipid counts, but individual protein/lipid sample pairing was not established in this summary.

James reported no observed body-weight or liver-weight effect in panels J/K and suggested placing those measurements in the supplement. The exact assignment of J versus K was not stated. A possible stronger phenotype after **six months of Western diet** was a hypothesis for a future challenge, not an agreed experiment or an available dataset. Absence of those weight effects does not establish absence of molecular effects or all phenotypes.

Figure references were corrected during the discussion to **Figure 3G** as the panel to repeat. A new bottom row, **H–L probably**, was provisional; confirm final lettering against the current manuscript artwork. The historical validation described loosely as an in vitro “cell line knockout” is documented in the repository as **siRNA knockdown in primary hepatocytes**.

### Input preparation and differential analysis

**Agreed input handling, 15:10–16:24:** use `lipidomics_cyp2c38ko.csv`. James said it combines positive/negative acquisition modes and removes duplicate identifications, using the previous cleaning procedure. The four accompanying `statistical results F-Ctrl/M-Ctrl vs F-KO/M-KO pos/neg mode.txt` files were for quick preliminary statistics; both speakers agreed they were unnecessary as primary analysis inputs. Recompute differential effects from measurements.

For provenance, the four files are `statistical results F-Ctrl vs F-KO neg mode.txt`, `statistical results F-Ctrl vs F-KO pos mode.txt`, `statistical results M-Ctrl vs M-KO neg mode.txt` and `statistical results M-Ctrl vs M-KO pos mode.txt`. Negative-mode files have 649 rows each and positive-mode files 929 each. Their mode-specific feature totals differ from the integrated CSV. Their generic `Statistical Value` has not been verified as a limma moderated t-statistic, and `fold change [Ctrl] vs [KO]` direction must not be assumed to match the planned KO-minus-WT contrast.

- [ ] Create explicit sample metadata from `Cyp-*_F-Ctrl`, `F-KO`, `M-Ctrl` and `M-KO`, retaining mouse ID, sex, genotype and any known technical batch. Verify the 8/8/7/9 group counts and whether any mice were excluded before export.
- [ ] Separate `Metabolite name`, annotation fields, biological measurements, internal standards and QC injections. Produce a numeric lipid × mouse matrix with stable feature IDs and an annotation cross-reference.
- [ ] Check intensity scale, missing/nonpositive values, normalization provenance, feature duplication, QC performance and sample PCA. Harmonize lipid names with the reference without silently merging different structural identifications.
- [ ] Fit the three agreed lipid contrasts: **pooled KO versus WT**, **female KO versus female WT**, and **male KO versus male WT**. Save complete tables with lipid ID, logFC, moderated t, P value and FDR, plus a t-statistic matrix.

**Implementation recommendation:** interpret the spoken “confound the sexes” as pooling the sexes for an overall genotype estimate, while retaining sex in the design. Because male KO/control counts differ and baseline lipid abundances can differ by sex, a naive KO-versus-control pooled model may mix genotype and sex composition. Use an explicit sex × genotype design; a transparent pooled contrast is `0.5 × [(F_KO − F_WT) + (M_KO − M_WT)]`, alongside the two sex-specific contrasts. A sex × genotype interaction is a useful secondary analysis, but was not specified in the meeting. Define `KO − WT` consistently so positive t means increased abundance after perturbation.

**Normalization decision to document:** the main FCG lipid loader applies log2, whereas the current Workflow 07 knockdown path does not explicitly log-transform after `clean_lipid_data()`. That helper reformats/filter names; it does not by itself establish an appropriate transformation for these new measurements. Confirm the export scale and select preprocessing before fitting, recording any difference from the existing validation analysis. Do not log-transform already logged data or import the annotation columns as measurements.

### Concordance analyses

**Agreed, 16:24–18:50:** no gene-expression profiling exists for the KO animals and none is required here. Prediction correlations come from the **earlier paired RNA/lipid cohort**, not from new KO RNA.

1. **Pooled prediction versus pooled KO effect:** extract the original global Cyp2c38–lipid coefficients and join to the pooled KO t-statistics by harmonized lipid identity. Reproduce the Figure 3G scatter: one point per shared lipid, reference correlation on x, observed KO t-statistic on y.
2. **Sex-specific effects against the pooled reference:** the discussion allowed examining female and male KO effects against the same pooled reference as an intermediate comparison. Retain these as secondary/contextual plots if useful.
3. **Sex-matched validation:** recompute Cyp2c38–lipid correlations separately in the original testes-bearing and ovary-bearing FCG mice, using the established paired hepatocyte expression/lipid input and fixed reference feature set. Compare male KO t-statistics with the testes-bearing reference, and female KO t-statistics with the ovary-bearing reference.
4. Export joined per-lipid data, overlap/missing-feature reports and concordance summaries for each comparison. Retain the original in vitro Cyp2c38 result as context; a direct new in vivo versus old in vitro effect comparison would be an additional analysis, not an explicit meeting requirement.

**Repository sign convention:** Workflow 07 computes `consistency = −(reference correlation × perturbation t-statistic)`. With `KO − WT`, concordant perturbations therefore tend to give a **negative slope** between the original gene–lipid correlation and observed t-statistic. Positive reference association predicts a lipid decrease when Cyp2c38 is reduced; negative association predicts an increase. Preserve the sign in labels and summaries. If using the existing decoupleR ULM against signed reference coefficients, concordance corresponds to a negative score with this contrast orientation; do not label a positive score as concordance without redefining the sign.

**Implementation recommendation:** use all eligible shared reference lipids for the main concordance calculation rather than selecting only significant KO responders. Report the number of shared lipids, effect association/slope, uncertainty and the signed consistency measure. Report discordance as well as concordance; the meeting explicitly left the outcome open. Correlation-based predictions and heterozygous perturbation effects need not agree perfectly.

## Work package 2 WDSD baseline conservation and Western diet perturbation

### Cohort and biological question

**Meeting facts, 00:03–02:44:** CF453 is the lowercase, unhyphenated `wdsd` dataset. It contains wild-type male/female animals on standard/Western diet, approximately **10 months old after 9 months of diet exposure**. The original manuscript cohort was described as **3 months old**. The new background was described as half C57BL/6 and half another subspecies; the exact strain/cross and subspecies were not named. Exact replication is therefore not expected.

James had already used the previous code to generate periportal/pericentral hepatocyte annotations but had **not performed detailed QC**. RNA had been available for some time; lipidomics had arrived approximately a month before the meeting. The central question is whether Figure 5 associations, including the partially validated **Mtmr11** network, are conserved and then altered by Western diet.

### Preparation required before network computation

- [ ] Inspect RNA metadata, assays and available layers; verify preprocessing, gene symbols, sex/diet labels and zone annotation. Check adequate contribution of nuclei to each pseudobulk and document any exclusions. Repository zonation markers are `Cyp2e1` high/`Cyp2f2` low for pericentral and the reverse for periportal; verify the supplied annotations used the intended rule and expression scale.
- [ ] Export paired periportal and pericentral gene × mouse matrices and a lipid × mouse matrix, with an explicit mouse-order check. Each mouse contributes two RNA zones but one bulk-liver lipid profile. Do not count its two zones as two independent animals.
- [ ] Apply a documented lipid import/QC procedure, excluding standards and QC samples from biological comparisons. Check PCA by sex and diet; review batch, outlier, detection and normalization information.
- [ ] Freeze the original Figure 5 context-specific genes, lipid species and edge list for the primary replication analysis. Audit which are measurable in every relevant cohort and save unmatched IDs and exclusions.

**Repository context:** Workflow 05 uses scaled zonation RNA and lipid abundances matched to the same mouse, builds correlations separately by gonad and zone, and includes inferred TF activities in its broader networks. Reuse its Pearson gene–lipid association logic for the meeting's primary comparisons. TF–gene replication could extend the analysis later, but should be labelled separately. Matching whole-liver lipid abundance to zonal RNA gives a zonation-associated expression–lipid relationship; it does not measure lipid abundance separately in the two zones.

### Stage A Assess conservation in standard diet

**Agreed, 02:44–03:58:** for the same gene–lipid pairs, compute correlations in the matching WDSD standard-diet stratum and compare them with original Figure 5 coefficients. Mtmr11 was the concrete example. Make scatterplots and calculate the association between the two vectors of coefficients, a **correlation of correlations**.

| Original reference | New comparison |
| --- | --- |
| Testes periportal | WDSD male SD periportal |
| Testes pericentral | WDSD male SD pericentral |
| Ovaries periportal | WDSD female SD periportal |
| Ovaries pericentral | WDSD female SD pericentral |

The illustrative values **0.8 or 0.9** discussed in the meeting were examples of strong conservation, not observed results or a predetermined pass threshold. Assess conservation honestly given age/background differences. If baseline conservation is weak, document that limitation before interpreting Western-diet departures from the original manuscript network.

**Implementation recommendations:** retain a primary panel based on the original displayed gene–lipid edges, and a complementary analysis of all eligible pairs in the original full context-specific export. Do not require an edge to exceed 0.65 again in the new data before evaluating it; doing so removes weakened/disappearing associations. Use a consistent common feature universe when directly comparing summary statistics, report both edge count and mouse count, and record coefficients that cannot be estimated because of missing or invariant features.

### Stage B Assess Western diet perturbation within WDSD

**Agreed, 03:58–05:26:** compare male SD versus male WD periportal coefficients as the worked example, and apply the same approach across sex/zone contexts. Start with scatterplots and identify the outliers: gene–lipid correlations that change most.

- [ ] Compute SD and WD coefficients for the same predefined pairs within male periportal, male pericentral, female periportal and female pericentral contexts.
- [ ] Plot `r_SD` against `r_WD` and annotate selected changing associations, including Mtmr11 where measurable.
- [ ] Export an edge-change table with sex, zone, gene, lipid, `r_SD`, `r_WD`, `Δr = r_WD − r_SD`, absolute change, sign change, and actual paired mouse counts.
- [ ] Rank and examine large changes; describe weakening, strengthening and sign reversal explicitly. Organize the strongest patterns by regulator, zone and lipid class without replacing species-level evidence with compressed class averages.
- [ ] Assemble a first Figure 6 draft: cohort/design context, standard-diet conservation, Western-diet comparison and selected changed associations. Exact panel number and lettering remain open.

**Implementation recommendations:** seven mice per WDSD group gives noisy individual correlation estimates. Use mouse-level sensitivity checks, such as leave-one-mouse-out estimates or an appropriate resampling procedure, before interpreting large `Δr` as reproducible rewiring. Correlation edges share genes/lipids and are dependent; a significance calculation treating thousands of edges as independent biological replicates would be misleading. Rank changes descriptively unless an explicit differential-association test and multiple-testing procedure have been defined.

### Optional integrative and composition analyses

**Discussed but not specified:** an analysis combining all four sex × diet groups, rather than only pairwise comparisons. A candidate implementation is a joint sex/diet interaction analysis of associations; a sex × diet lipid/RNA abundance model could provide complementary context. Changes in abundance are a different question from changes in correlation. No particular model was agreed.

James also suggested repeating **Figure 4 hepatocyte-type proportion analyses** and offered **Xenium** if needed. The supplied files are pseudobulks; subtype/nucleus proportions require annotated nucleus counts or appropriate per-mouse composition metadata, not simply the number of pseudobulk columns. Establish those inputs before scheduling composition analysis. Xenium is a possible experimental extension for James, not a prerequisite for the first network comparison draft.

## Work package 3 WD-only cross-cohort robustness

**Meeting facts, 06:05–09:08:** CF459 is `wd`, contains only Western-diet animals, has fewer replicates and was described as closer to the manuscript mice in age and genetic background. Exact age and background were not stated in the meeting. Both modalities are now available. **Dated metadata clarification from the user, 14 October 2026:** CF459 is **3 months old**, the original FCG cohort is **3 months old**, and CF453 is **10 months old**. These ages are supplied cohort descriptions, not values inferred from the raw exports. Diet formulation and CF459 diet-exposure duration, exact background, litter and processing/batch differences remain to be established. The meeting’s description of CF453 as approximately 10 months after 9 months of diet is preserved.

**Original agreed comparison:** reproduce the male WD network in CF459 and compare its coefficients with male WD CF453. **Explicit extension requested on 14 October:** compare CF459 also with the original paper’s male standard-diet cohort, and directly ask whether the 3-month WD network is closer to 3-month SD than to 10-month WD. Retain periportal and pericentral distinctions. Female comparisons remain an optional extension.

### Matched three-cohort comparison, specified for the amended analysis

1. Verify sex, WD designation and individual RNA/lipid pairing; parse the CF459 lipid table with semicolons and decimal commas. Use four CF459 males, seven CF453 WD males and the **22 testes-bearing FCG mice (11 XXT and 11 XYT)** in the original Figure 5 male networks. The FCG anchor is neither the five separate WT males nor only the XYT subgroup. Each animal supplies paired RNA and whole-liver lipid data; RNA zones are not independent animals.
2. Reconstruct the frozen original male Figure 5 coefficients from the original paired matrices and verify their numerical replay. Preserve the fixed genes and strict lipid identities. Use exactly the same eligible pairs for both anchors: **610 lipids; 25 periportal genes/15,250 full pairs, and 30 pericentral genes/18,300 full pairs**. All eligible original full pairs are primary. The original displayed subsets (**69 periportal, 202 pericentral**) are secondary and remain conditional on original FCG selection. Do not reselect edges by a new correlation threshold or remove unstable genes after viewing results.
3. For each zone and universe, report `cor(CF459, FCG)`, `cor(CF459, CF453)` and `cor(FCG, CF453)`, where each input vector contains the gene–lipid Pearson coefficients. Define the direct comparison as **`delta_r = cor(CF459, FCG) − cor(CF459, CF453)`**. A positive difference means more positive agreement with FCG on the same pairs; it does not imply strong absolute conservation. Plot both CF459-anchor comparisons on matched axes, showing edge counts and biological mouse counts.
4. Report native ULM scores and nominal ULM p-values alongside mouse-level checks. For the FCG anchor, reuse the same CF459 whole-mouse exact 24 RNA/lipid pairings and 4,999 empirical-marginal null-bootstrap indices as the original WP3 comparison. State the null hypothesis and valid/requested denominator; edge-vector p-values are descriptive because edges share animals and features. Record mathematically degenerate or undefined frozen-vector draws explicitly; never turn roundoff into a correlation or delete within-draw edges to recover a result.
5. Assess the **direct difference** with 999 paired three-cohort resamples: use the same CF459 paired mouse draw against both anchors, independently resample CF453 paired mice, and independently resample FCG paired mice within its 11 XXT/11 XYT strata. Export the seed and every index. Report how often the difference is positive among defined draws, together with undefined counts. This is a descriptive stability fraction, not a null p-value. Do not treat the two correlations as independent or use an independent-correlation Fisher test.
6. Check individual-mouse omissions in all three cohorts and the original CF453/CF459 joint omissions. Retain the frozen edge universe throughout. Include Mtmr11 as the predefined periportal example regardless of whether its ranking matches the full network; report gene-level exceptions.
7. Export matched per-edge coefficients, coverage/missing-feature ledgers, inference and resampling/omission tables under `results/wd_crosscohort_robustness/fcg_comparison/`. Preserve the initial WP3 numerical outputs and distinguish any new numerical-degeneracy correction from the historical analysis. Update this research plan and the main WP3 Markdown/PDF report in `../Meeting_Summaries`, archiving prior document versions.

There is **no within-CF459 SD-versus-WD contrast**, because no standard-diet CF459 arm is supplied. The new three-cohort comparison can establish measured relative similarity, conditional on the fixed genes, lipids and available mice. It cannot identify an age effect independently of diet, background, exposure, batch or other cohort differences. Agreement with one anchor, weak agreement with both anchors and a gene-specific exception are all valid outcomes; report magnitude and stability honestly.

**Reproducible entry:** [Workflow 13 plan](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/workflow_plans/13_wd_crosscohort_robustness.md) documents the implementation and test outputs. The dated results report records the actual findings; this section defines the expanded scientific question and completion criteria.

## Deliverables and completion criteria

The names below define required outputs; execution status and numerical findings are documented in the dated work-package reports. Keep new outputs in cohort-specific locations and preserve the existing FCG and siRNA reference results.

| Deliverable | Required contents |
| --- | --- |
| Dataset manifest and QC record | Original file names, biological/QC sample counts, annotation schema, RNA layers, sex/diet/genotype metadata, exclusions and preprocessing choices |
| WDSD paired sample manifest | Use the verified `cf453_id` mapping to record original lipid column, RNA profile, sex/diet, zone and final analysis order |
| Lipid identity crosswalk | Original annotation, harmonized ID, structural resolution, matched reference ID and unresolved mappings |
| Cyp2c38 differential tables | Pooled, female and male KO-minus-WT contrasts; logFC, moderated t, P value/FDR and actual n |
| Cyp2c38 concordance data and figures | Pooled and sex-matched reference coefficients joined to KO effects, signed consistency, overlap counts and comparison summaries |
| WDSD network coefficients | Mouse-paired coefficients per sex × diet × zone for a documented fixed feature universe |
| Baseline conservation figures | Four original-versus-WDSD SD context comparisons, coefficient-vector association and both edge/mouse counts |
| Western-diet edge-change report | Four within-WDSD comparisons, ranked `Δr`, sign changes, missing/invariant features and stability checks |
| CF459 three-cohort comparison | Matched original FCG SD/CF459 WD and CF453 WD/CF459 WD coefficients, direct delta r, full/displayed universes, 22/7/4 mouse counts, paired three-cohort stability, omissions and Mtmr11 exception assessment |
| Manuscript figure drafts | Figure 3 new validation row and first Figure 6 network comparisons; CF459 robustness including the same-age SD anchor as secondary material |
| Reproducibility notes | Code entry points, input provenance, sign conventions, software/configuration and limits on biological interpretation |

Completion requires faithful reporting regardless of whether conservation/concordance is strong, weak or absent. The meeting expressed hoped-for results but did not establish any new outcome.

## Owners dependencies and unresolved information

Ownership below follows the discussion; these are proposed work assignments, not additional commitments obtained after the meeting.

| Action or information | Suggested owner | What it enables |
| --- | --- | --- |
| Confirm strain/cross, diet formulation/exposure and batch information for CF453/CF459; cohort ages clarified by user on 14 October (CF459/FCG 3 months, CF453 10 months) | James / supplied user clarification | Accurate cohort descriptions and interpretation of replication differences |
| Confirm Cyp2c38 heterozygous genotype, WT control labels, age/background and exclusions | James | Correct KO design and sample accounting |
| Confirm normalization/cleaning provenance and remaining lipid/QC filtering | James with analytical processing records; Aurélien to implement | Consistent, auditable lipid preprocessing |
| Review pseudobulk layers, upstream QC and gene/lipid matching | Aurélien, with James for original QC | Analysis-ready matrices; avoids treating generated annotations as validated QC |
| Perform KO limma/concordance analysis, then WDSD baseline and diet comparisons | Aurélien | First Figure 3 and Figure 6 drafts |
| Provide Figure 3I proteomics validation and J/K weight source data/artwork | James, using Frank's proteomics pipeline | Complete experimental validation row and supplementary context |
| Review final figure numbering, CF459 inclusion and optional composition/Xenium work | James and Aurélien | Agreed manuscript presentation after initial results |

Other unresolved details: the exact manuscript artwork/legend version was not supplied; the transcript alone cannot fix H–L lettering or all Figure 4 analyses. Formal differential-correlation inference, uncertainty methods, feature-selection sensitivities and multiplicity handling remain implementation choices. Nucleus-level composition/QC inputs may need separate source files beyond the supplied pseudobulks.

**Dependency order:** intake/QC and sample manifests → Cyp2c38 differential analysis and existing-reference concordance; in parallel, WDSD QC and matrix alignment → standard-diet conservation → within-WDSD diet perturbation → CF459 robustness against CF453 WD and the original FCG SD anchor, with a direct matched-pair closeness comparison → joint manuscript figure review. Sex-specific Cyp2c38 predictions require the earlier paired FCG matrices, not WDSD data.

## Manuscript timing and documentation

James described a revision push **between now and February**. No precise analysis or submission deadline was agreed; February 2027 is the natural calendar inference from the recording date, not an explicitly stated year. The current Cell manuscript takes priority over the separate Western-diet RNA manuscript, which had been desk rejected at eLife and was not proceeding soon. Coordinate reuse/publication of those data accordingly. Reviewer comments were not attached, so this document does not invent a mapping to specific comment numbers.

Meeting summaries should be retained in `../Meeting_Summaries`, alongside `../original_datasets`, as requested directly by the user and discussed at 19:13–19:33. The closing London/Heidelberg travel discussion did not set analysis deadlines.

When implementation starts, read the relevant Workflow 03/04/05/07 plan files before changing their scripts, update those plans and cross-workflow documentation if dependencies/outputs change, use explicit tidyverse namespaces, and run each modified R script immediately as required by the repository instructions. The original summary did not change analysis code; dated work-package implementations and reports now document the executed analyses, including the 14 October WP3 extension.

## Evidence map

| Recording interval | Topic |
| --- | --- |
| 00:03–02:44 | CF453 design, age/background differences, existing zonation preparation, RNA/lipid availability and seven male SD mice |
| 02:44–03:58 | Same-edge correlation-vector conservation scatterplots; Mtmr11 example; hypothetical 0.8/0.9 |
| 03:58–05:26 | Within-WDSD SD/WD network comparisons, changing coefficients and possible joint analysis |
| 05:35–07:40 | CF453=`wdsd`, CF459=`wd`, WD-only limitations and evolving robustness-comparison scope |
| 08:04–10:24 | Proposed Figure 6, optional Figure 4/composition/Xenium, upload status and manuscript priority |
| 10:24–13:57 | In vivo Cyp2c38, Figure 3G/H–L, protein validation and heterozygous perturbation |
| 14:11–14:56 | Body/liver weights, young healthy animals and hypothetical six-month WD challenge |
| 15:02–16:24 | Suggested KO-first sequencing, combined/deduplicated lipid input, ignored preliminary statistics and column groups |
| 16:24–18:50 | Pooled reference, pooled/sex-specific KO contrasts, original sex-specific reference recomputation and no KO RNA requirement |
| 19:13–19:33 | Shared Dropbox meeting-summary destination |

Repository checks used the [README](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/README.md), agent guide, cross-workflow overview, [Workflow 03 plan](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/workflow_plans/03_fcg_lipidomic_analysis.md), [Workflow 04 plan](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/workflow_plans/04_mofa_integration.md), [Workflow 05 plan](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/workflow_plans/05_zonation_analysis.md), [Workflow 07 plan](/Users/dugourd/Library/CloudStorage/Dropbox/FCGL_FirstSubmission/gonads_snRNAseq_and_lipidomics_analysis_refactored/docs/workflow_plans/07_validation_experiments.md), their active scripts, the network CSVs and current RNA/lipid input metadata. The supplied new CSV headers and RDS metadata were inspected directly. Updated WDSD metadata and sample maps were rechecked against lipidomics on 6 October 2026; all 28 mouse IDs and group assignments agree. No external mechanistic claims or literature results are needed to support the decisions recorded here.
