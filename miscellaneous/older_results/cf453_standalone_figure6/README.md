# Standalone CF453 Figure 6 panel package

28 wild-type mice (seven per female/male × SD/WD group) contribute 56 paired hepatocyte-zone RNA profiles and one whole-liver lipid profile per mouse. The primary models retain 18,486 genes and 1,129 endogenous lipid features with analytical QC CV ≤30%.

The principal RNA result is diet-specific sex-by-zone association, with a separate test of diet modification:

| contrast | genes_tested | genes_q_005 |
| --- | --- | --- |
| sex_zone_SD | 18,486 | 536 |
| sex_zone_WD | 18,486 | 607 |
| sex_diet_zone | 18,486 | 194 |

Interaction contrasts are male-minus-female differences in periportal-minus-pericentral expression; the three-way contrast is WD minus SD of that interaction. All three share one BH family across genes. Context nominations also require a significant within-sex/diet zonation effect; memberships may overlap.

| diet | distinct_genes | context_memberships | shared_genes |
| --- | --- | --- | --- |
| SD | 529 | 677 | 148 |
| WD | 603 | 783 | 180 |

**No individual gene–lipid correlation passes BH q ≤0.05 among 2,750,244 primary tests (observed count: 0). Association panels are exploratory.**

| diet | tested_pairs | pairs_abs_r_065 | pairs_q_005 |
| --- | --- | --- | --- |
| SD | 1,375,122 | 184,782 | 0 |
| WD | 1,375,122 | 153,727 | 0 |

The original C any-edge rule assigns 100% of eligible lipids in every class and both diets to multiple contexts (shared). This saturated summary provides no class-specific spatial discrimination. It is retained transparently; use the complementary pair-density panels to inspect continuous context differences.

The complementary C_pair_density panels use identical sex/zone candidate-gene unions and lipid pairs in both diets. Across class/context/diet cells, 2.66%–31.05% of eligible pairs exceed |r| ≥ 0.65. These densities are descriptive and do not establish differential correlation or rewiring.

Representative genes are selected for descriptive association stability, without implying a metabolic or regulatory role:

| sex | zone | gene | stable_partner_count_across_diets |
| --- | --- | --- | --- |
| female | PC | St3gal6 | 537 |
| female | PP | Auts2 | 406 |
| male | PC | Unc5c | 575 |
| male | PP | Hspa12a | 456 |

Each representative uses the same selected partners in both diets, including coefficients below the display threshold in one diet. Their ranking and partner selection use the current cohort; examples are not independent validation.

The lipid-quality sensitivity universes contain 1,177 broadly eligible features, 1,129 primary QC-filtered features and 1,058 features after additionally requiring supplied Fill ≥0.5. Their class denominators are recalculated within each universe; these are membership-summary sensitivities rather than independent experiments.

The provenance audit confirms data = ln(1 + 10000 × supplied counts / profile library), and scale.data is a per-gene z-score of data. Maximum normalization replay errors were 4.441e-16 and 8.882e-16 for hepatocyte and zone objects, respectively. This establishes the stored transformation, not the original count-summing provenance.

This is an observational analysis of sex associations in one CF453 cohort; it does not isolate gonadal hormone or sex-chromosome effects. One whole-liver lipid vector is paired with both RNA zones per mouse: RNA-zone associations do not measure spatial lipid concentrations, enzyme activity, substrate specificity or flux. Nucleus contributions, original nuclear QC/gating, count aggregation, batch/litter and assignment provenance remain unverified. Seven mice per sex/diet limit correlation precision. Representative selection uses these same data, and omission/bootstrap checks condition on the fixed RNA and partner selections. Stable signs, large correlations and bootstrap ranges do not provide independent replication or confirm causal regulation.

## Panel files

30 separate panels are exported in PDF, SVG and PNG; [panel_manifest.csv](panel_manifest.csv) records dimensions, formats and matching source tables. The manuscript figure is to be assembled manually.

| panel | family | PDF | source_data |
| --- | --- | --- | --- |
| Fig6A_SD | A | [PDF](panels/Fig6A_SD.pdf) | [CSV](source_data/Fig6A_SD.csv) |
| Fig6A_WD | A | [PDF](panels/Fig6A_WD.pdf) | [CSV](source_data/Fig6A_WD.csv) |
| Fig6B_design | B | [PDF](panels/Fig6B_design.pdf) | [CSV](source_data/Fig6B_design.csv) |
| Fig6C_SD | C | [PDF](panels/Fig6C_SD.pdf) | [CSV](source_data/Fig6C_SD.csv) |
| Fig6C_pair_density_SD | C_pair_density | [PDF](panels/Fig6C_pair_density_SD.pdf) | [CSV](source_data/Fig6C_pair_density_SD.csv) |
| Fig6C_WD | C | [PDF](panels/Fig6C_WD.pdf) | [CSV](source_data/Fig6C_WD.csv) |
| Fig6C_pair_density_WD | C_pair_density | [PDF](panels/Fig6C_pair_density_WD.pdf) | [CSV](source_data/Fig6C_pair_density_WD.csv) |
| Fig6D_female_pericentral_SD | D | [PDF](panels/Fig6D_female_pericentral_SD.pdf) | [CSV](source_data/Fig6D_female_pericentral_SD.csv) |
| Fig6D_female_pericentral_WD | D | [PDF](panels/Fig6D_female_pericentral_WD.pdf) | [CSV](source_data/Fig6D_female_pericentral_WD.csv) |
| Fig6D_female_periportal_SD | D | [PDF](panels/Fig6D_female_periportal_SD.pdf) | [CSV](source_data/Fig6D_female_periportal_SD.csv) |
| Fig6D_female_periportal_WD | D | [PDF](panels/Fig6D_female_periportal_WD.pdf) | [CSV](source_data/Fig6D_female_periportal_WD.csv) |
| Fig6D_male_pericentral_SD | D | [PDF](panels/Fig6D_male_pericentral_SD.pdf) | [CSV](source_data/Fig6D_male_pericentral_SD.csv) |
| Fig6D_male_pericentral_WD | D | [PDF](panels/Fig6D_male_pericentral_WD.pdf) | [CSV](source_data/Fig6D_male_pericentral_WD.csv) |
| Fig6D_male_periportal_SD | D | [PDF](panels/Fig6D_male_periportal_SD.pdf) | [CSV](source_data/Fig6D_male_periportal_SD.csv) |
| Fig6D_male_periportal_WD | D | [PDF](panels/Fig6D_male_periportal_WD.pdf) | [CSV](source_data/Fig6D_male_periportal_WD.csv) |
| Fig6E_diet_zonation | E | [PDF](panels/Fig6E_diet_zonation.pdf) | [CSV](source_data/Fig6E_diet_zonation.csv) |
| Fig6F_female_pericentral | F | [PDF](panels/Fig6F_female_pericentral.pdf) | [CSV](source_data/Fig6F_female_pericentral.csv) |
| Fig6F_female_periportal | F | [PDF](panels/Fig6F_female_periportal.pdf) | [CSV](source_data/Fig6F_female_periportal.csv) |
| Fig6F_male_pericentral | F | [PDF](panels/Fig6F_male_pericentral.pdf) | [CSV](source_data/Fig6F_male_pericentral.csv) |
| Fig6F_male_periportal | F | [PDF](panels/Fig6F_male_periportal.pdf) | [CSV](source_data/Fig6F_male_periportal.csv) |
| FigS6_PCA_RNA | supplementary_QC | [PDF](panels/FigS6_PCA_RNA.pdf) | [CSV](source_data/FigS6_PCA_RNA.csv) |
| FigS6_PCA_Lipids | supplementary_QC | [PDF](panels/FigS6_PCA_Lipids.pdf) | [CSV](source_data/FigS6_PCA_Lipids.csv) |
| FigS6_lipid_sex_SD | supplementary_lipid | [PDF](panels/FigS6_lipid_sex_SD.pdf) | [CSV](source_data/FigS6_lipid_sex_SD.csv) |
| FigS6_lipid_sex_WD | supplementary_lipid | [PDF](panels/FigS6_lipid_sex_WD.pdf) | [CSV](source_data/FigS6_lipid_sex_WD.csv) |
| FigS6_lipid_diet_female | supplementary_lipid | [PDF](panels/FigS6_lipid_diet_female.pdf) | [CSV](source_data/FigS6_lipid_diet_female.csv) |
| FigS6_lipid_diet_male | supplementary_lipid | [PDF](panels/FigS6_lipid_diet_male.pdf) | [CSV](source_data/FigS6_lipid_diet_male.csv) |
| FigS6_lipid_sex_diet | supplementary_lipid | [PDF](panels/FigS6_lipid_sex_diet.pdf) | [CSV](source_data/FigS6_lipid_sex_diet.csv) |
| FigS6_display_threshold | supplementary_threshold | [PDF](panels/FigS6_display_threshold.pdf) | [CSV](source_data/FigS6_display_threshold.csv) |
| FigS6_lipid_quality | supplementary_quality | [PDF](panels/FigS6_lipid_quality.pdf) | [CSV](source_data/FigS6_lipid_quality.csv) |
| FigS6_zonation_markers | supplementary_QC | [PDF](panels/FigS6_zonation_markers.pdf) | [CSV](source_data/FigS6_zonation_markers.csv) |

## Reproducible records

[Methods](METHODS.md) and [panel legends](PANEL_LEGENDS.md) describe the tests and visual encodings. `tables/` contains complete gene tests, nomination/selection ledgers, lipid tests, all correlation coefficients, memberships and diet changes. `source_data/` matches individual panels. `analysis.rds` preserves the full fitted analysis; `bootstrap_indices.rds` preserves animal draw banks. `QC/`, `QA/provenance/`, independent replay checks under `QA/`, `provenance_sha256.csv` and `sessionInfo.txt` record input QC, numerical verification and provenance.

Paired bootstrap delta-r summaries for 24 selected gene–lipid partners are in [representative_delta_summary.csv](tables/representative_delta_summary.csv), with requested/valid draws and percentile ranges. These are descriptive fixed-selection quantities, not tests of differential correlation.
