# Liver zonation and gene-lipid associations under standard and Western diets

The CF453 cohort contains 28 approximately 10-month-old wild-type mice: seven females and seven males on each of standard diet (SD) and Western diet (WD). Each mouse contributes paired periportal/pericentral liver RNA and one whole-liver lipid measurement. Three panel groups connect selected zonation genes, relative lipid-class associations and exploratory correlation changes.

## RNA anchors and direct diet evidence

40 anchor selections represent 35 distinct genes. Selection uses anchor-diet RNA specificity and zonation magnitude; neither the SD–WD shift nor lipid correlations ranks the anchors.

| selection_diet | arrows | direct_threeway_q_005 |
| --- | --- | --- |
| SD | 20 | 4 |
| WD | 20 | 6 |

Arrows track the same genes across diets even if their destination classification changes. Solid arrows have direct three-way diet×sex×zone q≤0.05; other arrows remain visible and dashed. Changes in t-statistic coordinates combine effect and precision, and anchor selection can produce regression toward weaker destination estimates.

## Fixed selected programmes and class changes

| context | fixed_union_genes |
| --- | --- |
| female_PC | 7 |
| female_PP | 9 |
| male_PC | 10 |
| male_PP | 9 |

Both diets use the same per-context union of anchor genes and 1,129 eligible lipid features. Across the class/context cells, ΔE ranges from -27.13 to 72.96 percentage points. Defined bootstrap draws range from 991 to 993 of 999 per context; intervals and directional fractions condition on those defined draws.

0/120 conditional 95% percentile ranges exclude zero; 2 cells meet the descriptive ≥0.90 directional marker. The largest absolute displayed change is LPG in female PC (1 measured lipid). Small class denominators can generate extreme rates; marker stability is not class-level significance.

For class C, E is its fraction of gene–lipid pairs with |r|≥0.65 minus the corresponding fraction across all other lipid classes. Panel B shows ΔE=E_WD−E_SD. Positive ΔE means increased relative association-display density, not increased lipid abundance or spatial localization.

The explicitly selected pair/context/diet family contains 79,030 correlation tests; 0 pass its BH q≤0.05 threshold. This is conditional on the selected RNA programmes and is not a test of differential correlation.

## Selected gain/loss examples

| context | gene | n_edges | status |
| --- | --- | --- | --- |
| female_PP | Fmo3 | 6 | Exploratory selected example |
| female_PC | Cyp3a41a | 6 | Exploratory selected example |
| male_PP | Col27a1 | 6 | Exploratory selected example |
| male_PC | Osbpl3 | 6 | Exploratory selected example |

Network examples satisfy explicit strength, omission and bootstrap-stability gates. Gain/loss describes a threshold-defined change in observed |r|. The selected illustrations do not establish gene regulation or a significant rewiring event.

## Panel files

| panel | PDF | source |
| --- | --- | --- |
| Fig6A_anchor_gene_shifts | [PDF](panels/Fig6A_anchor_gene_shifts.pdf) | [CSV](source_data/Fig6A_anchor_gene_shifts.csv) |
| Fig6B_lipid_class_shift | [PDF](panels/Fig6B_lipid_class_shift.pdf) | [CSV](source_data/Fig6B_lipid_class_shift.csv) |
| Fig6C_network_gain_loss | [PDF](panels/Fig6C_network_gain_loss.pdf) | [CSV](source_data/Fig6C_network_gain_loss.csv) |

PDF/SVG/PNG versions and matching CSVs are listed in `panel_manifest.csv`. No extra or assembled figure is produced. [Methods](METHODS.md) and [legends](PANEL_LEGENDS.md) define all selection and uncertainty rules. Tables retain exclusions, exact denominators, correlation evidence, bootstrap valid counts and observed sensitivity summaries.

This is a standalone observational analysis of approximately 10-month-old CF453 wild-type mice, seven per sex/diet. The specificity gate and correlation thresholds define displayed programmes, not biological exclusivity, causal regulation or established rewiring. Whole-liver lipids are not spatial lipid measurements. Nucleus contributions, raw count aggregation, original nuclear QC/gating, batch/litter and assignment provenance remain unverified. Selection and bootstrap summaries use these same animals; fixed-selection stability is not independent validation or a calibrated post-selection significance test.
