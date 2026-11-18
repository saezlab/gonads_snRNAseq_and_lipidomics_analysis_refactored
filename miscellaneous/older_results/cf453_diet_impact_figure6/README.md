# CF453: direct diet-impact panels

28 mice (seven per female/male × SD/WD group), 18,486 modeled genes and 1,129 analytical-QC-filtered whole-liver lipid features. Four separate panels focus on direct diet interactions.

The fitted model compares WD with SD and tests whether that difference depends on sex, zone, or both. This tests diet impact directly; comparing separate SD and WD significance labels would not.

## Results

| test | genes_tested | q_005 |
| --- | --- | --- |
| Diet × sex, averaged over zones | 18,486 | 166 |
| Diet × zone, averaged over sexes | 18,486 | 350 |
| Diet × sex × zone | 18,486 | 137 |

The separate three-degree-of-freedom RNA omnibus test identifies 605 genes at BH q≤0.05. Its ranking supplies the 20 genes decomposed in panel B.

Diet × sex is supported for 267/1,129 lipids (23.65%); 228 interactions are positive and 39 negative. Positive means the male WD−SD response is more positive, or less negative, than the female response.

Panel C examples and their direct interaction estimates (natural-log normalized-expression units; intervals are nominal 95%):

| gene | contrast | effect | nominal_95_interval | primary_q |
| --- | --- | --- | --- | --- |
| Egfr | diet_sex | -2.107 | [-3.099, -1.115] | 0.0137 |
| Limch1 | diet_zone | 0.802 | [0.706, 0.898] | 9.32e-18 |
| Meg3 | diet_sex_zone | 0.878 | [0.594, 1.163] | 6.52e-05 |
| Cyp3a11 | diet_sex_zone | -0.733 | [-1.064, -0.403] | 0.00837 |

## Direct comparisons

Let Δ(s,z)=mean(WD,s,z)−mean(SD,s,z), with z=PP or PC. The direct contrasts are:

- Diet × sex: average over zones of Δ(male,z)−Δ(female,z).
- Diet × zone: average over sexes of Δ(s,PP)−Δ(s,PC).
- Three-way: [Δ(male,PP)−Δ(male,PC)]−[Δ(female,PP)−Δ(female,PC)].

Positive values follow these signed differences; they need not imply a positive diet effect in every group. The interaction is the test of a difference in response; significance in only one diet or sex is not such a test.

## Files and interpretation

| panel | PDF | source |
| --- | --- | --- |
| Fig6A_direct_RNA_interactions | [PDF](panels/Fig6A_direct_RNA_interactions.pdf) | [CSV](source_data/Fig6A_direct_RNA_interactions.csv) |
| Fig6B_conditional_diet_effects | [PDF](panels/Fig6B_conditional_diet_effects.pdf) | [CSV](source_data/Fig6B_conditional_diet_effects.csv) |
| Fig6C_mouse_examples | [PDF](panels/Fig6C_mouse_examples.pdf) | [CSV](source_data/Fig6C_mouse_examples.csv) |
| Fig6D_lipid_diet_effects | [PDF](panels/Fig6D_lipid_diet_effects.pdf) | [CSV](source_data/Fig6D_lipid_diet_effects.csv) |

PDF, SVG and PNG versions and matching source CSVs are listed in `panel_manifest.csv`. Assemble the final figure manually. [Methods](METHODS.md) and [legends](PANEL_LEGENDS.md) define the models, contrasts, uncertainty and selection.

These results describe diet-associated differences within CF453, with seven mice per sex/diet and paired RNA zones. They do not isolate hormonal or chromosome effects, and whole-liver lipids are not spatial lipid measurements. Source count aggregation, nucleus contributions/QC, batch/litter and assignment provenance remain unverified. Displayed genes are selected from the same data; nominal intervals do not account for that selection.
