# WP3 frozen CF459 inputs

Imported on **14 October 2026** for the male Western-diet cross-cohort robustness
comparison. These are byte-identical copies of the shared archive; no source was
removed or edited. `wp3_input_migration_manifest.csv` records original and local
paths, byte sizes and independently verified SHA256 hashes. Three files total
**13,118,772 bytes**.

| Repository path, relative to `data/` | Role |
|---|---|
| `metabolomics/wd/lipidomics_wd.csv` | CF459 WD-only lipidomics; 1,463 supplied feature rows and eight biological mice |
| `RNA/wd/wd_pseudobulks.rds` | Eight mice × five cell types; primary whole-hepatocyte matrix is the eight hepatocyte columns |
| `RNA/wd/wd_zonation_pseudobulks.rds` | Eight mice × two zones; primary comparison uses four male mice per zone |

The CF453 comparator uses the existing local files documented in
[WP1_WP2_INPUTS.md](WP1_WP2_INPUTS.md); there is no second import or alteration
of those inputs. WP3 reuses frozen Figure 5 feature selections and strict
formatting-only molecular matching. Historical WP1/WP2 results are preserved.

## Mouse map and scope

| CF459 IDs | Supplied RNA group | Recorded sex | Analysis diet | Role |
|---|---|---|---|---|
| `cf459-1`, `cf459-2`, `cf459-3`, `cf459-4` | `f-ctrl` | female | WD | Intake/sex-map QC and original scaling-roster diagnostics only |
| `cf459-13`, `cf459-14`, `cf459-15`, `cf459-16` | `m-ctrl` | male | WD | Primary cross-cohort comparison with seven CF453 male WD mice |

The 5 October meeting establishes CF459 as **Western diet only**. The word
`ctrl` in its RNA labels does not identify standard diet. Lipid headers identify
mice without independently encoding sex; the recorded sex map comes from RNA
metadata, validated against the fixed roster above in both Seurat objects.
All eight IDs match both modalities. Each mouse has one profile for every cell
type and exactly one Periportal and one Pericentral profile. Matrices are ordered
by explicit mouse IDs rather than by the differing lipid/RNA column orders.
There is no within-CF459 SD–WD comparison and no female network analysis in WP3.

## Lipid import and identity policy

The original CSV has a **UTF-8 BOM, semicolon delimiter, decimal commas,
scientific notation and quoted semicolons within annotations**. Read it using
`read.csv2(..., fileEncoding="UTF-8-BOM", check.names=FALSE)`; the supplied bytes
are kept intact. All 11,704 intensity cells are finite, positive numeric values
after correct parsing. The primary matrix is `log2` of supplied intensities;
this does not recover or independently validate upstream normalization or
imputation.

The intake retains all 1,463 annotation rows in an exclusion/multiplicity
ledger with unique `WD_0001`-style feature IDs. It reuses
`wp1_canonical_lipid()`: select the final pipe-delimited annotation and format
`CLASS structure` as `CLASS(structure)`. Chain order, slash/underscore, ether
prefixes and oxygenation are preserved. There are 1,103 pipe-delimited rows,
zero original/canonical-name duplicates and **244 duplicate Alignment-ID excess
rows**. Alignment IDs are not cross-cohort join keys. Molecular identity remains
strict and one-to-one; formatting agreement is not independent chemical proof.

Applying the existing class/standard/isotope policy retains **1,395 lipids**:
68 excluded class rows (CAR 20, DGGA 10, HBMP 12, NATau 1, PEtOH 4, PMeOH 9,
SL 8, SPB 3, VAE 1), with no internal-standard or labelled-isotope exclusions.
The exclusion-class vector is the same as in WP2, including classes absent
from this export. No feature or animal is excluded because of its network
association or apparent biological response.

**No analytical-QC injections, adduct, RT/m/z or Fill percentage are supplied.**
`qc_cv_percent` and `keep_qc30` are `NA`, with explicit availability flags false;
there is no invented QC-CV estimate, Fill sensitivity or substitute pass filter.

## Stored RNA layers and diagnostics

Both objects record Seurat **5.3.0**, one RNA Assay5 and the named layers
`counts`, `data`, `scale.data`. Access uses `SeuratObject::LayerData()` to retain
the Assay5 gene/profile names. Every layer has **20,123 genes**, and respectively
40 cell-type or 16 zonation profiles. All entries are finite; `counts` is sparse,
integer and nonnegative, `data` is dense/nonnegative and `scale.data` is
signed/dense. The primary matrices use the stored `scale.data` unchanged;
there is no re-normalization or re-scaling for the four male mice.

Row scaling the existing `data` across the **complete original object roster**
reproduces stored `scale.data` to maximum absolute differences of
9.77 × 10^-15 (40 cell-type profiles) and 1.24 × 10^-14 (16 zone profiles).
The cell-type diagnostic is performed before hepatocyte subsetting and the
zone diagnostic before male subsetting. This is a layer-fidelity diagnostic,
not independent normalization validation. Only `ScaleData.RNA` appears in
recorded command histories.

Aggregate supplied count-library ranges are 59,105–10,670,693 for cell-type
profiles and 305,248–7,192,075 for zonation profiles; detected-gene ranges are
8,576–16,810 and 11,567–16,323 respectively. These values are **not nucleus
counts**. Contributing nuclei, doublets, original QC thresholds, composition,
missing-value treatment and zonation-gating provenance remain unavailable.

All 16 aggregate zone profiles have the expected Cyp2e1/Cyp2f2 ordering and pass
the conventional aggregate expression threshold diagnostic. This cannot
validate original nucleus-level gates. Cyp2e1, Cyp2f2, Glul, Gls2, Hal, Sds,
Mtmr11 and Cyp2c38 are exported as per-profile/sex/zone diagnostics. The raw exports do not record exact ages, genetic background, diet exposure/
formulation, litter or batch details. On **14 October 2026**, the user clarified
cohort ages: **CF459 3 months; original FCG 3 months; CF453 10 months**. These
are supplied cohort descriptions; no age was inferred from network agreement.
CF453 was described in the 5 October meeting as approximately 10 months after
9 months of diet exposure. CF459 exposure duration and the other metadata
remain unestablished. The amended WP3 compares young CF459 WD against both
young FCG SD (22 testes-bearing FCG mice, 11 XXT/11 XYT) and older CF453 WD
(7 males); the new FCG anchor does not introduce a CF459 SD arm.

## Intake API and verification

From the repository root:

```r
source("scripts/support_functions.R")
source("scripts/core/wd_intake.R")
wd <- wd_load_inputs(root = ".")
```

`wd` returns aligned `metadata`, `lipid`, `lipid_supplied`, `lipid_ledger`,
zone-specific `rna`, `rna_data`, `rna_data_scaled`, `rna_counts`, whole-hepatocyte
matrices, profile/layer/marker diagnostics, processing limitations and local
`input_paths`. All eight mice are loaded for audit; statistical code explicitly
selects `sex == "male"` and `diet == "WD"` for the primary comparison.

Run `Rscript tests/test_wd_intake.R`. Checks include frozen imported hashes and
sizes, correct locale parsing, exclusions/unavailable-QC policy, original named
Seurat-layer replay, deliberately permuted matrix/metadata orders, and explicit
failures for duplicate/missing mice, changed sex labels, inconsistent profile
names and repeated or absent mouse-zone profiles.
