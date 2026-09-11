# Factor 1 MOON input preparation

These inputs use the saved Factor 1 in
`results/mofa/main_model_RNA_lipid_weight.RData`, verified against
`results/mofa/RNA/mofa_res_12factor.hdf5`. The model started with 12 factors and
retained 10. The saved Factor 1 loading column agrees with native HDF5 factor 1
to within 5.2e-15 for both hepatocyte RNA and lipids.

## Reproduce

Run from the repository root, in this order:

```sh
Rscript scripts/moon/prepare_factor1_inputs.R
python3 scripts/moon/audit_factor1_model.py
Rscript scripts/moon/map_factor1_inputs.R
```

The first step requires the already-installed R package `decoupleR`; the final
mapping step uses base R. The Python audit requires `numpy` and `h5py`. This run
used Python 3.12.13, numpy 2.5.2, and h5py 3.16.0. With `uv` available, an optional
dependency-managed equivalent for the middle command is:

```sh
uv run --python 3.12 --with numpy==2.5.2 --with h5py==3.16.0 python scripts/moon/audit_factor1_model.py
```

The exact installed-cache invocation used in this session was:

```sh
PYTHONPATH=/Users/dugourd/.cache/uv/archive-v0/3BxlphLcUkFq3Zbj:/Users/dugourd/.cache/uv/archive-v0/9zJ8qjSOH6odT0xd /Users/dugourd/.local/share/uv/python/cpython-3.12.13-macos-aarch64-none/bin/python3.12 scripts/moon/audit_factor1_model.py
```

Those cache locations are machine-specific and may be removed by cache cleanup;
the dependency-managed command does not rely on them. No model training is
performed by any preparation step. After preparation, the separate MOON runner
is `Rscript scripts/moon/run_factor1_moon.R`.

## Data and interpretation

`factor1_inputs.rds` contains:

- `tf_activity`: 367 named hepatocyte TF footprint scores, with `Gene__` IDs;
  353 are present in the selected PKN before reachability pruning.
- `lipid_weight`: 155 named lipid loadings mapped to PKN nodes: 149 explicit
  measured-species rule nodes and six unambiguous exact-resolution ChEBI nodes.
- `hepatocyte_gene_weight`: all 20,123 hepatocyte RNA loadings, with `Gene__` IDs.
  Nonzero loading is not evidence of detectable expression.
- `tf_regulon`: 37,262 signed target-regulation rows (`source`, `target`, `mor`),
  with the PKN identifier prefix. These reproduce the legacy TF score inference:
  CollecTRI with the original Ar/Esr1 regulons replaced by the repository's
  hormone-target table, and a minimum of 20 measured targets.
- `mapping`: all 1,254 lipid observations, including original and oriented
  weights, mapping type, resolution, exclusions, and compartment uncertainty.
- `orientation`: the joint Factor 1 sign convention; multiplier is +1.
- `tf_statistics`: full ULM scores, raw p-values, BH-adjusted p-values, measured
  target counts, saved-score comparison, and PKN presence.
- `original_tf_activity`, `original_lipid_weight`, and
  `original_hepatocyte_gene_weight`: preserved values before mapping/orientation.

TFs are inferred from their target-gene loading patterns, not from TF mRNA.
Recomputed ULM scores reproduce the saved scores with maximum error 2.26e-12.
The TF view's factor columns are lexicographically ordered, so all extraction is
by the name `Factor1`, never column position. These ULM p-values describe the
footprint calculation; they do not quantify uncertainty in the fitted MOFA
model, validate causal mechanisms, or account for later network selection.

All 1,254 lipid loadings are retained in `factor1_inputs_raw.rds` and
`factor1_lipid_weights_all.csv`. The prior network's 64 selected lipid features
are annotated but do not restrict this pilot. The mapped set has seven of those
64 features. Class ancestry, chemically compatible but more specific structures,
and protonation substitutions are not used as measurement identities. One lipid
has multiple exact ChEBI alternatives and is excluded; 1,098 have no accepted
exact node in the selected PKN. No two accepted features map to the same node,
and no accepted feature needs compartment copying. All accepted compartment
labels are `unknown`, consistent with the lack of measured localization.

Positive Factor 1 is already associated with testes: the sample-coordinate means
are -0.742291 for ovaries and +0.708550 for testes. No sign reversal was applied.
This common sign is used for RNA loadings, TF scores, and lipid loadings. It is
an axis orientation, not a direct testes-minus-ovaries differential effect.

Factor 1 explains 8.906% of hepatocyte RNA variance and 29.958% of lipid variance.
Its saved association with gonadal type is R² = 0.969836. These are distinct
quantities: the latter does not mean it explains 97% of total omics variance.
Other RNA views are preserved in `factor1_all_views_raw.csv`; the pilot uses
hepatocyte TFs and bulk-liver lipid measurements.

## Audit files

- `factor1_model_audit.json`, `factor1_model_alignment.csv`: fitted model
  dimensions, training settings, loading identity, and joint orientation.
- `factor1_sample_coordinates.csv`: all 43 matched samples and group labels.
- `factor1_variance_explained.csv`: per-view Factor 1 and total model variance.
- `factor1_existing_model_consistency.csv`: alignment to the existing fits
  initialized with 5–15 factors. Factor-coordinate correlations are 0.977–0.999
  for the alternative fits. This is descriptive consistency across factor-number
  choices, not repeated-seed, bootstrap, or donor-resampling stability evidence.
- `saved_factor_gonadal_enrichment.csv`: unchanged saved regression statistics.
- `input_source_manifest.csv`, `preparation_session_info.txt`: source hashes
  and the R software environment.

The 29,055 directed regulator–target pairs shared by the legacy TF regulon and
the new PKN have matching signs; no opposite-sign overlaps were found. Using the
same RNA loadings for footprint inference and TF–target coherence is an internal
consistency operation, not independent validation of the resulting mechanisms.
