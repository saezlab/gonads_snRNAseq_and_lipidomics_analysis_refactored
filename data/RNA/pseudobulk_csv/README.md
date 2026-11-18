# Public pseudobulk CSV inputs

This is the public CSV representation of seven aggregate-profile RNA inputs. There are no individual-nucleus matrices or Seurat objects in this directory. Every CSV file is strictly below 50,000,000 bytes; original serialized objects remain local and ignored by Git.

Each subdirectory has `manifest.json`, a CSV file index, numeric layer parts and profile metadata/identities. Start with the manifest rather than guessing filenames or part order. `counts`, normalized `data` and `scale.data` are supplied source values, not recalculated preprocessing.

- `current_fcg_zonation/`, `current_gdx_zonation/`, `current_wdsd_zonation/`: current expanded-marker inputs.
- `wdsd_hepatocytes/`: whole-hepatocyte diet-study profiles.
- `historical_wdsd_zonation/`: original zone definition for frozen comparisons.
- `wd_cell_types/`, `wd_zonation/`: cell-type and original zone profiles from the separate WD comparison cohort.

Current and historical zone definitions are different datasets. Source-object SHA256 identities remain distinct from hashes of the CSV files.

See [the input guide](../../../docs/PSEUDOBULK_CSV_INPUTS.md), [public file hashes](../../../docs/provenance/current_input_manifest.csv), [round-trip verification](../../../docs/provenance/pseudobulk_csv_roundtrip.json) and [release registry](../../../scripts/config/analysis_release.json).
