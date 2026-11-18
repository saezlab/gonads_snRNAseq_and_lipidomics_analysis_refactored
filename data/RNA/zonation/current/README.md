# Local source objects and compatibility CSVs

The public pseudobulk input format is now the seven verified CSV exports under `data/RNA/pseudobulk_csv/`. Every RNA CSV file must be strictly below 50,000,000 bytes. See `docs/PSEUDOBULK_CSV_INPUTS.md` and the release registry.

The three original expanded-marker Seurat RDS files remain here locally for provenance and optional re-export. They are ignored by Git and are not required by CSV-backed scientific defaults. Their identities are retained in `docs/provenance/current_input_manifest_20261116.csv` and the export manifests; the approved production objects are unchanged.

The small `fcg/` and `gdx/` wide scaled-expression CSVs/metadata remain compatibility inputs for the validated FCG engine. The three supplied preprocessing scripts under `james_scripts_20261104/` and candidate `Zonation_Genes.csv` remain public provenance. The candidate marker panel does not establish the exact final retained markers in each cohort.
