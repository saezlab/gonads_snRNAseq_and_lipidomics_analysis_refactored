# Collaborator Source Materials

This folder stores unrefactored collaborator-provided source material retained
for paper-submission provenance. These files are not part of the default
workflow rerun path unless a future task explicitly promotes them into a
normalized workflow.

## 2026-06-19 Missing-Material Bundle

The zip file labeled as MALDI, Xenium, and proteomics material was triaged as
follows:

- Attachment 1, MALDI: `maldi/Cleland_analysis_expt2_recal.ipynb`
  - Negative-mode MALDI-MSI liver zonation notebook.
  - Uses external AnnData, image-registration, Visium/FISH reference, and SCiLS
    inputs that are not committed in this repository.
  - Notebook outputs are local to its configured paths, including a `results/`
    directory.
- Attachments 2-5, Xenium: already represented by
  `scripts/workflows/10_xenium_analysis/00_infer_doublets.py`,
  `01_merge_doublets.py`, `02_compositional_analysis.py`, and
  `03_DGE_analysis.py`.
  - The archive copies were compared with the checked-in scripts and matched
    after ignoring CRLF line-ending differences.
- Attachment 6, proteomics: `proteomics/JamesCleland_P3225_V1.Rmd`
  - Original P3225 fractionated TMT proteomics R Markdown analysis.
  - It uses the collaborator/PCF working directory and local `metadata.csv`
    plus raw protein table files.
  - The active repository proteomics validation in workflow 02 uses the curated
    checked-in files under `data/proteomics/`; it does not source this Rmd.
