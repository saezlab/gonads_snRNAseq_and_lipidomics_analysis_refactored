#' Dataset Configurations
#'
#' This file defines all dataset configurations used in the FCG liver analysis project.
#' Each dataset is defined with its RNA, lipidomic, and metadata file paths.
#'
#' This centralizes dataset definitions from scripts like:
#' - seurat_to_df*.R
#' - pseudobulk_DA*.R
#' - lipidomic_DA*.R
#' - gdxhr_*.R

#' @description
#' DATASETS is a named list where each element defines a dataset configuration.
#' Each dataset configuration can include:
#' - rna: List with RNA data configuration (files, metadata, split info)
#' - lipid: List with lipidomic data configuration (file, metadata info)
#' - description: Human-readable description
#' - notes: Additional notes about the dataset
DATASETS <- list(

  # Four Core Genotypes (FCG) - Main Dataset
  # =========================================
  # This is the primary dataset for the study, containing ~466,986 high-quality nuclei
  # from 43 FCG livers + 5 wild-type males.
  # Samples are pseudobulked by cell type and mouse ID.
  fcg_main = list(
    rna = list(
      # RNA data is split into 3 files due to GitHub size constraints
      files = c("data/RNA/scaled_RNA_1.csv",
                "data/RNA/scaled_RNA_2.csv",
                "data/RNA/scaled_RNA_3.csv"),
      metadata = "data/RNA/meta_data_df.csv",
      split = TRUE,
      n_splits = 3
    ),
    lipid = list(
      file = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    overview = "data/Suppl_Table_1_OverviewData.xlsx",
    description = "Four Core Genotypes main dataset (43 FCG + 5 WT)",
    notes = "Primary dataset for gonadal hormone vs. sex chromosome analysis"
  ),

  # Four Core Genotypes with Wild-Type Controls
  # ============================================
  # Includes XYT wild-type controls to compare with FCG XYT
  fcg_with_wt = list(
    rna = list(
      source_rds = "data/RNA/fcgl_pseudobulks_gonad_withxytwt.rds",
      # This would need to be converted via seurat_to_df_withWT.R
      files = c("data/RNA/scaled_RNA_withWT_1.csv",
                "data/RNA/scaled_RNA_withWT_2.csv",
                "data/RNA/scaled_RNA_withWT_3.csv"),
      metadata = "data/RNA/meta_data_df_withWT.csv",
      split = TRUE,
      n_splits = 3
    ),
    lipid = list(
      file = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    overview = "data/Suppl_Table_1_OverviewData.xlsx",
    description = "FCG dataset including XYT-WT controls",
    notes = "Used to validate FCG model against wild-type XYT males"
  ),

  # Gonadectomy + Hormone Replacement (GDXHR) Dataset
  # ==================================================
  # Validation experiment: gonadectomy followed by hormone replacement
  # Used to confirm activational (rather than organizational) hormone effects
  gdxhr = list(
    rna = list(
      files = c("data/RNA/scaled_RNA_gdxhr_1.csv",
                "data/RNA/scaled_RNA_gdxhr_2.csv",
                "data/RNA/scaled_RNA_gdxhr_3.csv"),
      metadata = "data/RNA/meta_data_df_gdxhr.csv",
      split = TRUE,
      n_splits = 3,
      source_rds = "data/RNA/gdxhr_pseudobulks.rds"
    ),
    lipid = list(
      file = "data/metabolomics/lipidomics_gdxhr.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    overview = "data/Suppl_Table_1_OverviewData.xlsx",
    description = "Gonadectomy + hormone replacement validation dataset",
    notes = "Correlation with FCG: r = -0.48 (expected due to opposite directionality)"
  ),

  # Zonation-Specific Datasets
  # ===========================
  # Hepatocytes subdivided by zonation markers (periportal vs. pericentral)
  # Periportal: Cyp2f2+, Gls2+, Hal+, Sds+
  # Pericentral: Cyp2e1+, Glul+

  fcg_zonation = list(
    rna = list(
      files = c("data/RNA/scaled_RNA_zonation_1.csv",
                "data/RNA/scaled_RNA_zonation_2.csv",
                "data/RNA/scaled_RNA_zonation_3.csv"),
      metadata = "data/RNA/meta_data_df_zonation.csv",
      split = TRUE,
      n_splits = 3,
      source_rds = "data/RNA/fcg_zonation_pseudobulks.rds"
    ),
    lipid = list(
      file = "data/metabolomics/Lipid_data_file_JPC_cleaned.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    overview = "data/Suppl_Table_1_OverviewData.xlsx",
    description = "FCG hepatocytes subdivided by spatial zonation",
    notes = "Analyzes gonad-dependent effects in periportal vs. pericentral hepatocytes"
  ),

  gdxhr_zonation = list(
    rna = list(
      files = c("data/RNA/scaled_RNA_gdxhr_zonation_1.csv",
                "data/RNA/scaled_RNA_gdxhr_zonation_2.csv",
                "data/RNA/scaled_RNA_gdxhr_zonation_3.csv"),
      metadata = "data/RNA/meta_data_df_gdxhr_zonation.csv",
      split = TRUE,
      n_splits = 3,
      source_rds = "data/RNA/gdxhr_heps_zonationsubclustered_pseudobulk_group_20241028.rds"
    ),
    lipid = list(
      file = "data/metabolomics/lipidomics_gdxhr.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    overview = "data/Suppl_Table_1_OverviewData.xlsx",
    description = "GDXHR hepatocytes subdivided by spatial zonation",
    notes = "Validates zonation-specific gonad effects in GDXHR model"
  ),

  # Validation Experiments (siRNA Knockdowns)
  # ==========================================
  # Lipidomic data from siRNA knockdown experiments in primary hepatocytes
  # Used to validate gene-lipid associations identified in MOFA

  validation_expt1 = list(
    lipid = list(
      file = "data/metabolomics/lipidomics_gonad_validation.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    description = "Validation experiment 1: siRNA knockdowns",
    notes = "First batch of siRNA validation for Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9"
  ),

  validation_expt2 = list(
    lipid = list(
      file = "data/metabolomics/lipidomics_expt2.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    description = "Validation experiment 2: Additional siRNA knockdowns",
    notes = "Second batch including Sult3a1, Mtmr11, Abcb1a"
  ),

  validation_20250319 = list(
    lipid = list(
      file = "data/metabolomics/lipidomics_20250319transfections.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    description = "Transfection experiment: March 19, 2025",
    notes = "Zonation-specific gene validation"
  ),

  validation_20250325 = list(
    lipid = list(
      file = "data/metabolomics/lipidomics_20250325transfections.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    description = "Transfection experiment: March 25, 2025",
    notes = "Follow-up validation experiments"
  ),

  validation_20250423 = list(
    lipid = list(
      file = "data/metabolomics/lipidomics_20250423transfections.csv",
      log_transform = TRUE,
      clean_names = TRUE
    ),
    description = "Transfection experiment: April 23, 2025",
    notes = "Final validation round"
  ),

  # PNPLA3 Knockout Dataset
  # =======================
  # Side project examining PNPLA3 knockout effects
  pnpla3ko = list(
    rna = list(
      source_rds = "data/RNA/pnpla3ko_sci_pseudobulks.rds",
      files = c("data/RNA/scaled_RNA_pnpla3ko_1.csv",
                "data/RNA/scaled_RNA_pnpla3ko_2.csv",
                "data/RNA/scaled_RNA_pnpla3ko_3.csv"),
      metadata = "data/RNA/meta_data_df_pnpla3ko.csv",
      split = TRUE,
      n_splits = 3
    ),
    description = "PNPLA3 knockout comparison dataset",
    notes = "Side project not included in main manuscript"
  ),

  pnpla3ko_zonation = list(
    rna = list(
      source_rds = "data/RNA/pnpla3ko_sci_heps_zonationsubsetted_pseudobulks.rds",
      files = c("data/RNA/scaled_RNA_pnpla3ko_zonation_1.csv",
                "data/RNA/scaled_RNA_pnpla3ko_zonation_2.csv",
                "data/RNA/scaled_RNA_pnpla3ko_zonation_3.csv"),
      metadata = "data/RNA/meta_data_df_pnpla3ko_zonation.csv",
      split = TRUE,
      n_splits = 3
    ),
    description = "PNPLA3 knockout with zonation analysis",
    notes = "Zonation-specific analysis for PNPLA3 KO"
  )
)


#' Get Dataset Configuration
#'
#' Helper function to retrieve a specific dataset configuration.
#'
#' @param dataset_name Character. Name of the dataset (e.g., "fcg_main")
#'
#' @return List containing the dataset configuration
#'
#' @examples
#' fcg_config <- get_dataset_config("fcg_main")
#' print(fcg_config$description)
get_dataset_config <- function(dataset_name) {
  if (!dataset_name %in% names(DATASETS)) {
    stop("Dataset '", dataset_name, "' not found. Available datasets: ",
         paste(names(DATASETS), collapse = ", "))
  }
  return(DATASETS[[dataset_name]])
}


#' List Available Datasets
#'
#' Print a summary of all available datasets with descriptions.
#'
#' @examples
#' list_datasets()
list_datasets <- function() {
  cat("Available Datasets:\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  for (dataset_name in names(DATASETS)) {
    dataset <- DATASETS[[dataset_name]]
    cat("Name:       ", dataset_name, "\n")
    cat("Description:", dataset$description, "\n")
    if (!is.null(dataset$notes)) {
      cat("Notes:      ", dataset$notes, "\n")
    }
    if (!is.null(dataset$rna)) {
      cat("RNA:         Available\n")
    }
    if (!is.null(dataset$lipid)) {
      cat("Lipid:       Available\n")
    }
    cat("\n")
  }
}


#' Validate Dataset Configuration
#'
#' Checks if all files specified in a dataset configuration exist.
#'
#' @param dataset_name Character. Name of the dataset to validate.
#' @param check_split_only Logical. If TRUE, only check if split files can be created
#'                         from source RDS (don't require them to exist). Default: FALSE
#'
#' @return Logical. TRUE if all files exist, FALSE otherwise (with warnings).
#'
#' @examples
#' validate_dataset("fcg_main")
validate_dataset <- function(dataset_name, check_split_only = FALSE) {
  dataset <- get_dataset_config(dataset_name)
  all_exist <- TRUE

  cat("Validating dataset:", dataset_name, "\n")

  # Check RNA files
  if (!is.null(dataset$rna)) {
    if (!check_split_only) {
      for (file in dataset$rna$files) {
        if (!file.exists(file)) {
          warning("RNA file not found: ", file)
          all_exist <- FALSE
        } else {
          cat("  ✓ RNA file exists:", file, "\n")
        }
      }
      if (!file.exists(dataset$rna$metadata)) {
        warning("RNA metadata not found: ", dataset$rna$metadata)
        all_exist <- FALSE
      } else {
        cat("  ✓ RNA metadata exists:", dataset$rna$metadata, "\n")
      }
    } else if (!is.null(dataset$rna$source_rds)) {
      # Check if source RDS exists
      if (!file.exists(dataset$rna$source_rds)) {
        warning("Source RDS not found: ", dataset$rna$source_rds)
        all_exist <- FALSE
      } else {
        cat("  ✓ Source RDS exists:", dataset$rna$source_rds, "\n")
      }
    }
  }

  # Check lipid files
  if (!is.null(dataset$lipid)) {
    if (!file.exists(dataset$lipid$file)) {
      warning("Lipid file not found: ", dataset$lipid$file)
      all_exist <- FALSE
    } else {
      cat("  ✓ Lipid file exists:", dataset$lipid$file, "\n")
    }
  }

  # Check overview file
  if (!is.null(dataset$overview)) {
    if (!file.exists(dataset$overview)) {
      warning("Overview file not found: ", dataset$overview)
      all_exist <- FALSE
    } else {
      cat("  ✓ Overview file exists:", dataset$overview, "\n")
    }
  }

  if (all_exist) {
    cat("  ✓ All files validated successfully!\n")
  } else {
    cat("  ✗ Some files are missing (see warnings above)\n")
  }

  return(all_exist)
}
