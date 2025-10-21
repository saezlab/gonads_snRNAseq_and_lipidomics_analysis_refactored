#' Contrast Definitions
#'
#' This file defines all contrasts used in the FCG liver analysis project.
#' Contrasts specify which groups to compare in differential analysis.
#'
#' IMPORTANT: Understanding Contrast Vectors
#' ==========================================
#' Contrast vectors are numeric vectors that specify how to combine conditions.
#' The position in the vector corresponds to the order of unique(targets$condition).
#'
#' Examples:
#' - If unique(targets$condition) = c("Ovaries", "Testes")
#'   Then c(-1, 2) means: -1 * Ovaries + 2 * Testes = 2 * Testes - Ovaries
#'   This effectively computes Testes - Ovaries (the "2" is just a scaling factor)
#'
#' - If unique(targets$condition) = c("XXO", "XXT", "XYO", "XYT")
#'   Then c(1, -4) means: 1 * XXO + (-4) * XYT = XXO - XYT (scaled)
#'
#' The sign indicates direction:
#' - Positive coefficient: add this condition
#' - Negative coefficient: subtract this condition
#'
#' For complex contrasts:
#' - c(2, 2, -1, -3) means: (2*Cond1 + 2*Cond2 - 1*Cond3 - 3*Cond4)
#'   This averages Cond1 and Cond2, then subtracts average of Cond3 and Cond4
#'
#' See docs/CONTRAST_DEFINITIONS_EXPLAINED.md for detailed explanation.

#' @description
#' CONTRASTS is a named list where each element defines a contrast configuration.
#' Each contrast includes:
#' - name: Human-readable name
#' - group_variable: Column in metadata to use for grouping
#' - vector: Numeric contrast vector (see above for interpretation)
#' - description: What this contrast tests biologically
#' - colors: Color scheme for visualizations (optional)
CONTRASTS <- list(

  # Main Contrasts: Gonadal Hormones vs. Sex Chromosomes
  # =====================================================

  gonadal_type = list(
    name = "Testes vs Ovaries",
    group_variable = "gonadal_type",
    vector = c(-1, 2),
    description = "Effect of gonadal hormones (primary comparison in FCG model)",
    interpretation = "Positive logFC = higher in testes, negative logFC = higher in ovaries",
    colors = list(
      Ovaries = "red",
      Testes = "royalblue"
    ),
    notes = "This is the PRIMARY comparison. ~4000 DEGs in hepatocytes at FDR < 0.05"
  ),

  allosome_type = list(
    name = "XY vs XX",
    group_variable = "allosome_type",
    vector = c(-1, 2),
    description = "Effect of sex chromosomes (Y chromosome vs. second X)",
    interpretation = "Positive logFC = higher in XY, negative logFC = higher in XX",
    colors = list(
      XX = "red",
      XY = "royalblue"
    ),
    notes = "Minimal effect: few DEGs, mainly Y-linked genes (Ddx3y, Eif2s3y, Kdm5d) and X-escape genes"
  ),

  # Group-Specific Comparisons (Four Core Genotypes)
  # =================================================
  # FCG groups: XXO (XX + Ovaries), XXT (XX + Testes), XYO (XY + Ovaries), XYT (XY + Testes)

  XXT_vs_XXO = list(
    name = "XXT vs XXO",
    group_variable = "group",
    vector = c(1, -4),  # Assuming order: XXO, XXT, XYO, XYT
    description = "Effect of testes vs ovaries in XX background",
    interpretation = "Isolates gonadal hormone effect with sex chromosomes held constant (XX)",
    colors = list(
      XXO = "red",
      XXT = "royalblue"
    )
  ),

  XYT_vs_XYO = list(
    name = "XYT vs XYO",
    group_variable = "group",
    vector = c(3, -2),  # Assuming order: XXO, XXT, XYO, XYT
    description = "Effect of testes vs ovaries in XY background",
    interpretation = "Isolates gonadal hormone effect with sex chromosomes held constant (XY)",
    colors = list(
      XYO = "red",
      XYT = "royalblue"
    )
  ),

  XYO_vs_XXO = list(
    name = "XYO vs XXO",
    group_variable = "group",
    vector = c(1, -3),  # Assuming order: XXO, XXT, XYO, XYT
    description = "Effect of XY vs XX in ovaries background",
    interpretation = "Isolates sex chromosome effect with gonads held constant (ovaries)",
    colors = list(
      XXO = "brown",
      XYO = "cyan"
    )
  ),

  XYT_vs_XXT = list(
    name = "XYT vs XXT",
    group_variable = "group",
    vector = c(2, -4),  # Assuming order: XXO, XXT, XYO, XYT
    description = "Effect of XY vs XX in testes background",
    interpretation = "Isolates sex chromosome effect with gonads held constant (testes)",
    colors = list(
      XXT = "brown",
      XYT = "cyan"
    )
  ),

  # Gonadectomy + Hormone Replacement (GDXHR) Contrasts
  # ====================================================
  # Used for validation experiments

  gdxhr_treatment = list(
    name = "Treatment Effect (GDXHR)",
    group_variable = "treatment",
    vector = c(-1, 2),
    description = "Effect of hormone treatment after gonadectomy",
    interpretation = "Tests activational (not organizational) hormone effects",
    notes = "Correlation with FCG gonadal_type: r = -0.48 (opposite directionality expected)"
  ),

  # Zonation-Specific Contrasts
  # ============================
  # For hepatocytes subdivided by spatial zonation

  zonation_gonadal = list(
    name = "Testes vs Ovaries (Zonation-specific)",
    group_variable = "gonadal_type",
    vector = c(-1, 2),
    description = "Gonadal effect in zonation-subdivided hepatocytes",
    interpretation = "Reveals periportal vs pericentral differences in gonad response",
    notes = "43% of gonad-dependent genes are pericentral-biased, 18% periportal-biased"
  ),

  zonation_pericentral_vs_periportal = list(
    name = "Pericentral vs Periportal",
    group_variable = "zonation",
    vector = c(-1, 2),
    description = "Spatial zonation effect (independent of gonad)",
    interpretation = "Identifies zonation markers: Cyp2e1+ (pericentral), Cyp2f2+ (periportal)",
    colors = list(
      Periportal = "#8B4513",
      Pericentral = "#4169E1"
    )
  ),

  # Validation Experiment Contrasts
  # ================================

  validation_knockdown_vs_control = list(
    name = "siRNA Knockdown vs Control",
    group_variable = "treatment",
    vector = c(-1, 2),
    description = "Effect of gene knockdown on lipid abundance",
    interpretation = "Validates gene-lipid associations from MOFA",
    notes = "Applied to: Cyp2b9, Cyp2b13, Cyp2c38, Cyp2d9, Sult3a1, Mtmr11, Abcb1a"
  )
)


#' Get Contrast Configuration
#'
#' Helper function to retrieve a specific contrast configuration.
#'
#' @param contrast_name Character. Name of the contrast (e.g., "gonadal_type")
#'
#' @return List containing the contrast configuration
#'
#' @examples
#' gonadal_contrast <- get_contrast("gonadal_type")
#' print(gonadal_contrast$description)
get_contrast <- function(contrast_name) {
  if (!contrast_name %in% names(CONTRASTS)) {
    stop("Contrast '", contrast_name, "' not found. Available contrasts: ",
         paste(names(CONTRASTS), collapse = ", "))
  }
  return(CONTRASTS[[contrast_name]])
}


#' List Available Contrasts
#'
#' Print a summary of all available contrasts with descriptions.
#'
#' @examples
#' list_contrasts()
list_contrasts <- function() {
  cat("Available Contrasts:\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  # Group contrasts by category
  categories <- list(
    "Main Contrasts" = c("gonadal_type", "allosome_type"),
    "FCG Group Comparisons" = c("XXT_vs_XXO", "XYT_vs_XYO", "XYO_vs_XXO", "XYT_vs_XXT"),
    "Validation Experiments" = c("gdxhr_treatment", "validation_knockdown_vs_control"),
    "Zonation-Specific" = c("zonation_gonadal", "zonation_pericentral_vs_periportal")
  )

  for (category in names(categories)) {
    cat(category, "\n")
    cat(paste(rep("-", 80), collapse = ""), "\n")

    for (contrast_name in categories[[category]]) {
      if (contrast_name %in% names(CONTRASTS)) {
        contrast <- CONTRASTS[[contrast_name]]
        cat("  ", contrast_name, "\n")
        cat("    Name:       ", contrast$name, "\n")
        cat("    Variable:   ", contrast$group_variable, "\n")
        cat("    Vector:     ", paste(contrast$vector, collapse = ", "), "\n")
        cat("    Description:", contrast$description, "\n")
        if (!is.null(contrast$notes)) {
          cat("    Notes:      ", contrast$notes, "\n")
        }
        cat("\n")
      }
    }
  }

  cat(paste(rep("=", 80), collapse = ""), "\n")
}


#' Build Contrast List for Differential Analysis
#'
#' Creates a named list of contrast vectors suitable for run_differential_analysis().
#'
#' @param contrast_names Character vector. Names of contrasts to include.
#' @param use_formal_names Logical. Use formal contrast names (e.g., "Testes vs Ovaries")
#'                         instead of config names (e.g., "gonadal_type"). Default: TRUE
#'
#' @return Named list of contrast vectors.
#'
#' @examples
#' # Single contrast
#' comparisons <- build_contrast_list("gonadal_type")
#'
#' # Multiple contrasts
#' comparisons <- build_contrast_list(c("gonadal_type", "allosome_type"))
#'
#' # Use config names instead of formal names
#' comparisons <- build_contrast_list("gonadal_type", use_formal_names = FALSE)
build_contrast_list <- function(contrast_names, use_formal_names = TRUE) {

  comparisons <- list()

  for (contrast_name in contrast_names) {
    contrast <- get_contrast(contrast_name)

    if (use_formal_names) {
      # Use formal name (e.g., "Testes vs Ovaries")
      list_name <- contrast$name
    } else {
      # Use config name (e.g., "gonadal_type")
      list_name <- contrast_name
    }

    comparisons[[list_name]] <- contrast$vector
  }

  return(comparisons)
}


#' Get Colors for Contrast
#'
#' Retrieves color scheme for a specific contrast.
#'
#' @param contrast_name Character. Name of the contrast.
#'
#' @return Named vector of colors, or NULL if no colors defined.
#'
#' @examples
#' colors <- get_contrast_colors("gonadal_type")
#' # Returns: c(Ovaries = "red", Testes = "royalblue")
get_contrast_colors <- function(contrast_name) {
  contrast <- get_contrast(contrast_name)
  return(contrast$colors)
}


#' Validate Contrast Vector
#'
#' Checks if a contrast vector is valid for a given set of conditions.
#'
#' @param contrast_vector Numeric vector. The contrast to validate.
#' @param conditions Character vector. The conditions (from unique(targets$condition)).
#'
#' @return Logical. TRUE if valid, FALSE otherwise (with warning).
#'
#' @examples
#' validate_contrast_vector(c(-1, 2), c("Ovaries", "Testes"))  # TRUE
#' validate_contrast_vector(c(-1, 2, 3), c("Ovaries", "Testes"))  # FALSE (length mismatch)
validate_contrast_vector <- function(contrast_vector, conditions) {

  n_conditions <- length(conditions)
  n_vector <- length(contrast_vector)

  if (n_vector != n_conditions) {
    warning("Contrast vector length (", n_vector, ") does not match number of conditions (",
            n_conditions, "): ", paste(conditions, collapse = ", "))
    return(FALSE)
  }

  # Check that at least one positive and one negative coefficient exists
  if (all(contrast_vector >= 0) || all(contrast_vector <= 0)) {
    warning("Contrast vector should have both positive and negative coefficients")
    return(FALSE)
  }

  return(TRUE)
}
