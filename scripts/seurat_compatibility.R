################################################################################
# Seurat v4/v5 Compatibility Helper Functions
################################################################################
#
# This file provides compatibility functions for working with both Seurat v4
# (Assay class) and Seurat v5 (Assay5 class).
#
# Key differences:
#   Seurat v4: Data stored in @counts, @data, @scale.data slots
#   Seurat v5: Data stored in layers accessed via LayerData()
#
# Usage:
#   source("scripts/seurat_compatibility.R")
#   scaled_data <- get_layer_data(seurat_obj, assay = "RNA", layer = "scale.data")
#
################################################################################

#' Get layer data from Seurat object (compatible with v4 and v5)
#'
#' @param seurat_obj A Seurat object
#' @param assay Name of the assay (default: "RNA")
#' @param layer Layer to extract: "counts", "data", or "scale.data"
#' @return Matrix of expression data
get_layer_data <- function(seurat_obj, assay = "RNA", layer = "scale.data") {

  # Check if assay exists
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object")
  }

  # Determine Seurat version by checking assay class
  assay_obj <- seurat_obj@assays[[assay]]
  assay_class <- class(assay_obj)[1]

  if (assay_class == "Assay5") {
    # Seurat v5: Use LayerData()
    tryCatch({
      data <- LayerData(seurat_obj, assay = assay, layer = layer)
      return(data)
    }, error = function(e) {
      stop("Failed to extract layer '", layer, "' from Seurat v5 object: ", e$message)
    })

  } else if (assay_class == "Assay") {
    # Seurat v4: Access slots directly
    data <- switch(layer,
                   "counts" = assay_obj@counts,
                   "data" = assay_obj@data,
                   "scale.data" = assay_obj@scale.data,
                   stop("Unknown layer: ", layer))

    if (is.null(data) || length(data) == 0) {
      stop("Layer '", layer, "' is empty or NULL in Seurat v4 object")
    }

    return(data)

  } else {
    stop("Unknown assay class: ", assay_class)
  }
}


#' Check if layer exists in Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param assay Name of the assay (default: "RNA")
#' @param layer Layer to check: "counts", "data", or "scale.data"
#' @return Logical indicating if layer exists and has data
has_layer_data <- function(seurat_obj, assay = "RNA", layer = "scale.data") {

  # Check if assay exists
  if (!assay %in% names(seurat_obj@assays)) {
    return(FALSE)
  }

  # Determine Seurat version
  assay_obj <- seurat_obj@assays[[assay]]
  assay_class <- class(assay_obj)[1]

  if (assay_class == "Assay5") {
    # Seurat v5: Check if layer exists
    tryCatch({
      layers <- Layers(seurat_obj, assay = assay)
      return(layer %in% layers)
    }, error = function(e) {
      return(FALSE)
    })

  } else if (assay_class == "Assay") {
    # Seurat v4: Check if slot exists and has data
    data <- switch(layer,
                   "counts" = assay_obj@counts,
                   "data" = assay_obj@data,
                   "scale.data" = assay_obj@scale.data,
                   NULL)

    return(!is.null(data) && length(data) > 0)

  } else {
    return(FALSE)
  }
}


#' Get Seurat version information
#'
#' @param seurat_obj A Seurat object
#' @param assay Name of the assay (default: "RNA")
#' @return List with version info
get_seurat_version_info <- function(seurat_obj, assay = "RNA") {

  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object")
  }

  assay_obj <- seurat_obj@assays[[assay]]
  assay_class <- class(assay_obj)[1]

  info <- list(
    assay_class = assay_class,
    is_v5 = (assay_class == "Assay5"),
    is_v4 = (assay_class == "Assay")
  )

  # Get available data
  if (info$is_v5) {
    info$available_layers <- tryCatch(
      Layers(seurat_obj, assay = assay),
      error = function(e) character(0)
    )
  } else if (info$is_v4) {
    info$available_layers <- c(
      if (!is.null(assay_obj@counts) && length(assay_obj@counts) > 0) "counts",
      if (!is.null(assay_obj@data) && length(assay_obj@data) > 0) "data",
      if (!is.null(assay_obj@scale.data) && length(assay_obj@scale.data) > 0) "scale.data"
    )
  }

  return(info)
}


#' Print Seurat object information
#'
#' @param seurat_obj A Seurat object
print_seurat_info <- function(seurat_obj) {
  cat("Seurat object information:\n")
  cat("  Samples/Cells:", ncol(seurat_obj), "\n")
  cat("  Genes:", nrow(seurat_obj), "\n")
  cat("  Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n\n")

  for (assay_name in names(seurat_obj@assays)) {
    cat("  Assay:", assay_name, "\n")

    info <- get_seurat_version_info(seurat_obj, assay = assay_name)
    cat("    Class:", info$assay_class, "\n")
    cat("    Available layers:", paste(info$available_layers, collapse = ", "), "\n")
  }

  cat("\n")
}
