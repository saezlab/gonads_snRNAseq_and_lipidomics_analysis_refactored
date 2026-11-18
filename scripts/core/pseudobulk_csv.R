# Lossless CSV transport for sample-level pseudobulks; never renormalizes data.
# Reading CSV requires jsonlite/digest, not Seurat or individual-nucleus inputs.
pb_chr <- function(x) as.character(unlist(x, use.names = FALSE))
pb_sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
pb_matrix_sha256 <- function(x) digest::digest(
  writeBin(as.double(x), raw(), size = 8L, endian = "little"),
  algo = "sha256", serialize = FALSE)
# R's base decimal parser can differ by one ULP for 17-digit strings. jsonlite's
# JSON number parser preserves their binary64 round trip; parse columns in bulk.
pb_parse_decimal <- function(x) {
  out <- rep(NA_real_, length(x))
  special <- !is.na(x) & x %in% c("Inf", "-Inf", "NaN")
  out[special] <- as.numeric(x[special])
  keep <- !is.na(x) & !special
  if (any(keep)) out[keep] <- as.double(jsonlite::fromJSON(
    paste0("[", paste(x[keep], collapse = ","), "]")))
  out
}
pb_file_record <- function(path, root, role) list(
  file = substring(path, nchar(root) + 2L), role = role,
  bytes = unname(file.info(path)$size), sha256 = pb_sha256(path))

pb_manifest_path <- function(path) {
  if (dir.exists(path)) path <- file.path(path, "manifest.json")
  if (basename(path) == "manifest.csv") path <- file.path(dirname(path), "manifest.json")
  if (!file.exists(path) || basename(path) != "manifest.json")
    stop("Expected a pseudobulk directory, manifest.json or manifest.csv: ", path)
  path
}

pb_source_files <- function(source) {
  if (grepl("[.]rds$", source, ignore.case = TRUE)) return(source)
  path <- pb_manifest_path(source)
  m <- jsonlite::read_json(path, simplifyVector = FALSE)
  unique(c(path, vapply(m$files, function(f) file.path(dirname(path), f$file), character(1))))
}

pb_rds_components <- function(source) {
  if (!requireNamespace("SeuratObject", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE))
    stop("Explicit local RDS import needs SeuratObject and Matrix; CSV import does not")
  o <- readRDS(source)
  if (!inherits(o, "Seurat")) stop("Expected an aggregate-profile Seurat pseudobulk source")
  a <- methods::slot(o, "assays")[["RNA"]]
  if (!"layers" %in% methods::slotNames(a)) stop("Expected Assay5 pseudobulk layers")
  raw_layers <- methods::slot(a, "layers")
  required <- c("counts", "data", "scale.data")
  if (!all(required %in% names(raw_layers))) stop("Missing required RNA layer")
  fm <- as.matrix(methods::slot(a, "features"))
  cm <- as.matrix(methods::slot(a, "cells"))
  features <- rownames(fm); cells <- rownames(cm)
  classes <- lapply(raw_layers[required], function(x) as.character(class(x)))
  layers <- lapply(required, function(n) {
    x <- as.matrix(raw_layers[[n]])
    genes <- rownames(fm)[fm[, n]]; profiles <- rownames(cm)[cm[, n]]
    stopifnot(identical(dim(x), c(length(genes), length(profiles))))
    if (!is.null(rownames(x))) stopifnot(identical(rownames(x), genes))
    if (!is.null(colnames(x))) stopifnot(identical(colnames(x), profiles))
    dimnames(x) <- list(genes, profiles)
    x
  })
  names(layers) <- required
  md <- methods::slot(o, "meta.data")
  stopifnot(identical(rownames(md), cells), !anyDuplicated(features), !anyDuplicated(cells))
  ident <- methods::slot(o, "active.ident")
  identities <- as.character(ident); names(identities) <- names(ident)
  out <- list(layers = layers, meta.data = md, metadata = md, cells = cells,
    profile_ids = cells, features = features, identities = identities,
    identity_levels = levels(ident), identity_class = as.character(class(ident)),
    recorded_commands = as.character(SeuratObject::Command(o)),
    object_version = as.character(methods::slot(o, "version")),
    layer_classes = classes, source_paths = source, assay = "RNA",
    source = list(path = source, bytes = unname(file.info(source)$size),
      sha256 = pb_sha256(source), object_class = as.character(class(o))))
  class(out) <- "pseudobulk_components"
  out
}

pb_column_schema <- function(x, name) {
  if (is.list(x)) stop("List metadata columns are not supported: ", name)
  list(name = name, storage_type = typeof(x), class = as.character(class(x)),
    levels = if (is.factor(x)) levels(x) else NULL,
    ordered = is.ordered(x), timezone = attr(x, "tzone", exact = TRUE))
}

pb_column_text <- function(x, na_token) {
  if (is.factor(x) || is.character(x)) y <- as.character(x) else
    if (is.logical(x)) y <- as.character(x) else
      if (is.numeric(x) || inherits(x, c("Date", "POSIXct")))
        y <- sprintf("%.17g", as.numeric(x)) else stop("Unsupported metadata type: ", class(x)[1])
  y[is.na(x)] <- na_token
  y
}

pb_restore_column <- function(x, schema) {
  cl <- pb_chr(schema$class)
  if ("factor" %in% cl) return(factor(x, levels = pb_chr(schema$levels), ordered = isTRUE(schema$ordered)))
  type <- schema$storage_type
  out <- switch(type, character = x, integer = as.integer(x), double = pb_parse_decimal(x),
    logical = as.logical(x), stop("Unsupported metadata storage type: ", type))
  if ("Date" %in% cl) class(out) <- "Date"
  if ("POSIXct" %in% cl) {
    class(out) <- cl
    if (!is.null(schema$timezone)) attr(out, "tzone") <- pb_chr(schema$timezone)
  }
  out
}

pb_csv_read <- function(path, verify = TRUE) {
  path <- pb_manifest_path(path)
  m <- jsonlite::read_json(path, simplifyVector = FALSE)
  stopifnot(identical(m$format, "pseudobulk_csv"), identical(as.integer(m$format_version), 1L))
  root <- dirname(path)
  for (f in m$files) {
    if (startsWith(f$file, "/") || ".." %in% strsplit(f$file, "/", fixed = TRUE)[[1]])
      stop("Nonportable payload path")
    p <- file.path(root, f$file)
    stopifnot(file.exists(p), unname(file.info(p)$size) == f$bytes)
    if (grepl("[.]csv$", p)) stopifnot(file.info(p)$size < m$maximum_csv_bytes)
    if (verify) stopifnot(identical(pb_sha256(p), f$sha256))
  }
  layers <- lapply(m$layers, function(z) {
    pieces <- lapply(z$chunks, function(f) {
      tab <- utils::read.csv(file.path(root, f$file), check.names = FALSE,
        colClasses = "character",
        na.strings = m$na_token, stringsAsFactors = FALSE)
      stopifnot(identical(tab[[1]], pb_chr(z$row_ids)))
      mat <- do.call(cbind, lapply(tab[-1], pb_parse_decimal)); rownames(mat) <- tab[[1]]
      stopifnot(nrow(mat) == f$n_rows, ncol(mat) == f$n_columns,
        identical(colnames(mat), pb_chr(z$column_ids)[seq.int(f$column_start, f$column_end)]))
      mat
    })
    x <- do.call(cbind, pieces)
    stopifnot(identical(dim(x), as.integer(pb_chr(z$dimensions))),
      identical(rownames(x), pb_chr(z$row_ids)), identical(colnames(x), pb_chr(z$column_ids)))
    if (verify) stopifnot(identical(pb_matrix_sha256(x), z$content_sha256))
    x
  })
  md_raw <- utils::read.csv(file.path(root, m$metadata$file), colClasses = "character",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = m$na_token)
  stopifnot(identical(md_raw[[1]], pb_chr(m$profile_ids)))
  md <- as.data.frame(lapply(seq_along(m$metadata$columns), function(i)
    pb_restore_column(md_raw[[i + 1L]], m$metadata$columns[[i]])),
    stringsAsFactors = FALSE, check.names = FALSE)
  names(md) <- vapply(m$metadata$columns, function(z) z$name, character(1))
  rownames(md) <- md_raw[[1]]
  it <- utils::read.csv(file.path(root, m$identities$file), colClasses = "character",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = m$na_token)
  identities <- stats::setNames(it$identity, it$profile_id)
  stopifnot(identical(names(identities), pb_chr(m$identities$names)))
  out <- list(layers = layers, meta.data = md, metadata = md,
    cells = pb_chr(m$profile_ids), profile_ids = pb_chr(m$profile_ids),
    features = pb_chr(m$features), identities = identities,
    identity_levels = pb_chr(m$identities$levels), identity_class = pb_chr(m$identities$class),
    recorded_commands = pb_chr(m$recorded_commands), object_version = m$source$object_version,
    layer_classes = lapply(m$layers, function(z) pb_chr(z$original_class)),
    source_paths = pb_source_files(path), assay = m$assay, source = m$source, manifest = m)
  class(out) <- c("pseudobulk_csv", "pseudobulk_components")
  out
}

pb_read <- function(path, verify = TRUE) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) pb_rds_components(path) else pb_csv_read(path, verify)
}
pb_components <- pb_read

pb_csv_write <- function(source_rds, output, max_bytes = 50000000L) {
  if (file.exists(file.path(output, "manifest.json"))) stop("Preserve existing transport: ", output)
  if (!is.finite(max_bytes) || max_bytes < 500L || max_bytes > 50000000L) stop("Invalid CSV byte limit")
  x <- pb_rds_components(source_rds)
  stopifnot(ncol(x$layers$counts) <= 1000L, "orig.ident" %in% names(x$meta.data),
    any(c("sample_name", "cf453_id") %in% names(x$meta.data)),
    identical(as.character(x$meta.data$orig.ident), x$profile_ids),
    all(is.finite(x$layers$counts)), all(x$layers$counts >= 0),
    all(x$layers$counts == round(x$layers$counts)))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- sub("/+$", "", output)
  na_token <- paste0("__PB_MISSING_", substr(x$source$sha256, 1L, 20L), "__")
  if (any(vapply(x$meta.data, function(z) any(as.character(z) == na_token, na.rm = TRUE), logical(1))))
    stop("Missing-value token collides with metadata")
  quote_csv <- function(z) paste0('"', gsub('"', '""', z, fixed = TRUE), '"')
  records <- list(); layer_records <- list()
  for (name in names(x$layers)) {
    mat <- x$layers[[name]]
    stopifnot(is.matrix(mat), all(is.finite(mat)), !anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)))
    parts <- list(); number <- 1L
    write_part <- function(first, last) {
      file <- sprintf("%s_part%03d.csv", name, number)
      path <- file.path(output, file)
      ids <- seq.int(first, last)
      tab <- as.data.frame(lapply(ids, function(j) {
        if (name == "counts") sprintf("%.0f", mat[, j]) else sprintf("%.17g", mat[, j])
      }), stringsAsFactors = FALSE, check.names = FALSE)
      tab <- cbind(data.frame(.feature_id = rownames(mat), stringsAsFactors = FALSE), tab)
      con <- file(path, open = "wt", encoding = "UTF-8")
      writeLines(paste(quote_csv(c(".feature_id", colnames(mat)[ids])), collapse = ","), con)
      utils::write.table(tab, con, sep = ",", quote = 1L, row.names = FALSE,
        col.names = FALSE, qmethod = "double", na = na_token)
      close(con)
      if (file.info(path)$size >= max_bytes) {
        unlink(path)
        if (first == last) stop("One column exceeds CSV size limit; row splitting needed")
        middle <- floor((first + last) / 2)
        write_part(first, middle); write_part(middle + 1L, last)
        return(invisible(NULL))
      }
      rec <- pb_file_record(path, output, paste0("layer:", name))
      rec$column_start <- first; rec$column_end <- last
      rec$n_rows <- nrow(mat); rec$n_columns <- length(ids)
      parts[[length(parts) + 1L]] <<- rec
      number <<- number + 1L
    }
    write_part(1L, ncol(mat))
    layer_records[[name]] <- list(dimensions = dim(mat), row_ids = rownames(mat),
      column_ids = colnames(mat), original_class = x$layer_classes[[name]],
      storage_mode = typeof(mat), content_sha256 = pb_matrix_sha256(mat),
      content_encoding = "IEEE754-binary64-little-endian-column-major", chunks = parts)
    records <- c(records, parts)
  }
  metadata <- x$meta.data
  schemas <- lapply(seq_along(metadata), function(i) pb_column_schema(metadata[[i]], names(metadata)[i]))
  md_text <- as.data.frame(lapply(metadata, pb_column_text, na_token = na_token),
    stringsAsFactors = FALSE, check.names = FALSE)
  md_text <- cbind(data.frame(.profile_id = rownames(metadata), stringsAsFactors = FALSE), md_text)
  md_path <- file.path(output, "metadata_part001.csv")
  utils::write.table(md_text, md_path, sep = ",", row.names = FALSE, quote = TRUE,
    na = na_token, qmethod = "double", fileEncoding = "UTF-8")
  md_record <- pb_file_record(md_path, output, "metadata")
  md_record$n_rows <- nrow(metadata); md_record$n_columns <- ncol(metadata)
  md_record$columns <- schemas
  records <- c(records, list(md_record))
  id_path <- file.path(output, "identities_part001.csv")
  utils::write.table(data.frame(profile_id = names(x$identities), identity = unname(x$identities)),
    id_path, sep = ",", row.names = FALSE, quote = TRUE, na = na_token,
    qmethod = "double", fileEncoding = "UTF-8")
  id_record <- pb_file_record(id_path, output, "identities")
  id_record$names <- names(x$identities); id_record$levels <- x$identity_levels
  id_record$class <- x$identity_class
  records <- c(records, list(id_record))
  index <- do.call(rbind, lapply(records, function(f) data.frame(file = f$file, role = f$role,
    bytes = f$bytes, sha256 = f$sha256, source_rds_path = x$source$path,
    source_rds_sha256 = x$source$sha256, source_rds_bytes = x$source$bytes)))
  index_path <- file.path(output, "manifest.csv")
  utils::write.csv(index, index_path, row.names = FALSE, fileEncoding = "UTF-8")
  records <- c(records, list(pb_file_record(index_path, output, "file_index")))
  stopifnot(all(vapply(records, function(f) f$bytes < max_bytes, logical(1))))
  source <- x$source; source$object_version <- x$object_version
  m <- list(format = "pseudobulk_csv", format_version = 1L, assay = "RNA",
    profile_unit = "sample-level pseudobulk", individual_nuclei = FALSE,
    maximum_csv_bytes = as.integer(max_bytes), numeric_encoding = "17 significant decimal digits; counts integral",
    source = source, features = x$features, profile_ids = x$profile_ids,
    recorded_commands = x$recorded_commands, na_token = na_token,
    layers = layer_records, metadata = md_record, identities = id_record, files = records)
  manifest_path <- file.path(output, "manifest.json")
  jsonlite::write_json(m, manifest_path, pretty = TRUE, auto_unbox = TRUE, digits = NA, null = "null")
  imported <- pb_csv_read(manifest_path)
  stopifnot(identical(imported$meta.data, x$meta.data), identical(imported$identities, x$identities),
    identical(imported$identity_levels, x$identity_levels), identical(imported$recorded_commands, x$recorded_commands))
  for (name in names(x$layers)) stopifnot(identical(imported$layers[[name]], x$layers[[name]]))
  invisible(m)
}
