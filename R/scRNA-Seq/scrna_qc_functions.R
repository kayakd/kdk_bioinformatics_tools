# Purpose: Reusable, sparse-count Seurat QC and plotting. Remove only all-zero
#          features, calculate RNA/ADT metrics, and flag OSCA-style outliers.
#          No cells are removed unless qc_subset_cells() is explicitly called.
# Input:   Seurat objects, named lists of Seurat objects, Cell Ranger H5/MTX
#          matrices, or per-cell metadata tables. Use unnormalized RNA counts.
# Output:  Seurat objects, audit/threshold/retention tables, and ggplot objects.
# Author:  Dr. Koray Doğan Kaya
# Date:    2026-10-20
# Dependencies: Seurat >= 4, SeuratObject, Matrix, ggplot2; scuttle for MAD flags.
# References:
# https://bioconductor.org/books/release/OSCA.basic/quality-control.html
# https://bioconductor.org/books/release/OSCA.advanced/quality-control-redux.html
# https://satijalab.org/seurat/articles/seurat5_essential_commands.html
# Scope: In-memory base/Matrix counts, including disjoint Seurat v5 counts
#        layers. No normalization, integration, doublet calling, or ambient
#        RNA correction. Sourcing this file does not install or run anything.

qc_check_packages <- function(report = FALSE, h5 = FALSE) {
  packages <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "scuttle")
  if (report) packages <- c(packages, "rmarkdown", "knitr")
  if (h5) packages <- c(packages, "hdf5r")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". See README.md. Nothing was installed.", call. = FALSE)
  if (utils::packageVersion("Seurat") < "4.0.0") stop("Seurat >= 4 is required.")
  if (utils::packageVersion("ggplot2") < "3.4.0") stop("ggplot2 >= 3.4.0 is required.")
  invisible(TRUE)
}

.qc_counts_check <- function(x) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    stop("QC requires in-memory matrix/Matrix counts. On-disk backends are not supported.")
  }
  if (length(dim(x)) != 2L || any(dim(x) == 0L)) stop("Empty count matrix.")
  if (is.null(rownames(x)) || is.null(colnames(x)) ||
      anyNA(rownames(x)) || anyNA(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x)))) {
    stop("Count matrices need nonempty, unique feature and cell names.")
  }
  values <- if (inherits(x, "sparseMatrix") && "x" %in% methods::slotNames(x)) {
    methods::slot(x, "x")
  } else if (is.matrix(x)) x else {
    stop("Use a numeric sparse Matrix or a numeric base matrix.")
  }
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0) ||
      any(abs(values - round(values)) > 1e-6)) {
    stop("Expected finite, nonnegative, integer raw counts, not normalized data.")
  }
  invisible(TRUE)
}

# Access every counts layer, never just the first match in a merged v5 object.
.qc_layers <- function(object, assay = "RNA") {
  if (!inherits(object, "Seurat")) stop("Expected a Seurat object.")
  if (!assay %in% SeuratObject::Assays(object)) stop("Missing assay: ", assay)
  a <- object[[assay]]
  if (inherits(a, "Assay5")) {
    layers <- SeuratObject::Layers(a, search = NA)
    layers <- layers[grepl("^counts($|\\.)", layers)]
    if (!length(layers)) stop("No raw counts layers in assay ", assay)
    result <- lapply(layers, function(layer) SeuratObject::LayerData(a, layer = layer))
    names(result) <- layers
  } else {
    result <- list(counts = SeuratObject::GetAssayData(a, slot = "counts"))
  }
  invisible(lapply(result, .qc_counts_check))
  cells <- unlist(lapply(result, colnames), use.names = FALSE)
  if (anyDuplicated(cells)) {
    stop("Counts layers overlap in cell names in ", assay,
         ". Resolve overlapping layers before QC; counts will not be double-counted.")
  }
  result
}

# Remove zero-total features across ALL count layers of each requested assay.
# Cells and other assays are not subset by the RNA feature names.
qc_drop_zero_features <- function(object, assays = c("RNA", "ADT")) {
  original_cells <- colnames(object)
  present <- SeuratObject::Assays(object)
  if (is.null(object@misc$qc_zero_features)) object@misc$qc_zero_features <- list()
  for (assay in intersect(assays, present)) {
    layers <- .qc_layers(object, assay)
    keep <- unique(unlist(lapply(layers, function(m) {
      rownames(m)[Matrix::rowSums(m) > 0]
    }), use.names = FALSE))
    before <- rownames(object[[assay]])
    removed <- setdiff(before, keep)
    object@misc$qc_zero_features[[assay]] <- unique(c(
      object@misc$qc_zero_features[[assay]], removed))
    if (!length(keep)) {
      if (assay == "RNA" || assay == SeuratObject::DefaultAssay(object)) {
        stop("Assay ", assay, " contains no expressed features; inspect the input.")
      }
      warning("Removing entirely zero assay ", assay, "; all cells are retained.")
      object[[assay]] <- NULL
      object[[paste0("nCount_", assay)]] <- 0
      object[[paste0("nFeature_", assay)]] <- 0
    } else if (length(removed)) {
      object[[assay]] <- subset(object[[assay]], features = keep)
    }
  }
  if (!identical(original_cells, colnames(object))) stop("Unexpected cell loss during feature cleanup.")
  object
}

# Read one cell-called sample. Raw droplets require separate cell calling.
qc_read_10x <- function(path, sample_id, metadata = list(), keep_adt = TRUE) {
  if (length(sample_id) != 1L || is.na(sample_id) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", sample_id)) {
    stop("sample_id must be one nonempty filename-safe identifier.")
  }
  if (!file.exists(path)) stop("Input does not exist: ", path)
  if (grepl("(^|/)raw_feature_bc_matrix", path)) {
    stop("Use a cell-called filtered matrix here. Raw droplets need cell calling first.")
  }
  if (dir.exists(path)) {
    dat <- Seurat::Read10X(path, unique.features = TRUE)
  } else {
    if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Install hdf5r to read Cell Ranger H5 files.")
    dat <- Seurat::Read10X_h5(path, use.names = TRUE, unique.features = TRUE)
  }
  if (is.list(dat)) {
    if (!"Gene Expression" %in% names(dat)) {
      stop("No unambiguous 'Gene Expression' matrix. Available: ", paste(names(dat), collapse = ", "))
    }
    rna <- dat[["Gene Expression"]]
    adt <- dat[["Antibody Capture"]]
  } else {
    rna <- dat
    adt <- NULL
  }
  .qc_counts_check(rna)
  object <- SeuratObject::CreateSeuratObject(
    counts = rna, assay = "RNA", project = sample_id, min.cells = 0, min.features = 0)
  if (!identical(colnames(object), colnames(rna))) stop("Unexpected cell change on object creation.")
  if (keep_adt && !is.null(adt)) {
    .qc_counts_check(adt)
    if (!setequal(colnames(adt), colnames(rna))) stop("RNA and ADT cell barcodes differ.")
    adt <- adt[, colnames(rna), drop = FALSE]
    object[["ADT"]] <- if (inherits(object[["RNA"]], "Assay5")) {
      SeuratObject::CreateAssay5Object(counts = adt, min.cells = 0, min.features = 0)
    } else {
      SeuratObject::CreateAssayObject(counts = adt, min.cells = 0, min.features = 0)
    }
  }
  object$barcode_10x <- colnames(object)
  object$sample_id <- sample_id
  object <- SeuratObject::RenameCells(object, new.names = paste(sample_id, colnames(object), sep = "__"))
  if (length(metadata)) {
    if (is.null(names(metadata)) || anyDuplicated(names(metadata)) || any(!nzchar(names(metadata)))) {
      stop("metadata must be a named list of scalar values.")
    }
    if (any(names(metadata) %in% colnames(object[[]]))) stop("Metadata would overwrite an existing field.")
    for (nm in names(metadata)) {
      if (length(metadata[[nm]]) != 1L) stop("Each metadata value must be scalar: ", nm)
      object[[nm]] <- metadata[[nm]]
    }
  }
  object@misc$qc_input <- normalizePath(path, mustWork = TRUE)
  qc_drop_zero_features(object)
}

# Mitochondrial features must be identified, not silently assumed to be absent.
# Human/mouse symbols work with the default. For other IDs supply mito_features.
qc_add_metrics <- function(object, mito_pattern = "^[Mm][Tt]-", mito_features = NULL,
                           assay = "RNA", adt_assay = "ADT") {
  layers <- .qc_layers(object, assay)
  all_cells <- colnames(object)
  count_cells <- unlist(lapply(layers, colnames), use.names = FALSE)
  if (!setequal(all_cells, count_cells)) stop("Every object cell must have RNA counts for this QC report.")
  retained <- unique(unlist(lapply(layers, rownames), use.names = FALSE))
  known <- unique(c(retained, object@misc$qc_zero_features[[assay]]))
  if (is.null(mito_features)) {
    if (is.null(mito_pattern) || length(mito_pattern) != 1L || !nzchar(mito_pattern)) {
      stop("Supply mito_pattern or explicit mito_features.")
    }
    recognized <- grep(mito_pattern, known, value = TRUE)
  } else {
    mito_features <- unique(as.character(mito_features))
    if (anyNA(mito_features) || any(!nzchar(mito_features))) stop("Invalid mitochondrial feature list.")
    recognized <- intersect(mito_features, known)
  }
  if (!length(recognized)) {
    stop("No mitochondrial feature names matched. Inspect gene naming; for Ensembl IDs, ",
         "supply a mitochondrial feature list from the matching reference annotation. ",
         "Mitochondrial percentages were NOT set to zero.")
  }
  counts <- detected <- mito <- stats::setNames(rep(NA_real_, length(all_cells)), all_cells)
  for (m in layers) {
    cells <- colnames(m)
    counts[cells] <- Matrix::colSums(m)
    detected[cells] <- Matrix::colSums(m > 0)
    features <- intersect(recognized, rownames(m))
    mito[cells] <- if (length(features)) Matrix::colSums(m[features, , drop = FALSE]) else 0
  }
  percent <- rep(NA_real_, length(counts))
  valid <- counts > 0
  percent[valid] <- 100 * mito[valid] / counts[valid]
  metrics <- data.frame(row.names = all_cells)
  metrics[[paste0("nCount_", assay)]] <- unname(counts)
  metrics[[paste0("nFeature_", assay)]] <- unname(detected)
  metrics$mito_umis <- unname(mito)
  metrics$percent.mt <- percent
  metrics$qc_zero_rna <- unname(counts == 0)
  object <- SeuratObject::AddMetaData(object, metrics)
  if (adt_assay %in% SeuratObject::Assays(object)) {
    ac <- af <- stats::setNames(rep(NA_real_, length(all_cells)), all_cells)
    for (m in .qc_layers(object, adt_assay)) {
      ac[colnames(m)] <- Matrix::colSums(m)
      af[colnames(m)] <- Matrix::colSums(m > 0)
    }
    object[[paste0("nCount_", adt_assay)]] <- unname(ac)
    object[[paste0("nFeature_", adt_assay)]] <- unname(af)
  }
  object@misc$qc_mito_features <- recognized
  object@misc$qc_mito_detected <- intersect(recognized, retained)
  object
}

.qc_bind_rows <- function(tables) {
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(d) {
    for (nm in setdiff(columns, names(d))) d[[nm]] <- NA
    d[, columns, drop = FALSE]
  })
  result <- do.call(rbind, tables)
  rownames(result) <- NULL
  result
}

# Every plot accepts an object, named object list, or this metadata table.
qc_metadata <- function(x, sample_col = "sample_id") {
  if (inherits(x, "Seurat")) {
    d <- x[[]]
    d$cell_id <- rownames(d)
  } else if (is.data.frame(x)) {
    d <- x
    if (!"cell_id" %in% names(d)) d$cell_id <- rownames(d)
  } else if (is.list(x) && length(x) && all(vapply(x, function(z) inherits(z, "Seurat") || is.data.frame(z), logical(1)))) {
    d <- .qc_bind_rows(lapply(x, qc_metadata, sample_col = sample_col))
  } else stop("Expected a Seurat object, a nonempty list of objects/metadata tables, or a metadata data.frame.")
  if (!sample_col %in% names(d) || anyNA(d[[sample_col]]) || any(!nzchar(as.character(d[[sample_col]])))) {
    stop("Missing or invalid sample metadata: ", sample_col)
  }
  if (!nrow(d) || anyNA(d$cell_id) || anyDuplicated(d$cell_id)) {
    stop("Metadata must have cells with unique cell_id values. Prefix repeated barcodes with sample IDs.")
  }
  d[[sample_col]] <- factor(as.character(d[[sample_col]]), levels = unique(as.character(d[[sample_col]])))
  d
}

# OSCA diagnostics: lower log-count/log-feature tails, upper untransformed mt%.
# Small groups are not fitted. This safeguard is a toolkit choice, not an OSCA cutoff.
qc_flag_osca <- function(x, sample_col = "sample_id", nmads = 3, min_cells = 50) {
  if (!requireNamespace("scuttle", quietly = TRUE)) stop("Install Bioconductor package scuttle.")
  if (length(nmads) != 1L || !is.finite(nmads) || nmads <= 0) stop("nmads must be positive.")
  if (length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 3) stop("min_cells must be >= 3.")
  d <- qc_metadata(x, sample_col)
  spec <- data.frame(metric = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                     flag = c("qc_low_counts", "qc_low_features", "qc_high_mt_mad"),
                     direction = c("lower", "lower", "higher"),
                     log = c(TRUE, TRUE, FALSE), stringsAsFactors = FALSE)
  if (!all(spec$metric %in% names(d))) stop("Run qc_add_metrics() first.")
  for (flag in spec$flag) d[[flag]] <- NA
  thresholds <- list()
  for (sample in levels(d[[sample_col]])) {
    ix <- which(d[[sample_col]] == sample)
    for (j in seq_len(nrow(spec))) {
      values <- d[[spec$metric[j]]][ix]
      if (!is.numeric(values)) stop("QC metrics must be numeric.")
      valid <- is.finite(values)
      if (spec$log[j]) valid <- valid & values > 0
      lower <- higher <- med <- spread <- NA_real_
      status <- "too_few_evaluable_cells"
      if (sum(valid) >= min_cells) {
        # Use scuttle itself; do not reproduce its threshold implementation.
        flags <- scuttle::isOutlier(values[valid], nmads = nmads,
                                   type = spec$direction[j], log = spec$log[j])
        th <- attr(flags, "thresholds")
        if (is.null(th)) th <- attr(flags, "threshold")
        if (is.null(th) || length(th) != 2L) stop("Unexpected scuttle threshold format.")
        lower <- as.numeric(th[1]); higher <- as.numeric(th[2])
        d[[spec$flag[j]]][ix[valid]] <- as.logical(flags)
        transformed <- if (spec$log[j]) log2(values[valid]) else values[valid]
        med <- stats::median(transformed)
        spread <- stats::mad(transformed)
        status <- if (spread == 0) "zero_MAD_review" else "diagnostic_only"
      }
      # Zero RNA counts are independently flagged, but never discarded here.
      if (spec$log[j]) d[[spec$flag[j]]][ix[is.finite(values) & values == 0]] <- TRUE
      row <- data.frame(sample = sample, metric = spec$metric[j], direction = spec$direction[j],
                        log_transformed = spec$log[j], nmads = nmads, n_used = sum(valid),
                        median_on_fitted_scale = med, mad_on_fitted_scale = spread,
                        lower = lower, higher = higher,
                        threshold = if (spec$direction[j] == "lower") lower else higher,
                        status = status, stringsAsFactors = FALSE)
      names(row)[1] <- sample_col
      thresholds[[length(thresholds) + 1L]] <- row
    }
  }
  d$qc_osca_low_rna <- d$qc_low_counts | d$qc_low_features
  d$qc_osca_any <- d$qc_osca_low_rna | d$qc_high_mt_mad
  list(metrics = d, thresholds = do.call(rbind, thresholds))
}

qc_attach_metadata <- function(object, table, columns = grep("^qc_", names(table), value = TRUE)) {
  if (!"cell_id" %in% names(table) || anyDuplicated(table$cell_id)) stop("Unique cell_id is required.")
  pos <- match(colnames(object), table$cell_id)
  if (anyNA(pos)) stop("Missing cells in metadata table.")
  add <- table[pos, columns, drop = FALSE]
  rownames(add) <- colnames(object)
  SeuratObject::AddMetaData(object, add)
}

qc_summary <- function(x, sample_col = "sample_id") {
  d <- qc_metadata(x, sample_col)
  answer <- lapply(levels(d[[sample_col]]), function(s) {
    z <- d[d[[sample_col]] == s, , drop = FALSE]
    row <- data.frame(sample = s, cells = nrow(z), stringsAsFactors = FALSE)
    names(row)[1] <- sample_col
    for (metric in intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ADT", "nFeature_ADT"), names(z))) {
      v <- z[[metric]][is.finite(z[[metric]])]
      row[[paste0("median_", metric)]] <- if (length(v)) stats::median(v) else NA_real_
    }
    for (flag in grep("^qc_(zero_rna|low_counts|low_features|high_mt_mad|osca_any)$", names(z), value = TRUE)) {
      row[[paste0("n_", flag)]] <- sum(z[[flag]] %in% TRUE)
      row[[paste0("n_unscored_", flag)]] <- sum(is.na(z[[flag]]))
    }
    row
  })
  .qc_bind_rows(answer)
}

.qc_plot_data <- function(x, metrics, sample_col) {
  d <- qc_metadata(x, sample_col)
  for (metric in metrics) {
    if (!metric %in% names(d) || !is.numeric(d[[metric]])) stop("Missing/non-numeric metric: ", metric)
  }
  valid <- Reduce(`&`, lapply(d[metrics], is.finite))
  if (any(!valid)) warning(sum(!valid), " cells have undefined plot values; they remain in the metadata.")
  d <- d[valid, , drop = FALSE]
  if (!nrow(d)) stop("No finite values available for the requested plot.")
  d
}

.qc_plot_thresholds <- function(thresholds, metric, d, sample_col) {
  if (is.null(thresholds)) return(NULL)
  need <- c(sample_col, "metric", "threshold")
  if (!all(need %in% names(thresholds))) stop("Use the thresholds table from qc_flag_osca().")
  th <- thresholds[thresholds$metric == metric & is.finite(thresholds$threshold) &
                     thresholds[[sample_col]] %in% as.character(d[[sample_col]]), , drop = FALSE]
  th[[sample_col]] <- factor(th[[sample_col]], levels = levels(d[[sample_col]]))
  th
}

qc_plot_violin <- function(x, metric = "percent.mt", sample_col = "sample_id",
                           thresholds = NULL, log1p = FALSE) {
  d <- .qc_plot_data(x, metric, sample_col)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[sample_col]], y = .data[[metric]], fill = .data[[sample_col]])) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, alpha = 0.5) +
    ggplot2::geom_boxplot(width = 0.13, outlier.shape = NA, alpha = 0.7) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = metric, title = paste("QC distribution:", metric))
  if (log1p) p <- p + ggplot2::scale_y_continuous(trans = "log1p") +
    ggplot2::labs(y = paste0(metric, " (log1p axis; labels in original units)"))
  th <- .qc_plot_thresholds(thresholds, metric, d, sample_col)
  if (!is.null(th) && nrow(th)) p <- p +
    ggplot2::geom_point(data = th, ggplot2::aes(x = .data[[sample_col]], y = .data$threshold),
                        inherit.aes = FALSE, shape = 95, size = 10) +
    ggplot2::labs(subtitle = "Horizontal marks: sample-specific MAD diagnostic thresholds, not applied filters")
  p
}

qc_plot_histogram <- function(x, metric = "percent.mt", sample_col = "sample_id",
                              thresholds = NULL, bins = 60, log1p = FALSE) {
  d <- .qc_plot_data(x, metric, sample_col)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[metric]])) +
    ggplot2::geom_histogram(bins = bins) +
    ggplot2::facet_wrap(stats::reformulate(sample_col), scales = "free_y") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = metric, y = "Cells", title = paste("QC histogram:", metric))
  if (log1p) p <- p + ggplot2::scale_x_continuous(trans = "log1p") +
    ggplot2::labs(x = paste0(metric, " (log1p axis; labels in original units)"))
  th <- .qc_plot_thresholds(thresholds, metric, d, sample_col)
  if (!is.null(th) && nrow(th)) p <- p +
    ggplot2::geom_vline(data = th, ggplot2::aes(xintercept = .data$threshold), linetype = 2) +
    ggplot2::labs(subtitle = "Dashed lines: sample-specific MAD diagnostic thresholds, not applied filters")
  p
}

qc_plot_scatter <- function(x, x_metric = "nCount_RNA", y_metric = "percent.mt",
                            sample_col = "sample_id", colour_by = NULL,
                            log1p_x = TRUE, log1p_y = FALSE,
                            max_cells_per_sample = 30000, seed = 1) {
  if (length(max_cells_per_sample) != 1L || is.na(max_cells_per_sample) || max_cells_per_sample < 1) {
    stop("max_cells_per_sample must be positive, or Inf.")
  }
  d <- .qc_plot_data(x, c(x_metric, y_metric), sample_col)
  n_all <- nrow(d)
  if (!is.null(colour_by) && !colour_by %in% names(d)) stop("Missing colour_by column.")
  # Preserve the caller's random-number state.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  keep <- unlist(lapply(split(seq_len(nrow(d)), d[[sample_col]], drop = TRUE), function(ix) {
    if (length(ix) > max_cells_per_sample) ix[sample.int(length(ix), as.integer(max_cells_per_sample))] else ix
  }), use.names = FALSE)
  d <- d[keep, , drop = FALSE]
  mapping <- if (is.null(colour_by)) ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]]) else
    ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]], colour = .data[[colour_by]])
  p <- ggplot2::ggplot(d, mapping) + ggplot2::geom_point(size = 0.35, alpha = 0.35) +
    ggplot2::facet_wrap(stats::reformulate(sample_col)) + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = x_metric, y = y_metric, title = paste(y_metric, "versus", x_metric),
                  subtitle = paste(format(nrow(d), big.mark = ","), "of", format(n_all, big.mark = ","),
                                   "evaluable cells plotted; metrics and threshold counts use all cells"))
  if (log1p_x) p <- p + ggplot2::scale_x_continuous(trans = "log1p")
  if (log1p_y) p <- p + ggplot2::scale_y_continuous(trans = "log1p")
  p
}

# Exact hypothetical cell loss at each cutoff. Equality is retained.
# Undefined mt% (e.g. zero-RNA cells) is reported, never treated as zero.
qc_mt_sweep <- function(x, cutoffs = seq(0, 100, 1), sample_col = "sample_id") {
  d <- qc_metadata(x, sample_col)
  if (!"percent.mt" %in% names(d) || !is.numeric(d$percent.mt)) stop("Missing percent.mt.")
  if (!is.numeric(cutoffs) || !length(cutoffs) || any(!is.finite(cutoffs)) || any(cutoffs < 0 | cutoffs > 100)) {
    stop("cutoffs must be finite percentages from 0 through 100.")
  }
  cutoffs <- sort(unique(cutoffs))
  .qc_bind_rows(lapply(levels(d[[sample_col]]), function(s) {
    z <- d$percent.mt[d[[sample_col]] == s]
    valid <- is.finite(z)
    if (any(z[valid] < 0 | z[valid] > 100)) stop("percent.mt values must be within [0,100].")
    n <- sum(valid)
    discarded <- vapply(cutoffs, function(cutoff) sum(z[valid] > cutoff), numeric(1))
    ans <- data.frame(sample = s, cutoff = cutoffs, n_input = length(z),
                      n_evaluable = n, n_undefined_mt = sum(!valid),
                      n_discarded = discarded, n_retained = n - discarded,
                      percent_discarded = if (n > 0) 100 * discarded / n else NA_real_,
                      percent_retained = if (n > 0) 100 * (n - discarded) / n else NA_real_)
    names(ans)[1] <- sample_col
    ans
  }))
}

qc_plot_mt_sweep <- function(sweep, y = c("percent_discarded", "n_discarded"), sample_col = "sample_id") {
  y <- match.arg(y)
  if (!all(c(sample_col, "cutoff", y) %in% names(sweep))) stop("Use qc_mt_sweep() output.")
  sweep[[sample_col]] <- factor(sweep[[sample_col]], levels = unique(sweep[[sample_col]]))
  ggplot2::ggplot(sweep, ggplot2::aes(x = .data$cutoff, y = .data[[y]])) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(stats::reformulate(sample_col)) + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "Candidate maximum mitochondrial percentage",
                  y = if (y == "percent_discarded") "Cells discarded (% of evaluable cells)" else "Cells discarded (number)",
                  title = "Mitochondrial cutoff sensitivity: no filtering applied",
                  subtitle = "Hypothetical discard rule: percent.mt > cutoff. Undefined mt% is reported separately.")
}

# This writes a blank decision template, not inferred/automatic thresholds.
qc_threshold_template <- function(x, sample_col = "sample_id") {
  d <- qc_metadata(x, sample_col)
  out <- data.frame(sample = levels(d[[sample_col]]), min_counts = NA_real_,
                     min_features = NA_real_, max_counts = NA_real_,
                     max_features = NA_real_, max_mt = NA_real_)
  names(out)[1] <- sample_col
  out
}

# NA means "do not apply this gate". Boundary values are retained.
# No MAD flag is used unless the user explicitly copies its value into this table.
qc_manual_decisions <- function(x, thresholds, sample_col = "sample_id") {
  d <- qc_metadata(x, sample_col)
  gates <- c(min_counts = "nCount_RNA", min_features = "nFeature_RNA",
             max_counts = "nCount_RNA", max_features = "nFeature_RNA", max_mt = "percent.mt")
  if (!is.data.frame(thresholds) || !all(c(sample_col, names(gates)) %in% names(thresholds))) {
    stop("Use the complete table from qc_threshold_template().")
  }
  if (anyNA(thresholds[[sample_col]]) || anyDuplicated(thresholds[[sample_col]])) stop("Invalid threshold sample IDs.")
  pos <- match(as.character(d[[sample_col]]), as.character(thresholds[[sample_col]]))
  if (anyNA(pos)) stop("Every sample needs a threshold row; NA gates are allowed.")
  for (gate in names(gates)) {
    value <- thresholds[[gate]]
    if (!is.numeric(value) && !all(is.na(value))) stop("Threshold column must be numeric: ", gate)
    value <- as.numeric(value)
    if (any(!is.na(value) & (!is.finite(value) | value < 0))) stop("Invalid threshold: ", gate)
    if (gate == "max_mt" && any(value > 100, na.rm = TRUE)) stop("max_mt must be <= 100.")
    thresholds[[gate]] <- value
  }
  for (suffix in c("counts", "features")) {
    low <- thresholds[[paste0("min_", suffix)]]; high <- thresholds[[paste0("max_", suffix)]]
    if (any(low > high, na.rm = TRUE)) stop("Minimum exceeds maximum for ", suffix)
  }
  d$qc_manual_missing_metric <- FALSE
  for (gate in names(gates)) {
    if (!gates[[gate]] %in% names(d) || !is.numeric(d[[gates[[gate]]]])) stop("Missing metric: ", gates[[gate]])
    limit <- thresholds[[gate]][pos]
    value <- d[[gates[[gate]]]]
    active <- !is.na(limit)
    missing <- active & !is.finite(value)
    fail <- rep(FALSE, nrow(d))
    valid <- active & is.finite(value)
    fail[valid] <- if (startsWith(gate, "min_")) value[valid] < limit[valid] else value[valid] > limit[valid]
    d[[paste0("qc_manual_fail_", gate)]] <- fail
    d$qc_manual_missing_metric <- d$qc_manual_missing_metric | missing
  }
  reasons <- c(paste0("qc_manual_fail_", names(gates)), "qc_manual_missing_metric")
  d$qc_keep <- !Reduce(`|`, d[reasons])
  d$qc_manual_reason <- apply(d[reasons], 1, function(z) paste(sub("^qc_manual_", "", reasons[as.logical(z)]), collapse = ";"))
  d
}

# The ONLY cell-subsetting function. Call only after reviewing the report.
qc_subset_cells <- function(object, decisions, assays = c("RNA", "ADT")) {
  if (!all(c("cell_id", "qc_keep") %in% names(decisions)) ||
      !is.logical(decisions$qc_keep) || anyNA(decisions$qc_keep) || anyDuplicated(decisions$cell_id)) {
    stop("Use qc_manual_decisions() output with unique cell IDs and defined logical qc_keep.")
  }
  pos <- match(colnames(object), decisions$cell_id)
  if (anyNA(pos)) stop("Decisions are missing cells from this object.")
  keep <- colnames(object)[decisions$qc_keep[pos]]
  if (!length(keep)) stop("These decisions would remove every cell in this object. Review the thresholds.")
  object <- qc_attach_metadata(object, decisions, grep("^qc_", names(decisions), value = TRUE))
  object <- subset(object, cells = keep)
  object <- qc_drop_zero_features(object, assays)
  if (!is.null(object@misc$qc_review_settings)) {
    object@misc$qc_review_settings$no_cell_filters_applied <- FALSE
  }
  object@misc$qc_manual_filter <- list(input_cells = length(pos), retained_cells = length(keep),
                                      filtered_at = as.character(Sys.time()))
  object
}

qc_save_plots <- function(plots, out_dir, width = 10, height = 6, dpi = 160) {
  if (is.null(names(plots)) || any(!nzchar(names(plots))) || anyDuplicated(names(plots))) stop("Name each plot uniquely.")
  if (any(!grepl("^[A-Za-z0-9_.-]+$", names(plots)))) stop("Plot names must be filename-safe.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(plots)) {
    ggplot2::ggsave(filename = file.path(out_dir, paste0(name, ".png")),
                    plot = plots[[name]], width = width, height = height, units = "in", dpi = dpi)
  }
  invisible(file.path(out_dir, paste0(names(plots), ".png")))
}
