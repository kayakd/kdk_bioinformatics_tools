# Purpose: Import Cell Ranger RNA + ADT counts, remove all-zero features,
#          calculate QC metrics, and save one pre-QC Seurat object.
# Input: samples.csv and each sample's filtered_feature_bc_matrix.h5.
# Output: qc_raw/combined_raw.rds and qc_raw/cell_qc_metrics.csv.
# Author: Dr. Koray Doğan Kaya

library(Seurat)
library(Matrix)
options(Seurat.object.assay.version = "v5")

cellranger_base_dir <- "/EDIT/absolute/path/to/cellranger_results"
output_dir <- "qc_raw"
mito_pattern <- "^mt-"  # Mouse gene symbols; use ^MT- for human symbols.

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
samples <- read.csv("samples.csv", stringsAsFactors = FALSE)

# Relative paths are resolved against the alignment output base directory.
input_paths <- path.expand(samples$matrix_path)
relative <- !grepl("^(/|[A-Za-z]:[/\\\\])", input_paths)
input_paths[relative] <- file.path(cellranger_base_dir, input_paths[relative])

objects <- lapply(seq_len(nrow(samples)), function(i) {
  id <- samples$sample_id[i]
  counts <- Read10X_h5(input_paths[i])
  rna <- counts[["Gene Expression"]]
  adt <- counts[["Antibody Capture"]][, colnames(rna), drop = FALSE]

  # Calculate mt% from RNA only, before removing feature rows.
  mito <- grepl(mito_pattern, rownames(rna))
  percent_mt <- 100 * Matrix::colSums(rna[mito, , drop = FALSE]) /
    Matrix::colSums(rna)
  percent_mt[!is.finite(percent_mt) | !any(mito)] <- NA_real_

  # Only all-zero rows are removed. Every input barcode is retained.
  keep_rna <- Matrix::rowSums(rna) > 0
  keep_adt <- Matrix::rowSums(adt) > 0

  metadata <- samples[
    rep(i, ncol(rna)), setdiff(names(samples), "matrix_path"), drop = FALSE
  ]
  rownames(metadata) <- colnames(rna)
  metadata$barcode_10x <- colnames(rna)
  metadata$percent.mt <- unname(percent_mt)

  object <- CreateSeuratObject(
    counts = rna[keep_rna, , drop = FALSE],
    meta.data = metadata, project = id,
    min.cells = 0, min.features = 0
  )
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(
    counts = adt[keep_adt, , drop = FALSE],
    min.cells = 0, min.features = 0
  )
  object@misc$zero_features_removed <- list(
    RNA = rownames(rna)[!keep_rna], ADT = rownames(adt)[!keep_adt]
  )
  object@misc$mitochondrial_features <- rownames(rna)[mito]
  RenameCells(object, add.cell.id = id)
})
names(objects) <- samples$sample_id

# One combined object, with distinct RNA and ADT assays.
raw <- merge(objects[[1]], y = objects[-1], merge.data = FALSE,
             project = "GARP_CITEseq")
raw@misc$qc_import <- list(
  sample_manifest = samples,
  matrix_paths = input_paths,
  zero_features_removed = lapply(objects, function(x) x@misc$zero_features_removed),
  mitochondrial_features = lapply(objects, function(x) x@misc$mitochondrial_features),
  cell_filters_applied = FALSE,
  normalization_applied = FALSE,
  seurat_version = as.character(packageVersion("Seurat"))
)
rm(objects)
invisible(gc())

# Store one raw counts layer per assay, not separate sample layers.
raw <- JoinLayers(raw, assay = "RNA")
raw <- JoinLayers(raw, assay = "ADT")
raw$sample_id <- factor(raw$sample_id, levels = samples$sample_id)
Idents(raw) <- "sample_id"
DefaultAssay(raw) <- "RNA"

saveRDS(raw, file.path(output_dir, "combined_raw.rds"), compress = FALSE)
qc <- raw[[]]
qc$cell_id <- rownames(qc)
write.csv(qc, file.path(output_dir, "cell_qc_metrics.csv"), row.names = FALSE)

raw
