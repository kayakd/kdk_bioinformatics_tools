# Purpose:
#   Produce the initial RNA and CITE-seq QC figures from unnormalized data.
#   Plotting never filters cells, alters counts, or selects thresholds.
# Input:
#   A Seurat object with sample_id, nCount_RNA, nFeature_RNA, percent.mt,
#   and nCount_ADT in cell metadata. The sample column is configurable.
# Output:
#   A patchwork distribution figure and reusable ggplot scatter figures.
# Author: Dr. Koray Doğan Kaya
# Dependencies: ggplot2, patchwork, scales.

plot_qc_distributions <- function(
    object, sample_col = "sample_id",
    metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ADT")) {
  d <- object[[]]
  d$.sample <- factor(as.character(d[[sample_col]]),
                      levels = unique(as.character(d[[sample_col]])))

  # Counts use log1p axes; mitochondrial percentage stays linear.
  specifications <- list(
    nCount_RNA = list(label = "RNA UMIs", transform = "log1p",
                     breaks = c(0, 100, 1000, 10000, 100000, 1000000)),
    nFeature_RNA = list(label = "Detected RNA genes", transform = "log1p",
                       breaks = c(0, 100, 500, 1000, 3000, 10000)),
    percent.mt = list(label = "Mitochondrial RNA (%)", transform = "identity",
                      breaks = c(0, 25, 50, 75, 100)),
    nCount_ADT = list(label = "Antibody UMIs", transform = "log1p",
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000))
  )

  metrics <- intersect(metrics, names(d))
  panels <- lapply(metrics, function(metric) {
    specification <- specifications[[metric]]
    ggplot2::ggplot(d, ggplot2::aes(x = .sample, y = .data[[metric]])) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, na.rm = TRUE) +
      ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, na.rm = TRUE) +
      ggplot2::scale_y_continuous(
        trans = specification$transform,
        breaks = specification$breaks,
        labels = scales::label_number(scale_cut = scales::cut_short_scale())
      ) +
      ggplot2::labs(x = NULL, y = specification$label) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  })

  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "QC distributions by sample",
      subtitle = "No cell filtering. Count axes use log1p; tick labels show original counts.",
      caption = paste0(
        "Cells: ", format(nrow(d), big.mark = ",", scientific = FALSE),
        ". Undefined mt%: ", sum(!is.finite(d$percent.mt)),
        " (retained in the object, omitted only from the mt% distribution)."
      )
    )
}

plot_qc_relationship <- function(object, y = "nFeature_RNA",
                                 sample_col = "sample_id", ncol = 3) {
  d <- object[[]]
  d$.sample <- factor(as.character(d[[sample_col]]),
                      levels = unique(as.character(d[[sample_col]])))

  # Draw high-mitochondrial cells last so they are not hidden by overplotting.
  d <- d[order(d$percent.mt, na.last = FALSE), , drop = FALSE]
  y_labels <- c(nFeature_RNA = "Detected RNA genes", nCount_ADT = "Antibody UMIs")

  ggplot2::ggplot(
    d, ggplot2::aes(x = nCount_RNA, y = .data[[y]], colour = percent.mt)
  ) +
    ggplot2::geom_point(size = 0.25, alpha = 0.5) +
    ggplot2::facet_wrap(~ .sample, ncol = ncol) +
    ggplot2::scale_x_continuous(
      trans = "log1p", breaks = c(0, 100, 1000, 10000, 100000, 1000000),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::scale_y_continuous(
      trans = "log1p", breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
      labels = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::scale_colour_continuous(
      name = "Mitochondrial\nRNA (%)", limits = c(0, 100)
    ) +
    ggplot2::labs(
      x = "RNA UMIs", y = unname(y_labels[y]),
      title = paste(unname(y_labels[y]), "versus RNA UMIs"),
      subtitle = "All cells plotted; no subsampling. Log1p axes, original-count labels.",
      caption = "Colour represents mitochondrial RNA percentage, not a filtering decision."
    ) +
    ggplot2::theme_bw(base_size = 11)
}
