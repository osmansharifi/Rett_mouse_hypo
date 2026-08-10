#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

# Gene-level feature filtering for pseudobulk DEG.
# This script starts from pseudobulk Seurat objects that have already passed
# pseudobulk sample-level QC in pseudobulking_qc_pbsample.R.

DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/sample_level_view/qc"
FILE_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg/sample_level_view"
SAMPLE_QC_DIR <- file.path(FILE_DIR, "qc")

CLUSTER_RESOLUTIONS <- c("neuron_class", "cell_type_concise")

CPM_CUTOFF <- 1
MIN_PROP_SAMPLES <- 0.25
PROP_LEVELS_TO_PLOT <- c(0.25, 0.50)

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FILE_DIR, recursive = TRUE, showWarnings = FALSE)

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unknown", x)
}

get_counts_cpm <- function(counts) {
  lib_sizes <- Matrix::colSums(counts)
  keep_samples <- lib_sizes > 0
  counts <- counts[, keep_samples, drop = FALSE]
  lib_sizes <- lib_sizes[keep_samples]

  if (ncol(counts) == 0) {
    stop("No pseudobulk samples have non-zero library sizes.")
  }

  t(t(counts) / lib_sizes * 1e6)
}

get_feature_filter <- function(counts,
                               cpm_cutoff = CPM_CUTOFF,
                               min_prop_samples = MIN_PROP_SAMPLES) {
  cpm <- get_counts_cpm(counts)
  n_samples_cpm_ge_cutoff <- Matrix::rowSums(cpm >= cpm_cutoff)
  min_samples <- max(1, ceiling(ncol(cpm) * min_prop_samples))

  data.frame(
    gene = rownames(cpm),
    n_pseudobulk_samples = ncol(cpm),
    n_samples_cpm_ge_cutoff = n_samples_cpm_ge_cutoff,
    prop_samples_cpm_ge_cutoff = n_samples_cpm_ge_cutoff / ncol(cpm),
    cpm_cutoff = cpm_cutoff,
    min_prop_samples = min_prop_samples,
    min_samples = min_samples,
    keep_feature = n_samples_cpm_ge_cutoff >= min_samples,
    stringsAsFactors = FALSE
  )
}

get_feature_plot_df <- function(seurat_pb,
                                group_col,
                                cpm_cutoff = CPM_CUTOFF,
                                prop_levels = PROP_LEVELS_TO_PLOT) {
  group_levels <- sort(unique(as.character(seurat_pb[[group_col]][, 1])))

  bind_rows(lapply(group_levels, function(group_level) {
    cells_use <- colnames(seurat_pb)[as.character(seurat_pb[[group_col]][, 1]) == group_level]
    seurat_subset <- subset(seurat_pb, cells = cells_use)
    counts <- LayerData(seurat_subset, assay = "RNA", layer = "counts")
    cpm <- get_counts_cpm(counts)
    n_samples_cpm_ge_cutoff <- Matrix::rowSums(cpm >= cpm_cutoff)

    bind_rows(lapply(prop_levels, function(prop_level) {
      min_samples <- max(1, ceiling(ncol(cpm) * prop_level))
      data.frame(
        cluster_resolution = group_col,
        feature_subset = group_level,
        prop_label = paste0(prop_level * 100, "%"),
        cpm_cutoff = cpm_cutoff,
        min_prop_samples = prop_level,
        min_samples = min_samples,
        num_highly_expressed_genes = sum(n_samples_cpm_ge_cutoff >= min_samples),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

plot_feature_filter_summary <- function(plot_df, group_col) {
  plot_df$feature_subset <- factor(
    plot_df$feature_subset,
    levels = unique(plot_df$feature_subset)
  )
  plot_df$prop_label <- factor(
    plot_df$prop_label,
    levels = paste0(PROP_LEVELS_TO_PLOT * 100, "%")
  )

  ggplot(plot_df, aes(x = feature_subset, y = num_highly_expressed_genes)) +
    geom_col(fill = "grey35", width = 0.82) +
    facet_wrap(vars(prop_label), nrow = 1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = group_col,
      y = "num_highly_expressed_genes",
      title = "Genes with at least 1 CPM in 25% or 50% of pseudobulk samples"
    )
}

save_feature_filtered_subset <- function(seurat_pb,
                                         group_col,
                                         group_level,
                                         cpm_cutoff = CPM_CUTOFF,
                                         min_prop_samples = MIN_PROP_SAMPLES,
                                         object_dir = DATA_DIR,
                                         table_dir = FILE_DIR) {
  group_safe <- safe_filename(group_level)
  cells_use <- colnames(seurat_pb)[as.character(seurat_pb[[group_col]][, 1]) == group_level]
  seurat_subset <- subset(seurat_pb, cells = cells_use)
  counts <- LayerData(seurat_subset, assay = "RNA", layer = "counts")

  gene_filter <- get_feature_filter(
    counts = counts,
    cpm_cutoff = cpm_cutoff,
    min_prop_samples = min_prop_samples
  )

  keep_genes <- gene_filter$gene[gene_filter$keep_feature]
  if (length(keep_genes) == 0) {
    stop(
      "No genes passed feature filtering for ", group_col, " / ", group_level,
      ". Consider lowering CPM_CUTOFF or MIN_PROP_SAMPLES."
    )
  }

  seurat_subset_feature <- subset(seurat_subset, features = keep_genes)
  counts_feature <- LayerData(seurat_subset_feature, assay = "RNA", layer = "counts")
  meta_feature <- seurat_subset_feature[[]]

  prefix <- paste0("pseudobulk_", group_col, "_", group_safe, "_feature")

  saveRDS(
    seurat_subset_feature,
    file.path(object_dir, paste0(prefix, ".rds"))
  )
  write.csv(
    as.matrix(counts_feature),
    file.path(table_dir, paste0(prefix, "_counts.csv"))
  )
  write.csv(
    meta_feature,
    file.path(table_dir, paste0(prefix, "_metadata.csv")),
    row.names = TRUE
  )

  gene_filter[[group_col]] <- group_level
  gene_filter$n_features_before <- nrow(counts)
  gene_filter$n_features_after <- length(keep_genes)

  message(
    "  ", group_level, ": kept ", length(keep_genes), " / ", nrow(counts),
    " genes across ", ncol(seurat_subset), " pseudobulk samples"
  )

  gene_filter
}

process_cluster_resolution <- function(group_col) {
  message("Processing feature filter for: ", group_col)

  object_path <- file.path(SAMPLE_QC_DIR, paste0("pseudobulk_seurat_", group_col, ".rds"))
  if (!file.exists(object_path)) {
    stop("Sample-filtered pseudobulk Seurat object not found: ", object_path)
  }

  seurat_pb <- readRDS(object_path)
  if (!group_col %in% colnames(seurat_pb[[]])) {
    stop("Metadata column not found in pseudobulk object: ", group_col)
  }

  group_levels <- sort(unique(as.character(seurat_pb[[group_col]][, 1])))

  plot_df <- get_feature_plot_df(seurat_pb = seurat_pb, group_col = group_col)
  write.csv(
    plot_df,
    file.path(FILE_DIR, paste0("pb_feature_", group_col, ".csv")),
    row.names = FALSE
  )
  ggsave(
    file.path(FIG_DIR, paste0("qc_feature_", group_col, ".pdf")),
    plot = plot_feature_filter_summary(plot_df, group_col = group_col),
    height = 5,
    width = 10
  )

  gene_filter_list <- lapply(
    group_levels,
    function(group_level) {
      save_feature_filtered_subset(
        seurat_pb = seurat_pb,
        group_col = group_col,
        group_level = group_level
      )
    }
  )

  gene_filter_df <- bind_rows(gene_filter_list)
  names(gene_filter_df)[names(gene_filter_df) == group_col] <- "feature_subset"
  gene_filter_df$cluster_resolution <- group_col
  gene_filter_df$feature_subset <- factor(gene_filter_df$feature_subset, levels = group_levels)

  gene_filter_csv <- file.path(FILE_DIR, paste0("pseudobulk_gene_filter_", group_col, "_feature.csv"))
  write.csv(gene_filter_df, gene_filter_csv, row.names = FALSE)

  summary_df <- gene_filter_df %>%
    group_by(cluster_resolution, feature_subset) %>%
    summarise(
      n_pseudobulk_samples = first(n_pseudobulk_samples),
      cpm_cutoff = first(cpm_cutoff),
      min_prop_samples = first(min_prop_samples),
      min_samples = first(min_samples),
      n_features_before = first(n_features_before),
      n_features_after = first(n_features_after),
      .groups = "drop"
    )

  summary_csv <- file.path(FILE_DIR, paste0("pseudobulk_gene_filter_summary_", group_col, "_feature.csv"))
  write.csv(summary_df, summary_csv, row.names = FALSE)

  union_keep_genes <- unique(gene_filter_df$gene[gene_filter_df$keep_feature])
  if (length(union_keep_genes) == 0) {
    stop(
      "No genes passed feature filtering for ", group_col,
      ". Consider lowering CPM_CUTOFF or MIN_PROP_SAMPLES."
    )
  }

  seurat_pb_feature <- subset(seurat_pb, features = union_keep_genes)
  counts_feature <- LayerData(seurat_pb_feature, assay = "RNA", layer = "counts")

  saveRDS(
    seurat_pb_feature,
    file.path(DATA_DIR, paste0("pseudobulk_seurat_", group_col, "_feature.rds"))
  )
  write.csv(
    as.matrix(counts_feature),
    file.path(FILE_DIR, paste0("pseudobulk_counts_", group_col, "_feature.csv"))
  )
  write.csv(
    seurat_pb_feature[[]],
    file.path(FILE_DIR, paste0("pseudobulk_metadata_", group_col, "_feature.csv")),
    row.names = TRUE
  )

  message("Saved feature-filter gene table: ", gene_filter_csv)
  message("Saved feature-filter summary: ", summary_csv)
  message(
    "Saved union feature-filtered ", group_col, " object with ",
    length(union_keep_genes), " genes"
  )

  invisible(list(
    seurat = seurat_pb_feature,
    gene_filter = gene_filter_df,
    summary = summary_df,
    eda = plot_df
  ))
}

invisible(lapply(CLUSTER_RESOLUTIONS, process_cluster_resolution))
