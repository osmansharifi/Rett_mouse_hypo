#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
})

# Set path and load data
DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/sample_level_view/qc"
FILE_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg/sample_level_view/qc"
seu <- readRDS(file.path(DATA_DIR, "seu_final.rds"))

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FILE_DIR, recursive = TRUE, showWarnings = FALSE)

# Configuration
MIN_CELLS  <- 10
MIN_COUNTS <- 500
N_PCS_ASSOC <- 5
N_PCS_VAR   <- 30

# Columns to test against pseudobulk PCs when present in the pseudobulk metadata.
BASE_PCA_COVARIATES <- c(
  "sample",
  "label",
  "replicate",
  "animal_genotype",
  "animal_sex",
  "animal_timepoint",
  "psbulk_cells",
  "psbulk_cells_log",
  "psbulk_counts",
  "psbulk_counts_log"
)

clean_pb_component <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^[0-9]", x), paste0("g", x), x)
  gsub("_", "-", x)
}

make_pb_id <- function(sample, group) {
  paste(clean_pb_component(sample), clean_pb_component(group), sep = "_")
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unknown", x)
}

# 1. Generate pseduobulk samples --------------------------------------
aggregate_pseudobulk_object <- function(seu_obj, sample_col = "sample", group_col) {
  AggregateExpression(
    seu_obj,
    assays = "RNA",
    return.seurat = TRUE,
    group.by = c(sample_col, group_col)
  )
}

get_sample_metadata <- function(seu_obj, sample_col = "sample") {
  meta <- seu_obj[[]]
  sample_cols <- c(
    sample_col,
    intersect(BASE_PCA_COVARIATES, colnames(meta)),
    grep("^animal_", colnames(meta), value = TRUE)
  )
  sample_cols <- unique(sample_cols)

  meta %>%
    tibble::rownames_to_column("cell_id") %>%
    select(all_of(sample_cols)) %>%
    group_by(.data[[sample_col]]) %>%
    summarise(
      across(everything(), ~ {
        vals <- unique(na.omit(as.character(.x)))
        if (length(vals) == 0) NA_character_ else vals[1]
      }),
      .groups = "drop"
    )
}

# 2. Generate pseduobulk QC metrics --------------------------------------
get_pseudobulk_qc <- function(seu_obj, seurat_pb, sample_col = "sample", group_col) {
  pb_counts <- LayerData(seurat_pb, assay = "RNA", layer = "counts")

  cell_counts <- as.data.frame(table(seu_obj[[sample_col]][, 1], seu_obj[[group_col]][, 1])) %>%
    filter(Freq > 0) %>%
    mutate(pb_id = make_pb_id(Var1, Var2)) %>%
    rename(
      !!sample_col := Var1,
      !!group_col := Var2,
      psbulk_cells = Freq
    )

  sample_meta <- get_sample_metadata(seu_obj, sample_col = sample_col)

  data.frame(
    pb_id = colnames(pb_counts),
    psbulk_counts = Matrix::colSums(pb_counts),
    stringsAsFactors = FALSE
  ) %>%
    left_join(cell_counts, by = "pb_id") %>%
    left_join(sample_meta, by = sample_col) %>%
    mutate(
      psbulk_cells_log = log(psbulk_cells + 1),
      psbulk_counts_log = log(psbulk_counts + 1)
    )
}

add_pseudobulk_qc_metadata <- function(seurat_pb, qc_df, group_col) {
  qc_df$pb_id <- as.character(qc_df$pb_id)
  qc_df <- qc_df[match(colnames(seurat_pb), qc_df$pb_id), ]
  stopifnot(all(colnames(seurat_pb) == qc_df$pb_id))

  qc_df <- qc_df %>%
    mutate(across(where(is.character), as.factor))

  meta <- as.data.frame(qc_df)
  rownames(meta) <- as.character(meta$pb_id)

  seurat_pb <- AddMetaData(seurat_pb, metadata = meta)
  SetIdent(seurat_pb, value = group_col)
}

# 3. Visualize pseduobulk QC metrics --------------------------------------
plot_pseudobulk_qc <- function(qc_df,
                               group_col,
                               min_cells = MIN_CELLS,
                               min_counts = MIN_COUNTS,
                               facet = TRUE) {
  p <- ggplot(qc_df, aes(x = log10(psbulk_cells), y = log10(psbulk_counts))) +
    geom_point(aes(color = .data[[group_col]]), size = 2.5, alpha = 0.8) +
    geom_vline(xintercept = log10(min_cells), linetype = "dashed", color = "gray30") +
    geom_hline(yintercept = log10(min_counts), linetype = "dashed", color = "gray30") +
    theme_bw() +
    labs(
      x = expression(log[10] ~ "(cells)"),
      y = expression(log[10] ~ "(counts)"),
      title = paste("Pseudobulk QC:", group_col)
    )

  if (facet) {
    p <- p + facet_wrap(vars(.data[[group_col]])) + theme(legend.position = "none")
  }
  p
}

# 4. Filter pseduobulk QC samples --------------------------------------
filter_pseudobulk_object <- function(seurat_pb,
                                     group_col,
                                     min_cells = MIN_CELLS,
                                     min_counts = MIN_COUNTS,
                                     output_dir = ".") {
  seurat_pb$pass_qc <- seurat_pb$psbulk_cells >= min_cells &
    seurat_pb$psbulk_counts >= min_counts

  valid_ids <- colnames(seurat_pb)[seurat_pb$pass_qc]
  seurat_pb_filtered <- subset(seurat_pb, cells = valid_ids)
  qc_df <- seurat_pb[[]]

  detailed_log <- qc_df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n_samples_before = n(),
      n_samples_passed = sum(pass_qc),
      n_samples_filtered = sum(!pass_qc),
      min_cells_cutoff = min_cells,
      min_counts_cutoff = min_counts,
      .groups = "drop"
    )

  total_row <- data.frame(
    group_temp = "TOTAL_ALL_GROUPS",
    n_samples_before = ncol(seurat_pb),
    n_samples_passed = ncol(seurat_pb_filtered),
    n_samples_filtered = ncol(seurat_pb) - ncol(seurat_pb_filtered),
    min_cells_cutoff = min_cells,
    min_counts_cutoff = min_counts
  )
  colnames(total_row)[1] <- group_col

  detailed_log <- bind_rows(detailed_log, total_row)
  csv_filename <- file.path(output_dir, paste0("pseudobulk_filter_", group_col, ".csv"))
  write.csv(detailed_log, file = csv_filename, row.names = FALSE)
  message("Detailed pseudobulk filtering log saved to: ", csv_filename)

  seurat_pb_filtered
}

# 5. normalize and scale raw counts before computing PCA -------------------------------------- # nolint
preprocess_pseudobulk_pca <- function(seurat_pb, n_pcs = N_PCS_VAR) {
  DefaultAssay(seurat_pb) <- "RNA"
  seurat_pb <- FindVariableFeatures(seurat_pb)
  seurat_pb <- RunPCA(seurat_pb, npcs = n_pcs, verbose = FALSE)
  seurat_pb
}

# 6. Check PCA association with metadata --------------------------------------
check_metadata_pca <- function(object,
                               reduction = "pca",
                               n_pcs = N_PCS_ASSOC,
                               covariates) {
  pcs <- Embeddings(object, reduction = reduction)
  n_pcs <- min(n_pcs, ncol(pcs))
  pcs <- pcs[, seq_len(n_pcs), drop = FALSE]
  meta <- object[[]]
  covariates <- intersect(covariates, colnames(meta))

  assoc <- expand.grid(
    covariate = covariates,
    PC = colnames(pcs),
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      p_value = {
        y <- pcs[, PC]
        x <- meta[[covariate]]
        keep <- !is.na(y) & !is.na(x)
        y <- y[keep]
        x <- x[keep]
        if (length(unique(x)) < 2 || length(y) < 3) {
          NA_real_
        } else if (is.numeric(x)) {
          summary(lm(y ~ x))$coefficients[2, 4]
        } else {
          fit <- anova(lm(y ~ factor(x)))
          fit[1, "Pr(>F)"]
        }
      },
      score = {
        y <- pcs[, PC]
        x <- meta[[covariate]]
        keep <- !is.na(y) & !is.na(x)
        y <- y[keep]
        x <- x[keep]
        if (length(unique(x)) < 2 || length(y) < 3) {
          NA_real_
        } else if (is.numeric(x)) {
          suppressWarnings(cor(y, x, method = "spearman"))
        } else {
          fit <- anova(lm(y ~ factor(x)))
          as.numeric(fit[1, "F value"])
        }
      },
      variable_type = ifelse(is.numeric(meta[[covariate]]), "numeric", "categorical")
    ) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      neg_log10_padj = -log10(pmax(p_adj, .Machine$double.xmin))
    )

  assoc
}

plot_metadata_pca <- function(pca_assoc, padj_cutoff = 0.05) {
  sig_assoc <- pca_assoc %>% filter(!is.na(p_adj), p_adj < padj_cutoff)

  p <- ggplot(pca_assoc, aes(x = PC, y = covariate, fill = neg_log10_padj)) +
    geom_tile(color = "white") +
    scale_fill_gradient(
      low = "#F7F7F7",
      high = "#4B0082",
      na.value = "grey90",
      name = expression(-log[10] ~ "(padj)")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = "Principal component",
      y = "Metadata covariate",
      title = "Adjusted p-values from PCA metadata association tests"
    )

  if (nrow(sig_assoc) > 0) {
    p <- p +
      geom_point(
        data = sig_assoc,
        aes(x = PC, y = covariate, shape = "padj < 0.05", color = "padj < 0.05"),
        inherit.aes = FALSE,
        size = 2.1,
        stroke = 0.8
      ) +
      scale_shape_manual(name = NULL, values = c("padj < 0.05" = 8)) +
      scale_color_manual(name = NULL, values = c("padj < 0.05" = "red")) +
      guides(
        shape = guide_legend(override.aes = list(color = "red", size = 3)),
        color = "none"
      )
  }

  p
}

plot_pca_variance_rank <- function(object, reduction = "pca", n_pcs = N_PCS_VAR) {
  stdev <- object[[reduction]]@stdev
  n_pcs <- min(n_pcs, length(stdev))
  pct <- stdev^2 / sum(stdev^2)

  data.frame(
    ranking = seq_len(n_pcs),
    PC = paste0("PC", seq_len(n_pcs)),
    variance = pct[seq_len(n_pcs)]
  ) %>%
    ggplot(aes(x = ranking, y = variance, label = PC)) +
    geom_col(width = 0.7, fill = "grey75") +
    geom_text(angle = 90, hjust = -0.1, size = 3) +
    theme_bw() +
    labs(x = "ranking", y = "variance explained", title = "PCA variance explained")
}

save_pca_covariate_plots <- function(object,
                                     covariates,
                                     out_path,
                                     reduction = "pca") {
  coords <- as.data.frame(Embeddings(object, reduction = reduction)[, 1:2, drop = FALSE])
  colnames(coords) <- c("PC1", "PC2")
  meta <- object[[]]
  covariates <- intersect(covariates, colnames(meta))

  pdf(out_path, width = 10, height = 8)
  for (covariate in covariates) {
    plot_df <- bind_cols(coords, value = meta[[covariate]])
    p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = value)) +
      geom_point(size = 2.2, alpha = 0.9) +
      theme_bw() +
      labs(title = covariate, color = covariate)
    if (is.numeric(plot_df$value)) {
      p <- p + scale_color_viridis_c(option = "viridis", direction = 1)
    } else {
      plot_df$value <- as.factor(plot_df$value)
      p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = value)) +
        geom_point(size = 2.2, alpha = 0.9) +
        theme_bw() +
        labs(title = covariate, color = covariate)
    }
    print(p)
  }
  dev.off()
  message("PCA covariate plots saved to: ", out_path)
}

# 7. Export pseudobulk assay and metadata --------------------------------------
export_pseudobulk_assay <- function(seurat_pb, group_col, output_dir = ".") {
  counts <- LayerData(seurat_pb, assay = "RNA", layer = "counts")
  metadata <- seurat_pb[[]]

  counts_rds <- file.path(output_dir, paste0("pseudobulk_counts_", group_col, ".rds"))
  counts_csv <- file.path(output_dir, paste0("pseudobulk_counts_", group_col, ".csv"))
  meta_csv <- file.path(output_dir, paste0("pseudobulk_metadata_", group_col, ".csv"))
  object_rds <- file.path(output_dir, paste0("pseudobulk_seurat_", group_col, ".rds"))

  saveRDS(counts, counts_rds)
  write.csv(as.matrix(counts), counts_csv)
  write.csv(metadata, meta_csv, row.names = TRUE)
  saveRDS(seurat_pb, object_rds)

  message("Filtered pseudobulk counts saved to: ", counts_rds)
  message("Filtered pseudobulk counts CSV saved to: ", counts_csv)
  message("Filtered pseudobulk metadata saved to: ", meta_csv)
  message("Filtered pseudobulk Seurat object saved to: ", object_rds)
}


# 8. Run pseudobulk exploratory analysis --------------------------------------
run_pseudobulk_exploratory_analysis <- function(seurat_pb,
                                                group_col,
                                                fig_dir = ".",
                                                n_pcs_assoc = N_PCS_ASSOC,
                                                n_pcs_var = N_PCS_VAR,
                                                covariates = BASE_PCA_COVARIATES) {
  seurat_pb <- preprocess_pseudobulk_pca(seurat_pb, n_pcs = n_pcs_var)

  covariates <- unique(c(group_col, covariates))
  covariates <- intersect(covariates, colnames(seurat_pb[[]]))

  pca_assoc <- check_metadata_pca(
    object = seurat_pb,
    reduction = "pca",
    n_pcs = n_pcs_assoc,
    covariates = covariates
  )

  assoc_csv <- file.path(fig_dir, paste0("pseudobulk_pca_metadata_assoc_", group_col, ".csv"))
  write.csv(pca_assoc, assoc_csv, row.names = FALSE)

  p <- plot_metadata_pca(pca_assoc)
  ggsave(
    file.path(fig_dir, paste0("pseudobulk_pca_metadata_assoc_", group_col, ".pdf")),
    plot = p,
    height = 5.5,
    width = 8
  )

  p <- plot_pca_variance_rank(seurat_pb, n_pcs = n_pcs_var)
  ggsave(
    file.path(fig_dir, paste0("pseudobulk_pca_variance_", group_col, ".pdf")),
    plot = p,
    height = 4,
    width = 7
  )

  save_pca_covariate_plots(
    object = seurat_pb,
    covariates = covariates,
    out_path = file.path(fig_dir, paste0("pseudobulk_pca_covariates_", group_col, ".pdf"))
  )

  seurat_pb
}

# 9. Rerun PCA/covariate association within each filtered cell-type subset -------
# If all-data pseudobulk PCA is strongly associated with cell type, subset-level
# PCA helps check whether technical or biological covariate effects persist after
# removing between-cell-type variation.
run_celltype_subset_exploratory_analysis <- function(seurat_pb,
                                                     group_col,
                                                     fig_dir = ".",
                                                     n_pcs_assoc = N_PCS_ASSOC,
                                                     n_pcs_var = N_PCS_VAR,
                                                     covariates = BASE_PCA_COVARIATES) {
  subset_base_dir <- file.path(fig_dir, group_col)
  dir.create(subset_base_dir, recursive = TRUE, showWarnings = FALSE)

  subset_levels <- sort(unique(as.character(seurat_pb[[group_col]][, 1])))

  for (subset_level in subset_levels) {
    subset_safe <- safe_filename(subset_level)
    subset_dir <- file.path(subset_base_dir, subset_safe)
    dir.create(subset_dir, recursive = TRUE, showWarnings = FALSE)

    subset_cells <- colnames(seurat_pb)[as.character(seurat_pb[[group_col]][, 1]) == subset_level]
    seurat_subset <- subset(seurat_pb, cells = subset_cells)
    message("  Subset PCA: ", group_col, " / ", subset_level, " (", ncol(seurat_subset), " pseudobulk samples)")

    n_pcs_subset <- min(n_pcs_var, ncol(seurat_subset) - 1, nrow(seurat_subset) - 1)
    if (n_pcs_subset < 2) {
      message("    Skipped: fewer than 3 pseudobulk samples or features.")
      next
    }

    seurat_subset <- NormalizeData(
      seurat_subset,
      normalization.method = "LogNormalize",
      scale.factor = 1e6,
      verbose = FALSE
    )
    seurat_subset <- FindVariableFeatures(seurat_subset, verbose = FALSE)
    seurat_subset <- ScaleData(seurat_subset, verbose = FALSE)
    seurat_subset <- RunPCA(seurat_subset, npcs = n_pcs_subset, verbose = FALSE)

    subset_covariates <- setdiff(covariates, group_col)
    subset_covariates <- intersect(subset_covariates, colnames(seurat_subset[[]]))

    pca_assoc <- check_metadata_pca(
      object = seurat_subset,
      reduction = "pca",
      n_pcs = min(n_pcs_assoc, n_pcs_subset),
      covariates = subset_covariates
    )

    write.csv(
      pca_assoc,
      file.path(subset_dir, paste0("pca_metadata_assoc_", subset_safe, ".csv")),
      row.names = FALSE
    )

    p <- plot_metadata_pca(pca_assoc)
    p <- p + labs(title = paste(group_col, subset_level, "PCA metadata association", sep = " - "))
    ggsave(
      file.path(subset_dir, paste0("pca_metadata_assoc_", subset_safe, ".pdf")),
      plot = p,
      height = 5.5,
      width = 8
    )

    p <- plot_pca_variance_rank(seurat_subset, n_pcs = n_pcs_var)
    ggsave(
      file.path(subset_dir, paste0("pca_variance_", subset_safe, ".pdf")),
      plot = p,
      height = 4,
      width = 7
    )

    save_pca_covariate_plots(
      object = seurat_subset,
      covariates = subset_covariates,
      out_path = file.path(subset_dir, paste0("pca_covariates_", subset_safe, ".pdf"))
    )
  }

  invisible(NULL)
}

process_pseudobulk_level <- function(seu_obj,
                                     sample_col = "sample",
                                     group_col,
                                     fig_dir = FIG_DIR,
                                     output_dir = FILE_DIR,
                                     min_cells = MIN_CELLS,
                                     min_counts = MIN_COUNTS,
                                     n_pcs_assoc = N_PCS_ASSOC,
                                     n_pcs_var = N_PCS_VAR) {
  message("Processing pseudobulk level: ", group_col)

  seurat_pb <- aggregate_pseudobulk_object(
    seu_obj = seu_obj,
    sample_col = sample_col,
    group_col = group_col
  )

  qc_df <- get_pseudobulk_qc(
    seu_obj = seu_obj,
    seurat_pb = seurat_pb,
    sample_col = sample_col,
    group_col = group_col
  )

  p <- plot_pseudobulk_qc(
    qc_df,
    group_col = group_col,
    min_cells = min_cells,
    min_counts = min_counts
  )
  ggsave(
    file.path(fig_dir, paste0("pseudobulk_qc_", group_col, ".pdf")),
    plot = p,
    height = 8.5,
    width = 10
  )

  seurat_pb <- add_pseudobulk_qc_metadata(
    seurat_pb = seurat_pb,
    qc_df = qc_df,
    group_col = group_col
  )

  seurat_pb <- filter_pseudobulk_object(
    seurat_pb = seurat_pb,
    group_col = group_col,
    min_cells = min_cells,
    min_counts = min_counts,
    output_dir = output_dir
  )

  seurat_pb <- run_pseudobulk_exploratory_analysis(
    seurat_pb = seurat_pb,
    group_col = group_col,
    fig_dir = fig_dir,
    n_pcs_assoc = n_pcs_assoc,
    n_pcs_var = n_pcs_var
  )

  run_celltype_subset_exploratory_analysis(
    seurat_pb = seurat_pb,
    group_col = group_col,
    fig_dir = fig_dir,
    n_pcs_assoc = n_pcs_assoc,
    n_pcs_var = n_pcs_var
  )

  export_pseudobulk_assay(
    seurat_pb = seurat_pb,
    group_col = group_col,
    output_dir = output_dir
  )

  invisible(seurat_pb)
}

pb_broad <- process_pseudobulk_level(
  seu_obj = seu,
  sample_col = "sample",
  group_col = "neuron_class"
)

pb_fine <- process_pseudobulk_level(
  seu_obj = seu,
  sample_col = "sample",
  group_col = "cell_type_concise"
)
