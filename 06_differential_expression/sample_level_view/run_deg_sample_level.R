#!/usr/bin/env Rscript
# ==============================================================================
# run_deg_sample_level.R
#
# Pseudobulk sample-level DEG on feature-filtered pseudobulk objects.
#
# Flow:
#   1. Load cell-type-specific feature-filtered pseudobulk Seurat objects from
#      pseduobulking_qc_gene.R.
#   2. For each clustering resolution, cell type, sex, and timepoint, test
#      mut vs wt using design <- model.matrix(~ animal_genotype).
#   3. Run Limma-Voom, EdgeR-zingeR, and DESeq2 when their packages are
#      available.
#   4. Save EDA plots before DEG fitting where useful, including voom
#      mean-variance plots and DESeq2 dispersion plots.
#   5. Save DEG CSVs with shared DEG visualization columns.
#   6. Generate volcano plots and top-gene heatmaps in the same run.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(edgeR)
  library(limma)
  library(grid)
  library(circlize)
})

DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
BASE_FILE_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg/sample_level_view"
BASE_FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/sample_level_view"

CLUSTER_RESOLUTIONS <- c("neuron_class", "cell_type_concise")
METHODS <- c("limma_voomcc", "edger_zinger", "deseq2")

GENOTYPE_REF <- "wt"
GENOTYPE_TEST <- "mut"
MIN_SAMPLES_PER_GENOTYPE <- 2
DISCARD_CELL_TYPES <- list(
  cell_type_concise = c("Fibroblasts")
)

dir.create(BASE_FILE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(BASE_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

parse_cli_args <- function(default_methods, default_clusters) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(methods = default_methods, clusters = default_clusters, viz_only = FALSE)
  if (length(args) == 0) return(out)

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--viz_only") {
      out$viz_only <- TRUE
      i <- i + 1
    } else if (args[i] == "--method" && i < length(args)) {
      out$methods <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--cluster_col" && i < length(args)) {
      out$clusters <- args[i + 1]
      i <- i + 2
    } else {
      stop("Unknown or incomplete argument: ", args[i])
    }
  }

  bad_methods <- setdiff(out$methods, default_methods)
  bad_clusters <- setdiff(out$clusters, default_clusters)
  if (length(bad_methods) > 0) stop("Unsupported method(s): ", paste(bad_methods, collapse = ", "))
  if (length(bad_clusters) > 0) stop("Unsupported cluster_col(s): ", paste(bad_clusters, collapse = ", "))
  out
}

CLI <- parse_cli_args(METHODS, CLUSTER_RESOLUTIONS)

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unknown", x)
}

`%||%` <- function(a, b) {
  if (length(a) == 0 || is.na(a)) b else a
}

method_available <- function(method) {
  switch(method,
    limma_voomcc = requireNamespace("limma", quietly = TRUE) &&
      requireNamespace("edgeR", quietly = TRUE),
    edger_zinger = requireNamespace("edgeR", quietly = TRUE) &&
      requireNamespace("zingeR", quietly = TRUE),
    deseq2 = requireNamespace("DESeq2", quietly = TRUE),
    FALSE
  )
}

read_feature_subsets <- function(cluster_col) {
  summary_csv <- file.path(
    BASE_FILE_DIR,
    paste0("pseudobulk_gene_filter_summary_", cluster_col, "_feature.csv")
  )
  if (!file.exists(summary_csv)) {
    stop("Feature-filter summary CSV not found: ", summary_csv)
  }

  subsets <- read.csv(summary_csv, stringsAsFactors = FALSE)$feature_subset
  subsets <- sort(unique(as.character(subsets)))
  discard <- DISCARD_CELL_TYPES[[cluster_col]] %||% character()
  setdiff(subsets, discard)
}

load_feature_subset_object <- function(cluster_col, cell_type) {
  rds_path <- file.path(
    DATA_DIR,
    paste0("pseudobulk_", cluster_col, "_", safe_filename(cell_type), "_feature.rds")
  )
  if (!file.exists(rds_path)) {
    stop("Feature-filtered cell-type Seurat object not found: ", rds_path)
  }
  readRDS(rds_path)
}

get_counts_meta <- function(obj) {
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  meta <- obj[[]]
  meta <- meta[colnames(counts), , drop = FALSE]
  stopifnot(identical(rownames(meta), colnames(counts)))
  list(counts = counts, meta = meta)
}

make_design <- function(meta) {
  meta$animal_genotype <- factor(meta$animal_genotype, levels = c(GENOTYPE_REF, GENOTYPE_TEST))
  model.matrix(~ animal_genotype, data = meta)
}

coef_name <- function(design) {
  grep("^animal_genotype", colnames(design), value = TRUE)[1]
}

format_deg_table <- function(tab,
                             method,
                             cluster_col,
                             cell_type,
                             sex_i,
                             tp_i,
                             n_wt,
                             n_mut) {
  tab <- as.data.frame(tab)
  tab$gene <- rownames(tab)

  logfc_col <- intersect(c("logFC", "log2FoldChange"), colnames(tab))[1]
  p_col <- intersect(c("P.Value", "PValue", "pvalue"), colnames(tab))[1]
  padj_col <- intersect(c("adj.P.Val", "padjFilter", "FDR", "padj"), colnames(tab))[1]

  if (is.na(logfc_col) || is.na(p_col)) {
    stop("Could not identify logFC or p-value columns for method: ", method)
  }
  if (is.na(padj_col)) {
    tab$p_val_adj_tmp <- p.adjust(tab[[p_col]], method = "BH")
    padj_col <- "p_val_adj_tmp"
  }

  out <- tab %>%
    mutate(
      avg_log2FC = .data[[logfc_col]],
      p_val = .data[[p_col]],
      p_val_adj = .data[[padj_col]],
      animal_sex = sex_i,
      animal_timepoint = tp_i,
      method = method,
      cluster_col = cluster_col,
      n_wt_pseudobulk = n_wt,
      n_mut_pseudobulk = n_mut
    )
  out[[cluster_col]] <- cell_type
  out
}

save_voom_plot <- function(dge, design, out_path, quality_weights = FALSE) {
  pdf(out_path, width = 7, height = 6)
  if (quality_weights) {
    v <- limma::voomWithQualityWeights(dge, design = design, plot = TRUE)
  } else {
    v <- limma::voom(dge, design = design, plot = TRUE)
  }
  dev.off()
  v
}

make_biological_replicate_block <- function(meta) {
  required_cols <- c("animal_genotype", "animal_sex", "animal_timepoint")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Cannot define biological replicate block; missing metadata column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  interaction(
    meta$animal_genotype,
    meta$animal_sex,
    meta$animal_timepoint,
    drop = TRUE,
    sep = "_"
  )
}

run_limma_voomcc <- function(counts, meta, design, coef, eda_dir, prefix) {
  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)

  v <- save_voom_plot(
    dge = dge,
    design = design,
    out_path = file.path(eda_dir, paste0(prefix, "_voom_mean_variance.pdf")),
    quality_weights = FALSE
  )

  block <- make_biological_replicate_block(meta)
  use_block <- anyDuplicated(block) > 0

  if (use_block) {
    dupcor <- limma::duplicateCorrelation(v, design = design, block = block)
    fit <- limma::lmFit(v, design = design, block = block, correlation = dupcor$consensus)
  } else {
    dupcor <- NULL
    fit <- limma::lmFit(v, design = design)
  }

  fit <- limma::eBayes(fit, trend = TRUE)
  tab <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  attr(tab, "consensus_correlation") <- if (is.null(dupcor)) NA_real_ else dupcor$consensus
  attr(tab, "duplicate_correlation_block") <- if (use_block) "animal_genotype_animal_sex_animal_timepoint" else NA_character_
  tab
}

run_edger_zinger <- function(counts, design, coef, eda_dir, prefix) {
  dge <- edgeR::DGEList(counts = round(counts))
  dge <- edgeR::calcNormFactors(dge)

  weights <- zingeR::zeroWeightsLS(
    counts = dge$counts,
    design = design,
    maxit = 200,
    normalization = "TMM"
  )
  dge$weights <- weights
  dge <- edgeR::estimateDisp(dge, design)

  pdf(file.path(eda_dir, paste0(prefix, "_edger_zinger_bcv.pdf")), width = 7, height = 6)
  edgeR::plotBCV(dge)
  dev.off()

  fit <- edgeR::glmFit(dge, design)
  lrt <- zingeR::glmWeightedF(fit, coef = coef, independentFiltering = TRUE)
  edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
}

run_deseq2 <- function(counts, meta, design_formula, eda_dir, prefix) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = meta,
    design = design_formula
  )
  dds$animal_genotype <- stats::relevel(factor(dds$animal_genotype), ref = GENOTYPE_REF)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)

  pdf(file.path(eda_dir, paste0(prefix, "_deseq2_dispersion.pdf")), width = 7, height = 6)
  DESeq2::plotDispEsts(dds)
  dev.off()

  dds <- DESeq2::nbinomWaldTest(dds)
  result_name <- paste0("animal_genotype_", GENOTYPE_TEST, "_vs_", GENOTYPE_REF)
  if (!result_name %in% DESeq2::resultsNames(dds)) {
    stop("DESeq2 result name not found: ", result_name)
  }
  res <- as.data.frame(DESeq2::results(dds, name = result_name))
  res
}

run_one_stratum <- function(obj,
                            method,
                            cluster_col,
                            cell_type,
                            sex_i,
                            tp_i,
                            out_dir,
                            eda_dir) {
  cm <- get_counts_meta(obj)
  meta <- cm$meta
  keep <- meta$animal_sex == sex_i &
    meta$animal_timepoint == tp_i &
    meta[[cluster_col]] == cell_type &
    meta$animal_genotype %in% c(GENOTYPE_REF, GENOTYPE_TEST)

  meta_i <- droplevels(meta[keep, , drop = FALSE])
  counts_i <- cm$counts[, rownames(meta_i), drop = FALSE]
  genotype_tab <- table(meta_i$animal_genotype)
  n_wt <- unname(genotype_tab[GENOTYPE_REF] %||% 0)
  n_mut <- unname(genotype_tab[GENOTYPE_TEST] %||% 0)

  if (n_wt < MIN_SAMPLES_PER_GENOTYPE || n_mut < MIN_SAMPLES_PER_GENOTYPE) {
    message(
      "  Skipping ", cell_type, " / ", sex_i, " / ", tp_i,
      " (n_wt=", n_wt, ", n_mut=", n_mut, ")"
    )
    return(NULL)
  }

  counts_i <- counts_i[Matrix::rowSums(counts_i) > 0, , drop = FALSE]
  meta_i$animal_genotype <- factor(meta_i$animal_genotype, levels = c(GENOTYPE_REF, GENOTYPE_TEST))
  design <- make_design(meta_i)
  coef <- coef_name(design)
  if (is.na(coef)) {
    stop("Could not identify animal_genotype coefficient in design.")
  }

  prefix <- paste(method, safe_filename(sex_i), safe_filename(tp_i), safe_filename(cell_type), sep = "_")
  message("  Running ", method, ": ", sex_i, " / ", tp_i, " / ", cell_type)

  tab <- switch(method,
    limma_voomcc = run_limma_voomcc(counts_i, meta_i, design, coef, eda_dir, prefix),
    edger_zinger = run_edger_zinger(counts_i, design, coef, eda_dir, prefix),
    deseq2 = run_deseq2(counts_i, meta_i, ~ animal_genotype, eda_dir, prefix),
    stop("Unknown method: ", method)
  )

  deg <- format_deg_table(
    tab = tab,
    method = method,
    cluster_col = cluster_col,
    cell_type = cell_type,
    sex_i = sex_i,
    tp_i = tp_i,
    n_wt = n_wt,
    n_mut = n_mut
  )

  out_csv <- file.path(
    out_dir,
    sprintf("deg_%s_%s_%s.csv", sex_i, tp_i, safe_filename(cell_type))
  )
  write.csv(deg, out_csv, row.names = FALSE)
  message("    Saved: ", out_csv)
  deg
}

# ---------------------------------------------------------------------------- #
# Integrated visualization helpers                                               #
# ---------------------------------------------------------------------------- #
timepoint_label <- function(x) {
  paste0("P", gsub("[^0-9]", "", as.character(x)))
}

timepoint_levels_for <- function(x) {
  present <- unique(timepoint_label(x))
  preferred <- c("P30", "P60", "P120", "P150")
  c(
    preferred[preferred %in% present],
    setdiff(present[order(as.numeric(gsub("[^0-9]", "", present)))], preferred)
  )
}

format_sex_label <- function(x) {
  case_when(
    tolower(as.character(x)) %in% c("male", "m") ~ "Male",
    tolower(as.character(x)) %in% c("female", "f") ~ "Female",
    TRUE ~ as.character(x)
  )
}

sex_title <- function(x) {
  case_when(
    format_sex_label(x) == "Male" ~ "Males",
    format_sex_label(x) == "Female" ~ "Females",
    TRUE ~ format_sex_label(x)
  )
}

ditto_seq_colors <- c(
  "17" = "#FFBE2D",
  "18" = "#80C7EF",
  "19" = "#00F6B3",
  "21" = "#06A5FF",
  "22" = "#FF8320"
)

stepped_colors <- c(
  "#990F26FF", "#B33E52FF", "#CC7A88FF", "#E6B8BFFF",
  "#99600FFF", "#B3823EFF", "#CCAA7AFF", "#E6D2B8FF",
  "#54990FFF", "#78B33EFF", "#A3CC7AFF", "#CFE6B8FF",
  "#0F8299FF", "#3E9FB3FF", "#7ABECCFF", "#B8DEE6FF",
  "#3D0F99FF", "#653EB3FF", "#967ACCFF", "#C7B8E6FF",
  "#333333FF", "#666666FF", "#999999FF", "#CCCCCCFF"
)

sex_colors_map <- c("Male" = ditto_seq_colors[["21"]], "Female" = ditto_seq_colors[["22"]])
timepoint_colors_map <- c(
  "P30" = stepped_colors[19],
  "P60" = stepped_colors[18],
  "P120" = stepped_colors[17],
  "P150" = stepped_colors[17]
)

cell_type_concise_colors <- c(
  "GABA" = "#3B00FB",
  "GLU" = "#66B0FF",
  "Astrocytes" = "#683B79",
  "Oligodendrocytes" = "#D85FF7",
  "OPC" = "#1C7F93",
  "Immune" = "#7ED7D1",
  "Ependymal-like" = "#B5EFB5",
  "Fibroblasts" = "#822E1C",
  "Mural+Endothelial" = "#BDCDFF",
  "ParsTuber" = "#AAF400"
)

cell_type_concise_levels <- c(
  "GABA",
  "GLU",
  "Astrocytes",
  "Ependymal-like",
  "Immune",
  "Oligodendrocytes",
  "OPC",
  "Fibroblasts",
  "Mural+Endothelial",
  "ParsTuber"
)

neuron_class_colors <- c(
  "GABA" = "#782AB6",
  "GLU" = "#C075A6",
  "Non-neuronal" = "#F7E1A0"
)

neuron_class_levels <- c("GABA", "GLU", "Non-neuronal")

cluster_levels_for <- function(cluster_col, values) {
  preferred <- switch(
    cluster_col,
    cell_type_concise = cell_type_concise_levels,
    neuron_class = neuron_class_levels,
    character()
  )
  values <- as.character(values)
  c(preferred[preferred %in% values], setdiff(values, preferred))
}

cluster_colors_for <- function(cluster_col, levels_use) {
  color_map <- switch(
    cluster_col,
    cell_type_concise = cell_type_concise_colors,
    neuron_class = neuron_class_colors,
    NULL
  )
  colors <- unname(color_map[as.character(levels_use)])
  names(colors) <- as.character(levels_use)
  colors[is.na(colors)] <- "#8C8C8C"
  colors
}

strip_labels <- function(grob) {
  labels <- character()
  if (inherits(grob, "text")) labels <- c(labels, as.character(grob$label))
  if (!is.null(grob$grobs)) labels <- c(labels, unlist(lapply(grob$grobs, strip_labels), use.names = FALSE))
  if (!is.null(grob$children)) labels <- c(labels, unlist(lapply(grob$children, strip_labels), use.names = FALSE))
  labels[nzchar(labels)]
}

recolor_strip_rects <- function(grob, color) {
  if (inherits(grob, "rect")) {
    grob$gp$fill <- adjustcolor(color, alpha.f = 0.10)
    grob$gp$col <- color
    grob$gp$lwd <- 1.2
  }
  if (!is.null(grob$grobs)) grob$grobs <- lapply(grob$grobs, recolor_strip_rects, color = color)
  if (!is.null(grob$children)) {
    grob$children <- do.call(grid::gList, lapply(grob$children, recolor_strip_rects, color = color))
  }
  grob
}

lookup_strip_color <- function(color_map, label) {
  if (is.null(color_map) || is.na(label)) return(NULL)
  color <- unname(color_map[as.character(label)])
  if (length(color) == 0 || is.na(color)) return(NULL)
  color
}

color_facet_strips <- function(plot, top_colors = NULL, row_colors = NULL) {
  gt <- ggplotGrob(plot)
  strip_idx <- grep("^strip-", gt$layout$name)
  for (i in strip_idx) {
    strip_name <- gt$layout$name[i]
    label <- strip_labels(gt$grobs[[i]])[1]
    if (is.na(label)) next
    strip_color <- NULL
    if (grepl("^strip-t|^strip-b", strip_name)) {
      strip_color <- lookup_strip_color(top_colors, label)
    } else if (grepl("^strip-l|^strip-r", strip_name)) {
      strip_color <- lookup_strip_color(row_colors, label)
    }
    if (!is.null(strip_color) && !is.na(strip_color)) {
      gt$grobs[[i]] <- recolor_strip_rects(gt$grobs[[i]], strip_color)
    }
  }
  gt
}

save_colored_facet_pdf <- function(plot, out_path, width, height, top_colors = NULL, row_colors = NULL) {
  pdf(out_path, width = width, height = height)
  grid.newpage()
  grid.draw(color_facet_strips(plot, top_colors = top_colors, row_colors = row_colors))
  dev.off()
}

bgr_col_fun <- function(mat) {
  max_abs <- max(abs(mat), na.rm = TRUE)
  lim <- floor(max_abs)
  if (lim == 0) lim <- max_abs
  circlize::colorRamp2(c(-lim, 0, lim), c("#2166AC", "#EEEEEE", "#B2182B"))
}

add_volcano_label_layers <- function(plot, labels, sig_colors) {
  if (nrow(labels) == 0) return(plot)
  labels$sig <- as.character(labels$sig)

  label_specs <- list(
    "down in mut" = list(text = sig_colors[["down in mut"]], segment = sig_colors[["down in mut"]]),
    "up in mut" = list(text = sig_colors[["up in mut"]], segment = sig_colors[["up in mut"]]),
    "ns" = list(text = "black", segment = "grey50")
  )

  for (label_group in names(label_specs)) {
    label_df <- labels %>% filter(sig == label_group)
    if (nrow(label_df) == 0) next

    plot <- plot +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = gene),
        size = 2.3,
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.15,
        force = 2,
        force_pull = 0.05,
        max.time = 2,
        max.iter = 10000,
        seed = 123,
        min.segment.length = 0,
        segment.size = 0.25,
        segment.color = label_specs[[label_group]]$segment,
        color = label_specs[[label_group]]$text,
        show.legend = FALSE
      )
  }

  plot
}

plot_deg_volcanoes <- function(deg, fig_dir, cluster_col, method, top_n = 10) {
  cluster_col_name <- cluster_col

  deg <- deg %>%
    filter(!is.na(p_val_adj)) %>%
    mutate(
      neglog10padj = -log10(pmax(p_val_adj, 1e-300)),
      stratum = paste(animal_sex, animal_timepoint, .data[[cluster_col_name]], sep = " | "),
      sig = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "up in mut",
        p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "down in mut",
        TRUE ~ "ns"
      )
    )

  if (nrow(deg) == 0) return(invisible(NULL))

  top_n_genes <- deg %>%
    group_by(stratum) %>%
    slice_min(order_by = p_val_adj, n = top_n, with_ties = FALSE) %>%
    ungroup()

  sig_levels <- c("down in mut", "ns", "up in mut")
  sig_colors <- setNames(c("#2166AC", "grey80", "#B2182B"), sig_levels)

  bordered_theme <- theme_minimal(base_size = 8) +
    theme(
      strip.text = element_text(size = 7, face = "bold"),
      strip.background = element_rect(fill = "white", color = "grey50", linewidth = 0.6),
      strip.placement = "outside",
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.spacing = unit(0.6, "lines")
    )

  for (sx in unique(deg$animal_sex)) {
    deg_sx <- deg %>%
      filter(animal_sex == sx) %>%
      mutate(
        sig = factor(sig, levels = sig_levels),
        timepoint_display = factor(timepoint_label(animal_timepoint), levels = timepoint_levels_for(animal_timepoint)),
        cluster_display = factor(
          .data[[cluster_col_name]],
          levels = cluster_levels_for(cluster_col_name, unique(.data[[cluster_col_name]]))
        )
      )
    top_sx <- top_n_genes %>%
      filter(animal_sex == sx) %>%
      mutate(
        timepoint_display = factor(timepoint_label(animal_timepoint), levels = levels(deg_sx$timepoint_display)),
        cluster_display = factor(.data[[cluster_col_name]], levels = levels(deg_sx$cluster_display))
      )

    n_tp <- length(levels(droplevels(deg_sx$timepoint_display)))
    n_class <- length(levels(droplevels(deg_sx$cluster_display)))
    y_max_sx <- max(deg_sx$neglog10padj, na.rm = TRUE)
    timepoint_strip_colors <- timepoint_colors_map[levels(droplevels(deg_sx$timepoint_display))]
    if (any(is.na(timepoint_strip_colors))) timepoint_strip_colors[is.na(timepoint_strip_colors)] <- stepped_colors[23]
    cluster_strip_colors <- cluster_colors_for(cluster_col, levels(droplevels(deg_sx$cluster_display)))

    p_base <- ggplot(deg_sx, aes(x = avg_log2FC, y = neglog10padj, color = sig)) +
      geom_point(size = 0.6, alpha = 0.7) +
      scale_color_manual(values = sig_colors, name = NULL, drop = FALSE) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2, color = "grey50") +
      coord_cartesian(ylim = c(0, y_max_sx * 1.03)) +
      bordered_theme +
      labs(
        title = paste("Sample-level DEG volcano plots -", method, "-", format_sex_label(sx), "(genotype: mut vs wt)"),
        x = "avg log2FC (mut vs wt)",
        y = "-log10(adjusted p-value)"
      )
    p_base <- add_volcano_label_layers(p_base, top_sx, sig_colors)

    out_path <- file.path(fig_dir, sprintf("volcano_%s_rows_timepoint_cols_%s.pdf", sx, cluster_col))
    save_colored_facet_pdf(
      p_base + facet_grid(timepoint_display ~ cluster_display, scales = "free_x", switch = "y"),
      out_path,
      width = max(8, n_class * 3),
      height = max(5, n_tp * 3),
      top_colors = cluster_strip_colors,
      row_colors = timepoint_strip_colors
    )

    out_path <- file.path(fig_dir, sprintf("volcano_%s_rows_%s_cols_timepoint.pdf", sx, cluster_col))
    save_colored_facet_pdf(
      p_base + facet_grid(cluster_display ~ timepoint_display, scales = "free_x", switch = "y"),
      out_path,
      width = max(8, n_tp * 3),
      height = max(5, n_class * 3),
      top_colors = timepoint_strip_colors,
      row_colors = cluster_strip_colors
    )
  }
}

plot_deg_heatmaps <- function(deg, fig_dir, cluster_col, top_n = 5) {
  cluster_col_name <- cluster_col

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    message("ComplexHeatmap/circlize is not installed; skipping heatmaps.")
    return(invisible(NULL))
  }

  deg <- deg %>%
    filter(!is.na(p_val_adj)) %>%
    mutate(stratum = paste(animal_sex, animal_timepoint, .data[[cluster_col_name]], sep = " | "))
  top_n_genes <- deg %>%
    group_by(stratum) %>%
    slice_min(order_by = p_val_adj, n = top_n, with_ties = FALSE) %>%
    ungroup()

  for (sex_label in unique(deg$animal_sex)) {
    deg_sx <- deg %>% filter(animal_sex == sex_label)
    top_sx <- top_n_genes %>% filter(animal_sex == sex_label)
    heat_genes <- unique(top_sx$gene)
    if (length(heat_genes) == 0) next

    wide <- deg_sx %>%
      filter(gene %in% heat_genes) %>%
      select(gene, stratum, avg_log2FC) %>%
      pivot_wider(names_from = stratum, values_from = avg_log2FC, values_fill = 0)
    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$gene

    padj_wide <- deg_sx %>%
      filter(gene %in% heat_genes) %>%
      select(gene, stratum, p_val_adj) %>%
      pivot_wider(names_from = stratum, values_from = p_val_adj, values_fill = 1)
    padj_mat <- as.matrix(padj_wide[, -1])
    rownames(padj_mat) <- padj_wide$gene
    padj_mat <- padj_mat[rownames(mat), colnames(mat)]

    col_meta <- strsplit(colnames(mat), " \\| ") %>%
      do.call(rbind, .) %>%
      as.data.frame(stringsAsFactors = FALSE)
    colnames(col_meta) <- c("sex", "timepoint", cluster_col)
    col_meta$sex_display <- format_sex_label(col_meta$sex)
    col_meta$timepoint_display <- timepoint_label(col_meta$timepoint)
    col_meta[[cluster_col]] <- as.character(col_meta[[cluster_col]])

    tp_levels <- timepoint_levels_for(col_meta$timepoint)
    class_levels <- cluster_levels_for(cluster_col, unique(col_meta[[cluster_col]]))
    sex_levels <- unique(col_meta$sex_display)
    col_order <- order(
      factor(col_meta$timepoint_display, levels = tp_levels),
      factor(col_meta[[cluster_col]], levels = class_levels)
    )
    mat <- mat[, col_order, drop = FALSE]
    padj_mat <- padj_mat[, col_order, drop = FALSE]
    col_meta <- col_meta[col_order, , drop = FALSE]

    tp_colors <- setNames(unname(timepoint_colors_map[tp_levels]), tp_levels)
    if (any(is.na(tp_colors))) tp_colors[is.na(tp_colors)] <- "#8C8C8C"

    class_colors <- cluster_colors_for(cluster_col, class_levels)
    sex_colors <- setNames(unname(sex_colors_map[sex_levels]), sex_levels)
    if (any(is.na(sex_colors))) sex_colors[is.na(sex_colors)] <- "#8C8C8C"

    anno_df <- data.frame(
      `Cell Type` = factor(col_meta[[cluster_col]], levels = class_levels),
      `Time Point` = factor(col_meta$timepoint_display, levels = tp_levels),
      Sex = factor(col_meta$sex_display, levels = sex_levels),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    col_anno <- ComplexHeatmap::HeatmapAnnotation(
      df = anno_df,
      col = setNames(list(class_colors, tp_colors, sex_colors), colnames(anno_df)),
      annotation_name_side = "left",
      simple_anno_size = grid::unit(0.35, "cm")
    )

    sig_mat <- matrix(ifelse(padj_mat < 0.05, "*", ""), nrow = nrow(mat), dimnames = dimnames(mat))
    title_col <- lookup_strip_color(sex_colors_map, format_sex_label(sex_label)) %||% "#333333"
    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "avg log2FC\n(mut vs wt)",
      col = bgr_col_fun(mat),
      top_annotation = col_anno,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      border = title_col,
      row_names_gp = grid::gpar(fontsize = 10),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (sig_mat[i, j] != "") {
          grid::grid.text(
            "*",
            x,
            y - height * 0.12,
            gp = grid::gpar(fontsize = 15, fontface = "bold", col = "black")
          )
        }
      },
      heatmap_legend_param = list(
        direction = "vertical",
        title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9)
      ),
      column_title = sex_title(sex_label),
      column_title_gp = grid::gpar(col = title_col, fontface = "bold", fontsize = 12)
    )

    out_path <- file.path(fig_dir, sprintf("heatmap_top_genes_%s.pdf", sex_label))
    pdf(out_path, width = max(9, ncol(mat) * 0.35 + 3), height = max(6, nrow(mat) * 0.2 + 2))
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legend = TRUE
    )
    dev.off()
  }
}

run_method_resolution <- function(method, cluster_col) {
  if (!method_available(method)) {
    message("Skipping ", method, ": required package(s) not installed.")
    return(invisible(NULL))
  }

  out_dir <- file.path(BASE_FILE_DIR, method, cluster_col)
  fig_dir <- file.path(BASE_FIG_DIR, method, cluster_col)
  eda_dir <- file.path(fig_dir, "eda")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(eda_dir, recursive = TRUE, showWarnings = FALSE)

  message("========================================================")
  message("  Method      : ", method)
  message("  Cluster col : ", cluster_col)
  message("  Output dir  : ", out_dir)
  message("  Figure dir  : ", fig_dir)
  message("========================================================")

  combined_csv <- file.path(out_dir, sprintf("deg_all_%s_%s.csv", method, cluster_col))
  if (isTRUE(CLI$viz_only)) {
    if (!file.exists(combined_csv)) {
      message("Skipping viz-only ", method, " / ", cluster_col, ": combined CSV not found: ", combined_csv)
      return(invisible(NULL))
    }
    message("Viz-only mode: reading combined DEG table: ", combined_csv)
    deg_all <- read.csv(combined_csv, stringsAsFactors = FALSE)
    plot_deg_volcanoes(deg_all, fig_dir = fig_dir, cluster_col = cluster_col, method = method)
    plot_deg_heatmaps(deg_all, fig_dir = fig_dir, cluster_col = cluster_col)
    return(invisible(deg_all))
  }

  cell_types <- read_feature_subsets(cluster_col)
  results <- list()

  for (cell_type in cell_types) {
    obj <- load_feature_subset_object(cluster_col, cell_type)
    meta <- obj[[]]
    strata <- meta %>%
      distinct(animal_sex, animal_timepoint) %>%
      arrange(animal_sex, animal_timepoint)

    for (i in seq_len(nrow(strata))) {
      deg <- tryCatch(
        run_one_stratum(
          obj = obj,
          method = method,
          cluster_col = cluster_col,
          cell_type = cell_type,
          sex_i = strata$animal_sex[i],
          tp_i = strata$animal_timepoint[i],
          out_dir = out_dir,
          eda_dir = eda_dir
        ),
        error = function(e) {
          message("  Failed ", cell_type, " / ", strata$animal_sex[i], " / ", strata$animal_timepoint[i], ": ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(deg)) {
        results[[paste(cell_type, strata$animal_sex[i], strata$animal_timepoint[i], sep = "|")]] <- deg
      }
    }
  }

  if (length(results) == 0) {
    message("No DEG results for ", method, " / ", cluster_col)
    return(invisible(NULL))
  }

  deg_all <- bind_rows(results)
  write.csv(deg_all, combined_csv, row.names = FALSE)
  message("Combined DEG table saved: ", combined_csv)

  plot_deg_volcanoes(deg_all, fig_dir = fig_dir, cluster_col = cluster_col, method = method)
  plot_deg_heatmaps(deg_all, fig_dir = fig_dir, cluster_col = cluster_col)

  invisible(deg_all)
}

for (method in CLI$methods) {
  for (cluster_col in CLI$clusters) {
    run_method_resolution(method, cluster_col)
  }
}

message("Done.")
