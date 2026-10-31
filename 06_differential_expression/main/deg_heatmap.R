#!/usr/bin/env Rscript
# Heatmaps of pre-computed main/ DEG results.
# Named palettes, matrix construction and typography follow the supplied reference.
# One PDF per contrast; top_n genes per timepoint/cell type, ranked by BH p-value.
# Rscript main/viz_heatmap.R --method MAST_RE --cluster_col cell_type_concise
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
parser <- ArgumentParser(description = "Heatmaps of pairwise DEG contrasts")
parser$add_argument("--method", required = TRUE,
  choices = c("MAST", "MAST_RE", "voomLmFit", "voomLmFit_block", "wilcox"))
parser$add_argument("--cluster_col", default = "neuron_class",
  choices = c("neuron_class", "cell_type_concise"))
parser$add_argument("--deg_csv", default = NULL,
  help = "Explicit combined DEG CSV; otherwise resolved under results_dir")
parser$add_argument("--results_dir", default = "/quobyte/lasallegrp/Osman/shenyu/main/results")
parser$add_argument("--fig_dir", default = "/quobyte/lasallegrp/Osman/shenyu/main/figures",
  help = "Figure base directory; method/cluster_col appended")
parser$add_argument("--top_n", type = "integer", default = 5,
  help = "Top genes by adjusted p-value per contrast/timepoint/cell type")
parser$add_argument("--sig", action = "store_true",
  help = "Select top genes only among BH-adjusted p_val_adj < 0.05")
parser$add_argument("--main", action = "store_true",
  help = "Use avg_log2FC for every method instead of method-specific fold changes")
args <- parser$parse_args()
stopifnot(args$top_n >= 1L)
if (is.null(args$deg_csv)) {
  args$deg_csv <- file.path(args$results_dir, args$method, args$cluster_col,
    sprintf("deg_all_%s_%s.csv", args$method, args$cluster_col))
}
fig_dir <- file.path(args$fig_dir, args$method, args$cluster_col)
fc_col <- if (args$main) "avg_log2FC" else switch(args$method,
  MAST_RE = "avg_log2FC_viz",
  voomLmFit = "logFC",
  voomLmFit_block = "logFC",
  "avg_log2FC")
contrast_titles <- c(
  sex_wt = "Sex effect in wt (female - male)",
  sex_mut = "Sex effect in mut (female - male)",
  genotype_male = "Genotype effect in male (mut - wt)",
  genotype_female = "Genotype effect in female (mut - wt)"
)

timepoint_label <- function(x) {

  x <- as.character(x)
  x[x == "male_P120_female_P150"] <- "P120_P150"

  ifelse(
    x == "P120_P150",
    "P120_P150",
    paste0(
      "P",
      sub("^[Pp]", "", x)
    )
  )
}


timepoint_levels_for <- function(x) {

  present <- unique(
    timepoint_label(x)
  )

  preferred <- c(
    "P7",
    "P14",
    "P30",
    "P60",
    "P120_P150",
    "P120",
    "P150"
  )

  extra <- setdiff(
    present,
    preferred
  )

  if (length(extra)) {

    extra_num <- suppressWarnings(
      as.numeric(
        gsub("[^0-9]", "", extra)
      )
    )

    extra <- extra[
      order(
        extra_num,
        na.last = TRUE
      )
    ]
  }

  c(
    preferred[
      preferred %in% present
    ],
    extra
  )
}


stepped_colors <- c(
  "#990F26FF", "#B33E52FF", "#CC7A88FF", "#E6B8BFFF",
  "#99600FFF", "#B3823EFF", "#CCAA7AFF", "#E6D2B8FF",
  "#54990FFF", "#78B33EFF", "#A3CC7AFF", "#CFE6B8FF",
  "#0F8299FF", "#3E9FB3FF", "#7ABECCFF", "#B8DEE6FF",
  "#3D0F99FF", "#653EB3FF", "#967ACCFF", "#C7B8E6FF",
  "#333333FF", "#666666FF", "#999999FF", "#CCCCCCFF"
)


timepoint_colors_map <- c(
  P7 = stepped_colors[20],
  P14 = stepped_colors[19],
  P30 = stepped_colors[19],
  P60 = stepped_colors[18],
  P120 = stepped_colors[17],
  P150 = stepped_colors[17],
  P120_P150 = stepped_colors[17]
)


sex_colors_map <- c(
  Male = "#06A5FF",
  Female = "#FF8320"
)


# =============================================================================
# Cell type / neuron class definitions
# =============================================================================

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


cell_type_concise_colors <- c(
  GABA = "#3B00FB",
  GLU = "#66B0FF",
  Astrocytes = "#683B79",
  Oligodendrocytes = "#D85FF7",
  OPC = "#1C7F93",
  Immune = "#7ED7D1",
  `Ependymal-like` = "#B5EFB5",
  Fibroblasts = "#822E1C",
  `Mural+Endothelial` = "#BDCDFF",
  ParsTuber = "#AAF400"
)


neuron_class_levels <- c(
  "GABA",
  "GLU",
  "Non-neuronal"
)


neuron_class_colors <- c(
  GABA = "#782AB6",
  GLU = "#C075A6",
  `Non-neuronal` = "#F7E1A0"
)


cluster_levels_for <- function(cluster_col, values) {

  preferred <- if (
    cluster_col == "neuron_class"
  ) {
    neuron_class_levels
  } else {
    cell_type_concise_levels
  }

  values <- unique(
    as.character(values)
  )

  c(
    preferred[
      preferred %in% values
    ],
    setdiff(
      values,
      preferred
    )
  )
}


cluster_colors_for <- function(cluster_col, levels_use) {

  color_map <- if (
    cluster_col == "neuron_class"
  ) {
    neuron_class_colors
  } else {
    cell_type_concise_colors
  }

  ans <- unname(
    color_map[levels_use]
  )

  ans[is.na(ans)] <- "#8C8C8C"

  setNames(
    ans,
    levels_use
  )
}


get_tp_colors <- function(levels_use) {

  ans <- unname(
    timepoint_colors_map[
      levels_use
    ]
  )

  ans[is.na(ans)] <- "#8C8C8C"

  setNames(
    ans,
    levels_use
  )
}



heatmap_limit <- function(mat) {
  finite <- mat[is.finite(mat)]
  max_abs <- if (length(finite)) max(abs(finite)) else 0
  lim <- ceiling(max_abs)
  if (lim == 0) lim <- 1
  lim
}
bgr_col_fun <- function(mat, lim = heatmap_limit(mat)) {
  colorRamp2(c(-lim, 0, lim), c("#2166AC", "#EEEEEE", "#B2182B"))
}

# Keep the actual endpoints in the legend; format only their text labels.
legend_ticks <- function(lim) c(-lim, 0, lim)

select_heatmap_genes <- function(d, cluster_col, top_n, sig = FALSE) {
  top <- d %>%
    filter(is.finite(p_val_adj), is.finite(plot_log2FC), !sig | p_val_adj < 0.05) %>%
    group_by(across(all_of(c("animal_timepoint", cluster_col)))) %>%
    slice_min(p_val_adj, n = top_n, with_ties = FALSE) %>% ungroup()
  unique(top$gene)
}

# Compare the ranges of the actual selected heatmaps, including --sig filtering.
genotype_shared_limit <- function(deg, cluster_col, top_n, sig = FALSE) {
  limits <- vapply(c("genotype_male", "genotype_female"), function(ct) {
    d <- deg[deg$contrast == ct, , drop = FALSE]
    genes <- select_heatmap_genes(d, cluster_col, top_n, sig)
    if (!length(genes)) return(NA_real_)
    heatmap_limit(d$plot_log2FC[d$gene %in% genes & is.finite(d$p_val_adj)])
  }, numeric(1))
  available <- limits[is.finite(limits)]
  if (length(available)) min(available) else NULL
}

message("Input CSV: ", args$deg_csv)
message("Fold-change column: ", fc_col)
deg <- read.csv(args$deg_csv, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("gene", "p_val_adj", fc_col, "contrast", "animal_timepoint", args$cluster_col)
missing <- setdiff(required, names(deg))
if (length(missing)) stop("Input CSV is missing required columns: ", paste(missing, collapse = ", "))
if (!is.numeric(deg$p_val_adj) || !is.numeric(deg[[fc_col]])) {
  stop("Adjusted p-values and fold changes must be numeric")
}
unknown <- setdiff(unique(deg$contrast), names(contrast_titles))
if (length(unknown)) stop("Unrecognized contrast: ", paste(unknown, collapse = ", "))
deg$animal_timepoint <- as.character(deg$animal_timepoint)
deg$plot_log2FC <- deg[[fc_col]]
if (anyNA(deg[, c("gene", "contrast", "animal_timepoint", args$cluster_col)])) {
  stop("Missing gene, contrast, timepoint or cell type labels")
}
if (anyDuplicated(deg[, c("gene", "contrast", "animal_timepoint", args$cluster_col)])) {
  stop("Duplicate gene rows within a contrast/timepoint/cell type")
}
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

make_contrast_heatmap <- function(contrast_i, deg, fig_dir, cluster_col, top_n, sig = FALSE, scale_limit = NULL) {
  d <- deg[deg$contrast == contrast_i, , drop = FALSE]
  if (!nrow(d)) return(invisible(NULL))
  genes <- select_heatmap_genes(d, cluster_col, top_n, sig)
  if (!length(genes)) {
    message("No plottable genes for ", contrast_i, "; skipping")
    return(invisible(NULL))
  }
  # Match the reference: construct columns only from eligible selected results;
  # absent gene/stratum combinations are neutral for display (0 / padj 1).
  selected <- d %>%
    filter(gene %in% genes, is.finite(plot_log2FC), is.finite(p_val_adj)) %>%
    mutate(timepoint_display = timepoint_label(animal_timepoint),
      cluster_display = as.character(.data[[cluster_col]]),
      stratum = paste(timepoint_display, cluster_display, sep = " | "))
  wide <- selected %>% select(gene, stratum, plot_log2FC) %>%
    pivot_wider(names_from = stratum, values_from = plot_log2FC, values_fill = 0)
  pwide <- selected %>% select(gene, stratum, p_val_adj) %>%
    pivot_wider(names_from = stratum, values_from = p_val_adj, values_fill = 1)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$gene
  padj_mat <- as.matrix(pwide[, -1, drop = FALSE])
  rownames(padj_mat) <- pwide$gene
  padj_mat <- padj_mat[rownames(mat), colnames(mat), drop = FALSE]
  col_meta <- selected %>% distinct(stratum, timepoint_display, cluster_display)
  col_meta <- col_meta[match(colnames(mat), col_meta$stratum), , drop = FALSE]
  tp_levels <- timepoint_levels_for(col_meta$timepoint_display)
  class_levels <- cluster_levels_for(cluster_col, col_meta$cluster_display)
  ord <- order(match(col_meta$timepoint_display, tp_levels),
    match(col_meta$cluster_display, class_levels))
  mat <- mat[, ord, drop = FALSE]
  padj_mat <- padj_mat[, ord, drop = FALSE]
  col_meta <- col_meta[ord, , drop = FALSE]
  anno_df <- data.frame(
    `Cell Type` = factor(col_meta$cluster_display, levels = class_levels),
    `Time Point` = factor(col_meta$timepoint_display, levels = tp_levels),
    check.names = FALSE)
  anno_colors <- list(`Cell Type` = cluster_colors_for(cluster_col, class_levels),
    `Time Point` = get_tp_colors(tp_levels))
  title_col <- "#333333"
  if (contrast_i %in% c("genotype_male", "genotype_female")) {
    sex_label <- if (contrast_i == "genotype_male") "Male" else "Female"
    anno_df$Sex <- factor(rep(sex_label, nrow(col_meta)), levels = c("Male", "Female"))
    anno_colors$Sex <- sex_colors_map
    title_col <- unname(sex_colors_map[sex_label])
  }
  col_anno <- HeatmapAnnotation(df = anno_df, col = anno_colors,
    annotation_name_side = "left", simple_anno_size = unit(0.35, "cm"))
  sig_mat <- matrix(ifelse(is.finite(padj_mat) & padj_mat < 0.05, "*", ""),
    nrow = nrow(mat), dimnames = dimnames(mat))
  cluster_rows <- nrow(mat) > 1L
  if (is.null(scale_limit)) scale_limit <- heatmap_limit(mat)
  ht <- Heatmap(mat, name = if (fc_col == "avg_log2FC_viz") "log2_coef" else if (fc_col == "logFC") "logFC" else "avg log2FC",
    col = bgr_col_fun(mat, scale_limit),
    top_annotation = col_anno, cluster_columns = FALSE, cluster_rows = cluster_rows,
    show_column_names = FALSE, border = title_col,
    row_names_gp = gpar(fontsize = 10),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y - height * 0.12,
        gp = gpar(fontsize = 15, fontface = "bold"))
    },
    heatmap_legend_param = list(direction = "vertical",
      at = legend_ticks(scale_limit),
      labels = format(legend_ticks(scale_limit), trim = TRUE, digits = 3)),
    column_title = unname(contrast_titles[contrast_i]),
    column_title_gp = gpar(col = title_col, fontface = "bold", fontsize = 12))
  out_path <- file.path(fig_dir, paste0(if (args$main) "main_" else "",
    if (sig) "sig_" else "",
    sprintf("heatmap_top%d_%s.pdf", top_n, contrast_i)))
  pdf(out_path, width = max(9, ncol(mat) * 0.35 + 3), height = max(6, nrow(mat) * 0.2 + 2))
  tryCatch(draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
    merge_legend = TRUE), finally = dev.off())
  message("Saved: ", out_path)
}
shared_limit <- genotype_shared_limit(deg, args$cluster_col, args$top_n, args$sig)
if (!is.null(shared_limit)) message("Genotype male/female shared color range: [", -shared_limit, ", ", shared_limit, "]")
for (contrast_i in names(contrast_titles)) {
  make_contrast_heatmap(contrast_i, deg, fig_dir, args$cluster_col, args$top_n, sig = args$sig,
    scale_limit = if (contrast_i %in% c("genotype_male", "genotype_female")) shared_limit else NULL)
}
message("Done. Figures in: ", fig_dir)
