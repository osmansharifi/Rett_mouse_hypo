#!/usr/bin/env Rscript
# Rows = genes; columns = methods. One heatmap per contrast/timepoint/cell type.
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
parser <- ArgumentParser(description = "Compare gene fold changes across DEG methods")
parser$add_argument("--method", nargs = "+",
  default = c("wilcox", "MAST_original", "MAST_RE", "voomLmFit", "voomLmFit_block"),
  choices = c("MAST_original", "wilcox", "MAST_RE", "voomLmFit", "voomLmFit_block"),
  help = "Methods to compare, in column order; default: all five")
parser$add_argument("--cluster_col", default = "neuron_class",
  choices = c("neuron_class", "cell_type_concise"))
parser$add_argument("--results_dir", default = "/quobyte/lasallegrp/Osman/shenyu/main/results")
parser$add_argument("--fig_dir", default = "/quobyte/lasallegrp/Osman/shenyu/main/figures/logFC")
parser$add_argument("--top_n", type = "integer", default = 5,
  help = "Select the union of top N genes per method in each panel")
parser$add_argument("--sig", action = "store_true",
  help = "Restrict top-gene selection to BH p_val_adj < 0.05")
parser$add_argument("--genes", nargs = "+", default = NULL,
  help = "Explicit genes in row order; bypasses top_n and sig selection")
args <- parser$parse_args()
stopifnot(args$top_n >= 1L)
methods <- unique(args$method)
# MAST_original and MAST_RE use their model-based estimates in log2 units.
fc_columns <- c(wilcox = "avg_log2FC", MAST_original = "avg_log2FC_viz",
  MAST_RE = "avg_log2FC_viz", voomLmFit = "logFC", voomLmFit_block = "logFC")
contrast_titles <- c(sex_wt = "Sex effect in wt (female - male)",
  sex_mut = "Sex effect in mut (female - male)",
  genotype_male = "Genotype effect in male (mut - wt)",
  genotype_female = "Genotype effect in female (mut - wt)")
safe_filename <- function(x) gsub("[^A-Za-z0-9_-]+", "_", as.character(x))
timepoint_label <- function(x) {
  if (x == "male_P120_female_P150") return("Male P120 / Female P150")
  if (grepl("^[0-9]+$", x)) paste0("P", x) else x
}
select_panel_genes <- function(d, methods, top_n, sig, genes = NULL) {
  if (!is.null(genes)) return(unique(genes))
  unique(unlist(lapply(methods, function(m) {
    use <- d[d$method == m & is.finite(d$p_val_adj) & is.finite(d$fold_change), , drop = FALSE]
    if (sig) use <- use[use$p_val_adj < 0.05, , drop = FALSE]
    head(use$gene[order(use$p_val_adj, use$gene)], top_n)
  }), use.names = FALSE))
}
make_method_matrix <- function(d, genes, methods) {
  mat <- matrix(NA_real_, length(genes), length(methods), dimnames = list(genes, methods))
  use <- d[d$gene %in% genes & d$method %in% methods, , drop = FALSE]
  mat[cbind(match(use$gene, genes), match(use$method, methods))] <- use$fold_change
  mat[!is.finite(mat)] <- NA_real_
  mat
}
scale_limit <- function(mat) {
  values <- mat[is.finite(mat)]
  if (!length(values)) return(NA_real_)
  max(1, ceiling(max(abs(values))))
}
format_cells <- function(mat) {
  labels <- matrix("NA", nrow(mat), ncol(mat), dimnames = dimnames(mat))
  ok <- is.finite(mat)
  values <- round(mat[ok], 3)
  values[values == 0] <- 0
  labels[ok] <- sprintf("%.3f", values)
  labels
}

tables <- lapply(methods, function(m) {
  path <- file.path(args$results_dir, m, args$cluster_col,
    sprintf("deg_all_%s_%s.csv", m, args$cluster_col))
  if (!file.exists(path)) stop("Missing method CSV: ", path)
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  fc <- fc_columns[[m]]
  required <- c("gene", "contrast", "animal_timepoint", args$cluster_col, "p_val_adj", fc)
  missing <- setdiff(required, names(d))
  if (length(missing)) stop(path, ": missing columns ", paste(missing, collapse = ", "))
  if (!is.numeric(d[[fc]]) || !is.numeric(d$p_val_adj)) stop("Non-numeric statistics: ", path)
  ans <- data.frame(gene = d$gene, contrast = d$contrast,
    timepoint = as.character(d$animal_timepoint), cluster = as.character(d[[args$cluster_col]]),
    method = m, p_val_adj = d$p_val_adj, fold_change = d[[fc]], stringsAsFactors = FALSE)
  if (anyNA(ans[, c("gene", "contrast", "timepoint", "cluster")])) stop("Missing identifiers: ", path)
  if (anyDuplicated(ans[, c("gene", "contrast", "timepoint", "cluster")])) stop("Duplicate gene/stratum: ", path)
  message(m, ": ", fc, " from ", path)
  ans
})
deg <- bind_rows(tables)
strata <- unique(deg[, c("contrast", "timepoint", "cluster")])
panels <- lapply(seq_len(nrow(strata)), function(i) {
  st <- strata[i, ]
  d <- deg[deg$contrast == st$contrast & deg$timepoint == st$timepoint & deg$cluster == st$cluster, ]
  genes <- select_panel_genes(d, methods, args$top_n, args$sig, args$genes)
  mat <- make_method_matrix(d, genes, methods)
  list(mat = mat, limit = scale_limit(mat))
})
# Pair male/female genotype panels for the same cell type and stage.
# P120 male and P150 female correspond to the late stage in main.R.
stage <- strata$timepoint
stage[stage %in% c("120", "150", "P120", "P150")] <- "late"
for (i in seq_len(nrow(strata))) {
  if (!strata$contrast[i] %in% c("genotype_male", "genotype_female")) next
  pair <- which(strata$contrast %in% c("genotype_male", "genotype_female") &
    stage == stage[i] & strata$cluster == strata$cluster[i])
  limits <- vapply(panels[pair], function(p) p$limit, numeric(1))
  if (any(is.finite(limits))) for (j in pair) panels[[j]]$limit <- min(limits, na.rm = TRUE)
}
out_dir <- file.path(args$fig_dir, args$cluster_col)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(strata))) {
  mat <- panels[[i]]$mat
  if (!nrow(mat) || !any(is.finite(mat))) next
  st <- strata[i, ]
  lim <- panels[[i]]$limit
  labels <- format_cells(mat)
  title <- unname(contrast_titles[st$contrast])
  if (is.na(title)) title <- st$contrast
  title <- paste0(title, "\n", st$cluster, " | ", timepoint_label(st$timepoint))
  mode <- if (!is.null(args$genes)) "custom" else paste0(if (args$sig) "sig_" else "", "top", args$top_n)
  path <- file.path(out_dir, paste0("logFC_heatmap_", mode, "_", safe_filename(st$contrast),
    "_", safe_filename(st$timepoint), "_", safe_filename(st$cluster), ".pdf"))
  pdf(path, width = max(7, ncol(mat) * 1.2 + 3), height = max(4, nrow(mat) * 0.32 + 2.3))
  tryCatch({
    ht <- Heatmap(mat, name = "log2FC", col = colorRamp2(c(-lim, 0, lim),
      c("#00008B", "#FFFFFF", "#FF0000")), na_col = "#EEEEEE",
      cluster_rows = FALSE, cluster_columns = FALSE,
      row_names_side = "left", row_names_gp = gpar(fontsize = 11, fontface = "italic"),
      show_column_names = TRUE, column_names_side = "bottom", column_names_rot = 0,
      column_names_gp = gpar(fontsize = 10),
      rect_gp = gpar(col = "white", lwd = 1.5),
      column_title = title, column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      heatmap_legend_param = list(at = c(-lim, -lim/2, 0, lim/2, lim)),
      cell_fun = function(j, i, x, y, width, height, fill) {
        rgb <- col2rgb(fill)[, 1] / 255
        text_col <- if (sum(rgb * c(.299, .587, .114)) < .4) "white" else "black"
        grid.text(labels[i, j], x, y, gp = gpar(fontsize = 9, col = text_col))
      })
    draw(ht, heatmap_legend_side = "right")
  }, finally = dev.off())
  message("Saved: ", path)
}
