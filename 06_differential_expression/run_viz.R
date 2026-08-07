#!/usr/bin/env Rscript
# ==============================================================================
# run_viz.R (Visualization-only DEG script)
# Reads pre-computed DEG CSV produced by run_deg_cell_level.R and
# generates:
#   1. Volcano plots: bold-bordered panels, viridis colors
#   2. Heatmaps of top DEGs: one per sex, blue -> grey -> red color scale
#      (#2166AC / #EEEEEE / #B2182B), with timepoint/cluster headers
#
# Usage:
#   Rscript run_viz.R \
#     --method wilcox \
#     --cluster_col neuron_class
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---------------------------------------------------------------------------- #
# CLI arguments                                                                  #
# ---------------------------------------------------------------------------- #
parser <- ArgumentParser(description = "DEG visualization only (volcano + heatmap)")
parser$add_argument("--method", type = "character", required = TRUE,
  choices = c("wilcox", "MAST"),
  help = "Test method used for DEG (used to auto-construct file paths)")
parser$add_argument("--cluster_col", type = "character", default = "neuron_class",
  choices = c("neuron_class", "cell_type_concise"),
  metavar = "<col>",
  help = "Metadata column defining cluster identity [default: neuron_class]")
parser$add_argument("--deg_csv", type = "character", default = NULL,
  metavar = "<path>",
  help = "Explicit path to combined DEG CSV. Auto-resolved if not provided.")
parser$add_argument("--fig_dir", type = "character", default = NULL,
  metavar = "<path>",
  help = "Base directory for figures. Auto-resolved to 03_figures if not provided.")
parser$add_argument("--top_n", type = "integer", default = 5,
  metavar = "<int>",
  help = "Number of top genes (by adjusted p-value) per stratum for heatmap [default: 5]")
args <- parser$parse_args()

# ---------------------------------------------------------------------------- #
# Path resolution                                                                #
# ---------------------------------------------------------------------------- #
BASE_CSV_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg/cell_level_view"
BASE_FIG_DIR <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/cell_level_view"

# Resolve deg_csv
if (is.null(args$deg_csv)) {
  args$deg_csv <- file.path(
    BASE_CSV_DIR, 
    args$method, 
    args$cluster_col,
    sprintf("deg_all_%s_%s.csv", tolower(args$method), args$cluster_col)
  )
}

# Resolve fig_dir
base_fig <- if (is.null(args$fig_dir)) BASE_FIG_DIR else args$fig_dir
fig_dir <- file.path(base_fig, args$method, args$cluster_col)

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
message("========================================================")
message("  Method               : ", args$method)
message("  Cluster col          : ", args$cluster_col)
message("  Input CSV            : ", args$deg_csv)
message("  Figure output dir    : ", fig_dir)
message("========================================================")


# ---------------------------------------------------------------------------- #
# Load and prep DEG table                                                        #
# ---------------------------------------------------------------------------- #
message("Reading DEG table: ", args$deg_csv)
deg_percell <- read.csv(args$deg_csv, stringsAsFactors = FALSE)

# Validate expected columns
required_cols <- c("gene", "p_val_adj", "avg_log2FC", "animal_sex",
                   "animal_timepoint", args$cluster_col)
missing_cols <- setdiff(required_cols, colnames(deg_percell))
if (length(missing_cols) > 0) {
  stop("Input CSV is missing required columns: ", paste(missing_cols, collapse = ", "))
}

deg <- deg_percell %>%
  filter(!is.na(p_val_adj)) %>%
  mutate(
    neglog10padj = -log10(pmax(p_val_adj, 1e-300)),
    stratum = paste(animal_sex, animal_timepoint, .data[[args$cluster_col]], sep = " | "),
    sig = case_when(
      p_val_adj < 0.05 & avg_log2FC >  0.25 ~ "up in mut",
      p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "down in mut",
      TRUE ~ "ns"
    )
  )

top_n_genes <- deg %>%
  group_by(stratum) %>%
  slice_min(order_by = p_val_adj, n = args$top_n, with_ties = FALSE) %>%
  ungroup()

sig_levels <- c("down in mut", "ns", "up in mut")
sig_colors <- setNames(c("#2166AC", "grey80", "#B2182B"), sig_levels)

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
cluster_strip_color <- stepped_colors[22]

strip_labels <- function(grob) {
  labels <- character()
  if (inherits(grob, "text")) {
    labels <- c(labels, as.character(grob$label))
  }
  if (!is.null(grob$grobs)) {
    labels <- c(labels, unlist(lapply(grob$grobs, strip_labels), use.names = FALSE))
  }
  if (!is.null(grob$children)) {
    labels <- c(labels, unlist(lapply(grob$children, strip_labels), use.names = FALSE))
  }
  labels[nzchar(labels)]
}

recolor_strip_rects <- function(grob, color) {
  if (inherits(grob, "rect")) {
    grob$gp$fill <- adjustcolor(color, alpha.f = 0.10)
    grob$gp$col <- color
    grob$gp$lwd <- 1.2
  }
  if (!is.null(grob$grobs)) {
    grob$grobs <- lapply(grob$grobs, recolor_strip_rects, color = color)
  }
  if (!is.null(grob$children)) {
    grob$children <- do.call(grid::gList, lapply(grob$children, recolor_strip_rects, color = color))
  }
  grob
}

lookup_strip_color <- function(color_map, label) {
  if (is.null(color_map) || is.na(label)) {
    return(NULL)
  }
  color <- unname(color_map[as.character(label)])
  if (length(color) == 0 || is.na(color)) {
    return(NULL)
  }
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

bordered_theme <- theme_minimal(base_size = 8) +
  theme(
    strip.text       = element_text(size = 7, face = "bold"),
    strip.background = element_rect(fill = "white", color = "grey50", linewidth = 0.6),
    strip.placement  = "outside",
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing    = unit(0.6, "lines")
  )


# ---------------------------------------------------------------------------- #
# Volcano plots - one PDF per sex                                                #
# ---------------------------------------------------------------------------- #
message("Generating per-sex volcano plots ...")

for (sx in unique(deg$animal_sex)) {

  deg_sx <- deg %>% filter(animal_sex == sx) %>%
    mutate(
      sig = factor(sig, levels = sig_levels),
      timepoint_display = factor(
        timepoint_label(animal_timepoint),
        levels = timepoint_levels_for(animal_timepoint)
      ),
      cluster_display = factor(.data[[args$cluster_col]], levels = unique(.data[[args$cluster_col]]))
    )
  top_sx <- top_n_genes %>% filter(animal_sex == sx) %>%
    mutate(
      timepoint_display = factor(
        timepoint_label(animal_timepoint),
        levels = levels(deg_sx$timepoint_display)
      ),
      cluster_display = factor(.data[[args$cluster_col]], levels = levels(deg_sx$cluster_display))
    )

  n_tp    <- length(levels(droplevels(deg_sx$timepoint_display)))
  n_class <- length(levels(droplevels(deg_sx$cluster_display)))
  y_max_sx <- max(deg_sx$neglog10padj, na.rm = TRUE)
  timepoint_strip_colors <- timepoint_colors_map[levels(droplevels(deg_sx$timepoint_display))]
  if (any(is.na(timepoint_strip_colors))) {
    timepoint_strip_colors[is.na(timepoint_strip_colors)] <- stepped_colors[23]
  }
  cluster_strip_colors <- setNames(
    rep(cluster_strip_color, n_class),
    levels(droplevels(deg_sx$cluster_display))
  )

  p_volcano_base <- ggplot(deg_sx, aes(x = avg_log2FC, y = neglog10padj, color = sig)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_text_repel(
      data = top_sx, aes(label = gene), size = 2.3, max.overlaps = 15,
      color = "black", segment.size = 0.2
    ) +
    scale_color_manual(values = sig_colors, name = NULL, drop = FALSE) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", linewidth = 0.2, color = "grey50") +
    coord_cartesian(ylim = c(0, y_max_sx * 1.03)) +
    bordered_theme +
    labs(
      title = paste("DEG volcano plots -", format_sex_label(sx), "(genotype: mut vs wt)"),
      x     = "avg log2FC (mut vs wt)",
      y     = "-log10(adjusted p-value)"
    )

  p_volcano_timepoint_rows <- p_volcano_base +
    facet_grid(timepoint_display ~ cluster_display, scales = "free_x", switch = "y") +
    labs(subtitle = paste("per-cell; rows = timepoint, columns =", args$cluster_col))

  out_path <- file.path(fig_dir, sprintf("volcano_%s_rows_timepoint_cols_%s.pdf", sx, args$cluster_col))
  save_colored_facet_pdf(
    p_volcano_timepoint_rows,
    out_path,
    width = max(8, n_class * 3),
    height = max(5, n_tp * 3),
    top_colors = cluster_strip_colors,
    row_colors = timepoint_strip_colors
  )
  message("  Saved: ", out_path)

  p_volcano_cluster_rows <- p_volcano_base +
    facet_grid(cluster_display ~ timepoint_display, scales = "free_x", switch = "y") +
    labs(subtitle = paste("per-cell; rows =", args$cluster_col, ", columns = timepoint"))

  out_path <- file.path(fig_dir, sprintf("volcano_%s_rows_%s_cols_timepoint.pdf", sx, args$cluster_col))
  save_colored_facet_pdf(
    p_volcano_cluster_rows,
    out_path,
    width = max(8, n_tp * 3),
    height = max(5, n_class * 3),
    top_colors = timepoint_strip_colors,
    row_colors = cluster_strip_colors
  )
  message("  Saved: ", out_path)
}


# ---------------------------------------------------------------------------- #
# Heatmap color function - blue -> grey -> red                                   #
# ---------------------------------------------------------------------------- #
bgr_col_fun <- function(mat) {
  max_abs <- max(abs(mat), na.rm = TRUE)
  lim <- floor(max_abs)
  if (lim == 0) lim <- max_abs # fallback if all values < 1
  colorRamp2(
    c(-lim, 0, lim),
    c("#2166AC", "#EEEEEE", "#B2182B")   # blue -> grey -> red
  )
}


# ---------------------------------------------------------------------------- #
# Heatmaps - one per sex                                                         #
# ---------------------------------------------------------------------------- #
make_sex_heatmap <- function(sex_label, deg, top_n_genes, fig_dir, cluster_col) {

  deg_sx     <- deg %>% filter(animal_sex == sex_label)
  top_sx     <- top_n_genes %>% filter(animal_sex == sex_label)
  heat_genes <- unique(top_sx$gene)

  if (length(heat_genes) == 0) {
    message("  No genes to plot for sex=", sex_label, "; skipping heatmap.")
    return(invisible(NULL))
  }

  # ---- log2FC matrix --------------------------------------------------------
  wide <- deg_sx %>%
    filter(gene %in% heat_genes) %>%
    select(gene, stratum, avg_log2FC) %>%
    pivot_wider(names_from = stratum, values_from = avg_log2FC, values_fill = 0)

  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$gene

  # ---- adjusted p-value matrix ----------------------------------------------
  padj_wide <- deg_sx %>%
    filter(gene %in% heat_genes) %>%
    select(gene, stratum, p_val_adj) %>%
    pivot_wider(names_from = stratum, values_from = p_val_adj, values_fill = 1)

  padj_mat <- as.matrix(padj_wide[, -1])
  rownames(padj_mat) <- padj_wide$gene
  padj_mat <- padj_mat[rownames(mat), colnames(mat)]

  # ---- column metadata & annotation ----------------------------------------
  col_meta <- strsplit(colnames(mat), " \\| ") %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE)
  colnames(col_meta) <- c("sex", "timepoint", cluster_col)

  col_meta$sex_display <- format_sex_label(col_meta$sex)
  col_meta$timepoint_display <- timepoint_label(col_meta$timepoint)
  col_meta[[cluster_col]] <- as.character(col_meta[[cluster_col]])
  col_meta[[cluster_col]][is.na(col_meta[[cluster_col]]) | col_meta[[cluster_col]] == ""] <- "Unknown"
  col_meta$sex_display[is.na(col_meta$sex_display) | col_meta$sex_display == ""] <- "Unknown"
  col_meta$timepoint_display[is.na(col_meta$timepoint_display) | col_meta$timepoint_display == "P"] <- "Unknown"

  tp_levels    <- timepoint_levels_for(col_meta$timepoint)
  if ("Unknown" %in% col_meta$timepoint_display) {
    tp_levels <- c(tp_levels, "Unknown")
  }
  class_levels <- unique(col_meta[[cluster_col]])
  sex_levels   <- unique(col_meta$sex_display)

  col_order <- order(
    factor(col_meta$timepoint_display, levels = tp_levels),
    factor(col_meta[[cluster_col]], levels = class_levels)
  )
  mat <- mat[, col_order, drop = FALSE]
  padj_mat <- padj_mat[, col_order, drop = FALSE]
  col_meta <- col_meta[col_order, , drop = FALSE]
  
  tp_colors <- setNames(unname(timepoint_colors_map[tp_levels]), tp_levels)
  if (any(is.na(tp_colors))) {
    tp_colors[is.na(tp_colors)] <- "#8C8C8C"
  }

  # Polychrome 36 colors reversed (right to left)
  polychrome_36 <- c("#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", 
                     "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", 
                     "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", 
                     "#85660D", "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", 
                     "#F7E1A0", "#C075A6", "#782AB6", "#AAF400", "#BDCDFF", "#822E1C", 
                     "#B5EFB5", "#7ED7D1", "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", 
                     "#3B00FB")
  poly_cols <- rev(unname(polychrome_36))
  class_colors <- setNames(rep(poly_cols, length.out = length(class_levels)), class_levels)

  sex_colors <- setNames(unname(sex_colors_map[sex_levels]), sex_levels)
  if (any(is.na(sex_colors))) {
    sex_colors[is.na(sex_colors)] <- "#8C8C8C"
  }

  anno_df <- data.frame(
    `Cell Type` = factor(col_meta[[cluster_col]], levels = class_levels),
    `Time Point` = factor(col_meta$timepoint_display, levels = tp_levels),
    Sex = factor(col_meta$sex_display, levels = sex_levels),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  anno_col_list <- setNames(
    list(class_colors, tp_colors, sex_colors),
    colnames(anno_df)
  )

  col_anno <- HeatmapAnnotation(
    df = anno_df,
    col = anno_col_list,
    annotation_name_side = "left",
    simple_anno_size     = unit(0.35, "cm")
  )

  sig_mat <- matrix(
    ifelse(padj_mat < 0.05, "*", ""),
    nrow = nrow(mat), dimnames = dimnames(mat)
  )

  title_col <- lookup_strip_color(sex_colors_map, format_sex_label(sex_label))
  if (is.null(title_col)) {
    title_col <- "#333333"
  }

  # ---- draw heatmap ---------------------------------------------------------
  ht <- Heatmap(
    mat,
    name              = "avg log2FC\n(mut vs wt)",
    col               = bgr_col_fun(mat),        # <- blue / grey / red scale
    top_annotation    = col_anno,
    cluster_columns   = FALSE,
    show_column_names = FALSE,
    border            = title_col,
    row_names_gp      = gpar(fontsize = 7),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sig_mat[i, j], x, y, gp = gpar(fontsize = 9, col = "black"))
    },
    heatmap_legend_param = list(direction = "vertical"),
    column_title = sex_title(sex_label),
    column_title_gp = gpar(
      col = title_col,
      fontface = "bold",
      fontsize = 12
    )
  )

  out_path <- file.path(fig_dir, sprintf("heatmap_top_genes_%s.pdf", sex_label))
  pdf(out_path,
      width  = max(9, ncol(mat) * 0.35 + 3),
      height = max(6, nrow(mat) * 0.2  + 2))
  draw(ht,
       heatmap_legend_side    = "right",
       annotation_legend_side = "right",
       merge_legend           = TRUE)
  dev.off()
  message("  Saved: ", out_path)
}

message("Generating per-sex heatmaps ...")
for (sx in unique(deg$animal_sex)) {
  make_sex_heatmap(sx, deg, top_n_genes, fig_dir, args$cluster_col)
}

message("Done. Figures in: ", fig_dir)
