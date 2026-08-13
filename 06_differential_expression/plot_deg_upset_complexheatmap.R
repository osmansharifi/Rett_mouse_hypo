#!/usr/bin/env Rscript
# ==============================================================================
# plot_deg_upset_complexheatmap.R
#
# ComplexHeatmap UpSet plots for overlaps among significant DEGs (padj < 0.05).
#
# For each clustering resolution, this script writes one multi-page PDF. Each page
# is one cell type x postnatal timepoint, with male and female UpSet plots shown
# side by side. Re-run after long DEG jobs finish to refresh incomplete methods.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(grid)
})

BASE_RESULT_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg"
BASE_FIG_DIR    <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg"

OUT_RESULT_DIR <- file.path(BASE_RESULT_DIR, "deg_upset")
OUT_FIG_DIR    <- file.path(BASE_FIG_DIR, "deg_upset")
dir.create(OUT_RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

CLUSTER_RESOLUTIONS <- c("neuron_class", "cell_type_concise")
PADJ_CUTOFF <- 0.05

# Algorithm order in the UpSet rows and set-size bars.
METHODS <- tibble::tribble(
  ~method_label,           ~analysis_level, ~method_dir,
  "Wilcoxon",              "cell_level",    "wilcox",
  "MAST",                  "cell_level",    "MAST",
  "LimmaVoomCC_cell",      "cell_level",    "LimmaVoomCC",
  "LimmaVoomCC_sample",    "sample_level",  "limma_voomcc",
  "DESeq2",                "sample_level",  "deseq2",
  "EdgeR-zingeR",          "sample_level",  "edger_zinger"
)

METHOD_LEVELS <- METHODS$method_label

method_colors <- c(
  "Wilcoxon" = "#3BA6D0",
  "MAST" = "#2B6FAE",
  "LimmaVoomCC_sample" = "#3BA6D0",
  "LimmaVoomCC_cell" = "#2B6FAE",
  "DESeq2" = "#3BA6D0",
  "EdgeR-zingeR" = "#2B6FAE"
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

neuron_class_levels <- c("GABA", "GLU", "Non-neuronal")

format_sex_label <- function(x) {
  dplyr::case_when(
    tolower(as.character(x)) %in% c("male", "m") ~ "Male",
    tolower(as.character(x)) %in% c("female", "f") ~ "Female",
    TRUE ~ as.character(x)
  )
}

timepoint_label <- function(x) {
  paste0("P", gsub("[^0-9]", "", as.character(x)))
}

timepoint_num <- function(x) {
  as.integer(gsub("[^0-9]", "", as.character(x)))
}

cluster_levels_for <- function(cluster_col, values) {
  preferred <- switch(
    cluster_col,
    cell_type_concise = cell_type_concise_levels,
    neuron_class = neuron_class_levels,
    character()
  )
  values <- as.character(values)
  c(preferred[preferred %in% values], setdiff(sort(unique(values)), preferred))
}

result_dir_for <- function(analysis_level, method_dir, cluster_col) {
  view_dir <- switch(
    analysis_level,
    cell_level = "cell_level_view",
    sample_level = "sample_level_view",
    stop("Unsupported analysis level: ", analysis_level)
  )
  file.path(BASE_RESULT_DIR, view_dir, method_dir, cluster_col)
}

read_method_deg <- function(method_label, analysis_level, method_dir, cluster_col) {
  in_dir <- result_dir_for(analysis_level, method_dir, cluster_col)
  if (!dir.exists(in_dir)) {
    message("Missing result directory: ", in_dir)
    return(tibble())
  }

  csv_files <- list.files(in_dir, pattern = "^deg_.*\\.csv$", full.names = TRUE)
  csv_files <- csv_files[!grepl("^deg_all_", basename(csv_files))]
  if (length(csv_files) == 0) {
    message("No per-stratum DEG CSV files found: ", in_dir)
    return(tibble())
  }

  bind_rows(lapply(csv_files, function(path) {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE) %>%
      mutate(
        method_label = method_label,
        analysis_level = analysis_level,
        method_dir = method_dir,
        source_file = basename(path)
      )
  }))
}

read_cluster_deg <- function(cluster_col) {
  bind_rows(lapply(seq_len(nrow(METHODS)), function(i) {
    read_method_deg(
      method_label = METHODS$method_label[i],
      analysis_level = METHODS$analysis_level[i],
      method_dir = METHODS$method_dir[i],
      cluster_col = cluster_col
    )
  }))
}

make_sig_gene_sets <- function(sig_df, cluster_col, cluster_id, sex_display, timepoint) {
  sets <- lapply(METHOD_LEVELS, function(method_i) {
    sig_df %>%
      filter(
        method_label == method_i,
        .data[[cluster_col]] == cluster_id,
        sex_display == .env$sex_display,
        animal_timepoint_num == .env$timepoint
      ) %>%
      pull(gene) %>%
      unique() %>%
      as.character()
  })
  names(sets) <- METHOD_LEVELS
  sets
}

draw_empty_panel <- function(title, subtitle = "No significant DEGs available") {
  grid.rect(gp = gpar(fill = "white", col = NA))
  grid.text(title, x = unit(0.5, "npc"), y = unit(0.72, "npc"),
            gp = gpar(fontsize = 13, fontface = "bold"))
  grid.text(subtitle, x = unit(0.5, "npc"), y = unit(0.52, "npc"),
            gp = gpar(fontsize = 10, col = "grey40"))
}

draw_upset_panel <- function(gene_sets, title) {
  set_sizes <- vapply(gene_sets, length, integer(1))
  if (sum(set_sizes) == 0) {
    draw_empty_panel(title)
    return(invisible(NULL))
  }

  comb_mat <- ComplexHeatmap::make_comb_mat(gene_sets, mode = "intersect")
  comb_order <- order(ComplexHeatmap::comb_size(comb_mat), decreasing = TRUE)

  top_anno <- ComplexHeatmap::upset_top_annotation(
    comb_mat,
    add_numbers = TRUE,
    numbers_gp = gpar(fontsize = 10),
    gp = gpar(fill = "#333333"),
    height = unit(3.5, "cm"),
    axis_param = list(gp = gpar(fontsize = 9))
  )

  left_anno <- ComplexHeatmap::upset_left_annotation(
    comb_mat,
    add_numbers = TRUE,
    numbers_gp = gpar(fontsize = 10),
    gp = gpar(fill = unname(method_colors[METHOD_LEVELS])),
    width = unit(3.5, "cm"),
    axis_param = list(gp = gpar(fontsize = 9))
  )

  ht <- ComplexHeatmap::UpSet(
    comb_mat,
    comb_order = comb_order,
    set_order = METHOD_LEVELS,
    top_annotation = top_anno,
    left_annotation = left_anno,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 11),
    pt_size = unit(2.4, "mm"),
    lwd = 0.8,
    bg_col = "#F5DDDB",
    bg_pt_col = "#E0D6D1",
    comb_col = "#3A3A3A",
    column_title = title,
    column_title_gp = gpar(fontsize = 15, fontface = "bold")
  )

  ComplexHeatmap::draw(ht, newpage = FALSE)
  grid.text(
    "Intersection numDEGs",
    x = unit(0.03, "npc"),
    y = unit(0.78, "npc"),
    rot = 90,
    gp = gpar(fontsize = 11)
  )
}

plot_cluster_upset_pdf <- function(cluster_col, deg) {
  required_cols <- c("gene", "animal_sex", "animal_timepoint", cluster_col, "p_val_adj", "method_label")
  missing_cols <- setdiff(required_cols, colnames(deg))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) for ", cluster_col, ": ", paste(missing_cols, collapse = ", "))
  }

  sig_df <- deg %>%
    filter(!is.na(p_val_adj), p_val_adj < PADJ_CUTOFF) %>%
    mutate(
      sex_display = format_sex_label(animal_sex),
      animal_timepoint_num = timepoint_num(animal_timepoint),
      method_label = factor(method_label, levels = METHOD_LEVELS)
    )

  all_strata <- deg %>%
    mutate(animal_timepoint_num = timepoint_num(animal_timepoint)) %>%
    distinct(.data[[cluster_col]], animal_timepoint_num) %>%
    rename(cluster_id = all_of(cluster_col)) %>%
    arrange(
      factor(cluster_id, levels = cluster_levels_for(cluster_col, cluster_id)),
      animal_timepoint_num
    )

  if (nrow(all_strata) == 0) {
    message("No strata found for ", cluster_col)
    return(invisible(NULL))
  }

  out_pdf <- file.path(OUT_FIG_DIR, paste0("upset_numdeg_", cluster_col, ".pdf"))
  pdf(out_pdf, width = 16, height = 8.5, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)

  for (i in seq_len(nrow(all_strata))) {
    cluster_id <- as.character(all_strata$cluster_id[i])
    tp <- all_strata$animal_timepoint_num[i]

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(
      nrow = 2,
      ncol = 2,
      heights = unit.c(unit(0.55, "in"), unit(1, "null")),
      widths = unit.c(unit(1, "null"), unit(1, "null"))
    )))

    title <- paste(cluster_id, timepoint_label(tp), "numDEGs")
    grid.text(
      title,
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),
      gp = gpar(fontsize = 18, fontface = "bold")
    )

    for (j in seq_along(c("Male", "Female"))) {
      sex_i <- c("Male", "Female")[j]
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = j))
      panel_title <- paste(cluster_id, timepoint_label(tp), tolower(sex_i), "numDEGs")
      gene_sets <- make_sig_gene_sets(sig_df, cluster_col, cluster_id, sex_i, tp)
      draw_upset_panel(gene_sets, panel_title)
      upViewport()
    }

    upViewport()
  }

  message("Saved: ", out_pdf)
  invisible(out_pdf)
}

all_summaries <- list()

for (cluster_col in CLUSTER_RESOLUTIONS) {
  message("Reading DEG results for ", cluster_col)
  deg <- read_cluster_deg(cluster_col)

  if (nrow(deg) == 0) {
    message("No DEG rows available for ", cluster_col)
    next
  }

  summary_df <- deg %>%
    mutate(
      cluster_col = cluster_col,
      cluster_id = as.character(.data[[cluster_col]]),
      sex_display = format_sex_label(animal_sex),
      animal_timepoint_num = timepoint_num(animal_timepoint),
      is_sig = !is.na(p_val_adj) & p_val_adj < PADJ_CUTOFF
    ) %>%
    group_by(cluster_col, cluster_id, sex_display, animal_timepoint_num, method_label) %>%
    summarise(numDEGs = sum(is_sig), n_genes_tested = n(), .groups = "drop")

  all_summaries[[cluster_col]] <- summary_df
  plot_cluster_upset_pdf(cluster_col, deg)
}

summary_all <- bind_rows(all_summaries)
summary_csv <- file.path(OUT_RESULT_DIR, "upset_numdeg_set_sizes.csv")
write.csv(summary_all, summary_csv, row.names = FALSE)
message("Summary saved: ", summary_csv)
message("Done.")
