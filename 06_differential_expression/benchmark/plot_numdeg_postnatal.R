#!/usr/bin/env Rscript
# ==============================================================================
# plot_numdeg_postnatal.R
#
# Summarize significant DEG counts (padj < 0.05) across postnatal timepoints and
# plot male/female side-by-side line graphs for each DEG method and clustering
# resolution.
#
# Outputs:
#   /quobyte/lasallegrp/Osman/shenyu/04_results/deg/numdeg_postnatal/
#   /quobyte/lasallegrp/Osman/shenyu/03_figures/deg/numdeg_postnatal/
#
# The script reads per-stratum DEG CSV files, so incomplete long-running methods
# such as cell-level LimmaVoomCC can still be plotted from currently available
# strata. Re-run this script after more DEG CSVs finish to refresh the figures.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(grid)
})

BASE_RESULT_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg"
BASE_FIG_DIR    <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg"

OUT_RESULT_DIR <- file.path(BASE_RESULT_DIR, "numdeg_postnatal")
OUT_FIG_DIR    <- file.path(BASE_FIG_DIR, "numdeg_postnatal")
dir.create(OUT_RESULT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

CLUSTER_RESOLUTIONS <- c("neuron_class", "cell_type_concise")

METHODS <- tibble::tribble(
  ~analysis_level, ~method_dir,      ~method_label,
  "cell_level",    "wilcox",         "Wilcoxon",
  "cell_level",    "MAST",           "MAST",
  "sample_level",  "limma_voomcc",   "LimmaVoomCC_sample",
  "cell_level",    "LimmaVoomCC",    "LimmaVoomCC_cell",
  "sample_level",  "deseq2",         "DESeq2",
  "sample_level",  "edger_zinger",   "EdgeR-zingeR"
)

ditto_seq_colors <- c(
  "17" = "#FFBE2D",
  "18" = "#80C7EF",
  "19" = "#00F6B3",
  "21" = "#06A5FF",
  "22" = "#FF8320"
)

sex_colors_map <- c("Male" = ditto_seq_colors[["21"]], "Female" = ditto_seq_colors[["22"]])

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

format_sex_label <- function(x) {
  dplyr::case_when(
    tolower(as.character(x)) %in% c("male", "m") ~ "Male",
    tolower(as.character(x)) %in% c("female", "f") ~ "Female",
    TRUE ~ as.character(x)
  )
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

result_dir_for <- function(analysis_level, method_dir, cluster_col) {
  view_dir <- switch(
    analysis_level,
    cell_level = "cell_level_view",
    sample_level = "sample_level_view",
    stop("Unsupported analysis level: ", analysis_level)
  )
  file.path(BASE_RESULT_DIR, view_dir, method_dir, cluster_col)
}

read_method_deg <- function(analysis_level, method_dir, cluster_col) {
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
      mutate(source_file = basename(path))
  }))
}

summarize_numdeg <- function(deg, analysis_level, method_dir, method_label, cluster_col) {
  if (nrow(deg) == 0) return(tibble())

  required_cols <- c("animal_sex", "animal_timepoint", cluster_col, "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(deg))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s) for ", analysis_level, " / ", method_dir,
      " / ", cluster_col, ": ", paste(missing_cols, collapse = ", ")
    )
  }

  deg %>%
    mutate(
      analysis_level = analysis_level,
      method = method_dir,
      method_label = method_label,
      cluster_col = cluster_col,
      cluster_id = as.character(.data[[cluster_col]]),
      sex_display = format_sex_label(animal_sex),
      age_days = timepoint_num(animal_timepoint),
      is_sig = !is.na(p_val_adj) & p_val_adj < 0.05
    ) %>%
    group_by(analysis_level, method, method_label, cluster_col, sex_display, age_days, cluster_id) %>%
    summarise(numDEG = sum(is_sig), n_genes_tested = n(), .groups = "drop")
}

add_sex_border_layers <- function(plot) {
  for (sex_label in names(sex_colors_map)) {
    border_df <- data.frame(
      sex_display = sex_label,
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    )
    plot <- plot +
      geom_rect(
        data = border_df,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = sex_colors_map[[sex_label]],
        linewidth = 1.2
      )
  }
  plot
}

plot_numdeg_one <- function(numdeg, analysis_level, method_dir, method_label, cluster_col) {
  cluster_levels <- cluster_levels_for(cluster_col, unique(numdeg$cluster_id))
  cluster_colors <- cluster_colors_for(cluster_col, cluster_levels)

  plot_data <- numdeg %>%
    filter(
      analysis_level == .env$analysis_level,
      method == .env$method_dir,
      cluster_col == .env$cluster_col
    ) %>%
    mutate(
      cluster_id = factor(cluster_id, levels = cluster_levels),
      sex_display = factor(sex_display, levels = c("Male", "Female")),
      age_label = paste0("P", age_days)
    )

  plot_df <- plot_data %>%
    complete(
      sex_display,
      age_days,
      cluster_id,
      fill = list(numDEG = NA_integer_, n_genes_tested = NA_integer_)
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = age_days, y = numDEG, color = cluster_id, group = cluster_id)) +
    geom_line(data = plot_data, linewidth = 0.8) +
    geom_point(data = plot_data, size = 2) +
    facet_wrap(~sex_display, nrow = 1, scales = "free_y", drop = FALSE) +
    scale_color_manual(values = cluster_colors, limits = cluster_levels, drop = FALSE, name = "Cell Type") +
    scale_x_continuous(
      breaks = sort(unique(plot_df$age_days[!is.na(plot_df$age_days)])),
      labels = paste0("P", sort(unique(plot_df$age_days[!is.na(plot_df$age_days)])))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
      panel.border = element_blank(),
      strip.background = element_rect(fill = "white", color = "grey50", linewidth = 0.6),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0)
    ) +
    labs(
      title = paste(method_label, cluster_col, "numDEGs", sep = " - "),
      x = "Age",
      y = "numDEGs"
    )

  p <- add_sex_border_layers(p)

  out_dir <- file.path(OUT_FIG_DIR, analysis_level, cluster_col)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- paste("numdeg", analysis_level, method_dir, cluster_col, sep = "_")
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))

  ggsave(pdf_path, plot = p, width = 10.5, height = 4.2, useDingbats = FALSE)
  ggsave(png_path, plot = p, width = 10.5, height = 4.2, dpi = 300)
  message("Saved: ", pdf_path)
  message("Saved: ", png_path)

  invisible(p)
}

plot_numdeg_cluster_combined <- function(numdeg, cluster_col) {
  numdeg_cluster <- numdeg %>%
    filter(cluster_col == .env$cluster_col)
  if (nrow(numdeg_cluster) == 0) return(invisible(NULL))

  method_panels <- METHODS %>%
    mutate(
      method_key = paste(analysis_level, method_dir, sep = "__"),
      method_order = row_number()
    )

  cluster_levels <- cluster_levels_for(cluster_col, unique(numdeg_cluster$cluster_id))
  cluster_colors <- cluster_colors_for(cluster_col, cluster_levels)
  sex_levels <- c("Male", "Female")
  age_levels <- sort(unique(numdeg_cluster$age_days[!is.na(numdeg_cluster$age_days)]))

  panel_meta <- tidyr::expand_grid(
    method_key = method_panels$method_key,
    sex_display = sex_levels
  ) %>%
    left_join(method_panels, by = "method_key") %>%
    mutate(
      sex_display = factor(sex_display, levels = sex_levels),
      panel_title = paste(method_label, "Postnatal", tolower(as.character(sex_display)), "numDEGs")
    )

  panel_levels <- panel_meta %>%
    arrange(method_order, sex_display) %>%
    pull(panel_title) %>%
    as.character()

  plot_data <- numdeg_cluster %>%
    mutate(method_key = paste(analysis_level, method, sep = "__")) %>%
    select(method_key, sex_display, age_days, cluster_id, numDEG, n_genes_tested) %>%
    left_join(panel_meta %>% select(method_key, sex_display, method_order, method_label, panel_title),
      by = c("method_key", "sex_display")
    ) %>%
    mutate(
      cluster_id = factor(cluster_id, levels = cluster_levels),
      sex_display = factor(sex_display, levels = sex_levels),
      panel_title = factor(as.character(panel_title), levels = panel_levels)
    )

  plot_df <- numdeg_cluster %>%
    mutate(method_key = paste(analysis_level, method, sep = "__")) %>%
    select(method_key, sex_display, age_days, cluster_id, numDEG, n_genes_tested) %>%
    right_join(
      tidyr::expand_grid(
        method_key = method_panels$method_key,
        sex_display = sex_levels,
        age_days = age_levels,
        cluster_id = cluster_levels
      ),
      by = c("method_key", "sex_display", "age_days", "cluster_id")
    ) %>%
    left_join(panel_meta %>% select(method_key, sex_display, method_order, method_label, panel_title),
      by = c("method_key", "sex_display")
    ) %>%
    mutate(
      cluster_id = factor(cluster_id, levels = cluster_levels),
      sex_display = factor(sex_display, levels = sex_levels),
      panel_title = factor(as.character(panel_title), levels = panel_levels)
    )

  border_df <- panel_meta %>%
    mutate(
      panel_title = factor(as.character(panel_title), levels = panel_levels),
      border_color = unname(sex_colors_map[as.character(sex_display)]),
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    )

  p <- ggplot(plot_df, aes(x = age_days, y = numDEG, color = cluster_id, group = cluster_id)) +
    geom_line(data = plot_data, linewidth = 0.75) +
    geom_point(data = plot_data, size = 1.7) +
    facet_wrap(~panel_title, ncol = 2, scales = "free_y", drop = FALSE) +
    scale_color_manual(values = cluster_colors, limits = cluster_levels, drop = FALSE, name = "Cell Type") +
    scale_x_continuous(
      breaks = age_levels,
      labels = paste0("P", age_levels)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.10))) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.35, "cm"),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
      panel.border = element_blank(),
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(2.0, "lines"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", hjust = 0, size = 13)
    ) +
    labs(
      title = paste("Postnatal numDEGs by method -", cluster_col),
      x = "Age",
      y = "numDEGs"
    )

  for (sex_label in names(sex_colors_map)) {
    p <- p +
      geom_rect(
        data = border_df %>% filter(as.character(sex_display) == sex_label),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        inherit.aes = FALSE,
        fill = NA,
        color = sex_colors_map[[sex_label]],
        linewidth = 1.1
      )
  }

  out_dir <- file.path(OUT_FIG_DIR, "combined")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(out_dir, paste0("numdeg_combined_", cluster_col, ".pdf"))
  png_path <- file.path(out_dir, paste0("numdeg_combined_", cluster_col, ".png"))

  ggsave(pdf_path, plot = p, width = 13, height = 18, useDingbats = FALSE)
  ggsave(png_path, plot = p, width = 13, height = 18, dpi = 300)
  message("Saved combined figure: ", pdf_path)
  message("Saved combined figure: ", png_path)

  invisible(p)
}

all_numdeg <- list()

for (cluster_col in CLUSTER_RESOLUTIONS) {
  for (i in seq_len(nrow(METHODS))) {
    method_info <- METHODS[i, ]
    message("Reading ", method_info$analysis_level, " / ", method_info$method_dir, " / ", cluster_col)

    deg <- read_method_deg(
      analysis_level = method_info$analysis_level,
      method_dir = method_info$method_dir,
      cluster_col = cluster_col
    )

    summary_i <- summarize_numdeg(
      deg = deg,
      analysis_level = method_info$analysis_level,
      method_dir = method_info$method_dir,
      method_label = method_info$method_label,
      cluster_col = cluster_col
    )

    if (nrow(summary_i) == 0) {
      message("No plottable DEG rows for ", method_info$analysis_level, " / ", method_info$method_dir, " / ", cluster_col)
      next
    }

    all_numdeg[[paste(cluster_col, method_info$analysis_level, method_info$method_dir, sep = "|")]] <- summary_i
  }
}

numdeg_summary <- bind_rows(all_numdeg)
summary_csv <- file.path(OUT_RESULT_DIR, "numdeg_postnatal_summary.csv")
write.csv(numdeg_summary, summary_csv, row.names = FALSE)
message("Summary saved: ", summary_csv)

if (nrow(numdeg_summary) > 0) {
  for (cluster_col in CLUSTER_RESOLUTIONS) {
    plot_numdeg_cluster_combined(numdeg_summary, cluster_col)

    for (i in seq_len(nrow(METHODS))) {
      method_info <- METHODS[i, ]
      plot_numdeg_one(
        numdeg = numdeg_summary,
        analysis_level = method_info$analysis_level,
        method_dir = method_info$method_dir,
        method_label = method_info$method_label,
        cluster_col = cluster_col
      )
    }
  }
}

message("Done.")
