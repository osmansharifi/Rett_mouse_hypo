#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

args <- commandArgs(trailingOnly = TRUE)
input_rds <- if (length(args) >= 1) args[[1]] else "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects/seu_final.rds"
output_pdf <- if (length(args) >= 2) args[[2]] else "output/pdf/propCells_cell_type_concise.pdf"

if (!file.exists(input_rds)) stop("Input RDS does not exist: ", input_rds)
dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# This is the same cell-type mapping and order used by the DEG heatmaps.
cell_type_colors <- c(
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

cell_type_levels <- c(
  "GABA", "GLU", "Astrocytes", "Ependymal-like", "Immune",
  "Oligodendrocytes", "OPC", "Fibroblasts", "Mural+Endothelial", "ParsTuber"
)

pick_col <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

obj <- readRDS(input_rds)
md <- obj[[]]
rm(obj)
invisible(gc())

cell_col <- pick_col(md, "cell_type_concise", "cell_type_concise")
sample_col <- pick_col(md, c("sample", "sample_name", "orig.ident", "mouse"), "sample ID")
sex_col <- pick_col(md, c("animal_sex", "sex"), "sex")
genotype_col <- pick_col(md, c("animal_genotype", "genotype"), "genotype")
time_col <- pick_col(md, c("animal_timepoint", "timepoint", "age"), "timepoint")

plot_md <- md %>%
  transmute(
    sample_id = as.character(.data[[sample_col]]),
    cell_type = as.character(.data[[cell_col]]),
    sex = tolower(as.character(.data[[sex_col]])),
    genotype = tolower(as.character(.data[[genotype_col]])),
    age_days = as.numeric(gsub("[^0-9.]", "", as.character(.data[[time_col]])))
  ) %>%
  filter(
    !is.na(sample_id), !is.na(cell_type), !is.na(age_days),
    sex %in% c("male", "female"), genotype %in% c("wt", "mut")
  )

unknown_types <- setdiff(sort(unique(plot_md$cell_type)), names(cell_type_colors))
if (length(unknown_types)) {
  stop("cell_type_concise values missing from the heatmap mapping: ",
       paste(unknown_types, collapse = ", "))
}

# Compute proportions within each biological sample, then average samples in
# each sex x genotype x age group. Completing absent types with zero prevents
# accidental denominator or line omissions.
sample_meta <- plot_md %>%
  distinct(sample_id, sex, genotype, age_days) %>%
  count(sample_id, name = "n_metadata_rows")
if (any(sample_meta$n_metadata_rows != 1L)) {
  stop("At least one sample maps to multiple sex/genotype/age combinations.")
}

sample_props <- plot_md %>%
  count(sample_id, sex, genotype, age_days, cell_type, name = "n_cells") %>%
  group_by(sample_id, sex, genotype, age_days) %>%
  mutate(propCells = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  tidyr::complete(
    tidyr::nesting(sample_id, sex, genotype, age_days),
    cell_type = cell_type_levels,
    fill = list(n_cells = 0, propCells = 0)
  )

plot_df <- sample_props %>%
  group_by(sex, genotype, age_days, cell_type) %>%
  summarise(propCells = mean(propCells), n_samples = n(), .groups = "drop") %>%
  mutate(cell_type = factor(cell_type, levels = cell_type_levels))

panel_specs <- tibble::tribble(
  ~sex,     ~genotype, ~title,                    ~tag,
  "male",   "wt",      "propCells_WT_males",     "a",
  "male",   "mut",     "propCells_MUT_males",    "b",
  "female", "wt",      "propCells_WT_females",   "c",
  "female", "mut",     "propCells_MUT_females",  "d"
)

make_panel <- function(sex_use, genotype_use, title_use, show_legend = FALSE) {
  dat <- filter(plot_df, sex == sex_use, genotype == genotype_use)
  if (!nrow(dat)) stop("No data for ", genotype_use, " / ", sex_use)

  x_breaks <- sort(unique(dat$age_days))
  # Same sex-color mapping used by the heatmap annotations/titles.
  border_col <- if (sex_use == "male") "#06A5FF" else "#FF8320"

  ggplot(dat, aes(age_days, propCells, color = cell_type, group = cell_type)) +
    geom_line(linewidth = 1.5, lineend = "round") +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = cell_type_colors,
      limits = cell_type_levels,
      breaks = cell_type_levels,
      drop = FALSE,
      name = "Cell Type"
    ) +
    scale_x_continuous(breaks = x_breaks, minor_breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(title = title_use, x = "Age (days)", y = "propCells") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold", color = "#222222"),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = border_col, fill = NA, linewidth = 1.2),
      axis.line = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = grid::unit(0.38, "cm"),
      legend.key.width = grid::unit(0.55, "cm"),
      legend.title = element_text(face = "plain"),
      legend.text = element_text(size = 8.5),
      plot.margin = margin(8, 8, 8, 8)
    )
}

panels <- lapply(seq_len(nrow(panel_specs)), function(i) {
  make_panel(panel_specs$sex[i], panel_specs$genotype[i], panel_specs$title[i], show_legend = FALSE)
})

legend_plot <- make_panel("male", "wt", "", show_legend = TRUE)
legend_grob <- cowplot::get_legend(legend_plot)

tagged <- Map(
  function(p, tag) cowplot::ggdraw(p) +
    cowplot::draw_label(tag, x = 0, y = 1, hjust = -0.15, vjust = 1.25,
                        fontface = "bold", size = 13),
  panels, panel_specs$tag
)

panel_grid <- cowplot::plot_grid(plotlist = tagged, ncol = 2, align = "hv")
final_plot <- cowplot::plot_grid(panel_grid, legend_grob, nrow = 1, rel_widths = c(1, 0.22))

ggsave(output_pdf, final_plot, width = 14.5, height = 8.5, units = "in", device = cairo_pdf)
message("Saved: ", normalizePath(output_pdf, mustWork = FALSE))
