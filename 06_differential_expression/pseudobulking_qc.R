#!/usr/bin/env Rscript

library(Seurat)
library(dplyr)
library(ggplot2)

# Set path and load data
DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/sample_level_view/"
FILE_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg/sample_level_view/"
seu <- readRDS(file.path(DATA_DIR, "seu_metadata.rds"))

# Configuration
MIN_CELLS  <- 10
MIN_COUNTS <- 500

# 1. Generate pseduobulk samples --------------------------------------

# pseudobulk the counts based on donor-condition-celltype
pseu_broad <- as.data.frame(AggregateExpression(seu, assays = "RNA", return.seurat = F, group.by = c("sample", "neuron_class"))[["RNA"]])
pseu_fine <- as.data.frame(AggregateExpression(seu, assays = "RNA", return.seurat = F, group.by = c("sample", "cell_type_concise"))[["RNA"]])

# each 'cell' is a donor-condition-celltype pseudobulk profile
head(pseu_broad)
head(pseu_fine)


# 2. Generate pseduobulk QC metrics --------------------------------------
get_pseudobulk_qc <- function(seu_obj, pb_df, sample_col = "sample", group_col) {
  cell_counts <- as.data.frame(table(seu_obj[[sample_col]][,1], seu_obj[[group_col]][,1])) %>%
    filter(Freq > 0) %>%
    mutate(
      clean_sample = ifelse(grepl("^[0-9]", Var1), paste0("g", Var1), as.character(Var1)),
      clean_sample = gsub("_", "-", clean_sample),
      clean_group  = gsub("_", "-", Var2),
      pb_id = paste(clean_sample, clean_group, sep = "_")
    ) %>%
    rename(
      !!sample_col := Var1,
      !!group_col  := Var2,
      psbulk_cells = Freq
    )
  qc_df <- data.frame(
    pb_id = colnames(pb_df),
    psbulk_counts = colSums(pb_df)
  ) %>%
    left_join(cell_counts, by = "pb_id")
  return(qc_df)
}

qc_broad <- get_pseudobulk_qc(
  seu_obj = seu, 
  pb_df = pseu_broad, 
  sample_col = "sample", 
  group_col = "neuron_class"
)

qc_fine <- get_pseudobulk_qc(
  seu_obj = seu, 
  pb_df = pseu_fine, 
  sample_col = "sample", 
  group_col = "cell_type_concise"
)

head(qc_broad)
head(qc_fine)


# 3. Visualize pseduobulk QC metrics --------------------------------------
plot_pseudobulk_qc <- function(qc_df, group_col, min_cells = MIN_CELLS, min_counts = MIN_COUNTS, facet = TRUE) {
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
  return(p)
}

p <- plot_pseudobulk_qc(qc_fine, group_col = "cell_type_concise")
ggplot2::ggsave(file.path(FIG_DIR, "pseudobulk_qc_cell_type_concise.pdf"), plot = p, height = 8.5, width = 10)
p <- plot_pseudobulk_qc(qc_broad, group_col = "neuron_class")
ggplot2::ggsave(file.path(FIG_DIR, "pseudobulk_qc_neuron_class.pdf"), plot = p, height = 8.5, width = 10)


# 4. Filter pseduobulk QC samples --------------------------------------
filter_pseudobulk_samples <- function(pb_df, 
                                      qc_df, 
                                      group_col, 
                                      hue_col = NULL, 
                                      min_cells = MIN_CELLS, 
                                      min_counts = MIN_COUNTS,
                                      output_dir = ".") {
  
  # 1. Flag samples that pass quality control thresholds
  qc_df <- qc_df %>%
    mutate(pass_qc = psbulk_cells >= min_cells & psbulk_counts >= min_counts)
  
  # Extract valid sample IDs matching the thresholds
  valid_ids <- qc_df %>%
    filter(pass_qc) %>%
    pull(pb_id)
  
  # 2. Filter expression matrix and QC metadata
  pb_df_filtered <- pb_df[, valid_ids]
  qc_df_filtered <- qc_df %>% filter(pass_qc)
  
  # 3. Create a detailed per-cell-type summary log dataframe
  detailed_log <- qc_df %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n_samples_before  = n(),
      n_samples_passed  = sum(pass_qc),
      n_samples_filtered = sum(!pass_qc),
      min_cells_cutoff  = min_cells,
      min_counts_cutoff = min_counts,
      .groups = "drop"
    )
  
  # Add overall grand total summary row
  total_row <- data.frame(
    group_temp         = "TOTAL_ALL_GROUPS",
    n_samples_before   = ncol(pb_df),
    n_samples_passed   = ncol(pb_df_filtered),
    n_samples_filtered = ncol(pb_df) - ncol(pb_df_filtered),
    min_cells_cutoff   = min_cells,
    min_counts_cutoff  = min_counts
  )
  colnames(total_row)[1] <- group_col
  
  detailed_log <- bind_rows(detailed_log, total_row)
  
  # Export filtering summary to a CSV file
  csv_filename <- file.path(output_dir, paste0("pseudobulk_filter_", group_col, ".csv"))
  write.csv(detailed_log, file = csv_filename, row.names = FALSE)
  message(paste("Detailed cell-type filtering log saved to:", csv_filename))
  
  # Return filtered outputs
  return(list(
    filtered_counts = pb_df_filtered,
    filtered_qc     = qc_df_filtered
  ))
}
  
res_broad <- filter_pseudobulk_samples(
  pb_df     = pseu_broad,
  qc_df     = qc_broad,
  group_col = "neuron_class",
  hue_col   = "animal_genotype",
  output_dir = FILE_DIR # Saves to pseudobulk_filter_neuron_class.csv
)

# Run for fine resolution (cell_type_concise)
res_fine <- filter_pseudobulk_samples(
  pb_df     = pseu_fine,
  qc_df     = qc_fine,
  group_col = "cell_type_concise",
  hue_col   = "animal_genotype",
  output_dir = FILE_DIR # Saves to pseudobulk_filter_cell_type_concise.csv
)
