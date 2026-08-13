#!/usr/bin/env Rscript
# ==============================================================================
# run_deg_cell_level.R
#
# Cell-level DEG analysis.
# Supports: wilcox | MAST | LimmaVoomCC
# Loops over every animal_sex x animal_timepoint x <cluster_col> stratum.
# Comparison: mut (ident.1) vs wt (ident.2)
#
# Usage:
#   Rscript run_deg_cell_level.R \
#     --method      wilcox \
#     --cluster_col neuron_class \
#     [--seu        /path/to/seu_final.rds] \
#     [--file_dir   /path/to/output/] \
#     [--min_pct    0.1] \
#     [--min_cells  10]
#
#   Rscript run_deg_cell_level.R --method MAST --cluster_col cell_type_concise
#
# Output CSVs and figures are generated in the same run.
# ==============================================================================
suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(grid)
})
# ---------------------------------------------------------------------------- #
# CLI arguments                                                                  #
# ---------------------------------------------------------------------------- #
parser <- ArgumentParser(
  description = "Cell-level DEG via FindMarkers or LimmaVoomCC"
)
parser$add_argument("--method",
  type = "character", required = TRUE,
  choices = c("wilcox", "MAST", "LimmaVoomCC"),
  help = "DEG method"
)
parser$add_argument("--cluster_col",
  type = "character", required = TRUE,
  choices = c("neuron_class", "cell_type_concise"),
  help = "Metadata column defining cluster identity for stratification"
)
parser$add_argument("--seu",
  type = "character",
  default = "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects/seu_final.rds",
  metavar = "<path>",
  help = "Path to Seurat RDS object [default: quobyte path]"
)
parser$add_argument("--file_dir",
  type = "character",
  default = NULL, # auto-resolved below if not supplied
  metavar = "<path>",
  help = "Base output directory for DEG CSV files (method/cluster_col subfolders are always appended)"
)
parser$add_argument("--fig_dir",
  type = "character",
  default = NULL,
  metavar = "<path>",
  help = "Base output directory for DEG figures (method/cluster_col subfolders are always appended)"
)
parser$add_argument("--min_pct",
  type = "double", default = 0.25,
  metavar = "<float>",
  help = "min.pct threshold for FindMarkers [default: 0.25]"
)
parser$add_argument("--min_cells",
  type = "integer", default = 10,
  metavar = "<int>",
  help = "Minimum cells per genotype group per stratum [default: 10]"
)
parser$add_argument("--viz_only",
  action = "store_true", default = FALSE,
  help = "Skip DEG fitting and regenerate figures from the combined DEG CSV"
)
args <- parser$parse_args()
# ---------------------------------------------------------------------------- #
# Path resolution                                                                #
# Folder structure (ALWAYS applied, even if --file_dir is supplied manually):    #
#   <file_dir base>/                                                             #
#   └── <method>/          e.g. wilcox/ or MAST/                                #
#       └── <cluster_col>/ e.g. neuron_class/ or cell_type_concise/             #
#           ├── deg_<sex>_<tp>_<cluster>.csv   (per-stratum)                    #
#           └── deg_all_<method>_<cluster_col>.csv  (combined)                  #
# ---------------------------------------------------------------------------- #
BASE_FILE_DIR <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg/cell_level_view"
BASE_FIG_DIR <- "/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/cell_level_view"
base_dir <- if (is.null(args$file_dir)) BASE_FILE_DIR else args$file_dir
base_fig_dir <- if (is.null(args$fig_dir)) BASE_FIG_DIR else args$fig_dir
args$file_dir <- file.path(
  base_dir,
  args$method, # top-level: method name  (e.g. wilcox, MAST)
  args$cluster_col # subfolder: clustering granularity
)
args$fig_dir <- file.path(
  base_fig_dir,
  args$method,
  args$cluster_col
)
dir.create(args$file_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(args$fig_dir, showWarnings = FALSE, recursive = TRUE)
message("========================================================")
message("  Method               : ", args$method)
message("  Cluster col          : ", args$cluster_col)
message("  min.pct              : ", args$min_pct)
message("  min_cells            : ", args$min_cells)
message("  Output dir           : ", args$file_dir)
message("  Figure dir           : ", args$fig_dir)
message("========================================================")
# ---------------------------------------------------------------------------- #
# Helper                                                                         #
# ---------------------------------------------------------------------------- #
`%||%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a
min_cells_per_group <- args$min_cells # minimum cells per genotype group within a stratum

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "unknown", x)
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

run_limmavoomcc_deg <- function(sub, cluster_col, cluster_i, sex_i, tp_i, tab) {
  if (!requireNamespace("edgeR", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) {
    stop("LimmaVoomCC requires edgeR and limma.")
  }

  counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
  expressed_mut <- Matrix::rowMeans(counts[, sub$animal_genotype == "mut", drop = FALSE] > 0)
  expressed_wt <- Matrix::rowMeans(counts[, sub$animal_genotype == "wt", drop = FALSE] > 0)
  keep_genes <- pmax(expressed_mut, expressed_wt) >= args$min_pct
  counts <- counts[keep_genes & Matrix::rowSums(counts) > 0, , drop = FALSE]

  meta <- sub[[]]
  meta$animal_genotype <- factor(meta$animal_genotype, levels = c("wt", "mut"))
  design <- model.matrix(~ animal_genotype, data = meta)
  coef <- grep("^animal_genotype", colnames(design), value = TRUE)[1]

  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)

  eda_dir <- file.path(args$fig_dir, "eda")
  dir.create(eda_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- paste("LimmaVoomCC", safe_filename(sex_i), safe_filename(tp_i), safe_filename(cluster_i), sep = "_")
  pdf(file.path(eda_dir, paste0(prefix, "_voom_mean_variance.pdf")), width = 7, height = 6)
  v <- limma::voom(dge, design = design, plot = TRUE)
  dev.off()

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
  de <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none", adjust.method = "BH")
  de$gene <- rownames(de)
  de$avg_log2FC <- de$logFC
  de$p_val <- de$P.Value
  de$p_val_adj <- de$adj.P.Val
  de$pct_mut <- expressed_mut[de$gene]
  de$pct_wt <- expressed_wt[de$gene]
  de$max_pct <- pmax(de$pct_mut, de$pct_wt)
  de$animal_sex <- sex_i
  de$animal_timepoint <- tp_i
  de[[cluster_col]] <- cluster_i
  de$n_wt_cells <- unname(tab["wt"])
  de$n_mut_cells <- unname(tab["mut"])
  de$method <- "LimmaVoomCC"
  de$cluster_col <- cluster_col
  de$latent_vars_used <- if (is.null(dupcor)) NA_character_ else "animal_genotype_animal_sex_animal_timepoint_duplicateCorrelation"
  de
}
# ---------------------------------------------------------------------------- #
# FindMarkers DEG loop (sex x timepoint x cluster_col)                          #
# ---------------------------------------------------------------------------- #
run_findmarkers_deg <- function(seu,
                                cluster_col = "neuron_class",
                                test_method = "wilcox",
                                min_cells_per_group = 10,
                                min_pct = 0.25) {
  strata <- seu@meta.data %>%
    distinct(animal_sex, animal_timepoint, .data[[cluster_col]]) %>%
    rename(cluster_id = all_of(cluster_col)) %>%
    arrange(animal_sex, animal_timepoint, cluster_id)
  results_list <- list()
  for (i in seq_len(nrow(strata))) {
    sex_i <- strata$animal_sex[i]
    tp_i <- strata$animal_timepoint[i]
    cluster_i <- strata$cluster_id[i]
    # ---- subset to stratum --------------------------------------------------
    meta_df <- seu@meta.data
    cells_i <- rownames(meta_df[
      meta_df$animal_sex == sex_i &
        meta_df$animal_timepoint == tp_i &
        meta_df[[cluster_col]] == cluster_i,
    ])
    sub <- subset(seu, cells = cells_i)
    Idents(sub) <- sub$animal_genotype
    # ---- cell-count guard ---------------------------------------------------
    tab <- table(Idents(sub))
    if (length(tab) < 2 || any(tab < min_cells_per_group)) {
      message(sprintf(
        "[%d/%d] Skipping sex=%s timepoint=%s %s=%s  (n_wt=%s, n_mut=%s)",
        i, nrow(strata), sex_i, tp_i, cluster_col, cluster_i,
        tab["wt"] %||% 0, tab["mut"] %||% 0
      ))
      next
    }
    message(sprintf(
      "[%d/%d] Running [%s]: sex=%s timepoint=%s %s=%s  (n=%d cells)",
      i, nrow(strata), test_method, sex_i, tp_i, cluster_col, cluster_i, ncol(sub)
    ))
    # ---- MAST-only: CDR covariate (as scaled nFeature) ----------------------
    # CDR_i = nFeature_i / n_genes is a linear rescaling of nFeature_i by a
    # constant (n_genes is fixed within this stratum's `data` layer). Once you
    # z-score it, that constant cancels out, so scale(nFeature) == scale(CDR)
    # numerically. Computed from the SAME `data` layer used by FindMarkers
    # (non-zero status is identical between counts and log-normalized data),
    # and computed WITHIN this stratum's subset (not globally), matching the
    # standard MAST vignette recommendation.
    latent_vars_i <- NULL
    if (test_method == "MAST") {
      nfeat_i <- Matrix::colSums(GetAssayData(sub, layer = "data") > 0)
      sub$nFeature_scaled <- as.numeric(scale(nfeat_i))
      latent_vars_i <- "nFeature_scaled"
    }
    # ---- idempotency: skip if output CSV already exists --------------------
    out_csv <- file.path(
      args$file_dir,
      sprintf("deg_%s_%s_%s.csv", sex_i, tp_i, cluster_i)
    )
    if (file.exists(out_csv)) {
      de_cached <- read.csv(out_csv, stringsAsFactors = FALSE)
      n_sig_cached <- sum(de_cached$p_val_adj < 0.05, na.rm = TRUE)
      message(sprintf(
        "  Already done (cached): sex=%s timepoint=%s %s=%s [%s]  %d genes, %d sig (p_val_adj<0.05) -- %s",
        sex_i, tp_i, cluster_col, cluster_i, test_method, nrow(de_cached), n_sig_cached, basename(out_csv)
      ))
      results_list[[paste(sex_i, tp_i, cluster_i, sep = "_")]] <- de_cached
      next
    }
    # ---- DEG test -----------------------------------------------------------
    t_start <- Sys.time()
    de <- if (test_method == "LimmaVoomCC") {
      tryCatch(
        run_limmavoomcc_deg(sub, cluster_col, cluster_i, sex_i, tp_i, tab),
        error = function(e) {
          message("  LimmaVoomCC failed: ", conditionMessage(e))
          NULL
        }
      )
    } else {
      tryCatch(
        do.call(FindMarkers, c(
          list(
            object           = sub,
            ident.1          = "mut",
            ident.2          = "wt",
            test.use         = test_method,
            slot             = "data",
            min.pct          = min_pct,
            pseudocount.use  = 1,
            logfc.threshold  = 0
          ),
          if (!is.null(latent_vars_i)) list(latent.vars = latent_vars_i) else NULL
        )),
        error = function(e) {
          message("  FindMarkers failed: ", conditionMessage(e))
          NULL
        }
      )
    }
    if (is.null(de)) next
    # ---- annotate: stratum metadata -----------------------------------------
    if (test_method != "LimmaVoomCC") {
      de$gene <- rownames(de)
      de$animal_sex <- sex_i
      de$animal_timepoint <- tp_i
      de[[cluster_col]] <- cluster_i
      de$n_wt_cells <- unname(tab["wt"])
      de$n_mut_cells <- unname(tab["mut"])
      de$method <- test_method
      de$cluster_col <- cluster_col
      de$latent_vars_used <- if (is.null(latent_vars_i)) NA_character_ else latent_vars_i
      # pct.1 = fraction in ident.1 (mut); pct.2 = fraction in ident.2 (wt)
      de$pct_mut <- de$pct.1
      de$pct_wt <- de$pct.2
      de$max_pct <- pmax(de$pct.1, de$pct.2)
    }
    # ---- write per-stratum CSV ----------------------------------------------
    write.csv(de, out_csv, row.names = FALSE)
    elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    n_sig <- sum(de$p_val_adj < 0.05, na.rm = TRUE)
    message(sprintf(
      "  Done: sex=%s timepoint=%s %s=%s [%s]  %d genes tested, %d sig (p_val_adj<0.05), %.1fs elapsed",
      sex_i, tp_i, cluster_col, cluster_i, test_method, nrow(de), n_sig, elapsed_sec
    ))
    message("  Saved: ", out_csv)
    results_list[[paste(sex_i, tp_i, cluster_i, sep = "_")]] <- de
  }
  bind_rows(results_list)
}

# ---------------------------------------------------------------------------- #
# Visualizing Results                                                            #
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
  if (!is.finite(lim) || lim == 0) lim <- max_abs
  if (!is.finite(lim) || lim == 0) lim <- 1
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

save_deg_visualizations <- function(deg, fig_dir, cluster_col, method, top_n_volcano = 10, top_n_heatmap = 5) {
  cluster_col_name <- cluster_col

  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
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

  top_volcano_genes <- deg %>%
    group_by(stratum) %>%
    slice_min(order_by = p_val_adj, n = top_n_volcano, with_ties = FALSE) %>%
    ungroup()

  top_heatmap_genes <- deg %>%
    group_by(stratum) %>%
    slice_min(order_by = p_val_adj, n = top_n_heatmap, with_ties = FALSE) %>%
    ungroup()

  sig_levels <- c("down in mut", "ns", "up in mut")
  sig_colors <- setNames(c("#2166AC", "grey80", "#B2182B"), sig_levels)

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
    top_sx <- top_volcano_genes %>%
      filter(animal_sex == sx) %>%
      mutate(
        timepoint_display = factor(timepoint_label(animal_timepoint), levels = levels(deg_sx$timepoint_display)),
        cluster_display = factor(.data[[cluster_col_name]], levels = levels(deg_sx$cluster_display))
      )

    p <- ggplot(deg_sx, aes(x = avg_log2FC, y = neglog10padj, color = sig)) +
      geom_point(size = 0.6, alpha = 0.7) +
      scale_color_manual(values = sig_colors, name = NULL, drop = FALSE) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.2, color = "grey50") +
      theme_minimal(base_size = 8) +
      theme(
        strip.text = element_text(size = 7, face = "bold"),
        strip.background = element_rect(fill = "white", color = "grey50", linewidth = 0.6),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
      ) +
      labs(
        title = paste("Cell-level DEG volcano plots -", method, "-", format_sex_label(sx), "(genotype: mut vs wt)"),
        x = "avg log2FC (mut vs wt)",
        y = "-log10(adjusted p-value)"
      )
    p <- add_volcano_label_layers(p, top_sx, sig_colors)

    n_tp <- length(levels(droplevels(deg_sx$timepoint_display)))
    n_class <- length(levels(droplevels(deg_sx$cluster_display)))
    timepoint_strip_colors <- timepoint_colors_map[levels(droplevels(deg_sx$timepoint_display))]
    if (any(is.na(timepoint_strip_colors))) timepoint_strip_colors[is.na(timepoint_strip_colors)] <- "#8C8C8C"
    cluster_strip_colors <- cluster_colors_for(cluster_col, levels(droplevels(deg_sx$cluster_display)))

    save_colored_facet_pdf(
      p + facet_grid(timepoint_display ~ cluster_display, scales = "free_x", switch = "y"),
      file.path(fig_dir, sprintf("volcano_%s_rows_timepoint_cols_%s.pdf", sx, cluster_col)),
      width = max(8, n_class * 3),
      height = max(5, n_tp * 3),
      top_colors = cluster_strip_colors,
      row_colors = timepoint_strip_colors
    )
    save_colored_facet_pdf(
      p + facet_grid(cluster_display ~ timepoint_display, scales = "free_x", switch = "y"),
      file.path(fig_dir, sprintf("volcano_%s_rows_%s_cols_timepoint.pdf", sx, cluster_col)),
      width = max(8, n_tp * 3),
      height = max(5, n_class * 3),
      top_colors = timepoint_strip_colors,
      row_colors = cluster_strip_colors
    )
  }

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    message("ComplexHeatmap/circlize not installed; skipping heatmaps.")
    return(invisible(NULL))
  }

  for (sx in unique(deg$animal_sex)) {
    deg_sx <- deg %>% filter(animal_sex == sx)
    top_sx <- top_heatmap_genes %>% filter(animal_sex == sx)
    heat_genes <- unique(top_sx$gene)
    if (length(heat_genes) == 0) next

    wide <- deg_sx %>%
      filter(gene %in% heat_genes) %>%
      select(gene, stratum, avg_log2FC) %>%
      pivot_wider(names_from = stratum, values_from = avg_log2FC, values_fill = 0)
    if (nrow(wide) == 0) next
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
    title_col <- sex_colors_map[[format_sex_label(sx)]] %||% "#333333"
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
      column_title = sex_title(sx),
      column_title_gp = grid::gpar(col = title_col, fontface = "bold", fontsize = 12)
    )
    pdf(file.path(fig_dir, sprintf("heatmap_top_genes_%s.pdf", sx)),
      width = max(9, ncol(mat) * 0.35 + 3),
      height = max(6, nrow(mat) * 0.2 + 2)
    )
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legend = TRUE
    )
    dev.off()
  }

  invisible(NULL)
}

# ---------------------------------------------------------------------------- #
# Run DEG                                                                        #
# ---------------------------------------------------------------------------- #
combined_csv <- file.path(
  args$file_dir,
  sprintf("deg_all_%s_%s.csv", tolower(args$method), args$cluster_col)
)

if (isTRUE(args$viz_only)) {
  if (!file.exists(combined_csv)) {
    stop("Combined DEG CSV not found for viz-only mode: ", combined_csv)
  }
  message("Viz-only mode: reading combined DEG table: ", combined_csv)
  deg_all <- read.csv(combined_csv, stringsAsFactors = FALSE)
  save_deg_visualizations(
    deg = deg_all,
    fig_dir = args$fig_dir,
    cluster_col = args$cluster_col,
    method = args$method
  )
  message("Done.")
  quit(save = "no", status = 0)
}

message("Loading Seurat object: ", args$seu)
seu <- readRDS(args$seu)
message("Seurat object loaded.")
message("  Counts dim  : ", paste(dim(GetAssayData(seu, layer = "counts")), collapse = " x "))
message("  Data dim    : ", paste(dim(GetAssayData(seu, layer = "data")), collapse = " x "))
message("  Default idents (first 5): ", paste(head(Idents(seu), 5), collapse = ", "))

deg_all <- run_findmarkers_deg(
  seu                 = seu,
  cluster_col         = args$cluster_col,
  test_method         = args$method,
  min_cells_per_group = min_cells_per_group,
  min_pct             = args$min_pct
)
# ---- write combined CSV -----------------------------------------------------
write.csv(deg_all, combined_csv, row.names = FALSE)
message("Combined DEG table saved: ", combined_csv)
save_deg_visualizations(
  deg = deg_all,
  fig_dir = args$fig_dir,
  cluster_col = args$cluster_col,
  method = args$method
)
# ---- summary ----------------------------------------------------------------
strata_present <- deg_all %>%
  distinct(animal_sex, animal_timepoint, .data[[args$cluster_col]]) %>%
  arrange(animal_sex, animal_timepoint, .data[[args$cluster_col]])
message("Strata present in output (", nrow(strata_present), "):")
print(strata_present)
message("Done.")
