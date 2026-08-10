#!/usr/bin/env Rscript
# ==============================================================================
# run_deg_cell_level.R
#
# Cell-level DEG analysis using Seurat::FindMarkers.
# Supports two test methods:  wilcox | MAST
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
# Output CSVs (one combined + one per stratum) are written to --file_dir.
# Visualisation is handled separately by run_viz.R.
# ==============================================================================
suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(purrr)
})
# ---------------------------------------------------------------------------- #
# CLI arguments                                                                  #
# ---------------------------------------------------------------------------- #
parser <- ArgumentParser(
  description = "Cell-level DEG via FindMarkers (wilcox | MAST)"
)
parser$add_argument("--method",
  type = "character", required = TRUE,
  choices = c("wilcox", "MAST"),
  help = "Test method passed to FindMarkers test.use"
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
base_dir <- if (is.null(args$file_dir)) BASE_FILE_DIR else args$file_dir
args$file_dir <- file.path(
  base_dir,
  args$method, # top-level: method name  (e.g. wilcox, MAST)
  args$cluster_col # subfolder: clustering granularity
)
dir.create(args$file_dir, showWarnings = FALSE, recursive = TRUE)
message("========================================================")
message("  Method               : ", args$method)
message("  Cluster col          : ", args$cluster_col)
message("  min.pct              : ", args$min_pct)
message("  min_cells            : ", args$min_cells)
message("  Output dir           : ", args$file_dir)
message("========================================================")
# ---------------------------------------------------------------------------- #
# Helper                                                                         #
# ---------------------------------------------------------------------------- #
`%||%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a
min_cells_per_group <- args$min_cells # minimum cells per genotype group within a stratum
# ---------------------------------------------------------------------------- #
# Load Seurat object                                                              #
# ---------------------------------------------------------------------------- #
message("Loading Seurat object: ", args$seu)
seu <- readRDS(args$seu)
message("Seurat object loaded.")
message("  Counts dim  : ", paste(dim(GetAssayData(seu, layer = "counts")), collapse = " x "))
message("  Data dim    : ", paste(dim(GetAssayData(seu, layer = "data")), collapse = " x "))
message("  Default idents (first 5): ", paste(head(Idents(seu), 5), collapse = ", "))
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
      "[%d/%d] Running FindMarkers [%s]: sex=%s timepoint=%s %s=%s  (n=%d cells)",
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
    # ---- FindMarkers --------------------------------------------------------
    t_start <- Sys.time()
    de <- tryCatch(
      do.call(FindMarkers, c(
        list(
          object           = sub,
          ident.1          = "mut",
          ident.2          = "wt",
          test.use         = test_method,
          slot             = "data", # use normalised layer (never scaled.data for MAST)
          min.pct          = min_pct,
          pseudocount.use  = 1,
          logfc.threshold  = 0 # return all genes; filter downstream
          # p_val_adj: Bonferroni correction using all genes in the dataset
        ),
        if (!is.null(latent_vars_i)) list(latent.vars = latent_vars_i) else NULL
      )),
      error = function(e) {
        message("  FindMarkers failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(de)) next
    # ---- annotate: stratum metadata -----------------------------------------
    de$gene <- rownames(de)
    de$animal_sex <- sex_i
    de$animal_timepoint <- tp_i
    de[[cluster_col]] <- cluster_i
    de$n_wt_cells <- unname(tab["wt"])
    de$n_mut_cells <- unname(tab["mut"])
    de$method <- test_method
    de$cluster_col <- cluster_col
    de$latent_vars_used <- if (is.null(latent_vars_i)) NA_character_ else latent_vars_i
    # ---- annotate: expression prevalence (min.pct) --------------------------
    # pct.1 = fraction of mut cells expressing the gene
    # pct.2 = fraction of wt  cells expressing the gene
    # These are already returned by FindMarkers; rename for clarity.
    de$pct_mut <- de$pct.1 # proportion in ident.1 (mut)
    de$pct_wt <- de$pct.2 # proportion in ident.2 (wt)
    de$max_pct <- pmax(de$pct.1, de$pct.2) # max across both groups
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
# Run DEG                                                                        #
# ---------------------------------------------------------------------------- #
deg_all <- run_findmarkers_deg(
  seu                 = seu,
  cluster_col         = args$cluster_col,
  test_method         = args$method,
  min_cells_per_group = min_cells_per_group,
  min_pct             = args$min_pct
)
# ---- write combined CSV -----------------------------------------------------
combined_csv <- file.path(
  args$file_dir,
  sprintf("deg_all_%s_%s.csv", tolower(args$method), args$cluster_col)
)
write.csv(deg_all, combined_csv, row.names = FALSE)
message("Combined DEG table saved: ", combined_csv)
# ---- summary ----------------------------------------------------------------
strata_present <- deg_all %>%
  distinct(animal_sex, animal_timepoint, .data[[args$cluster_col]]) %>%
  arrange(animal_sex, animal_timepoint, .data[[args$cluster_col]])
message("Strata present in output (", nrow(strata_present), "):")
print(strata_present)
message("Done.")
