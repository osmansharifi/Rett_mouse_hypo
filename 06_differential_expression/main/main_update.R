#!/usr/bin/env Rscript
# Benchmark: each method fits each two-group contrast separately.
# Positive effects: female - male; mut - wt (same direction as original script).
# Example: Rscript main_update.R --method all --cluster_col neuron_class
suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Matrix)
  library(dplyr)
})
parser <- ArgumentParser(description = "Pairwise cell-level DEG benchmark")
parser$add_argument("--method", default = "all",
  choices = c("all", "wilcox", "MAST", "MAST_RE", "voomLmFit", "voomLmFit_block"))
parser$add_argument("--cluster_col", required = TRUE,
  choices = c("neuron_class", "cell_type_concise"))
parser$add_argument("--seu", default = "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects/seu_final.rds")
parser$add_argument("--file_dir", default = "/quobyte/lasallegrp/Osman/shenyu/main/results")
parser$add_argument("--fig_dir", default = "/quobyte/lasallegrp/Osman/shenyu/main/figures")
parser$add_argument("--min_pct", type = "double", default = 0.25)
parser$add_argument("--min_cells", type = "integer", default = 3)
parser$add_argument("--resume", action = "store_true",
  help = "Reuse existing comparison CSVs (same inputs/settings) and rebuild the combined CSV")
args <- parser$parse_args()
stopifnot(args$min_pct > 0, args$min_pct <= 1, args$min_cells >= 1)
methods <- if (args$method == "all") {
  c("wilcox", "MAST", "MAST_RE", "voomLmFit", "voomLmFit_block")
} else args$method
safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("^_+|_+$", "", x)
}
# Isolate both the full fit and the reduced-model LRT for each gene.
# ebayes=FALSE makes these fits independent across genes. Keep CDR computed
# on ALL input genes above; do not recompute it in the single-gene subset.
fit_mast_re_gene <- function(gene, sca, coef_i,
                             fit_fun = MAST::zlm, summary_fun = MAST::summary) {
  result <- data.frame(
    primerid = gene, gene = gene, p_val = NA_real_, coef = NA_real_,
    error_message = "", stringsAsFactors = FALSE
  )
  parsing_result <- FALSE
  tryCatch({
    # MAST 1.32 summary drops dimensions for a one-gene ZlmFit.
    # A second identical row avoids that bug; ebayes=FALSE keeps fits independent.
    # Only the original gene is extracted below; the copy never enters the CSV/BH.
    sca_gene <- sca[c(gene, gene), ]
    gene_ids <- c(gene, paste0(gene, "__summary_copy"))
    rownames(sca_gene) <- gene_ids
    SummarizedExperiment::rowData(sca_gene)$primerid <- gene_ids
    fit <- fit_fun(~ group + cngeneson + (1 | sample),
      sca = sca_gene, method = "glmer", ebayes = FALSE,
      silent = FALSE, parallel = FALSE)
    mast_summary <- summary_fun(fit, doLRT = coef_i)
    parsing_result <- TRUE
    dt <- if (is.data.frame(mast_summary)) {
      mast_summary
    } else if (is.list(mast_summary) && "datatable" %in% names(mast_summary)) {
      mast_summary[["datatable"]]
    } else if (isS4(mast_summary) &&
               "datatable" %in% methods::slotNames(mast_summary)) {
      methods::slot(mast_summary, "datatable")
    } else {
      message("Skipping result extraction for ", gene, ": unsupported summary format")
      return(result)
    }
    dt <- as.data.frame(dt)
    required <- c("primerid", "contrast", "component", "Pr(>Chisq)", "coef")
    if (!all(required %in% names(dt))) {
      message("Skipping result extraction for ", gene, ": missing summary columns")
      return(result)
    }

    # Missing output remains NA; do not manufacture errors from result checks.
    rows <- which(dt$primerid == gene & dt$contrast == coef_i)
    h <- rows[which(dt$component[rows] == "H")]
    fc <- rows[which(dt$component[rows] == "logFC")]
    if (length(h) && length(dt[["Pr(>Chisq)"]])) {
      result$p_val <- dt[["Pr(>Chisq)"]][h[1]]
    }
    if (length(fc) && length(dt$coef)) {
      result$coef <- dt$coef[fc[1]]
    }
    result
  }, error = function(e) {
    if (parsing_result) {
      message("Result extraction failed for ", gene, ": ", conditionMessage(e))
    } else {
      result$error_message <- conditionMessage(e)
    }
    result
  })
}

read_resume_result <- function(path) {
  tryCatch({
    de <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("gene", "p_val", "p_val_adj", "avg_log2FC")
    if (!nrow(de) || !all(required %in% names(de)) ||
        anyNA(de$gene) || anyDuplicated(de$gene)) stop("Invalid result CSV")
    if (!"error_message" %in% names(de)) de$error_message <- rep("", nrow(de))
    de$error_message <- as.character(de$error_message)
    de$error_message[is.na(de$error_message)] <- ""
    if ("animal_timepoint" %in% names(de)) de$animal_timepoint <- as.character(de$animal_timepoint)
    de
  }, error = function(e) {
    message("Recomputing unreadable/invalid CSV: ", path, " (", conditionMessage(e), ")")
    NULL
  })
}

contrasts <- data.frame(
  contrast = c("sex_wt", "sex_mut", "genotype_male", "genotype_female"),
  reference = c("wt_male", "mut_male", "wt_male", "wt_female"),
  comparison = c("wt_female", "mut_female", "mut_male", "mut_female")
)
seu <- readRDS(args$seu)
DefaultAssay(seu) <- "RNA"
meta_all <- seu[[]]
required_cols <- c("animal_genotype", "animal_sex", "animal_timepoint", "sample", args$cluster_col)
stopifnot(all(required_cols %in% colnames(meta_all)))
meta_all$group <- paste(meta_all$animal_genotype, meta_all$animal_sex, sep = "_")
meta_all$analysis_timepoint <- as.numeric(as.character(meta_all$animal_timepoint))
if (anyNA(meta_all$analysis_timepoint)) stop("Missing or non-numeric animal_timepoint")
# Preserve original metadata; male P120 and female P150 represent the late stage.
meta_all$analysis_timepoint[meta_all$animal_sex == "female" &
  meta_all$analysis_timepoint == 150] <- 120
strata <- meta_all %>%
  distinct(analysis_timepoint, .data[[args$cluster_col]]) %>%
  rename(cluster_id = all_of(args$cluster_col)) %>%
  arrange(analysis_timepoint, cluster_id)

for (method_i in methods) {
  out_dir <- file.path(args$file_dir, method_i, args$cluster_col)
  eda_dir <- file.path(args$fig_dir, method_i, args$cluster_col, "eda")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(eda_dir, recursive = TRUE, showWarnings = FALSE)
  all_results <- list()
  for (i in seq_len(nrow(strata))) {
    tp_i <- strata$analysis_timepoint[i]
    cluster_i <- strata$cluster_id[i]
    for (j in seq_len(nrow(contrasts))) {
      contrast_i <- contrasts$contrast[j]
      ref_i <- contrasts$reference[j]
      comp_i <- contrasts$comparison[j]
      cells_i <- rownames(meta_all)[which(
        meta_all$analysis_timepoint == tp_i &
          meta_all[[args$cluster_col]] == cluster_i &
          meta_all$group %in% c(ref_i, comp_i))]
      meta <- meta_all[cells_i, , drop = FALSE]
      meta$group <- factor(meta$group, levels = c(ref_i, comp_i))
      tab <- table(meta$group)
      original_timepoints <- unique(as.character(meta$animal_timepoint))
      timepoint_label <- if (length(original_timepoints) == 1L) {
        original_timepoints
      } else {
        age_labels <- unique(meta[, c("animal_sex", "animal_timepoint")])
        age_labels <- age_labels[order(match(age_labels$animal_sex, c("male", "female"))), ]
        paste(paste0(age_labels$animal_sex, "_P", age_labels$animal_timepoint), collapse = "_")
      }
      key <- paste(contrast_i, safe_filename(timepoint_label), safe_filename(cluster_i), sep = "_")
      if (any(tab < args$min_cells)) {
        message("Skipping ", method_i, " / ", key, ": insufficient cells")
        next
      }
      out_csv <- file.path(out_dir, paste0("deg_", key, ".csv"))
      if (args$resume && file.exists(out_csv)) {
        previous <- read_resume_result(out_csv)
        if (!is.null(previous)) {
          all_results[[key]] <- previous
          message("Reused: ", out_csv)
          next
        }
      }
      sub <- subset(seu, cells = cells_i)
      meta <- meta[colnames(sub), , drop = FALSE]
      sub$group <- meta$group
      Idents(sub) <- "group"
      counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
      expr <- GetAssayData(sub, assay = "RNA", layer = "data")
      # CDR is computed over ALL input genes BEFORE any gene filtering.
      # Scale within this pairwise fit; this is equivalent to scaled nFeature.
      if (method_i %in% c("MAST", "MAST_RE")) {
        cdr <- Matrix::colMeans(expr > 0)
        if (!is.finite(sd(cdr)) || sd(cdr) == 0) {
          stop("CDR has zero variance: ", key)
        }
        meta$cngeneson <- as.numeric(scale(cdr))
        sub$cngeneson <- meta$cngeneson
      }
      message("Running ", method_i, " / ", key, " (", comp_i, " - ", ref_i, ")")
      t_start <- Sys.time()
      consensus_correlation <- NA_real_
      # Shared filter: CPM from the complete, unfiltered library in each cell.
      if (any(Matrix::colSums(counts) <= 0)) stop("Zero-count cell: ", key)
      detected <- edgeR::cpm(counts, log = FALSE) >= 1
      pct_ref <- Matrix::rowMeans(detected[, meta$group == ref_i, drop = FALSE])
      pct_comp <- Matrix::rowMeans(detected[, meta$group == comp_i, drop = FALSE])
      keep <- pmax(pct_ref, pct_comp) >= args$min_pct
      genes_keep <- rownames(counts)[keep]
      rm(detected)
      if (!length(genes_keep)) {
        message("Skipping: no genes passed filtering")
        next
      }

      # wicoxon, MAST
      if (method_i %in% c("wilcox", "MAST")) {
        de <- FindMarkers(sub, ident.1 = comp_i, ident.2 = ref_i,
          test.use = method_i, slot = "data", features = genes_keep,
          min.pct = 0, min.diff.pct = -Inf, min.cells.feature = 0,
          pseudocount.use = 1, logfc.threshold = 0,
          latent.vars = if (method_i == "MAST") "cngeneson" else NULL)
        de$gene <- rownames(de)
        # Preserve FindMarkers Bonferroni values; use BH for benchmark plots.
        de$p_val_bonferroni <- de$p_val_adj
        de$p_val_adj <- p.adjust(de$p_val, method = "BH")
      } else {
        # MAST_RE
        if (method_i == "MAST_RE") {
          meta$sample <- factor(meta$sample)
          # CDR was computed above on all genes; densify only retained genes.
          expr_keep <- expr[genes_keep, , drop = FALSE]
          fd <- data.frame(primerid = rownames(expr_keep), row.names = rownames(expr_keep))
          sca <- MAST::FromMatrix(as.matrix(expr_keep), cData = meta, fData = fd)
          rm(expr_keep)
          coef_i <- paste0("group", comp_i)
          gene_results <- vector("list", length(genes_keep))
          for (g in seq_along(genes_keep)) {
            gene_i <- genes_keep[g]
            gene_results[[g]] <- fit_mast_re_gene(
              gene_i, sca, coef_i)
            if (nzchar(gene_results[[g]]$error_message)) {
              message("Failed ", method_i, " / ", key, " / ", gene_i,
                ": ",
                gene_results[[g]]$error_message)
            }
            if (g %% 100L == 0L || g == length(genes_keep)) {
              message("MAST_RE progress: ", g, "/", length(genes_keep), " genes / ", key)
            }
          }
          de <- bind_rows(gene_results)
          # Keep failed genes as NA; retain the full filtered testing family.
          de$p_val_adj <- p.adjust(de$p_val, method = "BH", n = length(genes_keep))
          de[["Pr(>Chisq)"]] <- de$p_val
          message("MAST_RE failed genes: ", sum(nzchar(de$error_message)),
            "/", nrow(de), " / ", key)
        } else {
          # voomLmFit with or without sample block
          design <- model.matrix(~ group, data = meta)
          coef_i <- 2
          dge <- edgeR::DGEList(counts = counts[keep, , drop = FALSE])
          dge <- edgeR::calcNormFactors(dge)
          pdf(file.path(eda_dir, paste0(key, "_mean_variance.pdf")), width = 7, height = 6)
          fit <- tryCatch({
            if (method_i == "voomLmFit_block") {
              edgeR::voomLmFit(dge, design = design, block = factor(meta$sample), plot = TRUE)
            } else {
              edgeR::voomLmFit(dge, design = design, plot = TRUE)
            }
          }, finally = dev.off())
          if (!is.null(fit$correlation)) consensus_correlation <- fit$correlation
          fit <- limma::eBayes(fit)
          de <- limma::topTable(fit, coef = coef_i, number = Inf, sort.by = "none", adjust.method = "BH")
          de$gene <- rownames(de)
          de$p_val <- de$P.Value
          de$p_val_adj <- de$adj.P.Val
        }
      }
      if (!nrow(de)) next
      # Common descriptive fold change; retain model-specific coef/logFC unchanged.
      seurat_fc <- Seurat::FoldChange(
        sub, ident.1 = comp_i, ident.2 = ref_i,
        assay = "RNA", slot = "data", features = de$gene,
        pseudocount.use = 1, base = 2, fc.name = "avg_log2FC"
      )
      de$avg_log2FC <- seurat_fc[de$gene, "avg_log2FC"]
      de$min_cpm <- 1
      de$pct_cpm1_comparison <- pct_comp[de$gene]
      de$pct_cpm1_reference <- pct_ref[de$gene]
      de$n_genes_filtered <- length(genes_keep)
      de$contrast <- contrast_i
      de$reference <- ref_i
      de$comparison <- comp_i
      de$animal_timepoint <- timepoint_label
      de[[args$cluster_col]] <- cluster_i
      de$method <- method_i
      de$cluster_col <- args$cluster_col
      de$n_reference_cells <- unname(tab[ref_i])
      de$n_comparison_cells <- unname(tab[comp_i])
      de$n_reference_samples <- length(unique(meta$sample[meta$group == ref_i]))
      de$n_comparison_samples <- length(unique(meta$sample[meta$group == comp_i]))
      de$min_pct <- args$min_pct
      de$latent_vars_used <- if (method_i %in% c("MAST", "MAST_RE")) "cngeneson" else NA_character_
      de$block <- if (method_i %in% c("MAST_RE", "voomLmFit_block")) "sample" else NA_character_
      if (method_i == "voomLmFit_block") {
        de$consensus_correlation <- consensus_correlation
      }
      de$elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      first_cols <- c("gene", "p_val", "p_val_adj", "avg_log2FC")
      model_fc_cols <- intersect(c("coef", "logFC"), names(de))
      first_cols <- c(first_cols, model_fc_cols)
      de <- de[, c(first_cols, setdiff(names(de), first_cols)), drop = FALSE]
      write.csv(de, out_csv, row.names = FALSE)
      all_results[[key]] <- de
      message("Saved: ", out_csv)
    }
  }
  write.csv(bind_rows(all_results),
    file.path(out_dir, paste0("deg_all_", method_i, "_", args$cluster_col, ".csv")), row.names = FALSE)
}
