#!/usr/bin/env Rscript

# Summarize the exact case/control input sizes recorded in completed DEG results.
# One row is emitted per analysis level, method, resolution, sex, timepoint,
# and cell group. Case = mut; control = wt.

result_root <- "/quobyte/lasallegrp/Osman/shenyu/04_results/deg"
output_file <- file.path(result_root, "samplesize_exp1and2.csv")

specs <- data.frame(
  analysis_level = c(rep("cell_level", 4), rep("sample_level", 6)),
  method = c(
    "wilcox", "wilcox", "LimmaVoomCC", "LimmaVoomCC",
    "deseq2", "deseq2", "edger_zinger", "edger_zinger",
    "limma_voomcc", "limma_voomcc"
  ),
  resolution = rep(c("cell_type_concise", "neuron_class"), 5),
  relative_file = c(
    "cell_level_view/wilcox/cell_type_concise/deg_all_wilcox_cell_type_concise.csv",
    "cell_level_view/wilcox/neuron_class/deg_all_wilcox_neuron_class.csv",
    "cell_level_view/LimmaVoomCC/cell_type_concise/deg_all_limmavoomcc_cell_type_concise.csv",
    "cell_level_view/LimmaVoomCC/neuron_class/deg_all_limmavoomcc_neuron_class.csv",
    "sample_level_view/deseq2/cell_type_concise/deg_all_deseq2_cell_type_concise.csv",
    "sample_level_view/deseq2/neuron_class/deg_all_deseq2_neuron_class.csv",
    "sample_level_view/edger_zinger/cell_type_concise/deg_all_edger_zinger_cell_type_concise.csv",
    "sample_level_view/edger_zinger/neuron_class/deg_all_edger_zinger_neuron_class.csv",
    "sample_level_view/limma_voomcc/cell_type_concise/deg_all_limma_voomcc_cell_type_concise.csv",
    "sample_level_view/limma_voomcc/neuron_class/deg_all_limma_voomcc_neuron_class.csv"
  ),
  stringsAsFactors = FALSE
)

extract_one <- function(i) {
  s <- specs[i, ]
  path <- file.path(result_root, s$relative_file)
  if (!file.exists(path)) stop("Missing DEG result: ", path)

  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  group_col <- s$resolution
  if (!group_col %in% names(x)) stop("Missing group column ", group_col, " in ", path)

  if (s$analysis_level == "cell_level") {
    wt_col <- "n_wt_cells"
    mut_col <- "n_mut_cells"
    unit <- "cells"
  } else {
    wt_col <- "n_wt_pseudobulk"
    mut_col <- "n_mut_pseudobulk"
    unit <- "pseudobulk_samples"
  }
  required <- c("animal_sex", "animal_timepoint", group_col, wt_col, mut_col)
  if (!all(required %in% names(x))) {
    stop("Missing required columns in ", path, ": ",
         paste(setdiff(required, names(x)), collapse = ", "))
  }

  z <- unique(x[required])
  names(z) <- c("sex", "timepoint", "cell_group", "control_n", "case_n")
  z$analysis_level <- s$analysis_level
  z$method <- s$method
  z$resolution <- s$resolution
  z$case <- "mut"
  z$control <- "wt"
  z$unit <- unit
  z$experiment <- ifelse(z$sex == "male", "exp1", "exp2")
  z$comparison <- paste(z$experiment, z$sex, z$timepoint, z$cell_group,
                        "mut_vs_wt", sep = "__")
  z[, c("experiment", "analysis_level", "method", "resolution", "sex",
        "timepoint", "cell_group", "comparison", "case", "case_n",
        "control", "control_n", "unit")]
}

raw <- do.call(rbind, lapply(seq_len(nrow(specs)), extract_one))

# Methods at the same analysis level use the same input for a given comparison.
# Collapse them to one comparison row, retaining all methods as provenance.
key <- c("experiment", "analysis_level", "resolution", "sex", "timepoint",
         "cell_group", "comparison", "case", "control", "unit")
size_check <- aggregate(cbind(case_n, control_n) ~ experiment + analysis_level +
                          resolution + sex + timepoint + cell_group + comparison +
                          case + control + unit,
                        data = raw, FUN = function(v) length(unique(v)))
if (any(size_check$case_n != 1 | size_check$control_n != 1)) {
  stop("Conflicting sample sizes found between DEG methods for the same comparison")
}

methods <- aggregate(method ~ experiment + analysis_level + resolution + sex +
                       timepoint + cell_group + comparison + case + control + unit,
                     data = raw,
                     FUN = function(v) paste(sort(unique(v)), collapse = ";"))
sizes <- unique(raw[c(key, "case_n", "control_n")])
out <- merge(methods, sizes, by = key, sort = FALSE)
out <- out[, c("experiment", "analysis_level", "method", "resolution", "sex",
               "timepoint", "cell_group", "comparison", "case", "case_n",
               "control", "control_n", "unit")]
out$timepoint_num <- suppressWarnings(as.numeric(out$timepoint))
out <- out[order(out$experiment, out$analysis_level, out$method,
                 out$resolution, out$timepoint_num, out$cell_group), ]
out$timepoint_num <- NULL
rownames(out) <- NULL

write.csv(out, output_file, row.names = FALSE, quote = TRUE)
message("Wrote ", nrow(out), " rows to ", output_file)
print(with(out, table(analysis_level, resolution)))
