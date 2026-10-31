#!/usr/bin/env Rscript
# Refresh all summary/evaluation CSVs from the existing neuron_class DEG results.
# This does not refit DEG models. Run from any working directory with:
#   Rscript /quobyte/lasallegrp/Osman/shenyu/main/evaluation.R

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot determine the location of evaluation.R")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
main_dir <- dirname(script_path)
results_dir <- file.path(main_dir, "results")
summary_dir <- file.path(results_dir, "robust_deg")
main_update <- file.path(main_dir, "main_update.R")
robust_script <- file.path(main_dir, "create_robust_deg.py")
num_deg_script <- file.path(main_dir, "create_num_deg_contrast.py")
methods <- c("wilcox", "MAST", "MAST_original", "MAST_RE",
  "voomLmFit", "voomLmFit_block")
required_inputs <- file.path(results_dir, methods, "neuron_class",
  paste0("deg_all_", methods, "_neuron_class.csv"))

run_checked <- function(command, args, label) {
  status <- system2(command, args = vapply(args, shQuote, character(1)))
  if (!identical(status, 0L)) {
    stop(label, " failed with exit status ", status)
  }
}

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
rscript <- Sys.which("Rscript")
python <- Sys.which("python3")
if (!nzchar(rscript)) stop("Rscript is not on PATH; activate/load the R 4.4.2 environment first")
if (!nzchar(python)) stop("python3 is not on PATH")
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing expected current-name input CSV(s):\n",
    paste(missing_inputs, collapse = "\n"))
}

# Refresh MAST coefficient visualization columns and logFC sign/NA summary tables.
run_checked(rscript, c(main_update, "--method", "all", "--cluster_col",
  "neuron_class", "--file_dir", results_dir, "--summary_only",
  "--summary_dir", summary_dir, "--skip_detailed_logfc"),
  "main_update.R summary step")

# Produce per-method significant DEG counts and cross-method robust DEG tables.
run_checked(python, num_deg_script, "num_deg_contrast step")
run_checked(python, robust_script, "robust DEG step")

message("Evaluation CSV refresh completed in: ", summary_dir)
message("Outputs: main_logFC_sign.csv, num_deg_contrast.csv, robust_deg.csv, ",
  "robust_deg_summary.csv, robust_deg_without_MAST_RE.csv, ",
  "robust_deg_summary_without_MAST_RE.csv")
