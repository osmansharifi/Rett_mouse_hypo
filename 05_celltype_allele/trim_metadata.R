#!/usr/bin/env Rscript
# Remove redundant / intermediate metadata columns. Prints kept vs dropped and
# writes a NEW file (original untouched) so nothing is lost if a column is wanted.
# Usage: Rscript 12_trim_metadata.R <in.rds> <out.rds>
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE); inf <- args[1]; outf <- args[2]
obj <- readRDS(inf)

# columns to DROP: duplicates + intermediate scaffolding
drop <- c(
  # duplicated animal metadata (keep the animal_* copies)
  "sex","genotype","time_point","animal_sample_name",
  # constants / recomputable / stale QC
  "old.ident","discard","nCount_RNA_normalized",
  # superseded per-cell C25 transfer (reconciled cols kept)
  "celltype","celltype_score","celltype_cluster",
  # reconciliation bookkeeping
  "hier_consistent","hier_confident","hier_level",
  # clustering sweep except the working resolution (0.6)
  "RNA_snn_res.0.2","RNA_snn_res.0.4","RNA_snn_res.0.8","RNA_snn_res.1"
)
# safety: never drop something that doesn't exist, and never drop core columns
core <- c("orig.ident","nCount_RNA","nFeature_RNA","sample","mouse",
          "cell_type","cell_type_c25","celltype_hier","mecp2_allele","mecp2_status")
drop <- setdiff(intersect(drop, colnames(obj[[]])), core)

cat("DROPPING", length(drop), "columns:\n"); print(drop)
for (d in drop) obj[[d]] <- NULL
cat("\nKEEPING", ncol(obj[[]]), "columns:\n"); print(colnames(obj[[]]))

saveRDS(obj, outf)
message("\nsaved ", outf)
