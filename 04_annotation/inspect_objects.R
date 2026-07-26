#!/usr/bin/env Rscript
# Pre-flight: report what is actually inside set_scaled.rds and hypoMap.rds
# so the clustering / annotation params are chosen from fact, not assumption.
# Usage: Rscript 05_inspect.R <query.rds> <reference.rds>
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly = TRUE)
qf <- args[1]; rf <- args[2]

report <- function(o, tag) {
  cat("\n=====", tag, "=====\n")
  cat("class:", class(o)[1], " assays:", paste(Assays(o), collapse=", "),
      " default:", DefaultAssay(o), "\n")
  cat("dims:", nrow(o), "features x", ncol(o), "cells\n")
  a <- DefaultAssay(o)
  cat("layers:", paste(Layers(o[[a]]), collapse=", "), "\n")
  cat("reductions:", paste(Reductions(o), collapse=", "), "\n")
  # is it normalized? (data layer present AND different from counts)
  L <- Layers(o[[a]])
  if ("data" %in% L && "counts" %in% L) {
    d <- LayerData(o, assay=a, layer="data")[1:min(50,nrow(o)), 1:min(50,ncol(o))]
    c_ <- LayerData(o, assay=a, layer="counts")[1:min(50,nrow(o)), 1:min(50,ncol(o))]
    cat("normalized:", !identical(as.numeric(d), as.numeric(c_)),
        " (max data value:", round(max(d),2), ")\n")
  } else cat("normalized: cannot tell (layers missing)\n")
  cat("scaled:", "scale.data" %in% L, "\n")
  cat("variable features:", length(VariableFeatures(o)), "\n")
  cat("\n-- metadata columns --\n"); print(colnames(o[[]]))
  cat("\n-- candidate label columns (character/factor, 2..1000 levels) --\n")
  md <- o[[]]
  for (cn in colnames(md)) {
    v <- md[[cn]]
    if (is.character(v) || is.factor(v)) {
      n <- length(unique(v))
      if (n >= 2 && n <= 1000) cat(sprintf("  %-28s %5d levels  e.g. %s\n", cn, n,
        paste(head(unique(as.character(v)),3), collapse=" | ")))
    }
  }
  invisible(NULL)
}

q <- readRDS(qf); report(q, "QUERY (set_scaled.rds)")
cat("\n-- query: cells per sample / timepoint / genotype --\n")
for (cn in c("sample","animal_timepoint","timepoint","animal_genotype","animal_sex")) {
  if (cn %in% colnames(q[[]])) { cat("\n", cn, ":\n"); print(table(q[[]][[cn]], useNA="ifany")) }
}
rm(q); gc()

r <- readRDS(rf); report(r, "REFERENCE (hypoMap.rds)")
cat("\nHypoMap annotation levels usually: C2/C7/C25/C66/C185/C286/C465 (+ _named).\n")
cat("Pick ONE for --label (C25_named or C66_named are good starting granularity).\n")
