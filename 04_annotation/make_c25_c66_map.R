#!/usr/bin/env Rscript
# Extract the HypoMap C25<->C66 nesting table used by hierarchical_labels.R.
# Usage: Rscript make_c25_c66_map.R <hypoMap.rds> <out.csv>
suppressPackageStartupMessages(library(Seurat))
a <- commandArgs(trailingOnly = TRUE)
r <- readRDS(a[1])
m <- unique(r[[]][, c("C25_named","C66_named")])
m <- m[order(m$C25_named, m$C66_named), ]
write.csv(m, a[2], row.names = FALSE)
cat("wrote", nrow(m), "C66->C25 pairs to", a[2], "\n")
