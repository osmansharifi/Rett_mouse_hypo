#!/usr/bin/env Rscript
# Paired pseudobulk DEG: MUT-active vs WT-active cells WITHIN het females.
#
# Design rationale: both alleles occur in the same animal (X-inactivation mosaic),
# so `~ mouse + allele` blocks on animal and removes animal/batch/age/dissection
# effects. Cells are aggregated to pseudobulk per (mouse x celltype x allele) first,
# because testing individual cells treats cells from one mouse as independent
# replicates and badly inflates significance.
#
# Usage: Rscript 10_pseudobulk_deg.R <obj.rds> <outdir> [min_cells] [min_pairs] [label_col]
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
inf <- args[1]; outdir <- args[2]
MIN_CELLS <- ifelse(length(args) >= 3, as.integer(args[3]), 10)   # cells per unit
MIN_PAIRS <- ifelse(length(args) >= 4, as.integer(args[4]), 4)    # animals w/ both alleles
LABEL     <- ifelse(length(args) >= 5, args[5], "celltype_hier")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE)

if (!requireNamespace("DESeq2", quietly = TRUE))
  stop("DESeq2 not installed. Run:\n  Rscript -e 'BiocManager::install(\"DESeq2\")'")
suppressPackageStartupMessages(library(DESeq2))

obj <- readRDS(inf)
md  <- obj[[]]
stopifnot(all(c("animal_sex","animal_genotype","mecp2_allele","mouse") %in% names(md)))

# ---- the mosaic subset: het females, cells with a confident allele call ----------
keep <- which(md$animal_sex == "female" & md$animal_genotype == "mut" &
              md$mecp2_allele %in% c("WT","MUT"))
message(sprintf("mosaic cells: %d across %d animals",
                length(keep), length(unique(md$mouse[keep]))))
md_k <- md[keep, ]
md_k$celltype <- as.character(md_k[[LABEL]])
md_k$unit <- paste(md_k$mouse, md_k$celltype, md_k$mecp2_allele, sep = "|")

# ---- feasibility: which cell types have enough paired units? ---------------------
cnt <- as.data.table(md_k)[, .(n = .N), by = .(mouse, celltype, allele = mecp2_allele)]
fwrite(cnt[order(celltype, mouse, allele)], file.path(outdir, "unit_cell_counts.csv"))

ok_units <- cnt[n >= MIN_CELLS]
pairs <- ok_units[, .(alleles = uniqueN(allele)), by = .(celltype, mouse)][alleles == 2]
viable <- pairs[, .(n_pairs = .N), by = celltype][n_pairs >= MIN_PAIRS][order(-n_pairs)]
cat("\n=== feasibility (units need >=", MIN_CELLS, "cells; celltypes need >=",
    MIN_PAIRS, "paired animals) ===\n"); print(viable)
if (nrow(viable) == 0)
  stop("no cell type has enough paired units. Lower --min_cells/--min_pairs, or ",
       "rerun with label_col = 'class_broad' or 'neuron_class'.")

# ---- aggregate raw counts to pseudobulk ------------------------------------------
cnts <- LayerData(obj, assay = "RNA", layer = "counts")[, keep, drop = FALSE]
u    <- factor(md_k$unit)
message("aggregating ", ncol(cnts), " cells into ", nlevels(u), " pseudobulk units")
# sparse cell -> unit indicator, so aggregation is a single matrix product
ind  <- sparse.model.matrix(~ 0 + u)
colnames(ind) <- levels(u)
pb   <- as.matrix(cnts %*% ind)                       # genes x units

meta <- data.table(unit = colnames(pb))
meta[, c("mouse","celltype","allele") := tstrsplit(unit, "|", fixed = TRUE)]
meta <- merge(meta, cnt, by = c("mouse","celltype","allele"), all.x = TRUE)
setkey(meta, unit); meta <- meta[colnames(pb)]

# ---- per-celltype paired DESeq2 --------------------------------------------------
all_res <- list()
for (ct in viable$celltype) {
  keep_pairs <- pairs[celltype == ct, mouse]
  sel <- which(meta$celltype == ct & meta$mouse %in% keep_pairs & meta$n >= MIN_CELLS)
  m   <- meta[sel]; M <- pb[, sel, drop = FALSE]
  # both alleles must survive for every retained mouse
  good <- m[, .(k = uniqueN(allele)), by = mouse][k == 2, mouse]
  sel2 <- which(m$mouse %in% good); m <- m[sel2]; M <- M[, sel2, drop = FALSE]
  if (length(unique(m$mouse)) < MIN_PAIRS) next

  cd <- data.frame(mouse  = factor(m$mouse),
                   allele = factor(m$allele, levels = c("WT","MUT")),
                   row.names = m$unit)
  M <- M[rowSums(M) > 0, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(round(M), cd, design = ~ mouse + allele)
  dds <- dds[rowSums(counts(dds) >= 5) >= ncol(dds)/2, ]   # light expression filter
  dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) NULL)
  if (is.null(dds)) { message("DESeq failed for ", ct); next }

  r <- as.data.frame(results(dds, name = "allele_MUT_vs_WT"))
  r$gene <- rownames(r); r$celltype <- ct
  r$n_pairs <- length(unique(m$mouse)); r$n_units <- ncol(M)
  r <- r[order(r$padj), ]
  fwrite(r, file.path(outdir, paste0("DEG_", gsub("[^A-Za-z0-9]+","_",ct), ".csv")))
  all_res[[ct]] <- r

  sig <- sum(r$padj < 0.05, na.rm = TRUE)
  cat(sprintf("%-32s pairs=%d units=%d genes=%d  sig(padj<0.05)=%d\n",
              ct, length(unique(m$mouse)), ncol(M), nrow(r), sig))

  v <- r[!is.na(r$padj), ]
  pdf(file.path(outdir, "figures", paste0("volcano_", gsub("[^A-Za-z0-9]+","_",ct), ".pdf")),
      width = 7, height = 6)
  print(ggplot(v, aes(log2FoldChange, -log10(padj))) +
        geom_point(aes(colour = padj < 0.05), size = .8, alpha = .6) +
        scale_colour_manual(values = c("grey70","firebrick")) +
        geom_hline(yintercept = -log10(0.05), linetype = 2) + theme_classic() +
        labs(title = paste0(ct, "  (MUT-active vs WT-active, ", sig, " sig)"),
             subtitle = "positive log2FC = higher in MUT-active cells"))
  dev.off()
}

if (length(all_res)) {
  comb <- rbindlist(all_res, fill = TRUE)
  fwrite(comb, file.path(outdir, "DEG_all_celltypes.csv"))
  cat("\n=== significant genes per cell type (padj<0.05) ===\n")
  print(comb[!is.na(padj) & padj < 0.05, .N, by = celltype][order(-N)])
  message("\nresults written to ", outdir)
} else message("no cell type produced results")

cat("\nNOTE: allele calls require a Mecp2 read over the start codon, so tested cells\n",
    "are enriched for Mecp2 expressors relative to the tissue. The MUT-vs-WT contrast\n",
    "is internally symmetric, but absolute expression statements inherit that selection.\n")
