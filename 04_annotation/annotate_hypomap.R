#!/usr/bin/env Rscript
# Reference-based cell type labeling against HypoMap (Steuernagel et al. 2022).
# Usage: Rscript 07_annotate.R <query.rds> <hypoMap.rds> <out.rds> <figdir> [label_col] [dims]
#   label_col : HypoMap metadata column to transfer (default C25_named; run 05_inspect.R to list)
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
qf <- args[1]; rf <- args[2]; outf <- args[3]; figdir <- args[4]
label <- ifelse(length(args) >= 5, args[5], "C25_named")
dims  <- 1:ifelse(length(args) >= 6, as.integer(args[6]), 30)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fig <- function(n) file.path(figdir, n)

obj <- readRDS(qf); DefaultAssay(obj) <- "RNA"
ref <- readRDS(rf)
if (inherits(try(Version(ref), silent=TRUE), "try-error") ||
    package_version(Version(ref)) < package_version("5.0.0")) {
  message("updating reference object to current Seurat"); ref <- UpdateSeuratObject(ref)
}
DefaultAssay(ref) <- if ("RNA" %in% Assays(ref)) "RNA" else DefaultAssay(ref)
message(sprintf("query %d cells | reference %d cells", ncol(obj), ncol(ref)))

if (!label %in% colnames(ref[[]])) {
  cat("available reference metadata columns:\n"); print(colnames(ref[[]]))
  stop("label column '", label, "' not in reference. Pick one from the list above.")
}
ref[[label]][,1] <- as.character(ref[[label]][,1])
message("transferring '", label, "' (", length(unique(ref[[label]][,1])), " classes)")

# reference must be normalized + have a PCA to project onto
if (!"data" %in% Layers(ref[["RNA"]])) ref <- NormalizeData(ref, verbose = FALSE)
if (!"pca" %in% Reductions(ref)) {
  message("reference has no PCA; computing (this is the slow step)")
  ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
  ref <- ScaleData(ref, verbose = FALSE)
  ref <- RunPCA(ref, npcs = max(dims), verbose = FALSE)
}

anchors <- FindTransferAnchors(reference = ref, query = obj,
             normalization.method = "LogNormalize", reference.reduction = "pca",
             dims = dims, verbose = TRUE)
pred <- TransferData(anchorset = anchors, refdata = ref[[label]][,1],
                     dims = dims, verbose = TRUE)
obj$celltype        <- pred$predicted.id
obj$celltype_score  <- pred$prediction.score.max
message("labeled. median prediction score: ", round(median(obj$celltype_score), 3))

# ------------------------------- figures ------------------------------------------
pdf(fig("09_UMAP_celltype.pdf"), width = 14, height = 9)
print(DimPlot(obj, reduction="umap", group.by="celltype", label=TRUE, repel=TRUE,
              raster=TRUE, label.size=3) + ggtitle(paste("HypoMap:", label)) +
      theme(legend.text=element_text(size=6)))
dev.off()

pdf(fig("10_prediction_score.pdf"), width = 13, height = 5)
print((FeaturePlot(obj, "celltype_score", reduction="umap", raster=TRUE) +
         ggtitle("max prediction score")) |
      (ggplot(obj[[]], aes(celltype_score)) + geom_histogram(bins=50) +
         geom_vline(xintercept=0.5, linetype=2, colour="red") + theme_classic() +
         labs(title="score distribution (low = unreliable label)")))
dev.off()

# cluster x celltype concordance: a clean diagonal means clusters map to real types
ct <- as.data.frame(table(cluster = obj$seurat_clusters, celltype = obj$celltype)) %>%
  group_by(cluster) %>% mutate(frac = Freq/sum(Freq))
pdf(fig("11_cluster_vs_celltype.pdf"), width = 13, height = 9)
print(ggplot(ct, aes(cluster, celltype, fill = frac)) + geom_tile() +
      scale_fill_viridis_c() + theme_classic() +
      theme(axis.text.x = element_text(angle=90, size=6), axis.text.y = element_text(size=5)) +
      labs(title="cluster vs transferred label (fraction of cluster)"))
dev.off()

# consensus label per cluster: usually cleaner than per-cell for downstream work
lab <- ct %>% group_by(cluster) %>% slice_max(frac, n=1, with_ties=FALSE)
obj$celltype_cluster <- lab$celltype[match(obj$seurat_clusters, lab$cluster)]
pdf(fig("12_UMAP_celltype_consensus.pdf"), width = 14, height = 9)
print(DimPlot(obj, reduction="umap", group.by="celltype_cluster", label=TRUE, repel=TRUE,
              raster=TRUE, label.size=3) + ggtitle("consensus label per cluster") +
      theme(legend.text=element_text(size=6)))
dev.off()

# INDEPENDENT VALIDATION: canonical markers must agree with the transferred labels
mk <- c("Snap25","Syt1","Rbfox3","Slc17a6","Slc32a1","Gad1","Gad2",  # neurons / E / I
        "Aqp4","Gja1","Agt","Plp1","Mbp","Mog","Pdgfra","Cspg4",     # astro / oligo / OPC
        "Cx3cr1","P2ry12","Ctss","Cldn5","Flt1","Ttr","Folr1",       # microglia / endo / ependymal
        "Agrp","Pomc","Hcrt","Pmch","Avp","Oxt","Gal","Th","Mecp2")  # hypothalamic neuropeptides
mk <- mk[mk %in% rownames(obj)]
pdf(fig("13_marker_dotplot.pdf"), width = 15, height = 10)
print(DotPlot(obj, features = mk, group.by = "celltype_cluster") + RotatedAxis() +
      theme(axis.text.y = element_text(size=6), axis.text.x = element_text(size=7)) +
      labs(title="canonical markers vs transferred labels (independent check)"))
dev.off()

pdf(fig("14_marker_featureplots.pdf"), width = 15, height = 12)
print(FeaturePlot(obj, features = head(mk, 12), reduction="umap", raster=TRUE, ncol=4))
dev.off()

# composition across the experimental design
comp_vars <- intersect(c("animal_timepoint","timepoint","animal_genotype","animal_sex"),
                       colnames(obj[[]]))
if (length(comp_vars)) {
  pdf(fig("15_composition_by_design.pdf"), width = 14, height = 5*length(comp_vars))
  print(wrap_plots(lapply(comp_vars, function(v) {
    d <- as.data.frame(table(celltype = obj$celltype_cluster, grp = obj[[]][[v]])) %>%
      group_by(grp) %>% mutate(frac = Freq/sum(Freq))
    ggplot(d, aes(grp, frac, fill = celltype)) + geom_col() + theme_classic() +
      theme(legend.text=element_text(size=5)) + labs(x=v, y="fraction", title=v)
  }), ncol = 1))
  dev.off()
}

saveRDS(obj, outf)
write.csv(as.data.frame(table(obj$celltype_cluster)), file.path(figdir, "celltype_counts.csv"))
message("saved ", outf, "; figures in ", figdir)
cat("\ncells per consensus cell type:\n"); print(sort(table(obj$celltype_cluster), decreasing=TRUE))
cat("\nlow-confidence cells (score < 0.5): ",
    sum(obj$celltype_score < 0.5), " of ", ncol(obj), "\n")
