#!/usr/bin/env Rscript
# PCA -> (optional Harmony) -> neighbors -> clustering -> UMAP, with a PDF at each step.
# Usage: Rscript 06_cluster.R <in.rds> <out.rds> <figdir> [batch_var] [npcs] [integrate]
#   batch_var : metadata column to correct over (default "sample"; "none" to skip)
#   npcs      : PCs to compute (default 50)
#   integrate : "harmony" or "none" (default "harmony")
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})

args   <- commandArgs(trailingOnly = TRUE)
inf    <- args[1]; outf <- args[2]; figdir <- args[3]
batch  <- ifelse(length(args) >= 4, args[4], "sample")
npcs   <- ifelse(length(args) >= 5, as.integer(args[5]), 50)
integ  <- ifelse(length(args) >= 6, args[6], "auto")   # auto|harmony|none
RES    <- c(0.2, 0.4, 0.6, 0.8, 1.0)     # clustering resolutions to sweep
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fig <- function(n) file.path(figdir, n)

obj <- readRDS(inf)
DefaultAssay(obj) <- "RNA"
message(sprintf("loaded %d features x %d cells", nrow(obj), ncol(obj)))

# ---- guard: label transfer downstream needs LOG-NORMALIZED data, not just scaled ----
L <- Layers(obj[["RNA"]])
if (!"data" %in% L) stop("no 'data' layer: run NormalizeData() before this script")
if ("counts" %in% L) {
  s <- 1:min(50, nrow(obj)); k <- 1:min(50, ncol(obj))
  if (identical(as.numeric(LayerData(obj, layer="data")[s,k]),
                as.numeric(LayerData(obj, layer="counts")[s,k])))
    stop("'data' layer == 'counts': data is NOT normalized. Run NormalizeData() first.")
}
if (!"scale.data" %in% L) { message("no scale.data; scaling now"); obj <- ScaleData(obj) }
if (length(VariableFeatures(obj)) == 0) {
  message("no variable features; running FindVariableFeatures")
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
}

# ---------------- QC snapshot (shows whether debris was ever filtered) -------------
qc <- intersect(c("nFeature_RNA","nCount_RNA","percent.mt"), colnames(obj[[]]))
if (length(qc)) {
  pdf(fig("01_QC_violin.pdf"), width = 14, height = 5)
  print(VlnPlot(obj, features = qc, group.by = "sample", pt.size = 0, ncol = length(qc)) &
        theme(axis.text.x = element_text(size = 6, angle = 90)))
  dev.off()
  cat("QC summary:\n"); print(summary(obj[[]][, qc, drop = FALSE]))
}

# ---------------------------------- PCA -------------------------------------------
obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
pdf(fig("02_PCA_elbow.pdf"), width = 7, height = 5); print(ElbowPlot(obj, ndims = npcs)); dev.off()
pdf(fig("03_PCA_scatter.pdf"), width = 12, height = 5)
print((DimPlot(obj, reduction="pca", group.by="sample") + NoLegend() + ggtitle("PCA by sample")) |
      (DimPlot(obj, reduction="pca", group.by=if("animal_timepoint" %in% colnames(obj[[]]))
               "animal_timepoint" else "sample") + ggtitle("PCA by timepoint")))
dev.off()
pdf(fig("04_PCA_loadings.pdf"), width = 10, height = 12)
print(VizDimLoadings(obj, dims = 1:6, reduction = "pca", ncol = 3)); dev.off()

# variance-explained: a defensible npc choice instead of eyeballing the elbow
sdev <- Stdev(obj, reduction = "pca"); pv <- sdev^2/sum(sdev^2)
n_use <- max(10, min(npcs, which(cumsum(pv) >= 0.90)[1]))
if (is.na(n_use)) n_use <- 30
message("using ", n_use, " PCs (>=90% variance or elbow floor)")

# ------------------------- integration over batch ---------------------------------
# Detect integration that already exists so we do not correct twice (over-correction
# can flatten real biology). Any of these means the object was already integrated.
known_int <- c("harmony","integrated.cca","integrated.rpca","integrated.mnn",
               "integrated.scvi","scvi","mnn","cca","rpca")
have_int  <- intersect(tolower(Reductions(obj)), known_int)
has_intassay <- "integrated" %in% Assays(obj)
cat("\n-- integration status --\n")
cat("reductions:", paste(Reductions(obj), collapse=", "), "\n")
cat("assays    :", paste(Assays(obj), collapse=", "), "\n")
cat("detected integrated reduction:", ifelse(length(have_int), have_int[1], "none"), "\n")
cat("integrated assay present     :", has_intassay, "\n")

red <- "pca"
if (length(have_int) > 0) {
  red <- Reductions(obj)[tolower(Reductions(obj)) %in% have_int][1]
  message("object is ALREADY integrated; using existing reduction '", red,
          "' and skipping Harmony")
} else if (has_intassay) {
  message("v4-style 'integrated' assay found; using its PCA (no new integration)")
} else if (integ %in% c("auto","harmony") && batch != "none") {
  if (!requireNamespace("harmony", quietly = TRUE))
    stop("harmony not installed. Run:\n  Rscript -e 'install.packages(\"harmony\", repos=\"https://cloud.r-project.org\")'")
  suppressPackageStartupMessages(library(harmony))
  obj <- RunHarmony(obj, group.by.vars = batch, reduction.use = "pca",
                    reduction.save = "harmony", dims.use = 1:n_use, plot_convergence = FALSE)
  red <- "harmony"
  message("no prior integration found; ran Harmony on '", batch, "'")
} else message("no integration (reduction = pca)")

# --------------------- neighbors, clustering sweep, UMAP --------------------------
obj <- FindNeighbors(obj, reduction = red, dims = 1:n_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = RES, verbose = FALSE)
obj <- RunUMAP(obj, reduction = red, dims = 1:n_use, verbose = FALSE)

rescols <- grep("^RNA_snn_res\\.", colnames(obj[[]]), value = TRUE)
pdf(fig("05_UMAP_resolution_sweep.pdf"), width = 16, height = 10)
print(wrap_plots(lapply(rescols, function(r)
  DimPlot(obj, reduction="umap", group.by=r, label=TRUE, raster=TRUE) +
    ggtitle(r) + NoLegend()), ncol = 3))
dev.off()

# pick a working resolution (0.6 unless absent) and set identities
pick <- if ("RNA_snn_res.0.6" %in% rescols) "RNA_snn_res.0.6" else rescols[1]
Idents(obj) <- obj[[pick]][,1]; obj$seurat_clusters <- Idents(obj)
message("working resolution: ", pick, " -> ", nlevels(Idents(obj)), " clusters")

covars <- intersect(c("sample","animal_timepoint","timepoint","animal_genotype",
                      "animal_sex","mecp2_allele_label"), colnames(obj[[]]))
pdf(fig("06_UMAP_covariates.pdf"), width = 15, height = 5*ceiling((length(covars)+1)/3))
print(wrap_plots(c(
  list(DimPlot(obj, reduction="umap", label=TRUE, raster=TRUE) + ggtitle("clusters") + NoLegend()),
  lapply(covars, function(c_) DimPlot(obj, reduction="umap", group.by=c_, raster=TRUE) +
           ggtitle(c_) + theme(legend.text = element_text(size=6)))), ncol = 3))
dev.off()

# experimental design: timepoints are NOT crossed with sex
#   males P30/P60/P120, females P30/P60/P150 -> late timepoint is sex-specific,
#   so stratify by sex rather than pooling across it.
tp <- intersect(c("animal_timepoint","timepoint"), colnames(obj[[]]))[1]
if (!is.na(tp) && "animal_sex" %in% colnames(obj[[]])) {
  cat("\n-- design (cells per sex x timepoint) --\n")
  print(table(sex = obj[[]][["animal_sex"]], timepoint = obj[[]][[tp]]))
  pdf(fig("06b_UMAP_by_sex_timepoint.pdf"), width = 14, height = 8)
  print(DimPlot(obj, reduction="umap", group.by=tp, split.by="animal_sex",
                raster=TRUE) + ggtitle("timepoint, split by sex"))
  dev.off()
}

# batch-mixing check: does any cluster come from a single sample?
comp <- as.data.frame(table(cluster = Idents(obj), sample = obj$sample)) %>%
  group_by(cluster) %>% mutate(frac = Freq/sum(Freq))
pdf(fig("07_cluster_composition.pdf"), width = 12, height = 6)
print(ggplot(comp, aes(cluster, frac, fill = sample)) + geom_col() +
      theme_classic() + theme(legend.position="none",
        axis.text.x = element_text(angle=90, size=6)) +
      labs(y="fraction of cluster", title="Cluster composition by sample (flat = well mixed)"))
dev.off()
worst <- comp %>% group_by(cluster) %>% summarise(top = max(frac)) %>% arrange(desc(top))
cat("\nclusters most dominated by one sample (possible batch/debris):\n"); print(head(worst, 10))

pdf(fig("08_QC_on_umap.pdf"), width = 15, height = 5)
print(FeaturePlot(obj, features = qc, reduction="umap", raster=TRUE, ncol=length(qc)))
dev.off()

saveRDS(obj, outf)
message("saved ", outf, "; figures in ", figdir)
