#!/usr/bin/env Rscript

# Environment Setup and Dependencies
library(Seurat)
library(scCustomize)
DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures"

# Set path and load data
seu <- readRDS(file.path(DATA_DIR, "seu_scaled.rds"))
head(seu@assays) # Layers: counts, data, scale.data


# 1. PCA ----------------------------------------

# execute dimensionality reduction on the reduced feature space
seu <- RunPCA(seu) # Default: npcs = 50, features= VariableFeatures(seu) for the scaled Assay

# scree graph plots of PoV vs k, stop at “elbow”
p <- ElbowPlot(seu, ndims = 50, reduction = "pca")
ggplot2::ggsave(file.path(FIG_DIR, "ElbowPlot_50pc.pdf"), plot = p, height = 8.5, width = 10)

# calculate Proportion of Variance (PoV) explained
pov <- (seu[["pca"]]@stdev^2) / seu[["pca"]]@misc$total.variance
cum_pov <- cumsum(pov)
df_pov <- data.frame(PC = 1:length(cum_pov), cum_pov = cum_pov)
write.csv(df_pov, file = file.path(FIG_DIR, "pov_50pc.csv"), row.names = FALSE)
p <- ggplot(df_pov, aes(x = PC, y = cum_pov)) +
  geom_line(color = "blue") +
  geom_point(color = "blue", shape = 16) +
  geom_hline(yintercept = 0.9, color = "red", linetype = "dashed") +
  labs(x = "Principal Component (k)", y = "PoV Explained", title = "Proportion of Variance (PoV)") +
  theme_bw() 
ggplot2::ggsave(file.path(FIG_DIR, "pov_50pc.pdf"), plot = p, height = 8.5, width = 10)

# Examine and visualize PCA results a few different ways
# Visualize top genes associated with reduction components
p <- VizDimLoadings((seu), dims = 1:2, reduction = "pca") 
ggplot2::ggsave(file.path(FIG_DIR, "topgenes_2pc.pdf"), plot = p, height = 8.5, width = 10)

Idents(seu) <- "old.ident"
p <- DimPlot(seu, reduction = "pca") + NoLegend()

p <- DimHeatmap(seu, dims = 1:20, cells = 500, balanced = TRUE)
ggplot2::ggsave(file.path(FIG_DIR, "heatmap_20pc.pdf"), plot = p, height = 8.5, width = 10)


# Inspecting quality control metrics ----------------------------------------
p <- FeaturePlot(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), reduction = "pca", ncol = 3, raster=FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "qc_pca_barcodeqc.pdf"), plot = p, height = 8.5, width = 10)
p <- FeaturePlot(seu, features = c("S.Score", "G2M.Score"), reduction = "pca", ncol = 2, raster=FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "qc_pca_cellcycleregression.pdf"), plot = p, height = 8.5, width = 10)

p1 <- DimPlot(seu, group.by = "Phase", reduction = "pca", raster=FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "qc_pca_phase.pdf"), plot = p1, height = 8.5, width = 10)

p2 <- DimPlot(seu, group.by = "sex", reduction = "pca", raster=FALSE)
p3 <- DimPlot(seu, group.by = "genotype", reduction = "pca", raster=FALSE)
p4 <- DimPlot(seu, group.by = "time_point", reduction = "pca", raster=FALSE)
p <- p2+p3+p4
ggplot2::ggsave(file.path(FIG_DIR, "qc_pca_metadata.pdf"), plot = p, height = 8.5, width = 10)