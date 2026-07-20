#!/usr/bin/env Rscript

# Environment Setup and Dependencies
library(SingleCellExperiment)
library(Seurat)
library(scater)
library(ggplot2)

DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures"

# Set path and load data
seu <- readRDS(file.path(DATA_DIR, "seu_filtered_batch.rds"))
seu # 20811 features across 100182 samples within 1 assay 


# 1. Calculate Feature-level QC Metrics with two ecosystems as cross-validation ----------------------------------------

# SeuratObject
gene_cell_counts <- rowSums(GetAssayData(seu, layer = "counts") > 0)
summary(gene_cell_counts) # min = 2

# SingleCellExperiment
sce <- as.SingleCellExperiment(seu) 
feature <- perFeatureQCMetrics(sce)
head(feature)
# mean: the mean counts for each feature.
# detected: the percentage of observations above threshold (default = 0).


# 2. Define Feature Filtering Threshold ------------------------------------------------

# filter out genes that are not detected in at least 20 cells
cutoff_cell_num <- 20
feature_removed <- names(gene_cell_counts[gene_cell_counts < cutoff_cell_num])
length(feature_removed) # 2080
head(feature_removed)

feature_threshold <- (cutoff_cell_num / ncol(seu)) * 100
feature_filtered <- rownames(feature)[feature$detected < feature_threshold]
length(feature_filtered) # 2080
head(feature_filtered)

# Plot Gene Prevalence Statistics
setwd("/quobyte/lasallegrp/Osman/shenyu/03_figures")
p <- ggplot(feature, aes(x = detected)) + geom_histogram() + 
  labs(x = "Gene Detection Rate per Cell (%)", y = "Number of Genes", 
       title = "Distribution of Gene Prevalence Across Cells",
       subtitle = paste0(
         "Red dashed line indicates quality control threshold (", round(feature_threshold, 3), "%)\n",
         "Filtered out ", length(feature_filtered), " genes detected in fewer than ", cutoff_cell_num, " cells"
       )) + 
  geom_vline(xintercept = feature_threshold, color = "red")
ggplot2::ggsave(file.path(FIG_DIR, "Histogram_GenePrevalence.pdf"), plot = p, height = 8.5, width = 10)


# 3. Execute Feature Filtering  ------------------------------------------------
identical(feature_removed,feature_filtered) # TRUE
seu_filtered <- subset(seu, features = rownames(feature)[feature$detected >= feature_threshold])
seu_filtered # 18731 features across 100182 samples within 1 assay 

saveRDS(seu_filtered, file = file.path(DATA_DIR, "seu_filtered_feature.rds"))