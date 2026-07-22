#!/usr/bin/env Rscript

# Environment Setup and Dependencies
library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)
library(biomaRt)
DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures"

# Set path and load data
seu <- readRDS(file.path(DATA_DIR, "seu_filtered_feature.rds"))
head(seu@assays) # Layers: counts


# 1. Shifted logarithm normalization ----------------------------------------

# inspect the distribution of the raw counts which we already calculated during quality control
QC_Histogram(seu, features = "nCount_RNA")
seu <- NormalizeData(seu) # Default: normalization.method = "LogNormalize", scale.factor = 10000
head(seu@assays) # Layers: counts, data

# inspect the distribution of the Shifted logarithm counts reflecting 
seu$nCount_RNA_normalized <- colSums(GetAssayData(seu, layer = "data"))
p <- QC_Histogram(seu, features = "nCount_RNA_normalized", plot_title = "Shifted logarithm")
ggplot2::ggsave(file.path(FIG_DIR, "Histogram_nCount_RNA_normalized.pdf"), plot = p, height = 8.5, width = 10)

saveRDS(seu, file = file.path(DATA_DIR, "seu_norm.rds"))


# 2. Cell Cycle Scoring ----------------------------------------

# Seurat's built-in human cell cycle markers from Tirosh et al, 2015
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# create compiled cell cycle gene lists for mouse homologues 

# gprofiler2::gorth for Orthology search
# https://github.com/satijalab/seurat/issues/2493
m.s.genes = gorth(cc.genes$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
m.g2m.genes = gorth(cc.genes$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
head(m.s.genes) 
head(m.g2m.genes)

s_lost   <- setdiff(cc.genes$s.genes, gorth(cc.genes$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$input)
g2m_lost <- setdiff(cc.genes$g2m.genes, gorth(cc.genes$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$input)
cat("S phase lost (", length(s_lost), "):\n", paste(s_lost, collapse = ", "), "\n\n")
cat("G2M phase lost (", length(g2m_lost), "):\n", paste(g2m_lost, collapse = ", "), "\n\n")

# score each gene for cell cycle phase
# https://github.com/satijalab/seurat/issues/5880
seu <- CellCycleScoring(seu, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
head(seu[[]])
table(seu$Phase)
summary(seu$S.Score)
summary(seu$G2M.Score)

# Visualize the distribution of S and G2M cell cycle markers expression
c("Pcna", "Top2a", "Mcm6", "Mki67") %in% rownames(seu)
p <- RidgePlot(seu, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2) # layer = "data"
ggplot2::ggsave(file.path(FIG_DIR, "RidgePlot_CellCycle.pdf"), plot = p, height = 8.5, width = 10)

# Evaluate the biological relevance of cell cycle scores
p <- QC_Histogram(seu, features = "S.Score")
ggplot2::ggsave(file.path(FIG_DIR, "histogram_S.Score.pdf"), plot = p, height = 8.5, width = 10)
p <- QC_Histogram(seu, features = "G2M.Score")
ggplot2::ggsave(file.path(FIG_DIR, "histogram_G2M.Score.pdf"), plot = p, height = 8.5, width = 10)
p <- FeatureScatter(seu, feature1 = "S.Score", feature2 = "Mcm6")
ggplot2::ggsave(file.path(FIG_DIR, "scatter_S.Score_Mcm6.pdf"), plot = p, height = 8.5, width = 10)

# Phase does not show clear separation in PCA space
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = c(m.s.genes, m.g2m.genes)) # group.by = "Phase"
p <- DimPlot(seu)
ggplot2::ggsave(file.path(FIG_DIR, "pca_exprmarker_phase.pdf"), plot = p, height = 8.5, width = 10)

sdev <- seu[["pca"]]@stdev
var_explained <- (sdev^2) / sum(sdev^2)
var_explained[1:2] # 0.03353326 0.02749839


# 3. Feature Selection  ----------------------------------------

# HVG by Seurat: Genes with a higher-than-expected variance
nfeature = 2000
seu <- FindVariableFeatures(seu, nfeatures = nfeature) # selection.method = "vst" by default

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu), 10)
hvg <- HVFInfo(seu) %>% 
  tibble::rownames_to_column(var = "Gene")
hvg$is_HVG <- hvg$Gene %in% VariableFeatures(seu)
hvg <- hvg %>% 
  dplyr::arrange(desc(variance.standardized))
write.csv(hvg, file = file.path(FIG_DIR, glue::glue("HVG_variance_{nfeature}.csv")), row.names = FALSE)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p <- plot1 + plot2
ggplot2::ggsave(file.path(FIG_DIR, glue::glue("HVG_{nfeature}.pdf")), plot = p, height = 8.5, width = 10)


# 4. Standardization & Regression ----------------------------------------
seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = VariableFeatures(seu))


# Save Scaled Data -------------------------------------------------------
seu[["pca"]] <- NULL
Reductions(seu)
saveRDS(seu, file = file.path(DATA_DIR, "seu_scaled.rds"))


# biomaRt::getLDS ----------------------------------------
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
m.s.genes <- convertHumanGeneList(s.genes) # HTTP 500 Internal Server Error.
m.g2m.genes <- convertHumanGeneList(g2m.genes)