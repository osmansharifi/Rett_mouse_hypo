#!/usr/bin/env Rscript

# Environment Setup and Dependencies
library(Seurat)
library(scCustomize)
library(scater)
library(ggplot2)
library(SingleCellExperiment)
library(glue)
library(dplyr)
library(tibble)

# Set path and load data
setwd("/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects")
seu <- readRDS("seu_metadata.rds")

## metadata: study design
head(seu)
unique(seu$orig.ident) 
table(unique(seu@meta.data[, c("orig.ident", "sex", "time_point")])[, c("sex", "time_point")])


setwd("/quobyte/lasallegrp/Osman/shenyu/03_figures")
### 1. Pre-filtering Exploratory Analysis --------------------------------------
summary(seu$nFeature_RNA)
summary(seu$nCount_RNA)
summary(seu$percent.mt)

# Sample Sequencing Depth
depth_df <- seu@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Total_Depth = sum(nCount_RNA))
p_depth <- ggplot(depth_df, aes(x = reorder(orig.ident, -Total_Depth), y = Total_Depth)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Total Sequencing Depth per Sample",
       x = "Sample",
       y = "Total UMI Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggplot2::ggsave("Total_Sequencing_Depth.pdf", plot = p_depth, height = 8.5, width = 10)

# Plot distributions (scCustomize)
Idents(seu) <- "scRTT_hypo"
qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
lapply(qc_features, function(feature) {
  p <- QC_Histogram(seu, features = feature)
  file_name <- glue("Histogram_{feature}.pdf")
  ggplot2::ggsave(file_name, plot = p, height = 8.5, width = 10)
  message(glue("Successfully saved: {file_name}"))
})

# Violin plot
p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
ggplot2::ggsave("violinplot_RNA.pdf", plot = p, height = 8.5, width = 10)

# Pre-filter Dependencies (Scatter)
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = NULL, raster=FALSE) # r = .07
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = NULL, raster=FALSE) # r = .98
p_scatter <- plot1 + plot2
ggplot2::ggsave("featurescatter_RNA.pdf", plot = p_scatter, height = 8.5, width = 10)


### 2. Barcode Filtering (Adaptive & Hard Thresholds) --------------------------

# Explore MAD-based thresholds
qc.nCount_RNA <- isOutlier(seu$nCount_RNA, log = TRUE, type = "both", batch = seu$orig.ident, nmads = 5)
umi_lower_limit <- attr(qc.nCount_RNA, "thresholds")["lower", ]
umi_upper_limit <- attr(qc.nCount_RNA, "thresholds")["higher", ]

qc.nFeature_RNA <- isOutlier(seu$nFeature_RNA, log = TRUE, type = "both", batch = seu$orig.ident, nmads = 5)
gene_lower_limit <- attr(qc.nFeature_RNA, "thresholds")["lower", ]
gene_upper_limit <- attr(qc.nFeature_RNA, "thresholds")["higher", ]

# Execute Filtering
seu_filtered <- subset(seu, 
              nCount_RNA >= umi_lower_limit & nCount_RNA <= umi_upper_limit & 
                nFeature_RNA >= gene_lower_limit & nFeature_RNA <= gene_upper_limit & 
                percent.mt <= 1)
seu$discard <- qc.nCount_RNA | qc.nFeature_RNA | (seu$percent.mt > 1)
head(seu)

# Document cell filtration statistics
thresh_df <- data.frame(
  orig.ident = names(gene_lower_limit),
  Gene_Lower_Limit = as.numeric(gene_lower_limit),
  Gene_Upper_Limit = as.numeric(gene_upper_limit),
  UMI_Lower_Limit  = as.numeric(umi_lower_limit),
  UMI_Upper_Limit  = as.numeric(umi_upper_limit),
  stringsAsFactors = FALSE
)
qc_summary <- seu@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    Discard_Mito_Count = sum(percent.mt > 1, na.rm = TRUE),
    Discard_UMI_Count = sum(nCount_RNA < umi_lower_limit[orig.ident] | nCount_RNA > umi_upper_limit[orig.ident], na.rm = TRUE),
    Discard_Gene_Count = sum(nFeature_RNA < gene_lower_limit[orig.ident] | nFeature_RNA > gene_upper_limit[orig.ident], na.rm = TRUE),
    Total_Discard_Cells = sum(discard, na.rm = TRUE)
  )
final_summary <- depth_df %>%
  left_join(thresh_df, by = "orig.ident") %>%
  left_join(qc_summary, by = "orig.ident") %>%
  arrange(desc(Total_Depth))
print(head(final_summary))
sum(final_summary$Total_Discard_Cells) # 193
write.csv(final_summary, "SampleSequencingDepth_QCSummary.csv", row.names = FALSE)


### 3. Post-filtering QC Re-assessment -------------------------

# Violin plot
p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "discard", raster = FALSE)
ggplot2::ggsave("violinplot_RNA_filtered.pdf", plot = p, height = 8.5, width = 10)

# Global threshold among 28 samples
global_gene_low  <- min(gene_lower_limit, na.rm = TRUE) # 48.99
global_gene_high <- max(gene_upper_limit, na.rm = TRUE) # 75170.59
global_umi_low   <- min(umi_lower_limit, na.rm = TRUE) # 27.54
global_umi_high  <- max(umi_upper_limit, na.rm = TRUE) # 298097.92

# UMI vs Gene with threshold lines
p1 <- QC_Plot_UMIvsGene(seu, 
                        x_axis_label = "UMIs per Nucleus", y_axis_label = "Genes per Nucleus",
                        group.by = "discard",
                        raster = FALSE) +
  ggplot2::labs(title = "Genes vs. UMIs per Nucleus")
ggplot2::ggsave("QC_Plot_UMIvsGene_filtered.pdf", plot = p1, height = 8.5, width = 10)

# percent.mt vs Gene with threshold lines
p2 <- QC_Plot_GenevsFeature(seu, 
                            feature1 = "percent.mt", 
                            high_cutoff_feature = 1, 
                            group.by = "discard",
                            y_axis_label = "Genes per Nucleus", x_axis_label = "percent.mt per Nucleus",
                            raster = FALSE) # # r = -.09
ggplot2::ggsave("QC_Plot_percent.mtvsGene_raw.pdf", plot = p2, height = 8.5, width = 10)
p2 <- QC_Plot_GenevsFeature(seu_filtered, 
                            feature1 = "percent.mt", 
                            high_cutoff_feature = 1, 
                            y_axis_label = "Genes per Nucleus", x_axis_label = "percent.mt per Nucleus",
                            raster = FALSE) # r = -.13
ggplot2::ggsave("QC_Plot_percent.mtvsGene_filtered.pdf", plot = p2, height = 8.5, width = 10) 

# Mito histogram post-filtering
p3 <- QC_Histogram(seu, features = "percent.mt", low_cutoff = 1)
ggplot2::ggsave("QC_Histogram_percent.mt_filtered.pdf", plot = p3, height = 8.5, width = 10)


### 4. Sex difference in Median number of Genes and UMIs per Nucleus -------------------------------------------
p_sex_gene <- Plot_Median_Genes(seu_filtered, group.by = "sex")
ggplot2::ggsave("Plot_Median_Genes_sex_filtered.pdf", plot = p_sex_gene, height = 8.5, width = 10)

p_sex_umi <- Plot_Median_UMIs(seu_filtered, group.by = "sex")
ggplot2::ggsave("Plot_Median_UMIs_sex_filtered.pdf", plot = p_sex_umi, height = 8.5, width = 10)


### 5. Save Cleaned Data -------------------------------------------------------
setwd("/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects")
saveRDS(seu_filtered, file = "seu_filtered_batch.rds")