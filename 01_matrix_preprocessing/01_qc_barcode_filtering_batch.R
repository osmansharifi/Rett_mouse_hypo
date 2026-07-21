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
DATA_DIR <- "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
FIG_DIR  <- "/quobyte/lasallegrp/Osman/shenyu/03_figures"
seu <- readRDS(file.path(DATA_DIR, "seu_metadata.rds"))

# metadata: study design
head(seu)
unique(seu$orig.ident) 
table(unique(seu@meta.data[, c("orig.ident", "sex", "time_point")])[, c("sex", "time_point")])


### 1. Pre-filtering Exploratory Analysis --------------------------------------
summary(seu$nFeature_RNA)
summary(seu$nCount_RNA)
summary(seu$percent.mt)

# Sample Total UMI Count Distribution
df <- seu@meta.data %>%
  group_by(orig.ident) %>%
  summarise(Total_UMI = sum(nCount_RNA))
p <- ggplot(df, aes(x = reorder(orig.ident, -Total_UMI), y = Total_UMI)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Total UMI Count per Sample",
       x = "Sample",
       y = "Total UMI Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggplot2::ggsave(file.path(FIG_DIR, "Sample_UMI_Count.pdf"), plot = p, height = 8.5, width = 10)

# Plot QC Covariate Distributions (scCustomize)
Idents(seu) <- "scRTT_hypo"
qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
lapply(qc_features, function(feature) {
  p <- QC_Histogram(seu, features = feature)
  file_name <- glue("Histogram_{feature}.pdf")
  ggplot2::ggsave(file_name, plot = p, height = 8.5, width = 10)
  message(glue("Successfully saved: {file_name}"))
})

# Violin plot (Seurat)
p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "violinplot_RNA.pdf"), plot = p, height = 8.5, width = 10)

# Pre-filter Dependencies (Seurat)
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = NULL, raster=FALSE) # r = -.07
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = NULL, raster=FALSE) # r = .98
p_scatter <- plot1 + plot2
ggplot2::ggsave(file.path(FIG_DIR, "featurescatter_RNA.pdf"), plot = p_scatter, height = 8.5, width = 10)


### 2. Barcode Filtering (Adaptive & Hard Thresholds) --------------------------

# Explore MAD-based thresholds (scater)
qc.nCount_RNA <- isOutlier(seu$nCount_RNA, log = TRUE, type = "both", batch = seu$orig.ident, nmads = 5)
umi_lower_limit <- attr(qc.nCount_RNA, "thresholds")["lower", ]
umi_upper_limit <- attr(qc.nCount_RNA, "thresholds")["higher", ]

qc.nFeature_RNA <- isOutlier(seu$nFeature_RNA, log = TRUE, type = "both", batch = seu$orig.ident, nmads = 5)
gene_lower_limit <- attr(qc.nFeature_RNA, "thresholds")["lower", ]
gene_upper_limit <- attr(qc.nFeature_RNA, "thresholds")["higher", ]

# Execute Filtering
seu$discard <- qc.nCount_RNA | qc.nFeature_RNA | (seu$percent.mt > 1)
seu_filtered <- subset(seu, subset = !discard)
head(seu)

# Global threshold among 28 samples
global_gene_low  <- min(gene_lower_limit, na.rm = TRUE) # 48.99
global_gene_high <- max(gene_upper_limit, na.rm = TRUE) # 75170.59
global_umi_low   <- min(umi_lower_limit, na.rm = TRUE) # 27.54
global_umi_high  <- max(umi_upper_limit, na.rm = TRUE) # 298097.92

global_thresh_df <- data.frame(
  Metric = c("global_gene_low", "global_gene_high", "global_umi_low", "global_umi_high"),
  Global_Value = c(
    min(gene_lower_limit, na.rm = TRUE),
    max(gene_upper_limit, na.rm = TRUE),
    min(umi_lower_limit, na.rm = TRUE),
    max(umi_upper_limit, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
print(global_thresh_df)
write.csv(global_thresh_df, file = file.path(FIG_DIR, "barcode_qc_globalthreshold.csv"), row.names = FALSE)

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
final_summary <- df %>%
  left_join(thresh_df, by = "orig.ident") %>%
  left_join(qc_summary, by = "orig.ident") %>%
  arrange(desc(Total_UMI))
write.csv(final_summary, file.path(FIG_DIR, "barcode_qc_summary.csv"), row.names = FALSE)
print(head(final_summary))
sum(final_summary$Total_Discard_Cells) # 193


### 3. Post-filtering QC Re-assessment -------------------------

# Violin plot (Seurat)
p <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "discard", raster = FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "violinplot_RNA_filtered.pdf"), plot = p, height = 8.5, width = 10)

# UMI vs Gene (scCustomize)
p1 <- QC_Plot_UMIvsGene(seu, 
                        x_axis_label = "UMIs per Nucleus", y_axis_label = "Genes per Nucleus",
                        group.by = "discard",
                        raster = FALSE) +
  ggplot2::labs(title = "Genes vs. UMIs per Nucleus")
ggplot2::ggsave(file.path(FIG_DIR, "QC_Plot_UMIvsGene_raw.pdf"), plot = p1, height = 8.5, width = 10)

p1 <- QC_Plot_UMIvsGene(seu_filtered, 
                        x_axis_label = "UMIs per Nucleus", y_axis_label = "Genes per Nucleus",
                        group.by = "discard",
                        raster = FALSE) +
  ggplot2::labs(title = "Genes vs. UMIs per Nucleus")
ggplot2::ggsave(file.path(FIG_DIR, "QC_Plot_UMIvsGene_filtered.pdf"), plot = p1, height = 8.5, width = 10)

# percent.mt vs Gene (scCustomize)
# r = -.09
p2 <- QC_Plot_GenevsFeature(seu, 
                            feature1 = "percent.mt", 
                            high_cutoff_feature = 1, 
                            group.by = "discard",
                            y_axis_label = "Genes per Nucleus", x_axis_label = "percent.mt per Nucleus",
                            raster = FALSE) 
ggplot2::ggsave(file.path(FIG_DIR, "QC_Plot_percent.mtvsGene_raw.pdf"), plot = p2, height = 8.5, width = 10)

# r = -.12
p2 <- QC_Plot_GenevsFeature(seu_filtered, 
                            feature1 = "percent.mt", 
                            high_cutoff_feature = 1, 
                            group.by = "discard",
                            y_axis_label = "Genes per Nucleus", x_axis_label = "percent.mt per Nucleus",
                            raster = FALSE) 
ggplot2::ggsave(file.path(FIG_DIR, "QC_Plot_percent.mtvsGene_filtered.pdf"), plot = p2, height = 8.5, width = 10) 


### 4. Sex difference in Median number of Genes and UMIs per Nucleus -------------------------------------------
p_sex_gene <- Plot_Median_Genes(seu_filtered, group.by = "sex")
ggplot2::ggsave(file.path(FIG_DIR, "Plot_Median_Genes_sex_filtered.pdf"), plot = p_sex_gene, height = 8.5, width = 10)

p_sex_umi <- Plot_Median_UMIs(seu_filtered, group.by = "sex")
ggplot2::ggsave(file.path(FIG_DIR, "Plot_Median_UMIs_sex_filtered.pdf"), plot = p_sex_umi, height = 8.5, width = 10)


### 5. Save Cleaned Data -------------------------------------------------------
Idents(seu_filtered)
saveRDS(seu_filtered, file = file.path(DATA_DIR,"seu_filtered_batch.rds"))