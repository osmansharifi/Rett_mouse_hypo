library(Seurat)
library(biomaRt)
library(gprofiler2) # ‘0.2.4’
packageVersion("Seurat") # ‘5.3.0’

setwd("/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects")
seu <- readRDS("mecp2_H_master.rds")
attributes(seu) # Assay (v5) data with 20811 features for 100375 cells
# seu <- CreateSeuratObject(seu, min.cells = 0, min.features = 0)
gene_cell_counts <- rowSums(GetAssayData(seu, layer = "counts") > 0)
min(gene_cell_counts) # min.cells = 3
min(seu$nFeature_RNA) # min.features = 225
head(seu)

## orig.ident: metadata mapping
setwd("/quobyte/lasallegrp/Osman/shenyu/01_raw")
hypo_sample_metadata <- read.csv("hypo_sample_metadata.csv")
head(hypo_sample_metadata)
match_key <- paste0(gsub("-", "_", hypo_sample_metadata$mouse), "_H")
map_dict <- setNames(hypo_sample_metadata$sample_name, match_key)
new_idents <- map_dict[as.character(seu$orig.ident)]
new_idents <- unname(new_idents)  
seu$orig.ident <- new_idents
Idents(seu) <- "orig.ident"
old_prefix <- paste0(gsub("-", "_", hypo_sample_metadata$mouse), "_H") 
new_prefix <- hypo_sample_metadata$sample_name
old_cell_names <- colnames(seu)
new_cell_names <- old_cell_names
for (i in seq_along(old_prefix)) {
  new_cell_names <- gsub(old_prefix[i], new_prefix[i], new_cell_names)
}
seu <- RenameCells(seu, new.names = new_cell_names)

## sex
seu$sex <- ifelse(grepl("_M_", seu$orig.ident), "M",
                  ifelse(grepl("_F_", seu$orig.ident), "F", NA))
head(seu)

## genotype
seu$genotype <- sub("_.*", "", seu$orig.ident)
head(seu)

## time_point
seu$time_point <- regmatches(seu$orig.ident, regexpr("P\\d+", seu$orig.ident))
head(seu)
seu$time_point <- factor(
  seu$time_point, 
  levels = c("P30", "P60", "P120", "P150")
)

## number of genes detected per UMI: this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
# merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

## percent.mt
seu [["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
head(seu)

setwd("/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects")
saveRDS(seu, file = "seu_metadata.rds")