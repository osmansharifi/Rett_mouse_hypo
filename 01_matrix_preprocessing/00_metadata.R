library(Seurat)
library(biomaRt)
library(gprofiler2) # ‘0.2.4’
packageVersion("Seurat") # ‘5.3.0’

getwd() # "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
load("hypo_adult_unfiltered.RData")
seu <- CreateSeuratObject(counts = seu, project = "Rett_hypothalamus", min.cells = 0, min.features = 0)
gene_cell_counts <- rowSums(GetAssayData(seu, layer = "counts") > 0)
min(gene_cell_counts) # 0
min(seu$nFeature_RNA) # 0
head(seu)

## orig.ident: metadata mapping
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

## cell cycle
seu <- NormalizeData(seu)
head(seu@assays)
# A list of cell cycle markers from Tirosh et al, 2015
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# Convert Seurat's built-in human cell cycle gene lists to their mouse homologues 
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
# m.s.genes <- convertHumanGeneList(s.genes)
# m.g2m.genes <- convertHumanGeneList(g2m.genes)

# gprofiler2::gorth for Orthology search
# https://github.com/satijalab/seurat/issues/2493
m.s.genes = gorth(cc.genes$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
m.g2m.genes = gorth(cc.genes$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
head(m.s.genes) # "Mcm5" "Pcna" "Tyms" "Fen1" "Mcm2" "Mcm4"
head(m.g2m.genes) # "Hmgb2"  "Cdk1"   "Nusap1" "Birc5"  "Tpx2"   "Top2a"

seu <- CellCycleScoring(seu, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
head(seu)

setwd("/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects")
saveRDS(seu, file = "seu_metadata.rds") # seu with metadata & the normalized data layer