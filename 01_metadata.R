library(Seurat)
library(biomaRt)
library(gprofiler2) # ‘0.2.4’
packageVersion("Seurat") # ‘5.3.0’

getwd() # "/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
load(all_hypo_soupx.RData)
# all_hypo_soupx <- CreateSeuratObject(counts = all_hypo_soupx, project = "all_hypo_soupx", min.cells = 3, min.features = 200)
head(all_hypo_soupx)
gene_cell_counts <- rowSums(GetAssayData(all_hypo_soupx, layer = "counts") > 0)
min(gene_cell_counts) # 10
min(all_hypo_soupx$nFeature_RNA) # 200


## orig.ident: sample name
sample <- sub("^(([^_]+_){3}[^_]+).*", "\\1", row.names(all_hypo_soupx@meta.data))
all_hypo_soupx$orig.ident <- sample
head(all_hypo_soupx)

## sex
all_hypo_soupx$sex <- ifelse(grepl("_M_", all_hypo_soupx$orig.ident), "M",
                             ifelse(grepl("_F_", all_hypo_soupx$orig.ident), "F", NA))
head(all_hypo_soupx)

## genotype
all_hypo_soupx$genotype <- sub("_.*", "", all_hypo_soupx$orig.ident)
head(all_hypo_soupx)

## time_point
all_hypo_soupx$time_point <- regmatches(all_hypo_soupx$orig.ident, regexpr("P\\d+", all_hypo_soupx$orig.ident))
head(all_hypo_soupx)
all_hypo_soupx$time_point <- factor(
  all_hypo_soupx$time_point, 
  levels = c("P30", "P60", "P120", "P150")
)

## number of genes detected per UMI: this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
# merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

## percent.mt
all_hypo_soupx [["percent.mt"]] <- PercentageFeatureSet(all_hypo_soupx, pattern = "^mt-")
head(all_hypo_soupx)

## cell cycle
head(all_hypo_soupx@assays)
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

all_hypo_soupx <- CellCycleScoring(all_hypo_soupx, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
head(all_hypo_soupx)

saveRDS(all_hypo_soupx, file = "all_hypo_soupx_metadata.rds")