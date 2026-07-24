#!/usr/bin/env Rscript
# Build ONE analysis-ready `cell_type` column for DEG + mosaicism.
#
# Stage A (structural, always): collapse every C66 label to its C25 parent. This
#   removes confidence artifacts -- e.g. "Ermn.Oligodendrocytes" (confident C66) vs
#   "C25-19: Oligodendrocytes" (same cells, just below the C66 score cutoff) are the
#   same biology split by how sure the transfer was, not by cell identity.
#
# Stage B (adaptive, count-driven): a group survives only if enough animals carry
#   enough cells in BOTH alleles in the limiting stratum (het females with allele
#   calls). Survivors keep their C25 identity; the rest pool into GLU-other /
#   GABA-other / Other-glia / Other-nonneuronal.
#
# Usage: Rscript 11_make_celltype.R <obj.rds> <c25c66_map.csv> <out.rds> <outdir> \
#                                   [min_cells] [min_pairs]
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
inf <- args[1]; mapf <- args[2]; outf <- args[3]; outdir <- args[4]
MIN_CELLS <- ifelse(length(args) >= 5, as.integer(args[5]), 10)
MIN_PAIRS <- ifelse(length(args) >= 6, as.integer(args[6]), 4)
MIN_KEEP  <- ifelse(length(args) >= 7, as.integer(args[7]), 25)  # cells to stay distinct
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(inf); md <- obj[[]]; cn <- colnames(obj)
stopifnot("celltype_hier" %in% names(md))

# ---------- Stage A: C66 -> C25 parent, then strip the "C25-n: " prefix -----------
hmap <- fread(mapf)
setnames(hmap, c("C25_named","C66_named"), c("c25","c66"))
parent_of <- setNames(hmap$c25, hmap$c66)

hier   <- as.character(md$celltype_hier)
c25lab <- ifelse(hier %in% names(parent_of), unname(parent_of[hier]), hier)
tidy   <- function(x) sub("^C\\d+-\\d+:\\s*", "", x)      # "C25-2: GLU-2" -> "GLU-2"
base_type <- tidy(c25lab)

cat("=== Stage A: C66 collapsed to C25 parent ===\n")
cat("distinct labels:", length(unique(hier)), "->", length(unique(base_type)), "\n")

# ---------- Stage B: keep a group distinct if it has enough cells ----------------
# NOTE: group identity and DEG power are different questions. Here we only decide
# whether a group is big enough to stand on its own (mosaicism needs far fewer cells
# than DEG). The DEG script applies its own stricter per-animal paired gate, so
# keeping more groups here costs nothing and preserves resolution.
lim <- which(md$animal_sex == "female" & md$animal_genotype == "mut" &
             md$mecp2_allele %in% c("WT","MUT"))
sz <- as.data.table(list(group = base_type[lim]))[, .N, by = group][order(-N)]
keep_set <- sz[N >= MIN_KEEP, group]

cat("\n=== Stage B: allele-called cells per group (het females) ===\n"); print(sz)
cat("\nkept distinct (>=", MIN_KEEP, "cells):", paste(keep_set, collapse = ", "), "\n")

# pooling: keep the major glial classes SEPARATE (they are transcriptionally
# distinct); only genuinely small/rare groups pool together.
GLIA <- c("Astrocytes","Oligodendrocytes","OPC","Immune","Tanycytes","Ependymal-like")
pool_of <- function(g) {
  ifelse(grepl("^GLU",  g), "GLU-other",
  ifelse(grepl("^GABA", g), "GABA-other",
  ifelse(g %in% GLIA,       "Other-glia", "Other-nonneuronal")))
}
cell_type <- ifelse(base_type %in% keep_set, base_type, pool_of(base_type))

# report which FINAL groups clear the stricter paired-DEG bar (informational only)
dtl <- data.table(group = cell_type[lim], mouse = as.character(md$mouse[lim]),
                  allele = as.character(md$mecp2_allele[lim]))
pu  <- dtl[, .N, by = .(group, mouse, allele)]
pr  <- pu[N >= MIN_CELLS][, .(k = uniqueN(allele)), by = .(group, mouse)][k == 2]
vb  <- pr[, .(n_paired_animals = .N), by = group][order(-n_paired_animals)]
cat("\n=== of the final groups, DEG-viable (>=", MIN_CELLS,
    "cells both alleles in >=", MIN_PAIRS, "animals) ===\n")
print(vb[n_paired_animals >= MIN_PAIRS])
cat("\n(below the bar -- usable for mosaicism, not for per-group DEG)\n")
print(vb[n_paired_animals < MIN_PAIRS])

# ---------- report + write --------------------------------------------------------
mapping <- unique(data.table(celltype_hier = hier, c25 = base_type, cell_type = cell_type))
fwrite(mapping[order(cell_type, c25)], file.path(outdir, "celltype_mapping.csv"))

obj <- AddMetaData(obj, data.frame(row.names = cn,
                                   cell_type_c25 = base_type,
                                   cell_type     = cell_type,
                                   stringsAsFactors = FALSE))
saveRDS(obj, outf)

cat("\n=== final cell_type: all cells ===\n")
print(sort(table(cell_type), decreasing = TRUE))
cat("\n=== final cell_type x allele (het females, the DEG stratum) ===\n")
tb <- table(cell_type = cell_type[lim], allele = md$mecp2_allele[lim])
print(tb[order(-rowSums(tb)), , drop = FALSE])
cat("\n=== per-animal units (cells per mouse x cell_type x allele) ===\n")
dt2 <- data.table(cell_type = cell_type[lim], mouse = as.character(md$mouse[lim]),
                  allele = as.character(md$mecp2_allele[lim]))[, .N, by = .(cell_type, mouse, allele)]
fwrite(dt2[order(cell_type, mouse, allele)], file.path(outdir, "unit_counts_by_celltype.csv"))
print(dt2[, .(animals_with_both = uniqueN(mouse[N >= MIN_CELLS])), by = cell_type][order(-animals_with_both)])

pdf(file.path(outdir, "22_celltype_allele_counts.pdf"), width = 11, height = 6)
d <- as.data.frame(tb)
print(ggplot(d, aes(reorder(cell_type, -Freq), Freq, fill = allele)) +
      geom_col(position = "dodge") + theme_classic() +
      geom_hline(yintercept = MIN_CELLS * MIN_PAIRS, linetype = 2, colour = "grey40") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = NULL, y = "allele-called cells (het females)",
           title = "final cell_type by allele; dashed = rough viability floor"))
dev.off()

message("\nsaved ", outf, "\nmapping + counts in ", outdir)
