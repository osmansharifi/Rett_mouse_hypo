#!/usr/bin/env Rscript
# Build nested cell-type labels: keep a cell's C66 label only where it is BOTH
# confident (score >= cutoff) AND consistent with its C25 parent in the HypoMap
# hierarchy; otherwise fall back to the (coarser, more reliable) C25 label.
#
# Produces these metadata columns on the output object:
#   celltype_c25        C25 transferred label            (+ score_c25)
#   celltype_c66        C66 transferred label            (+ score_c66)
#   neuron_class        E/I collapse: GLU / GABA / Non-neuronal   (coarse, trivially confident)
#   class_broad         glia kept as-is, all neurons -> Neuron    (backbone tier)
#   celltype_hier       per-cell "deepest confident + consistent" label  <-- main label
#   hier_level          which level celltype_hier came from: C66 | C25
#
# Usage: Rscript 08_hier_labels.R <c25_obj.rds> <c66_obj.rds> <c25c66_map.csv> \
#                                 <out.rds> <figdir> [c66_cutoff]
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
f25 <- args[1]; f66 <- args[2]; mapf <- args[3]
outf <- args[4]; figdir <- args[5]
cutoff <- ifelse(length(args) >= 6, as.numeric(args[6]), NA)   # NA -> data-driven below
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fig <- function(n) file.path(figdir, n)

o25 <- readRDS(f25); o66 <- readRDS(f66)
stopifnot(setequal(colnames(o25), colnames(o66)))   # same cells (order handled by name-match)

# ---- pull the two transfers onto one object (use the C66 object as the base) ----
obj <- o66
# align C25 fields onto obj by CELL NAME (build named vectors so Seurat matches
# on names, not position; avoids the "no cell overlap" quirk in v5 $<- ).
i <- match(colnames(obj), colnames(o25))
stopifnot(!anyNA(i))                              # every obj cell must exist in o25
v_ct <- as.character(o25$celltype)[i];        names(v_ct) <- colnames(obj)
v_sc <- as.numeric(o25$celltype_score)[i];    names(v_sc) <- colnames(obj)
obj <- AddMetaData(obj, v_ct, col.name = "celltype_c25")
obj <- AddMetaData(obj, v_sc, col.name = "score_c25")
obj$celltype_c66 <- as.character(obj$celltype)
obj$score_c66    <- as.numeric(obj$celltype_score)
rm(o25, o66); gc()

# ---- data-driven C66 cutoff if not supplied: the antimode between the two peaks,
#      searched in the plausible 0.3-0.7 band; falls back to 0.5 ----
if (is.na(cutoff)) {
  d <- density(obj$score_c66[!is.na(obj$score_c66)], n = 512)
  lo <- 0.3; hi <- 0.9
  win <- which(d$x > lo & d$x < hi)
  cutoff <- if (length(win)) round(d$x[win][which.min(d$y[win])], 2) else 0.5
  if (cutoff <= lo + 0.02 || cutoff >= hi - 0.02)
    message("WARNING: antimode hit the search boundary (", cutoff,
            ") - no clear trough. Inspect 16_c66_score_cutoff.pdf and consider ",
            "passing an explicit cutoff as the 6th argument.")
}
message("C66 confidence cutoff: ", cutoff)

# ---- reference hierarchy: which C66 types nest under which C25 class ----
hmap <- fread(mapf)                       # cols: C25_named, C66_named
setnames(hmap, c("C25_named","C66_named"), c("c25_parent","c66_child"))
parent_of <- setNames(hmap$c25_parent, hmap$c66_child)   # C66 -> its true C25 parent

# ---- consistency: does the cell's predicted C66 nest under its predicted C25? ----
# Compute on PLAIN VECTORS (pulled out of the object), then write finished columns
# back once via a named data.frame -> AddMetaData. Avoids the v5 $<- "no cell
# overlap" error that named-index lookups (parent_of[...]) trigger.
cells    <- colnames(obj)
ct_c25   <- as.character(obj$celltype_c25)
ct_c66   <- as.character(obj$celltype_c66)
sc_c66   <- as.numeric(obj$score_c66)
# PRIMARY (authoritative): does the HypoMap table say this C66 type's parent is
# exactly the C25 class we predicted for this cell? Full-string comparison.
parent_map <- unname(parent_of[ct_c66])
primary_ok <- !is.na(parent_map) & parent_map == ct_c25

# SECONDARY guard, NEURONS ONLY: the "...GABA-1"/"...GLU-2" suffix must agree with
# the parent's class token. Non-neuronal labels (Astrocytes, OPC, Immune, ...) carry
# no such token, so the guard is simply not applicable to them and must NOT fail them.
ei_token <- function(x) {
  x <- as.character(x); out <- rep(NA_character_, length(x))
  hit <- regexpr("(GLU|GABA)-[0-9]+$", x)
  ok <- hit > 0
  out[ok] <- substring(x[ok], hit[ok])
  out
}
tok_parent <- ei_token(parent_map)
tok_c66    <- ei_token(ct_c66)
suffix_applicable <- !is.na(tok_parent) & !is.na(tok_c66)
suffix_ok <- !suffix_applicable | (tok_parent == tok_c66)   # TRUE when not applicable

consistent <- primary_ok & suffix_ok
confident  <- !is.na(sc_c66) & sc_c66 >= cutoff
n_drift <- sum(suffix_applicable & !suffix_ok, na.rm = TRUE)
if (n_drift > 0) message("NOTE: ", n_drift,
  " neuronal cells where map-parent and name-suffix disagree -> treated inconsistent")
cat("consistency detail: primary_ok =", sum(primary_ok),
    " | suffix guard applicable =", sum(suffix_applicable),
    " | final consistent =", sum(consistent), "\n")

# ---- the main label: deepest that is confident AND consistent, else C25 ----
keep_c66     <- confident & consistent
celltype_hier<- ifelse(keep_c66, ct_c66, ct_c25)
hier_level   <- ifelse(keep_c66, "C66", "C25")

# ---- coarse tiers ----
neuron_class <- ifelse(grepl("GLU",  ct_c25), "GLU",
                ifelse(grepl("GABA", ct_c25), "GABA", "Non-neuronal"))
class_broad  <- ifelse(neuron_class %in% c("GLU","GABA"), "Neuron", ct_c25)

# ---- write all derived columns back at once, keyed by cell name ----
add <- data.frame(row.names       = cells,
                  hier_consistent = consistent,
                  hier_confident  = confident,
                  celltype_hier   = celltype_hier,
                  hier_level      = hier_level,
                  neuron_class    = neuron_class,
                  class_broad     = class_broad,
                  stringsAsFactors = FALSE)
obj <- AddMetaData(obj, add)

# ---------------------------------- figures ---------------------------------------
pdf(fig("16_c66_score_cutoff.pdf"), width = 8, height = 5)
print(ggplot(obj[[]], aes(score_c66)) + geom_histogram(bins = 60) +
      geom_vline(xintercept = cutoff, colour = "red", linetype = 2) +
      theme_classic() + labs(title = paste("C66 score, cutoff =", cutoff)))
dev.off()

# where did the final label come from, and why C66 was dropped
lvl <- obj[[]] %>% mutate(reason = case_when(
          hier_level == "C66" ~ "kept C66",
          !hier_confident & !hier_consistent ~ "C25: low score & inconsistent",
          !hier_confident ~ "C25: low score",
          !hier_consistent ~ "C25: inconsistent w/ parent",
          TRUE ~ "C25: other"))
pdf(fig("17_hier_level_breakdown.pdf"), width = 8, height = 5)
print(ggplot(lvl, aes(reason, fill = reason)) + geom_bar() + coord_flip() +
      theme_classic() + NoLegend() + labs(title = "why each cell's final label level"))
dev.off()

if ("umap" %in% Reductions(obj)) {
  pdf(fig("18_UMAP_hier_label.pdf"), width = 15, height = 9)
  print(DimPlot(obj, group.by = "celltype_hier", label = TRUE, repel = TRUE,
                raster = TRUE, label.size = 2.5) + NoLegend() +
        ggtitle("celltype_hier (deepest confident + consistent)"))
  dev.off()
  pdf(fig("19_UMAP_label_tiers.pdf"), width = 18, height = 6)
  print((DimPlot(obj, group.by="class_broad",  label=TRUE, raster=TRUE) + NoLegend() + ggtitle("class_broad")) |
        (DimPlot(obj, group.by="neuron_class", raster=TRUE) + ggtitle("neuron_class (E/I)")) |
        (DimPlot(obj, group.by="hier_level",   raster=TRUE) + ggtitle("label level")))
  dev.off()
}

saveRDS(obj, outf)
cat("\n=== label-level summary ===\n"); print(table(obj$hier_level, useNA="ifany"))
cat("\n=== consistency of C66 calls ===\n")
print(table(confident = obj$hier_confident, consistent = obj$hier_consistent))
cat("\n=== cells per final (hier) label ===\n")
print(sort(table(obj$celltype_hier), decreasing = TRUE))
cat("\n=== neuron_class ===\n"); print(table(obj$neuron_class))
message("\nsaved ", outf)
