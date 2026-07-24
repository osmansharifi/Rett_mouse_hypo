#!/usr/bin/env Rscript
# Attach per-cell Mecp2 allele calls to the annotated object and report which DEG
# comparisons the cell counts can actually support.
#
# Join key: the object carries `sample` (19_XXXX_H) and its cell names end in the
# 10x barcode, so key = <sample>_<barcode>, which matches the allele CSVs directly.
#
# Adds:
#   mecp2_wt_umi / mecp2_mut_umi / mecp2_snp_umi   raw SNP-covering UMI counts
#   mecp2_vaf_mut                                   mut / (wt+mut) at the SNP
#   mecp2_snp_purity                                dominant-allele share
#   mecp2_allele                                    WT | MUT | conflicted | NA(no SNP read)
#   mecp2_expr                                      total Mecp2 UMIs (from counts layer)
#   mecp2_allele_expr_wt / _mut                     cell's FULL Mecp2 attributed to its allele
#   sample_ambient_vaf                              per-sample soup MUT fraction (contamination context)
#   deg_group                                       convenience label for modelling
#
# Usage: Rscript 09_attach_alleles.R <obj.rds> <alleledir> <out.rds> <figdir> [min_purity]
suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
inf <- args[1]; alleledir <- args[2]; outf <- args[3]; figdir <- args[4]
BASE <- args[5]                                   # CellRanger dir (for RAW Mecp2)
MIN_PURITY <- ifelse(length(args) >= 6, as.numeric(args[6]), 0.75)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
fig <- function(n) file.path(figdir, n)

obj <- readRDS(inf)
message(sprintf("object: %d cells", ncol(obj)))
stopifnot("sample" %in% colnames(obj[[]]))

# ---------------- read allele tables (cells + AMBIENT soup rows) -------------------
files <- list.files(alleledir, pattern = "\\.mecp2_allele\\.csv$", full.names = TRUE)
message("allele files: ", length(files))
al_all <- rbindlist(lapply(files, function(f) {
  d <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !"barcode" %in% names(d)) return(NULL)
  d
}), fill = TRUE)

# per-sample ambient (soup) VAF: the contamination baseline for interpreting calls
amb <- al_all[barcode == "AMBIENT",
              .(soup_wt = wt_umi, soup_mut = mut_umi), by = sample]
amb[, sample_ambient_vaf := fifelse(soup_wt + soup_mut > 0,
                                    soup_mut / (soup_wt + soup_mut), NA_real_)]

cells <- al_all[barcode != "AMBIENT"]
cells[, ckey := paste0(sample, "_", barcode)]

# ---------------- build the same key on the object --------------------------------
cn      <- colnames(obj)
barcode <- sub(".*_", "", cn)              # trailing 10x barcode, e.g. AAACCTGAGCAGCCTC-1
objkey  <- paste0(as.character(obj$sample), "_", barcode)
i <- match(objkey, cells$ckey)
message(sprintf("matched %d of %d cells to an allele record (%.1f%%)",
                sum(!is.na(i)), length(i), 100 * mean(!is.na(i))))
if (sum(!is.na(i)) == 0)
  stop("no cells matched. Check that obj$sample looks like '19_0136_H' and cell names end in the barcode.")

wt  <- cells$wt_umi[i];  wt[is.na(wt)]   <- 0
mut <- cells$mut_umi[i]; mut[is.na(mut)] <- 0
snp_tot <- wt + mut
purity  <- ifelse(snp_tot > 0, pmax(wt, mut) / snp_tot, NA_real_)
vaf     <- ifelse(snp_tot > 0, mut / snp_tot, NA_real_)

allele <- rep(NA_character_, length(cn))                       # no SNP-covering read
allele[snp_tot > 0 & mut > wt  & purity >= MIN_PURITY] <- "MUT"
allele[snp_tot > 0 & wt  > mut & purity >= MIN_PURITY] <- "WT"
allele[snp_tot > 0 & is.na(allele)]                    <- "conflicted"

# ---------------- Mecp2 expression: RAW (uncorrected) + as-stored --------------------
# The object's counts layer is SoupX-corrected (it descends from the SoupX merge), so
# pull RAW Mecp2 straight from the CellRanger filtered matrices as well.
cnts <- LayerData(obj, assay = "RNA", layer = "counts")
mecp2_stored <- if ("Mecp2" %in% rownames(cnts)) as.numeric(cnts["Mecp2", ]) else rep(NA_real_, length(cn))

message("reading RAW Mecp2 from CellRanger matrices ...")
raw_list <- lapply(sort(unique(as.character(obj$sample))), function(s) {
  d <- file.path(BASE, s, "outs/filtered_feature_bc_matrix")
  if (!dir.exists(d)) { warning("missing: ", d); return(NULL) }
  m <- Read10X(d); if (is.list(m)) m <- m[["Gene Expression"]]
  if (!"Mecp2" %in% rownames(m)) return(NULL)
  data.table(ckey = paste0(s, "_", colnames(m)), raw = as.numeric(m["Mecp2", ]))
})
rawdt <- rbindlist(raw_list[!sapply(raw_list, is.null)])
mecp2_raw <- rawdt$raw[match(objkey, rawdt$ckey)]
message(sprintf("raw Mecp2 matched for %d of %d cells", sum(!is.na(mecp2_raw)), length(cn)))

# propagate on RAW counts: one SNP read fixes the allele, so the cell's WHOLE raw
# Mecp2 load (incl. UMIs that never spanned the codon) is attributed to that allele
expr_wt  <- ifelse(allele == "WT",  mecp2_raw, ifelse(allele == "MUT", 0, NA))
expr_mut <- ifelse(allele == "MUT", mecp2_raw, ifelse(allele == "WT",  0, NA))

# explicit four-level status so unlabeled-but-expressing cells are visible, not lost
status <- ifelse(!is.na(allele) & allele %in% c("WT","MUT"), paste0("labeled_", allele),
          ifelse(!is.na(allele) & allele == "conflicted",    "conflicted",
          ifelse(!is.na(mecp2_raw) & mecp2_raw > 0,          "expressed_unlabeled",
                                                             "not_expressed")))

amb_vaf <- amb$sample_ambient_vaf[match(as.character(obj$sample), amb$sample)]

add <- data.frame(row.names        = cn,
                  mecp2_wt_umi     = wt,
                  mecp2_mut_umi    = mut,
                  mecp2_snp_umi    = snp_tot,
                  mecp2_vaf_mut    = vaf,
                  mecp2_snp_purity = purity,
                  mecp2_allele     = allele,
                  mecp2_raw_umi    = mecp2_raw,
                  mecp2_stored_umi = mecp2_stored,
                  mecp2_status     = status,
                  mecp2_allele_expr_wt  = expr_wt,
                  mecp2_allele_expr_mut = expr_mut,
                  sample_ambient_vaf    = amb_vaf,
                  stringsAsFactors = FALSE)
obj <- AddMetaData(obj, add)

# convenience grouping for modelling: sex_genotype_allele (allele only where called)
md <- obj[[]]
sx <- if ("animal_sex" %in% names(md)) md$animal_sex else NA
gt <- if ("animal_genotype" %in% names(md)) md$animal_genotype else NA
obj <- AddMetaData(obj, data.frame(row.names = cn,
  deg_group = ifelse(is.na(allele), paste(sx, gt, "unassigned", sep = "_"),
                                     paste(sx, gt, allele, sep = "_")),
  stringsAsFactors = FALSE))

saveRDS(obj, outf)
message("saved ", outf)

# ================= DEG FEASIBILITY REPORT =========================================
md <- obj[[]]
lab <- if ("celltype_hier" %in% names(md)) "celltype_hier" else "celltype_c25"

cat("\n=== allele calls overall ===\n")
print(table(md$mecp2_allele, useNA = "ifany"))
cat("\n=== Mecp2 status (raw, uncorrected) ===\n"); print(table(md$mecp2_status, useNA="ifany"))
cat("\nMecp2+ cells (raw>0):", sum(md$mecp2_raw_umi > 0, na.rm = TRUE),
    "| allele-callable:", sum(md$mecp2_raw_umi > 0 & md$mecp2_allele %in% c("WT","MUT"), na.rm = TRUE),
    "| Mecp2+ but NOT callable:", sum(md$mecp2_status == "expressed_unlabeled", na.rm = TRUE), "\n")
cat("\npropagated Mecp2 UMIs attributed to an allele:",
    sum(md$mecp2_allele_expr_wt, md$mecp2_allele_expr_mut, na.rm = TRUE),
    "(vs", sum(md$mecp2_snp_umi, na.rm=TRUE), "UMIs actually spanning the SNP)\n")

cat("\n=== allele x sex x genotype ===\n")
print(ftable(table(sex = md$animal_sex, genotype = md$animal_genotype,
                   allele = factor(md$mecp2_allele, exclude = NULL))))

# The key within-animal comparison: MUT-active vs WT-active cells in het (mut) FEMALES
fem <- md[which(md$animal_sex == "female" & md$animal_genotype == "mut" &
                md$mecp2_allele %in% c("WT","MUT")), ]
cat("\n=== het-female mosaic comparison: cells per cell type x allele ===\n")
if (nrow(fem)) {
  tb <- table(celltype = fem[[lab]], allele = fem$mecp2_allele)
  tb <- tb[rowSums(tb) > 0, , drop = FALSE]
  print(tb[order(-rowSums(tb)), , drop = FALSE])
  cat("\ncell types with >=20 cells in BOTH alleles (DEG-viable):\n")
  ok <- rownames(tb)[apply(tb, 1, function(r) all(r >= 20))]
  print(if (length(ok)) ok else "none - counts too thin for per-celltype DEG")
  cat("\nper-animal breakdown (needed for pseudobulk):\n")
  print(table(mouse = fem$mouse, allele = fem$mecp2_allele))
} else cat("no het-female allele-called cells found\n")

pdf(fig("20_allele_by_celltype.pdf"), width = 13, height = 7)
if (nrow(fem)) {
  d <- as.data.frame(table(celltype = fem[[lab]], allele = fem$mecp2_allele))
  print(ggplot(d, aes(reorder(celltype, -Freq), Freq, fill = allele)) +
        geom_col(position = "dodge") + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
        labs(x = NULL, y = "cells", title = "het females: allele-called cells per cell type"))
}
dev.off()

pdf(fig("21_ambient_vs_cells.pdf"), width = 8, height = 5)
ds <- md %>% filter(!is.na(mecp2_allele), mecp2_allele %in% c("WT","MUT")) %>%
  group_by(sample) %>%
  summarise(cell_vaf = mean(mecp2_allele == "MUT"),
            ambient  = first(sample_ambient_vaf), n = n())
print(ggplot(ds, aes(ambient, cell_vaf, size = n)) + geom_point(alpha = .7) +
      geom_abline(linetype = 2) + theme_classic() +
      labs(x = "ambient (soup) MUT fraction", y = "per-cell MUT fraction",
           title = "cell calls vs soup baseline (points far above the line = possible contamination)"))
dev.off()

cat("\nNOTE: for DEG, aggregate to pseudobulk per (animal x celltype x allele) before\n",
    "testing. Single-cell tests treat cells from one animal as independent and will\n",
    "inflate significance. Cell types failing the >=20-per-allele check above should be\n",
    "merged upward (use class_broad / neuron_class) rather than tested individually.\n")
