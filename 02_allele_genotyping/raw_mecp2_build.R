#!/usr/bin/env Rscript
# RAW Mecp2 build -- no SoupX, no QC filtering, every called cell kept -- PLUS
# within-cell allele propagation: once a SNP-covering read fixes a cell's allele,
# attribute ALL of that cell's raw Mecp2 UMIs to that allele (biallelic Mecp2 is
# biologically rare, so one read labels the whole cell).
# Usage: Rscript 04_raw_mecp2.R <base> <alleledir> <csv> <out_rds> <out_csv>
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})

MIN_PURITY <- 0.75   # dominant-allele share at the SNP needed to label a cell;
                     # below this (e.g. genuine ~50/50) the cell is 'conflicted'.
                     # Set 0.5 for pure majority vote, 1.0 to require unanimity.

args <- commandArgs(trailingOnly = TRUE)
base <- args[1]; alleledir <- args[2]; csvpath <- args[3]
out_rds <- args[4]; out_csv <- args[5]

samples <- basename(Sys.glob(file.path(base, "*_H")))
message("raw build + allele propagation from ", length(samples), " samples")

objs <- lapply(samples, function(s) {
  m <- Read10X(file.path(base, s, "outs/filtered_feature_bc_matrix"))
  if (is.list(m)) m <- m[["Gene Expression"]]
  o <- CreateSeuratObject(counts = m, project = s, min.cells = 0, min.features = 0)
  o$sample <- s
  mecp2 <- if ("Mecp2" %in% rownames(o))
             as.numeric(LayerData(o, layer = "counts")["Mecp2", ]) else 0
  o$mecp2_raw <- mecp2; o$mecp2_expressed <- mecp2 > 0
  o
})
merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = samples, project = "Mecp2_H_raw")
merged <- JoinLayers(merged)

# --- raw allelic UMI counts at the SNP (soup-inclusive) ---
al <- rbindlist(lapply(samples, function(s) {
  f <- file.path(alleledir, paste0(s, ".mecp2_allele.csv"))
  if (file.exists(f)) { d <- tryCatch(fread(f), error = function(e) NULL)
    if (!is.null(d) && "barcode" %in% names(d)) d[barcode != "AMBIENT"] else NULL } else NULL
}), fill = TRUE)
al[, cellkey := paste0(sample, "_", barcode)]
idx <- match(colnames(merged), al$cellkey)
merged$mecp2_wt_umi  <- al$wt_umi[idx]
merged$mecp2_mut_umi <- al$mut_umi[idx]
merged$mecp2_vaf_mut <- suppressWarnings(as.numeric(al$vaf_mut[idx]))
merged$mecp2_allele  <- al$call[idx]

# --- animal metadata ---
meta <- fread(csvpath)
merged$mouse <- gsub("_", "-", sub("_H$", "", merged$sample))
mm <- meta[match(merged$mouse, meta$mouse)]
for (col in c("genotype","sex","timepoint","region","actual_age","sample_name"))
  if (col %in% names(meta)) merged[[paste0("animal_", col)]] <- mm[[col]]

# --- within-cell allele propagation ---
wt <- ifelse(is.na(merged$mecp2_wt_umi),  0, merged$mecp2_wt_umi)
mut<- ifelse(is.na(merged$mecp2_mut_umi), 0, merged$mecp2_mut_umi)
snp_tot <- wt + mut
purity  <- ifelse(snp_tot > 0, pmax(wt, mut) / snp_tot, NA_real_)
dom <- ifelse(mut > wt, "MUT", ifelse(wt > mut, "WT", "tie"))
label <- rep("unlabeled", ncol(merged))                       # no SNP read
label[snp_tot > 0 & dom == "MUT" & purity >= MIN_PURITY] <- "MUT"
label[snp_tot > 0 & dom == "WT"  & purity >= MIN_PURITY] <- "WT"
label[snp_tot > 0 & (dom == "tie" | purity < MIN_PURITY)] <- "conflicted"  # doublet/biallelic
merged$mecp2_snp_umi      <- snp_tot                           # label evidence depth
merged$mecp2_snp_purity   <- purity
merged$mecp2_allele_label <- label
# attribute the cell's FULL raw Mecp2 to its allele (0 to the other; NA if unassignable)
merged$mecp2_wt_expr  <- ifelse(label=="WT",  merged$mecp2_raw, ifelse(label=="MUT", 0, NA))
merged$mecp2_mut_expr <- ifelse(label=="MUT", merged$mecp2_raw, ifelse(label=="WT",  0, NA))

saveRDS(merged, out_rds)
message(sprintf("saved %s: %d cells, %d Mecp2+; labeled WT=%d MUT=%d conflicted=%d",
        out_rds, ncol(merged), sum(merged$mecp2_expressed),
        sum(label=="WT"), sum(label=="MUT"), sum(label=="conflicted")))

# --- per-cell table (every Mecp2+ cell) ---
md <- as.data.table(merged[[]], keep.rownames = "cell")
pos <- md[mecp2_raw > 0, .(cell, sample, mouse, animal_genotype, animal_sex,
            mecp2_raw, mecp2_wt_umi, mecp2_mut_umi, mecp2_snp_umi, mecp2_snp_purity,
            mecp2_allele_label, mecp2_wt_expr, mecp2_mut_expr)]
fwrite(pos, out_csv)
message("wrote per-cell Mecp2+ table: ", out_csv, " (", nrow(pos), " cells)")

# --- per-sample summary: SNP-read pseudobulk + propagated allele expression ---
raw_all <- rbindlist(lapply(samples, function(s) {
  f <- file.path(alleledir, paste0(s, ".mecp2_allele.csv"))
  if (file.exists(f)) tryCatch(fread(f), error = function(e) NULL) else NULL }), fill = TRUE)
snp <- raw_all[, .(snp_wt = sum(wt_umi[barcode!="AMBIENT"]),
                   snp_mut= sum(mut_umi[barcode!="AMBIENT"]),
                   soup_wt= sum(wt_umi[barcode=="AMBIENT"]),
                   soup_mut=sum(mut_umi[barcode=="AMBIENT"])), by = sample]
prop <- md[, .(mecp2_pos_cells = sum(mecp2_raw > 0),
               wt_cells  = sum(mecp2_allele_label=="WT"),
               mut_cells = sum(mecp2_allele_label=="MUT"),
               conflicted= sum(mecp2_allele_label=="conflicted"),
               wt_expr_umi  = sum(mecp2_wt_expr,  na.rm=TRUE),   # propagated expression
               mut_expr_umi = sum(mecp2_mut_expr, na.rm=TRUE)), by = sample]
prop[, frac_mut_expr := mut_expr_umi/(wt_expr_umi+mut_expr_umi)]
psb <- merge(snp, prop, by = "sample", all = TRUE)
persample <- sub("\\.csv$", "_persample.csv", out_csv)
fwrite(psb[order(sample)], persample)
message("wrote per-sample summary: ", persample)
print(psb[order(sample)])
