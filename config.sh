#!/bin/bash
# Central configuration -- edit these paths for your environment, everything sources this.
# ---------------------------------------------------------------------------------------
# Conda env holding R (Seurat/SoupX/harmony/DESeq2), samtools, python (see env/).
export PIPELINE_ENV="${PIPELINE_ENV:-/path/to/envs/rett_hypo}"

# CellRanger output root: one <SAMPLE>/outs/ dir per sample (samples end in _H).
export CELLRANGER_BASE="${CELLRANGER_BASE:-/path/to/cellranger}"

# Working dirs.
export SCRIPTS="${SCRIPTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export OBJDIR="${OBJDIR:-/path/to/seurat_objects}"      # .rds objects live here
export ALLELEDIR="${ALLELEDIR:-$OBJDIR/allele}"          # per-sample allele CSVs
export LOGDIR="${LOGDIR:-$SCRIPTS/logs}"

# Reference files.
export SAMPLE_META="${SAMPLE_META:-/path/to/hypo_sample_metadata.csv}"
export HYPOMAP="${HYPOMAP:-$OBJDIR/hypoMap.rds}"
export C25C66_MAP="${C25C66_MAP:-$OBJDIR/hypomap_c25_c66_map.csv}"

# Mecp2-e1 start codon (mm10 / GRCm38), Mecp2 is minus-strand so WT=T, MUT=A on + strand.
export SNP_CHR="${SNP_CHR:-X}"
export SNP_POS="${SNP_POS:-74085586}"
export SNP_WT="${SNP_WT:-T}"
export SNP_MUT="${SNP_MUT:-A}"

# SLURM defaults.
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-your_account}"
export SLURM_PARTITION="${SLURM_PARTITION:-high}"
export SLURM_QOS="${SLURM_QOS:-your_qos}"

mkdir -p "$LOGDIR" "$ALLELEDIR"
