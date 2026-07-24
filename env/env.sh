#!/bin/bash
# Activate the self-contained pipeline environment (R + Seurat + SoupX + harmony +
# DESeq2, samtools, python). Edit CONDA_SH / PIPELINE_ENV for your system.
CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
[ -f "$CONDA_SH" ] || { echo "FATAL: conda.sh not found at $CONDA_SH" >&2; exit 1; }
source "$CONDA_SH"
conda activate "${PIPELINE_ENV:?set PIPELINE_ENV (see config.sh)}"
command -v Rscript >/dev/null || { echo "FATAL: R not found in $PIPELINE_ENV" >&2; exit 1; }
command -v samtools >/dev/null || { echo "FATAL: samtools not found in $PIPELINE_ENV" >&2; exit 1; }
