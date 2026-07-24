#!/bin/bash
#SBATCH --job-name=attach_alleles
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Partition/account/qos come from config.sh; pass with --export or set via sbatch flags.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config.sh"
source "$HERE/env/env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

Rscript "$HERE/05_celltype_allele/attach_alleles.R" \
  "$OBJDIR/seu_annotated_hier.rds" "$ALLELEDIR" \
  "$OBJDIR/seu_annotated_allele.rds" "$OBJDIR/figures/allele" "$CELLRANGER_BASE" 0.75
