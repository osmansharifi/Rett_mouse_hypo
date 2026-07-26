#!/bin/bash
#SBATCH --job-name=raw_mecp2
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Partition/account/qos come from config.sh; pass with --export or set via sbatch flags.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config.sh"
source "$HERE/env/env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

Rscript "$HERE/02_allele_genotyping/raw_mecp2_build.R" \
  "$CELLRANGER_BASE" "$ALLELEDIR" "$SAMPLE_META" "$OBJDIR/mecp2_raw.rds" "$OBJDIR/mecp2_raw_percell.csv"
