#!/bin/bash
#SBATCH --job-name=cluster
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Partition/account/qos come from config.sh; pass with --export or set via sbatch flags.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config.sh"
source "$HERE/env/env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

Rscript "$HERE/03_clustering/cluster.R" \
  "$OBJDIR/seu_scaled.rds" "$OBJDIR/seu_clustered.rds" "$OBJDIR/figures/clustering" sample 50 auto
