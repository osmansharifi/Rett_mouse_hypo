#!/bin/bash
#SBATCH --job-name=hier_labels
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Partition/account/qos come from config.sh; pass with --export or set via sbatch flags.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config.sh"
source "$HERE/env/env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

Rscript "$HERE/04_annotation/hierarchical_labels.R" \
  "$OBJDIR/seu_annotated_C25_named.rds" "$OBJDIR/seu_annotated_C66_named.rds" \
  "$C25C66_MAP" "$OBJDIR/seu_annotated_hier.rds" "$OBJDIR/figures/hier_labels"
