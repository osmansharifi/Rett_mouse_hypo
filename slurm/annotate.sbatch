#!/bin/bash
#SBATCH --job-name=annotate
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Partition/account/qos come from config.sh; pass with --export or set via sbatch flags.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config.sh"
source "$HERE/env/env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"

# $1 = C25_named or C66_named ; writes seu_annotated_<lvl>.rds
LVL="${1:-C25_named}"
Rscript "$HERE/04_annotation/annotate_hypomap.R" \
  "$OBJDIR/seu_clustered.rds" "$HYPOMAP" \
  "$OBJDIR/seu_annotated_${LVL}.rds" "$OBJDIR/figures/annotation_${LVL}" "$LVL" 30
