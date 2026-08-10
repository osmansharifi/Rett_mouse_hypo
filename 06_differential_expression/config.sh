#!/bin/bash
# ==============================================================================
# config.sh
# Cluster/project-specific settings, sourced by job scripts. Keeping this
# separate from the .slurm file means the job script itself doesn't need to
# change if the account, partition, or base paths change later.
# ==============================================================================

# --- SLURM submission defaults (pass at submit time, e.g.:
#     sbatch --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" run_deg_cell_level.slurm ) ---
export SLURM_ACCOUNT="lasallegrp"
export SLURM_PARTITION="high"

# --- project paths ------------------------------------------------------------
export OBJDIR="/quobyte/lasallegrp/Osman/shenyu/02_seurat_objects"
export RESULTS_DIR="/quobyte/lasallegrp/Osman/shenyu/04_results/deg/cell_level_view"
export LOG_DIR="/quobyte/lasallegrp/Osman/shenyu/05_logs"

mkdir -p "$LOG_DIR"
