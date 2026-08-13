#!/bin/bash
#SBATCH --job-name=deg_cell_level
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=4-04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lotwang@ucdavis.edu
#SBATCH --array=0-5
# Partition/account come from config.sh; pass at submit time, e.g.:
#   sbatch --account="$SLURM_ACCOUNT" --partition="$SLURM_PARTITION" run_deg_cell_level.slurm
# (source config.sh in your shell first so those vars are set before sbatch runs)

set -euo pipefail

# HERE = the directory this script lives in (adjust the path if you move this
# script into a subfolder, mirroring your mentor's 06_differential_expression/ layout)
HERE="/quobyte/lasallegrp/Osman/shenyu/deg"
cd "$HERE"
mkdir -p logs

source "$HERE/config.sh"
source "$HERE/env/env.sh"

SCRIPT="$HERE/cell_level_view/run_deg_cell_level.R"
if [[ ! -f "$SCRIPT" ]]; then
  SCRIPT="$HERE/deg_cell_level/run_deg_cell_level.R"
fi
if [[ ! -f "$SCRIPT" ]]; then
  SCRIPT="$HERE/run_deg_cell_level.R"
fi

# --- combination map (index -> method, cluster_col) --------------------------
METHODS=(wilcox wilcox MAST MAST LimmaVoomCC LimmaVoomCC)
CLUSTER_COLS=(neuron_class cell_type_concise neuron_class cell_type_concise neuron_class cell_type_concise)

METHOD=${METHODS[$SLURM_ARRAY_TASK_ID]}
CLUSTER_COL=${CLUSTER_COLS[$SLURM_ARRAY_TASK_ID]}

echo "Task ${SLURM_ARRAY_TASK_ID}: method=${METHOD} cluster_col=${CLUSTER_COL}"
echo "Using script: ${SCRIPT}"

Rscript "$SCRIPT" \
  --seu "$OBJDIR/seu_final.rds" \
  --file_dir "$RESULTS_DIR" \
  --fig_dir "${FIGURES_DIR:-/quobyte/lasallegrp/Osman/shenyu/03_figures/deg/cell_level_view}" \
  --method "${METHOD}" \
  --cluster_col "${CLUSTER_COL}" \
  --min_pct 0.1 \
  --min_cells 10 \
