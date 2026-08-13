#!/bin/bash
#SBATCH --job-name=deg_sample_level
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=4-04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lotwang@ucdavis.edu
#SBATCH --array=0-5

set -euo pipefail

# Submit from the repository root or adjust HERE to the deployed DEG directory.
HERE="/quobyte/lasallegrp/Osman/shenyu/deg"
cd "$HERE"
mkdir -p logs

source "$HERE/config.sh"
source "$HERE/env/env.sh"

SCRIPT="$HERE/sample_level_view/run_deg_sample_level.R"
if [[ ! -f "$SCRIPT" ]]; then
  SCRIPT="$HERE/deg_sample_level/run_deg_sample_level.R"
fi
if [[ ! -f "$SCRIPT" ]]; then
  SCRIPT="$HERE/run_deg_sample_level.R"
fi

METHODS=(limma_voomcc limma_voomcc edger_zinger edger_zinger deseq2 deseq2)
CLUSTER_COLS=(neuron_class cell_type_concise neuron_class cell_type_concise neuron_class cell_type_concise)

METHOD=${METHODS[$SLURM_ARRAY_TASK_ID]}
CLUSTER_COL=${CLUSTER_COLS[$SLURM_ARRAY_TASK_ID]}

echo "Task ${SLURM_ARRAY_TASK_ID}: method=${METHOD} cluster_col=${CLUSTER_COL}"
echo "Using script: ${SCRIPT}"

Rscript "$SCRIPT" \
  --method "${METHOD}" \
  --cluster_col "${CLUSTER_COL}"
