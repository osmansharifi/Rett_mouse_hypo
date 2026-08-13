#!/bin/bash
# Submit sample-level and cell-level DEG Slurm arrays.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE"

source "$ROOT/config.sh"

mkdir -p "$ROOT/logs" "$HERE/logs"

SAMPLE_LEVEL_DIR="$HERE/sample_level_view"
CELL_LEVEL_DIR="$HERE/cell_level_view"
if [[ ! -d "$SAMPLE_LEVEL_DIR" ]]; then
  SAMPLE_LEVEL_DIR="$HERE/deg_sample_level"
fi
if [[ ! -d "$CELL_LEVEL_DIR" ]]; then
  CELL_LEVEL_DIR="$HERE/deg_cell_level"
fi

SBATCH_ARGS=(--partition="$SLURM_PARTITION")
if [[ -n "${SLURM_ACCOUNT:-}" ]]; then
  SBATCH_ARGS+=(--account="$SLURM_ACCOUNT")
fi
if [[ -n "${SLURM_QOS:-}" ]]; then
  SBATCH_ARGS+=(--qos="$SLURM_QOS")
fi

echo "Submitting sample-level DEG array ..."
sbatch "${SBATCH_ARGS[@]}" "$SAMPLE_LEVEL_DIR/run_deg_sample_level.slurm"

echo "Submitting cell-level DEG array ..."
sbatch "${SBATCH_ARGS[@]}" "$CELL_LEVEL_DIR/run_deg_cell_level.slurm"

echo "Submitted DEG jobs."
