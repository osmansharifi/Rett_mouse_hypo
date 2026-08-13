#!/bin/bash
# ==============================================================================
# env/env.sh
# ==============================================================================

module load conda/base/latest
module load R/4.4.2

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

