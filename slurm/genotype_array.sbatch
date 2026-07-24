#!/bin/bash
#SBATCH --job-name=genotype_array
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Partition/account/qos come from config.sh; pass with --export or set via sbatch flags.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/config.sh"
source "$HERE/env/env.sh"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# genotype every _H sample at the SNP (array: one task per sample)
mapfile -t SAMPLES < <(ls -d "$CELLRANGER_BASE"/*_H | xargs -n1 basename | sort)
s="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
python "$HERE/02_allele_genotyping/genotype_mecp2.py" \
  --bam "$CELLRANGER_BASE/$s/outs/possorted_genome_bam.bam" \
  --sample "$s" --chrom "$SNP_CHR" --pos "$SNP_POS" --wt "$SNP_WT" --mut "$SNP_MUT" \
  --out "$ALLELEDIR/${s}.mecp2_allele.csv"
