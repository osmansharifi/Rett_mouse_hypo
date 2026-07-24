#!/bin/bash
# End-to-end submission with SLURM dependencies. Edit config.sh first.
# Assumes 01_matrix_preprocessing/ has already produced seu_scaled.rds.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/config.sh"
SB="sbatch --parsable --account=$SLURM_ACCOUNT --partition=$SLURM_PARTITION --qos=$SLURM_QOS"
N=$(ls -d "$CELLRANGER_BASE"/*_H | wc -l)

# --- allele track (independent of clustering) ---
g=$($SB --array=0-$((N-1)) slurm/genotype_array.sbatch)
r=$($SB --dependency=afterok:$g slurm/raw_mecp2.sbatch)

# --- expression track ---
c=$($SB slurm/cluster.sbatch)
a25=$($SB --dependency=afterok:$c slurm/annotate.sbatch C25_named)
a66=$($SB --dependency=afterok:$c slurm/annotate.sbatch C66_named)
h=$($SB --dependency=afterok:$a25:$a66 slurm/hier_labels.sbatch)

# --- convergence + outputs ---
at=$($SB --dependency=afterok:$h:$g slurm/attach_alleles.sbatch)
ct=$($SB --dependency=afterok:$at slurm/make_celltype.sbatch)
$SB --dependency=afterok:$ct slurm/pseudobulk_deg.sbatch

echo "submitted: genotype=$g raw=$r cluster=$c annotate=$a25,$a66 hier=$h attach=$at celltype=$ct"
