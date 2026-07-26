# 02 · Allele genotyping

Barcode-aware genotyping of the *Mecp2* e1 start codon and a raw (uncorrected)
*Mecp2* expression object.

- `genotype_mecp2.py` — `samtools mpileup --output-extra CB,UB` at the SNP; per
  `(barcode,UMI)` dedup; WT/MUT tally per barcode + an `AMBIENT` soup row.
  Run once per sample (see `slurm/genotype_array.sbatch`).
- `raw_mecp2_build.R` — raw *Mecp2* counts from CellRanger matrices (no SoupX),
  every cell kept, with within-cell allele propagation.
