# 05 · Cell type + allele integration

- `attach_alleles.R` — join per-cell WT/MUT calls + raw *Mecp2* + propagated
  allele expression + ambient baseline onto the annotated object. Prints a DEG
  feasibility report.
- `make_celltype.R` — collapse C66->C25, then a count-aware grouping into one
  analysis-ready `cell_type` column (glia kept distinct). Adds `cell_type_c25`.
- `trim_metadata.R` — drop redundant/intermediate metadata columns.
