# 04 · Annotation

Cell-type labels transferred from HypoMap.

- `make_c25_c66_map.R` — derive the C25<->C66 hierarchy table (run once).
- `inspect_objects.R` — report layers/labels of query + reference before transfer.
- `annotate_hypomap.R` — `FindTransferAnchors`/`TransferData` at one level
  (`C25_named` or `C66_named`); validation figures incl. marker dotplot.
- `hierarchical_labels.R` — keep C66 only where confident **and** consistent with
  the C25 parent, else fall back to C25. Also emits coarse tiers.
