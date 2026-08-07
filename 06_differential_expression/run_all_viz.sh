#!/bin/bash
# ==============================================================================
# run_all_viz.sh
# Run DEG visualizations for all methods and cluster granularities
# ==============================================================================

# Load system R environment (the most effective and conflict-free method confirmed earlier)
module load R/4.4.2

# Define all combinations to iterate over
METHODS=("wilcox" "MAST")
CLUSTERS=("neuron_class" "cell_type_concise")

for method in "${METHODS[@]}"; do
  for cluster in "${CLUSTERS[@]}"; do
    echo ""
    echo "=========================================================="
    echo "▶ Processing: Method = ${method}, Cluster Level = ${cluster}"
    echo "=========================================================="
    
    Rscript run_viz.R --method "$method" --cluster_col "$cluster"
    
  done
done

echo ""
echo "🎉 All visualization tasks are complete! Figures are saved in: /quobyte/lasallegrp/Osman/shenyu/03_figures/deg/cell_level_view/"
