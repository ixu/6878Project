#!/bin/sh
classes="const";
#cell_types="GM12878 HELA K562 H1";
cell_types="GM12878";
#types="vanilla matchnodes matchdists matchboth";
#nums="100 500 1000 2000";
types="matchdists";             # Only 1 FG/BG match-type, to simplify to get Iris started.
nums="1000";                    # Only one value for the number of interactions to be selected.
sets="1 2 3 4 5"

clusters = "0 1"

for cluster in $clusters; do
  bsub -P compbiofolk -q compbio-week -J alignment_${cluster} \
       -o out/output_alignment_${cluster} -n 1 \
       -R "rusage[mem=16]" python Main.py ./Input/mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv 5FPKMNormalizedAligned_${cluster} ${cluster}
done;
