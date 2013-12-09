#!/bin/sh
clusters="0 1"

for cluster in $clusters; do
  bsub -P compbiofolk -q compbio-week -J alignment_${cluster} \
       -o out/output_alignment_${cluster} -n 8 \
       -R "rusage[mem=16]" -R "span[hosts=1]" python Main.py ./Input/GSE675.csv GSE675_${cluster} ${cluster}
done;
