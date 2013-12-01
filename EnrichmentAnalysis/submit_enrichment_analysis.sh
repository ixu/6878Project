#!/bin/sh
clusters="Affinity 10 20 30";
for cluster in $clusters; do
  bsub -P compbiofolk -q compbio-week -J enrichment_analysis \
       -o out/enrichmentAnalysis_${cluster}.out -n 1\
       -R "rusage[mem=8]" R CMD BATCH --no-save --no-restore \
       "--args cluster='${cluster}'" \
          enrichmentAnalysis.R Rout/enrichmentAnalysis.Rout
done;
