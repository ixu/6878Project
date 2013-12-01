#!/bin/sh
clusters="80 5FPKMAffinity 5FPKMAgg5 5FPKMAgg12";
for cluster in $clusters; do
  bsub -P compbiofolk -q compbio-week -J enrichment_analysis_${cluster} \
       -o out/enrichmentAnalysis_${cluster}.out -n 1 -R "rusage[mem=8]" \
        R CMD BATCH --no-save --no-restore \
       "--args cluster='${cluster}'" \
          enrichmentAnalysis.R Rout/enrichmentAnalysis${cluster}.Rout
done;
