#setwd('~/Dropbox/2013-2014/6.878/6878Project/EnrichmentAnalysis')
args=(commandArgs(TRUE));
if (length(args)==0) {
  stop("No arguments supplied.");
} else {
  eval(parse(text=args[[1]])); # parse first argument: cluster
}
load('human-GO.RData');
file <- paste('input/SetGeneMatrix',cluster,'.csv',sep="");
print(file);
setGeneMatrix <- read.table(file, header = FALSE, sep = "\t",row.names = 1);
source('fisherTestGeneSets.R');
geneAnnotationMatrix <-as.matrix(GOdata);
setGeneMatrix <- t(as.matrix(setGeneMatrix));
print(dim(setGeneMatrix))
backgroundIDset <- colnames(setGeneMatrix);
minSetSize <- 0;
maxSetSize <- Inf;
results <- fisherTestGeneSets(setGeneMatrix, geneAnnotationMatrix);
print(dim(results$pValueMatrix))
save(results, file=paste('output/results',cluster,'.RData',sep=""))
write.table(t(results$pValueMatrix),file=paste('output/results',cluster,'.csv',sep=","))