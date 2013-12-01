args=(commandArgs(TRUE));
if (length(args)==0) {
  stop("No arguments supplied.");
} else {
  eval(parse(text=args[[1]])); # parse first argument: cluster
}
load('human-GO.RData');
file <- paste('input/SetGeneMatrix',cluster,'.csv',sep="");
setGeneMatrix <- read.table(file, header = FALSE, sep = "\t",row.names = 1);
source('fisherTestGeneSets.R');
geneAnnotationMatrix <-as.matrix(GOdata);
setGeneMatrix <- t(as.matrix(setGeneMatrix));
#backgroundIDset <- colnames(setGeneMatrix);
#minSetSize <- 0;
#maxSetSize <- Inf;
results <- fisherTestGeneSets(setGeneMatrix, geneAnnotationMatrix);
save(results, file=paste('output/results',cluster,'.RData',sep=""))
write.csv(t(results$pValueMatrix),file=paste('output/results',cluster,'.csv',sep="")) 

