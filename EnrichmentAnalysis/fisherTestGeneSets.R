library(Matrix);

fisherTestGeneSets <- function(setGeneMatrix, geneAnnotationMatrix, backgroundIDset=colnames(setGeneMatrix), minSetSize=0, maxSetSize=Inf) {
#setGeneMatrix -- matrix of gene sets (rows) and gene memberships (columns)
#backgroundIDset -- set of genes from which the columns of setGeneMatrix is drawn
#geneAnnotationMatrix is a GxA matrix, where G is the # of genes, A is the # of annotations, geneAnnotationMatrix[gg,aa]==1 if gene gg belongs to annotation group aa

#note that first, we will intersect the columns of setGeneMatrix with the backgroundIDset, because we are measuring enrichment of drawing the setGeneMatrix from backgroundIDset, so doesn't make sense if the the columns of setGeneMatrix is not a subset.
#note however that the backgroundIDset may not be a strict subset of geneAnnotationMatrix -- we will add extra rows to geneAnnotationMatrix to make it a subset -- just because there is no annotation, does not mean we should consider the gene set smaller than it actually is.


source('./mapNamesToSuperset.Matrix.R');

#first, only consider columns of setGeneMatrix that are in the background.  Note typically we use the whole (i.e. GO) geneAnnotationMatrix as the background gene set.  ideally, the background set (and the geneAnnotationMatrix) would cover the entire genome.  but i think almost all genes in geneAnnotationMatrix had at least one annotation (in each of BP, CC, and BF), so it's a biased set.  so it's ok that we're reducing our 'foreground set' to those that are actually also in geneAnnotationMatrix (though it should also technically include those genes that have no annotation, but since the background set all have at least 1 annotation, it is ok).
commongenes = intersect(colnames(setGeneMatrix), backgroundIDset);
#return if there are no genes left
if (length(commongenes)==0) {
	return;
}
setGeneMatrix=setGeneMatrix[,commongenes,drop=F];

if (length(setdiff(backgroundIDset,colnames(setGeneMatrix)))>0) {
	setGeneMatrix=mapNamesToSuperset.Matrix(setGeneMatrix,newColnames=backgroundIDset);
}

geneAnnotationMatrix=geneAnnotationMatrix[intersect(backgroundIDset,rownames(geneAnnotationMatrix)),,drop=F];

#next, pad geneAnnotationMatrix with any genes in the background set that are not in the annotation matrix
missingGenes=setdiff(backgroundIDset, rownames(geneAnnotationMatrix));
if (length(missingGenes)>0) {
	geneAnnotationMatrix=mapNamesToSuperset.Matrix(geneAnnotationMatrix,newRownames=backgroundIDset,newVal=0);
}

#restrict enrichment testing to those gene groups with at least minSetSize members, but no bigger than maxSetSize 
goodix = which(colSums(geneAnnotationMatrix)>=minSetSize & colSums(geneAnnotationMatrix) <= maxSetSize);
geneAnnotationMatrix = geneAnnotationMatrix[,goodix,drop=F];


#store pvalues
pValueMatrix=matrix(1,nrow(setGeneMatrix),ncol(geneAnnotationMatrix));
rownames(pValueMatrix)=rownames(setGeneMatrix);
colnames(pValueMatrix)=colnames(geneAnnotationMatrix);

#this is for the Fisher's test

QQ = setGeneMatrix %*% geneAnnotationMatrix;
MM = colSums(geneAnnotationMatrix);
NN = nrow(geneAnnotationMatrix)-MM;
KK = rowSums(setGeneMatrix);

for (ii in 1:nrow(pValueMatrix)) {
	for (jj in which(QQ[ii,]>0)) {
		pValueMatrix[ii,jj] = phyper(q=(QQ[ii,jj]-1),m=MM[jj],n=NN[jj],k=KK[ii],lower.tail=FALSE);
	}
}

geneSetSizes=cbind(rowSums(setGeneMatrix));
annotationSizes=rbind(colSums(geneAnnotationMatrix));
rownames(geneSetSizes)=rownames(setGeneMatrix);
colnames(annotationSizes)=colnames(geneAnnotationMatrix)
returnlist=list(pValueMatrix=pValueMatrix,geneSetSizes=geneSetSizes,annotationSizes=annotationSizes);
returnlist
}
