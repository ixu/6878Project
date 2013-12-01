#useful for when we have a network with certain rows and columns and now are putting in the same dimensions as other networks, but we specifically want to change only either the rows or columns (i.e. when some networks are directed regulatory ones, others undirected) -- otherwise we could use adjustNetworkDims
mapNamesToSuperset.Matrix <- function(myMatrix, newRownames=NULL, newColnames=NULL,newVal=NA) {

	#don't do anything if nothing to map to
	if (is.null(newRownames) && is.null(newColnames)) {
		return(myMatrix);
	}

	if (is.null(newRownames)) {newRownames=rownames(myMatrix);};
	if (is.null(newColnames)) {newColnames=colnames(myMatrix);};
	
	#names we are mapping to must be superset
	stopifnot(length(setdiff(rownames(myMatrix),newRownames))==0 && length(setdiff(colnames(myMatrix),newColnames)) == 0);


	additionalRowNames=setdiff(newRownames,rownames(myMatrix));
	additionalColNames=setdiff(newColnames,colnames(myMatrix));

	currM=nrow(myMatrix);
	currN=ncol(myMatrix);

	#may not be adding rows if we are only adding columns
	if (length(additionalRowNames)>0) {
		newnet=rbind(myMatrix,matrix(newVal,length(additionalRowNames),ncol(myMatrix)));
		rownames(newnet)[(currM+1):nrow(newnet)]=additionalRowNames;
	}
	

	#may not be adding columns if we are only adding rows
	if (length(additionalColNames)>0) {
		newnet=cbind(newnet,matrix(newVal,nrow(newnet),length(additionalColNames)));
		colnames(newnet)[(currN+1):ncol(newnet)]=additionalColNames;
	};
	newnet=newnet[newRownames,newColnames];

#	newnet=matrix(newVal,length(newRownames),length(newColnames));
#	rownames(newnet)=newRownames;
#	colnames(newnet)=newColnames;

#	newnet[rownames(myMatrix),colnames(myMatrix)]=myMatrix;

	return(newnet);		
}
