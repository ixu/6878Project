User Interface : 
---------------
	index.html
	d3.layout.cloud.js
	clusters.js
	clusters.css

to execute 
	from the project directory python -m SimpleHTTPServer 8888
	and then go to  localhost:8888

Backend Python code. : 
------------------
python Main.py ./Input/mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv 5FPKMNormalizedAligned 0 Here 0 stands for agglomerative\Heirarchical clustering and 1 would be affinity.
Main file Main.py : Performs alignment , followed by clustering and annotation. 
Input : Look at samples in the input directory to get an idea of the expected timecourse matrix. Basically
gene names in the first columns , time\day values in header row and expression values in the matrix.
Output
Look at samples in the Output directory.


Folders
PreprocessAndUtils
	1. Cufflinks_ExpressionAbundance.sh : top hat output to abundane using cufflinks
	2. PerformingQuantileNormalizationInR.txt : steps to perform quantile normalization
	3. FPKM filtering : filtering the o\p above
	4. Log2 ration with a given time point and vector normalization.
	5. InputUtils : helpers to read the various csv files.
	6. Helpers to output results.

Alignment 
	1. seqalign.py : Contains implementation of the needleman wunsh algorithm for sequence alignment
	2. AlignmentUtils.py : Wrapper code to create a distance matrix of alignment scores given a distance matrix.

Clustering
	1. Clustering.py : contains our implementation of agglomerative\heirarchical clustering.
	2. AffinityPropogationClustering.py : Contains our implementation (buggy) and usage of sklearns affinity prop

Annotation
	1.OMIMAnnotation.py : Given the omim.txt , creates a dictionary of gene symbols and corresponding text from it.
	2.Annotation.py : Given a cluster ,
		a. creates an OMIM summary 
		b. converts to ensemble ids using ensemble apis
		c. Fetches david chart and table reports using davids apis

Input 
	1. Contains the various data files that acted as our input. Only 1 of them is used at a time , these merely represent the number
	   of filtering we ended up doing
Output
	1. Output directories contain the clusters as a csv file in addition they contain a CDT and GCT format output to view the clusters
	in other tools like treeview


Python Packages Required
General : numpy
For annotation : httplib2 , json , suds
For Affinity Propogration clustering : sklearn