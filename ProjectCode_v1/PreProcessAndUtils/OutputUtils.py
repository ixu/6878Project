import csv
import os
#### writeScoreMatrix:
#### 	INPUTS:	ScoreMatrix - matrix of scores
####            Output file name
def writeScoreMatrix(ScoreMatrix,outputFileName):
    scorefile = open(outputFileName,"wb")
    writer = csv.writer(scorefile)
    header = [ i for i in xrange(len(ScoreMatrix[0]))]
    writer.writerow(header)
    for i in xrange(len(ScoreMatrix)):
        writer.writerow(ScoreMatrix[i])
    scorefile.close()
    return


#### writeScoreMatrix with sequences:
def writeScoresWithSequences(ScoreMatrix,outputFileName,SignalSeqs):
    f = open(outputFileName,"wb")
    f.write("score,i,j,seqi,seqj\n")
    for i in range(len(SignalSeqs)):
        for j in range(i + 1, len(SignalSeqs)):
            f.write(str(ScoreMatrix[i][j])+","+str(i)+","+str(j)+",'"+''.join(SignalSeqs[i])+"','"+''.join(SignalSeqs[j])+"'\n")
    f.close()
    return

# write clusters to files
def writeClusters(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory):
 for i in xrange(len(clusters)):
        writeCluster(i,clusters[i],Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory)
 WriteClustersAsCDT(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory)
 WriteClustersAsGCT(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory)
 ConvertCDTToSetGeneMatrixFile(directory)
 return

def writeCluster(clusterNumber,cluster,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    scorefile = open(directory+"/Cluster_"+str(clusterNumber)+".csv","wb")
    writer = csv.writer(scorefile)
    header =[0]+timeline
    writer.writerow(header)
    for j in xrange(len(cluster)):
        row = [ExprNoToExprIdentifier[cluster[j]]]+Signals[cluster[j]]
        writer.writerow(row)
    scorefile.close()
    return

def WriteClustersAsCDT(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    scorefile = open(directory+"/Cluster_ALL.CDT","wb")
    writer = csv.writer(scorefile,dialect='excel-tab')
    header =['day','ClusterNum','GWEIGHT']
    header =header+timeline
    writer.writerow(header)
    secondRow = ['EWEIGHT','','']
    for i in xrange(len(timeline)):
        secondRow.append('1.0')
    writer.writerow(secondRow)
    for i in xrange(len(clusters)):
        for j in xrange(len(clusters[i])):
            row = [ExprNoToExprIdentifier[clusters[i][j]],str(i),1.0]+Signals[clusters[i][j]]
            writer.writerow(row)
    scorefile.close()
    return

def WriteClustersAsGCT(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    scorefile = open(directory+"/Cluster_ALL.gct","wb")
    writer = csv.writer(scorefile,dialect='excel-tab')
    firstRow = ['#1.2']
    writer.writerow(firstRow)
    secondRow = [len(Signals),len(timeline)]
    writer.writerow(secondRow)
    header =['NAME','DescriptionClusterNo']
    header =header+timeline
    writer.writerow(header)
    for i in xrange(len(clusters)):
        for j in xrange(len(clusters[i])):
            row = [ExprNoToExprIdentifier[clusters[i][j]],str(i)]+Signals[clusters[i][j]]
            writer.writerow(row)
    scorefile.close()
    return

def ConvertCDTToSetGeneMatrixFile(outputDir):
    print " converting :",outputDir+"/Cluster_ALL.CDT"
    cdtFile = open(outputDir+"/Cluster_ALL.CDT","rb")
    reader = csv.reader(cdtFile,dialect='excel-tab')
    GeneSymbolToClusterIdDict = {}
    maxClusterNo = 0
    for row in reader:
        if(row[0] == "EWEIGHT" or row[0] == "day"): continue
        geneSymbol = row[0].split(":")[0]
        clusterNumber = int(row[1])
        GeneSymbolToClusterIdDict[geneSymbol] = clusterNumber
        if(clusterNumber > maxClusterNo): maxClusterNo = clusterNumber
    cdtFile.close()
    setGeneMatrixFile = open(outputDir+"/SetGeneMatrix.csv","wb")
    writer = csv.writer(setGeneMatrixFile,dialect='excel-tab')
    holder = [str(0) for i in xrange(maxClusterNo+2)]
    for geneSymbol in GeneSymbolToClusterIdDict.keys():
        holder[0] = geneSymbol
        holder[GeneSymbolToClusterIdDict[geneSymbol]+1]=str(1)
        writer.writerow(holder)
        holder[0] = ""
        holder[GeneSymbolToClusterIdDict[geneSymbol]+1]=str(0)
    setGeneMatrixFile.close()
    return
