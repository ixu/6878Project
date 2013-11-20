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
 return

def writeCluster(clusterNumber,cluster,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    scorefile = open(directory+"/Cluster_"+str(clusterNumber)+".csv","wb")
    writer = csv.writer(scorefile)
    header =[0]+timeline
    writer.writerow(header)
    for clusterNumber in xrange(len(cluster)):
        row = [ExprNoToExprIdentifier[cluster[clusterNumber]]]+Signals[cluster[clusterNumber]]
        writer.writerow(row)
    scorefile.close()
    return