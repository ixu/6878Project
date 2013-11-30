
import sys
import csv
from copy import deepcopy
from sklearn.cluster import AffinityPropagation
from numpy import matrix
import  numpy as np
import os
import PreProcessAndUtils

def HeirarchicalClusterFast(scoreMatrix,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,OutputDirectory):
    #1. Compute a thresholdDistance from the scoreMatrix
    minAverageDistanceScoreThreshold = getScoreThreshold(scoreMatrix)
    
    #2. start with as many clusters as there are sequences
    clusters = [[i] for i in xrange(len(scoreMatrix[0]))]
    AvgDistanceMatrix = initPairWiseAvgDistanceMatrix(scoreMatrix)

    #3. While the averageClusterDistance is less than threshold or number or clusters = 2
    currentaverageDistanceScore = sys.maxint;
    #while (currentaverageDistanceScore > minAverageDistanceScoreThreshold) & (len(clusters)>1):
    while (len(clusters)>1):
        print len(clusters)
        length = len(clusters)
        #3.1 get two clusters with max pair wise currentaverageDistanceScore
        a,b,averageDistanceScore,AvgDistanceMatrix = getClustersWithMaxPairWiseAverageDistanceScoreFast(clusters,scoreMatrix,AvgDistanceMatrix)

        #3.2 merge the two clusters
        newCluster = clusters[a]  + clusters[b]
        #3.3 if the averageClusterDistance > thresholdDistance
        #   3.4 return the cluster
        if(averageDistanceScore < minAverageDistanceScoreThreshold):
            PreProcessAndUtils.OutputUtils.writeClusters(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,"./Output/"+OutputDirectory+"/Heirarchical_NoOfClusters_"+str(len(clusters)))
         #   return clusters
        #3.5 Add the merged cluster into the clusterlist
        clusters[a] = newCluster
        clusters.pop(b)

    #4. return clusters
    return clusters

VERY_LOW = -999999
def initPairWiseAvgDistanceMatrix(ScoreMatrix):
    ADM = matrix(ScoreMatrix)
    for i in xrange(len(ScoreMatrix[0])):
        ADM[i,i] = VERY_LOW
    return ADM

def getClustersWithMaxPairWiseAverageDistanceScoreFast(clusters,ScoreMatrix,ADM):
    rowLen = len(ADM[0:])
    maxindex = np.argmax(ADM)
    i = maxindex/rowLen
    j = maxindex%rowLen
    averageDistanceScore = ADM[i,j]

    # Correct I and J
    #for c in xrange(rowLen):
    #    ADM[i,c] = (ADM[i,c] + ADM[j,c])/2.0
    #    ADM[c,i] = ADM[i,c]

    newCluster = clusters[i]  + clusters[j]
    noOfClusters = len(clusters)
    for c in xrange(noOfClusters):
            cl1len = len(clusters[c])
            cl2len = len(newCluster)
            sumOfDistances = 0
            for c1 in xrange(cl1len):
                for c2 in xrange(cl2len):
                    sumOfDistances = sumOfDistances + ScoreMatrix[clusters[c][c1]][newCluster[c2]]
            pairwiseAvgDistance = sumOfDistances/(float)(cl1len*cl2len)
            ADM[i,c] = pairwiseAvgDistance
            ADM[c,i] = pairwiseAvgDistance

    ADM[i,i] = VERY_LOW
    ADM = np.delete(ADM,j,axis=0)
    ADM = np.delete(ADM,j,axis=1)
    newrowLen = len(ADM[0:])
    #ADM[j,:] = VERY_LOW
    #ADM[:,j] = VERY_LOW
    
    return i,j,averageDistanceScore,ADM
def cluster(scoreMatrix):
    #1. Compute a thresholdDistance from the scoreMatrix
    minAverageDistanceScoreThreshold = getScoreThreshold(scoreMatrix)

    #2. start with as many clusters as there are sequences
    clusters = [[i] for i in xrange(len(scoreMatrix[0]))]
    
    #3. While the averageClusterDistance is less than threshold or number or clusters = 2
    currentaverageDistanceScore = sys.maxint;
    while (currentaverageDistanceScore > minAverageDistanceScoreThreshold) & (len(clusters)>1):
        #3.1 get two clusters with max pair wise currentaverageDistanceScore
        a,b,averageDistanceScore = getClustersWithMaxPairWiseAverageDistanceScore(clusters,scoreMatrix)

        #3.2 merge the two clusters
        newCluster = clusters[a]  + clusters[b]

        #3.3 if the averageClusterDistance > thresholdDistance
        #   3.4 return the cluster
        if(averageDistanceScore < minAverageDistanceScoreThreshold):
            return clusters
        #3.5 Add the merged cluster into the clusterlist
        clusters[a] = newCluster
        clusters.pop(b)

    #4. return clusters
    return clusters

def getClustersWithMaxPairWiseAverageDistanceScore(clusters,ScoreMatrix):
     #The distance between two clusters C1 and C2 is defined to be the average of the distances between 
    # each pair of sequences in C1 and C2 

    noOfClusters = len(clusters)
    maxAvgDistPair = 0,0
    maxAvgDistance = -sys.maxint
    returnResult =0,0,0
    for i in xrange(noOfClusters):
        for j in xrange(i+1,noOfClusters):
            cl1len = len(clusters[i])
            cl2len = len(clusters[j])
            sumOfDistances = 0
            for c1 in xrange(cl1len):
                for c2 in xrange(cl2len):
                    sumOfDistances = sumOfDistances + ScoreMatrix[clusters[i][c1]][clusters[j][c2]]
            pairwiseAvgDistance = sumOfDistances/(float)(cl1len*cl2len)
            if( pairwiseAvgDistance > maxAvgDistance):
                maxAvgDistance = pairwiseAvgDistance
                maxAvgDistPair = i,j
                returnResult = i,j,pairwiseAvgDistance
    return returnResult

def getScoreThreshold(scoreMatrix):
    maxValue = -sys.maxint 
    minValue = sys.maxint
    for i in xrange(len(scoreMatrix)):
        if(max(scoreMatrix[i])>maxValue): maxValue = max(scoreMatrix[i])
        if(min(scoreMatrix[i])<minValue): minValue = min(scoreMatrix[i])
    return maxValue - ((maxValue-minValue)*0.40) # upper 40 percentile of the score selected for clustering
