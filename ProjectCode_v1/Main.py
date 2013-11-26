import sys
from PreProcessAndUtils import InputUtils
from PreProcessAndUtils import OutputUtils
from Alignment import AlignmentUtils
from Clustering import Clustering
from Clustering import AffinityPropogationClustering
from Annotation import Annotation
import os

UPGMA=0
AFFINITY_PROP=1


def AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectory,ClusterType):
    if not os.path.exists("./Output/"+OutputDirectory):
        os.makedirs("./Output/"+OutputDirectory)
    # 2. Create an alignment matrix between the sequences.
    print "  Aligning signals and creating a score matrix ..."
    scores = AlignmentUtils.getAlignmentScoreMatrixWithLookup(SignalSeqs,S,gap_pen,timeline)
    print " -- Scores generated."
    print " -- Writing to score matrix to AlignmentScores.csv .."
    OutputUtils.writeScoreMatrix(scores,"./Output/"+OutputDirectory+"/AlignmentScoreMatrix.csv")
    print " -- Writing to score and sequences to AlignmentScoresWithSequences.csv .."
    OutputUtils.writeScoresWithSequences(scores,"./Output/"+OutputDirectory+"/AlignmentScoresWithSequences.csv",SignalSeqs)
    
    # 3. Cluster the signals based on the scoring matrix.
    print "  clustering signals based on Score ..."
    if(ClusterType == AFFINITY_PROP):
        print "  Performing Affinity Prop ..."
        clusters = AffinityPropogationClustering.AffinityPropCluster(scores)
    else:
        print "  Performing Heirarchial ..."
        clusters = Clustering.cluster(scores)
    print " -- No of Clusters generated = ",len(clusters)
    print " -- Writing clusters to ",OutputDirectory," directory"
    OutputUtils.writeClusters(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,"./Output/"+OutputDirectory)
    print " Done. "

def TrialClusteringWithVariousScoreMatrix(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName):
    S = [
     # R D S
     [3, -3, -3], # R
     [-3, 3, -3], # D
     [-3, -3, 3]  # S
     ]
    gap_pen = 1
    OutputDirectoryName = "StrictMatchScoreMatrix"
    AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)

    S = [
     # R D S
     [6, -6, -6], # R
     [-6, 6, -6], # D
     [-6, -6, 6]  # S
     ]
    gap_pen = 3
    OutputDirectoryName = "VeryStrictMatchScoreMatrix"
    AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)

    S = [
     # R D S
     [3, 3, -3], # R
     [3, 3, -3], # D
     [-3, -3, 3]  # S
     ]
    gap_pen = 1
    OutputDirectoryName = "MatchEqualsAntiMatchScoreMatrix"
    AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)
    return

def main():
    if len(sys.argv) < 4:
        print "you must call program as:  "
        print "   python Main.py expressionFile.csv outputFile clusterType(0-UPGMA,1-Affinity) eg ./Input/mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv 5FPKMNormalizedAlignedAffinityPropogation 1 or ./Input/mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv 5FPKMNormalizedAligned 0"
        return
    OutputDirectoryName = sys.argv[2]

    ClusterType = int(sys.argv[3])
    # 1. Reads an expression time series from a csv file
    print "1. Reading expression time series ..."
    Signals,SignalSeqs,ExprNoToExprIdentifier,timeline = InputUtils.readExpressionTimeSeries(sys.argv[1])
    print " -- Read ",len(Signals)," signals ."

    #2. Perform Alignment and  Clustering
    print "2. Perform Alignment and  Clustering ..."
    S = [
     # R D S
     [3, -1, -3], # R
     [-1, 3, -3], # D
     [-3, -3, 1]  # S
     ]
    gap_pen = 1
    #OutputDirectoryName = "5FPKMNormalizedAlignedAffinityPropogation"
    #OutputDirectoryName = "5FPKMNormalizedAlignedScoreMatrix"
    AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName,ClusterType)

    #3. Annotation
    print "2. Perform Annotation ..."
    OutputDirectory = "./Output/"+OutputDirectoryName
    Annotation.Annotateclusters(OutputDirectory)
    
    return

main()