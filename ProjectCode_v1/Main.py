import sys
from PreProcessAndUtils import InputUtils
from PreProcessAndUtils import OutputUtils
from Alignment import AlignmentUtils
from Clustering import Clustering

def AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectory):
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
    clusters = Clustering.cluster(scores)
    print " -- No of Clusters generated = ",len(clusters)
    print " -- Writing clusters to ",OutputDirectory," directory"
    OutputUtils.writeClusters(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,"./Output/"+OutputDirectory)
    print " Done. "

def main():
    if len(sys.argv) < 2:
        print "you must call program as:  "
        print "   python Main.py expressionFile.csv eg Main.csv ./Input/mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv "
        return

    
    # 1. Reads an expression time series from a csv file
    print "1. Reading expression time series ..."
    Signals,SignalSeqs,ExprNoToExprIdentifier,timeline = InputUtils.readExpressionTimeSeries(sys.argv[1])
    print " -- Read ",len(Signals)," signals ."

   
    
    S = [
     # R D S
     [3, -1, -3], # R
     [-1, 3, -3], # D
     [-3, -3, 1]  # S
     ]
    gap_pen = 1
    OutputDirectoryName = "AlignedScoreMatrix"
    AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)

    #S = [
    # # R D S
    # [3, -3, -3], # R
    # [-3, 3, -3], # D
    # [-3, -3, 3]  # S
    # ]
    #gap_pen = 1
    #OutputDirectoryName = "StrictMatchScoreMatrix"
    #AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)

    #S = [
    # # R D S
    # [6, -6, -6], # R
    # [-6, 6, -6], # D
    # [-6, -6, 6]  # S
    # ]
    #gap_pen = 3
    #OutputDirectoryName = "VeryStrictMatchScoreMatrix"
    #AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)

    #S = [
    # # R D S
    # [3, 3, -3], # R
    # [3, 3, -3], # D
    # [-3, -3, 3]  # S
    # ]
    #gap_pen = 1
    #OutputDirectoryName = "MatchEqualsAntiMatchScoreMatrix"
    #AlignAndCluster(Signals,SignalSeqs,ExprNoToExprIdentifier,S,gap_pen,timeline,OutputDirectoryName)
    return

main()