import sys
from PreProcessAndUtils import InputUtils
from PreProcessAndUtils import OutputUtils
from Alignment import AlignmentUtils
from Clustering import Clustering

def main():
    if len(sys.argv) < 2:
        print "you must call program as:  "
        print "   python Main.py expressionFile.csv eg Main.csv ./Input/mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv "
        return

    
    # 1. Reads an expression time series from a csv file
    print "1. Reading expression time series ..."
    Signals,SignalSeqs,ExprNoToExprIdentifier,timeline = InputUtils.readExpressionTimeSeries(sys.argv[1])
    print " -- Read ",len(Signals)," signals ."

    # 2. Create an alignment matrix between the sequences.
    print "2. Aligning signals and creating a score matrix ..."
    S = [
     # R D S
     [3, -1, -3], # R
     [-1, 3, -3], # D
     [-3, -3, 1]  # S
     ]
    gap_pen = 1
    scores = AlignmentUtils.getAlignmentScoreMatrixWithLookup(SignalSeqs,S,gap_pen,timeline)
    print " -- Scores generated."
    print " -- Writing to score matrix to AlignmentScores.csv .."
    OutputUtils.writeScoreMatrix(scores,"./Output/AlignmentScoreMatrix.csv")
    print " -- Writing to score and sequences to AlignmentScoresWithSequences.csv .."
    OutputUtils.writeScoresWithSequences(scores,"./Output/AlignmentScoresWithSequences.csv",SignalSeqs)
    
    # 2. Cluster the signals based on the scoring matrix.
    print "3.  clustering signals based on Score ..."
    clusters = Clustering.cluster(scores)
    print " -- No of Clusters generated = ",len(clusters)
    print " -- Writing clusters to ./Output directory"
    OutputUtils.writeClusters(clusters,Signals,SignalSeqs,ExprNoToExprIdentifier,timeline,"./Output")
    print " Done. "
    return

main()