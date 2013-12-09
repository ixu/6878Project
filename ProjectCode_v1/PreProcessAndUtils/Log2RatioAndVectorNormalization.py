import InputUtils
import OutputUtils
import sys
import math

def main():
    if len(sys.argv) < 3:
        print "you must call program as:  "
        print "   pythonFPKMFiltering.py expressionRawFileName normColumn OutputFileName eg .\input\mRNAExpression_FilteredBy5FPKM_QuantileNormalized.csv  7 .\input\mRNAExpression_FilteredBy5FPKM_QuantileNormalizedLog2Filtered.csv "
        return
    rawExpressionFileName = sys.argv[1]
    Signals,SignalSeqs,ExprNoToExprIdentifier,timeline = InputUtils.readExpressionTimeSeries(rawExpressionFileName)
    normColumn = int(sys.argv[2])-1
    mRNAExpressionDatacorrectedAndNormalized = {}
    for key in ExprNoToExprIdentifier:
        identifier  = ExprNoToExprIdentifier[key]
        normValue = Signals[key][normColumn]
        mRNAExpressionDatacorrectedAndNormalized[identifier] = Signals[key]
        for i in xrange(len(Signals[key])):
            newValue = math.log((Signals[key][i]+0.000000000001)/(normValue+0.000000000001),2)
            mRNAExpressionDatacorrectedAndNormalized[identifier][i] = newValue
    
    # vector normalization.....
    for j in xrange(len(mRNAExpressionDatacorrectedAndNormalized[ExprNoToExprIdentifier[0]])):
        sumSq = 0.0
        for key in ExprNoToExprIdentifier:
            i  = ExprNoToExprIdentifier[key]
            sumSq = sumSq + mRNAExpressionDatacorrectedAndNormalized[i][j]*mRNAExpressionDatacorrectedAndNormalized[i][j]
        sqSumSq = math.sqrt(sumSq)
        if(sqSumSq == 0):
            sqSumSq = sys.float_info.epsilon
        for key in ExprNoToExprIdentifier:
            i  = ExprNoToExprIdentifier[key]
            mRNAExpressionDatacorrectedAndNormalized[i][j] = mRNAExpressionDatacorrectedAndNormalized[i][j]/sqSumSq


    OutputFileName = sys.argv[3]
    OutputUtils.writeRNAExpressionasCSV(OutputFileName,mRNAExpressionDatacorrectedAndNormalized)
    return
main()