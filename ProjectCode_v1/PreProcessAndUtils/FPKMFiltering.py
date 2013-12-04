import InputUtils
import OutputUtils
import sys

def main():
    if len(sys.argv) < 3:
        print "you must call program as:  "
        print "   pythonFPKMFiltering.py expressionRawFileName FilteringThreshold OutputFileName eg Input\mRNAExpression_Raw.csv 0.5 Input\mRNAExpression_FilteredBy0.5FPKM.csv"
        return
    rawExpressionFileName = sys.argv[1]
    Signals,SignalSeqs,ExprNoToExprIdentifier,timeline = InputUtils.readExpressionTimeSeries(rawExpressionFileName)
    FilterFPKMValue = float(sys.argv[2])
    mRNAExpressionDataFiltered = {}
    for key in ExprNoToExprIdentifier:
        identifier  = ExprNoToExprIdentifier[key]
        for value in Signals[key]:
                if(value > FilterFPKMValue):
                    mRNAExpressionDataFiltered[identifier] = Signals[key]
                    break
    OutputFileName = sys.argv[3]
    OutputUtils.writeRNAExpressionasCSV(OutputFileName,mRNAExpressionDataFiltered)
    return
main()