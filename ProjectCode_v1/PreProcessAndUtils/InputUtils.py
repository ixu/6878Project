import csv
#import seqalign
import timeit
thresholdPercent = 0.10

#def read_file(filename):
#	f = open(filename)
#	signals = getSignals(f.readlines()[:10])
#	f.close()
#	scores = [[0] * len(signals) for i in range(len(signals))]
#	for i in range(len(signals)):
#		for j in range(i + 1, len(signals)):
#			scores[i][j] = seqalign.seqalignDP(toString(signals[i]),toString(signals[j]))
#	return scores

def toString(signal):
    floatSignal = toFloat(signal)
    thresholdValue = (max(floatSignal)-min(floatSignal))*thresholdPercent
    string = [getChar((floatSignal[i + 1] - floatSignal[i]),thresholdValue)  for i in xrange(len(signal)-1)]
    return ['S'] + string


def toFloat(signal):
    return [float(value) for value in signal]

def getChar(delta,threshold):
	if delta >= threshold:
		return 'R' # rise
	elif delta <= -1 * threshold:
		return 'D' # drop
	return 'S' # steady

def getSignals(lines):
	return [line.strip().split(',') for line in lines]

# Reads an Expression sequence file as a csv matrix
#    First row specifies timepoints in days for now ( 1, 4, 5....etc)
#    First column specifies Identifiers : gene , protein etc
# example
# 0,1,4,21,116,185,186,255,289,290,292,294,297,301,307,311,322,329,369,380,400
# POLR2E:NM_002695,14.9422196803,16.3586229044,0.0,21.4938814525,18.5816896932,18.136534988,17.8763498023,24.2978732055,22.1762345548,22.5305949443,30.5495531776,30.9540226599,25.1319306705,17.2693724361,36.0807441314,48.9203151113,40.424099104,39.4981743925,38.9747630512,30.8782176975
# MIR3179-1:NR_036140_2,0.0,0.0,0.0,0.0,19.1162126737,21.8722911461,0.0,3.0121016638,0.0,15.3486787946,0.0,0.0,0.0,0.0,0.0,0.0,15.7971888665,5.4176044866,0.0,0.0
# NPW:NM_001099456,10.6080365108,7.9065508811,0.0,3.9936847143,3.4689824865,3.3243029263,9.1902126449,3.8102066883,3.4207679421,4.199802789,3.7636620944,3.0873540925,4.5812543872,1.4491448367,4.9263915936,4.9433420343,3.6486359278,4.3001887484,3.043022313,4.8214587228
#
# Returns
#       TimeLine []
#       Signals [][] 14.9422196803,16.3586229044,0.0,21.4938814525..
#       SignalSeq [][] 'R','S'...
#       GeneIdToSignalId Dictionary {}
def readExpressionTimeSeries(expressionTimeSeriesFileName):
    f = open(expressionTimeSeriesFileName)
    lines = getSignals(f.readlines())

     # read Timeline
    timeline = [ float(timeValue) for timeValue in lines[0][1:] ]

    # read signal
    ExprNoToExprIdentifier = {}
    Signals= []
    SignalSeqs = []
    signalId = 0
    for line in lines[1:]:
        if(len(line) < 2): continue
        ExprNoToExprIdentifier[signalId] = line[0]
        signalId += 1
        SignalSeqs.append(toString(line[1:]))
        Signals.append(toFloat(line[1:]))
    
    f.close()

    return Signals,SignalSeqs,ExprNoToExprIdentifier,timeline

# reads every alternate line , usefull for affy type arrays 22000 long
def readExpressionTimeSeriesTrimmed(expressionTimeSeriesFileName):
    f = open(expressionTimeSeriesFileName)
    lines = getSignals(f.readlines())
    noOflines = len(lines)
    
     # read Timeline
    timeline = [ float(timeValue) for timeValue in lines[0][1:] ]

    # read signal
    ExprNoToExprIdentifier = {}
    Signals= []
    SignalSeqs = []
    signalId = 0
    alternate = False
    for line in lines[1:]:
        if(len(line) < 2): continue
        if(alternate == True): 
            alternate = False
            continue
        else: alternate = True
        ExprNoToExprIdentifier[signalId] = line[0]
        signalId += 1
        SignalSeqs.append(toString(line[1:]))
        Signals.append(toFloat(line[1:]))
    
    f.close()

    return Signals,SignalSeqs,ExprNoToExprIdentifier,timeline

def readScore(ScoreFile):
    f = open(ScoreFile,'rb')
    reader = csv.reader(f)
    scores = []
    for row in reader:
        scoreValues = [ float(i) for i in row]
        scores.append(scoreValues)
    scores = scores[1:]
    return scores