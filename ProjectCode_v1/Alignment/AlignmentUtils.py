import seqalign


def getAlignmentScoreMatrix(signalSeqs,S,gap_pen,timeline):
    scores = [[0] * len(signalSeqs) for i in range(len(signalSeqs))]
    for i in range(len(signalSeqs)):
        for j in range(i + 1, len(signalSeqs)):
            scores[i][j] = seqalign.seqalignDP(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline)
            scores[j][i] = scores[i][j]
    return scores

def getAlignmentScoreMatrixWithLookup(signalSeqs,S,gap_pen,timeline):
    ScoreDictionary = {}
    lookupCount = 0
    scores = [[0] * len(signalSeqs) for i in range(len(signalSeqs))]
    for i in range(len(signalSeqs)):
        for j in range(i + 1, len(signalSeqs)):
            doesScoreExist,lookupScore,lookupCount = lookUpScore(signalSeqs[i],signalSeqs[j],ScoreDictionary,lookupCount)
            if(doesScoreExist == True): scores[i][j] = lookupScore
            else:
                scores[i][j] = seqalign.seqalignDP(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline)
                addToLookupScore(signalSeqs[i],signalSeqs[j],scores[i][j],ScoreDictionary)
            scores[j][i] = scores[i][j]
    print " -- Successful lookups :",lookupCount
    return scores



def lookUpScore(signalSeq1,signalSeq2,ScoreDictionary,lookupCount):
    key = ''.join(signalSeq1)+''.join(signalSeq2)
    reverseKey = ''.join(signalSeq2)+''.join(signalSeq1)
    if(ScoreDictionary.has_key(key)): return True,ScoreDictionary[key],lookupCount+1
    elif(ScoreDictionary.has_key(reverseKey)): return True,ScoreDictionary[reverseKey],lookupCount+1
    return False,0,lookupCount

def addToLookupScore(signalSeq1,signalSeq2,score,ScoreDictionary):
    key = ''.join(signalSeq1)+''.join(signalSeq2)
    ScoreDictionary[key] = score
    return