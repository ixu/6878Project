import seqalign
import timeit
base_idx = { 'R' : 0, 'D' : 1, 'S' : 2}
def getAlignmentScoreMatrix(signalSeqs,S,gap_pen,timeline):
    scores = [[0] * len(signalSeqs) for i in range(len(signalSeqs))]
    startFull = timeit.default_timer()
    gapPenLut = seqalign.createGapPenaltyLUT(len(signalSeqs[0]),gap_pen,timeline)
    for i in range(len(signalSeqs)):
        if(i%100 == 0): print i
        for j in range(i + 1, len(signalSeqs)):
            #scores[i][j] = seqalign.seqalignDP(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline)
            scores[i][j] = seqalign.seqalignDPFast(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline,gapPenLut)
            scores[j][i] = scores[i][j]
    endFull = timeit.default_timer()
    print "Time for creating matrix : ",endFull - startFull
    return scores
    


def getAlignmentScoreMatrixWithLookup(signalSeqs,S,gap_pen,timeline):
    ScoreDictionary = {}
    lookupCount = 0
    startFull = timeit.default_timer()
    scores = [[0] * len(signalSeqs) for i in range(len(signalSeqs))]

    gapPenLut = seqalign.createGapPenaltyLUT(len(signalSeqs[0]),gap_pen,timeline)
    for i in range(len(signalSeqs)):
        if(i%100==0): print i
        for j in range(i + 1, len(signalSeqs)):
            doesScoreExist,lookupScore,lookupCount = lookUpScore(signalSeqs[i],signalSeqs[j],ScoreDictionary,lookupCount)
            if(doesScoreExist == True): scores[i][j] = lookupScore
            else:
                scores[i][j] = seqalign.seqalignDP(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline)
                #scores[i][j] = seqalign.seqalignDPFast(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline,gapPenLut)
                addToLookupScore(signalSeqs[i],signalSeqs[j],scores[i][j],ScoreDictionary)
            scores[j][i] = scores[i][j]
    print " -- Successful lookups :",lookupCount
    endFull = timeit.default_timer()
    print "Time for creating matrix : ",endFull - startFull
    return scores

def getAlignmentScoreMatrixWithLookupPlusSignalKeys(signalSeqs,S,gap_pen,timeline):
    ScoreDictionary = {}
    lookupCount = 0
    startFull = timeit.default_timer()
    scores = [[0] * len(signalSeqs) for i in range(len(signalSeqs))]
    seqKeys = {}
    for i in xrange(len(signalSeqs)):
        signalKey = 0
        for j in xrange(len(signalSeqs[i])):
            signalKey = signalKey*3+base_idx[signalSeqs[i][j]]
        seqKeys[i] = signalKey

    gapPenLut = seqalign.createGapPenaltyLUT(len(signalSeqs[0]),gap_pen,timeline)
    for i in range(len(signalSeqs)):
        if(i%100==0): print i
        for j in range(i + 1, len(signalSeqs)):
            doesScoreExist,lookupScore,lookupCount = lookUpScoreWithKeys(signalSeqs[i],signalSeqs[j],ScoreDictionary,lookupCount,seqKeys,i,j)
            if(doesScoreExist == True): scores[i][j] = lookupScore
            else:
                scores[i][j] = seqalign.seqalignDP(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline)
                #scores[i][j] = seqalign.seqalignDPFast(signalSeqs[i],signalSeqs[j],S,gap_pen,timeline,gapPenLut)
                addToLookupScoreWithKeys(signalSeqs[i],signalSeqs[j],scores[i][j],ScoreDictionary,seqKeys,i,j)
            scores[j][i] = scores[i][j]
    print " -- Successful lookups :",lookupCount
    endFull = timeit.default_timer()
    print "Time for creating matrix : ",endFull - startFull
    return scores

def lookUpScoreWithKeys(signalSeq1,signalSeq2,ScoreDictionary,lookupCount,seqKeys,i,j):
    key = (seqKeys[i],seqKeys[j])
    reverseKey = (seqKeys[j],seqKeys[i])
    if(ScoreDictionary.has_key(key)): return True,ScoreDictionary[key],lookupCount+1
    elif(ScoreDictionary.has_key(reverseKey)): return True,ScoreDictionary[reverseKey],lookupCount+1
    return False,0,lookupCount

def addToLookupScoreWithKeys(signalSeq1,signalSeq2,score,ScoreDictionary,seqKeys,i,j):
    key = (seqKeys[i],seqKeys[j])
    reverseKey = (seqKeys[j],seqKeys[i])
    ScoreDictionary[key] = score
    ScoreDictionary[reverseKey] = score

    return

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