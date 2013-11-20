#!/usr/bin/env python

import sys

base_idx = { 'R' : 0, 'D' : 1, 'S' : 2}
PTR_NONE, PTR_GAP1, PTR_GAP2, PTR_BASE = 0, 1, 2, 3
S = [
     # R D S
     [3, -1, -3], # R
     [-1, 3, -3], # D
     [-3, -3, 1]  # S
     ]
gap_pen = 1

def seqalignDP(seq1,seq2,S,gap_pen,timelines):
    """return the score of the optimal Needleman-Wunsch alignment for seq1 and seq2
        Note: gap_pen should be positive (it is subtracted)
        """

    subst_matrix=S
    F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
    
    # initialize dynamic programming table for Needleman-Wunsch alignment (Durbin p.20)
    for i in range(1,len(seq1)+1):
        F[i][0] = F[i-1][0] - getGapPenalty(i-1,i,gap_pen,timelines)
    for j in range(1,len(seq2)+1):
        F[0][j] = F[0][j-1] - getGapPenalty(j-1,j,gap_pen,timelines)
    
    
    
    # YOUR CODE HERE
    # Fill in the dynamic programming tables F and TB, starting at [1][1]
    # Hints: The first row and first column of the table F[i][0] and F[0][j] are dummies
    #        (see for illustration Durbin p.21, Figure 2.5, but be careful what you
    #         think of as rows and what you think of as columns)
    #        Hence, the bases corresponding to F[i][j] are actually seq1[i-1] and seq2[j-1].
    #        Use the dictionary base_idx to convert from the character to an index to
    #         look up entries of the substitution matrix.
    #        To get started, you can complete and run the algorithm filling in only F,
    #         and then figure out how to do TB.
    
    
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            seq1Base = base_idx[seq1[i - 1]]
            seq2Base = base_idx[seq2[j - 1]]
            score = S[seq1Base][seq2Base]
            max_choice = max(F[i][j - 1] - getGapPenalty(i-1,j-2,gap_pen,timelines), F[i - 1][j - 1] + score, F[i - 1][j] - getGapPenalty(i-2,j-1,gap_pen,timelines))
            F[i][j] = max_choice

    return F[len(seq1)][len(seq2)]

def getGapPenalty(i,j,gap_pen,timelines):
    if((i>(len(timelines)-1)) | (j>(len(timelines)-1))): return 1000 # last row and last column
    timeLineThreshold = (max(timelines)-min(timelines))/len(timelines)
    if(abs(timelines[j]-timelines[i]) > timeLineThreshold):
        return 1000 # basically we penalize with a huge gap penalty if the distance is too big in time.
    else:
        return gap_pen

def traceback(seq1,seq2,TB):
    s1 = ""
    s2 = ""
    
    i = len(seq1)
    j = len(seq2)
    
    while TB[i][j] != PTR_NONE:
        if TB[i][j] == PTR_BASE:
            s1 = seq1[i-1] + s1
            s2 = seq2[j-1] + s2
            i=i-1
            j=j-1
        elif TB[i][j] == PTR_GAP1:
            s1 = '-' + s1
            s2 = seq2[j-1] + s2
            j=j-1
        elif TB[i][j] == PTR_GAP2:
            s1 = seq1[i-1] + s1
            s2 = '-' + s2
            i=i-1
        else: assert False
    
    return s1,s2

def readSeq(filename):
    """reads in a FASTA sequence"""
    
    stream = open(filename)
    seq = []
    
    for line in stream:
        if line.startswith(">"):
            continue
        seq.append(line.rstrip())
    
    return "".join(seq)




def main():
    # parse commandline
    if len(sys.argv) < 3:
        print "you must call program as: python ps1-seqalign.py <FASTA 1> <FASTA 2>"
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    
    seq1 = readSeq(file1)
    seq2 = readSeq(file2)
    
    score, F, TB = seqalignDP(seq1,seq2,S,gap_pen)
    
    print >> sys.stderr, score
    
    s1, s2 = traceback(seq1,seq2,TB)
    print s1
    print s2

if __name__ == "__main__":
    main()
