from copy import deepcopy
from sklearn.cluster import AffinityPropagation
from numpy import matrix

# there are two versions of affinity clustering in this code ,
# the implemented version has some bugs so had to use one from sklearn.

def AffinityPropCluster(S):

    SM = matrix(S)
    af = AffinityPropagation(preference=-50,affinity="precomputed",copy=True).fit(SM)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    #Initialization: Responsibilities are initialized to the similarity between data points (r(i; j)   S[i][j]).
    #Availabilities are initialized to 0 (a(i; j)   0).
    n_clusters_ = len(cluster_centers_indices)
    print "no of clusters found byAffy prop",n_clusters_

    clusterDict = {}
    
    for i in xrange(len(S[0])):
        if(clusterDict.has_key(labels[i])): clusterDict[labels[i]].append(i)
        else:clusterDict[labels[i]] = [i]
    clusters = []
    for key in clusterDict:
        clusters.append(clusterDict[key])
    
    return clusters

# Some bug in this code , hence had to use an external affinity clustering 
def AffinityPropClusterInternal(S):
    print "Indices :",cluster_centers_indices
    dampner = 0.5
    r = deepcopy(S)
    N = len(S[0])
    a = [[0 for j in xrange(N)] for i in xrange(N)]
    exemplarConstantForITers=0
    PreviousNoOfExemplars = -1
    while exemplarConstantForITers<5:
        # responsibility update r(i,j) = S[i][j] - maxj j1 (a(i,j1) + s(i,j1) ) where j1 != j
        for i in xrange(N):
            for j in xrange(N):
                maxAplusS = 0
                for j1 in xrange(N):
                    if((j1!=j) & ((a[i][j1]+S[i][j1])>maxAplusS)): maxAplusS = a[i][j1]+S[i][j1]
                oldrij = r[i][j]
                r[i][j] = (dampner)*(S[i][j] - maxAplusS)+(1-dampner)*oldrij

        # Availability update rule a(i,j) = min (0,r(i,j)+sum(max(0,ri1,j) where i1 is not i or j) 
        for i in xrange(N):
            for j in xrange(N):
                sumRij = 0
                for i1 in xrange(N):
                    if((i1!=j) & (i1 != i)):
                        sumRij = sumRij + (r[i1][j] if (r[i1][j]>0) else 0)
                oldaij = a[i][j]
                a[i][j] =(r[j][j]+sumRij) if (r[j][j]+sumRij)<0 else 0
                a[i][j] = (dampner)*a[i][j] +(1-dampner)*oldaij

        # self availability update rule a(j,j) = sum(max(0,ri1,j) where i1 is not j)
        for j in xrange(N):
            sumRij = 0
            for i1 in xrange(N):
                if(i1 != j): sumRij = sumRij + (r[i1][j] if (r[i1][j]>0) else 0)
            oldaij = a[i][j]
            a[j][j] = sumRij
            a[i][j] = (dampner)*a[i][j] +(1-dampner)*oldaij
        # count no of exemplars maximizes a(i; j) + r(i; j).
        exemplars = []
        for i in xrange(N):
            exemplar = -1
            max = -999999
            for j in xrange(N):
                if((a[i][j]+r[i][j])>max):
                    max = a[i][j]+r[i][j]
                    exemplar = j
            exemplars.append(exemplar)
        exemplarSet = set(exemplars)
        if(len(exemplarSet) != PreviousNoOfExemplars): 
            PreviousNoOfExemplars = len(exemplarSet)
            exemplarConstantForITers = 0
        else:
            exemplarConstantForITers += 1
        print " No Of Exemplars ",len(exemplarSet)
    #Identifying exemplars: We identify the exemplar of point i by nding the point j that maximizes
    #a(i; j) + r(i; j). If i = j, then point i is the exemplar.
    clusterDict = {}
    
    for i in xrange(N):
            exemplar = -1
            max = -999999
            for j in xrange(N):
                if((a[i][j]+r[i][j])>max):
                    max = a[i][j]+r[i][j]
                    exemplar = j
            if(clusterDict.has_key(exemplar)): clusterDict[exemplar].append(i)
            else:
                clusterDict[exemplar] = [i]
    clusters = []
    for key in clusterDict:
        clusters.append(clusterDict[key])
    
    return clusters