import glob
import os
import httplib2
import json
import logging
import traceback as tb
import suds.metrics as metrics
from tests import *
from suds import *
from suds.client import Client
from datetime import datetime
import OMIMAnnotation
import sys
sys.path.append('../')
 
ClusterFilePattern="Cluster*.csv"
errors = 0
geneClusters = {}
geneSymbolDict = {}
geneDescriptions = {}

# Input : directory containg the cluster files of the format Cluster_x*.csv
def Geneclusters(clusterFolderName):
    print " Geneclusters ..."
    #1. for each cluster file.
    geneClusters["name"] = clusterFolderName.split("/")[-1]
    geneClusters["children"] = []
    geneClusters["type"] = "root"

    for clusterFileName in glob.glob(clusterFolderName+"/"+ClusterFilePattern):
        clusterName = (clusterFileName.split("/")[-1]).split(".")[0]
        print "     Gene Cluster  :",clusterName
        #2.     Parse Gene list.
        print "     "
        GeneSymbolList = readGeneSymbols(clusterFileName)
        print " NoOfGGeneSymbs",len(GeneSymbolList)
        #3. Convert to ensemble ids
        EnsemblIdList = getEnsemblIds(GeneSymbolList)
        print " NoOfEnsemble Symbols",len(EnsemblIdList)
        print "     Writing david file ..."
        #4. Generate reports from David.
        getDavidReports(GeneSymbolList,clusterName,clusterFolderName)
        '''
        try:
            getDavidReports(EnsemblIdList,clusterName,clusterFolderName)
        except:
            print "Exception generating report for ",clusterName
        '''
    f = open(clusterFolderName + "/geneCluster.json", 'w')
    jsonarray = json.dumps(geneClusters)
    f.write(jsonarray)
    f.close()

    f = open(clusterFolderName + "/geneDescriptions.txt", 'w')
    for description in geneDescriptions:
        f.write(description + "," + geneDescriptions[description] + "\n")
    f.close()
    return

def getEnsemblIds(GeneSymbolList):
    EnsemblIdList = []
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    
    for geneSymbol in GeneSymbolList:
        ext = "/xrefs/symbol/homo_sapiens/"+geneSymbol+"?"
        resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
        try:
            decoded = json.loads(content)
            EnsemblIdList.append(decoded[0]['id'])
            geneSymbolDict[decoded[0]['id']] = geneSymbol
        except:
            print "couldntConvert :",geneSymbol
    return EnsemblIdList
    

def getDavidReports(EnsemblIdList,clusterName,clusterFolderName):
    errors = 0
    setup_logging()
    logging.getLogger('suds.client').setLevel(logging.ERROR)
    
    url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
        
    client = Client(url)
    
    #authenticate user email 
    client.service.authenticate('ixu@mit.edu')
    
  
    inputIds = ""
    for i in xrange(len(EnsemblIdList)-1):
        inputIds = inputIds+EnsemblIdList[i]+","
    inputIds = inputIds+EnsemblIdList[len(EnsemblIdList)-1]
    idType='AFFYMETRIX_3PRIME_IVT_ID'
    listName = 'testList'
    listType = 0
    print "     Creating list of ",len(EnsemblIdList)," symbols"
    print client.service.addList(inputIds, idType, listName, listType)
    #category = 'BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'
    #no = client.service.setCategories(category)
    #print no

    print " ... get Gene Clusters"
    overlap=3
    initialSeed = 3
    finalSeed = 3
    linkage = 0.5
    kappa = 90
    geneClusterReportFile = open(clusterFolderName+"/"+clusterName+".GeneCluster","wb")
    geneClusterReport = client.service.getGeneClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)
    #geneClusterReport = client.service.getGeneClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)
    #geneClusterReportFile.write('Category,Term,Count,%,Pvalue,Genes,List Total,Pop Hits,Pop Total,Fold Enrichment,Bonferroni,Benjamini,FDR\n')
    
    clusters = {}
    clusters["name"] = clusterName
    clusters["children"] = []
    clusters["type"] = "cluster"
    scores = {}
    for row in geneClusterReport:
        try:
            rowDict = dict(row)
            score = rowDict['score']
            name = rowDict['name']
            scores[name] = score
            group = {}
            group["name"] = name
            group["children"] = []
            group["type"] = "group"
            group["score"] = score
            for listRecord in rowDict['listRecords']:
                listRecordDict = dict(listRecord)
                description = listRecordDict['name']
                geneId = listRecordDict['values'][0]
                geneName = geneSymbolDict[geneId]
                geneDescriptions[geneName] = description
                group["children"].append({"name": geneName, "description" : description, "type": "gene"})
            clusters["children"].append(group)
        except:
            print geneClusterReport
    if len(clusters["children"]) > 0:
        geneClusters["children"].append(clusters)
        #print rowDict['name']
    '''
        categoryName = str(rowDict['categoryName'])
        termName = str(rowDict['termName'])
        listHits = str(rowDict['listHits'])
        percent = str(rowDict['percent'])
        ease = str(rowDict['ease'])
        Genes = str(rowDict['geneIds'])
        listTotals = str(rowDict['listTotals'])
        popHits = str(rowDict['popHits'])
        popTotals = str(rowDict['popTotals'])
        foldEnrichment = str(rowDict['foldEnrichment'])
        bonferroni = str(rowDict['bonferroni'])
        benjamini = str(rowDict['benjamini'])
        FDR = str(rowDict['afdr'])
        rowList = [categoryName,termName,listHits,percent,ease,Genes,listTotals,popHits,popTotals,foldEnrichment,bonferroni,benjamini,FDR]
        chartReportFile.write(','.join(rowList)+'\r\n')'''
    geneClusterReportFile.close()

    #summaryReport = client.service.getSummaryReport()
    #listReport = client.service.getListReport()

    return 

def readGeneSymbols(clusterFileName):
    f = open(clusterFileName,"rb")
    GeneSymbolList = []
    for line in f.readlines():
        tokens = line.split(',');
        if(tokens[0] == "0"): continue
        token = tokens[0].split(":")[0]# splitting tokens:transcript
        token = token.replace("\"","")
        token = token.replace("\'","")
        GeneSymbolList.append(token) 
    return GeneSymbolList

def setup_logging():
    if sys.version_info < (2, 5):
        fmt = '%(asctime)s [%(levelname)s] @%(filename)s:%(lineno)d\n%(message)s\n'
    else:
        fmt = '%(asctime)s [%(levelname)s] %(funcName)s() @%(filename)s:%(lineno)d\n%(message)s\n'
    logging.basicConfig(level=logging.INFO, format=fmt)

Geneclusters('../Output/GSE675_0/Heirarchical_NoOfClusters_9')

        
    