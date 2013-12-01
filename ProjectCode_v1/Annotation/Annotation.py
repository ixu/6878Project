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
 
ClusterFilePattern="Cluster*.csv"
errors = 0


# Input : directory containg the cluster files of the format Cluster_x*.csv
def Annotateclusters(clusterFolderName):
    print " Annotateclusters ..."
    print "     Getting the OMIM Dictionary ..."
    OMIMDict = OMIMAnnotation.GetOMIMDictionary()
    #1. for each cluster file.
    for clusterFileName in glob.glob(clusterFolderName+"/"+ClusterFilePattern):
        clusterName = (os.path.basename(clusterFileName)).split(".")[0]
        print "     Annotating Cluster  :",clusterName
        #2.     Parse Gene list and write omim content.
        print "     "
        GeneSymbolList = readGeneSymbols(clusterFileName)
        print " NoOfGGeneSymbs",len(GeneSymbolList)
        print "     Writing omim file ..."
        writeOMIMFile(GeneSymbolList,OMIMDict,clusterName,clusterFolderName)
        #3. Convert to ensemble ids
        EnsemblIdList = getEnsemblIds(GeneSymbolList)
        print " NoOfEnsemble Symbols",len(EnsemblIdList)
        print "     Writing david file ..."
        #4. Generate reports from David.
        try:
            getDavidReports(EnsemblIdList,clusterName,clusterFolderName)
        except:
            print "Exception generating report for ",clusterName
    return

def writeOMIMFile(GeneSymbolList,OMIMDict,clusterName,clusterFolderName):
    omimFile = open(clusterFolderName+"/"+clusterName+"_OMIMSummary.txt","wb")
    for geneSymbol in GeneSymbolList:
        if(OMIMDict.has_key(geneSymbol)):
            omimFile.write(OMIMDict[geneSymbol])
    omimFile.close()
    return

GeneSymbolToEnsembIDDict = {}

def getEnsemblIds(GeneSymbolList):
    EnsemblIdList = []
    http = httplib2.Http(".cache")
    server = "http://beta.rest.ensembl.org"
    
    for geneSymbol in GeneSymbolList:
        if(GeneSymbolToEnsembIDDict.has_key(geneSymbol)):
            if(GeneSymbolToEnsembIDDict[geneSymbol] != "couldntConvert"):
                EnsemblIdList.append(GeneSymbolToEnsembIDDict[geneSymbol])
        else:
            ext = "/xrefs/symbol/homo_sapiens/"+geneSymbol+"?"
            resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
            try:
                decoded = json.loads(content)
                EnsemblIdList.append(decoded[0]['id'])
                GeneSymbolToEnsembIDDict[geneSymbol]=decoded[0]['id']
            except:
                print "couldntConvert :",geneSymbol
                GeneSymbolToEnsembIDDict[geneSymbol]="couldntConvert :"
    return EnsemblIdList
    

def getDavidReports(EnsemblIdList,clusterName,clusterFolderName):
    errors = 0
    setup_logging()
    logging.getLogger('suds.client').setLevel(logging.ERROR)
    
    url = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl'
        
    client = Client(url)
    
    #authenticate user email 
    client.service.authenticate('dalesh@mit.edu')
    
  
    inputIds = ""
    for i in xrange(len(EnsemblIdList)-1):
        inputIds = inputIds+EnsemblIdList[i]+","
    inputIds = inputIds+EnsemblIdList[len(EnsemblIdList)-1]
    idType='ENSEMBL_GENE_ID'
    listName = 'testList'
    listType = 0
    print "     Creating list of ",len(EnsemblIdList)," symbols"
    print client.service.addList(inputIds, idType, listName, listType)
    category = 'BBID,BIOCARTA,COG_ONTOLOGY,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'
    no = client.service.setCategories(category)
    species = client.service.getCurrentSpecies()
    #allspecies = client.service.getSpecies ()
    #client.service.setCurrentSpecies(string)
    #getChartReport
    print " ... get Chart Report"
    thd=0.1
    count = 2
    chartReportFile = open(clusterFolderName+"/"+clusterName+".ChartReport","wb")
    chartReport = client.service.getChartReport(thd, count)
    chartReportFile.write('Category,Term,Count,%,Pvalue,Genes,List Total,Pop Hits,Pop Total,Fold Enrichment,Bonferroni,Benjamini,FDR\n')
    for row in chartReport:
        rowDict = dict(row)
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
        chartReportFile.write(','.join(rowList)+'\r\n')
    chartReportFile.close()
    tableReportFile = open(clusterFolderName+"/"+clusterName+".TableReport","wb")
    tableReport = client.service.getTableReport()
    tableReportFile.write('name,category,terms')
    for record in tableReport:
        for annot in record['annotationRecords']:
            tableReportFile.write(record['name']+","+annot['category']+str(annot['terms'])+"\r\n")
    tableReportFile.close()
    summaryReport = client.service.getSummaryReport()
    listReport = client.service.getListReport()

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


        
    