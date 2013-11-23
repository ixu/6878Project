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

ClusterFilePattern="Cluster*.csv"
errors = 0


# Input : directory containg the cluster files of the format Cluster_x*.csv
def Annotateclusters(clusterFolderName):
    print " Annotateclusters ..."

    #1. for each cluster file.
    for clusterFileName in glob.glob(clusterFolderName+"/"+ClusterFilePattern):
        clusterName = (clusterFileName.split("\\")[1]).split(".")[0]
        #2.     Parse Gene list
        print "     "
 
        GeneSymbolList = readGeneSymbols(clusterFileName)
        print " NoOfGGeneSymbs",len(GeneSymbolList)
        #3. Convert to ensemble ids
        EnsemblIdList = getEnsemblIds(GeneSymbolList)
        print " NoOfEnsemble Symbols",len(EnsemblIdList)
        #4. Generate reports from David.
        try:
            getDavidReports(EnsemblIdList,clusterName,clusterFolderName)
        except:
            print "Exception generating report for ",clusterName
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
        chartReportFile.write(','.join(rowList)+'\n')
    chartReportFile.close()
    tableReportFile = open(clusterFolderName+"/"+clusterName+".TableReport","wb")
    tableReport = client.service.getTableReport()
    tableReportFile.write('name,category,terms')
    for record in tableReport:
        for annot in record['annotationRecords']:
            tableReportFile.write(record['name']+","+annot['category']+str(annot['terms'])+"\n")
    tableReportFile.close()

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


        
    