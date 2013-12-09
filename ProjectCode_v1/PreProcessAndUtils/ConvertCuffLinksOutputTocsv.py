import sys
import os

TRANSCRIPT_FILE_NAME = "transcripts.gtf"
ISOFORMS_FILE_NAME = "isoforms.fpkm_tracking"
GENES_FILE_NAME = "genes.fpkm_tracking"
INDEX_PREFIX = "GSM8185"
FilterFPKMValue = 0.000
# given a directory search for all transcripts.gtf files
# read each transcripts.gtf file and create a matrix of genes vs fpkm values.



#########################################################################
##### getTranscriptFilesDictionary
#########################################################################
# searches for transcripts.gtf and creates a dictionary with the GSM** 
#   number as the key and filename as value.
# Input : directory to search
# output : dictionary as described above.
#########################################################################
def getFilesDictionary(inputDirectory,fileNameType):
    print " Processing ",inputDirectory,"...."
    transcriptFilesDictionary = {}
    for root, dirs, files in os.walk(inputDirectory):
        for file in files:
            if file.endswith(fileNameType):
                fileName = os.path.join(root, file)
                pos = fileName.find(INDEX_PREFIX)
                indexKey = fileName[pos:pos+9]
                transcriptFilesDictionary[indexKey] = fileName
    return transcriptFilesDictionary

#########################################################################
##### fetchmRNAExpressionForFile
#########################################################################
# runs through each line looking for transcript entries and fills the rna expression
# data with FPKM values
#########################################################################
def fetchmRNAExpressionForTranscriptFile(fileName,mRNAExpressionData):
    stream = open(fileName)
    for line in stream:
        tokens = line.split()
        if(tokens[2] == 'transcript'):
            geneId = tokens[9][1:len(tokens[9])-2]
            transcriptId=tokens[11][1:len(tokens[11])-2]
            geneTranscriptId = geneId+":"+transcriptId
            FPKMValue = float(tokens[13][1:len(tokens[13])-2])
            if(mRNAExpressionData.has_key(geneTranscriptId)):
                mRNAExpressionData[geneTranscriptId].append(FPKMValue)
            else:
                mRNAExpressionData[geneTranscriptId] = [FPKMValue]
    stream.close()
    return

def fetchmRNAExpressionForIsoformFile(fileName,mRNAExpressionData):
    stream = open(fileName)
    for line in stream:
        tokens = line.split()
        if(tokens[0] != 'tracking_id'):
            geneId = tokens[3]
            transcriptId=tokens[0]
            geneTranscriptId = geneId+":"+transcriptId
            FPKMValue = float(tokens[9])
            if(mRNAExpressionData.has_key(geneTranscriptId)):
                mRNAExpressionData[geneTranscriptId].append(FPKMValue)
            else:
                mRNAExpressionData[geneTranscriptId] = [FPKMValue]
    stream.close()
    print " No Of Keys ",len(mRNAExpressionData.keys())
    return

def fetchmRNAExpressionForGenesFile(fileName,mRNAExpressionData):
    stream = open(fileName)
    for line in stream:
        tokens = line.split()
        if(tokens[0] != 'tracking_id'):
            geneId = tokens[0]
            geneTranscriptId = geneId
            FPKMValue = float(tokens[9])
            if(mRNAExpressionData.has_key(geneTranscriptId)):
                if(len(mRNAExpressionData[geneTranscriptId])<20):
                   mRNAExpressionData[geneTranscriptId].append(FPKMValue)
            else:
                mRNAExpressionData[geneTranscriptId] = [FPKMValue]
    stream.close()
    print " No Of Keys ",len(mRNAExpressionData.keys())
    return
#########################################################################
##### writeRNAExpressionasCSV
#########################################################################
# runs through the dictionary and writes it to a file as a , seperated file.
#########################################################################
def writeRNAExpressionasCSV(outputFileName,mRNAExpressionData):
    g = open(outputFileName,"w")
    # write the header
    g.write("GeneID:TranscriptId")
    for columnHeader in mRNAExpressionData["HEADER"]:
        if(columnHeader != "GeneID:TranscriptId"):
            g.write(","+columnHeader)
    g.write("\n")
    for key in mRNAExpressionData:
        if key != "HEADER":
            if key == "" or key == " ": print "mess exists"
            g.write(key)
            for value in mRNAExpressionData[key]:
                g.write(","+str(value))
        if key != "HEADER":
            g.write("\n")
    g.close()
    return

def ProcessTranscriptFiles(inputDirectory):
    transcriptFilesDictionary = getFilesDictionary(inputDirectory,TRANSCRIPT_FILE_NAME)
    sortedKeys = []
    for key in transcriptFilesDictionary:
        sortedKeys.append(key)
    sortedKeys.sort()
    mRNAExpressionData = {}
    mRNAExpressionData["HEADER"] = ["GeneID:TranscriptId"]
    for key in sortedKeys:
        mRNAExpressionData["HEADER"].append(key)
        print " Processing " + key + "..."
        fetchmRNAExpressionForTranscriptFile(transcriptFilesDictionary[key],mRNAExpressionData)
    prev  = ""
    # if atleast 1 value is greater than the filterValue include that gene:transcript
    mRNAExpressionDataFiltered = {}
    for key in mRNAExpressionData:
        if key == 'GeneID:TranscriptId':
            mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
        else:
            for value in mRNAExpressionData[key]:
                if(value > FilterFPKMValue):
                    mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
                    break
        prev = key
    writeRNAExpressionasCSV("mRNAExpressionData_transcript.csv",mRNAExpressionData)
    writeRNAExpressionasCSV("mRNAExpressionDataFiltered_transcript_"+str(FilterFPKMValue)+".csv",mRNAExpressionDataFiltered)

def ProcessIsoFormsFiles(inputDirectory):
    transcriptFilesDictionary = getFilesDictionary(inputDirectory,ISOFORMS_FILE_NAME)
    sortedKeys = []
    for key in transcriptFilesDictionary:
        sortedKeys.append(key)
    sortedKeys.sort()
    mRNAExpressionData = {}
    mRNAExpressionData["HEADER"] = ["GeneID:TranscriptId"]
    for key in sortedKeys:
        mRNAExpressionData["HEADER"].append(key)
        print " Processing " + key + "..."
        fetchmRNAExpressionForIsoformFile(transcriptFilesDictionary[key],mRNAExpressionData)
    
    # if atleast 1 value is greater than the filterValue include that gene:transcript
    mRNAExpressionDataFiltered = {}
    for key in mRNAExpressionData:
        if key == 'GeneID:TranscriptId':
            mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
        else:
            for value in mRNAExpressionData[key]:
                if(value > FilterFPKMValue):
                    mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
                    break
    print " No Of Keys afterFiltering",len(mRNAExpressionDataFiltered.keys())
    writeRNAExpressionasCSV("mRNAExpressionData_isoforms.csv",mRNAExpressionData)
    writeRNAExpressionasCSV("mRNAExpressionDataFiltered_isoforms_"+str(FilterFPKMValue)+".csv",mRNAExpressionDataFiltered)
    return

def ProcessGenesFiles(inputDirectory):
    transcriptFilesDictionary = getFilesDictionary(inputDirectory,GENES_FILE_NAME)
    sortedKeys = []
    for key in transcriptFilesDictionary:
        sortedKeys.append(key)
    sortedKeys.sort()
    mRNAExpressionData = {}
    mRNAExpressionData["HEADER"] = ["GeneID:TranscriptId"]
    for key in sortedKeys:
        mRNAExpressionData["HEADER"].append(key)
        print " Processing " + key + "..."
        fetchmRNAExpressionForGenesFile(transcriptFilesDictionary[key],mRNAExpressionData)
    
    # if atleast 1 value is greater than the filterValue include that gene:transcript
    mRNAExpressionDataFiltered = {}
    for key in mRNAExpressionData:
        if key == 'GeneID:TranscriptId':
            mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
        else:
            for value in mRNAExpressionData[key]:
                if(value > FilterFPKMValue):
                    mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
                    break
    writeRNAExpressionasCSV("mRNAExpressionData_genes.csv",mRNAExpressionData)
    writeRNAExpressionasCSV("mRNAExpressionDataFiltered_genes_"+str(FilterFPKMValue)+".csv",mRNAExpressionDataFiltered)
   

def ProcessIsoFormsFilesFilteringTrials(inputDirectory):
    transcriptFilesDictionary = getFilesDictionary(inputDirectory,ISOFORMS_FILE_NAME)
    sortedKeys = []
    for key in transcriptFilesDictionary:
        sortedKeys.append(key)
    sortedKeys.sort()
    mRNAExpressionData = {}
    mRNAExpressionData["HEADER"] = ["GeneID:TranscriptId"]
    for key in sortedKeys:
        mRNAExpressionData["HEADER"].append(key)
        print " Processing " + key + "..."
        fetchmRNAExpressionForIsoformFile(transcriptFilesDictionary[key],mRNAExpressionData)
    ranges = [0.005,0.006,0.007,0.01,0.05,0.5]
    for range in ranges:
        FilterFPKMValue = range
        # if atleast 1 value is greater than the filterValue include that gene:transcript
        mRNAExpressionDataFiltered = {}
        for key in mRNAExpressionData:
           if key == 'GeneID:TranscriptId':
               mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
           else:
               for value in mRNAExpressionData[key]:
                   if(value > FilterFPKMValue):
                      mRNAExpressionDataFiltered[key] = mRNAExpressionData[key]
                      break
        print " No Of Keys afterFiltering for Range",range," : ",len(mRNAExpressionDataFiltered.keys())
    return
    


def main():
    if len(sys.argv) < 2:
        print "you must call program as:  "
        print "   python ConvertTranscriptGTFToCsvMatrix.py directoryContainingGTFfiles"

    inputDirectory = sys.argv[1]
    #ProcessIsoFormsFilesFilteringTrials(inputDirectory)
    ProcessTranscriptFiles(inputDirectory)
    ProcessIsoFormsFiles(inputDirectory)
    ProcessGenesFiles(inputDirectory)
    return
    

main()
