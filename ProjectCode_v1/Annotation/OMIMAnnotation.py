
def GetOMIMDictionary():
    FILE=open("Annotation\omim.txt","r")
    control1=0
    OMIMDict = {}
    GeneKey = "NONE"
    GeneText = ""
    lines = FILE.readlines()
    lineIndex = -1
    for line in lines:
        line=line.strip()
        lineIndex += 1
        if control1==1:
            GeneText = GeneText + line + "\t"
            control1=0
        #change abbreviations into their expanded form
        elif line=="*RECORD*":
            if GeneKey != "NONE": # write previous record
                OMIMDict[GeneKey] = GeneText
                GeneText = ""
        elif line=="*FIELD* NO":
            control1=1
        elif line=="*FIELD* TI":
            GeneText = GeneText + "TITLE "
            title = lines[lineIndex+1]
            GeneKey = title[title.find(';')+1:].strip()
        elif line=="*FIELD* TX":
            GeneText = GeneText + "DESCRIPTION "
        elif line=="*FIELD* CD":
            GeneText = GeneText + "CREATED "
        elif line=="*FIELD* CN":
            GeneText = GeneText + "CONTRIBUTORS "
        elif line=="*FIELD* ED":
            GeneText = GeneText + "EDITED "
        elif line=="*FIELD* SA":
            GeneText = GeneText + "ADDITIONAL REFERENCES "
        elif line=="*FIELD* CS":
            GeneText = GeneText + "CLINICAL SYMPTOMS "
        elif line=="*FIELD* MN":
            GeneText = GeneText + "MINIMIM "
        elif line=="*FIELD* RF":
            GeneText = GeneText + "REFERENCE "
        elif line=="*FIELD* AV":
            GeneText = GeneText + "ALLELIC VARIANT "
        else:
            GeneText = GeneText +line+ "\r\n"
    return OMIMDict


