#!/bin/bash
FILES="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32874/suppl/GSE32874_GSM818578_unk17.bam ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32874/suppl/GSE32874_GSM818579_unk18.bam ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32874/suppl/GSE32874_GSM818580_unk19.bam ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32874/suppl/GSE32874_GSM818581_unk20.bam ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32874/suppl/GSE32874_GSM818582_unk21.bam"
echo "begin"
dirExt="Dir"
for f in $FILES
do
 echo " downloading $f ....."
 wget $f
 shortName=$(basename $f)
 echo " Processing cufflinks on $shortName ..."
 cp /home/dalesh/data/oldgenes.gtf .
./cufflinks -G ./oldgenes.gtf $shortName
 rm -rf  *.bam
 dirName=$shortName$dirExt
 echo "making $dirName"
 mkdir $dirName
 echo " Copying files...."
 cp *.gtf $dirName
 cp *.fpkm* $dirName
 rm -rf *.gtf
 rm -rf *.fpkm*
 # do something on $f
done
