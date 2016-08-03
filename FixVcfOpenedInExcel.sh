#!/bin/bash

#Necesary cleanup script summary
cat InputFile | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' > OutputFile

#Below are the individual steps takn...

#This is needed to remove the commas from calls as "G,T" in ALT columns.
awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | 

#This was used to remove columns behond the typical VCF.  It would remove comments
#Is not necessary
awk -F '\t' 'BEGIN{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10;}' | 

#Fixed extra tabs and mis-commas in header
#Is not necessary
sed 's:\"\(.*\)\"\">\".*:\1\">:g' | 

#Used to ensure correctness of first header
#Is not necessary
sed 's:\(##fileformat=VCFv4.1\).*:\1:g' | 

#Used to remove double commas
#Is not necessary
sed 's/\"\"/\"/g' | 

#This is necessary to fix the header indicator ##
sed 's/\"##/##/'