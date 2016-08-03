#!/bin/sh

# Working directory must be data folder

#directoryList=`ls`

#for i in $directoryList; do 
for i in `cat myfile`; do
	echo $i
	vcf=`ls ${i}*/BWAmem-GATK/*mem.vcf`
	stats=`ls ${i}*/BWAmem-GATK/QualityValues/*stat*`
	
	# vcf 
	ac1count=`grep -v "#" $vcf | grep -c "AC=1;"` 
	ac2count=`grep -v "#" $vcf | grep -c "AC=2"`	
	allSNPcount=`grep -v "#" $vcf | awk ' $5 != "." {print $0}' | grep -c ".*"`
	passSNPcount=`grep -v "#" $vcf | egrep "AC=2;A" | awk -v Q="$QUAL" '$6 > 150' | grep -c ".*"`
	nocoverage=`awk ' $8 == "." {print $0}' $vcf | grep -c ".*"`
	mappingzero=`grep -v "#" $vcf | grep -c "MQ=0.00"` 	
	totalpositions=`grep -v "#" $vcf | grep -c ".*"`

	# stats
	refcoverage=`grep "Reference with coverage" $stats | tr -d [:alpha:] | tr -d [:blank:] | sed 's/://g'`
	unmappedcontigs=`grep -A 1 "Unmapped contig count" $stats | tr -d "\n" | tr -d [:alpha:] | tr -d [:blank:] | sed 's/://g'`
	avecoverage=`grep "Average coverage" $stats | tr -d [:alpha:] | tr -d [:blank:] | sed 's/://g'`
	
	# print
	echo "ac1count $ac1count "
	echo "ac2count $ac2count "
	echo "allSNPcount $allSNPcount "
	echo "passSNPcount $passSNPcount "
	echo "nocoverage $nocoverage "
	echo "mappingzero $mappingzero"
	echo "totalpositions $totalpositions "
	echo "refcoverage $refcoverage "
	echo "unmappedcontigs $unmappedcontigs "
	echo "avecoverage ${avecoverage}X"

done

# Created 2015-02-19, stuber
