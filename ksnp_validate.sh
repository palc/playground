#!/bin/sh

# Summary...
# Use AlleleCounts.ML.NodeLabel.tre
# Then ClusterInfo.ML
# Then kSNP VCF to get positions
# Cross reference those positions to GATK VCFs to get quality info

### !!! VERIFY VCF REFERENCES EQUAL BETWEEN ALIGNMENT AND KSNP" !!! ###

# location of vcfs
vcf_dir="/scratch/analysis_staging/ksnp/h37_vcfs"

# current working directory must be kSNP analysis run directory containing files
	# ClusterInfo.ML
	# *vcf file
# kSNP script must have output vcf

# Argument requirement, $1 = node number

if [ -z ${1} ]; then
	echo "provide arg1"
	echo "set arg1 to node number"
	exit 1
fi

num=$1
echo "getting context sequence"
# get context sequence for node
echo "finding $num positions"
contextseq=($(grep "Node.${num}" ClusterInfo.ML | awk '{print $2}'))

rm Node.${num}.positions*

echo "getting positions"
# from context sequence get the positions from VCF output by kSNP
count=1
for c in "${contextseq[@]}"; do 
	myposition=`grep "$c" *vcf | awk '{print $2}'`  
	echo "*** $count using context: $c position found: $myposition" | tee -a Node.${num}.positions.detail Node.${num}.positions
	awk -v p="${myposition}" '$2 == p {print FILENAME, $2, $3, $4, $5, $6, $7, $8}' $vcf_dir/* | sed "s:${vcf_dir}/::" | column -t | tee -a Node.${num}.positions.detail
	echo "----------------------------" | tee -a Node.${num}.positions.detail
	echo "" | tee -a Node.${num}.positions.detail
	let count=count+1
done

# extracting context sequence
# sed 's/.*using context: \(.*\) position found: .*/\1/' Node.78.positions
# 2016-07-29 stuber
