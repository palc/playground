#!/bin/sh

infile = $1

read1=`grep -A 2 "fastq.gz file sizes" $infile | tail -1`
read2=`grep -A 2 "fastq.gz file sizes" $infile | tail -2 | head -1`

ave_coverage=`grep "Average.*coverage" $infile | awk '{print $NF}'`
ref_with_coverage=`grep "eference.*coverage" $infile | awk '{print $NF}'`

mean_length=`grep -A 1 "Mean_Read_Length" $infile | tail -1`
quality_snps=`grep -A 1 "SNPs of AC2" $infile | tail -1`

printf "$read1\t$read2\t$ave_coverage\t$ref_with_coverage\t$mean_length\t$quality_snps\n"


# tstuber 2017-04-19
