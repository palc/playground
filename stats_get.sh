#!/bin/sh

infile=$1
<<<<<<< HEAD

name=`echo $infile | sed 's/\.stats\.txt//' | sed 's/stats-//' | sed 's/\.txt//'`
=======
>>>>>>> c400466b98d232f527c3b1e64878dcc86dc65404

name=`echo $infile | sed 's/.stats.txt//'`

read1=`grep -A 2 "^fastq.gz file sizes:" $infile | tail -1`
read2=`grep -A 2 "^fastq.gz file sizes:" $infile | tail -2 | head -1`

ave_coverage=`grep "Average.*coverage" $infile | awk '{print $NF}'`
ref_with_coverage=`grep "eference.*coverage" $infile | awk '{print $NF}'`

mean_length=`grep -A 1 "Mean_Read_Length" $infile | tail -1`
quality_snps=`grep -A 1 "SNPs of AC2" $infile | tail -1`

<<<<<<< HEAD
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$name" "$read1" "$read2" "$ave_coverage" "$ref_with_coverage" "$mean_length" "$quality_snps"
=======
echo  "$name $read1 $read2 $ave_coverage $ref_with_coverage $mean_length $quality_snps"

>>>>>>> c400466b98d232f527c3b1e64878dcc86dc65404

# tstuber 2017-04-19
