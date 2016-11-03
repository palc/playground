#!/usr/bin/env bash

help () {
echo 'usage: $$ assemblathon_reformat_stats.sh <scaffold file.fasta>'
}

if [ -e ! $1 ]; then
    help 
    exit 1
fi

name=`echo $1 | sed 's/[ .].*//'`

r1=`find . -wholename "*_R1*gz" | head -1`
file_r1=`basename $r1`
size_r1=`ls -lh $r1 | awk '{print $5}'`

r2=`find . -wholename "*_R2*gz" | head -1`
file_r2=`basename $r2`
size_r2=`ls -lh $r2 | awk '{print $5}'`

assemblathon_stats.pl $1 | sed -e 's/^[ \t]*//' | sed -e 's/ \{2,\}/\t/g' > ${name}-assemblathon_stats.txt

number_of_scaffolds1=`grep -v ">" ${name}-assemblathon_stats.txt | grep "Number of scaffolds" | awk -F "\t" '{print $1}'`
number_of_scaffolds2=`grep -v ">" ${name}-assemblathon_stats.txt | grep "Number of scaffolds" | awk -F "\t" '{print $2}'`

total_size1=`grep "Total size of scaffolds" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $1}'`
total_size2=`grep "Total size of scaffolds" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $2}'`

largest1=`grep "Longest scaffold" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $1}'`
largest2=`grep "Longest scaffold" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $2}'`

large_scaffolds1=`grep "Number of scaffolds > 1K nt" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $1}'`
large_scaffolds2=`grep "Number of scaffolds > 1K nt" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $2}'`

n501=`grep "N50 scaffold length" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $1}'`
n502=`grep "N50 scaffold length" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $2}'`

l501=`grep "L50 scaffold count" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $1}'`
l502=`grep "L50 scaffold count" ${name}-assemblathon_stats.txt | awk -F "\t" '{print $2}'` 

printf "Sample\tFile_R1\tSize_R1\tFile_R2\tSize_R2\t$number_of_scaffolds1\t$total_size1\t$largest1\t$large_scaffolds1\t$n501\t$l501\n" > ${name}-assemblathon_reformat_stats.txt
printf "$name\t$file_r1\t$size_r1\t$file_r2\t$size_r2\t$number_of_scaffolds2\t$total_size2\t$largest2\t$large_scaffolds2\t$n502\t$l502\n" >> ${name}-assemblathon_reformat_stats.txt

# create 2016-11-03 stuber
