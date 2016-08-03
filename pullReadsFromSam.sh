#!/bin/sh

# To pull mapped reads from a bam file, input bam file as argument $1.
# Must have bam files and corresponding fastq files present in working direcotry
# Change the header type below

NR_CPUS=60 # Computer cores to use when analyzing

echo "This is the argument: $1"
name=`echo $1 | sed 's/\..*//' | sed 's/_.*//'`

samtools view -h -F4 $1 | awk '{print $1}' | grep "^M" >> $name.mappedReadHeaders

mkdir readsfound
cd readsfound
split -l 399 ../$name.mappedReadHeaders
for i in `ls`; do
        (pos=`cat ${i} | tr "\n" "|" | sed 's/|$//'`
        egrep -A3 -h ${pos} ../*.fastq > $i) &
        let count+=1
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
cat [a-z][a-z][a-z] >> ${name}.fastq
rm [a-z][a-z][a-z]
egrep -A3 "1:N:" ${name}.fastq | grep -v '^--$' > ${name}_found_R1.fastq
egrep -A3 "2:N:" ${name}.fastq | grep -v '^--$' > ${name}_found_R2.fastq

grep -v '^--$' ${name}.fastq | perl -pe 's|@|>|;s|.*||s if $.%4==3 || $.%4==0;close $ARGV if eof' > ${name}.fasta





