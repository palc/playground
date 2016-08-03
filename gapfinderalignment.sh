#!/bin/sh

picard='/usr/local/bin/picard-tools-1.117/'
gatk='/usr/local/bin/GenomeAnalysisTK/GenomeAnalysisTK.jar'
igvtools='/usr/local/bin/IGVTools/igvtools.jar'

ref=$1
forReads=$2
revReads=$3

r=`echo $ref | sed 's/\..*//'`
n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`

if [[ $ref == *.fasta ]]; then
	echo "Reference: $ref"
else
	echo "Bad argument! Argument 1 must be a fasta file ending in .fasta"
	echo 'Run using: $ quickassembly.sh reference.fasta sample_R1.fastq(.gz) sample_R2.fastq(.gz)' 
	exit 1
fi

if [[ $forReads == *.fastq || $forReads == *.fastq.gz ]]; then
	echo "Read 1: $forReads"
else
	echo "Bad argument! Argument 2 must be a fastq file ending in .fastq or .fastq.gz"
	echo 'Run using: $ quickassembly.sh reference.fasta sample_R1.fastq(.gz) sample_R2.fastq(.gz)'
	exit 1
fi

if [[ $revReads == *.fastq || $revReads == *.fastq.gz ]]; then
	echo "Read 2: $revReads"
else
	echo "Bad argument! Argument 3 must be a fastq file ending in .fastq or .fastq.gz"
	echo 'Run using: $ quickassembly.sh reference.fasta sample_R1.fastq(.gz) sample_R2.fastq(.gz)'
fi

samtools faidx $ref
java -Xmx4g -jar ${picard}CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${r}.dict

if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
    echo "Index and dict are present, continue script"
    else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx4g -jar ${picard}CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${r}.dict
        if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
        fi
fi

# See echo comments
echo "***bwa index $r"
bwa index $ref

echo "***Making Sam file"
bwa mem -M -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > $n.sam

####### unmapped reads #######
#Bam with mapped and unmapped reads
#samtools view -bh -T $ref $n.sam > $n.all.bam
#Strip off the unmapped reads
#samtools view -h -f4 $n.all.bam > $n.unmappedReads.sam
#Create fastqs of unmapped reads to assemble
#java -Xmx6g -jar ${picard}SamToFastq.jar INPUT=$n.unmappedReads.sam FASTQ=${n}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-unmapped_R2.fastq
#rm $n.all.bam
#rm $n.unmappedReads.sam
#(abyss-pe name=${n}_abyss k=64 in="${n}-unmapped_R1.fastq ${n}-unmapped_R1.fastq" && blast-contigs.sh ./*-8.fa) &

#mkdir ../unmappedReads
#mv ${n}-unmapped_R1.fastq ../unmappedReads
#mv ${n}-unmapped_R2.fastq ../unmappedReads
#mv ${n}_abyss-3.fa ../unmappedReads
#mv ${n}_abyss-8.fa ../unmappedReads
#mv ${n}_abyss-stats ../unmappedReads
#mv *coverage* ../unmappedReads
#rm *abyss*
#(blast-contigs.sh ../unmappedReads/*-8.fa) &
######################

echo "***Making Bam file"
samtools view -bh -F4 -T $ref $n.sam > $n.raw.bam

echo "***Sorting Bam"
samtools sort $n.raw.bam $n.sorted
echo "***Indexing Bam"
samtools index $n.sorted.bam

# After the first awk $1 -> position $2 -> REF $3 -> ALT, $4 -> QUAL, $5 -> DP, $6 -> 0 or [1-9] telling if ref as coverage, all zeros no coverage.
java -Xmx4g -jar ${gatk} -R $ref -T UnifiedGenotyper -out_mode EMIT_ALL_SITES -I $n.sorted.bam -o $n.vcf -nct 8

java -Xmx4g -jar ${gatk} -R $ref -T UnifiedGenotyper -I $n.sorted.bam -o ${n}-SNPs.vcf -nct 8

grep -v "#" $n.vcf | awk '{print $2, $4, $5, $6, $8, $10}' | sed 's:\(A[C,N]=[1,2]\).*DP=\([0-9]*\).*:\1 \2:' | awk '{if ($4 !~ /[0-9]+/) print "N"; else if ( $6 < 2 ) print "N"; else if ( $4 > 700 && $3 ~ /[A,T,C,G]/ && $5 ~ /AC=2/ || $4 > 1800 && $3 ~ /[A,T,C,G]/ && $5 ~ /AC=1/) print $3; else print $2}' | sed 's/[ATGC]*,[ATGC]*/N/g' | tr -d "\n" | sed "s/^/>${n};/" | tr ";" "\n" > ${n}-consensus.fasta


wait

#rm *.fasta.amb
#rm *.fasta.ann
#rm *.fasta.bwt
#rm *.dict
#rm *.fasta.pac
#rm *.fasta.sa
rm $n.raw.bam
rm $n.sam

# ~~ Author: Tod Stuber Created: 2014-10-30
