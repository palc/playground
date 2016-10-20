#!/bin/sh

#ref=$1
#forReads=$2
#revReads=$3
echo "File info."
ls -lh *
echo ""

picard=`which picard.jar`
if [[ -z $picard ]]; then
    echo "picard.jar not in PATH"
    echo "picard version >1.14"
    echo "Add picard.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

gatk=`which GenomeAnalysisTK.jar`
if [[ -z $gatk ]]; then
    echo "GenomeAnalysisTK.jar is not in PATH"
    echo "Add GenomeAnalysisTK.jar to PATH"
    echo "See line: $LINENO"
    exit 1
fi

# Grab reads and reference and place them in variables
ref=`ls | grep .fasta`
echo "Reference Input:  $ref"

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

###
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
java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict

if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
    echo "Index and dict are present, continue script"
    else
    sleep 5
    echo "Either index or dict for reference is missing, try making again"
    samtools faidx $ref
    java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
        if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
        read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
        fi
fi

# See echo comments
echo "***bwa index $r"
bwa index $ref

echo "***Making Sam file"
#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
bwa mem -M -B 8 -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > $n.sam

samtools view -bh -T $ref $n.sam > $n.all.bam
#Strip off the unmapped reads
#samtools view -h -f4 $n.all.bam > $n.unmappedReads.sam

#Create fastqs of unmapped reads to assemble
#java -Xmx6g -jar ${picard}SamToFastq.jar INPUT=$n.unmappedReads.sam FASTQ=${n}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-unmapped_R2.fastq
#rm $n.all.bam
#rm $n.unmappedReads.sam
#(abyss-pe name=${n}_abyss k=64 in="${n}-unmapped_R1.fastq ${n}-unmapped_R1.fastq" && blast-contigs.sh ./*-8.fa) &
#wait
#mkdir ./unmappedReads
#wait
#mv ${n}-unmapped_R1.fastq ./unmappedReads
#mv ${n}-unmapped_R2.fastq ./unmappedReads
#mv ${n}_abyss-3.fa ./unmappedReads
#mv ${n}_abyss-8.fa ./unmappedReads
#mv ${n}_abyss-stats ./unmappedReads
#mv *coverage* ./unmappedReads
#rm *abyss*
#(blast-contigs.sh ../unmappedReads/*-8.fa) &
######################

#Keep mapped reads 
samtools view -h -F4 $n.all.bam > $n.mappedReads.sam
#java -Xmx6g -jar ${picard} SamToFastq INPUT=$n.mappedReads.sam FASTQ=${n}-mapped_R1.fastq SECOND_END_FASTQ=${n}-mapped_R2.fastq
#rm $n.all.bam
#abyss-pe name=${n}_abyss k=64 in="${n}-mapped_R1.fastq ${n}-mapped_R1.fastq" && blast-contigs.sh ./*-3.fa
#wait
samtools view -bh -F4 -T $ref  $n.mappedReads.sam > $n.raw.bam
#samtools sort $n.raw.bam -o $n.sorted
#samtools index $n.sorted.bam

rm *sam

#echo "***Making Bam file"
#samtools view -bh -F4 -T $ref $n.sam > $n.raw.bam

echo "***Sorting Bam"
samtools sort $n.raw.bam -o $n.sorted.bam
echo "***Indexing Bam"
samtools index $n.sorted.bam

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picard} MarkDuplicates INPUT=$n.sorted.bam OUTPUT=$n.dup.bam METRICS_FILE=$n.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $n.dup.bam"
samtools index $n.dup.bam

#Number of nucleotides in reference with coverage
echo "*** Bamtools is getting coverage..."
bamtools coverage -in ${n}.sorted.bam | awk -v x=${n} 'BEGIN{OFS="\t"}{print x, $2, $3}' >> ${n}-coveragefile

#Mean depth of coverage
meancov=`awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' ${n}-coveragefile`

#Length of reference
countNTs=`awk 'END{print NR}' ${n}-coveragefile`

#count positions with coverage
covCount=`awk '{ if ($3 != 0) count++ } END { print count }' ${n}-coveragefile`

declare -i x=${covCount}
declare -i y=${countNTs}

#Percent of reference with coverage
perc=`awk -v x=$x -v y=$y 'BEGIN { print(x/y)*100}'`
echo "percent of genome with coverage: $perc"
echo "average coverage depth: $meancov"

echo "$n" > stats.txt
echo "Percent of genome with coverage: $perc" >> stats.txt
echo "Average coverage depth: $meancov" >> stats.txt
echo "" >> stats.txt

printf "%-20s %11.2f%% %'10dX\n" ${n} $perc $meancov >> ${n}.findbest

java -Xmx4g -jar ${gatk} -R $ref -T UnifiedGenotyper -I $n.sorted.bam -o ${n}.UG.vcf -nct 8

# make reference guided contig
java -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${n}.readreference.fasta -V ${n}.UG.vcf

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
