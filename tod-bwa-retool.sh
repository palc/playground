#!/bin/sh

###############################################################
alias pause='read -p "$LINENO Enter"'
mystart=$(date)

root=$(pwd)
prog="/usr/local/bin"
genome="/home/shared/mycobacterium/tbc/snppipeline/tbbov/NC_002945.fasta"
dbSNP="/home/shared/mycobacterium/tbc/snppipeline/tbbov/HighestQualitySNPs.vcf.gz"

#Number of real cpu (not number from hyper-threading).

NR_CPUS=6
cpu1=$(($(nproc)-2))
cpu2=$(($(nproc)/NR_CPUS))
mem="-Xmx80g"

for i in *.fastq*; do 
    n=`echo $i | sed 's/[._].*//'`
    echo "n is : $n"
    mkdir -p $n
    mv $i $n/
done
filecount=$(ls -d */ | wc -l)

initial_alignment () {

#Create required directory structure

echo "*** Data Structure"
#######################
#                     #
#   Data Structure    #
#                     #
#######################
 
mkdir logs
mkdir zips
mkdir trimmed_reads
mkdir aligned
mkdir unmapped
 
cp *gz zips 

echo "*** Trimming"
#######################
#                     #
#      Trimming       #
#                     #
#######################
 
forReads=`ls | grep _R1`
echo "Forward Reads: $forReads"

revReads=`ls | grep _R2`
echo "Reverse Reads: $revReads"

#sequence nomenclature:
# 2014-SEQ-0729_S5_L001_R1_001.fastq.gz
 
#Trim the reads with bbmap tool kit (bbduk plugin)
#about twice as fast as trimmomatic
start=$(date +%s)
strain=$(echo $revReads | sed 's/_.*//' | sed 's/\..*//') 
echo -e "Quality trimming sample "$strain""
  
    bbduk.sh "$mem" \
    in1="$forReads" \
    in2="$revReads" \
    ref="$prog"/bbmap/resources/nextera.fa.gz \
    ktrim=r k=23 mink=11 hdist=1 tbo tbe \
    qtrim=lr trimq=5 \
    minlen=64 \
    out1=trimmed_reads/${strain}_Trimmed_1P.fastq.gz \
    out2=trimmed_reads/${strain}_Trimmed_2P.fastq.gz \
    pigz=t \
    unpigz=t 2>&1 | tee ${root}/${strain}/logs/logs.txt 

rm ${forReads} ${revReads}

echo "*** Reference Genome"
#######################
#                     #
#  Reference Genome   #
#                     #
#######################

cd ${root}/${strain}/aligned
cp $genome ./
ref=$(ls | grep .fasta)
echo "Reference Input:  $ref"
#Retrieves reference name and name from sorted BAM file name
r=$(echo $ref | sed 's/\..*//')
n=${strain}
echo "***Reference naming convention:  $r"
echo "***Isolate naming convention:  $n"

#ref genome indexing
echo "***Reference Genome $r is being indexed"
samtools faidx $ref
java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict


echo "***Alignment"
#######################
#                     #
#      Alignment      #
#                     #
#######################

#For each sample:
# 2014-SEQ-0729_Trimmed_1P.fq.gz - for paired forward reads
# 2014-SEQ-0729_Trimmed_2P.fq.gz - for paired reverse reads
start=$(date +%s)
  
echo -e "Aligning $strain sequencing data on reference genome $(basename "$genome")"

# See echo comments
echo "***bwa index $r"
bwa index $ref

# -t sets the number of threads/cores
# -r ST	 Specify the read group in a format like ‘@RG\tID:foo\tSM:bar’ Needed for GATK
#adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
echo "***Making Sam file"
bwa mem -M -t 16 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref ${root}/${strain}/trimmed_reads/${strain}_Trimmed_1P.fastq.gz ${root}/${strain}/trimmed_reads/${strain}_Trimmed_2P.fastq.gz > $n.sam

echo "***Making Bam file"
samtools view -bh -F4 -T $ref $n.sam > $n.raw.bam

echo "***Sorting Bam"
samtools sort $n.raw.bam $n.sorted
echo "***Indexing Bam"
samtools index $n.sorted.bam
# Remove duplicate molecules

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picard} MarkDuplicates INPUT=$n.sorted.bam OUTPUT=$n.dup.bam METRICS_FILE=$n.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index $n.dup.bam"
samtools index $n.dup.bam

#remove temporary files
rm $n.sam*
rm $n.raw.bam*
rm $n.sorted.bam*

end=$(date +%s)
elapsed=$(($end - $start))
#log
printf "Trimmed reads aligned to $(basename "$genome") in %dh:%dm:%ds\n" \
$(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a ${root}/${strain}/logs/logs.txt

}

variant_calling () {

cd ${root}/${strain}/aligned
echo "Variant Calling"
strain=$(basename *dup.bam | sed 's/[_.].*//')
echo "$strain"
pwd

#######################
#                     #
#   Variant Calling   #
#                     #
#######################

#Collect Depth of coverage info
echo "***Collect Depth of Coverage"
	java -jar ${gatk} -T DepthOfCoverage -R $genome -I $n.dup.bam -o ${strain}.coverage -omitIntervals --omitLocusTable --omitPerSampleStats -nt 8

echo ${strain}

echo "***HaplotypeCaller, aka calling SNPs"
#-allowNonUniqueKmersInRef
java -Xmx4g -jar ${gatk} -R $genome -T HaplotypeCaller -I $n.dup.bam -o ${strain}.hapreadyAll.vcf -bamout ${strain}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef
java -Xmx4g -jar ${igvtools} index ${strain}.hapreadyAll.vcf

echo "******Awk VCF leaving just SNPs******"
awk '/#/ || $4 ~ /^[ATGC]$/ && $5 ~ /^[ATGC]$/ {print $0}' ${strain}.hapreadyAll.vcf > ${strain}.hapreadyOnlySNPs.vcf

#Split header lines from position calls
grep "#" ${strain}.hapreadyOnlySNPs.vcf > ${strain}.header
grep -v "#" ${strain}.hapreadyOnlySNPs.vcf > ${strain}.body

#SNP positons that will be used
awk '{print $1 "%" $2}' ${strain}.body > ${strain}.calledSNPpositions

#Zero coverage positions
awk 'BEGIN {FS="[:\t]"} $3 == 0 {print $1 "%" $2}' ${strain}.coverage > ${strain}.zeroCoveragePositions
#Remove zero coverage positions that will be are in "${strain}".hapreadyOnlySNPs.vcf
cat ${strain}.calledSNPpositions ${strain}.zeroCoveragePositions | sort | uniq -d > ${strain}.duplicates
cat ${strain}.zeroCoveragePositions ${strain}.duplicates | sort | uniq -u > ${strain}.keepTheseZeroCovPositions
zeroposition=`grep -c ".*" ${strain}.keepTheseZeroCovPositions`
refsize=`wc -m $genome | awk '{print $1}'`

#Fromat "${strain}".keepTheseZeroCovPositions to VCF
sed 's/%/ /' ${strain}.keepTheseZeroCovPositions | awk 'BEGIN{OFS="\t"}{print $1, $2, ".", ".", ".", ".", ".", ".", "GT", "./."}' > ${strain}.vcfFormated
cat ${strain}.body ${strain}.vcfFormated | awk 'BEGIN{OFS="\t"}{if ($4 == ".") print $1, $2, $3, "N", $5, $6, $7, $8, $9, $10; else print $0}' > ${strain}.SNPsMapzeroNoHeader.vcf
cat ${strain}.header ${strain}.SNPsMapzeroNoHeader.vcf > ${strain}.unsortSNPsZeroCoverage.vcf
java -Xmx4g -jar ${igvtools} sort ${strain}.unsortSNPsZeroCoverage.vcf ${strain}.SNPsZeroCoverage.vcf
java -Xmx4g -jar ${igvtools} index ${strain}.SNPsZeroCoverage.vcf


rm ${strain}.duplicates
rm ${strain}.body
rm ${strain}.calledSNPpositions
rm ${strain}.allsites.vcf*
rm ${strain}.hapreadyOnlySNPs.vcf
rm ${strain}.header
rm ${strain}.keepTheseZeroCovPositions
rm ${strain}.SNPsMapzeroNoHeader.vcf
rm ${strain}.unsortSNPsZeroCoverage.vcf
rm ${strain}.vcfFormated
rm ${strain}.zeroCoveragePositions
rm igv.log
cp $genome ./

end=$(date +%s)
elapsed=$(($end - $start))
printf 'Variant calling done in %dh:%dm:%ds\n' $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a ${root}/${strain}/logs/logs.txt

}

#counter=1
#
#cd $root
#for f in */; do 
#    cd $root
#    cd $f
#    baseDir=$(pwd)
#    echo $f
#    echo "*** File ${f}, $counter of $filecount"
#    initial_alignment
#    ((counter++))
#done

wait

echo "     *********************************"
echo "     ***                           ***"
echo "     ***      variant_calling      ***"
echo "     ***                           ***"
echo "     *********************************"

NR_CPUS=12
cd $root
for f in */; do
    (cd $root
    cd $f
    echo $f
    initial_alignment
    variant_calling) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

wait

myend=$(date)
echo "started: $mystart"
echo "ended: $myend"

###############################################################
