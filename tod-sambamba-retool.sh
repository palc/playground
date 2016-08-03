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
  
bbmap.sh "$mem" \
    in1="trimmed_reads/${strain}_Trimmed_1P.fastq.gz" \
    in2="trimmed_reads/${strain}_Trimmed_2P.fastq.gz" \
    rgid="$strain" \
    rgcn="OLC" \
    rglb="NexteraXT" \
    rgpl="ILLUMINA" \
    rgsm="$strain" \
    outm=aligned/${strain}_mapped.bam \
    outu=unmapped/${strain}_unmapped.fastq.gz \
    ref="$genome" \
    #question
    path=$HOME/ref/NCBI_RefSeq_DEC232015 \
    pairedonly=t \
    pigz=t \
    unpigz=t 2>&1 | tee -a ${root}/${strain}/logs/logs.txt

#Sort and index VCF
sambamba sort -t "$cpu1" -o aligned/${strain}_mapped_sorted.bam aligned/${strain}_mapped.bam
 
#remove duplicates
sambamba markdup -r -t "$cpu1" aligned/${strain}_mapped_sorted.bam aligned/${strain}_mapped_sorted_nodup.bam 2>&1 | tee -a ${root}/${strain}/logs/logs.txt
 
#remove temporary files
rm -r trimmed_reads
rm aligned/${strain}_mapped.bam*
rm aligned/${strain}_mapped_sorted.bam*

#Create dictionary file for the reference genome (needed by gatk)
if [ -e "${genome%.*}".dict ]; then #check if indexing already done
  echo -e "Reference genome $(basename "$genome") already indexed. Skipping this step."
else
  echo -e "Indexing reference genome $(basename "$genome") with Picard tools"
  java "$mem" -jar "$picard" CreateSequenceDictionary \
  R="$genome" \
  O="${genome%.*}".dict #have to strip oof the ".fasta" extension and replace it by ".dict"
fi
}

variant_calling () {

cd aligned
echo "Variant Calling"
strain=$(basename *_mapped_sorted_nodup.bam | sed 's/_.*//')
echo "$strain"
pwd

#######################
#                     #
#   Variant Calling   #
#                     #
#######################
 

#Collect Depth of coverage info
echo "***Collect Depth of Coverage"
	java -jar ${gatk} -T DepthOfCoverage -R $genome -I ${strain}_mapped_sorted_nodup.bam -o ${strain}.coverage -omitIntervals --omitLocusTable --omitPerSampleStats -nt 8

echo ${strain}

echo "***HaplotypeCaller, aka calling SNPs"
#-allowNonUniqueKmersInRef
java -Xmx4g -jar ${gatk} -R $genome -T HaplotypeCaller -I ${strain}_mapped_sorted_nodup.bam -o ${strain}.hapreadyAll.vcf -bamout ${strain}.bamout.bam -dontUseSoftClippedBases -allowNonUniqueKmersInRef
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

counter=1

cd $root
for f in */; do 
    cd $root
    cd $f
    baseDir=$(pwd)
    echo $f
    echo "*** File ${f}, $counter of $filecount"
    initial_alignment
    ((counter++))
done

echo "     *********************************"
echo "     ***                           ***"
echo "     ***        bam_preproc        ***"
echo "     ***                           ***"
echo "     *********************************"
cd $root

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
    variant_calling) &
let count+=1
[[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

wait

myend=$(date)
echo "started: $mystart"
echo "ended: $myend"

###############################################################
