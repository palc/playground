#!/bin/sh

starttime=$(date +%s)
alias pause='read -p "$LINENO Enter"'

NR_CPUS=40

help () {
    printf "\n\n"
    printf "### Working directory containing FASTQ/s and Reference with .fasta extension\n"
    printf "Usage: $$ simple_align.sh\n"
    exit 1
}

# zip FASTQs
if [[ -e $(ls *fastq | head -1) ]]; then
    pigz *fastq
fi

ref=`ls | grep .fasta`
echo "Reference Input:  $ref"
if [[ ! -s $ref ]]; then
    printf "\n### No reference provided with .fasta extension\n\n"
    help
fi
r=`echo $ref | sed 's/\..*//'`

# check number of fastq.gz
# 1 == single read, 2 == paired, >2 == error
filecount=`ls *gz | wc -l`
if [ $filecount -gt 2 ]; then 
    help
elif [ $filecount -eq 2 ]; then
    printf "paired reads found\n"
    readtype=paired
else
    printf "single read found\n"
    readtype=single
fi

# save original reads
mkdir original_reads
mv *gz original_reads
# Make alias in working directory to zip files
ln -s original_reads/*gz ./

# get sample name
strain=$(echo *fastq.gz | head -1 | sed 's/[._].*//')
printf "Strain name $strain\n\n"
echo "Quality trimming sample $strain"

# trim reads
echo "Trimming reads"
if [ $filecount -eq 2 ]; then
    read1=$(ls *gz | head -1)
    read2=$(ls *gz | tail -1)
    echo "Forward Reads to be trimmed: $read1"
    echo "Reverse Reads to be trimmed: $read2"

    bbduk.sh -Xmx80g \
        in1="$read1" \
        in2="$read2" \
        ref="/usr/local/bin/bbmap/resources/nextera.fa.gz" \
        ktrim=r k=23 mink=11 hdist=1 \
        qtrim=lr trimq=5 \
        minlen=36 \
        out1=${strain}_Trimmed_R1.fastq.gz \
        out2=${strain}_Trimmed_R2.fastq.gz \
        stats=trim_stats.txt \
        qchist=qc_by_base.txt \
        threads=auto \
        showspeed=f
    rm ${read1} ${read2}
    read1=$(ls *gz | head -1)
    echo "Forward Reads to be used after trimmed: $read1"
    read2=$(ls *gz | tail -1)
    echo "Reverse Reads to be used after trimmed:: $read2"
elif [ $filecount -eq 1 ]; then
    read1=`ls *fastq*`

    bbduk.sh -Xmx80g \
        in=$read1 \
        ref="/usr/local/bin/bbmap/resources/nextera_flu.fa.gz" \
        trimq=20 \
        minlength=75 \
        minavgquality=20 \
        chastityfilter=t \
        maxns=2 \
        out=${strain}_Trimmed.fastq.gz \
        overwrite=t \
        stats=stats.txt \
        qchist=qc_by_base.txt \
        threads=auto \
        showspeed=f &>> trimming_stats.txt

    rm $read1
    read1=${strain}_Trimmed.fastq.gz
    echo "Forward Reads to be used after trimmed: $read1"
fi

### alignment
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

echo "***bwa index $r"
bwa index $ref

if [ $filecount -eq 2 ]; then
echo "***Making sam file from paired reads"
    #adding -B 8 will require reads to have few mismatches to align to reference.  -B 1 will allow more mismatch per read.
    bwa mem -M -B 8 -t 16 -R @RG"\t"ID:"$strain""\t"PL:ILLUMINA"\t"PU:"$strain"_RG1_UNIT1"\t"LB:"$strain"_LIB1"\t"SM:"$strain" $ref ${read1} ${read2} > $strain.sam
    rm ${read1} ${read2}
elif [ $filecount -eq 1 ]; then
    echo "***Making sam file from single read"
    bwa mem -M -B 1 -t 10 -T 20 -P -a -R @RG"\t"ID:"$strain""\t"PL:ILLUMINA"\t"PU:"$strain"_RG1_UNIT1"\t"LB:"$strain"_LIB1"\t"SM:"$strain" $ref ${read1} > $strain.sam
    rm ${read1}
fi

samtools view -bh -T $ref $strain.sam > $strain.all.bam
samtools view -h -F4 $strain.all.bam > $strain.mappedReads.sam
samtools view -bh -F4 -T $ref  $strain.mappedReads.sam > $strain.raw.bam
rm $strain.sam
echo "***Sorting Bam"
samtools sort $strain.raw.bam -o $strain.sorted.bam
echo "***Indexing Bam"
samtools index $strain.sorted.bam

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picard} MarkDuplicates INPUT=$strain.sorted.bam OUTPUT=$strain.dup.bam METRICS_FILE=$strain.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true
samtools view -h -F4 $strain.dup.bam > $strain.dedupmappedReads.sam
java -Xmx4g -jar ${picard} SamToFastq INPUT=$strain.dedupmappedReads.sam FASTQ=${strain}-dedup-filtered_R1.fastq SECOND_END_FASTQ=${strain}-dedup-filtered_R2.fastq

echo "***Index $strain.dup.bam"
samtools index $strain.dup.bam

#Number of nucleotides in reference with coverage
echo "*** Bamtools is getting coverage..."
bamtools coverage -in ${strain}.sorted.bam | awk -v x=${strain} 'BEGIN{OFS="\t"}{print x, $2, $3}' >> ${strain}-coveragefile

#Mean depth of coverage
meancov=`awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' ${strain}-coveragefile`

java -Xmx4g -jar ${gatk} -R $ref -T UnifiedGenotyper -I $strain.sorted.bam -o ${strain}.UG.vcf -nct 8

echo '##fileformat=VCFv4.1' > $strain-highqualitysnps.vcf
echo '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1 rsIDs' | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' >> $strain-highqualitysnps.vcf
grep -v '^#' ${strain}.UG.vcf | awk 'BEGIN{OFS="\t"} $6 > 500 {print $1, $2, ".", $4, ".",  ".", ".", ".", ".", ".", "."}' >> $strain-highqualitysnps.vcf

# file clean up
rm ${strain}.all.bam
rm ${strain}.dedupmappedReads.sam
rm ${strain}.dup.bam*
rm ${strain}.FilteredReads.xls
rm ${strain}.raw.bam
rm ${strain}.mappedReads.sam

printf "\nMean Coverage: $meancov\n\n"
endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))
printf "***DONE\n"
# created 2015-10-26, Tod Stuber
