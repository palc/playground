#!/bin/sh

starttime=$(date +%s)
alias pause='read -p "$LINENO Enter"'

if [[ -e *fastq ]]; then
    pigz *fastq
fi

# check number of fastq.gz
# 1 == single read, 2 == paired, >2 == error
filecount=`ls *gz | wc -l`
if [ $filecount -gt 2 ]; then 
    printf ">2 fastq.gz in working directory\n"
    printf "Anymore than 2 reads script exits\n\n"
    exit 1
elif [ $filecount -eq 2 ]; then
    printf "paired reads found\n"
    readtype=paired
else
    printf "single read found\n"
    readtype=single
fi

#save original reads
mkdir original_reads
mv *gz original_reads
# Make alias in working directory to zip files
ln -s original_reads/*gz ./

strain=$(echo *fastq* | head -1 | sed 's/_.*//' | sed 's/\..*//')

pause
echo "Trimming reads"
if [ $filecount -eq 2 ]; then 
    read1=$(ls *gz | head -1)
    read2=$(ls *gz | tail -1)

    read1=`ls | grep _R1`
    echo "Forward Reads to be trimmed: $read1"
    read2=`ls | grep _R2`
    echo "Reverse Reads to be trimmed:: $read2"
    
    #Trim the reads with bbmap tool kit (bbduk plugin)
    #about twice as fast as trimmomatic

    echo -e "Quality trimming sample "$strain""

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
    
    rm "$read1"
    rm "$read2"
    
    read1=`ls | grep _R1`
    echo "Forward Reads to be used after trimmed: $read1"
    read2=`ls | grep _R2`
    echo "Reverse Reads to be used after trimmed:: $read2"

elif [ $filecount -eq 1 ]; then
    #Trim the reads with bbmap tool kit (bbduk plugin)
    #about twice as fast as trimmomatic

    printf "Quality trimming sample $strain\n"
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
    read1=`ls | grep _R1`
    echo "Forward Reads to be used after trimmed: $read1"
fi
pause

# assemble trimmed reads
printf "SPAdes running...\n"
if [[ $readtype == paired ]]; then
    spades.py -k ${kmer} --careful -1 ${read1} -2 ${read2} -o spades_output &> /dev/null
    file="./spades_output/scaffolds.fasta"
else
    spades.py -k ${kmer} -s ${read1} -o spades_output &> /dev/null
    file="./spades_output/scaffolds.fasta"
fi

if [[ ! -s $file ]]; then
    printf "*** error\n*** error\n*** error\n"
    printf "A proper assembled file did not complete\n"
    printf "Exiting script \n\n"
    exit 1
fi

endtime=`date +%s`
runtime=$((endtime-starttime))
#totaltime=`date -u -d @${runtime} +"%T"`
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) >> sectiontime
print "***DONE\n"
# created 2015-10-26, Tod Stuber
