#!/bin/sh

###############################################################
alias pause='read -p "$LINENO Enter"'
mystart=$(date)

root=$(pwd)
prog="/usr/local/bin"

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

trim_reads () {

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

cp $forReads $revReads ./zips

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
ref="$prog"/bbmap/resources/truseq.fa.gz \
ktrim=r k=23 mink=11 hdist=1 tbo tbe \
qtrim=lr trimq=5 \
minlen=64 \
out1=trimmed_reads/${strain}_Trimmed_R1.fastq.gz \
out2=trimmed_reads/${strain}_Trimmed_R2.fastq.gz \
pigz=t \
unpigz=t 2>&1 | tee ${root}/${strain}/logs/logs.txt

rm ${forReads} ${revReads}
end=$(date +%s)
elapsed=$(($end - $start))
printf 'Trimming finished in %dh:%dm:%ds\n' \
$(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a ${root}/${strain}/logs/logs.txt

}

counter=1
cd $root
for f in */; do
    cd $root
    cd $f
    echo $f
    echo "*** File ${f}, $counter of $filecount"
    trim_reads
    ((counter++))
done


# tstuber 2016-05-05
