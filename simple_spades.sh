#!/bin/sh

starttime=$(date +%s)
alias pause='read -p "$LINENO Enter"'

# unzip a single FASTQ if needed
if [[ -e $(ls *fastq | head -1) ]]; then
    in_fastq=`ls *fastq | head -1`
elif [[ -e $(ls *gz | head -1) ]]; then
    pigz -d $(ls *gz | head -1)
    in_fastq=`ls *fastq | head -1`
else
    printf "### ERROR no FASTQs available\n\n"
    exit 1
fi

####################### 
# perl script to choose optimal kmer
cat >./get_size.pl <<'EOL'
#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;

# Output kmer file
open (my $outkmer, '>', "kmer.txt") or die "$!";

my $infile = $ARGV[0];
my $ave_length;

my $inseq;
my $frag_size_total=0;
my $counter=0;
$inseq = Bio::SeqIO->new(-file => $infile, -format => "fastq");
while ((my $seq_obj = $inseq->next_seq) && ($counter <= 100)) {
    $counter++;
    my $length = $seq_obj->length;
    $frag_size_total = $frag_size_total + $length;
}
$ave_length = $frag_size_total / 100;
print "$ave_length\n";

my $kmer = 0;
if ($ave_length > 130) {
    $kmer = 127
} else {
    $ave_length=(int $ave_length);
    if ($ave_length % 2 == 0) {
        print "average is even number\n";
        $kmer = $ave_length - 7;
    }  else {
        print "average is odd number\n";
        $kmer = $ave_length - 8;
    }
}

print $outkmer "kmer: $kmer";
EOL

chmod 755 ./get_size.pl
#######################
./get_size.pl ${in_fastq}
kmer=`awk '{print $2}' kmer.txt`
rm get_size.pl kmer.txt
printf "kmer used: $kmer\n"

# zip file if needed
if [[ -e $(ls *fastq | head -1) ]]; then
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

# get sample name
strain=$(echo *fastq* | head -1 | sed 's/[._].*//')

# trim reads
echo "Trimming reads"
if [ $filecount -eq 2 ]; then 

    read1=$(ls *gz | head -1)
    read2=$(ls *gz | tail -1)
    echo "Forward Reads to be trimmed: $read1"
    echo "Reverse Reads to be trimmed: $read2"

    echo "Quality trimming sample $strain"

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
    
read1=`ls | grep _R1`
    echo "Forward Reads to be used after trimmed: $read1"
    read2=`ls | grep _R2`
    echo "Reverse Reads to be used after trimmed:: $read2"

elif [ $filecount -eq 1 ]; then

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
    read1=${strain}_Trimmed.fastq.gz
    echo "Forward Reads to be used after trimmed: $read1"
fi

# assemble trimmed reads
printf "SPAdes running...\n"
if [[ $readtype == paired ]]; then
    spades.py -k ${kmer} --careful -1 ${read1} -2 ${read2} -o spades_output &> /dev/null
    file="./spades_output/scaffolds.fasta"
    rm ${read1} ${read2}
else
    spades.py -k ${kmer} -s ${read1} -o spades_output &> /dev/null
    file="./spades_output/scaffolds.fasta"
    rm ${read1}
fi

# error if scaffolds.fasta not made
if [[ ! -s $file ]]; then
    printf "*** error\n*** error\n*** error\n"
    printf "A proper assembled file did not complete\n"
    printf "Exiting script \n\n"
    exit 1
fi

# save needed files
cp spades_output/scaffolds.fasta ${strain}.scaffolds.fasta
cp spades_output/params.txt ./
cp spades_output/spades.log ./
cp spades_output/warnings.log ./
rm -r spades_output/

endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))
printf "***DONE\n"
# created 2015-10-26, Tod Stuber
