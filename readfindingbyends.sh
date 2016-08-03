#!/bin/sh

#Usage: ~$ readfindingbyends.sh <fastq_file_in> <fasta_seed_file_in> <kmer_length_(optional)>

if [[ -z $1 || -z $2 ]]; then
    echo "Usage: ~$ readfindingbyends.sh <fastq_file_in> <fasta_seed_file_in> <kmer_length_(optional)>"
    exit
fi

root=$(pwd)

#for debug
alias pause='read -p "$LINENO Enter"'

# set kmer_value to 20 if arg3 not given
new_kmer=$3
kmer_value=${new_kmer:-20}
echo "kmer_value: $kmer_value"

#put fasta and fastq file to variable
fasta_file=${root}/$2
fastq_file=${root}/$1

#specify header in fastq to avoid conflict
head_descript=$(grep -m 1 "^@" $fastq_file | sed 's/\(.....\).*/\1/')

endgrabber () {
# Create "here-document"

cat >./endtemp.py <<EOL
#!/usr/bin/python

import Bio
from Bio import SeqIO
from sys import argv 
          
script, prekmer, filename = argv

kmer = int(prekmer) + 1
count = 0
for record in SeqIO.parse(filename, "fasta"):
    myseq=record.seq
    print("%s" % myseq[1:int(kmer)])
    print("%s" % myseq[1:int(kmer)].reverse_complement())
    print("%s" % myseq[-int(kmer):-1]) 
    print("%s" % myseq[-int(kmer):-1].reverse_complement())   

EOL

chmod 755 ./endtemp.py

# $1 must be an interger, see script usage
# default kmer_value = 20
./endtemp.py $kmer_value $fasta_file

rm endtemp.py

}
#####

mkdir iteration0
cd iteration0
count=1
while [ $count -lt 40 ]; do 
    echo "kmer_value: $kmer_value"
    echo "fasta_file: $fasta_file"

    endgrabber > kmerinputlist

    # get headers that include the kmers
    knumber=$(wc -l kmerinputlist | awk '{print $1}')
    echo "kmer number $knumber: "
    echo "head_descript: $head_descript"

    while read l; do 
        grep -B 1 "$l" $fastq_file
    done < kmerinputlist | grep -v '^--$' | grep "${head_descript}" | sort | uniq > headers

    #merge headers with previous header
    #get a full list of headers
    if [ $count -lt 2 ]; then 
        echo "First iteration"

        cat headers | sort | uniq > newheaderlist
    else
        echo "current iteration $count"
        echo "using past iteration headers from ${past_iteration}/headers"

        cp ${past_iteration}/headers ./oldheaders
        cat headers oldheaders | sort | uniq > newheaderlist
    fi

    while read l; do 
        grep -A 3 "$l" $fastq_file
    done < newheaderlist > newreads.fastq

    reads_found=$(grep -c '^+$' newreads.fastq)
    
    echo "SPAdes is running"
    spades.py -k 55 -s newreads.fastq -o spades_output 2> /dev/null

    cp spades_output/scaffolds.fasta ./
    scaffold_number=$(grep -c '^>' scaffolds.fasta)
    scaffold_size=$(ls -lh scaffolds.fasta | awk '{print $5}')
    dir=$(pwd)
    fasta_file=${dir}/scaffolds.fasta
   
    echo "Interation${count} used $knumber kmers, found $reads_found reads, which assembled $scaffold_number scaffolds with a total file size of $scaffold_size" 
    echo "Interation${count} used $knumber kmers, found $reads_found reads, which assembled $scaffold_number scaffolds with a total file size of $scaffold_size" >> ${root}/stats

    past_iteration=$(pwd)
    count=$[$count+1]
    mkdir ${root}/iteration${count}
    cd ${root}/iteration${count}

done








# stuber 2016-04-20
