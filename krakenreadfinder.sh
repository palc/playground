#!/bin/sh

# Script to be used to isolate reads in fastq that have been identified using Kraken

# Set $1 as search term
# Working directory must have kraken report and output file, and fastq
# Fastqs can be zipped, single read or paired

# Usage example:  krakenreadfinder.sh brucell

# Unzip files if needed, else put std error to null
# find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip 2> /dev/null

# If needing to run Kraken...
# Paired reads
# find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P 60 gunzip; krakenDatabase="/home/shared/databases/kraken/std/"; sampleName=`ls *.fastq | head -1 | sed 's/_.*//' | sed 's/\..*//'`; kraken --db ${krakenDatabase} --threads 60 --paired *fastq* > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-report.txt

# Single reads
# krakenDatabase="/home/shared/databases/kraken/std/"; sampleName=`ls *.fastq | head -1 | sed 's/_.*//' | sed 's/\..*//'`; kraken --db ${krakenDatabase} --threads 60  *fastq* > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-report.txt

# Set flags
# flag -s Search term
# flag -n Number range
# flag -a run SPAdes
# flag -b run BLAST

NR_CPUS=$(($(nproc) - 5))

alias pause='read -p "$LINENO Enter"'

sflag=
nflag=
eflag=
aflag=
while getopts 'snab' OPTION; do
    case $OPTION in
        s) sflag=1
        ;;
        n) nflag=1
        ;;
        a) aflag=1
        ;;
        b) bflag=1
        ;;
        ?) echo "Invalid option: -$OPTARG" >&2
        ;;
    esac
done
shift $(($OPTIND - 1))


if [ "$sflag" ]; then
	# check that search argument was provided
	searcharg="$@"	
	echo ""
	echo "Argument given:  $searcharg"
elif [ "$nflag" ]; then 
	number1="$1"
	number2="$2"
	echo "The rows taxon numbers will be taken from in Kraken report"
	echo "$number1 - $number2"
else
	echo ""
	echo "Supply a flag and argument"
	echo ""
	echo "-n flag, Number usage: krakenreadfinder.sh -n 1180 1188"
	echo "-n is line numbers, not taxon numbers"
	echo "so to find all unclassified reads, krakenreadfinder.sh -n 1 1"
	echo "-s flag, Search usage: krakenreadfinder.sh -s brucell"
        echo "-a flag, run Assembly using SPAdes"
	echo "-b flag, run BLAST search on SPAdes output, -a must be called with -b"
	exit 1

fi
fastqfound=`ls *fastq* | head -1`

# check that needed files are in working directory
if [ -f *report*.txt -a -f *output*.txt -a $fastqfound ]; then
	echo ""
	echo "All files needed are present"
	echo ""
else
	echo ""
	echo "Check files in working directory"
	echo "Missing Kraken report, Kraken output or fastq file"
	echo "Looking for *-report*.txt, *-output*.txt and *fastq in working directory"
	echo ""
	exit 1
fi

sampleName=`ls *report*.txt | sed 's/[_-]report*.txt//'`

if [ "$sflag" ]; then
	echo "FINDINGS..."
	grep -i "$searcharg" *report*.txt | awk '{print $5}' > taxon
	cat taxon
elif [ "$nflag" ]; then
	echo "FINDINGS..."
	awk -v x=$number1 -v y=$number2 'NR==x,NR==y {print $5}' *report*.txt > taxon
	cat taxon
else
 	echo "provide flage"
	echo "ERROR"
	echo ""
	exit 1	
fi



# get taxon from search term provided as arguement
#grep -i "$searcharg" *-report.txt | awk '{print $5}' > taxon
idsfound=`grep -c ".*" taxon`
echo "The number of uniq taxon ids found using $searcharg:  $idsfound"
echo ""
echo "Getting read headers..."
echo ""

# get fastq headers from taxon identification
for i in `cat taxon`; do awk -v id=$i '$3 == id {print $2}' *-output*.txt; done >> headertemp

# REVERSE SEARCH, MAKE NEW FASTQ FOR ANYTHING THAT DOES NOT MATCH
#for i in `cat taxon`; do awk -v id=$i '$3 != id {print $2}' *-output.txt; done >> headertemp

# verify no duplication
sort < headertemp | uniq > headers
rm headertemp

headersfound=`grep -c ".*" headers`
echo "Fastq reads found using \"${searcharg}\" as search term:  $headersfound"
echo ""

echo "Unzipping read files..."
echo ""

# Unzip files if needed, else put std error to null
mkdir original_zips; cp *gz original_zips
forReads=`ls *_R1*`
revReads=`ls *_R2*`

find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip 2> /dev/null

echo "Gathering reads..."
echo ""

# isolate reads from original fastqs using identified headers
if [ -f *_R2* ]; then
    grep -F -A3 -h -f ./headers *_R1*.fastq | grep -v '^--$' > ${sampleName}-${searcharg}Reads_R1.fastq
    grep -F -A3 -h -f ./headers *_R2*.fastq | grep -v '^--$' > ${sampleName}-${searcharg}Reads_R2.fastq
else
    grep -F -A3 -h -f ./headers *fastq | grep -v '^--$' > ${sampleName}-${searcharg}Reads.fastq    
fi

rm ${forReads%.gz}
rm ${revReads%.gz}
rm headers
rm taxon

if [ "$aflag" ]; then
	echo "Assemblying isolated reads using SPAdes"
	spades.py -t 32 -k 127 --careful -1 ${sampleName}-${searcharg}Reads_R1.fastq -2 ${sampleName}-${searcharg}Reads_R2.fastq -o spades_output
else
	echo "No assembly done"
fi

if [ "$bflag" ]; then
	file="./spades_output/scaffolds.fasta"
	echo "BLAST scaffolds.fasta"
    echo "$file was reformated placing sequence on a sigle line"
        awk '{if ($0 ~ />/ ) {print ">"} else print $0}' $file | tr -d "\n" | sed -e 's:>:\n>\n:g' | awk '{if ($0 ~  />/ ) {print ">contigs-" x++} else print $0}' | grep -v "^$" > ${sampleName}-contigsoriginal.fa
	ls -lh ${sampleName}-contigsoriginal.fa
    blast-contigs.sh ${sampleName}-contigsoriginal.fa
	cat BLAST-summ*
    mkdir blast_results; mv BLAST-* blast_results
    mv ${sampleName}-contigsoriginal.fa blast_results
    rm *xml
    pigz ${sampleName}-Reads_R?.fastq
else
	echo "No BLAST done"
fi	

echo ""
echo "FINISHED"
echo ""

# If few reads transform fastq to fasta and identify with:  blast-contig.sh outfile.fasta
# awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' infile.fastq > outfile.fasta
# Fastq to fasta convert

# If many reads assemble and blast using:  assemble_verify.sh

# Created 2015-04-14, Stuber

