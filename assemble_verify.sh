#!/bin/sh

starttime=`date +%s`

alias pause='read -p "$LINENO Enter"'

help () {

printf "\n\n\n"
printf "Input must be single or paired FASTQs\n"
printf "Working directory must have single sample\n"

printf "\n -a assembler\n -k kmer to use in assembler\n -s BLAST search terms to verify contigs\n -m Kraken database\n -v no email\n\n"

printf 'Full usage: ~$ assembly_verify.sh -a abyss -k 61 -s "myco bruc" -m host\n'
printf "Assembler options: spades, abyss, ray\n"
printf "kmer defaults--> spades:127 abyss:64 ray:51\n\n"
printf "Search term (-s) usage:\n"
printf "    ~$ assembly_verify.sh -a Ray -s virus\n"
printf "    ...or for multiple search terms:\n"
printf '    Usage: ~$ assembly_verify.sh -a Ray -s "myco bruc virus"\n\n'
printf "FASTA file can be input.  Extention must be .fasta\n"
printf 'A search term is still needed. ~$ assembly_verify.sh -s "myco"\n\n'
printf "Available Kraken databases: jhu-flu, host, std\n"
printf "Defaults sending an email to tod.p.stuber@usda.gov\n\n\n"
exit 1
}

################################################################################################################################
#|||||||||||||||||||||||||||||||||||||||||| FUNCTION:  Create binned histo of contig sizes |||||||||||||||||||||||||||||||||||||
################################################################################################################################
function histo () {

# Create "here-document" to prevent a dependent file.
cat >./histo.py <<EOL
#!/usr/bin/env python

from sys import argv
from Bio import SeqIO
import pylab

script, input = argv

print (input)

sizes = [len(rec) for rec in SeqIO.parse(input, "fasta")]
len(sizes), min(sizes), max(sizes)

pylab.hist(sizes, bins=100)
pylab.title("%i $n sequences\nLengths %i to %i" \
            % (len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Sequence length (bp)\n100 bins along x-axis")
pylab.ylabel("Contig count")
pylab.savefig('./${n}_fig.pdf')

EOL

chmod 755 ./histo.py

./histo.py $file 1>&2 histo_std_out_error

#rm ./histo.py

}

################################################################################################################################
#|||||||||||||||||||||||||||||||||||||||||| FUNCTION:  Create binned histo of contig sizes |||||||||||||||||||||||||||||||||||||
################################################################################################################################
function contig_select () {

# Create "here-document" to prevent a dependent file.
cat >./contigs.py <<EOL
#!/usr/bin/env python

import Bio
from Bio import SeqIO
from sys import argv
script, input = argv

count_all=0
count_small=0
count_passing=0
passing_contigs = []

records = SeqIO.parse(input, "fasta")

for record in records:
    count_all += 1
    if len(record.seq) < 300:
        count_small += 1
    else:
        count_passing += 1
        passing_contigs.append(record)

print ("    < 300bp: %s" % (count_small))
print ('    => 300bp: %s' % (count_passing))

SeqIO.write(passing_contigs, "./$n-gt300_contigs-ver.fasta", "fasta")

EOL

chmod 755 ./contigs.py

./contigs.py verifiedreads-${n}.fasta 

rm ./contigs.py

}

################################################################################################################################

# System variable "NTDB" must be set to the path of the nt database
# System variable "NR_CPUS" must be set for Number of Cores to use

NR_CPUS=$(($(nproc) - 5))

# Set flags
hflag=
assembler=
searchterm=
kmer=
krakendb=
vflag=
while getopts ':ha:s:k:m:v' OPTION; do
    case $OPTION in
        h) hflag=1
        ;;
        a) assembler=$OPTARG
        ;;
        s) searchterm=$OPTARG
        ;;
        k) kmer=$OPTARG
        ;;
        m) krakendb=$OPTARG
        ;;
        v) vflag=1
        ;;
        ?) printf "\n\nIncorrect argument given"; help
        ;; 
    esac
done
shift $(($OPTIND - 1))

##

if [ "$hflag" ]; then
    help
    exit 1
fi

echo "Database: ${NTDB}"
echo "Number of cores being used: ${NR_CPUS}"
echo "Program used: $assembler"
echo "search/s term provide: $searchterm"
echo "kmer size for assembler: $kmer"

echo "Number of cores being used: ${NR_CPUS}"
log1="Program used: $assembler"
log2= "Search terms provide: $searchterm"


if [ -z $NTDB ]; then 
 	echo "System variable NTDB must be set to the path of the nt database"
	echo ""
	exit 1
fi

if [ -z $NR_CPUS ]; then 
 	echo "System variable NR_CPUS must be set for Number of Cores to use"
	echo ""
	exit 1
fi

# zip fastq/s if not already
if [[ -s $(ls *fastq | head -1) ]]; then 
    printf "Files being zipped\n\n"
    gzip *fastq
fi

ls -lh *gz | awk '{print $5, $9}' > fileinfo  

# check number of fastq.gz
# 1 == single read, 2 == paired, >2 == error
filecount=`ls *gz | wc -l`
if [ $filecount -gt 2 ]; then 
    printf ">2 fastq.gz in working directory\n"
    printf "Anymore than 2 reads script exits\n\n"
    exit 1 
fi

if [ $filecount -eq 1 ]; then
    read1=$(ls *gz)
    readtype=single
    
    #save original reads
    mkdir original_reads
    mv *gz original_reads
    # Make alias in working directory to zip files
    ls ../zips/*.fastq* | while read file; do ln -s $file; done
    printf  "Single end read: %s \n\n" $read1
fi

if [ $filecount -eq 2 ]; then 
    read1=$( ls *gz | head -1)
    read2=$(ls *gz | tail -1)
    readtype=paired
    
    #save original reads
    mkdir original_reads
    cp *gz original_reads
    # Make alias in working directory to zip files
    ls ../zips/*.fastq* | while read file; do ln -s $file; done
    printf "Paired reads: %s, %s \n\n" $read1 $read2
    echo "Trimming reads"
    echo "*** Trimming"
    #######################
    #                     #
    #      Trimming       #
    #                     #
    #######################
    
    read1=`ls | grep _R1`
    echo "Forward Reads to be trimmed: $read1"
    
    read2=`ls | grep _R2`
    echo "Reverse Reads to be trimmed:: $read2"
    
    #Trim the reads with bbmap tool kit (bbduk plugin)
    #about twice as fast as trimmomatic
    
    strain=$(echo $read2 | sed 's/_.*//' | sed 's/\..*//')
    echo -e "Quality trimming sample "$strain""
    
        bbduk.sh -Xmx80g \
        in1="$read1" \
        in2="$read2" \
        ref="/usr/local/bin/bbmap/resources/nextera.fa.gz" \
        ktrim=r k=23 mink=11 hdist=1 \
        qtrim=lr trimq=5 \
        minlen=36 \
        out1=trimmed_reads/${strain}_Trimmed_R1.fastq.gz \
        out2=trimmed_reads/${strain}_Trimmed_R2.fastq.gz \
        stats=trim_stats.txt \
        qchist=qc_by_base.txt \
        threads=auto \
        showspeed=f
    
    mv -v trimmed_reads/${strain}_Trimmed_R1.fastq.gz ./
    mv -v trimmed_reads/${strain}_Trimmed_R2.fastq.gz ./
    rm -r trimmed_reads
    rm "$read1"
    rm "$read2"
    
    read1=`ls | grep _R1`
    echo "Forward Reads to be used after trimmed: $read1"
    
    read2=`ls | grep _R2`
    echo "Reverse Reads to be used after trimmed:: $read2"
fi

log4="$readtype read/s"

# Get name of isolate to use for naming output files
n=`echo $read1 | sed 's/[_\.].*//'`
echo "***Isolate naming convention:  $n"

if [[ -z $assembler ]]; then
    printf "\n\n Assembler and BLAST not ran \n\n"
    log9=`printf "\n\n Assembler \n"` 
    file=$(ls *fasta | head -1)
    n=`echo $file | sed 's/[_\.].*//'`
    echo "***Isolate naming convention:  $n"
else    
    echo "$LINENO"
    # Check for dependent programs
    check1=`which spades.py`
    check2=`which blastn`
    check3=`which Ray`
    check4=`which abyss-pe`
    
    
    if [ -z $check1 ]; then 
    	echo "Program spades.py cannot be found"
    	echo 'Check program is installed and in $PATH'
    	exit 1
    fi
    
    if [ -z $check2 ]; then 
    	echo "Program blastn cannot be found"
    	echo 'Check program is installed and in $PATH'
    	exit 1
    fi
    
    if [ -z $check3 ]; then
        echo "Program Ray cannot be found"
        echo 'Check program is installed and in $PATH'
        exit 1
    fi
    
    if [ -z $check4 ]; then
        echo "Program abyss-pe cannot be found"
        echo 'Check program is installed and in $PATH'
        exit 1
    fi
    
    # Called assembler is ran
    # kmer is defaulted if not given as argument
    # organize assembly output
    # set "file" to be BLAST
    if [[ $assembler == spades ]]; then
        printf "`date` \nspades is running... \n\n" 
        if [[ -z $kmer ]] || [[ $kmer -gt 127 ]]; then
            kmer=127
        fi
        if [[ $readtype == paired ]]; then
            spades.py -k ${kmer} --careful -1 ${read1} -2 ${read2} -o spades_output &> /dev/null
            file="./spades_output/scaffolds.fasta"
        else
            spades.py -k ${kmer} -s ${read1} -o spades_output &> /dev/null
            file="./spades_output/scaffolds.fasta"
        fi
    fi
    
    if [[ $assembler == abyss ]]; then 
        printf "`date` \nabyss is running... \n\n"
        if [[ -z $kmer ]] || [[ $kmer -gt 64 ]]; then
            kmer=64
        fi
        if [[ $readtype == paired ]]; then
            abyss-pe name=Abyss k=${kmer} in="${read1} ${read2}" &> abysslog
            mkdir abyss_assembly
            mv *Abyss-* abysslog coverage.hist abyss_assembly
            file="./abyss_assembly/Abyss-8.fa"
        else
            abyss-pe name=Abyss k=${kmer} in="${read1}" &> abysslog
            mkdir abyss_assembly
            mv *Abyss-* abysslog coverage.hist abyss_assembly
            file="./abyss_assembly/Abyss-3.fa"   
        fi
    fi
    
    if [[ $assembler == ray ]]; then 
        printf "`date` \nray is running... \n\n"
        # kmer input to Ray must be odd
        if [[ -z $kmer ]] || [[ $kmer -gt 160 ]]; then
            kmer=51
        fi
        if [ $((kmer%2)) -eq 0 ]; then ((kmer+=1)); echo "kmer $kmer"; else echo "kmer: $kmer"; fi
        gunzip *fastq.gz
        if [[ $readtype == paired ]]; then
            Ray -o assembly-ray -p ${read1%.gz} ${read2%.gz} &> raylog
            file="./assembly-ray/Scaffolds.fasta"
        else
            Ray -o assembly-ray -s ${read1%.gz} &> raylog &> raylog
            file="./assembly-ray/Scaffolds.fasta" 
        fi
        gzip  *fastq
    fi
    
    histo
    printf '%s is done without error\n\n' $assembler

fi
echo "$LINENO"  

if [[ ! -s $file ]]; then
    printf "*** error\n*** error\n*** error\n"
    printf "A proper assembled file did not complete\n"
    printf "Exiting script \n\n"
    exit 1
fi
echo $file
echo "$LINENO"

if [[ ! -s $kmer  ]]; then
    log3="Kmer size: $kmer"
fi

# Make histograph of contig sizes
echo "file $file"
    
################################################################################################################################
#|||||||||||||||||||||||||||||||||||||||||| FUNCTION:  PARSE BLAST XML FROM VELVET CONTIGS |||||||||||||||||||||||||||||||||||||
################################################################################################################################
function parseXML () {

# Create "here-document" to prevent a dependent file.
cat >./instantparse.py <<EOL
#!/usr/bin/env python

import xml.etree.ElementTree as ET
from sys import argv

script, input = argv

mytree = ET.parse(input)
myroot = mytree.getroot()

for Blast_iteration in myroot.findall('BlastOutput_iterations'):
    for Iteration in Blast_iteration.findall('Iteration'):
        queryID = Iteration.find('Iteration_query-def').text.split("_")
        #length = queryID[-3]
        #coverage = queryID[-1]
        #totalseq = int(length) * float(coverage)
        for Iteration_hits in Iteration.findall('Iteration_hits'):
            for Hit in Iteration_hits.findall('Hit'):
                Hit_def = Hit.find('Hit_def').text
                Hit_def_remove = Hit_def.replace("PREDICTED: ", "")
                Hit_def_split = Hit_def_remove.split(' ')
                Hit_accession = Hit.find('Hit_accession').text
                Hit_len = Hit.find('Hit_len').text
                for Hit_hsps in Hit.findall('Hit_hsps'):
                    for Hsp in Hit_hsps.findall('Hsp'):
                        Hsp_bitscore = Hsp.find('Hsp_bit-score').text
                        Hsp_evalue = Hsp.find('Hsp_evalue').text
                print (queryID, Hit_def_split[0], Hit_def_split[1], Hit_accession, Hit_len, Hsp_bitscore, Hsp_evalue, Hit_def)

EOL

chmod 755 ./instantparse.py

./instantparse.py ${n}_contig-blastResults-5.xml | sed 's/\[.*\]//g' | sort -k1,1 -k2,2 | awk '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' | awk '{k=$1; a[k]++; b[k]=$0}; END{for (k in a) print a[k], b[k]}' | sort -rnk1,1 | awk 'BEGIN{print "n", "acc", "length", "score", "e-value", "ID"} {print $0}' | column -t > BLAST-summary-${n}-contigIDs.txt

./instantparse.py ${n}_contig-blastResults-5.xml > BLAST-all-${n}-contigIDs.txt

rm ./instantparse.py

}
################################################################################################################################

if [[ -z $searchterm ]]; then 
    printf "Reads were not verified with search term\n"
    printf "Provide search term if wishing to verify reads\n"
    help
else
    h=`grep -c ">" ${file}`
    a=`grep -c ".*" ${file}`
    echo "h: $h a: $a"
    singlelinetest=`expr $a / $h`
    echo "singlelingtest: $singlelinetest"
    
    contig_count=`grep -c ">" ${file}`  
    log5=`printf "Total contig count: %'.0f\n" $contig_count`

    if [[ $singlelinetest == 2 ]]; then
            echo "scaffolds file ( $file ) is formatted correctly"
            cat $file > ${n}-contigsoriginal.fa
    else
            # If contigs going in are not all on the same line then newlines are removed to put fasta on a single line, but the name header name is changed.
            echo "contig file ( $file ) was reformated placing sequence on a sigle line"
            awk '{if ($0 ~ />/ ) {print ">"} else print $0}' $file | tr -d "\n" | sed -e 's:>:\n>\n:g' | awk '{if ($0 ~  />/ ) {print ">contigs-" x++} else print $0}' | grep -v "^$" > ${n}-contigsoriginal.fa
    fi
    
    # Remove everything behond the space.  Meant for ">name otherinfo" to just ">name"
    sed 's/ .*$//' ${n}-contigsoriginal.fa > ${n}-contigs.fa
    
    # Cut contigs down to a length no longer than 450 bases to for faster BLAST results.  Remove contigs < 40 bases.
    awk '{ if ($0 ~ /^>/ ) { print $0 } else if (length($0) > 600 ) {print substr($0, 50, 500) } else {print $0}}' ${n}-contigs.fa | grep -B1 '.\{39\}' | grep -v "^--$" > ${n}-contigs2.fa
    echo "##### Blasting contig file #####"
    echo "Start Time: `date`"
    
    blastContigs=`grep -c ">" ${n}-contigs2.fa`
    echo "$blastContigs BLASTed contigs"
    
    blastn -query ${n}-contigs2.fa -db ${NTDB} -task blastn -num_threads ${NR_CPUS} -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1
    parseXML
    wait
    
    totalreads=`grep -v ">" ${n}-contigs.fa | tr -d "\n" | wc -m`
    rm ${n}-contigsoriginal.fa
    mv ${n}-contigs2.fa BLAST-${n}-INFILE.fasta
    
    # Pass all arguements as search terms.
    echo $searchterm >> search
    
    grep -F -i -h -f search BLAST-all-${n}-contigIDs.txt | awk '{print $1}' | sed 's/\[//' | sed 's/\]//' | tr -d "\'" > cherrypickedheader.txt
    
    mkdir blast-output
    mv BLAST-${n}-INFILE.fasta ./blast-output
    mv *contigIDs.txt ./blast-output
    mv ${n}_contig-blastResults-5.xml ./blast-output
    mv ${n}-contigs.fa ./blast-output
    
    echo ""
    echo "Search terms:"
    cat search
    echo ""
    grep -F -i -A1 -h -w -f cherrypickedheader.txt ./blast-output/${n}-contigs.fa | grep -v "^--$" > verifiedreads-${n}.fasta
    echo "Headers found"
    cat cherrypickedlabeledheader.txt
    
    foundreads=`grep -v ">" verifiedreads-${n}.fasta | tr -d "\n" | wc -m`
    
    printf "Total sequence %'.0f\n" $totalreads 
    printf  "Verified %s sequence: %'.0f\n" $searchterm $foundreads 
    
    log6=`printf "\nTotal sequence %'.0f\n" $totalreads` 
    log7=`printf  "Verified %s sequence: %'.0f\n" $searchterm $foundreads` 
    perc=`awk -v x=$foundreads -v y=$totalreads 'BEGIN { print (x / y)*100 }'`
    log8=`printf "Percent of total: %s\n" $perc`
    echo "Percent of total: $perc"
    
    echo ""
    echo "cherrypick.sh has completed"
    echo "see file:  verifiedreads-${n}.fasta"
    echo ""
    
    email_list="tod.p.stuber@usda.gov"
    
    rm cherrypickedheader.txt
    rm search
fi

# Kraken

printf "\nRunning Kraken\n"

if [[ -z $krakendb ]]; then
    printf "\nKraken not ran\n\n"
else
    if [[ $krakendb == jhu_flu ]]; then 
        krakenDatabase="/home/shared/databases/kraken/flu_jhu/fludb_20150820_with_hosts"
    elif [[ $krakendb == std ]]; then 
        krakenDatabase="/home/shared/databases/kraken/std/"
    elif [[ $krakendb == host ]]; then
        krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
    else
        printf "\nIncorrect argument! Provide first argument: jhu_flu, std, host\n\n"
        exit 1
    fi
    
    if [[ $readtype == paired ]]; then
        kraken --db ${krakenDatabase} --threads ${NR_CPUS} --paired ${read1} ${read2} > ${n}-output-kraken.txt && kraken-report --db ${krakenDatabase} ${n}-output-kraken.txt > ${n}-report-kraken.txt
    else
        echo "kraken single"
        kraken --db ${krakenDatabase} --threads ${NR_CPUS} ${read1} > ${n}-output-kraken.txt && kraken-report --db ${krakenDatabase} ${n}-output-kraken.txt > ${n}-report-kraken.txt
    fi
    
    if [[ $krakendb == jhu_flu ]]; then
    
        echo "------> Building Krona Graph... using JHU kraken2krona.sh"
        date
    
        kraken2krona.sh -i ${n}-output-kraken.txt -k ${krakenDatabase} -o ${n}-jhu-output.txt -r ${n}-jhu-Krona_id_graphic.html
    
    elif [[ $krakendb == std ]] || [[ $krakendb == host ]]; then
        cut -f2,3 ${n}-output-kraken.txt > ${n}-kronaInput.txt
            ktImportTaxonomy -a ${n}-kronaInput.txt
            mv taxonomy.krona.html ${n}-taxonomy.krona.html
            mv taxonomy.krona.html.files ${n}-taxonomy.krona.html.files
    else
        printf "\nIncorrect argument! Provide first argument: jhu_flu, std, host\n\n"
        exit 1
    fi
    
    mkdir kraken
    mv *kronaInput.txt *output-kraken.txt *report-kraken.txt *taxonomy.krona.html *taxonomy.krona.html.files kraken
fi

# Removed Trimmed reads
rm ${strain}_Trimmed_R*

endtime=`date +%s`
runtime=$((endtime-starttime))

printf 'Runtime: %dh:%dm:%ds\n\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) > summary-verReads.txt
cat fileinfo >> summary-verReads.txt
rm fileinfo
printf "\nAssembly --> \n" >> summary-verReads.txt
echo "$log1" >> summary-verReads.txt
echo "$log3" >> summary-verReads.txt
echo "$log4" >> summary-verReads.txt 
echo "$log2" >> summary-verReads.txt 
echo "$log5" >> summary-verReads.txt 
echo "$log6" >> summary-verReads.txt
echo "$log7" >> summary-verReads.txt
echo "$log8" >> summary-verReads.txt
if [ "$log9" ]; then
    echo "$log9" >> summary-verReads.txt
fi

contigs_verified=`grep -c ">" verifiedreads-${n}.fasta`
printf "\nVerified %s contigs: %'.0f\n" $searchterm $contigs_verified >> summary-verReads.txt 

contig_select >> summary-verReads.txt

printf "\nKraken -->\n" >> summary-verReads.txt
printf "database: $krakendb\n" >> summary-verReads.txt

printf "summary\n" >> summary-verReads.txt
awk '$4 == "D" || $4 == "U" {print $1, $2, $4, $5, $6}' ./kraken/*-report-kraken.txt | column -t >> summary-verReads.txt 
printf "\n\n" >> summary-verReads.txt

enscript ./blast-output/BLAST-summary-${n}-contigIDs.txt -B -j -r -f "Courier5" -o - | ps2pdf - ./blast_summary.pdf

#fix making histo graph and add back "${n}_fig.pdf" as attachment
if [ ! "$vflag" ]; then 
    echo "Sending email..."
    (cat summary-verReads.txt; pwd) |  mutt -s "Sameple: ${n}" -a ./blast_summary.pdf ./kraken/${n}-*.html -- "tod.p.stuber@usda.gov"
fi

rm ./blast_summary.pdf

# created 2015-01-07, Tod Stuber 
