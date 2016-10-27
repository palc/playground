#!/bin/sh

starttime=$(date +%s)
alias pause='read -p "$LINENO Enter"'

NR_CPUS=40

help () {
    printf "\n\n"
    printf "### ncorrect argument! Provide first argument: jhu_flu, std, host\n\n"
    printf "### Example: [prompt]$ idkraken.sh jhu_flu <usda email address>\n\n"
    printf '### email address or shorthand: "tod" "jess" or "na"\n'
    printf "### jhu_flu -- flu optimized, individual segment identification\n"
    printf "### std -- Kraken's standard bacteria and virus database\n"
    printf "### host -- Includes human, horse, cow, pig, chicken, babesia, bacteria and virus genomes\n"
    exit 1
}

# check arg1 for database type
if [[ $1 == jhu_flu ]]; then
    krakenDatabase="/home/shared/databases/kraken/flu_jhu/fludb_20150820_with_hosts"
elif [[ $1 == std ]]; then
    krakenDatabase="/home/shared/databases/kraken/std/"
elif [[ $1 == host ]]; then
    krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
else
    help
fi

# check arg2 for email
if [[ $2 == tod ]]; then
    email="tod.p.stuber@usda.gov"
elif [[ $2 == jess ]]; then
    email="Jessica.A.Hicks@aphis.usda.gov"
elif [[ $2 == na ]]; then
    printf "No email being sent\n"
elif [[ -n $2 ]]; then
    email="$2"
else
    help
fi

# zip FASTQs
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
cp *gz original_reads
# Make alias in working directory to zip files

# Kraken requires unzip FASTQs
pigz -d *fastq

# get sample name
strain=$(echo *fastq* | head -1 | sed 's/\..*//')
echo "Quality trimming sample $strain"

# Kraken
date
printf "Kraken running...\n\n"
if [ $filecount -eq 2 ]; then 
    read1=$(ls *fastq | head -1)
    read2=$(ls *fastq | tail -1)
    echo "Forward Reads: $read1"
    echo "Reverse Reads: $read2"
    kraken --db ${krakenDatabase} --threads ${NR_CPUS} --paired ${read1} ${read2} > $strain-outputkraken.txt && kraken-report --db ${krakenDatabase} $strain-outputkraken.txt > $strain-reportkraken.txt
    rm ${read1} ${read2}
elif [ $filecount -eq 1 ]; then
    read1=`ls *fastq`
    echo "Forward Reads: $read1"
    kraken --db ${krakenDatabase} --threads ${NR_CPUS} ${read1} > $strain-outputkraken.txt && kraken-report --db ${krakenDatabase} $strain-outputkraken.txt > $strain-reportkraken.txt
    rm ${read1}
fi

# error if scaffolds.fasta not made
if [[ ! -s $strain-reportkraken.txt ]]; then
    printf "*** error\n*** error\n*** error\n"
    printf "A proper assembled file did not complete\n"
    printf "Exiting script \n\n"
    exit 1
fi

date
printf "*** Krona transforming Kraken output to graph\n"

if [[ $1 == jhu_flu ]]; then
    printf "------> Building Krona Graph... using JHU kraken2krona.sh\n"
    date
    kraken2krona.sh -i $sampleName-kraken-output.txt -k ${krakenDatabase} -o $sampleName-jhu-output.txt -r $sampleName-jhu-Krona_id_graphic.html
else
    printf "------> Building Krona Graph... using standard taxonomy\n"
    cut -f2,3 $sampleName-kraken-output.txt > $sampleName-kronaInput.txt
    # removed -a
    ktImportTaxonomy $sampleName-kronaInput.txt
    mv taxonomy.krona.html $sampleName-taxonomy.krona.html
    mv taxonomy.krona.html.files $sampleName-taxonomy.krona.html.files
fi

# save needed files


# send email
if [[ $2 == na ]]; then
    printf "No email is being sent\n"
else
    printf "Sent email to $2\n"
    echo "idkraken.sh $1 has completed" | mutt -s "Sameple: $sampleName" -a $sampleName-*.html -- ${email}
fi

endtime=`date +%s`
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60))
printf "***DONE\n"
# created 2015-10-26, Tod Stuber
