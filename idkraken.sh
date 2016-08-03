#!/bin/sh

# working directory must contain fastq files (can be zipped) to be identified.  
# files can be paired or single read

NR_CPUS=20

if [[ $1 == jhu_flu ]]; then 

	krakenDatabase="/home/shared/databases/kraken/flu_jhu/fludb_20150820_with_hosts"

elif [[ $1 == std ]]; then 

	krakenDatabase="/home/shared/databases/kraken/std/"

elif [[ $1 == host ]]; then
	
	krakenDatabase="/home/shared/databases/kraken/host-bac-vir"
else

    echo ""
    echo ""
    echo "Incorrect argument! Provide first argument: jhu_flu, std, host"
    echo ""
    echo "Example: [prompt]$ idkraken.sh jhu_flu <usda email address>"
    echo ""
    echo "jhu_flu -- flu optimized, individual segment identification"
    echo "std -- Kraken's standard bacteria and virus database"
    echo "host -- Includes human, horse, cow, pig, chicken, babesia, bacteria and virus genomes"
    exit 1

fi 

if [[ $2 == tod ]]; then 
	email="tod.p.stuber@usda.gov"
elif [[ -n $2 ]]; then
	email="$2"
else
	echo ""
	echo ""
	echo "Incorrect argument! Provide second argument: email name or address, options: tod, or <usda email address>"
	echo ""
    	echo "Example: [prompt]$ idkraken.sh jhu_flu <usda email address>"
	echo ""
	exit 1
fi

mkdir originalreads
cp *fastq* originalreads

# Unzip files if needed, else put std error to null
find . -maxdepth 1 -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip 2> /dev/null

#Set sample name- sample name is the name of the fastq files minus any identifying information from the sequencer
sampleName=`ls *.fastq | head -1 | sed 's/_.*//' | sed 's/\..*//'`

root=`pwd`

#Establish Read Files
if [ -f *_R2* ]; then
    echo "R2 paired end read file present"
    export sampleType="paired"
    forFile=`ls | grep _R1`
    forReads="$root/$forFile"
    echo "Forward Reads:  $forReads"
    revFile=`ls | grep _R2`
    revReads="$root/$revFile"
    echo "Reverse Reads:  $revReads"
else
    echo "Just a single read present"
    export sampleType="single"
    forFile=`ls | grep fastq`
    forReads="$root/$forFile"
    echo "Forward Reads:  $forReads"
fi

echo "Kraken database selected is: $krakenDatabase"
date
echo "*** Kraken is finding reads"

#Run Kraken
if [[ $sampleType == "paired" ]]; then
	kraken --db ${krakenDatabase} --threads ${NR_CPUS} --paired *fastq* > $sampleName-kraken-output.txt && kraken-report --db ${krakenDatabase} $sampleName-kraken-output.txt > $sampleName-kraken-report.txt
else
	kraken --db ${krakenDatabase} --threads ${NR_CPUS}  $forReads > $sampleName-kraken-output.txt && kraken-report --db ${krakenDatabase} $sampleName-kraken-output.txt > $sampleName-kraken-report.txt
fi

date
echo "*** Krona transforming Kraken output to graph"

if [[ $1 == jhu_flu ]]; then

	echo "------> Building Krona Graph... using JHU kraken2krona.sh"
	date

	kraken2krona.sh -i $sampleName-kraken-output.txt -k ${krakenDatabase} -o $sampleName-jhu-output.txt -r $sampleName-jhu-Krona_id_graphic.html

elif [[ $1 == std ]] || [[ $1 == host ]]; then
	cut -f2,3 $sampleName-kraken-output.txt > $sampleName-kronaInput.txt
        ktImportTaxonomy -a $sampleName-kronaInput.txt
        mv taxonomy.krona.html $sampleName-taxonomy.krona.html
        mv taxonomy.krona.html.files $sampleName-taxonomy.krona.html.files
else

    echo ""
    echo "Incorrect argument! Provide argument: jhu_flu, std"
    echo ""
    echo "Example: [prompt]$ idkraken.sh jhu_flu"
    echo ""
    exit 1

fi

rm ./*fastq

echo "Sending email..."
echo "idkraken.sh $1 has completed" | mutt -s "Sameple: $sampleName" -a $sampleName-*.html -- "tod.p.stuber@usda.gov" 
echo ""
echo "*** Finished running Kraken on $sampleName"
echo ""

# created 2015-06-29, tstuber

