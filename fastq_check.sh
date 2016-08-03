#!/bin/sh

echo "" > /data/id/pedv/reports/fastq_check_run.txt

find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P 60 gunzip

for infile in *fastq*; do 

# Run fastqc on file
fastqc $infile

singleheader=$(head -1 $infile | grep "^@")

# Check that there is a header in file
if [[ -z $singleheader ]]; then 
	echo "There is an error in the $infile fastq header"
	exit 1
else
	echo "$infile header: $singleheader"
fi

# Check platform
if [[ $(echo $singleheader | sed 's/[^:]//g' | awk '{print length}') == 2 ]]; then
	platform="IonTorrent"
	echo "Ion Torrent fastq"
elif [[  $(echo $singleheader | sed 's/[^:]//g' | awk '{print length}') == 9 ]]; then
	platform="Illumina2"
		if [[ $(echo $singleheader | awk '{print substr($0, 1, 2) }' | sed 's/@//') == H ]]; then
        		echo "Illumina Hiseq fastq"
		elif [[ $(echo $singleheader | awk '{print substr($0, 1, 2) }' | sed 's/@//') == M ]]; then
        		echo "Illumina Miseq fastq"
		else
        		echo "Unable to determine if Miseq or Hiseq data"
        		echo "$singleheader"
		fi
else
	echo "Unable to determine platform"
	echo "$singleheader"
fi

#Parse headers
if [[ $platform == IonTorrent ]]; then 
	echo "Ion Torrent is the platform"
elif [[ $platform == Illumina2 ]]; then
	
	instrumentname=$(grep -B 2 "^+$" $infile | grep "^@" | awk 'BEGIN{FS = "[: ]"; OFS = "\t"} {k=$1; b[k]=$1}; END{for (k in b) print b[k]}'| sed 's/@//')
	
	echo "***File name: $infile"
	echo "***File name: $infile"  >> /data/id/pedv/reports/log_all.txt
	echo "--> File name: $infile"  >> /data/id/pedv/reports/fastq_check_run.txt	
	for i in `echo $instrumentname`; do
		instrumentnumber=$(grep $i $infile | awk 'BEGIN{FS = "[: ]"; OFS = "\t"} {k=$1; a[k]++}; END{ print a[k]}' & )
		runnumbers=$(grep $i $infile | awk 'BEGIN{FS = "[: ]"; OFS = "\t"} {k=$2; a[k]++; b[k]=$2}; END{for (k in a) print a[k] " reads from run " b[k]}' & )	
		readnumbers=$(grep $i $infile | awk 'BEGIN{FS = "[: ]"; OFS = "\t"} {k=$8; a[k]++; b[k]=$8}; END{for (k in a) print  a[k] " reads from R" b[k]}' & )
		indexseqs=$(grep $i $infile | awk 'BEGIN{FS = "[: ]"; OFS = "\t"} {k=$10; a[k]++; b[k]=$10}; END{for (k in a) print a[k] " reads have an index of " b[k]}' & )
	wait
	echo "*Instrument: $i"
	echo "$runnumbers"
	echo "Read numbers $readnumbers"
	echo "$indexseqs"
	echo ""		

	echo "*Instrument: $i" >> /data/id/pedv/reports/log_all.txt
        echo "$runnumbers" >> /data/id/pedv/reports/log_all.txt
        echo "Read numbers $readnumbers" >> /data/id/pedv/reports/log_all.txt
        echo "$indexseqs" >> /data/id/pedv/reports/log_all.txt
	date  >> /data/id/pedv/reports/log_all.txt
        echo "" >> /data/id/pedv/reports/log_all.txt

	echo "" >> /data/id/pedv/reports/fastq_check_run.txt
	echo "*Ran on instrument: $i" >> /data/id/pedv/reports/fastq_check_run.txt
        echo "Read length: $(head -2 $infile | grep -v "^@" | wc | awk '{print $3}')" >> /data/id/pedv/reports/fastq_check_run.txt
	echo "$runnumbers" >> /data/id/pedv/reports/fastq_check_run.txt
        echo "$readnumbers" >> /data/id/pedv/reports/fastq_check_run.txt
        echo "$indexseqs" >> /data/id/pedv/reports/fastq_check_run.txt

	#echo "$infile, Instrument: $i, $runnumbers, $readnumbers, $indexseqs" >> /data/id/pedv/reports/fastq_check_run.txt
	echo "$infile, Instrument: $i, $runnumbers, $readnumbers, $indexseqs" >> /data/id/pedv/reports/${i}.txt 

	done
fi

echo "" >> /data/id/pedv/reports/fastq_check_run.txt
echo "*************************************************" >> /data/id/pedv/reports/fastq_check_run.txt
echo "" >> /data/id/pedv/reports/fastq_check_run.txt


done

date  >> /data/id/pedv/reports/fastq_check_run.txt

email_list="tod.p.stuber@aphis.usda.gov"
#email_list="tod.p.stuber@aphis.usda.gov suelee.robbe-austerman@aphis.usda.gov"

find . -name "*_fastqc.html" >> list

cat /data/id/pedv/reports/fastq_check_run.txt | mutt -s "fastq check" -a `cat list` -- $email_list

# ~~end
# Created by: Tod Stuber, 2014-11-15
