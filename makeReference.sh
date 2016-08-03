#!/bin/sh

#Script to make a reference.  Selected isolates that had the largest genomes making the assumption that the largest genome is the most intact "oldest" representation.  Interestingly the largest genomes turned out to be the oldest isolates that have been sequenced.


tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."

# Coverage below this value will be changed to N
Ncov=1

tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9]\{4,6\}\).*/\1/' #Only tb Number
dropEXT='s/\(.*\)\..*/\1/' #Just drop the extention from the file name.

mkdir ./orginalFastqs
cp * ./orginalFastqs
gunzip *.gz

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

#echo "Starting read counts---" > log.txt
loopNumber=1
accListNumber=0
fastaNumber=1
oldfilesize=10
root=`pwd`
echo  "$root"
logfile=${root}/log.txt
matchfile=${root}/match.txt
targetlist=${root}/targetList.txt
newContigsBlastNT=${root}/newContigsBlastNT
prinseq=${root}/prinseq.txt
echo "Start Time" > $logfile
date >> $logfile
echo "Start Time" > $matchfile
date >> $matchfile
echo "PcntCov MatchingReads RefCoverage RefSize RefName" | awk '{printf "%-10s|%-13s|%-13s|%-13s|%s\n", $1, $2, $3, $4, $5}' >> $matchfile

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi

m=`basename ${forReads}`; n=`echo $m | sed 's/_.*//' | sed 's/\..*//'`
echo "***Sample naming convention:  ${n}"
###########################################

################################
# PARSE BLAST XML FROM VELVET CONTIGS
function parseXML () {
# Parse the XML file.

echo "Cnt ConLn Coverage NTs Accession Length BitScore evalue id id id - - - - - - -" | awk '{ printf "%-12s| %-12s| %-14s| %-12s| %-16s| %-10s| %-11s| %-14s| %-22s| %-22s| %-4s %-4s %-4s %-4s %-4s %-4s %-4s %-4s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}' > _${n}contigIDs_of_unassembled_reads.txt

parseXML-Blast5.py ${n}_contig-blastResults-5.xml | sort -k11,11 -k12,12 -k7,7 | awk '{ printf "%-12s| %-12s| %-14s| %-12s| %-16s| %-10s| %-11s| %-14s| %-22s| %-22s| %-4s %-4s %-4s %-4s %-4s %-4s %-4s %-4s\n", NR, $1, $2, $3, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' >> _${n}contigIDs_of_unassembled_reads.txt
}
################################

####### START FASTA LOOP #######
function fastaLoop () {
m=`basename ${ref}`; nameref=`echo $m | sed 's/\..*//'`
echo "nameref: $nameref"
echo "n: $n"
echo "ref: $ref"
/Users/Shared/_programs/bwa-0.7.5a/bwa index $ref
samtools faidx $ref
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${nameref}.dict

if [ -s ${ref}.fai ] && [ -s ${nameref}.dict ]; then
echo "Index and dict are present, continue script"
else
sleep 5
echo "Either index or dict for reference is missing, try making again"
samtools faidx $ref
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${nameref}.dict
if [ -s ${ref}.fai ] && [ -s ${nameref}.dict ]; then
read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
fi
fi

# Using -B 1 achieved a slightly higher reference coverage
/Users/Shared/_programs/bwa-0.7.5a/bwa mem -M -B 1 -t 5 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > ${n}-${acc}.sam
samtools view -bh -T $ref ${n}-${acc}.sam > ${n}.raw.bam
echo "***Sorting Bam"
samtools sort ${n}.raw.bam ${n}-${acc}.sorted
echo "***Indexing Bam"
samtools index ${n}-${acc}.sorted.bam

java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -out_mode EMIT_ALL_SITES -I ${n}-${acc}.sorted.bam -o ${n}-${acc}.ALLready-mem.vcf -nt 8

totCount=`samtools view -c ${n}-${acc}.sorted.bam`
mapCount=`samtools view -c -F 4 ${n}-${acc}.sorted.bam`
unmapCount=`samtools view -c -f 4 ${n}-${acc}.sorted.bam`
echo "Total Reads: ${totCount}, Unmapped Reads: ${unmapCount}, Mapped Reads: ${mapCount}" >> $logfile

countNTs=`grep -v "#" ${n}-${acc}.ALLready-mem.vcf | grep -c "."`
covCount=`grep -v "#" ${n}-${acc}.ALLready-mem.vcf | awk '{print $2, $4, $5, $6, $8}' | sed 's/\(.*\) .*DP=\([0-9]*\).*/\1 \2/g' | awk '{if ( $5 != "." ) print $0}' | grep -c "."`
echo "Reference NTs with coverage: ${covCount}, NT in reference: $countNTs" >> $logfile

if [ $mapCount -gt 0 ]; then
echo "In THEN line $LINENO"
declare -i x=${covCount}
declare -i y=${countNTs}
perc=`awk -v x=$x -v y=$y 'BEGIN { print (x / y)*100 }'`
printf "%-10s|%-13s|%-13s|%-13s --%s\n" "${perc}%" "$mapCount" "${covCount}" "${countNTs}" "$refName" >> $matchfile
echo "                                                         Total Reads: ${totCount}, Unmapped Reads: ${unmapCount}, Mapped Reads: ${mapCount}" >> $matchfile
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> $matchfile

fi

echo "fastaNumber: $fastaNumber"
echo "fastaNumber: $fastaNumber" >> $logfile
grep -v "#" *.ALLready-mem.vcf | awk '{if ($5 == "." || $6 < 300) print $4; else print $5}' | tr -d "\n" | sed "s/^/>${acc}-${fastaNumber};/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > ${acc}-${fastaNumber}_composite.fasta
rm ${n}-${acc}.sam
rm ${n}.raw.bam
#rm ${n}_${acc}.unmappedReads.sam
rm *.dict
rm ${ref}.fai
mkdir ${n}-${loopNumber}-${acc}-${fastaNumber}
cp *.bam* ./${n}-${loopNumber}-${acc}-${fastaNumber}
#mv ${forReads} ./${n}-${loopNumber}-${acc}-${fastaNumber}
#mv ${n}_${acc}.coverage ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${ref} ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${ref}.amb ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${ref}.ann ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${ref}.bwt ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${ref}.pac ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${ref}.sa ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${n}-${acc}.ALLready-mem.vcf ./${n}-${loopNumber}-${acc}-${fastaNumber}
mv ${n}-${acc}.ALLready-mem.vcf.idx ./${n}-${loopNumber}-${acc}-${fastaNumber}
ref=${acc}-${fastaNumber}_composite.fasta
echo "REFERENCE: ${ref}"
newfilesize=`wc $ref | awk '{print $3}'`
echo "fastaNumber: $fastaNumber"
echo "oldfileszie: $oldfilesize"
echo "newfilesize: $newfilesize"

if [ $mapCount == 0 ]; then
echo "Zero mapped fragment to reference, Ditching Reference."
echo $mapCount
pwd
elif [ "$newfilesize" == "$oldfilesize" ]; then
echo "IN THE THEN STATEMENT"
# If the number of fragments aligning to the target is greater than 10 fragments then send to virusContigs.fa to be BLASTed with NT

echo $mapCount
pwd
if [ $mapCount -gt 10 ]; then
cat ${acc}-${fastaNumber}_composite.fasta >> virusContigs.fa
pwd
grep -v "#" ./${n}-${loopNumber}-${acc}-${fastaNumber}/${n}-${acc}.ALLready-mem.vcf | awk '{print $2, $4, $5, $6, $8}' | sed 's/\(.*\) .*DP=\([0-9]*\).*/\1 \2/g' | awk '{if ( $5 == "." ) print "-"; else if ( $5 == 1 ) print "N"; else if ( $3 ~ /[A-Z]/ ) print $3; else print $2}' | tr -d "\n" | sed "s/^/>${n}-${acc};/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > ./${n}-${loopNumber}-${acc}-${fastaNumber}/${n}-${acc}_bestSeqfromAlignment.fasta

else
rm -r ./${n}-${loopNumber}-${acc}*/
fi
# Place unmapped reads into R1 and R2 fastq files.
samtools view -h -f4 ${n}-${acc}.sorted.bam > ${n}_${acc}.unmappedReads.sam
java -Xmx8g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${n}_${acc}.unmappedReads.sam FASTQ=${n}-${loopNumber}_${acc}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-${loopNumber}-${acc}-unmapped_R2.fastq
rm ${n}_${acc}.unmappedReads.sam
rm *.bam*
rm *.fasta
rm ${forReads}
rm ${revReads}
rm ${acc}-${fastaNumber}_composite.fasta
fastaNumber=$(( $fastaNumber + 1))

else
echo "JUMPED INTO THE ELSE STATEMENT"
oldfilesize=$newfilesize
rm *.bam*
fastaNumber=$(( $fastaNumber + 1))
fastaLoop
echo "Just completed if statement, iteration: ${fastaNumber}"
fi

}
####### END FASTA LOOP #######

############################################################################
############################################################################
############################################################################
############################################################################
############################# Start of Script ##############################
############################################################################
############################################################################
############################################################################
############################################################################

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi

m=`basename ${forReads}`; n=`echo $m | sed 's/_.*//' | sed 's/\..*//'`
echo "***Sample naming convention:  ${n}"

###########################################

#Make file of initial targets.

for i in `cat acc*`; do
grep ">" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${i}.fasta >> $targetlist
done
echo "_____________________" >> $targetlist
###########################################
# Improve virus selection
#read -p "--> Line Number: $LINENO, Press Enter to continue"
acclist=`cat acc*`
echo acc* >> $logfile
#echo "Accesion list: ${acclist}"

for i in ${acclist}; do
acc=$i
echo "Number ${loopNumber}, ${acc}" >> $logfile
echo "Number ${loopNumber}, ${acc}"

ls /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb > /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt
p=`grep "${acc}" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt`

if [[ -z "$p" ]]
then
echo "Downloading from NCBI"
writeFile="${root}/${acc}.fasta"
echo "This is the file to write to:  $writeFile"
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
wait
if [ -s $writeFile ]; then
echo "Downloaded from NCBI, Good to continue."
else
echo "Try downloading again"
sleep 10
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
sleep 5
if [ ! -s $writeFile ]; then
read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
fi
fi

cp $writeFile /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/
else
echo "File is local"
cp /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${p} ./
writeFile=${p}
fi

grep ">" $writeFile >> $logfile
refName=`grep ">" $writeFile`

ref=${writeFile}

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi

###########################################
fastaLoop
###########################################

loopNumber=$(( $loopNumber + 1))
echo "loopNumber: $loopNumber"
echo "____________________________________" >> $logfile
date >> $logfile
done
read -p "$LINENO Finished????"

# Remove Ns from the virusContig file before Blastings
#sed 's/NNNNN//g' virusContigs.fa > virusContigs.temp
#rm virusContigs.fa
#mv virusContigs.temp virusContigs.fa

# Count reads
contigNum=`grep -c ">" virusContigs.fa`

echo "$contigNum contigs assembled"
echo "$contigNum new contigs assembled" >> $logfile

# Used megablast for longer reads versus blastn which is better for shorter reads
blastn -query virusContigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/db/nt -task megablast -num_threads 20 -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1

function newContigsB5 () {
# Parse the XML file.

echo "Cnt BuiltWith id id Accession Length BitScore evalue id id id - - - - - - -" | awk '{ printf "%-10s| %-15s| %-15s| %-15s| %-14s| %-8s| %-9s| %-9s| %-17s| %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' > ${n}-${loopNumber}_newContigs_BlastNT.txt

newContigs-parseXML-Blast5.py ${n}_contig-blastResults-5.xml | awk '{ printf "%-10s| %-15s| %-15s| %-15s| %-14s| %-8s| %-9s| %-9s| %-17s| %s %s %s %s %s\n", NR, $1, $2, $3, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' >> ${n}-${loopNumber}_newContigs_BlastNT.txt
}

newContigsB5
cat ${n}-${loopNumber}_newContigs_BlastNT.txt >> $newContigsBlastNT
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Second Search Started of Virus Database #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> $logfile
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Second Search Started of Virus Database #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> $matchfile

############################################################################
############################################################################
# Using the unmapped reads Blast virus database again an look for more hits.
############################################################################
############################################################################

mkdir First_Search_VirusDB
mv ./*-*-* ./First_Search_VirusDB/
mv ./accessions.txt ./First_Search_VirusDB
mv ./virusBLAST ./First_Search_VirusDB
mv ./First_Search_VirusDB/*.fastq ./

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi

m=`basename ${forReads}`; n=`echo $m | sed 's/_.*//' | sed 's/\..*//'`
echo "***Sample naming convention:  ${n}"

### VIRUS DATABASE TARGET ###
if [ $virusData == yes ]; then

echo "Aligning Fragments"
velveth ./ 31 -fastq -short -separate ${forReads} ${revReads}
velvetg ./ -ins_length 500 -cov_cutoff auto -min_contig_lgth 25 -unused_reads yes

echo "Blast contigs against virus database"
blastn -query contigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/virus_db/virus -task blastn -num_threads 5 -outfmt 5 -out blastResults-5.xml -max_target_seqs 1
# Count reads

contigNum=`grep -c ">" contigs.fa`
UnusedNum=`grep -c ">" UnusedReads.fa`

echo "Contigs used to Blast virus database" >> $logfile
echo "$contigNum contigs assembled"
echo "$contigNum contigs assembled" >> $logfile
echo "$contigNum contigs assembled" >> $matchfile
echo "$UnusedNum UnusedReads"
echo "$UnusedNum UnusedReads" >> $logfile
echo "$UnusedNum UnusedReads" >> $matchfile

printf "\n" >> accessions.txt
parseXML-Blast5.py blastResults-5.xml | sort -nk9,9 | awk '{print $4}' | sed 's/.*\(NC_[0-9]*\).*/\1/g' | tail -50 | sort | uniq -du >> accessions.txt

#Clean folder
mkdir virusBLAST
mv contigs.fa ./virusBLAST
mv blastResults-5.xml ./virusBLAST
mv stats.txt ./virusBLAST
mv UnusedReads.fa ./virusBLAST

rm Graph2
rm LastGraph
rm Log
rm PreGraph
rm Roadmaps
rm Sequences

fi
###########################################

#Make file of initial targets.
for i in `cat acc*`; do
grep ">" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${i}.fasta >> $targetlist
done
echo "_____________________" >> $targetlist
###########################################
# Improve virus selection
#read -p "--> Line Number: $LINENO, Press Enter to continue"
acclist=`cat acc*`
echo acc* >> $logfile
#echo "Accesion list: ${acclist}"

for i in ${acclist}; do
acc=$i
echo "Number ${loopNumber}, ${acc}" >> $logfile
echo "Number ${loopNumber}, ${acc}"

ls /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb > /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt
p=`grep "${acc}" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt`

if [[ -z "$p" ]]
then
echo "Downloading from NCBI"
writeFile="${root}/${acc}.fasta"
echo "This is the file to write to:  $writeFile"
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
wait
if [ -s $writeFile ]; then
echo "Downloaded from NCBI, Good to continue."
else
echo "Try downloading again"
sleep 10
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
sleep 5
if [ ! -s $writeFile ]; then
read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
fi
fi


cp $writeFile /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/
else
echo "File is local"
cp /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${p} ./
writeFile=${p}
fi

grep ">" $writeFile >> $logfile
refName=`grep ">" $writeFile`

ref=${writeFile}

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi

###########################################
fastaLoop
###########################################

loopNumber=$(( $loopNumber + 1))
echo "loopNumber: $loopNumber"
echo "____________________________________" >> $logfile
date >> $logfile
done

# Remove Ns from the virusContig file before Blastings
#sed 's/NNNNN//g' virusContigs.fa > virusContigs.temp
#rm virusContigs.fa
#mv virusContigs.temp virusContigs.fa

# Count reads
contigNum=`grep -c ">" virusContigs.fa`

echo "$contigNum contigs assembled"
echo "$contigNum new contigs assembled" >> $logfile

# Used megablast for longer reads versus blastn which is better for shorter reads
blastn -query virusContigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/db/nt -task megablast -num_threads 20 -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1

newContigsB5
echo "_____________________" >> $newContigsBlastNT
cat ${n}-${loopNumber}_newContigs_BlastNT.txt >> $newContigsBlastNT

mkdir Second_Search_VirusDB
mv ./*-*-* ./Second_Search_VirusDB
mv ./accessions.txt ./Second_Search_VirusDB
mv ./virusBLAST ./Second_Search_VirusDB
mv ./Second_Search_VirusDB/*.fastq ./

echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* Searching NT Database Using Remaining Contigs #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> $logfile
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* Searching NT Database Using Remaining Contigs #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> $matchfile

######################################
# Create the First NT Accession List #
######################################

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"
pwd

# Align fragments that have not assembled to selected accessions
velveth ./ 31 -fastq -short -separate ${forReads} ${revReads}
velvetg ./ -ins_length 500 -cov_cutoff auto -min_contig_lgth 100 -unused_reads yes

# Count reads
contigNum=`grep -c ">" contigs.fa`
UnusedNum=`grep -c ">" UnusedReads.fa`

echo "$contigNum contigs assembled"
echo "$contigNum contigs assembled" >> $logfile
echo "$contigNum contigs assembled" >> $matchfile
echo "$UnusedNum UnusedReads"
echo "$UnusedNum UnusedReads" >> $logfile
echo "$UnusedNum UnusedReads" >> $matchfile

#echo "##### Blasting first 600 lines of contig file #####"
#head -600 contigs.fa > headcontigs.fa
blastn -query contigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/db/nt -task blastn -num_threads 20 -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1

accListNumber=$(( $accListNumber + 1))
echo "accListNumber: $accListNumber"
#read -p "--> Line: $LINENO What is being parsed in the XML file.  Is it doing what I thing or is an error coming back?"
parseXML
wait
#read -p "--> Line: $LINENO After the parseXML.  Is there a file or did it throw an error?"
#This will output a accession list of 100 top scored accessions ordered by genome size largest to smallest
sed 1d *contigIDs_of_unassembled_reads.txt | awk ' $15 == 0 { print $0 }' | sort -nk13,13 | awk ' { if (a[$9]++ == 0) print $0; }' "$@" | tail -100 | sort -rnk11,11 | awk '{print $9}' > accession${accListNumber}.txt
#read -p "--> Line: $LINENO What is in the Accession list if one was made?"

rm Graph2
rm LastGraph
rm Log
rm PreGraph
rm Roadmaps
rm Sequences
mv stats.txt assembly-stats.txt

mkdir ./assembly
mv *xml ./assembly
mv assembly-stats.txt ./assembly
mv contigs.fa ./assembly
mv headcontigs.fa ./assembly
mv UnusedReads.fa ./assembly

#######################################################################################################################
############################################### BEGIN LOOP ACCESSIONS #################################################
#######################################################################################################################

#Started at the function's end

function loopAccessions () {
# Align reads to genomes listed in Accession file.
# Use the unmapped reads to loop back as working down list.
#read -p "--> Line Number: $LINENO, Press Enter to continue"
acclist=`cat acc*`
echo acc* >> $logfile
#echo "Accesion list: ${acclist}"

#######################
for i in ${acclist}; do
acc=$i
echo "Number ${loopNumber}, ${acc}" >> $logfile
echo "Number ${loopNumber}, ${acc}"

ls /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb > /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt
p=`grep "${acc}" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt`

if [[ -z "$p" ]]
then
echo "Downloading from NCBI"
writeFile="${root}/${acc}.fasta"
echo "This is the file to write to:  $writeFile"
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
wait
if [ -s $writeFile ]; then
echo "Downloaded from NCBI, Good to continue."
else
echo "Try downloading again"
sleep 10
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
sleep 5
if [ ! -s $writeFile ]; then
read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
fi
fi

cp $writeFile /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/
else
echo "File is local"
cp /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${p} ./
writeFile=${p}
fi

grep ">" $writeFile >> $logfile
refName=`grep ">" $writeFile`

# Grab reads and reference and place them in variables
ref=${writeFile}

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi

m=`basename ${ref}`; nameref=`echo $m | sed 's/\..*//'`
echo "nameref: $nameref"
echo "n: $n"
echo "ref: $ref"
/Users/Shared/_programs/bwa-0.7.5a/bwa index $ref
samtools faidx $ref
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${nameref}.dict

if [ -s ${ref}.fai ] && [ -s ${nameref}.dict ]; then
echo "Index and dict are present, continue script"
else
sleep 5
echo "Either index or dict for reference is missing, try making again"
samtools faidx $ref
java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.100/CreateSequenceDictionary.jar REFERENCE=${ref} OUTPUT=${nameref}.dict
if [ -s ${ref}.fai ] && [ -s ${nameref}.dict ]; then
read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
fi
fi

echo "Current working directory: `pwd`"

# Using -B 1 achieved a slightly higher reference coverage
/Users/Shared/_programs/bwa-0.7.5a/bwa mem -M -B 1 -t 5 -R @RG"\t"ID:"$n""\t"PL:ILLUMINA"\t"PU:"$n"_RG1_UNIT1"\t"LB:"$n"_LIB1"\t"SM:"$n" $ref $forReads $revReads > ${n}-${acc}.sam

samtools view -bh -T $ref ${n}-${acc}.sam > ${n}.raw.bam
echo "***Sorting Bam"
samtools sort ${n}.raw.bam ${n}-${acc}.sorted
echo "***Indexing Bam"
samtools index ${n}-${acc}.sorted.bam

java -Xmx2g -jar /Users/Shared/_programs/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -out_mode EMIT_ALL_SITES -I ${n}-${acc}.sorted.bam -o ${n}-${acc}.ALLready-mem.vcf -nt 8

totCount=`samtools view -c ${n}-${acc}.sorted.bam`
mapCount=`samtools view -c -F 4 ${n}-${acc}.sorted.bam`
unmapCount=`samtools view -c -f 4 ${n}-${acc}.sorted.bam`
#percentMap=`expr $mapCount / $totCount`
#echo "Total Reads: ${totCount}, Unmapped Reads: ${unmapCount}, % Mapped: ${percentMap}%, Mapped Reads: ${mapCount}" >> $logfile
echo "Total Reads: ${totCount}, Unmapped Reads: ${unmapCount}, Mapped Reads: ${mapCount}" >> $logfile

#countNTs=`grep -v ">" ${ref} | grep -o [ATGC] | wc -l`
#bamtools coverage -in ${n}-${acc}.sorted.bam > ${n}_${acc}.coverage 2> myerror.txt
#covCount=`awk '{if ($3 > 0) ++b } END {print b}' ${n}_${acc}.coverage`
#echo "Reference NTs with coverage: ${covCount}, NT in reference: $countNTs" >> $logfile

countNTs=`grep -v "#" ${n}-${acc}.ALLready-mem.vcf | grep -c "."`
covCount=`grep -v "#" ${n}-${acc}.ALLready-mem.vcf | awk '{print $2, $4, $5, $6, $8}' | sed 's/\(.*\) .*DP=\([0-9]*\).*/\1 \2/g' | awk '{if ( $5 != "." ) print $0}' | grep -c "."`
echo "Reference NTs with coverage: ${covCount}, NT in reference: $countNTs" >> $logfile

declare -i x=${covCount}
declare -i y=${countNTs}
perc=`awk -v x=$x -v y=$y 'BEGIN { print (x / y)*100 }'`

echo "Total Reads: ${totCount}, Unmapped Reads: ${unmapCount}, % Mapped: ${perc}%, Mapped Reads: ${mapCount}" >> $logfile

if [ $mapCount -gt 0 ]; then
echo "Jumped into THEN statement"
pwd

printf "%-10s|%-13s|%-13s|%-13s --%s\n" "${perc}%" "$mapCount" "${covCount}" "${countNTs}" "$refName" >> $matchfile
echo "                                                         Total Reads: ${totCount}, Unmapped Reads: ${unmapCount}, Mapped Reads: ${mapCount}" >> $matchfile
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> $matchfile

grep -v "#" ${n}-${acc}.ALLready-mem.vcf | awk '{print $2, $4, $5, $6, $8}' | sed 's/\(.*\) .*DP=\([0-9]*\).*/\1 \2/g' | awk '{if ( $5 == "." ) print "-"; else if ( $5 == 1 ) print "N"; else if ( $3 ~ /[A-Z]/ ) print $3; else print $2}' | tr -d "\n" | sed "s/^/>${n}-${acc};/" | tr ";" "\n" | sed 's/[A-Z],[A-Z]/N/g' > ${n}-${acc}_bestSeqfromAlignment.fasta

rm ${n}-${acc}.sam
rm ${n}.raw.bam

mkdir ${n}-${loopNumber}-${acc}
mv ./${n}-${acc}.ALLready-mem.vcf ./${n}-${loopNumber}-${acc}/${n}-${acc}.ALLready-mem.vcf
mv ./${n}-${acc}.ALLready-mem.vcf.idx ./${n}-${loopNumber}-${acc}/${n}-${acc}.ALLready-mem.vcf

mv ./${n}-${acc}_bestSeqfromAlignment.fasta ./${n}-${loopNumber}-${acc}/${n}-${acc}_bestSeqfromAlignment.fasta
mv ${forReads} ./${n}-${loopNumber}-${acc}
mv ${revReads} ./${n}-${loopNumber}-${acc}
#mv ${n}_${acc}.coverage ./${n}-${loopNumber}-${acc}
mv ${ref} ./${n}-${loopNumber}-${acc}
mv ${ref}.amb ./${n}-${loopNumber}-${acc}
mv ${ref}.ann ./${n}-${loopNumber}-${acc}
mv ${ref}.bwt ./${n}-${loopNumber}-${acc}
mv ${ref}.pac ./${n}-${loopNumber}-${acc}
mv ${ref}.sa ./${n}-${loopNumber}-${acc}

# Place unmapped reads into R1 and R2 fastq files.
samtools view -h -f4 ${n}-${acc}.sorted.bam > ${n}_${acc}.unmappedReads.sam
java -Xmx8g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${n}_${acc}.unmappedReads.sam FASTQ=${n}-${loopNumber}_${acc}-unmapped_R1.fastq SECOND_END_FASTQ=${n}-${loopNumber}-${acc}-unmapped_R2.fastq
rm *.dict
rm *.fasta.fai
mv ${n}_${acc}.unmappedReads.sam ./${n}-${loopNumber}-${acc}/${n}_${acc}.unmappedReads.sam
mv *.bam* ./${n}-${loopNumber}-${acc}

else
echo "Jumped into ELSE statement"
pwd

rm ./${n}-${acc}.ALLready-mem.vcf
rm ./${n}-${acc}.ALLready-mem.vcf.idx
rm ${n}-${acc}.sam
rm ${n}.raw.bam
rm ${n}_${acc}.unmappedReads.sam
rm *.bam*
rm ${ref}
rm ${ref}.amb
rm ${ref}.ann
rm ${ref}.bwt
rm ${ref}.pac
rm ${ref}.sa
rm *.dict
rm *.fasta.fai

fi

loopNumber=$(( $loopNumber + 1))
echo "loopNumber: $loopNumber"
echo "____________________________________" >> $logfile
date >> $logfile
done
#######################
pwd
forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi


# Align fragments that have not assembled to selected accessions
echo "Aligning Fragments"
velveth ./ 31 -fastq -short -separate ${forReads} ${revReads}
velvetg ./ -ins_length 500 -cov_cutoff auto -min_contig_lgth 100 -unused_reads yes

# Count reads
contigNum=`grep -c ">" contigs.fa`
UnusedNum=`grep -c ">" UnusedReads.fa`


echo "$contigNum contigs assembled"
echo "$contigNum contigs assembled" >> $logfile
echo "$contigNum contigs assembled" >> $matchfile
echo "$UnusedNum UnusedReads"
echo "$UnusedNum UnusedReads" >> $logfile
echo "$UnusedNum UnusedReads" >> $matchfile

if [ $contigNum -lt 500 ] || [ $accListNumber -gt 3 ]; then
rm accession*.txt
# Local blast on contig file return only the highest match
echo "Less than 500 contigs or accession loop ended"
echo "Blastn on leftover contigs: $contigNum"
blastn -query contigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/db/nt -task blastn -num_threads 20 -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1
mv contigs.fa usedendContigFile.txt
#Function to parse XML
parseXML

###########################################
###########################################
###########################################
pwd
rm Graph2
rm LastGraph
rm Log
rm PreGraph
rm Roadmaps
rm Sequences
rm stats.txt
rm *xml
rm contigs.fa

mkdir ./Final_alignments
mv ./${n}-*-*/ ./Final_alignments
mv acc* ./Final_alignments

echo "Blasting Unused Reads against virus database: $UnusedNum Reads"
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Unused Reads against virus database #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> $matchfile
blastn -query UnusedReads.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/virus_db/virus -task blastn -num_threads 20 -outfmt 5 -out blastResults-5.xml -max_target_seqs 1

parseXML-Blast5.py blastResults-5.xml > UnusedReadsBlast.txt

awk '$11 {a[$11]++} END { for (i in a) {printf "%s\t\t%s\n",i , a[i] }}' UnusedReadsBlast.txt | awk '$2 > 10 {print $0}' | sed 's/.*\(NC_[0-9]*\).*/\1/g' > accession-UnusedReads.txt


###########################################

#Make file of initial targets.
echo "From Unused Reads" >> $targetlist
for i in `cat acc*`; do
grep ">" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${i}.fasta >> $targetlist
done
echo "_____________________" >> $targetlist

###########################################
# Improve virus selection
#read -p "--> Line Number: $LINENO, Press Enter to continue"
acclist=`cat acc*`
echo acc* >> $logfile
#echo "Accesion list: ${acclist}"

for i in ${acclist}; do
acc=$i
echo "Number ${loopNumber}, ${acc}" >> $logfile
echo "Number ${loopNumber}, ${acc}"

ls /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb > /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt
p=`grep "${acc}" /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/list.txt`

if [[ -z "$p" ]]
then
echo "Downloading from NCBI"
writeFile="${root}/${acc}.fasta"
echo "This is the file to write to:  $writeFile"
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
wait
if [ -s $writeFile ]; then
echo "Downloaded from NCBI, Good to continue."
else
echo "Try downloading again"
sleep 10
echo "Running python script to grab fasta file"
python -u /Users/Shared/_programs/Python/GetFASTAbyGI.py $acc $writeFile
sleep 5
if [ ! -s $writeFile ]; then
read -p "Fasta file ${acc} failed to download from NCBI, Manually download if desired to proceed.  Line: $LINENO"
fi
fi

cp $writeFile /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/
else
echo "File is local"
cp /Users/Shared/_programs/ncbi-blast-2.2.27+/mydb/${p} ./
writeFile=${p}
fi

grep ">" $writeFile >> $logfile
refName=`grep ">" $writeFile`

ref=${writeFile}

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

if [ -s $forReads ]; then
echo "At $LINENO R1 read file exists"
else
"Aborting script at line number $LINENO, R1 read file is empty"
exit 1 #error status
fi


####################
fastaLoop
wait
#echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" >> $matchfile
####################

loopNumber=$(( $loopNumber + 1))
echo "loopNumber: $loopNumber"
echo "____________________________________" >> $logfile
date >> $logfile

done

# Remove Ns from the virusContig file before Blastings
#sed 's/NNNNN//g' virusContigs.fa > virusContigs.temp
#rm virusContigs.fa
#mv virusContigs.temp virusContigs.fa

# Count reads
contigNum=`grep -c ">" virusContigs.fa`

echo "$contigNum contigs assembled"
echo "$contigNum new contigs assembled" >> $logfile

# Used megablast for longer reads versus blastn which is better for shorter reads
blastn -query virusContigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/db/nt -task megablast -num_threads 20 -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1

#function newContigsB5 () {
# Parse the XML file.

#echo "Cnt BuiltWith id id Accession Length BitScore evalue id id id - - - - - - -" | awk '{ printf "%-10s| %-15s| %-15s| %-15s| %-14s| %-8s| %-9s| %-9s| %-17s| %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' > ${n}-${loopNumber}_newContigs_BlastNT.txt

#newContigs-parseXML-Blast5.py ${n}_contig-blastResults-5.xml | awk '{ printf "%-10s| %-15s| %-15s| %-15s| %-14s| %-8s| %-9s| %-9s| %-17s| %s %s %s %s %s\n", NR, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' >> ${n}-${loopNumber}_newContigs_BlastNT.txt
#}

newContigsB5
echo "_____________________" >> $newContigsBlastNT
echo "Alingments from the Unused Reads of Assembles -->" >> $newContigsBlastNT
cat ${n}-${loopNumber}_newContigs_BlastNT.txt >> $newContigsBlastNT

###########################################
###########################################
###########################################
#Clean folder
mv ${root}/log.txt ${root}/_${n}alignment_log.txt

mkdir ./Final_alignments-UnusedReads
mv ${n}-*-*/ ./Final_alignments-UnusedReads
mv usedendContigFile.txt ./Final_alignments-UnusedReads
mv UnusedReadsBlast.txt ./Final_alignments-UnusedReads
mv UnusedReads.fa ./Final_alignments-UnusedReads
mv blastResults.xml ./Final_alignments-UnusedReads
mv accession-UnusedReads.txt ./Final_alignments-UnusedReads
mv ${n}* ./Final_alignments-UnusedReads

echo "${accListNumber} Accession Lists Used" >> $matchfile
cat $matchfile > ${root}/matching-fragments.txt
cp _${n}contigIDs_of_unassembled_reads.txt ${root}/FINAL_${n}contigIDs_of_unassembled_reads.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" > ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Filtered Reads  #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
cat $prinseq >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Virus Target List #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
cat $targetlist >> ${root}/_${n}-SummaryFile.txt
printf "\n" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Blast NT of Fasta Composites | Gaps filed with Ref Sequence #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
cat $newContigsBlastNT >> ${root}/_${n}-SummaryFile.txt
printf "\n" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*# Aligned Reads #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
cat ${root}/matching-fragments.txt >> ${root}/_${n}-SummaryFile.txt
printf "\n" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* Blast NT Resutls of Final Contig Assembly #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
echo "#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#" >> ${root}/_${n}-SummaryFile.txt
cat ${root}/FINAL_${n}contigIDs_of_unassembled_reads.txt >> ${root}/_${n}-SummaryFile.txt
echo "End Time" >> ${root}/_${n}-SummaryFile.txt
date >> ${root}/_${n}-SummaryFile.txt
awk 'sub("$", "\r")' ${root}/_${n}-SummaryFile.txt > PC-Notebook_${n}-SummaryFile.txt

# E-mail results
basen=`basename $0`
( cat PC-Notebook_${n}-SummaryFile.txt; uuencode PC-Notebook_${n}-SummaryFile.txt PC-Notebook_${n}-SummaryFile.txt ) | mail -s "${n} $basen completed | Open with Notebook Program" tod.p.stuber@aphis.usda.gov

echo "#*#*#*#*#*#*#*# DONE #*#*#*#*#*#*#*#"

else
#echo "##### Blasting first 600 lines of contig file #####"
#head -600 contigs.fa > headcontigs.fa
blastn -query contigs.fa -db /Users/Shared/_programs/ncbi-blast-2.2.27+/db/nt -task blastn -num_threads 20 -outfmt 5 -out ${n}_contig-blastResults-5.xml -max_target_seqs 1
rm contigs.fa
rm headcontigs.fa
#Function to parse XML
parseXML

mv ${n}_contig-blastResults-5.xml XMLfileUsed-${loopNumber}.txt

mkdir NT_Search-${accListNumber}
mv ./*-*-* ./NT_Search-${accListNumber}/
mv ./accession*.txt ./NT_Search-${accListNumber}
mv ./XMLfileUsed* ./NT_Search-${accListNumber}
rm Graph2
rm LastGraph
rm Log
rm PreGraph
rm Roadmaps
rm Sequences
rm stats.txt
rm -r ./assembly/
rm UnusedReads.fa
mv ./NT_Search-${accListNumber}/*.fastq ./
mv ./NT_Search-${accListNumber}/*contigIDs_of_unassembled_reads.txt ./

#This will output a accession list of 100 top scored accessions ordered by genome size largest to smallest
sed 1d _${n}contigIDs_of_unassembled_reads.txt | awk ' $15 == 0 { print $0 }' | sort -nk13,13 | awk ' { if (a[$9]++ == 0) print $0; }' "$@" | tail -100 | sort -rnk11,11 | awk '{print $9}' > ./accession${accListNumber}.txt
mv *contigIDs_of_unassembled_reads.txt ./NT_Search-${accListNumber}

wait
echo "#*#*#*#*#*#*#*# NEW ACCESSION LIST CREATED #*#*#*#*#*#*#*#"
echo "#*#*#*#*#*#*#*# NEW ACCESSION LIST CREATED, ACCESSION LIST $accListNumber #*#*#*#*#*#*#*#" >> $logfile

accListNumber=$(( $accListNumber + 1))
echo "accListNumber: $accListNumber"

loopAccessions

fi
}

##### Start Script #####
loopAccessions

#
#  Created by Stuber, Tod P - APHIS on 1/18/2014.
#