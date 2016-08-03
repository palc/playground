#!/bin/sh

# Arguments
urls=`ls -d -1 $PWD/** | grep "url"`
downloadDir=`pwd`

#################################################################################
# Set variables:

# Sed searches put into variables
tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9,FM]\{4,6\}\).*/\1/' #Only tb Number, *laboratory specific*
dropEXT='s/\(.*\)\..*/\1/' #Just drop the extention from the file

NR_CPUS=18 # Computer cores to use when analyzing

#################################################################################

function abyss_run () {
    
    # Place R1 Reads into variable
    forReads=`ls ./Zips | grep _R1`
    echo "Forward Reads:  $forReads"
    # Place R2 Reads into variable
    revReads=`ls ./Zips | grep _R2`
    echo "Reverse Reads:  $revReads"

    # Get name of isolate to use for naming output files
    n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`
    echo "***Isolate naming convention:  $n"

    # De Novo assembly using ABySS
    abyss-pe name=${n}_abyss k=31 in="./Zips/$forReads ./Zips/$revReads"

    # Copy contig file to set directory
    #cp ${n}_abyss-3.fa $abyssFile

    # Clean-up output
    # Move reads to their own directory
    rm coverage.hist
    mkdir ./${n}_abyss
    mv ${n}_abyss* ./${n}_abyss
    mv ${n}_abyss* ./${n}_abyss
    mv coverage.hist ./${n}_abyss
    unitigsFileSize=`ls -lL ./${n}_abyss | grep "_abyss-3.fa" | awk '{print $5"-bytes"}'`
    unitigs=`grep "_abyss-unitigs.fa" ./${n}_abyss/*abyss-stats`
    scaffolds=`grep "_abyss-scaffolds.fa" ./${n}_abyss/*abyss-stats`
    if [ -e $unitigsFileSize ]; then
        echo "$n ABySS did not assemble reads" >> ${downloadDir}/downloadSummary.txt
        else
        echo "$n unitigs file size: $unitigsFileSize" >> ${downloadDir}/downloadSummary.txt
        echo "$n unitigs summary: $unitigs" >> ${downloadDir}/downloadSummary.txt
        echo "$n scaffolds summary: $scaffolds" >> ${downloadDir}/downloadSummary.txt
    fi
    fileSizes_R1=`ls -lh ./Zips | grep "R1" | awk '{print $5}'`
    fileSizes_R2=`ls -lh ./Zips | grep "R2" | awk '{print $5}'`

    echo "$n R1 zipped file size: $fileSizes_R1" >> ${downloadDir}/downloadSummary.txt
    echo "$n R2 zipped file size: $fileSizes_R2" >> ${downloadDir}/downloadSummary.txt
    echo "*** Done ***"

}

#################################################################################

function spoligoSpacerFinder () {
    
    echo "**********************START**********************"

    spacer01="TGATCCAGAGCCGGCGACCCTCTAT|ATAGAGGGTCGCCGGCTCTGGATCA"
    spacer02="CAAAAGCTGTCGCCCAA|TTGGGCGACAGCTTTTG"
    spacer03="CCGTGCTTCCAGTGATCGCCTTCTA|TAGAAGGCGATCACTGGAAGCACGG"
    spacer04="ACGTCATACGCCGACCAATCATCAG|CTGATGATTGGTCGGCGTATGACGT"
    spacer05="TTTTCTGACCACTTGTGCGGGATTA|TAATCCCGCACAAGTGGTCAGAAAA"
    spacer06="CGTCGTCATTTCCGGCTTCAATTTC|GAAATTGAAGCCGGAAATGACGACG"
    spacer07="GAGGAGAGCGAGTACTCGGGGCTGC|GCAGCCCCGAGTACTCGCTCTCCTC"
    spacer08="CGTGAAACCGCCCCCAGCCTCGCCG|CGGCGAGGCTGGGGGCGGTTTCACG"
    spacer09="ACTCGGAATCCCATGTGCTGACAGC|GCTGTCAGCACATGGGATTCCGAGT"
    spacer10="TCGACACCCGCTCTAGTTGACTTCC|GGAAGTCAACTAGAGCGGGTGTCGA"
    spacer11="GTGAGCAACGGCGGCGGCAACCTGG|CCAGGTTGCCGCCGCCGTTGCTCAC"
    spacer12="ATATCTGCTGCCCGCCCGGGGAGAT|ATCTCCCCGGGCGGGCAGCAGATAT"
    spacer13="GACCATCATTGCCATTCCCTCTCCC|GGGAGAGGGAATGGCAATGATGGTC"
    spacer14="GGTGTGATGCGGATGGTCGGCTCGG|CCGAGCCGACCATCCGCATCACACC"
    spacer15="CTTGAATAACGCGCAGTGAATTTCG|CGAAATTCACTGCGCGTTATTCAAG"
    spacer16="CGAGTTCCCGTCAGCGTCGTAAATC|GATTTACGACGCTGACGGGAACTCG"
    spacer17="GCGCCGGCCCGCGCGGATGACTCCG|CGGAGTCATCCGCGCGGGCCGGCGC"
    spacer18="CATGGACCCGGGCGAGCTGCAGATG|CATCTGCAGCTCGCCCGGGTCCATG"
    spacer19="TAACTGGCTTGGCGCTGATCCTGGT|ACCAGGATCAGCGCCAAGCCAGTTA"
    spacer20="TTGACCTCGCCAGGAGAGAAGATCA|TGATCTTCTCTCCTGGCGAGGTCAA"
    spacer21="TCGATGTCGATGTCCCAATCGTCGA|TCGACGATTGGGACATCGACATCGA"
    spacer22="ACCGCAGACGGCACGATTGAGACAA|TTGTCTCAATCGTGCCGTCTGCGGT"
    spacer23="AGCATCGCTGATGCGGTCCAGCTCG|CGAGCTGGACCGCATCAGCGATGCT"
    spacer24="CCGCCTGCTGGGTGAGACGTGCTCG|CGAGCACGTCTCACCCAGCAGGCGG"
    spacer25="GATCAGCGACCACCGCACCCTGTCA|TGACAGGGTGCGGTGGTCGCTGATC"
    spacer26="CTTCAGCACCACCATCATCCGGCGC|GCGCCGGATGATGGTGGTGCTGAAG"
    spacer27="GGATTCGTGATCTCTTCCCGCGGAT|ATCCGCGGGAAGAGATCACGAATCC"
    spacer28="TGCCCCGGCGTTTAGCGATCACAAC|GTTGTGATCGCTAAACGCCGGGGCA"
    spacer29="AAATACAGGCTCCACGACACGACCA|TGGTCGTGTCGTGGAGCCTGTATTT"
    spacer30="GGTTGCCCCGCGCCCTTTTCCAGCC|GGCTGGAAAAGGGCGCGGGGCAACC"
    spacer31="TCAGACAGGTTCGCGTCGATCAAGT|ACTTGATCGACGCGAACCTGTCTGA"
    spacer32="GACCAAATAGGTATCGGCGTGTTCA|TGAACACGCCGATACCTATTTGGTC"
    spacer33="GACATGACGGCGGTGCCGCACTTGA|TCAAGTGCGGCACCGCCGTCATGTC"
    spacer34="AAGTCACCTCGCCCACACCGTCGAA|TTCGACGGTGTGGGCGAGGTGACTT"
    spacer35="TCCGTACGCTCGAAACGCTTCCAAC|GTTGGAAGCGTTTCGAGCGTACGGA"
    spacer36="CGAAATCCAGCACCACATCCGCAGC|GCTGCGGATGTGGTGCTGGATTTCG"
    spacer37="CGCGAACTCGTCCACAGTCCCCCTT|AAGGGGGACTGTGGACGAGTTCGCG"
    spacer38="CGTGGATGGCGGATGCGTTGTGCGC|GCGCACAACGCATCCGCCATCCACG"
    spacer39="GACGATGGCCAGTAAATCGGCGTGG|CCACGCCGATTTACTGGCCATCGTC"
    spacer40="CGCCATCTGTGCCTCATACAGGTCC|GGACCTGTATGAGGCACAGATGGCG"
    spacer41="GGAGCTTTCCGGCTTCTATCAGGTA|TACCTGATAGAAGCCGGAAAGCTCC"
    spacer42="ATGGTGGGACATGGACGAGCGCGAC|GTCGCGCTCGTCCATGTCCCACCAT"
    spacer43="CGCAGAATCGCACCGGGTGCGGGAG|CTCCCGCACCCGGTGCGATTCTGCG"

    # Starting working directory must be BWA-GATK folder with included /Zips file containing 2 zipped fastq files.
    echo "directory"
    # Make fastq directory
    pwd
    mkdir ./fastq
    cp ./Zips/*.fastq.gz ./fastq

    echo "starting to unzip files"
    # Unzip files
    gunzip ./fastq/*.fastq.gz
    echo "finished unzipping files"

    spoligoforReads=`ls ./fastq | grep _R1`
    echo "Forward Reads:  $spoligoforReads"

    spoligorevReads=`ls ./fastq | grep _R2`
    echo "Reverse Reads:  $spoligorevReads"

    n=`echo $spoligorevReads | sed $tbNumberV | sed $tbNumberW`

    sp01=`egrep $spacer01 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    echo "Spacer 1 for $n finished at `date`"
    sp02=`egrep $spacer02 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp03=`egrep $spacer03 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp04=`egrep $spacer04 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp05=`egrep $spacer05 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    echo "Spacer 5 for $n finished at `date`"
    sp06=`egrep $spacer06 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp07=`egrep $spacer07 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp08=`egrep $spacer08 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp09=`egrep $spacer09 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp10=`egrep $spacer10 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp11=`egrep $spacer11 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp12=`egrep $spacer12 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp13=`egrep $spacer13 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp14=`egrep $spacer14 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp15=`egrep $spacer15 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp16=`egrep $spacer16 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp17=`egrep $spacer17 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp18=`egrep $spacer18 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp19=`egrep $spacer19 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp20=`egrep $spacer20 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp21=`egrep $spacer21 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp22=`egrep $spacer22 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp23=`egrep $spacer23 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp24=`egrep $spacer24 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp25=`egrep $spacer25 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp26=`egrep $spacer26 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp27=`egrep $spacer27 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp28=`egrep $spacer28 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp29=`egrep $spacer29 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp30=`egrep $spacer30 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp31=`egrep $spacer31 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp32=`egrep $spacer32 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp33=`egrep $spacer33 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp34=`egrep $spacer34 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp35=`egrep $spacer35 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp36=`egrep $spacer36 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp37=`egrep $spacer37 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp38=`egrep $spacer38 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp39=`egrep $spacer39 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp40=`egrep $spacer40 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp41=`egrep $spacer41 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp42=`egrep $spacer42 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`
    sp43=`egrep $spacer43 ./fastq/$spoligoforReads ./fastq/$spoligorevReads | wc -l`

    echo "spacer01	spacer02	spacer03	spacer04	spacer05	spacer06	spacer07	spacer08	spacer09	spacer10	spacer11	spacer12	spacer13	spacer14	spacer15	spacer16	spacer17	spacer18	spacer19	spacer20	spacer21	spacer22	spacer23	spacer24	spacer25	spacer26	spacer27	spacer28	spacer29	spacer30	spacer31	spacer32	spacer33	spacer34	spacer35	spacer36	spacer37	spacer38	spacer39	spacer40	spacer41	spacer42	spacer43" > $n.spacer.txt

    echo "$sp01	$sp02	$sp03	$sp04	$sp05	$sp06	$sp07	$sp08	$sp09	$sp10	$sp11	$sp12	$sp13	$sp14	$sp15	$sp16	$sp17	$sp18	$sp19	$sp20	$sp21	$sp22	$sp23	$sp24	$sp25	$sp26	$sp27	$sp28	$sp29	$sp30	$sp31	$sp32	$sp33	$sp34	$sp35	$sp36	$sp37	$sp38	$sp39	$sp40	$sp41	$sp42	$sp43" >> $n.spacer.txt

    cat $n.spacer.txt | awk 'NR==2 {for(i=1;i<=NF;i++) if ($i >= 5) print 1; else print 0}' | tr -cd "[:print:]" | fold -w3 > $n.myspacers

    chmod 755 ./$n.myspacers

    mybinaries=`cat ./$n.myspacers`

    for i in $mybinaries; do
    if [ $i == 000 ]
    then
        echo "0" >> $n.octalcode.txt
    elif [ $i == 001 ]
    then
        echo "1" >> $n.octalcode.txt
    elif [ $i == 010 ]
        then
        echo "2" >> $n.octalcode.txt
    elif [ $i == 011 ]
        then
        echo "3" >> $n.octalcode.txt
    elif [ $i == 100 ]
        then
        echo "4" >> $n.octalcode.txt
    elif [ $i == 101 ]
        then
        echo "5" >> $n.octalcode.txt
    elif [ $i == 110 ]
        then
        echo "6" >> $n.octalcode.txt
    elif [ $i == 111 ]
        then
        echo "7" >> $n.octalcode.txt
    elif [ $i == 0 ]
        then
        echo "0" >> $n.octalcode.txt
    elif [ $i == 1 ]
        then
        # Changed 4 to 1 on 2013-10-25
        echo "1" >> $n.octalcode.txt
    else
        echo "***Error***" >> $n.octalcode.txt
    fi
    done

    WGSpoligo=`cat $n.octalcode.txt | tr -cd "[:print:]"`

    # Add infor to spoligoCheck.txt
    echo "$n -----> $WGSpoligo *****************************" >> ${downloadDir}/downloadSummary.txt
   
    # Add infor to spoligoCheck_all.txt
    echo "<----- $n ----->" >> /Users/Shared/_WGS/spoligoCheck_all.txt
    echo "WGSpoligo:	$WGSpoligo" >> /Users/Shared/_WGS/spoligoCheck_all.txt

}

#################################################################################
################################## START ########################################
#################################################################################

links=`sed 's/ftp:\/\/ftp-trace.ncbi.nlm.nih.gov\///g' ${urls}`

for l in $links; do
    echo "Downloading: $l"
/Users/tstuberadmin/Applications/Aspera\ Connect.app/Contents/Resources/ascp -i /Users/tstuberadmin/Applications/Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh -k1 -Tr -l100m anonftp@ftp-private.ncbi.nlm.nih.gov:/${l} ${downloadDir}
done

for i in ${downloadDir}/*.sra; do
    fastq-dump --split-files --gzip $i
done
wait

ls > listOfFiles.txt
mail -s "***Check file names at $downloadDir $0" tod.p.stuber@aphis.usda.gov < listOfFiles.txt
rm listOfFiles.txt

read -p "Check file names then press enter"
#
#for i in *_1.fastq.gz; do
#    mv $i ${i%_1.fastq.gz}_R1.fastq.gz
#done
#
#for i in *_2.fastq.gz; do
#    mv $i ${i%_2.fastq.gz}_R2.fastq.gz
#done
#
#for i in *_3.fastq.gz; do
#    mv $i ${i%_3.fastq.gz}_R3.fastq.gz
#done
#
#for i in *_4.fastq.gz; do
#    mv $i ${i%_4.fastq.gz}_R4.fastq.gz
#done

for l in $links; do base=`basename $l`; sample=${base%.sra}; echo "The sample name is: $sample"; wget -O ${downloadDir}/${sample}.csv "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= ${sample}"; done

cat ${downloadDir}/*.csv | sort | uniq -ud > ${downloadDir}/metadata
awk 'BEGIN{FS=","}{print $28"-"$1}' ${downloadDir}/metadata | sed 's/Mycobacterium /M/g' | sed 's/tuberculosis/tb/g' >${downloadDir}/genotypeNames.txt
read -p "Paused at line ${LINENO}.  Does the genotypeNames.txt look ready?"

rm *.csv

for i in *.fastq.gz; do
    echo "******************** Changing the file names ********************"
    echo "Original File: $i"
    searchName=`echo $i | sed $tbNumberV | sed $tbNumberW`
    echo "searchName: $searchName"
    # Direct script to text file containing a list of the correct labels to use.
    # The file must be a txt file.
    p=`grep "$searchName" "${downloadDir}/genotypeNames.txt"`
    echo "This is what was found in tag file: $p"
    newName=`echo $p | awk '{print $1}' | tr -d "[:space:]"` # Captured the new name
    echo "newName: $newName"
    newFileName=`echo $i | sed "s/${searchName}/${newName}/"`
    echo "This will be the new file name: $newFileName"
    if [ -e $newFileName ]; then
        echo "File not renamed"
        else
        mv $i $newFileName
    fi
done

# Group samples into their own directories
for i in *.fastq.gz; do
    n=`echo $i | sed $tbNumberV | sed $tbNumberW`
    echo "n is : $n"
    mkdir -p $n
    mv $i $n/
done
currentDir=`pwd`
read -p "Paused at line ${LINENO}. "
for dir in */; do
    echo "dir: $dir"
    cd ${dir}
    n=`echo *_R1* | sed $tbNumberV | sed $tbNumberW`
    mkdir ./Zips
    cat *_R1* *_R3* > temp
    rm *_R1*
    rm *_R3*
    mv temp ./Zips/${n}_R1.fastq.gz
    cat *_R2* *_R4* >temp
    rm *_R2*
    rm *_R4*
    mv temp ./Zips/${n}_R2.fastq.gz
    cd $currentDir
done

for dir in */; do
    (cd ${dir}
     abyss_run
     spoligoSpacerFinder
    cd $currentDir) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait

cat ${downloadDir}/downloadSummary.txt | sort -nk1,1 > ${downloadDir}/downloadSummarysorted.txt

mail -s "$0 completed at $downloadDir" tod.p.stuber@aphis.usda.gov < ${downloadDir}/downloadSummarysorted.txt

#
#  Created by Stuber, Tod P - APHIS on 5/17/2014.
#
