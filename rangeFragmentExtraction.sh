#!/bin/sh

# rangeFragmentExtraction.sh
# Script is used to grab the fragments within a specified range of a bam file, perform a de novo on those fragments and output the results to a specified folder.

# Set working directory to _TB-Data
# Set variables
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9]\{4,6\}\).*/\1/' #Only tb Number

#Targets##
target1F=409392
target1R=412831
##
target2F=824329
target2R=825675
##
target3F=1018193
target3R=1021431
##
target4F=1664165
target4R=1666890
##
target5F=1899833
target5R=1902252
##
target6F=2446534
target6R=2450610
##
target7F=2784921
target7R=2787667
##
target8F=2972494
target8R=2975291
##
target9F=3392326
target9R=3396520
##
target10F=3765427
target10R=3768335
##
target11F=3850395
target11R=3852844
###

echo "**************** START ****************"
# Set the working directory to the _TB-Data
cd /Volumes/Data_HD/Mycobacterium/_TB-Data

mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target1F}-${target1R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target2F}-${target2R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target3F}-${target3R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target4F}-${target4R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target5F}-${target5R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target6F}-${target6R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target7F}-${target7R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target8F}-${target8R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target9F}-${target9R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target10F}-${target10R}"
mkdir "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target11F}-${target11R}"

# Make a list of the directorys, minus the misfit folder.
allDirectories=`ls`
allDirectories=($allDirectories)
allDirectories=${allDirectories[@]//_misf*}

# Iterate through directories
for d in $allDirectories; do
    echo "Current directory: $d"
    cd ./$d/BWA-GATK
    mkdir ./fragments

    # Get name
    name=`echo *.ready.bam | sed 's/\..*//'`

    # Convert the binary Bam file to a readable Sam file
    echo "Coverting ${name}.ready.bam to ${name}.ready.sam"
    samtools view -h ${name}.ready.bam > ./fragments/${name}.ready.sam

    ################################## This part will need to be repeated-1 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target1F} -v y=${target1R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
# Changed to -3 from -7.  file -7 would error sometimes.  Probably to few fragments
    abyss-pe name=${name}_${target1F}-${target1R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target1F}-${target1R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam

    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-2 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target2F} -v y=${target2R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target2F}-${target2R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target2F}-${target2R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-3 ##################################


    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target3F} -v y=${target3R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target3F}-${target3R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target3F}-${target3R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-4 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target4F} -v y=${target4R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target4F}-${target4R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target4F}-${target4R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-5 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target5F} -v y=${target5R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target5F}-${target5R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target5F}-${target5R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-6 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target6F} -v y=${target6R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target6F}-${target6R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target6F}-${target6R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-7 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target7F} -v y=${target7R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target7F}-${target7R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target7F}-${target7R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-8 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target8F} -v y=${target8R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target8F}-${target8R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target8F}-${target8R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-9 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target9F} -v y=${target9R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target9F}-${target9R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target9F}-${target9R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-10 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target10F} -v y=${target10R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target10F}-${target10R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target10F}-${target10R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments/*.fastq
    rm -rf ./fragments/*_range.sam
    rm -rf ./fragments/*_range.bam
    rm -rf ./fragments/*_rangeShuf.bam
    ##################################           End of Repeat           ##################################
    ################################## This part will need to be repeated-11 ##################################
    

    cd ./fragments

    echo "Capturing fragments within range"
    head -20 ${name}.ready.sam | grep "^@" > ${name}_range.sam; awk -v x=${target11F} -v y=${target11R} '$4 >= x && $4 <= y {print $0}' ${name}.ready.sam >> ${name}_range.sam

    echo "Converting .sam to .bam"
    samtools view -bh -T /Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta ${name}_range.sam > ${name}_range.bam

    echo "Shuffling Reads"
    htscmd bamshuf -uOn 128 ${name}_range.bam tmp > ${name}_rangeShuf.bam

    echo "Converting to .fastq"
    java -Xmx2g -jar /Users/Shared/_programs/picard-tools-1.94/SamToFastq.jar INPUT=${name}_rangeShuf.bam FASTQ=${name}_R1.fastq SECOND_END_FASTQ=${name}_R2.fastq

    cd ..
    mkdir abyss
    cd abyss

    echo "Running abyss"
    abyss-pe name=${name}_${target11F}-${target11R} k=60 in="../fragments/*.fastq"
    cp *-3.fa "/Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target11F}-${target11R}"

    # Clean-up files and directories
    cd ..
    rm -rf ./abyss

    rm -rf ./fragments

    ##################################           End of Repeat           ##################################

    cd /Volumes/Data_HD/Mycobacterium/_TB-Data

done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target1F}-${target1R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target2F}-${target2R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target3F}-${target3R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target4F}-${target4R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target5F}-${target5R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target6F}-${target6R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target7F}-${target7R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target8F}-${target8R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target9F}-${target9R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target10F}-${target10R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

cd /Volumes/Data_HD/Mycobacterium/rangeFragmentExtraction/${target11F}-${target11R}

for i in *.fa; do
echo "$i"
getbase=`basename "$i"`
number=`echo $getbase | sed 's/_.*//'`
sed "s/>.*/>$number/g" $i > ${getbase}sta
rm $i
done

echo "****************************** END ******************************"
pwd



#
#  Created by Stuber, Tod P - APHIS on 7/26/2013.
#
