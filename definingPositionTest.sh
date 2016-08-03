#!/bin/sh

#################################################################################
#   Set variables:
reference="/Volumes/Data_HD/Mycobacterium/Go_To_File/NC_002945.fasta"
vcfFilter=yes #(yes or no)
groupFilteredSNPs="/Volumes/Data_HD/Mycobacterium/Mark_files/bovisGroups/" # File containing regions to avoid applied to all vcfs
groupFilter=yes #(yes or no)
#coverageFiles="/Volumes/Data_HD/Mycobacterium/coverageFiles"
coverageFiles="/Volumes/Data_HD/Mycobacterium/coverageFiles"
#   Sed searches
tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."
tbNumberOnly='s/.*\([0-9]\{2\}-[0-9,FM]\{4,6\}\).*/\1/' #Only tb Number
dropEXT='s/\(.*\)\..*/\1/' #Just drop the extention from the file
fulDir=$PWD
# Number of Computer cores
NR_CPUS=18
genotypingcodes="/Users/tstuberadmin/Desktop/Untitled.txt"
DefiningSNPs="/Volumes/Data_HD/Mycobacterium/DefiningSNPsGroupDesignations.txt"
#################################################################################

for i in *.txt; do
    mv $i ${i%.txt}.vcf
done

for i in *.vcf; do
    echo "******************** Naming convention ********************"
    echo "Original File: $i"
    base=`basename "$i"`
    searchName=`echo $base | sed $tbNumberV | sed $tbNumberW | sed 's/V//'`
    echo "searchName: $searchName"
    # Direct script to text file containing a list of the correct labels to use.
    # The file must be a txt file.
    p=`grep "$searchName" "$genotypingcodes"`
    echo "This is what was found in tag file: $p"
    newName=`echo $p | awk '{print $1}' | tr -d "[:space:]"` # Captured the new name
    n=`echo $base | sed $tbNumberV | sed $tbNumberW`
    noExtention=`echo $base | sed $dropEXT`
    VALtest=`echo $i | grep "VAL"`
    echo "VALtest: $VALtest"
    #Check if a name was found in the tag file.  If no name was found, keep original name, make note in log and cp file to unnamed folder.
    if [[ -z "$p" ]]; then # new name was NOT found
        if [[ -z "$VALtest" ]]; then
            name=$searchName
            echo "n is $n"
            echo "--> $name was not found in the FileMaker database" >> log
            mv $i ${name}.vcf
            echo "A"
        else
            name=${searchName}-Val
            mv $i ${name}.vcf
            echo "B"
        fi
    else # New name WAS found
        if [[ -z "$VALtest" ]]; then
            name=$newName
            mv $i ${name}.vcf
            echo "C"
        else
            name=${newName}-Val
            echo "newName is $name"
            mv $i ${name}.vcf
            echo "D"
        fi
    fi
done

list=`awk '{print $2}' ${DefiningSNPs}`
for pos in $list; do
    echo "" >> position_quality.txt
    cluster=`grep "$pos" ${DefiningSNPs} | awk '{print $2, $4"-"$1}'`
    echo "******************* $cluster ********************"
    echo "******************* $cluster ********************" >> position_quality.txt
    for i in *.vcf; do
        m=`basename "$i"`; n=`echo $m | sed $tbNumberV | sed $tbNumberW`
        awk -v x=$pos '$2 ~ "^"x"$" { print FILENAME, $2, $6, $8}' $i | awk 'BEGIN {FS=";"} {print $1, $4, $5}' | sed 's/=-/=/g' | sed 's/\(.*DP=[0-9]\{1,5\}\).*/\1/g' | sed 's/BaseQRankSum=[0-9]\.[0-9]\{1,4\}//g' | awk 'BEGIN {OFS="\t"} {print $1, $3","$4","$5}' >> position_quality.txt
    done
done

function AConeCallPosition () {

positionList=`awk ' { print $2 }' "${DefiningSNPs}" | awk ' NF > 0 '`

echo "AConeCallPosition.sh is running"
echo "*********************************************************************" >> log
echo "Possible Mixed Isolates" >> log
echo "Defining SNPs that are called as AC=1" >> log
echo "" >> log
for i in *.vcf; do
(echo "Finding possible AC1 matches for $i"; for pos in $positionList; do awk -v x=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ "^"x"$" ) print FILENAME, "Pos:", $2, "QUAL:", $6, $8 }' $i; done | grep "AC=1;A" | awk 'BEGIN {FS=";"} {print $1, $2}' >> log) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done
wait
sleep 2

echo "*********************************************************************" >> log
}


