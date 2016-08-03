#!/bin/sh

# buildExclusionList.sh

# This script will compare positions and build exclusion list of low quality SNPs
# Working directory must contain unfiltered VCFs

# Create directory with vcf with coverage above 90X
# To determine coverage use file listing TB numbers with known high coverage

#################################################################################
# Set all position file for reference
reference="/Volumes/Data_HD/Mycobacterium/1315_realignment/1315positions.temp"

#################################################################################


# Do not build exclusion list with isolates with < 75X coverage.
`>dpositions.temp`
`>exclusionList.temp`
`>positions.temp`
`>exclusionTracking.temp`
`>allVCFs.temp`

#echo ${`date`} >> exclusionTracking.temp

# Gather all positions from all VCFs.  VCFs should not be filtered.
# uniq -d positions for normal QUAL and AC2 filter
for i in *; do cat $i | awk '{print $2}' | grep "[0-9]\{4,7\}" | grep -v "[A-Za-z]"; done >> positions.temp
cat positions.temp | sort | uniq -d > dpositions.temp
cat positions.temp | sort | uniq -u > upositions.temp

# uniq -u positions to save for higher QUAL selection
# Iterate through each position grepping positions from VCFs.
# After grepping position if any one position is of low quality the position is placed into a exclusion list.
# There must be at least two SNPs present to be called

echo "******************************** START SNPs that are present 2 or more times ********************************"

dlist=`cat dpositions.temp`

for a in *.vcf; do
grep -v "#" $a >> allVCFs.temp

done

for l in $dlist; do
    # Reset grepPostions.temp file to 0 bytes
    `>grepPositions.temp`
    `>ACtest.temp`

    awk -v x=$l '$2 == x { print $2, $6, $8 }' allVCFs.temp >> grepPositions.temp
    
    #cat grepPostions.temp
    # Filter SNPs on QUAL value
    # Filter SNPs on AC2 value
    # Filter SNPs with DP < 5
        
    # Parse QUAL value
    
    min=`awk 'BEGIN { min = 10000 }; $2 < min {  min = $2 }; END  { printf "%.0f\n", min }' grepPositions.temp`
    echo "Position: $l"
    echo "QUAL min: $min"
    # Parse AC1
    ACtest=`awk '$3 ~ /^AC=1/ { len=length ($3)}; END { print len }' grepPositions.temp`
    echo "ACtest: $ACtest"
    # Parse Depth of coverage
    leastDP=`cat grepPositions.temp | awk '{print $3}' | sed 's/.*DP=\([0-9]\{1,5\}\).*/\1/g' | sort -n | head -1 | awk '{printf "%.0f\n", $1}'`
    echo "Least DP: $leastDP"

    if (( min < 650 ))
        then
        echo "Filtered, Low QUAL"
        # Pass position to exclusion list
        echo $l >> exclusionList.temp
        # Track exclusion
        echo "$l excluded because of QUAL < 650, call: $min QUAL---" >> exclusionTracking.temp
    elif (( ACtest > 10 ))
        then
        echo "Filtered, AC1"
        echo "ACtest: $ACtest"
        echo $l >> exclusionList.temp
        echo "$l excluded because of AC=1, call: $ACtest AC1***" >> exclusionTracking.temp
    elif (( leastDP < 8 ))
        then
        echo "Filtered, DP < 10"
        echo "Least DP: $leastDP"
        echo $l >> exclusionList.temp
        echo "$l excluded because of DP < 8, call: $leastDP DP+++" >> exclusionTracking.temp

    else
        echo "Good SNP"    
    fi
done

echo "******************************** START SNPs that are only once ********************************"
# If only one SNP is present the QUAL value must be of higher value.
ulist=`cat upositions.temp`

for u in $ulist; do
    # Reset grepPostions.temp file to 0 bytes
    `>grepPositions.temp`
    `>ACtest.temp`
    
#grep "$u" allVCFs.temp | awk '{print $2, $6, $8}' >> grepPositions.temp
    awk -v x=$u '$2 == x { print $2, $6, $8 }' allVCFs.temp >> grepPositions.temp

    #cat grepPostions.temp
    # Filter SNPs on QUAL value
    # Filter SNPs on AC2 value
    # Filter SNPs with DP < 5
    
    # Parse QUAL value
    min=`awk 'BEGIN { min = 10000 }; $2 < min {  min = $2 }; END  { printf "%.0f\n", min }' grepPositions.temp`
    echo "Position: $u"
    echo "QUAL min: $min***"
    # Parse AC1
    ACtest=`awk '$3 ~ /^AC=1/ { len=length ($3) }; END { print len }' grepPositions.temp`
    echo "ACtest: $ACtest***"
    # Parse Depth of coverage
    leastDP=`cat grepPositions.temp | awk '{print $3}' | sed 's/.*DP=\([0-9]\{1,5\}\).*/\1/g' | sort -n | head -1 | awk '{printf "%.0f\n", $1}'`
    echo "Least DP: $leastDP***"

    if (( min < 1000 ))
        then
        echo "Filtered, Low QUAL"
        # Pass position to exclusion list
        echo $u >> exclusionList.temp
        # Track exclusion
        echo "$u excluded because of QUAL < 1000, call: $min QUAL---***" >> exclusionTracking.temp
    elif (( ACtest > 10 ))
        then
        echo "Filtered, AC1"
        echo "ACtest: $ACtest"
        echo $u >> exclusionList.temp
        echo "$u excluded because of AC=1, call: $ACtest AC1******" >> exclusionTracking.temp
    elif (( leastDP < 10 ))
        then
        echo "Filtered, DP < 10"
        echo "Least DP: $leastDP"
        echo $u >> exclusionList.temp
        echo "$u excluded because of DP < 10, call: $leastDP DP+++***" >> exclusionTracking.temp
    
    else
        echo "Good SNP"
    fi
done

echo "Number of positions in Exclusion list"
cat "exclusionList.temp" | wc -l
pwd
cp "$0" "$PWD"

echo "Making full exclusion list"
cat positions.temp | sort | uniq > positionsuniq.temp

cat positionsuniq.temp exclusionList.temp | sort | uniq -u > usableSNPs.temp

echo "The number of usable SNPs:"
cat usableSNPs.temp | wc -l

cat $reference usableSNPs.temp | sort | uniq -u > allExclusions.txt

echo "The number of SNPs that will not be considered in the analysis:"
cat allExclusions.txt | wc -l

#
#  Created by Stuber, Tod P - APHIS on 05/07/2013.
#
