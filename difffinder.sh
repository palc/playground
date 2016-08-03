#!/bin/sh

# have 2 vcfs is a working directory
# vcfs must be from same reference

vcfone=`ls -1 *vcf | head -1`
vcftwo=`ls -1 *vcf | tail -1`

echo ""
echo "Comparing"
echo "VCF 1: $vcfone"
echo "VCF 2: $vcftwo"

awk ' $0 !~ /^#/ && $6 > 300 || $6 == "." {print $1 "-" $2} ' $vcfone > onepositions
awk ' $0 !~ /^#/ && $6 > 300 || $6 == "." {print $1 "-" $2} ' $vcftwo > twopositions

diff onepositions twopositions | grep "<" | sed 's/< //' | tr '|-' '|\t' > inone
diff onepositions twopositions | grep ">" | sed 's/> //' | tr '|-' '|\t' > intwo

echo "These are uniq to $vcfone" > uniqtoone
while read l; do
	ref=`echo $l | awk '{print $1}'`
	pos=`echo $l | awk '{print $2}'`
	awk -v r=$ref -v p=$pos 'BEGIN{OFS="\t"} $1 == r && $2 == p {print $0}' $vcfone >> uniqtoone 
done < inone
echo "" >> uniqtoone

echo "These are uniq to $vcftwo" > uniqtotwo
while read l; do
        ref=`echo $l | awk '{print $1}'`
        pos=`echo $l | awk '{print $2}'`
        awk -v r=$ref -v p=$pos 'BEGIN{OFS="\t"} $1 == r && $2 == p {print $0}' $vcftwo >> uniqtotwo
done < intwo

cat uniqtoone uniqtotwo > ${vcfone%.vcf}-${vcftwo%.vcf}-comparison.tab 

echo ""
# created 2015-12-02 stuber
