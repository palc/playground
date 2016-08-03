#!/bin/sh

mkdir findseq
cp *.gz findseq
cd findseq

# Unzip files
find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip

forReads=`ls | grep _R1`

revReads=`ls | grep _R2`

n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'` #grab name minus the .vcf

#############

onemismatch () {
patt=($forward)
for ((i=0; i<${#patt[0]}; i++)); do
patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
echo "${patt[*]}" | tr " " "\n" > searchpattern

patt=($reverse)
for ((i=0; i<${#patt[0]}; i++)); do
patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
echo "${patt[*]}" | tr " " "\n" >> searchpattern
}

##############

#IS900F
forward="CCTTTCTTGAAGGGTGTTCG"
reverse="CGAACACCCTTCAAGAAAGG"
onemismatch
IS900F=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#IS900R
forward="CCACCAGATCGGAACGTC"
reverse="GACGTTCCGATCTGGTGG"
onemismatch
IS900R=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#F57F
forward="CCCGATAGCTTTCCTCTCCT"
reverse="AGGAGAGGAAAGCTATCGGG"
onemismatch
F57F=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#F57R
forward="GATCTCAGACAGTGGCAGGTG"
reverse="CACCTGCCACTGTCTGAGATC"
onemismatch
F57R=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#IS1245F
forward="GGTGAGCGGATCACTCAAG"
reverse="CTTGAGTGATCCGCTCACC"
onemismatch
IS1245F=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#IS1245R
forward="GAATCCGCAGTTCCAGGTC"
reverse="GACCTGGAACTGCGGATTC"
onemismatch
IS1245R=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#IS1311F
forward="GCTGGACGCATTACGCAATG"
reverse="CATTGCGTAATGCGTCCAGC"
onemismatch
IS1311F=`cat $forReads $revReads | egrep -h -c -f searchpattern`

#IS1311R
forward="CGCAACTCCAAATCGCCAG"
reverse="CTGGCGATTTGGAGTTGCG"
onemismatch
IS1311R=`cat $forReads $revReads | egrep -h -c -f searchpattern`

rm searchpattern
cd ..
rm -r ./findseq

echo "$n *************"
echo "IS900F: $IS900F"
echo "IS900R: $IS900R"
echo "F57F: $F57F"
echo "F57R: $F57R"
echo "IS1245F: $IS1245F"
echo "IS1245R: $IS1245R"
echo "IS1311F: $IS1311F"
echo "IS1311R: $IS1311R"

