#!/bin/sh

#################################################################################
#mkdir starting_files
#cp * starting_files

NR_CPUS=30

echo "unzipping files..."
# Unzip files if needed, else put std error to null
find . -name "*gz" -type f -print0 | xargs -0 -n 1 -P $NR_CPUS gunzip 2> /dev/null

onemismatch () {
patt=($1)
for ((i=0; i<${#patt[0]}; i++)); do
    patt+=( "${patt[0]:0:i}.${patt[0]:i+1}" )
done
regex=$(IFS='|'; echo "${patt[*]}")

}

for i in *_R1*; do

#(
echo "**********************ID FILES START**********************"

Bsbv1_5_1F="AGTTCCTCTCGCCGCTTTCA"
Bsbv1_5_1P="CGAAAACCGCGCGCGCCTG"
Bsbv1_5_1R="AAAGCCGCCGCGCAATG"
Bsbv4_8_1F="GCTTCACGCCGCGCTC"
Bsbv4_8_1Falt="CCCACCAGCGTCCGGTAG"
Bsbv4_8_1P="AGCGTCCGGTAGCCCTTGGGC"
Bsbv4_8_1Palt="TCTCGGCCTGTGATGCGCC"
Bsbv4_8_1R="GTCCCTGCGTTCGCAGAATG"
Bsbv4_9_1F="TCGATATCTGGGCGCACCG"
Bsbv4_9_1P="CGAGGATCGCATTGCCGGTGAC"
Bsbv4_9_1R="AGGGCAGGTCTATTCACGCG"
Bsbv4_9_2F="GTTATCGTCCACCACCGGCT"
Bsbv4_9_2P="CCGATCGCCCCAACCTGTACACC"
Bsbv4_9_2R="CGGTGCGCCCAGATATCGA"

targets=(Bsbv1_5_1F Bsbv1_5_1P Bsbv1_5_1R Bsbv4_8_1F Bsbv4_8_1Falt Bsbv4_8_1P Bsbv4_8_1Palt Bsbv4_8_1R Bsbv4_9_1F Bsbv4_9_1P Bsbv4_9_1R Bsbv4_9_2F Bsbv4_9_2P Bsbv4_9_2R)

forReads=`echo $i`                                                                           
echo "Forward Reads:  $forReads"                                                             

name=`echo $i`                                                                               
n=`echo $name | sed 's/\..*//'`                                                              
echo "$n"
echo "$n" > $n.findings

for t in ${targets[@]}; do
	(
#echo "target: $t"
	#echo "opened target: ${!t}"
	
	onemismatch ${!t}
	mismatch=`egrep -c $regex $forReads`
	findings=`egrep -c ${!t} $forReads`
	echo "${t}  $findings $mismatch"
) &
    let count+=1
    [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done | sort -k1,1 >> $n.findings
wait
done
paste *findings > columns.txt

#
#  Created by Stuber, Tod P - APHIS on 02/04/2015
#
