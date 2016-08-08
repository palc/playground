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

Bov2_2F="CGGATGGCGGCATCGATC"
Bov2_2P="CGGTGAGAACCCGTGCCGAA"
Bov2_2R="GGACCTATATCGGCATTGCAGC"
Bov4_1Falt="CTGGAAACGGCGGAAGTCAC"
Bov4_1P="CGCAACGTGTGATGGCGCTGA"
Bov4_1Ralt="CGGTTGCGGCGACAATCTC"
Bmel2_1F="TGCTGGTGGAAATCGTGGAAC"
Bmel2_1P="CTTGCCCTCGTCAATGGTGGCCG"
Bmel2_1R="ACGTAGCGTTTCCATCCAGGT"
Bs1_2F="TGCGTATCCGTAACCGCATTC"
Bs1_2P="CGAGGGAAATTGTTCTGGACGCGG"
Bs1_2R="TTTAGTAGCCCCCAACGCG"
Bs2_3F="GAAACTGCTCCAGGGGCCT"
Bs2_3P="AGTTTGATGATGGTCCGCGTCAGG"
Bs2_3R="CCACGCGGAACCTCTGTTTG"
Bs2_3Ralt="CACGCGGAACCTCTGTTTGT"
Bsbv3_1_1F="AAAGTTTAATCCGGCAGCGGTG"
Bsbv3_1_1P="CGTTTTGCGGGCGGTTGGCAC"
Bsbv3_1_1R="TGCGCATTCAGGCTTAATGTTTGT"
Bsbv3_1Ralt="GCGCATTCAGGCTTAATGTTTGTTAC"
Bsbv4_8_1Falt="CCCACCAGCGTCCGGTAG"
Bsbv4_8_1Palt="TCTCGGCCTGTGATGCGCC"
Bsbv4_8_1R="GTCCCTGCGTTCGCAGAATG"
Bsbv4_9_1Falt="TCGATATCCTGGGCGCACC"
Bsbv4_9_1P="CGAGGATCGCATTGCCGGTGAC"
Bsbv4_9_1R="AGGGCAGGTCTATTCACGCG"

targets=(Bov2_2F Bov2_2P Bov2_2R Bov4_1Falt Bov4_1P Bov4_1Ralt Bmel2_1F Bmel2_1P Bmel2_1R Bs1_2F Bs1_2P Bs1_2R Bs2_3F Bs2_3P Bs2_3R Bs2_3Ralt Bsbv3_1_1F Bsbv3_1_1P Bsbv3_1_1R Bsbv3_1Ralt Bsbv4_8_1Falt Bsbv4_8_1Palt Bsbv4_8_1R Bsbv4_9_1Falt Bsbv4_9_1P Bsbv4_9_1R)

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

# clean columns.txt
# sed 's/Bmel2_1F  //g' | sed 's/Bmel2_1P  //g' | sed 's/Bmel2_1R  //g' | sed 's/Bov2_2F  //g' | sed 's/Bov2_2P  //g' | sed 's/Bov2_2R  //g' | sed 's/Bov4_1Falt  //g' | sed 's/Bov4_1P  //g' | sed 's/Bov4_1Ralt  //g' | sed 's/Bs1_2F  //g' | sed 's/Bs1_2P  //g' | sed 's/Bs1_2R  //g' | sed 's/Bs2_3F  //g' | sed 's/Bs2_3P  //g' | sed 's/Bs2_3R  //g' | sed 's/Bs2_3Ralt  //g' | sed 's/Bsbv3_1_1F  //g' | sed 's/Bsbv3_1_1P  //g' | sed 's/Bsbv3_1_1R  //g' | sed 's/Bsbv3_1Ralt  //g' | sed 's/Bsbv4_8_1Falt  //g' | sed 's/Bsbv4_8_1Palt  //g' | sed 's/Bsbv4_8_1R  //g' | sed 's/Bsbv4_9_1Falt  //g' | sed 's/Bsbv4_9_1P  //g' | sed 's/Bsbv4_9_1R  //g' 

#
#  Created by Stuber, Tod P - APHIS on 02/04/2015
#
