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

Indel-1-2="CGCTGCTGGAAGGCTTCTATGGCAGCGCGTGAGGTTGCGCCGATCTTGCCGTCATAATAGCCGAGCGCCTTGAGGTGCGTTTGCAGTTCC"
Indel-12-1="CCTGCGGCCGGTGAACAGGCCGCAGCCTTCGTGAATGGGGAATGGGGCTATCGCCGTTTTCTCAACCATTGCGCCGATCTGGCCGTGCGGGGTCGATGC"
Indel-14-1="TTTCATGCGCTCCATCTATATCAGCGTTGATGCGATCTCCGGCGTTCCAAACTACCTGCATGTCGTGGGGCTGACTGATCTGATGCGGTGTCTCATTATCGATGCGACATCGCACGATAGTATCCCGGCGCCGGAATACCC"
Indel-16-1="GAACAGTCAGACAAGGACCGCCATTTTTACATAAGTCATTGTAATTGCTTGGTAAAATAATGGTGCCCGGAGGCGGATTCGAACCACCGACACGCG"
Indel-24-1="GCACATTGAGCGCACTTTTCCATCACCGGCGTAAATTCATAGCCATTCCATTCATCGCCGCCAATCAGGATACGCTCGCGATAATCATTGAAATCACCACAAATCAGCCAGC"
Indel-26-1="CGGGAATGCCAGCATGACGTCCATGATGCGCATGGCGACCGTATCGACGCCACCGCCCATATAGCCGGAAAACACGCCGATGGCGATAACGAAGCCAACG"
Indel-28-2="ATCGAGATCCTCCACCTTGCGCATTTTTCGCCCTTGGCGGAAAATGCTTCCAACCACGTTCATGCATGATTGGTCACATAAGGGCGGAATCG"
Indel-6-2="GGCCTTTTCGATGAAGGAGCGGATCGATAATCTGGGCGCGGAAATGAGTTTTTCCGCTTCGCAGGATTCGGTTTCCGGCGGCGTCCGTATGCTGGCCGAAAATCGTGATGCCGTGACCG"

targets=(Indel-1-2 Indel-12-1 Indel-14-1 Indel-16-1 Indel-24-1 Indel-26-1 Indel-28-2 Indel-6-2)

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
