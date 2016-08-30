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

Bmel4_2F="GAACGTTGAATTCGGCCTTATCGG"
Bmel4_2P="CCTTGTCAGGCGGCATGAAGCAG"
Bmel4_2R="ATAGGCGCGCGCGATGG"
Bsbv1_5_1F="AGTTCCTCTCGCCGCTTTCA"
Bsbv1_5_1P="CGAAAACCGCGCGCGCCTG"
Bsbv1_5_1R="AAAGCCGCCGCGCAATG"
Bsbv3_2_1F="GATATCAACAAGGACAGCGAAACCT"
Bsbv3_2_1altF="TCAACAAGGACAGCGAAACCTAC"
Bsbv3_2_1P="CAGGAAGCCCCGTAAGCGCCTT"
Bsbv3_2_1R="GCAAGACCGCTTTGCATTCTG"
Bsbv3_3_1F="AATGAACGGAAGCGGGCAAC"
Bsbv3_3_1P="CGAGCGCGCCATCAGAAACGC"
Bsbv3_3_1R="GCCGATTCATAGATTCGCTGAGC"
Bsbv3_3_2F="CCATCAGAAACGCATTTGCTCAG"
Bsbv3_3_2P="CGCCCATGAGCGCGCCATG"
Bsbv3_3_2R="TTTCAGTTTTTCGCGGCGCT"
Bsbv4_6_2altF="GTCCTCGCGCAAACTGTCAA"
Bsbv4_6_2altP="GCACCGGCGCTTCGATTGGC"
Bsbv4_6_2altR="CGCCATTCTCCAGACGCTG"
Bsbv4_9_1Falt="TCGATATCCTGGGCGCACC"
Bsbv4_9_1P="CGAGGATCGCATTGCCGGTGAC"
Bsbv4_9_1R="AGGGCAGGTCTATTCACGCG"

targets=(Bmel4_2F Bmel4_2P Bmel4_2R Bsbv1_5_1F Bsbv1_5_1P Bsbv1_5_1R Bsbv3_2_1F Bsbv3_2_1altF Bsbv3_2_1P Bsbv3_2_1R Bsbv3_3_1F Bsbv3_3_1P Bsbv3_3_1R Bsbv3_3_2F Bsbv3_3_2P Bsbv3_3_2R Bsbv4_6_2altF Bsbv4_6_2altP Bsbv4_6_2altR Bsbv4_9_1Falt Bsbv4_9_1P Bsbv4_9_1R)

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
