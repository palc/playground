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

p16M1="ACGCCTTCCAGCGGTGCATCATGTTTTTCA"
p16M2="ACCGCAACCGCGCTTCCGCCGATGTGAAGA"
p16M3="CTTTAATCTGAAAAGAGATGAATATGAGGC"
p16M4="TCACGACCGGCGGACGCGAAACCATTGCGG"
p16M6="GCTGAGCGGGCCTTGCAGGTGGCCAAGCAG"
p16M7="ATCACGATCGCCGTGCGGTAGCCCAATGCC"
p16M9="CAGCCCCAACCCCATCAACCATAAAATAGA"
p16M0="TCCCGCCGCCATGCCGCCGAAAGTCGCCGT"
p63B1="CCTGTTTAAAAGAATCGTCGGAACCGCTCT"
p63B2="CGCTTCGCGGCTGCGGGCGGGCGGGTAACG"
p63B4="AATTCACGGAACGTCACAGATTCAAAACG"
p63B5="GGTGCCGATGTTTTCGGCAGATGCCCAAGG"
p63B6="TCGAGCCAAAAGTGCGAAGCGGTTTTGCGT"
p63B7="GCCAAAATCCTCGGCGGCCATTTTTGCCGC"
p63B8="CCGAACCTGAGCGCGGCACTCCTCCTCCCA"
p63B9="GATAGTTTGAACCAGTCATTATTACCTCCA"
p63B0="GCGATACGCCGCTCGGACCCGAAAAGATTG"
Ether1="GCCCCTGCGCCAACGCCTGCTCAAGCGCGA"
Ether2="CGAAATCGTGGTGAAGGACGGGACCGAACC"
Ether4="ATCCAGGAGGAAGAACGGCAGCAGCGGTTA"
Ether6="TATCTCGATGTATGGGCCGATATTTTCCCT"
Ether7="GATATAGGCAGGCACGAAGCCGGGGAAGGC"
Ether8="GCATAGCGGCTCGGCCAATATCGTCACCGA"
Ether9="GGCACCGTTTGGAAACACCGACCGTCGCGT"
Ether0="CTCGGCTGGTTTTCGCTGCCGGTCGCCATC"

targets=(p16M1 p16M2 p16M3 p16M4 p16M6 p16M7 p16M9 p16M0 p63B1 p63B2 p63B4 p63B5 p63B6 p63B7 p63B8 p63B9 p63B0 Ether1 Ether2 Ether4 Ether6 Ether7 Ether8 Ether9 Ether0)

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
