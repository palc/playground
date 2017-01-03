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

BAP6119="ACATTCCAGATAATACACCCG"
BAP6120="ATTGGAGCACCTAGTAACCC"
BAP6121="CTTAAAGTAACACTCGCTATTGC"
BAP6122="TTTGATTTCCCTTGGGATAGC"
BAP6123="TCCTTATCTGACATTGAAATCG"
BAP6124="CTAGACATCTGGTGGTTGCG"
BAP6125="TTTCCATAGATTAGCAATGCCG"
BAP6126="CTTTATTTGGTCTTTATATATACC"
BAP6127="CCTATATTTATATCTCCTCCCC"
BAP6128="CTAATATATAAACCATCCAACGC"
BAP6129="AGATTGCATGGCGAAATGGC"
BAP6130="CAATCCTCGTAAGACCCCC"
BAP6132="AATGAAGGTTTAAAAGAGATAGC"
BAP6133="GAGAGTTACAAAAATGATCGGC"
BAP6134="TCCTGGTTCATATATAGGTAGG"
BAP7039="AATATCTTTATAATTATACTCTCCC"
BAP7213="TGCAGGCGAGAGTTGATAAACCATC"
BAP7214="CAAAGATTGGTTCCAAATCTGAATGGA"
BAP7292="TCTTTATAATTATACTCTCCCAAGG"
BAP7293="AATGAAGGTTTAAAAGAGATAGCTGGAG"

targets=(BAP6119 BAP6120 BAP6121 BAP6122 BAP6123 BAP6124 BAP6125 BAP6126 BAP6127 BAP6128 BAP6129 BAP6130 BAP6132 BAP6133 BAP6134 BAP7039 BAP7213 BAP7214 BAP7292 BAP7293)

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
