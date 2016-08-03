#!/bin/sh

root=`pwd`
filelocation="/bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script1"
file=$1
mkdir -p ${root}/snp_script_output

for i in `cat ${file}`; do 
    
    echo  "working on $i"
    mkdir ${root}/${i}
    
    cp ${filelocation}/${i}*/zips/*_R1*gz ${root}/${i}
    
    cd ${root}/${i}
    forReads=$(ls | grep _R1)
    echo "Forward Reads:  $forReads"
    pigz -d $forReads
    ../scripts/gil2.sh ${forReads%.gz}
    cp SNP-${forReads%.gz}.txt ${root}/snp_script_output/ 
    cd ${root}
    rm -r ${root}/${i}
done

echo "script has completed"

# 2016-03-31
