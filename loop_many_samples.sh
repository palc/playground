#!/bin/sh

COUNTER=17
files=`ls *gz | wc -l`

root=`pwd`

while [ $files -gt 2 ]; do
    mkdir $COUNTER 
    mv `ls *gz | head -30` ${COUNTER}
    cd ${COUNTER}     
    packagefastqs.sh
    currentdir=`pwd`; for f in *; do cd $currentdir; echo $f; cd ./$f; processZips.sh H37Rv & done; wait 
    cd ${root}
    files=`ls *gz | wc -l`
    echo "files left to do: $files"
    let COUNTER=COUNTER+1
    echo "batch $COUNTER"
done

emailme
# tstuber 2016-11-30
