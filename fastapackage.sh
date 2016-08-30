#!/bin/sh

# Get rid of spaces
for i in *; do 
    name=`echo "$i" | sed 's/ /_/g'`
    mv "$i" $name
done

# Change extensions to .fasta"
for i in *.fa; do
    mv "$i" ${i%.fa}.fasta
done

for i in *.fna; do 
    mv "$i" ${i%.fna}.fasta
done

for i in *.fas; do
    mv "$i" ${i}ta
done

for i in *.fasta; do 
    mkdir ${i%.fasta}
    mv $i ${i%.fasta}
done

echo 'currentdir=`pwd`; for f in *; do cd $currentdir; echo $f; cd ./$f; sixteen_s.sh *fasta 16s/rpoB ;& done; wait; echo "DONE $PWD" > tempfile; cat tempfile | mutt -s "DONE" -- tod.p.stuber@aphis.usda.gov; rm tempfile'

# stuber 2016-05-24
