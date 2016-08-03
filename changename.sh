#!/bin/sh


for i in *.fastq.gz; do

new=""
old=`echo $i | sed 's/_.*//'`
new=`awk -v x=$old ' $1 == x {print $2}' namelist`


echo "old name: $old"
echo "new name: $new"

fixname=`echo $i | sed "s/$old/$new/"`
echo "Fixed name: $fixname"

read -p "Press enter if okay, Prese ctrl-c if not"

mv $i $fixname

done



