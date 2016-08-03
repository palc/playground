#!/bin/sh

#Clean Tag file
#Run script as: clean_tag.sh <tag_file.txt>

cat $1 | sed 's/\*//g' | sed 's/(/_/g' | sed 's/)/_/g' | sed 's/ /_/g' | sed 's/-_/_/g' | sed 's/_-/_/g' | sed 's/,/_/g' | sed 's#/#_#g' | sed 's#\\#_#g' | sed 's/__/_/g' | sed 's/__/_/g' | sed 's/__/_/g' | sed 's/-$//g' | sed 's/_$//g' |awk 'BEGIN {OFS="\t"}{gsub("_$","",$1)}1' > outfile.txt

rm $1
mv outfile.txt $1

#
#  Created by Stuber, Tod P - APHIS on 4/22/2013.
#
