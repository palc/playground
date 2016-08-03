#!/bin/sh


name=`echo $1 | sed 's/[._]//'`

if [[ -z $1 ]]; then
    printf "\n\tsupply a single argument\n"
    printf "\targument must be text file\n\n"
    exit 1
fi

enscript $1 -B -j -r -f "Courier5" -o - | ps2pdf - $name-txt_to.pdf


# stuber 2016-06-11
