#!/bin/sh


for i in *.fastq*; do n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`; echo "n is : $n"; mkdir -p $n; mv $i $n/; done

echo 'currentdir=`pwd`; for f in *; do cd $currentdir; echo $f; cd ./$f; `PROGRAM_CALL` ;& done; wait; echo "DONE $PWD" > tempfile; cat tempfile | mutt -s "DONE" -- tod.p.stuber@aphis.usda.gov; rm tempfile'

# created 2015-06-29
