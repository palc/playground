#!/bin/sh

# Sed searches put into variables
tbNumberV='s/_.*//' #Remove all charaters at and beyond "_"
tbNumberW='s/\..*//' #Remove all charaters at and beyond "."

echo "Checking for duplicate VCFs."

for i in *; do
    getbase=`basename "$i"`
    number=`echo $getbase | sed $tbNumberV | sed $tbNumberW`
    echo $number >> list
done
duplist=`sort list | uniq -d`
rm list
dupNumberSize=`echo $duplist | wc | awk '{print $3}'`
if [ $dupNumberSize -gt 4 ]
then
    echo "These are duplicated VCFs."
    echo "Must remove duplication, and restart script."
    echo $duplist
    exit 1 # Error status
else
    echo "Good! No duplicate VCFs present"
fi

