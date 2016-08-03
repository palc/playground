#!/bin/sh

# Download xml at: http://www.ncbi.nlm.nih.gov
# Select assemblies
# "Send to:" File as xml

# Check if the last 4 characters of first argument end in ".xml"
if [ ${1: -4} != ".xml" ]; then
    echo "Script takes on argument"
    echo "Argument must have .xml extenstion"
    exit 1
fi


grep 'ftp://ftp.ncbi.nlm.nih.gov' $1 | sed 's/.*\(ftp:\/\/ftp.ncbi.nlm.nih.gov.*\)&lt;\/FtpPath.*/\1/' | sed 's:\(.*\)/\(.*\):\1/\2//\2_genomic.fna.gz:' > temp_list


printf "\n\nThere are `wc -l temp_list` files to download"
printf "\nWould you like to continue?\n"
read -p "Press Enter to continue, or ctrl-c to stop"
while read l; do wget $l; done < temp_list

sleep 5
gunzip *gz

printf "Files downloaded:\n"
grep -m 1 ">" *fna

for i in *fna; do 
    name=`grep -m 1 ">" $i | sed 's/\(.\{50\}\).*/\1/' | sed 's/[. ]/_/g' | sed 's/>//' | sed 's/__/_/' | sed 's/_$//' `
    mv "$i" ${name}.fasta
done

printf "\nFiles have been renamed:\n"
ls -lh *fasta

printf "\n\n\nScript has completed downloading and renaming files.\n\n"


# stuber 2016-05-25
