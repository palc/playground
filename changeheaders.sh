#!/bin/sh

if [ -z $1 ]; then
	echo ""
	echo "	usage: changeheaders.sh <fasta file> <updated names - 2 column, tab delimited, updated name in column 2>"
	echo ""
	exit
fi

if [ -z $2 ]; then
	echo ""
	echo "  usage: changeheaders.sh <fasta file> <updated names - 2 column, tab delimited, updated name in column 2>"
	echo ""
	exit    
fi

fflag=
hflag=
while getopts 'fh' OPTION; do
case $OPTION in
	f) fflag=1
	;;
	h) hflag=1
	;;
esac
done
shift $(($OPTIND - 1))

if [ "$fflag" ]; then
	echo "changing names based on file name"

fi

if [ "$hflag" ]; then
	echo "changing names based on header name"


while read -r line; do
	name=`echo $line | sed 's/_.*//' | sed 's/>//'`
	newname=`grep "$name" names`
	string=`echo $line | grep -v '^>'`
	header=`echo $line | grep '^>'`
	#echo $line | sed 's/_.*//' | sed 's/>//'
	if [ -n $header ]; then
		name=`echo $line | sed 's/_.*//' | sed 's/>//'`
		newname=`grep "$name" names`
		if [ -n $newname ]; then
			echo "$newname"
		else
			echo "$header"
		fi
	fi

	if [ -n $string ]; then
		echo "$string"
	fi
done < $1

fi


#2016-01-19 tstuber
