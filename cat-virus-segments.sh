#!/bin/sh

# cat all 8 virus segments

if [ -z $@ ]; then
        echo "No arguments given, Example: cat-virus-segments.sh all or cat-virus-segments.sh 3 5 8"
        echo ""
        exit 1
fi 

if [ $1 == all ]; then

for i in *fasta_H5N2.txt; do 
	type="all"
	echo $i
	dos2unix $i
	name=`echo $i | sed 's/_fasta_H5N2.txt//'`
	header=`head -1 $i | sed 's/>Seq1 //' | sed 's/segment.*//'`
	echo "name: $name"
	echo "header: $header"
	mkdir -p cat-${type}
	echo ">${name} $header" > ./cat-${type}/${name}_cat.fasta
	grep -v ">" $i | tr -d "\n" >> ./cat-${type}/${name}_cat.fasta
done

else

for i in *fasta_H5N2.txt; do
	echo $i
        dos2unix $i
        name=`echo $i | sed 's/_fasta_H5N2.txt//'`
        header=`head -1 $i | sed 's/>Seq1 //' | sed 's/segment.*//'`
        echo "name: $name"
        echo "header: $header"
        segments=`echo $@ | sed 's/ /_/g'`
        type="$segments"
	echo "segments: $segments"
	mkdir -p cat-${type}	
	echo ">${name} $header" > ./cat-${type}/${name}_cat_segments-${segments}.fasta
	
	for s in $@; do 
		echo $s
	grep -A 1 ">Seq${s}" $i | grep -v ">" | tr -d "\n" >> ./cat-${type}/${name}_cat_segments-${segments}.fasta
	done

done
fi

echo "samples are in cat directory"

