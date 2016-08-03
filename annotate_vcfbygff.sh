#!/bin/sh

if [[ -z $1 || -z $2 ]]; then
	echo "usage: ~$ annotate_vcfbygff.sh <in gff> <in vcf>"
	exit
fi

awk '$3 == "gene" {print $0}' $1 | awk '{print $4, $5, $9}' > list.genes

while read l; do
	echo $l | awk '{for(i=$1;i<=$2;i++) print i, $3}'
done < list.genes | sed -e 's/\([0-9]*\).*;\(Name=.*\);gbkey.*\(gene_biotype=.*\);\(locus_tag=.*\)/\1   \2;\3;\4/' > expand.gene

#Split header lines from position calls
grep "#" $2 > header
grep -v "#" $2 > body

#http://unix.stackexchange.com/questions/113898/how-to-merge-two-files-based-on-the-matching-of-two-columns

awk 'BEGIN{OFS="\t"}NR==FNR {h[$1] = $2; next} {print $1,$2,h[$2],$4,$5,$6,$7,$8,$9,$10}' expand.gene body | awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1' > body2

cat header body2 > ${2%.vcf}-annotated.vcf

rm body
rm body2
rm header
rm expand.gene
rm list.genes
