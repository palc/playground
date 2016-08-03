#!/bin/sh

name=`echo $1 | sed 's/\..*//'`

echo '##fileformat=VCFv4.1' > $name-highqualitysnps.vcf                   
echo '#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1 rsIDs' | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' >> $name-highqualitysnps.vcf
grep -v '^#' $1 | awk 'BEGIN{OFS="\t"} $6 > 500 {print $1, $2, ".", $4, ".",  ".", ".", ".", ".", ".", "."}' >> $name-highqualitysnps.vcf 

echo "$name-highqualitysnps.vcf is created"
echo ""

# created 2015-07-04, tstuber
