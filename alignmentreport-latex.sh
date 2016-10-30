#!/bin/sh

# usage: working directory must have reference (fasta), and paired FASTQ read files.

alias pause='read -p "$LINENO Enter"'
picard='/usr/local/bin/picard-tools-1.141/picard.jar'

forReads=`ls | grep _R1`; sampleName=`echo $forReads | sed 's/[._].*//'`
ref=`ls | grep .fasta`; r=`echo $ref | sed 's/\..*//'`

root=`pwd`
echo "" > $root/${sampleName}.summaryfile
summaryfile="$root/${sampleName}.summaryfile"
mytex="$root/${sampleName}.tex"

#######################################################################################
#|||||||||||||||||||||| Function to make R graph of coverage ||||||||||||||||||||||||||
#######################################################################################
function plotR () {
cat > ./plotR.r << EOL
#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(scales)

arg <- commandArgs(trailingOnly=TRUE)

data <- read.csv(arg[1], header=FALSE, sep="\t")
names(data) <- c("Species", "position", "coverage")

pdf("myplot.pdf", width=20, height=6)

#ggplot(data, aes(x=position, y=log(coverage), colour=species, group=species)) + geom_point(size=2.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=10)))

#bp <- ggplot(data, aes(x=position, y=log10(coverage), colour=species, group=species)) + geom_point(size=2.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=10)))

#bp + scale_y_continuous(breaks=seq(0, 3.0, 0.5))

ggplot(data, aes(x=position, y=coverage, colour=Species, group=Species)) + geom_point(size=1.0) + ggtitle(arg[2]) + scale_colour_brewer(palette="Set1")+ theme_bw() + guides(colour = guide_legend(override.aes = list(size=8)))+ scale_y_log10() + theme(axis.title.x = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=0,face="plain"), axis.title.y = element_text(colour="grey20",size=20,angle=90,hjust=.5,vjust=.5,face="plain"), title = element_text(colour="grey20",size=20,angle=0,hjust=.5,vjust=.5,face="plain")) + ylab("Depth of Coverage") + xlab("Position")

dev.off()
EOL

chmod 755 ./plotR.r

./plotR.r $1 $2
rm ./plotR.r
}
#######################################################################################

#here-document
#latex preamble
cat << EOL > $mytex
\documentclass[a4paper,10pt]{article}
\usepackage[margin=0.5in]{geometry}
\usepackage{graphicx}
\usepackage[table]{xcolor}
\usepackage{floatrow}
\usepackage{float}
\floatsetup[table]{capposition=top}
\usepackage{caption}
\captionsetup{labelformat=empty,justification=justified,singlelinecheck=false}
\usepackage{helvet}
\renewcommand{\familydefault}{\sfdefault}
\usepackage{lastpage}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}

\cfoot{Appendix --  page \thepage\ of \pageref{LastPage}}

\begin{document}

\includegraphics[scale=0.2]{/home/tstuber/report_doc/usdalogo.png}


\today

\vspace{5mm}

\textbf{Alignment Report:  ${sampleName}}

\vspace{1mm}

\textbf{Reference: $r}

\vspace{5mm}

EOL
#######################################################################################

# flag -m will email just "M"e
# flag -b will turn off muliple for starts for "B"ug finding
# flag -k will run Kraken
# flag -e flag used when running script from idemail.sh
bflag=
mflag=
kflag=
eflag=
while getopts 'bmke' OPTION; do
	case $OPTION in
		b) bflag=1
		;;
		m) mflag=1
		;;
		k) kflag=1
		;;
		e) eflag=1
		;;
		?) echo "Invalid option: -$OPTARG" >&2
		;;
	esac
done
shift $(($OPTIND - 1))

# Grab reads and reference and place them in variables
ref=`ls | grep .fasta`
echo "Reference Input:  $ref"
if [ -z $ref  ]; then
	echo "No reference file in working directory.  Required format: referencename.fasta"
	exit 1
fi

forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
if [ -z $forReads ]; then
	echo "No forward read in working directory.  Required format: samplename*_R1*.fastq.gz"
    exit 1
fi

sampleName=`echo $forReads | sed 's/[._].*//'`
echo "" > $root/${sampleName}.summaryfile
summaryfile="$root/${sampleName}.summaryfile"
mytex="$root/${sampleName}.tex"

revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"
if [ -z $revReads ]; then
	echo "No reverse read in working directory.  Required format: samplename*_R2*.fastq.gz"
fi

if [ -f *R2* ]; then
	forFileSize=`ls -lh $forReads | awk '{print $5}'`
	revFileSize=`ls -lh $revReads | awk '{print $5}'`
	forCount=`zgrep -c '^+$' $forReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
	revCount=`zgrep -c '^+$' $revReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
	declare -i x=${forCount}
	declare -i y=${revCount}
	echo "forCount: $forCount"
else
	forFileSize=`ls -lh $forReads | awk '{print $5}'`
        forCount=`zgrep -c '^+$' $forReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
        declare -i x=${forCount}
        echo "forCount: $forCount"
fi

# START KRAKEN
if [ "$kflag" ]; then
	echo "`date` -------> Kraken has started"
	krakenDatabase="/home/shared/databases/kraken/std/"
	echo "Kraken database selected is: $krakenDatabase"
	kraken --db ${krakenDatabase} --threads ${NR_CPUS} --paired *fastq* > $sampleName-output.txt && kraken-report --db ${krakenDatabase} $sampleName-output.txt > $sampleName-kraken_report.txt

	echo "*** Krona transforming Kraken output to graph"

	# Run Krona
	cut -f2,3 $sampleName-output.txt > $sampleName-kronaInput.txt;
	/usr/local/bin/ktImportTaxonomy $sampleName-kronaInput.txt;
	mv taxonomy.krona.html $sampleName-Krona_identification_graphic.html;
	mv taxonomy.krona.html.files $sampleName-taxonomy.krona.html.files

	# Set variables and paths
	output=`ls *-output.txt`
	report=`ls *kraken_report.txt`
	
	printf "%s, %s file size, %s reads\n" ${forReads} ${forFileSize} ${forCount}
	printf "%s, %s file size, %s reads\n" ${revReads} ${revFileSize} ${revCount}

	#Section of results summary that calculates number of reads per type of organism (ex: ssRNA virus)
	echo "Summary of Kraken Findings"
	cRead=`grep -c "^C" $output`
	uRead=`grep -c "^U" $output`
	let allReads=cRead+uRead
	echo "allReads: $allReads"

echo "`date` -------> Kraken has completed"
# END KRAKEN
fi

# COLLECT READ STATS
echo "`date` -------> Collecting read stats"
echo "R1 file size: ${forFileSize}, read count: $forCount" >> $summaryfile
echo "R2 file size: ${revFileSize}, read count: $revCount" >> $summaryfile

forsize=`ls -lh $forReads | awk '{print $5}'`
revsize=`ls -lh $revReads | awk '{print $5}'`

echo "" >> ${mytex}.filestats
echo "\vspace{5mm}" >> ${mytex}.filestats
echo "" >> ${mytex}.filestats

echo "\begin{table}[H]" >> ${mytex}.filestats
echo "\begin{tabular}{ l | p{7cm} | p{7cm} }" >> ${mytex}.filestats
echo "\hline" >> ${mytex}.filestats
echo "$sampleName & R1 file & R2 file \\\\ " >> ${mytex}.filestats
echo "\hline \hline" >> ${mytex}.filestats
echo "read count & $forCount & $revCount \\\\ " >> ${mytex}.filestats
echo "\hline" >> ${mytex}.filestats
echo "file size & $forsize & $revsize \\\\ " >> ${mytex}.filestats
echo "\hline" >> ${mytex}.filestats
echo "\end{tabular}" >> ${mytex}.filestats
echo "\caption{\textbf{File stats}}" >> ${mytex}.filestats
echo "\end{table}" >> ${mytex}.filestats

#START ALIGNMENT
echo "`date` -------> Alignment Started"
r=`echo $ref | sed 's/\..*//'`
n=`echo $revReads | sed 's/_.*//' | sed 's/\..*//'`

samtools faidx $ref
java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict

if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
	echo "Index and dict are present, continue script"
else
	sleep 5
	echo "Either index or dict for reference is missing, try making again"
	samtools faidx $ref
	java -Xmx4g -jar ${picard} CreateSequenceDictionary REFERENCE=${ref} OUTPUT=${r}.dict
	if [ -s ${ref}.fai ] && [ -s ${r}.dict ]; then
		read -p "--> Script has been paused.  Must fix.  No reference index and/or dict file present. Press Enter to continue.  Line $LINENO"
	fi
fi

# See echo comments
echo "***bwa index $r"
bwa index $ref

echo "***Making Sam file"
if [ -f *R2* ]; then
     echo "Read 2 file is present"
     # Using -B 1 achieved a slightly higher reference coverage, M allows compatibility with Picard
     bwa mem -M -B 1 -t 5 -R @RG"\t"ID:"${sampleName}""\t"PL:ILLUMINA"\t"PU:"${sampleName}"_RG1_UNIT1"\t"LB:"${sampleName}"_LIB1"\t"SM:"${sampleName}" $ref $forReads $revReads > ${sampleName}.sam
else
     echo "Read 1 file is present"
     bwa mem -M -t 5 -R @RG"\t"ID:"${sampleName}""\t"PL:ILLUMINA"\t"PU:"${sampleName}"_RG1_UNIT1"\t"LB:"${sampleName}"_LIB1"\t"SM:"${sampleName}" $ref $forReads > ${sampleName}.sam
fi

#F4 just keeps the mapped reads in the bam file
samtools view -bh -F4 -T $ref ${sampleName}.sam > ${sampleName}.mapped.bam
echo "Sorting Bam"
samtools sort ${sampleName}.mapped.bam -o ${sampleName}.mapped.sorted.bam
echo "****Indexing Bam"
samtools index ${sampleName}.mapped.sorted.bam

rm ${sampleName}.sam
rm ${sampleName}.mapped.bam

echo "***Marking Duplicates"
java -Xmx4g -jar  ${picard} MarkDuplicates INPUT=${sampleName}.mapped.sorted.bam OUTPUT=${sampleName}.dup.bam METRICS_FILE=${sampleName}.FilteredReads.xls ASSUME_SORTED=true REMOVE_DUPLICATES=true

echo "***Index ${sampleName}.dup.bam"
samtools index ${sampleName}.dup.bam

#Mean Quality by Cycle
echo "***Mean Quality by Cycle"
java -Xmx4g -jar ${picard} CollectMultipleMetrics REFERENCE_SEQUENCE=$ref INPUT=${sampleName}.dup.bam OUTPUT=$n.Quality_by_cycle PROGRAM=MeanQualityByCycle ASSUME_SORTED=true

#Collect Alignment Summary Metrics
echo "***Collect Alignment Summary Metrics"
java -Xmx4g -jar ${picard} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=$ref INPUT=${sampleName}.dup.bam OUTPUT=$n.AlignmentMetrics ASSUME_SORTED=true

echo "***Bamtools is getting coverage"
bamtools coverage -in ${sampleName}.dup.bam | awk -v x=${sampleName} 'BEGIN{OFS="\t"}{print x, $2, $3}' >> ${sampleName}-coveragefile
plotR ${sampleName}-coveragefile ${sampleName}

sed -n 7,8p ${sampleName}.FilteredReads.xls | awk '{print $2}' >> ${summaryfile}
sed -n 7,8p ${sampleName}.FilteredReads.xls | awk '{print $3}' >> ${summaryfile}
sed -n 7,8p ${sampleName}.FilteredReads.xls | awk '{print $8}' >> ${summaryfile}
readcount=`sed -n 8p ${sampleName}.FilteredReads.xls | awk '{print $3}'`

echo 'Mean_Insert_Size  Standard_Deviation:' >> ${summaryfile}
awk 'BEGIN {OFS="\t"} { print $5,$6 }' ${sampleName}.Quality_by_cycle.insert_size_metrics | awk 'FNR == 8 {print $0}' >> ${summaryfile}

echo 'Mean_Read_Length:' >> ${summaryfile}
awk 'BEGIN {OFS="\t"} { print $16 }' ${sampleName}.AlignmentMetrics | awk 'FNR == 10 {print $0}' >> ${summaryfile}

echo "" >> ${summaryfile}
echo "Alignment stats (reference guided):" >> ${summaryfile}
printf "%-45s %10s %11s %10s\n" "reference used" "genome size" "percent cov" "ave depth" >> ${summaryfile}

echo "***Bamtools running"
aveCoverage=`awk '{sum+=$3} END { print sum/NR"X"}' ${sampleName}-coveragefile`
echo "Average depth of coverage: $aveCoverage"

#genome coverage
#Length of reference
countNTs=`grep -v ">" $ref | tr -d "\n" | wc | awk '{print $3}'`
covCount=`awk '{ if ($3 != 0) count++ } END { print count }' ${sampleName}-coveragefile`
echo "covCount $covCount"

declare -i x=${covCount}
declare -i y=${countNTs}

echo "----------------------> ${orgref} covCount: $x"
echo "----------------------> ${orgref} countNTs: $y"

#Percent of reference with coverage
perc=`awk -v x=$x -v y=$y 'BEGIN { print(x/y)*100}'`
echo "perc: $perc"

LC_NUMERIC=en_US
printf "%-45s %'10d %11.2f%% %s\n" ${sampleName} $countNTs $perc $aveCoverage >> ${summaryfile}.pre

cat ${summaryfile}.pre >> ${summaryfile}

echo "\begin{table}[H]" >> ${mytex}.alignmentstats
echo "\begin{tabular}{ l | l | l | l }" >> ${mytex}.alignmentstats
echo "\hline" >> ${mytex}.alignmentstats
echo "reference used & genome size & percent cov & ave depth \\\\" >> ${mytex}.alignmentstats
echo "\hline" >> ${mytex}.alignmentstats
echo "\hline" >> ${mytex}.alignmentstats
awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4}' ${summaryfile}.pre | sort -k1,1 | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\ \\hline:' | sed 's/[_]/\\_/g' | sed 's/[%]/\\%/g' >> ${mytex}.alignmentstats
echo "\end{tabular}" >> ${mytex}.alignmentstats
echo "\caption{\textbf{Alignment stats (reference guided)}}" >> ${mytex}.alignmentstats
echo "\end{table}" >> ${mytex}.alignmentstats

java -Xmx4g -jar ${gatk} -R $ref -T UnifiedGenotyper -I ${sampleName}.dup.bam -o ${sampleName}.UG.vcf -nct 8

# make reference guided contig
java -jar ${gatk} -T FastaAlternateReferenceMaker -R $ref -o ${sampleName}.readreference.fasta -V ${sampleName}.UG.vcf

echo "SNP in ${sampleName}.UG.vcf:" >> ${summaryfile}
egrep -v "#" ${sampleName}.UG.vcf | grep -c ".*" >> ${summaryfile}
lsnps=`egrep -v "#" ${sampleName}.UG.vcf | grep -c ".*"`

echo "SNPs of AC2 and QUAL > 300:" >> ${summaryfile}
egrep -v "#" ${sampleName}.UG.vcf | egrep "AC=2" | awk '$6 > 300' | grep -c ".*" >> ${summaryfile}
hsnps=`egrep -v "#" ${sampleName}.UG.vcf | egrep "AC=2" | awk '$6 > 300' | grep -c ".*"`

echo "\begin{table}[H]" >> ${mytex}.snpstats
echo "\begin{tabular}{ l | l }" >> ${mytex}.snpstats
echo "\hline" >> ${mytex}.snpstats
echo "All SNPs & High Quality SNPs \\\\ " | sed "s/$n[._]//g" | sed 's/_/\\_/g' >> ${mytex}.snpstats
echo "\hline \hline" >> ${mytex}.snpstats
echo "$lsnps & $hsnps \\\\ " >> ${mytex}.snpstats
echo "\hline" >> ${mytex}.snpstats
echo "\end{tabular}" >> ${mytex}.snpstats
echo "\caption{\textbf{SNP stats}}" >> ${mytex}.snpstats
echo "\end{table}" >> ${mytex}.snpstats

# Complete Latex document
cat ${mytex}.filestats >> $mytex

echo "\vspace{5mm}" >> $mytex
echo "" >> $mytex

echo "\vspace{5mm}" >> $mytex
echo "" >> $mytex

echo "\begin{figure}[H]" >> $mytex
echo "\begin{flushleft}" >> $mytex
echo "\textbf{Coverage Graph}\par\medskip" >> $mytex
echo "\end{flushleft}" >> $mytex

echo "\includegraphics[width=450pt]{myplot.pdf}" >> $mytex
echo "" >> $mytex
echo "\end{figure}" >> $mytex
echo "" >> $mytex

cat ${mytex}.alignmentstats >> $mytex

cat ${mytex}.snpstats >> $mytex

echo "\end{document}" >> $mytex
echo "" >> $mytex

pdflatex $mytex
pdflatex $mytex




















































# created 2015-01-26 stuber
