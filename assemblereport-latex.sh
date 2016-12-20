#!/bin/sh

alias pause='read -p "$LINENO Enter"'

######
function blastcontig () {
grep -v "^$" $1 > file
h=`grep -c ">" file`
a=`grep -c ".*" file`
echo "h: $h a: $a"
singlelinetest=`expr $a / $h`
echo "singlelingtest: $singlelinetest"

if [[ $singlelinetest == 2 ]]; then
        echo "contig file ( $1 ) is formatted correctly"
        cat $1 > ${n}-contigsoriginal.fa
else
        # If contigs going in are not all on the same line then newlines are removed to put fasta on a single line.
        echo "contig file ( $1 ) was reformated placing sequence on a sigle line"
        #header names are changed
        #awk '{if ($0 ~ />/ ) {print ">"} else print $0}' $1 | tr -d "\n" | sed -e 's:>:\n>\n:g' | awk '{if ($0 ~  />/ ) {print ">contigs-" x++} else print $0}' | grep -v "^$" > ${n}-contigsoriginal.fa
        # header are not changed
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' $1 | grep -v '^$' > ${n}-contigsoriginal.fa
fi
rm file
###

# Remove everything behond the space.  Meant for ">name otherinfo" to just ">name"
sed 's/ .*$//' ${n}-contigsoriginal.fa > ${n}-contigs.fa

# Cut down contigs and this also should put all fasta contigs onto a single line.
#sed -e 's/$/=/g' ${n}-contigs.fa | tr -d "\n" | tr "=" "\n" | grep -v "^$" | awk '{ if ($0 ~ /^>/ ) { print $0 } else if (length($O) > 600 ) {print substr($0, 50, 500) } else {print $0}}' > ${n}-contigs3.fa
cat ${n}-contigs.fa > ${n}-contigs3.fa
echo "##### Blasting contig file #####"
echo "Start Time: `date`"

# removed "-word_size 11", didn't help in this case.  Caused more identification/contig
echo "short BLAST"
blastn -query ${n}-contigs3.fa -db /data/BLAST/db/nt -num_threads 40 -out ${n}-consensus-max1-nt.txt -max_target_seqs 1 -outfmt "6 saccver stitle"

awk 'BEGIN{OFS="\t"}{k=$1; a[k]++; b[k]=$0}; END{for (k in a) print a[k], b[k]}' ${n}-consensus-max1-nt.txt | sort -rnk1,1 > $n.pre.tex

echo "" >> ${mystart}/${n}.tex
echo "\begin{longtable}{ l | l | p{13cm} }" >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
echo "n & accession & identification \\\\" >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
tr "\t" "&" < $n.pre.tex | sed 's/&/ & /g' | sed 's:$: \\\\:' | sed 's/_/\\_/g' >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
echo "\end{longtable}" >> ${mystart}/${n}.tex
echo "\end{document}" >> ${mystart}/${n}.tex


rm $n.pre.tex
rm ${n}-contigs*fa

}

######   ######   START   ######   ######

mystart=`pwd`

mkdir ./zips
mv *fastq* ./zips

cd ./zips
# Place R1 Reads into variable
forReads=`ls | grep _R1`
echo "Forward Reads:  $forReads"
# Place R2 Reads into variable
revReads=`ls | grep _R2`
echo "Reverse Reads:  $revReads"

# Get name of isolate to use for naming output files
n=`echo $forReads | sed 's/[._].*//'`
echo "***Isolate naming convention:  $n"

#here-document
#latex preamble
cat << EOL > ${mystart}/${n}.tex
\documentclass[a4paper,11pt]{article}
\usepackage[margin=0.5in]{geometry}
\usepackage{graphicx}
\usepackage[table]{xcolor}
\usepackage{longtable}

\renewcommand{\thepage}{Appendix --  page \arabic{page}}

\begin{document}

\includegraphics[scale=0.2]{/home/tstuber/report_doc/usdalogo.png}


\today

\vspace{5mm}
\textbf{Whole Genome Sequencing Report:  ${n}} 

EOL

forcount=`zgrep -c "^+$" $forReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
revcount=`zgrep -c "^+$" $revReads | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'`
forsize=`ls -lh $forReads | awk '{print $5}'`
revsize=`ls -lh $revReads | awk '{print $5}'`

echo "" >> ${mystart}/${n}.tex
echo "\vspace{5mm}" >> ${mystart}/${n}.tex
echo "\textbf{File Stats}" >> ${mystart}/${n}.tex
echo "\vspace{2mm}" >> ${mystart}/${n}.tex
echo "" >> ${mystart}/${n}.tex

echo "" >> ${mystart}/${n}.tex
echo "\begin{tabular}{ l | p{7cm} | p{7cm} }" >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
echo "file name & $forReads & $revReads \\\\ " | sed "s/$n[._]//g" | sed 's/_/\\_/g' >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
echo "read count & $forcount & $revcount \\\\ " >> ${mystart}/${n}.tex
echo "file size & $forsize & $revsize \\\\ " >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
echo "\end{tabular}" >> ${mystart}/${n}.tex

# De Novo assembly using ABySS
abyss-pe name=${n}_abyss k=64 in="$forReads $revReads"

echo "" >> ${mystart}/${n}.tex
echo "\vspace{5mm}" >> ${mystart}/${n}.tex
echo "\textbf{Assembly}" >> ${mystart}/${n}.tex
echo "\vspace{2mm}" >> ${mystart}/${n}.tex
echo "" >> ${mystart}/${n}.tex

# Clean-up output
mkdir ../${n}_abyss-blast
cp ${n}_abyss-3.fa ../${n}_abyss-blast/${n}_abyss-3.fasta
cp ${n}_abyss-8.fa ../${n}_abyss-blast/${n}_abyss-8.fasta
cp ${n}_abyss-stats.tab ../${n}_abyss-blast

echo "" >> ${mystart}/${n}.tex
echo "\begin{tabular}{ l | l | l | l | l | l | l | p{7cm} }" >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex

awk 'BEGIN{OFS="\t"} {print $1, $2, $4, $6, $8, $9, $10, $11}' ${n}_abyss-stats.tab | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\:' | sed 's/_/\\_/g' | sed 's/name \\\\$/name \\\\ \\hline/' | sed "s/$n\\\_abyss-unitigs.fa/assembly/"  | sed "s/$n\\\_abyss-contigs.fa/assembly using paired ends/"  | sed "s/$n\\\_abyss-scaffolds.fa/contig scaffolds/" >> ${mystart}/${n}.tex
echo "\hline" >> ${mystart}/${n}.tex
echo "\end{tabular}" >> ${mystart}/${n}.tex

rm `ls ${n}_abyss* | grep -v "${n}_abyss-blast"`
rm coverage.hist

echo "" >> ${mystart}/${n}.tex
echo "\vspace{5mm}" >> ${mystart}/${n}.tex
echo "\textbf{Identification}" >> ${mystart}/${n}.tex
echo "\vspace{2mm}" >> ${mystart}/${n}.tex
echo "" >> ${mystart}/${n}.tex

#BLAST reads
cd ../${n}_abyss-blast

if [ -e ${n}_abyss-8.fasta ]; then
        blastcontig *-8.fasta
else
        blastcontig *-3.fasta
fi

cd ${mystart}
pdflatex ${mystart}/${n}.tex
rm ${n}.aux
rm ${n}.log

echo ""
echo "*** Done ***"
echo ""

# created 2015-10-26, stuber
