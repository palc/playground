#!/bin/sh

alias pause='read -p "$LINENO Enter"'

###########################################################################################################
#||||||||||||||||||||||||||||||||| FUNCTION:  Finds contig and cuts to size  |||||||||||||||||||||||||
###########################################################################################################
function find_contig_cut () {

# Create "here-document" to prevent a dependent file.
cat >./find_contig_cut.py <<EOL
#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO
from sys import argv

passing_contigs = []
fasta_file = sys.argv[1]
header = sys.argv[2]

print (fasta_file)
print (header)

records = SeqIO.parse(fasta_file, "fasta")
for record in records:
     if record.name == header:
        if $sstart > $send:
            record = record[$send:$sstart]
            passing_contigs.append(record.reverse_complement())
        else:
            passing_contigs.append(record[$sstart:$send])

SeqIO.write(passing_contigs, "./finding$counter.fasta", "fasta")


EOL

chmod 755 ./find_contig_cut.py

echo $assembled_genome_file

./find_contig_cut.py "$assembled_genome_file" "$header" 

rm ./find_contig_cut.py

}

if [[ $2 == 16s ]]; then 
    sixteen_s_seed="/home/shared/myseq16s.fasta"
    BLASTdatabase="/data/BLAST/db/16SMicrobial"
    #BLASTdatabase="/data/BLAST/db/nt"
    #BLASTdatabase="/data/BLAST/db/other_genomic"
elif [[ $2 == rpoB ]]; then
    sixteen_s_seed="/home/shared/myrpoB1.fasta"
    BLASTdatabase="/data/BLAST/db/nt"
    #BLASTdatabase="/data/BLAST/db/other_genomic"
    #BLASTdatabase="/data/BLAST/db/Representative_Genomes"
elif [[ $2 == 16slpsn ]]; then
    sixteen_s_seed="/home/shared/myseq16s.fasta"
    BLASTdatabase="/data/BLAST/db/myco_16s_db/16slpsn" 
else 
    printf "\nProvide a proper #2 argument: either 16s, rpoB, or 16slpsn\n\n"
    exit 1
fi 



assembled_genome_file=$1
root=`pwd`
printf "\n File for seeding: $sixteen_s_seed \n\n"
printf "BLAST database: $BLASTdatabase\n\n"


name=`echo $assembled_genome_file | sed -e 's/_/-/g' -e 's/[.].*//'`

# Test if sequence is sanger or wgs
# Determine via file size

file_size=`ls -l $assembled_genome_file | awk '{print $5}'`
if [[ $file_size -lt 3000 ]]; then
    echo "file is sanger"
    file_type="Sanger"
else
    echo "file is wgs"
    file_type="WGS"
fi

function initial_blast () {
# Make input sequence (sanger or wgs) into BLAST database
# BLAST a known 16s/rpoB/etc against the database made
# Output findings
makeblastdb -in $assembled_genome_file -dbtype nucl -out $name.makeblast
blastn -query $sixteen_s_seed -db ./${name}.makeblast -num_threads 40 -out rpoB_found -outfmt "6 sseqid slen sstart send length"
}

initial_blast

if [[ ! -s rpoB_found ]]; then
    echo "File is empty, testing other rpoB"
    sixteen_s_seed="/home/shared/myrpoB2.fasta"
    initial_blast
fi


counter=1
# For each sequence that was found matching 16s/rpoB/etc.
while read l; do 
    cd ${root}
    echo "line: $l"
    header=$(echo $l | awk '{print $1}')
    printf "Found header: %s\n" $header
    sstart=$(echo $l | awk '{print $3}')
    printf "position sstart: %s\n" $sstart 
    send=$(echo $l | awk '{print $4}')
    printf "position send: %s\n" $send

    # Use python here-doc to extract 16s/rpoB/etc. sequence
    # Output is finding$counter.fasta
    find_contig_cut
    # Move the found seuence to new directory to isolate next steps in loop
    mkdir ./finding${counter}
    mv ./finding${counter}.fasta ./finding${counter}   
    cd ./finding${counter}
    
    # Output format as table file
    blastn -query ./finding${counter}.fasta -db $BLASTdatabase -word_size 11 -num_threads 40 -out table_blast${counter} -outfmt '6 sacc qlen slen gaps length qcovs evalue bitscore pident mismatch stitle' -num_alignments 15
    avium_check=$(head -4 table_blast${counter} | awk 'BEGIN{FS="\t"} {print $11}' | awk '{print $2}' | grep "avium" | sort -u) 
    if [[ $avium_check == avium ]] && [[ $2 == rpoB ]]; then
        printf "\n  avium found \n  avium specific BLAST search preformed\n\n"
        blastn -query ./finding${counter}.fasta -db /data/BLAST/avium_db/avium_db -out table_avium${counter} -outfmt '6 sacc qlen slen gaps length qcovs evalue bitscore pident mismatch stitle'
    awk 'BEGIN{FS=OFS="\t"}{gsub(/gi.*gb\|/, "", $1);print}' table_avium${counter} | awk 'BEGIN{FS=OFS="\t"}{gsub(/gi.*\|/, "", $11);print}' | sed 's/|//g' > table_avium${counter}.temp; mv table_avium${counter}.temp table_avium${counter} 
    fi

    # Ouput format as fasta file
    blastn -query ./finding${counter}.fasta -db $BLASTdatabase -num_threads 40 -out fasta_blast${counter} -outfmt '6 sacc stitle mismatch sseq' -num_alignments 50
    # Sort the table by mismatches
    sort -nk10,10 table_blast${counter} > table_blast.temp; mv table_blast.temp table_blast${counter}
    # Table header
    printf "sacc\tqlen\tslen\tgaps\tlength\tqcovs\tevalue\tbitscore\tpident\tmismatch\tstitle\n" > header
    # Table formated and ready for awk to latex 
    cat header table_blast${counter} > table_formated_blast${counter}
    cat header table_avium${counter} > table_formated_avium${counter}
    # Format the fasta_blast as a properly formated fasta file
    # Add limit characters to 50 with still keeping mismatch number
    sort -uk1,1 < fasta_blast${counter} | sed 's/\(.*\)\(\t[A-Z-][A-Z-][A-Z-][A-Z-][A-Z-][A-Z-][A-Z-][A-Z-][A-Z-][A-Z-][A-Z-].*\)/\1 break \2/' | sed -e 's/\t/_/g' -e 's/ /_/g' -e 's/[.]/_/g' -e 's/__/-/g' | sed -e 's/^/>/' -e 's/_break-/\n/' | sed ' /^>/ s/Mycobacterium/M/g' | sed 's/\(.\{50\}\).*_\([0-9]*\)/\1_\2/' | sed -e 's/__/_/g' -e 's/_$//' -e 's/[()]//g' > formated_fasta_blast${counter}.fasta 
    
    sed "s/>.*/>${name}-subject${counter}/" finding${counter}.fasta > ${name}-subject${counter}.fasta
    
    cat ${name}-subject${counter}.fasta >> formated_fasta_blast${counter}.fasta

    
#clustalw2 -OUTFILE=alignment${counter}.fasta -OUTPUT=FASTA -INFILE=formated_fasta_blast${counter}.fasta  -STATS=stats${counter}-clustalw2.txt -SEED=100 -KIMURA -BOOTLABELS=node

    # Find fasta with most mismatches
    highest_mismatch_number=`grep ">" formated_fasta_blast${counter}.fasta | sed -e 's/^.*_\(.*$\)/\1/' | sort -n | grep -v "subject" | tail -1`
    highest_mismatch_name=`grep ">" formated_fasta_blast${counter}.fasta | sed 's/..........\(.*\)/\1/' | grep "$highest_mismatch_number\$" | tail -l`
    highest_mismatch_fullname=`grep ">" formated_fasta_blast${counter}.fasta | grep "$highest_mismatch_name" | sed 's/>//' | tail -1`

    echo "highest_mismatch_number: $highest_mismatch_number" 
    echo "highest_mismatch_name: $highest_mismatch_name"
    echo "highest_mismatch_fullname: $highest_mismatch_fullname"

    sed -i 's/$/_______________/' formated_fasta_blast${counter}.fasta

    mafft --auto formated_fasta_blast${counter}.fasta > alignment${counter}.fasta
    #kalign -input formated_fasta_blast${counter}.fasta -output alignment${counter}.fasta -format fasta
    #remove illegal characters
    sed -i 's/[:;,]/_/g' alignment${counter}.fasta

    #raxmlHPC-SSE3 -s alignment${counter}.fasta -n raxml-out.tre -m GTRCAT -p 12345
    raxmlHPC-SSE3 -f a -s alignment${counter}.fasta -p 12345 -x 12345 -# 100 -m GTRCAT -n raxml-out${counter}.tre 
    #raxmlHPC-SSE3 -f b -t ref -z tree -m GTRCAT -s alg -n ${d}
    
    nw_reroot RAxML_bipartitions.raxml-out${counter}.tre ${highest_mismatch_fullname}_______________ > rooted_mytree${counter}.tre
    
    printf "fill:red;font-style:italic L ${name}-subject${counter}_______________" > mycss.css
#nw_display -s -S -w 1500 -v 20 -b 'opacity:0' -l 'font-size:16;font-family:serif;font-style:italic' -d 'stroke-width:1;stroke:blue' -c 'mycss.css'

    # "&" showing up in accession number causing svg creation to error
    sed -i 's/&//g' rooted_mytree${counter}.tre
    cat rooted_mytree${counter}.tre | nw_display -s -S -w 1300 -t -v 30 -i 'opacity:0' -b 'opacity:0' -l 'font-size:18;font-family:serif;font-style:italic' -d 'stroke-width:1;stroke:blue' -c 'mycss.css' - > tree${counter}.svg

    inkscape -f tree${counter}.svg -A tree${counter}.pdf

    printf "Date ran: `date`\n"
    printf "Sample name: $name\n"
    printf "Sample type: $file_type\n\n"

    enscript table_formated_blast${counter} -B -j -r -f "Courier5" -o - | ps2pdf - table_blast${counter}.pdf

#here-document
#latex preamble
cat << EOL > ${name}_${counter}.tex 
\documentclass[a4paper,11pt]{article}
\usepackage[left=2cm,right=2cm,top=3cm,bottom=3cm,headsep=\dimexpr3cm-72pt\relax,headheight=72pt]{geometry}
\usepackage{graphicx}
\usepackage[table]{xcolor}
\usepackage{caption}
\usepackage{helvet}
\usepackage{seqsplit}
\usepackage{fancyhdr}
\usepackage{lastpage}
\pagestyle{fancy}
\thispagestyle{fancy}

\renewcommand{\headrulewidth}{1pt}

\lhead{\includegraphics[scale=0.2]{/home/tstuber/report_doc/usdalogo.png}}

\cfoot{Appendix --  page \thepage\ of \pageref{LastPage}}

\begin{document}

\vspace{15mm}

\today

\vspace{5mm}

\textbf{${2} Report:  ${name}}

\vspace{5mm} 

\textbf{Source: ${file_type}}

\vspace{5mm} 

EOL

    ###
    echo "" >> ${name}_${counter}.tex
    echo "\normalsize{${2} sequence for ${name}}" >> ${name}_${counter}.tex
    echo "" >> ${name}_${counter}.tex
    echo "\tiny{\seqsplit{$(grep -v ">" ./finding${counter}.fasta | tr -d "\n")}}" >> ${name}_${counter}.tex
   
    echo "\begin{table}[ht]" >> ${name}_${counter}.tex
    echo "\caption*{Top 15 hits sort by mismatches}" >> ${name}_${counter}.tex
    echo "\centering" >> ${name}_${counter}.tex
    echo "\resizebox{\textwidth}{!}{\begin{tabular}{ | l | l | l | l | l | l | l | l | l | l | l | }" >> ${name}_${counter}.tex
    echo "\hline" >> ${name}_${counter}.tex
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' table_formated_blast${counter} | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\:' | sed 's/_/\\_/g' | sed 's/name \\\\$/name \\\\ \\hline/' >> ${name}_${counter}.tex
    echo "\hline" >> ${name}_${counter}.tex
    echo "\end{tabular}}" >> ${name}_${counter}.tex
    echo "\end{table}" >> ${name}_${counter}.tex

    if [[ $avium_check == avium ]] && [[ $2 == rpoB ]]  ; then
        echo "\begin{table}[ht]" >> ${name}_${counter}.tex
        echo "\caption*{M. avium specific BLAST results}" >> ${name}_${counter}.tex
        echo "\centering" >> ${name}_${counter}.tex
        echo "\resizebox{\textwidth}{!}{\begin{tabular}{ | l | l | l | l | l | l | l | l | l | l | l | }" >> ${name}_${counter}.tex
        echo "\hline" >> ${name}_${counter}.tex
        awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' table_formated_avium${counter} | tr "\t" "&" | sed 's/&/ & /g' | sed 's:$: \\\\:' | sed 's/_/\\_/g' | sed 's/name \\\\$/name \\\\ \\hline/' >> ${name}_${counter}.tex    
        echo "\hline" >> ${name}_${counter}.tex
        echo "\end{tabular}}" >> ${name}_${counter}.tex
        echo "\end{table}" >> ${name}_${counter}.tex    
    fi
    echo ""  >> ${name}_${counter}.tex
    echo "\vspace{5mm}" >> ${name}_${counter}.tex
    echo "" >> ${name}_${counter}.tex
    echo "\normalsize{Link to column header description:}" >> ${name}_${counter}.tex 
    echo "" >> ${name}_${counter}.tex
    echo "\normalsize{http://www.metagenomics.wiki/tools/blast/blastn-output-format-6}" >> ${name}_${counter}.tex
    echo "" >> ${name}_${counter}.tex 
    echo "\begin{figure}[ht]" >> ${name}_${counter}.tex
    echo "\centering" >> ${name}_${counter}.tex
    echo "\caption*{Maximum Likelihood Tree, Top 50 hits}" >> ${name}_${counter}.tex
    echo "\includegraphics[scale=0.45]{tree${counter}.pdf}" >> ${name}_${counter}.tex
    echo "\end{figure}" >> ${name}_${counter}.tex

    echo "\end{document}" >> ${name}_${counter}.tex
    echo "" >> ${name}_${counter}.tex

    ####

    pdflatex -interaction=nonstopmode ${name}_${counter}.tex
    pdflatex -interaction=nonstopmode ${name}_${counter}.tex

    mv "${name}_${counter}.pdf" "${name}_${counter}_${2}-report.pdf"

    ((counter++))
done < rpoB_found

cd ${root}

printf "Date ran: `date`\nSample name: $name\nSample type: $file_type\n\n" | mutt -a $(find . -name "${name}*report.pdf") -s "$name $2 analysis" -- tod.p.stuber@usda.gov

printf "\nDone\n\n\n"
# tstuber 2016-05-10
