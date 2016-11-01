#!/usr/bin/env perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;

my $file = $ARGV[0];
if (not defined $file) {
    die "### Need FASTQ as arg1\n"
}

(my $file_name = $file) =~ s/\..*//;

if ( -B $file ) {
    `pigz -d $file`;
    $file =~ s/.gz//;
}

my $frag_size_total=0;
my $counter=0;
my $ave_length;
my $inblast = "$file_name" . "_to_blast.fasta";

my $inseq = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $outseq = Bio::SeqIO->new(-file => ">$inblast", -fomat => 'fasta');
while (my $seq_obj = $inseq->next_seq) {
    $counter++;
    my $length = $seq_obj->length;
    $frag_size_total = $frag_size_total + $length;
    if ( $length > 2000 ){
        my $subsequence=$seq_obj->trunc(50,1950);
        $outseq->write_seq($subsequence);
    }else{
        $outseq->write_seq($seq_obj);
    }
}

print "\n$counter contigs\n";
$frag_size_total =~ s/(?<=\d)(?=(?:\d\d\d)+\b)/,/g; 
print "Total bases: $frag_size_total\n\n";

my $blast_order = "qlen slen pident mismatch evalue bitscore stitle saccver qseqid";
my $outblast = "$file_name" . "_id.txt";
my $outblast_uniq = "$file_name" . "_id_uniq.txt";
# word_size needs to be set to 11 to return hits for all contigs
`blastn -query $inblast -db /data/BLAST/db/nt -word_size 11 -num_threads 40 -out $outblast -max_target_seqs 1 -outfmt "6 $blast_order"`;

# BLAST full contig length
#`blastn -query $file -db /data/BLAST/db/nt -word_size 11 -num_threads 40 -out $outblast -max_target_seqs 1 -outfmt "6 $blast_order"`;

# know bug in -max_target_seqs.  Parameter used within blastn algorithm causing multiple hits to sometimes be returned
`awk -F'\t' '!seen[\$9]++' $outblast > $outblast_uniq`;

# place blast out file into array.  each line and element in list
open ( my $blast_handle, '<', $outblast_uniq) or die "Can't open $outblast";
chomp(my @lines = <$blast_handle>);
close $blast_handle;

# put accession numbers into list
my @accessions;
foreach my $line (@lines) {
    my @tab = split(/\t/, $line);
    push @accessions, $tab[-2];
}

# count the number each accession occures, put into hash
my %id_count;
my %acc_counts;
my $lenght_of_contig;
$acc_counts{$_}++ for @accessions;
# for each key in hash
while (my ($key, $value) = each %acc_counts) {
    my $isolate_name;
    # iterate through lines to get isolate name matching accession
    my $count_total_length=0;
    foreach my $line (@lines) {
        my @tab = split(/\t/, $line);    
        if ($key eq $tab[-2]){
            # right to variable, when for multi times will just overwrite
            $isolate_name = $tab[-3];
            # get the total bases identified with uniq accession
            # only works with spades output with length in header
            ($lenght_of_contig = $tab[-1]) =~ s/.*_length_(.*)_cov.*/$1/;
            print "lenght_of_contig: $lenght_of_contig\n";
            $count_total_length = $count_total_length + $lenght_of_contig;
            print "count_total_length: $count_total_length\n";
        }
    }
    # value is the number accession was counted in BLAST output
    # put isolate id and count value into hash for sorting below
    $count_total_length =~ s/(?<=\d)(?=(?:\d\d\d)+\b)/,/g;
    $id_count{"$isolate_name" . "_" . "$count_total_length"} = $value;
}

#while (my ($key, $value) = each %acc_counts) {
#    my $count_total_length=0;
#    my $lenght_of_contig;
#    # iterate through lines to get isolate name matching accession
#    foreach my $line (@lines) {
#        my @tab = split(/\t/, $line);
#        if ($key eq $tab[-2]){
#            # right to variable, when for multi times will just overwrite
#            ($lenght_of_contig = $tab[-1]) =~ s/.*_length_(.*)_cov.*/$1/;
#            print "lenght_of_contig: $lenght_of_contig\n";
#            $count_total_length = $count_total_length + $lenght_of_contig; 
#            print "count_total_length: $count_total_length\n";
#        }
#    }
#    # value is the number accession was counted in BLAST output
#    # put isolate id and count value into hash for sorting below
#}

# output sort
# sorted on count number with id information
# reverse with: {$id_count{$b} <=> $id_count{$a}}
foreach my $name (sort {$id_count{$a} <=> $id_count{$b}} keys %id_count) {
    print "$id_count{$name} --> $name\n";
}

# check that blast out elements == in reads, ie number of expected ids output
my $array_size = scalar @lines;
print "The number of top hits: $counter\n";
if ( $counter == $array_size ){
        print "\nPASS -> Records in, match records out\n\n";
    }else{
        print "\n### WARNING read input: $counter uniq BLAST output: $array_size\n\n"
    }

print "\n$counter contigs\n";
print "Total bases: $frag_size_total\n\n";

print "\nDONE\n\n";
# 2016-10-27 stuber
