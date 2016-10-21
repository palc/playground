#!/usr/bin/env perl

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;

my $infasta = $ARGV[0];
if (not defined $infasta) {
    die "ERROR!!! --> An input file was not provided.  Provide FASTA!"
}
print "# --> BEGIN\n\n";
print "Working on: $infasta\n";
(my $file_name = $infasta) =~ s/.fasta//;

open (my $log, '>', "log_" . $file_name . ".txt") or die "$!";

my $outfile_name = $file_name . "-hformat.fasta";
open(my $out_handle, ">", $outfile_name)
    or die "Could not open file $outfile_name";

my $infastafirst_obj = Bio::SeqIO->new(-file => $infasta, -format => "fasta");
while (my $fasta_obj = $infastafirst_obj->next_seq){
    my $header = $fasta_obj->desc;
    $header =~/(\(.*\)\))/;
    my $main_desc = $1;
    $main_desc =~ s/ /_/g;
    $main_desc =~ s/\(/_/;
    $main_desc =~ s/\(/_/;
    $main_desc =~ s/\)/_/g;
    $main_desc =~ s/_-/-/g;
    $main_desc =~ s:/:-:g;
    $main_desc =~ s/__/_/;
    $main_desc =~ s/__/_/;
    $main_desc =~ s/^_//;
    $main_desc =~ s/_$//;
    print $out_handle ">$main_desc\n";
    my $sequence = $fasta_obj->seq;
    print $out_handle "$sequence\n";
}
close $out_handle;

my $total_count=0;
my $total_length=0;
print "outfile_name: $outfile_name\n";
my $infasta_obj = Bio::SeqIO->new(-file => $outfile_name, -format => "fasta");
while (my $fasta_obj = $infasta_obj->next_seq){
    $total_count++;
    my $length = $fasta_obj->length;
#    print "length: $length\n";
    $total_length = $total_length + $length;
#    print "total_length: $total_length\n";
}

print "total_count: $total_count\n";
print $log "Total FASTA sequences in file: $total_count\n";
my $average_length = $total_length/$total_count;
print "\nAverage FASTA length: $average_length\n\n";
print $log "Average FASTA length: $average_length\n\n";
my $minus_number=90;
my $select_size = $average_length - $minus_number;
#print "Removing FASTAs < $select_size\n";
print "If sequences $minus_number less than the average of $average_length they are removed and listed below...\n";
print $log "If sequences $minus_number less than the average of $average_length they are removed and listed below...\n";

$infasta_obj = Bio::SeqIO->new(-file => $outfile_name, -format => "fasta");
my $select_file = "$file_name" . "_selected.fasta";
my $seqout_obj = Bio::SeqIO->new(-file => ">$select_file", -fomat => 'fasta');
while ( my $seq_obj = $infasta_obj->next_seq ) {
    # print the sequence
    my $scaffold_length = $seq_obj->length;
    my $scaffold_display_display_id = $seq_obj->display_id;
#    print "scaffold_length: $scaffold_length\n";
#    print "display_id: $scaffold_display_display_id\n";
    if ( $scaffold_length > $select_size ) {
        $seqout_obj->write_seq($seq_obj);
    }else{
        print "### CAUTION $scaffold_length bases < threshold of $select_size bases. REMOVED --> $scaffold_display_display_id\n\n";
        print $log "### CAUTION $scaffold_length bases < threshold of $select_size bases. REMOVED --> $scaffold_display_display_id\n\n";
    }
}

$seqout_obj->close();

###
print "T-Coffee running...\n\n";
`t_coffee -in=$select_file -mode=regular -output=fasta_aln -case=upper -run_name=$select_file -clean_seq_name=1 &> /dev/null`;
###

$total_count=0;
my $add_all_lengths=0;
my $length;
my $check_position=5;
my $check_character="-";
my $in_aln="$select_file" . "_aln";

# get aligned fasta file lengths
my $in_fasta_obj = Bio::SeqIO->new(-file => $in_aln, -format => "fasta");
while (my $fasta_obj = $in_fasta_obj->next_seq){
    $total_count++;
    $length = $fasta_obj->length;
    $add_all_lengths = $length + $add_all_lengths;
}

# check that the average length is the same as the last length in alignment file
my $length_ave=$add_all_lengths/$total_count;
if ($length == $length_ave) {
    print "--> Can continue all incoming lengths the same\n";
    print $log "--> Can continue all incoming lengths the same\n";
} else {
    print $log "\n####FASTA sizes are not all the same therefore not in alignment\n\n###SCRIPT DIED AND EXITED";
    die "\n####FASTA sizes are not all the same therefore not in alignment\n\n"
}

print "\nTotal contig count: \t $total_count\n";
print $log "\nTotal contig count: \t $total_count\n";
print "Alignment file incoming FASTA sizes: $length_ave\n";
print $log "Alignment file incoming FASTA sizes: $length_ave\n";

# find the position to cut at 5' end
sub cut_forward_ends {
    my $in_fasta_obj = Bio::SeqIO->new(-file => $in_aln, -format => "fasta");
    while (my $fasta_obj = $in_fasta_obj->next_seq){
        #print "check_character: $check_character\n";
        #print "check_position: $check_position\n\n";
        $check_character=$fasta_obj->subseq($check_position,$check_position);
        if ($check_character eq "-"){
            return $check_position++;
        }
    }
}

my $forward_cut_position;

while ($check_character eq "-"){
    cut_forward_ends;
    $forward_cut_position=$check_position;
}

$check_position=$length - 5;

# find position to cut at 3' end
sub cut_reverse_ends {
    my $in_fasta_obj = Bio::SeqIO->new(-file => $in_aln, -format => "fasta");
    while (my $fasta_obj = $in_fasta_obj->next_seq){
        $check_character=$fasta_obj->subseq($check_position,$check_position);
        if ($check_character eq "-"){
            return $check_position--;
        }
    }
}

$check_character="-";
my $reverse_cut_position;

while ($check_character eq "-"){
    cut_reverse_ends;
    $reverse_cut_position=$check_position;
}

print "\n***forward_cut_position: $forward_cut_position\n";
print $log "\n***forward_cut_position: $forward_cut_position\n";
print "***reverse_cut_position: $reverse_cut_position\n\n";
print $log "***reverse_cut_position: $reverse_cut_position\n\n";

# write out
my $out_alignment_file = "$file_name" . "-trimmed" . ".fasta_aln";
my $outseq_obj = Bio::SeqIO->new(-file => ">$out_alignment_file", -format => 'fasta');

$in_fasta_obj = Bio::SeqIO->new(-file => $in_aln, -format => "fasta");
while (my $fasta_obj = $in_fasta_obj->next_seq){
    my $subsequence=$fasta_obj->trunc($forward_cut_position,$reverse_cut_position);
    $outseq_obj->write_seq($subsequence);
}

# reset variables
$total_count=0;
$length=0;
$add_all_lengths=0;

# check aligned fasta file lengths
$in_fasta_obj = Bio::SeqIO->new(-file => $out_alignment_file, -format => "fasta");
while (my $fasta_obj = $in_fasta_obj->next_seq){
    $total_count++;
    $length = $fasta_obj->length;
    $add_all_lengths = $length + $add_all_lengths;
}

# check that the average length is the same as the last length in alignment file
$length_ave=$add_all_lengths/$total_count;
if ($length == $length_ave) {
    print "--> Can continue all outgoing lengths the same\n";
    print $log "--> Can continue all outgoing lengths the same\n";
} else {
    print $log "\n####FASTA sizes are not all the same therefore not in alignment\n\nSCRIPT DIED AND EXITED-";
    die "\n####FASTA sizes are not all the same therefore not in alignment\n\n"
}

print "\nTotal contig count: \t $total_count\n";
print $log "\nTotal contig count: \t $total_count\n";
print "Trimmed FASTA sizes in alignment file: $length_ave\n\n";
print $log "Trimmed FASTA sizes in alignment file: $length_ave\n\n";

my $tree_file="$file_name" . ".tre";

print "tree file: $tree_file\n";
print "out_alignment_file: $out_alignment_file\n";

print "RAxML running...\n";
`raxmlHPC-SSE3 -f a -s $out_alignment_file -p 12345 -x 12345 -# 100 -m GTRCAT -n $tree_file`;

`rm *dnd *reduced RAxML_bootstrap* RAxML_info* RAxML_bipartitions*`;

print "\n### DONE\n\n";
print $log "\n### DONE\n";
# tstuber 2016-10-20
