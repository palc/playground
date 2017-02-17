#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;

# Usage: $$ tree_builder.pl <in.fasta>
# in.fasta contains list of genomes to aligned and build tree
my $infasta = $ARGV[0];
if (not defined $infasta) {
    die "ERROR!!! --> An input file was not provided.  Provide FASTA!"
}
(my $file_name = $infasta) =~ s/.fasta//;

print "# --> BEGIN\n\n";

# Output log file
open (my $log, '>', "log_" . $file_name . ".txt") or die "$!";

# Create hash
# desc and display_id placed into key
# sequence placed into value
my %sequences;
my $seqio = Bio::SeqIO->new(-file => $infasta, -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
    my $main_desc;
    my $desc  = $seqobj->desc;
    my $header_name;
    print "desc: $desc\n";
    my $display_id  = $seqobj->display_id;
    print "display_id: $display_id\n";
    my $id  = $display_id . " " . $desc;    # there's your key
    $display_id =~ s/_.*//;
    if ($id =~ /\(/ ){
        #print "there is a parethesis in header\n";
        $id =~/(\(.*\)\))/;
        $main_desc = $1;
        $main_desc =~ s/ /_/g;
        $main_desc =~ s/\(/_/;
        $main_desc =~ s/\(/_/;
        $main_desc =~ s/\)/_/g;
        $main_desc =~ s/_-/-/g;
        #$main_desc =~ s:/:-:g;
        $main_desc =~ s/__/_/;
        $main_desc =~ s/__/_/;
        $main_desc =~ s/^_//;
        $main_desc =~ s/_$//;
        #print "$main_desc\n";    
}else{
        #print "there is NOT a parethesis\n";
        $main_desc = $id;
        $main_desc =~ s/ /_/g;
        $main_desc =~ s/__/_/;
        $main_desc =~ s/__/_/;
        $main_desc =~ s/^_//;
        $main_desc =~ s/_$//;
        #print "$main_desc\n";
    }
    my $seq = $seqobj->seq;                 # and there's your value
    $header_name = $display_id; # . "_" . $main_desc;
    $sequences{$header_name} = $seq;
}

# Get average FASTA length
my $total_count=0;
my $total_length=0;

my $size = keys %sequences;
my $total_seq = 0;
print "Number of sequences: $size\n";
while (my ($key, $value) = each %sequences) {
    $total_seq = $total_seq + length $value;
}

print "total_count: $size\n";
print $log "Total FASTA sequences in file: $size\n";

print "$total_seq\n";
my $average_size = $total_seq / $size;
print "average size: $average_size\n";
print $log "Average FASTA length: $average_size\n\n";

# Remove genomes < $select_size
my $minus_number=9000;
my $select_size = $average_size - $minus_number;

print "If sequences $minus_number less than the average of $average_size they are removed and listed below...\n";
print $log "If sequences $minus_number less than the average of $average_size they are removed and listed below...\n";

my $select_file = $file_name . "_selected.fasta";
open (my $out_file, '>', $select_file) or die "$!";

while (my ($key, $value) = each %sequences) {
    my $value_length = length $value;
    if ( $value_length > $select_size ) {
        print $out_file ">$key\n$value\n";
    }else{
        print "### CAUTION $value_length bases < threshold of $select_size bases. REMOVED --> $key\n\n";
        print $log "### CAUTION $value_length bases < threshold of $select_size bases. REMOVED --> $key\n\n";
    }
}

$out_file->close();

###
# T-Coffee to align sequences
print "T-Coffee running...\n\n";
`t_coffee -in=$select_file -mode=regular -output=fasta_aln -case=upper -run_name=$select_file -clean_seq_name=1 &> /dev/null`;
###

# Remove "-" ends...
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

# gapped "-" ends removed
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

# run aligned, trimmed "-" gapped ends FASTAs in RAxML
# -f a --> conduct a rapid Bootstrap analysis and search for the best-scoring ML tree in one single program run
# -s --> input file
# -p --> random number seed for the parsimony inferences, allows one to reproduce your results
# -x --> integer number (random seed) and turn on rapid bootstrapping, used in place of -b
# -m --> model, GTR: Generalised time-reversible
# -n --> out file
# bootstrap with --> -b 123476 -N 100 (using -f a option... much faster)
#   -b --> random seed
#   -N --> number of alternative runs on distinct starting trees
print "RAxML running...\n";
`raxmlHPC-PTHREADS-AVX2 -s $out_alignment_file -f a -x 12345 -T 50 -p 12345 -N 100 -m GTRCAT -n $tree_file`;

# remove unneeded files
`rm *dnd *reduced RAxML_bootstrap* RAxML_info* RAxML_bipartitions*`;

my $svgname = $file_name . ".svg";
my $pdfname = $file_name . ".pdf";
my $raxmlname = "RAxML_bestTree." . $tree_file;
print "raxmlname: $raxmlname\n";

# create svg and pdf trees
`cat $raxmlname | nw_display -s -S -w 1300 -t -v 30 -i 'opacity:0' -b 'opacity:0' -l 'font-size:18;font-family:serif;font-style:italic' -d 'stroke-width:1;stroke:blue' - > $svgname && inkscape -f $svgname -A $pdfname`;
             
print "\n### DONE\n\n";
print $log "\n### DONE\n";
# tstuber 2016-10-20
