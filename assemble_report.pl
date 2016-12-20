#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;
use Data::Dumper qw(Dumper);
use Cwd;

# Usage:
    # virus_segment_finder.pl < full path to fastq files >
    # Files can be zipped or unzipped, single or paired reads
    # Paired files must be labed R1 and R2

print "----- START -----\n\n";



# one argument expected , which must be path to fastq files 
#my $path = $ARGV[0];
my $path = getcwd;
if (not defined $path) {
    die "ERROR!!! --> Argument to directory containing FASTQ/s was not provided"
}

# place a "/" at the end of path if one not provided
$path =~ s|/?$|/|;
#print "Path to files: $path\n";

# Determine if single or paired read
my @files = < ${path}*fastq* >;
my $read_type;
my $filename;
# check if array is empty
if (@files) {
    foreach (@files) {
        # get basename
        my $basename = basename($_);
        $filename=$basename;
        chomp $filename;
        # try to find R[1-2] from basename
        if ($basename =~ /R[1-2]{1}/g) {
            $read_type = "paired";
        } else {
            $read_type = "single";
        }
    }
} else {
    die "ERROR!!! FASTQ not found at $path Look";
}

print "Read type found: $read_type\n\n";

# make zipped and unzipped reads available
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);
my @files = < ${path}*fastq* >;
foreach my $file (@files) {
    # -B  File is a "binary" file
    # -T  File is an ASCII or UTF-8 text file
    # -e File exists
    # Test if FASTQs are zipped or not.  Make zip and unzip available
    # Binary and text file not present
    # strip .gz
    my $no_gz = substr($file, 0, -3);
        if ( -B $file && ! ( -e $no_gz )) {
        (my $file_noext = $file) =~ s/\.[^.]+$//;
        gunzip $file => $file_noext or die "ERROR!!! gunzip failed for $file";
        }
    my $with_gz = $file . ".gz";
        # Text file and binary not present
    if ( -T $file && ! ( -e $with_gz )) {
        gzip $file => "${file}.gz" or die "gzip failed: $GzipError";
        }
}

# allow to be globally available
# paired reads
my $input_R1_zip;
my $input_R2_zip;
my $input_R1_unzip;
my $input_R2_unzip;
# single reads
my $input_zip;
my $input_unzip;

if ($read_type eq "paired") {
    # place zip and unzipped file into variables
    my @input_R1_zip = < ${path}*R1*fastq.gz >;
    my @input_R2_zip = < ${path}*R2*fastq.gz >;
    my @input_R1_unzip = < ${path}*R1*fastq >;
    my @input_R2_unzip = < ${path}*R2*fastq >;

    foreach (@input_R1_zip) {
        $input_R1_zip = $_;
        }
    foreach (@input_R2_zip) {
        $input_R2_zip = $_;
        }
    foreach (@input_R1_unzip) {
        $input_R1_unzip = $_;
        }
    foreach (@input_R2_unzip) {
        $input_R2_unzip = $_;
        }

    print "Zipped Files:\n";
    print "$input_R1_zip\n";
    print "$input_R2_zip\n\n";

    print "Unzipped Files:\n";
    print "$input_R1_unzip\n";
    print "$input_R2_unzip\n\n";
} else {
    
    my @input_zip = < ${path}*fastq.gz >;
    my @input_unzip = < ${path}*fastq >;
    
    foreach (@input_zip) {
        $input_zip = $_;
    }
    foreach (@input_unzip) {
        $input_unzip = $_;
    }
    
    print "Zipped File:\n";
    print "$input_zip\n";
    
    print "Unzipped File:\n";
    print "$input_unzip\n\n";
}

if ($read_type eq "paired") {
    print "RUNNING SPADES\n\n";
    `spades.py -t 32 -k 45,47,49,51,53,57,59,61,63,65,77,99,127 --careful -1 $input_R1_zip -2 $input_R2_zip -o ./`;
} else {
    print "RUNNING SPADES\n\n";
    `spades.py -t 32 -k 45,47,49,51,53,57,59,61,63,65,77,99,127 --mismatch-correction -s $input_zip -o ./`;
}






























# 2016-12-20 tstuber
