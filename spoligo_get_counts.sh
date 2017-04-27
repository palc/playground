#!/bin/sh

infile=$1

f19=`grep -c "AAGTCACCTCGCCCACACCGTCGAA" $infile`
r19=`grep -c "TTCGACGGTGTGGGCGAGGTGACTT" $infile`

f34=`grep -c "AAGTCACCTCGCCCACACCGTCGAA" $infile`
r34=`grep -c "TTCGACGGTGTGGGCGAGGTGACTT" $infile`

f36=`grep -c "CGAAATCCAGCACCACATCCGCAGC" $infile`
r36=`grep -c "GCTGCGGATGTGGTGCTGGATTTCG" $infile`

f19a1=`agrep -1 -c "AAGTCACCTCGCCCACACCGTCGAA" $infile`
r19a1=`agrep -1 -c "TTCGACGGTGTGGGCGAGGTGACTT" $infile`

f34a1=`agrep -1 -c "AAGTCACCTCGCCCACACCGTCGAA" $infile`
r34a1=`agrep -1 -c "TTCGACGGTGTGGGCGAGGTGACTT" $infile`

f36a1=`agrep -1 -c "CGAAATCCAGCACCACATCCGCAGC" $infile`
r36a1=`agrep -1 -c "GCTGCGGATGTGGTGCTGGATTTCG" $infile`

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$infile" "spacercount: $f19" "spacercount: $r19" "spacercount: $f34" "spacercount: $r34" "spacercount: $f36" "spacercount: $r36" "spacercount: $f19a1" "spacercount: $r19a1" "spacercount: $f34a1" "spacercount: $r34a1" "spacercount: $f36a1" "spacercount: $r36a1"


# tstuber 2017-04-19
