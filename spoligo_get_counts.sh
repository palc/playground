#!/bin/sh

infile = $1

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

fprint "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$infile" "f19" "r19" "f34" "r34" "f36" "r36" "f19a1" "r19a1" "f34a1" "r34a1" "f36a1" "r36a1"

# tstuber 2017-04-19
