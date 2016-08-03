#!/bin/bash

# LICENSE AND DISCLAIMER
#
# Copyright (c) 2015 The Johns Hopkins University/Applied Physics Laboratory
#
# This software was developed at The Johns Hopkins University/Applied Physics Laboratory ("JHU/APL") that is the author thereof under the "work made for hire" provisions of the copyright law. Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation (the "Software"), to use the Software without restriction, including without limitation the rights to copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit others to do so, subject to the following conditions:
#
# 1. This LICENSE AND DISCLAIMER, including the copyright notice, shall be included in all copies of the Software, including copies of substantial portions of the Software;
#
# 2. JHU/APL assumes no obligation to provide support of any kind with regard to the Software. This includes no obligation to provide assistance in using the Software nor to provide updated versions of the Software; and
#
# 3. THE SOFTWARE AND ITS DOCUMENTATION ARE PROVIDED AS IS AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES WHATSOEVER. ALL WARRANTIES INCLUDING, BUT NOT LIMITED TO, PERFORMANCE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NONINFRINGEMENT ARE HEREBY DISCLAIMED. USERS ASSUME THE ENTIRE RISK AND LIABILITY OF USING THE SOFTWARE. USERS ARE ADVISED TO TEST THE SOFTWARE THOROUGHLY BEFORE RELYING ON IT. IN NO EVENT SHALL THE JOHNS HOPKINS UNIVERSITY BE LIABLE FOR ANY DAMAGES WHATSOEVER, INCLUDING, WITHOUT LIMITATION, ANY LOST PROFITS, LOST SAVINGS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF THE USE OR INABILITY TO USE THE SOFTWARE.”

#---------------------------------------------------------------------------------------------------
# This script converts Kraken's tabular output to a standard taxonomic format usable by tools such as Krona
#      by: Thomas Mehoke
# written: November 7, 2014
#---------------------------------------------------------------------------------------------------

export PATH=/data/apps/bin/surpi/:/data/apps/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

runtime=$(date +"%Y%m%d%H%M%S%N")

usage()
{
cat << EOF
usage: $0 -i <input.kraken> -k </path/to/kraken-db> -o <output.txt> -r <output.html>


OPTIONS:
   -h      show this message
   -i      output file from Kraken
   -k      Kraken database folder (path)
   -u      ignore unclassified results (taxid = 0)
   -o      file to place text output file
   -r      file to place Krona HTML file
   -p      path to Krona install directory (default: /data/apps/src/KronaTools-2.4)
   -d      threads to use for parallel processing
            The default is one less than the available processors.
            On this computer ($(hostname)) there are $(nproc) threads,
             so the default is $(($(nproc) - 1)).
   -w      working directory (default: /tmp)
EOF
}

# set default values
ignore="false"
workdir="/tmp"
THREADS=$(($(nproc) - 1))
kronafile=""
#stuberKRONA_PATH="/data/apps/src/KronaTools-2.4"
KRONA_PATH="/usr/local/bin/KronaTools-2.5"
outputfile=""
localkrona="true"

# parse input arguments
while getopts "hi:k:uw:d:o:r:p:" OPTION
do
	case $OPTION in
		h) usage; exit 1 ;;
		i) input=$OPTARG ;;
		k) DB=$OPTARG ;;
		u) ignore="true" ;;
		w) workdir=$OPTARG ;;
		d) THREADS=$OPTARG ;;
		o) outputfile=$OPTARG ;;
		r) kronafile=$OPTARG ;;
		p) KRONA_PATH=$OPTARG ;;
		?) usage; exit ;;
	esac
done

# if necessary arguments are not present, display usage info and exit
if [[ -z "$input" ]]; then
	echo "Specify a Kraken output file with -i" >&2
	usage
	exit 2
fi
if [[ -z "$DB" ]]; then
	echo "Select a kraken database with -k" >&2
	usage
	exit 2
fi
kronascript="$KRONA_PATH/src/krona-2.0.js"
if ! [[ -s "$kronascript" ]]; then
	echo "Warning: $KRONA_PATH/src/krona-2.0.js does not exist"
	echo "Will attempt to fetch from Krona server"
fi

# Need to add format checking to make sure kraken output is properly formatted
# Need to add check that Kraken DB contains names.dmp and nodes.dmp

#===================================================================================================
# function declarations

make_tax_string() {

	# parse input arguments
	count="$1"
	tax="$2"

	# pull time so temp files don't overwrite each other
	runtime_sub=$(date +"%Y%m%d%H%M%S%N")

	# pull out name and parent taxid for this taxon
	name=$(egrep "^$tax"$'\t' "$DB/taxonomy/names.dmp" | grep "scientific name" | cut -f3)
	parent=$(egrep "^$tax"$'\t' "$DB/taxonomy/nodes.dmp" | cut -f3)

	# deal with taxa near the root
	if [[ $tax -eq 0 ]]; then
		echo -e "$count\tUnclassified" >> "$OUTPUT"
	elif [[ $parent -eq 1 ]]; then
		if [[ "$ignore" == "true" ]]; then
			echo -e "$count\t$name" >> "$OUTPUT"
		elif [[ "$name" != "root" ]]; then
			echo -e "$count\troot\t$name" >> "$OUTPUT"
		else
			echo -e "$count\troot" >> "$OUTPUT"
		fi
	# for all others, create full taxonomic string
	else
		while [[ $parent -gt 1 ]]; do
			parentname=$(egrep "^$parent"$'\t' "$DB/taxonomy/names.dmp" | grep "scientific name" | cut -f3)
			name="$parentname\t$name"
			echo -e "$name" > "$tempdir/krona.temp-$runtime_sub"
			parent=$(egrep "^$parent"$'\t' "$DB/taxonomy/nodes.dmp" | cut -f3)
		done
		if [[ "$ignore" == "true" ]]; then
			echo -e "$count\t$(cat "$tempdir/krona.temp-$runtime_sub")" >> "$OUTPUT"
		else
			echo -e "$count\troot\t$(cat "$tempdir/krona.temp-$runtime_sub")" >> "$OUTPUT"
		fi
	fi
}

#===================================================================================================

tempdir="$workdir/kraken2krona-$runtime"
mkdir -p "$tempdir"

# download kronascript if it is not provided
if ! [[ -s "$kronascript" ]]; then
	kronascript="$tempdir/krona-2.0.js"
	wget -q http://krona.sourceforge.net/src/krona-2.0.js -O "$kronascript" || localkrona="false"
fi

tax_counts="$tempdir/tax_counts"
OUTPUT="$tempdir/output"

# create space-separated count of taxonomic ranks
if [[ "$ignore" == "true" ]]; then
	cut -f3 "$input" | sort -T "$tempdir" | uniq -c | sed 's/^ *//' | sort -T "$tempdir" -k1nr,1 | gawk -F" " '{if($2 > 1){print $0}}' > "$tax_counts"
else
	cut -f3 "$input" | sort -T "$tempdir" | uniq -c | sed 's/^ *//' | sort -T "$tempdir" -k1nr,1 > "$tax_counts"
fi

# loop through all taxonomic levels
export DB
export tempdir
export ignore
export OUTPUT
export -f make_tax_string
parallel -a "$tax_counts" -d $'\n' --colsep ' ' -n 1 -P "$THREADS" make_tax_string {1} {2}

# generate krona plot if specified by -k flag
if [[ -n "$kronafile" ]]; then
	export TERM=xterm-256color
	ktImportText -o "$kronafile" -n $(basename "$DB") "$OUTPUT"

	if [[ "$localkrona" == "true" ]]; then
		# add javascript directly to the HTML file so it can be read anywhere without an internet connection
		sed -i '/<script id="notfound">.*/d' "$kronafile"
		sed -i 's@<script src="http://krona.sourceforge.net/src/krona-2.0.js"></script>@<script>\n  </script>@' "$kronafile"
		sed -i "/<script>/ r $kronascript" "$kronafile"

		# make file viewable without an internet connection (I'm not sure why this line causes errors, but commenting it seems harmless)
		sed -i "s/^\(\s*hiddenPattern = context.createPattern(image, 'repeat');\)/\/\/\1/" "$kronafile"

		# unselect the collapse box
		sed -i "s/\(collapse = kronaElement.getAttribute('collapse') == '\)true';/\1false';/" "$kronafile"
	fi
fi

# output to STDOUT
if [[ -z "$outputfile" ]]; then
	sort -T "$tempdir" -k1nr,1 "$OUTPUT"
# or output to output file specified by -o argument
else
	sort -T "$tempdir" -k1nr,1 "$OUTPUT" > "$outputfile"
fi

rm -rf "$tempdir"

#~~eof~~#
