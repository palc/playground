#/bin/sh

QUAL=300 # Minimum quality for calling a SNP
NR_CPUS=40

egrep -vh "#" ./starting_files/*vcf | egrep "AC=2;A" | awk -v Q="$QUAL" ' $6 > Q {print $2}' | sort | uniq > prepositionlist

for n  in `cat prepositionlist`; do

(positioncount=`awk -v n=$n ' $2 == n {count++} END {print count}' ./starting_files/*vcf`
echo "position count: $positioncount"        
if [ $positioncount -gt 2 ]; then
        echo "$n" >> positionlist
fi) &
        let count+=1
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done

for p in `cat positionlist`; do

(maxqual=`awk -v p=$p 'BEGIN{max=0} $2 == p {if ($6>max) max=$6} END {print max}' ./starting_files/*vcf | sed 's/\..*//'`

        maxmap=`awk -v p=$p ' $2 == p {print $8}' ./starting_files/*vcf | sed 's/.*;MQ=\(.*\);MQ.*/\1/' | sed 's/;MQ.*//' | awk 'BEGIN{max=0}{if ($1>max) max=$1} END {print max}' | sed 's/\..*//'`

        if [ $maxqual -lt 800  ] || [ $maxmap -lt 52  ]; then
                echo "maxqual $maxqual" >> filterpositiondetail
                echo "maxmap $maxmap" >> filterpositiondetail
                echo "position $p" >> filterpositiondetail
                echo ""  >> filterpositiondetail
                echo "$p" >> ${d}-filterposition

                echo "maxqual $maxqual"
                echo "maxmap $maxmap"
                echo "position $p"
                echo ""
        fi) &
        let count+=1
        [[ $((count%NR_CPUS)) -eq 0 ]] && wait
done


# for i in *-filterposition; do name=`echo $i | sed 's/-filterposition//'`; echo $name; echo "$name" > temp; cat $i >> temp; mv temp ${i}.txt; done
# paste *txt >> file.txt
# Created 2015-04-15, Stuber
