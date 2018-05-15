WORK_DIR=/data/feifei/Eimeria/mapping2/
cd ${WORK_DIR}

declare -a arr=("eimeria" "gallus")
declare -a hours=("0" "4" "4" "24" "24" "48" "48" "72" "72" "0" "48" "48" "72" "72" "48" "48" "72" "72" "2" "2" "2" "2" "12" "12" "12" "12" "12" "12" "0")
declare -a conditions=("control" "infected" "control" "infected" "control" "infected" "control" "infected" "control" "control" \
                       "infected" "control" "infected" "control" "infected" "control" "infected" "control" "infected" "infected" \
                       "control" "control" "infected" "infected" "infected" "control" "control" "control" "infected")
x=0
for j in `seq 1 33`
do
    if [ $j != 10 -a $j != 16 -a $j != 17 -a $j != 22 ] ; then
        # (10, 16, 17, 22 are empty)
        x=$(expr $x + 1)
    	# print identifier & hour & condition
    	printf "${j}\t${hours[x-1]}\t${conditions[x-1]}\t"
    	prefix=star_${j}
	
    	# print the # reads sequenced
    	read1=${j}_S*_L00*_R1_001.fastq.gz
    	n_total_reads=$(($(zcat ${read1} | wc -l) / 4))
    	printf "%d\t" ${n_total_reads}
	
    	# Print # reads mapped to data
    	line=$(head -n 1 ${prefix}.flagstat)
    	n_mapped_reads=$(sed 's/^[^0-9]*\([0-9]\+\).*$/\1/' <<< "$line")
    	printf "%d\t" ${n_mapped_reads}
	
    	printf "%0.2f%%\t" $(bc <<< "scale=2; 100*$n_mapped_reads/($n_total_reads*2)")
    	# Print # reads mapped to features
    	n_mapped_reads_feature=$(head -n-5 ${prefix}.counts| cut -f 2 | awk '{s+=$0} END {printf "%d", s}')
    	printf "%d\t" ${n_mapped_reads_feature}
    	printf "%0.2f%%\t" $(bc <<< "scale=2; 100*$n_mapped_reads_feature/($n_total_reads*2)")

    	for ((i=0;i<${#arr[@]};++i))
    	do
    		n_mapped_reads_part=$(head -n-5 ${prefix}.${arr[i]}.counts| cut -f 2 | awk '{s+=$0} END {printf "%d", s}')
    		printf "%d\t" ${n_mapped_reads_part} 
    		printf "%9.2f%%\t" $(bc <<< "scale=2; 100*$n_mapped_reads_part/$n_mapped_reads_feature")
    	done
    	printf "\n"
    fi
done