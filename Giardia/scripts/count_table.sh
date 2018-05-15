WORK_DIR=$(pwd)

READ1_FILES="${WORK_DIR}/data/*R1*.gz"

printf "Identifier\tTotal Reads\tMappped Reads\tMapped Reads p\tMapped Reads in Features\tMapped Reads in Features p\n"
for READ1_FILE in ${READ1_FILES}
do
	prefix=$(echo $READ1_FILE| cut -d'_' -f 1)
	prefix=${prefix##*/}

	printf "%s\t" $prefix
	
	
	
	# print the # reads sequenced
	n_total_reads=$(($(zcat ${READ1_FILE} | wc -l) / 4))
	printf "%d\t" ${n_total_reads}
	
	# Print # reads mapped to data
	line=$(head -n 1 ${WORK_DIR}/mapping/${prefix}.flagstat)
	n_mapped_reads=$(sed 's/^[^0-9]*\([0-9]\+\).*$/\1/' <<< "$line")
	printf "%d\t" ${n_mapped_reads}
	printf "%0.2f%%\t" $(bc <<< "scale=2; 100*$n_mapped_reads/($n_total_reads*2)")

	# Print # reads mapped to features
	count_file="${WORK_DIR}/mapping/$prefix.counts"	
	n_mapped_reads_feature=$(head -n-5 $count_file| cut -f 2 | awk '{s+=$0} END {printf "%d", s}')
	printf "%d\t" ${n_mapped_reads_feature}
	printf "%0.2f%%\t" $(bc <<< "scale=2; 100*$n_mapped_reads_feature/($n_total_reads*2)")

	printf "\n"
done