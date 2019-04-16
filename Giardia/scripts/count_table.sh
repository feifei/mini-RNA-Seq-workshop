#!/bin/bash

WORK_DIR=$(dirname $PWD)
MAP_DIR="${WORK_DIR}/mapping2"
COUNT_FILES="${MAP_DIR}/*_ReadsPerGene.out.tab"


printf "Identifier\tSequenced_read_pairs\tMapped_read_pairs\tMapped_read_p\tMapped_reads_in_features\tMapped_reads_in_features_p\n"


for COUNT_FILE in $COUNT_FILES 
do
	prefix=$(echo $COUNT_FILE| cut -d'_' -f 1)
	prefix=${prefix##*/}
	
	# print identifier
	printf "${prefix}\t"

	# print the # reads sequenced
	n_total_reads=$(cat ${prefix}_Log.final.out |grep "Number of input reads" |cut -d'|' -f2)
	printf "%d\t" ${n_total_reads}
	
	# Print # reads mapped to data	
	n_mapped_reads=$(cat ${prefix}_Log.final.out |grep "Uniquely mapped reads number" |cut -d'|' -f2)
	printf "%d\t" ${n_mapped_reads}
	printf "%0.2f%%\t" $(bc <<< "scale=2; 100*$n_mapped_reads/($n_total_reads)")
	
	
	# Count from ReadsPerGene.out.tab
	awk 'NR>4 {print $1 "\t" $2}' ${prefix}_ReadsPerGene.out.tab >| ${prefix}.star.counts
	
	# Print # reads mapped to features
	n_mapped_reads_feature=$(cat ${prefix}.star.counts| cut -f 2 | awk '{s+=$0} END {printf "%d", s}')
	printf "%d\t" ${n_mapped_reads_feature}
	printf "%0.2f%%\n" $(bc <<< "scale=2; 100*$n_mapped_reads_feature/($n_total_reads*2)")


done

