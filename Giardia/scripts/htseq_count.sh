WORK_DIR=$(pwd)
GTF_FILE="${WORK_DIR}/data/wb.gtf"
MAP_DIR="${WORK_DIR}/mapping"
BAM_FILES="${MAP_DIR}/*.sorted.bam"

for BAM_FILE in ${BAM_FILES}
do
	prefix=$(echo $BAM_FILE| cut -d'.' -f 1)
	prefix=${prefix##*/}
	htseq-count -f bam -s no ${BAM_FILE}  ${GTF_FILE} >| ${prefix}.counts
done


