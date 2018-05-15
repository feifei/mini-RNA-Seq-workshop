GTF_FILE="/data/feifei/Others/Jingyi/data/wb.gtf"
WORK_DIR="/data/feifei/Others/Jingyi/mapping/"
BAM_FILES="${WORK_DIR}*.sorted.bam"

for BAM_FILE in ${BAM_FILES}
do
	prefix=$(echo $BAM_FILE| cut -d'.' -f 1)
	prefix=${prefix##*/}
	htseq-count -f bam -s no ${BAM_FILE}  ${GTF_FILE} >| ${prefix}.counts
done


