WORK_DIR="/data/feifei/Others/Jingyi/"
GENOME="${WORK_DIR}data/GiardiaDB-5.0_GintestinalisAssemblageA_Genome.fasta"
GTF_FILE="${WORK_DIR}data/wb.gtf"
MAPPING_DIR="${WORK_DIR}mapping/"

echo $GENOME
STAR --runMode genomeGenerate --genomeDir ${MAPPING_DIR}genome --genomeFastaFiles $GENOME --sjdbGTFfile ${GTF_FILE} --runThreadN 16 --sjdbOverhang 125


READ1_FILES="${WORK_DIR}data/*R1*.gz"
for READ1 in ${READ1_FILES} 
	do
		prefix=$(echo $READ1| cut -d'_' -f 1)
		prefix=${prefix##*/}
		READ2=$(echo ${READ1/R1/R2})
		echo ${prefix}
		echo $READ2
		STAR --genomeDir ${MAPPING_DIR}genome --readFilesIn $READ1 $READ2 --readFilesCommand gunzip -c --runThreadN 16 --outFileNamePrefix ${MAPPING_DIR}${prefix}_ --quantMode GeneCounts
		samtools view -u ${prefix}_Aligned.out.sam | \
		samtools sort -T tmp -o $prefix.sorted.bam -
		samtools index $prefix.sorted.bam
		samtools flagstat $prefix.sorted.bam > $prefix.flagstat
	done
