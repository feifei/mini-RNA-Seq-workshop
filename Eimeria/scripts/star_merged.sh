WORK_DIR=$(pwd)
EIMERIA_GENOME="ToxoDB-29_EtenellaHoughton_Genome.fasta"
GALLUS_GENOME="Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa"

cat ${WORK_DIR}/${EIMERIA_GENOME} ${WORK_DIR}/${GALLUS_GENOME} >| ${WORK_DIR}/merged.fasta

STAR --runMode genomeGenerate --genomeDir ${WORK_DIR}/genomes --genomeFastaFiles ${WORK_DIR}/merged.fasta \
	--runThreadN 8 --sjdbGTFfile ${WORK_DIR}/merged.gtf --sjdbOverhang 125 --genomeChrBinNbits 15

for j in `seq 1 33`
do
    prefix=star_${j}
    READS=${WORK_DIR}/${j}_S*_L00*_R*_001.fastq.gz
    STAR --genomeDir ${WORK_DIR}/genomes --readFilesIn $READS --readFilesCommand gunzip -c \
			--runThreadN 16 --outFileNamePrefix ${WORK_DIR}/${prefix}_ --quantMode GeneCounts
    samtools view -u ${prefix}_Aligned.out.sam | \
    samtools sort -T tmp -o $prefix.sorted.bam -
    samtools index $prefix.sorted.bam
    samtools flagstat $prefix.sorted.bam > $prefix.flagstat
done
