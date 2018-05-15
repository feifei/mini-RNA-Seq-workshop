# declare -a tags=("star_eimeria" "star_gallus")
# declare -a gtfs=("../data/ToxoDB-29_EtenellaHoughton.gtf" "../data/Gallus_gallus.Gallus_gallus-5.0.86.gtf")
#
# for ((i=0;i<${#tags[@]};++i))
# do
# 	for j in `seq 1 10`
# 	do
# 		prefix=${tags[i]}_$j
# 		htseq-count -f bam -s no ${prefix}.sorted.bam ${gtfs[i]} >| ${prefix}.counts
# 	done
# done


for j in `seq 1 33`
do
	prefix="star_$j"
	htseq-count -f bam -s no ${prefix}.sorted.bam merged.gtf >| ${prefix}.counts
	awk '/ETH/{print}' ${prefix}.counts >${prefix}.eimeria.counts
	awk '/ENSGAL/{print}' ${prefix}.counts >${prefix}.gallus.counts	
done


