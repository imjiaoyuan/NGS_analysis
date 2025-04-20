while read sample
do
	echo "samtools rmdup ./${sample}.q30.bam  ./${sample}.last.bam && samtools index ./${sample}.last.bam && samtools flagstat ./${sample}.last.bam > ./${sample}.last.flagstat &&  picard MarkDuplicates  I=./${sample}.q30.bam   O=./${sample}.picard.bam  M=./${sample}.picard_marked_dup_metrics.txt REMOVE_DUPLICATES=true && samtools index ./${sample}.picard.bam && samtools flagstat ./${sample}.picard.bam > ./${sample}.picard.flagstat  " >> rmdup.sh
done
