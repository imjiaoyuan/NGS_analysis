while read sample
do
	echo "samtools view -b -q 30 -f 2 ./${sample}.uniq.bam >./${sample}.q30.bam" >>uniq-to-q30.sh
done
