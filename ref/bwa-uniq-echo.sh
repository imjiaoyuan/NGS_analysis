while read sample
do
	echo "bwa mem -t 8  -M -k 16   /mnt/may1nov1/u5011/genome/pto_hap1/pto_A/bwa_index/pto_1_ref   ./${sample}.clean_read1.fq.gz    ${sample}.clean_read2.fq.gz   | samtools sort -@ 8 -o ./${sample}.bam - &&  samtools flagstat  ./${sample}.bam  >./${sample}.bam.flagstat  &&  samtools  view  -h  ./${sample}.bam | grep -v -e \"MW376760.1\" -e \"MZ675536.1\" | samtools view -o ./${sample}.uniq.bam" >>bwa-uniq.sh
done
