while read sample
do
	echo " bamCoverage  --bam  ./${sample}.last.bam  -o  ./${sample}.last-10.bw  --binSize 10 --effectiveGenomeSize  352383277  --normalizeUsing RPKM  --smoothLength 50  && bamCoverage --bam  ./${sample}.last.bam  -o  ./${sample}.last-50.bw  --binSize 50 --effectiveGenomeSize  352383277  --normalizeUsing RPKM  --smoothLength 150 && bamCoverage  --bam  ./${sample}.picard.bam  -o ./${sample}.picard-10.bw --binSize 10 --effectiveGenomeSize  352383277  --normalizeUsing RPKM  --smoothLength 50 && bamCoverage  --bam  ./${sample}.picard.bam  -o ./${sample}.picard-50.bw --binSize 50 --effectiveGenomeSize  352383277  --normalizeUsing RPKM  --smoothLength 150 " >> bam-to-bw.sh	
done
