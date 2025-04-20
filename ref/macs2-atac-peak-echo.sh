while read sample
do
	echo "macs2 callpeak  -g 352383277  -t  ${sample}.bam  -f  BAMPE --nomodel --extsize 200 --shift -100 --pvalue 0.05   -B  -n  ${sample}" >> macs2-atac-peak.sh
done
