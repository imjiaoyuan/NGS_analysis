while read sample
do
	echo "macs2 callpeak  -g 352383277  -t  ${sample}.bam  -f  BAMPE --nomodel --extsize 200 --shift -100 --pvalue 0.01   -B  -n  ${sample}" >> macs2-atac-peak-0.01.sh
done
