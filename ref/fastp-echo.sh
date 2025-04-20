while read line
do
        arr=($line)
        fq1=${arr[0]}
        fq2=${arr[1]}
        sample=${arr[2]}
        echo " fastp --in1 $fq1  --in2 $fq2  --out1 ./${sample}.clean_read1.fq.gz --out2  ${sample}.clean_read2.fq.gz  --detect_adapter_for_pe   --html ./${sample}.html  --json ${sample}.json "  >> fastp.sh
done

