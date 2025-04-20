#!/bin/bash
#SBATCH --job-name=RNAseq_fastp
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/rna-seq

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/rna-seq/logs"
RAW_DIR="/mnt/may1nov1/u5023/pal/rna-seq/rawdata"
CLEAN_DIR="/mnt/may1nov1/u5023/pal/rna-seq/cleandata"
REPORT_DIR="/mnt/may1nov1/u5023/pal/rna-seq/report"
TMP_DIR="/mnt/may1nov1/u5023/pal/rna-seq/tmp"

mkdir -p ${LOG_DIR} ${CLEAN_DIR} ${REPORT_DIR} ${TMP_DIR}

eval "$(${CONDA_EXE} shell.bash hook)" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1
conda activate ${CONDA_ENV} >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "错误:conda环境激活失败" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log
    exit 1
fi

FASTQ_FILES_R1=($(ls ${RAW_DIR}/*_R1.fq.gz 2>> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log | sort -V))
FASTQ_FILES_R2=($(ls ${RAW_DIR}/*_R2.fq.gz 2>> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log | sort -V))

if [ ${#FASTQ_FILES_R1[@]} -eq 0 ]; then
    echo "错误:未找到R1端fastq文件" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log
    exit 1
elif [ ${#FASTQ_FILES_R2[@]} -eq 0 ]; then
    echo "错误:未找到R2端fastq文件" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log
    exit 1
elif [ ${#FASTQ_FILES_R1[@]} -ne ${#FASTQ_FILES_R2[@]} ]; then
    echo "错误:R1/R2文件数量不匹配" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log
    exit 1
fi

echo "找到 ${#FASTQ_FILES_R1[@]} 对双端fastq文件" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log

for ((i=0; i<${#FASTQ_FILES_R1[@]}; i++)); do
    R1=${FASTQ_FILES_R1[$i]}
    R2=${FASTQ_FILES_R2[$i]}
    SAMPLE_NAME=$(basename ${R1} _R1.fq.gz)
    
    echo "正在处理样本: ${SAMPLE_NAME}" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log
    
    fastp -i ${R1} -I ${R2} \
          -o ${CLEAN_DIR}/${SAMPLE_NAME}_R1.fq.gz \
          -O ${CLEAN_DIR}/${SAMPLE_NAME}_R2.fq.gz \
          --detect_adapter_for_pe \
          -w 4 \
          --compression 9 \
          -h ${TMP_DIR}/${SAMPLE_NAME}.html \
          -j ${TMP_DIR}/${SAMPLE_NAME}.json >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1
    
    fastqc -o ${TMP_DIR} -t 4 ${CLEAN_DIR}/${SAMPLE_NAME}_R1.fq.gz >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1
done

multiqc ${TMP_DIR} -o ${REPORT_DIR} -n "RNAseq_post_QC_report" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1

mkdir -p ${REPORT_DIR}/fastp_reports
mv ${TMP_DIR}/*.html ${REPORT_DIR}/fastp_reports/ >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1
mv ${TMP_DIR}/*.json ${REPORT_DIR}/fastp_reports/ >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1

find ${TMP_DIR} -maxdepth 1 -type f -name "*.html" -o -name "*.json" -o -name "*_fastqc*" -delete >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log 2>&1

echo "处理完成" >> ${LOG_DIR}/RNAseq_fastp_${SLURM_JOB_ID}.log