#!/bin/bash
#SBATCH --job-name=RNAseq_QC
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -D /mnt/may1nov1/u5023/pal/rna-seq

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"

LOG_DIR="/mnt/may1nov1/u5023/pal/rna-seq/logs"
RAW_DIR="/mnt/may1nov1/u5023/pal/rna-seq/cleandata"
REPORT_DIR="/mnt/may1nov1/u5023/pal/rna-seq/report"
TMP_DIR="/mnt/may1nov1/u5023/pal/rna-seq/tmp"

mkdir -p ${LOG_DIR} ${REPORT_DIR} ${TMP_DIR}

eval "$(${CONDA_EXE} shell.bash hook)" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log 2>&1
conda activate ${CONDA_ENV} >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log 2>&1

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "错误:conda环境激活失败" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log
    exit 1
fi

FASTQ_FILES=($(ls ${RAW_DIR}/*.fq.gz 2>> ${LOG_DIR}/job_${SLURM_JOB_ID}.log))

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "错误:未找到任何fastq.gz文件" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log
    exit 1
fi

echo "找到 ${#FASTQ_FILES[@]} 个fastq文件需要处理" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log

for fastq in "${FASTQ_FILES[@]}"; do
    echo "正在处理: $(basename ${fastq})" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log
    fastqc -o ${TMP_DIR} -t ${SLURM_CPUS_PER_TASK} ${fastq} >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log 2>&1
done

wait

echo "合并质控报告..." >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log
multiqc ${TMP_DIR} -o ${REPORT_DIR} -n "RNAseq_QC_report" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log 2>&1
rm -rf ${TMP_DIR}/* >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log 2>&1

echo "质控分析完成报告保存在: ${REPORT_DIR}" >> ${LOG_DIR}/job_${SLURM_JOB_ID}.log
