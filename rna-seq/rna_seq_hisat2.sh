#!/bin/bash
#SBATCH --job-name=RNAseq_hisat2
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/rna-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/rna-seq/logs/RNAseq_hisat2_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/rna-seq/logs/RNAseq_hisat2_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/rna-seq/logs"
CLEAN_DIR="/mnt/may1nov1/u5023/pal/rna-seq/cleandata"
SAM_DIR="/mnt/may1nov1/u5023/pal/rna-seq/sam_files"
REPORT_DIR="/mnt/may1nov1/u5023/pal/rna-seq/report"
TMP_DIR="/mnt/may1nov1/u5023/pal/rna-seq/tmp"
GENOME_FASTA="/mnt/may1nov1/u5023/pal/Palv.hap2.chr.genome.fa"
INDEX_PREFIX="/mnt/may1nov1/u5023/pal/rna-seq/Palv.hap2.chr.genome"

mkdir -p ${LOG_DIR} ${CLEAN_DIR} ${SAM_DIR} ${REPORT_DIR} ${TMP_DIR}

eval "$(${CONDA_EXE} shell.bash hook)" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log 2>&1
conda activate ${CONDA_ENV} >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log 2>&1

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "错误:conda环境激活失败" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
fi

# Create HISAT2 index if it doesn't exist
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "索引不存在，正在创建HISAT2索引..." >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
    hisat2-build -p 4 ${GENOME_FASTA} ${INDEX_PREFIX} >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log 2>&1
    echo "索引创建完成" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
else
    echo "索引已存在，跳过创建步骤" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
fi

# 从cleandata目录获取样本文件
FASTQ_FILES_R1=($(ls ${CLEAN_DIR}/*_clean_R1.fq.gz 2>> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log | sort -V))
FASTQ_FILES_R2=($(ls ${CLEAN_DIR}/*_clean_R2.fq.gz 2>> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log | sort -V))

# Check file counts
if [ ${#FASTQ_FILES_R1[@]} -eq 0 ]; then
    echo "错误:未找到R1端fastq文件" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
elif [ ${#FASTQ_FILES_R2[@]} -eq 0 ]; then
    echo "错误:未找到R2端fastq文件" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
elif [ ${#FASTQ_FILES_R1[@]} -ne ${#FASTQ_FILES_R2[@]} ]; then
    echo "错误:R1/R2文件数量不匹配" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
fi

echo "找到 ${#FASTQ_FILES_R1[@]} 对双端fastq文件" >> ${LOG_DIR}/RNAseq_hisat2_${SLURM_JOB_ID}.log

# Process samples in parallel using GNU Parallel
# Note: Since each hisat2 uses 4 threads (-c 4), we set -j 1 to avoid oversubscription
# If you have more cores available, you could adjust these numbers accordingly
parallel --link -j 1 "
    SAMPLE_NAME=\$(basename {1} _clean_R1.fq.gz);
    SAMPLE_LOG=\"${LOG_DIR}/RNAseq_hisat2_\${SLURM_JOB_ID}_\${SAMPLE_NAME}.log\";

    echo \"正在处理样本: \${SAMPLE_NAME}\" > \${SAMPLE_LOG};
    echo \"内存使用情况:\" >> \${SAMPLE_LOG};
    free -h >> \${SAMPLE_LOG};

    hisat2 -x ${INDEX_PREFIX} \
           -1 {1} \
           -2 {2} \
           -S ${SAM_DIR}/\${SAMPLE_NAME}.sam \
           --threads 4 >> \${SAMPLE_LOG} 2>&1;

    echo \"样本 \${SAMPLE_NAME} 处理完成\" >> \${SAMPLE_LOG}
" ::: ${FASTQ_FILES_R1[@]} ::: ${FASTQ_FILES_R2[@]}
