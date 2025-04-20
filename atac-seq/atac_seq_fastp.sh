#!/bin/bash
#SBATCH --job-name=ATACseq_fastp
#SBATCH -c 4
#SBATCH --mem=32G  # 增加内存以防止OOM
#SBATCH -D /mnt/may1nov1/u5023/pal/atac-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_fastp_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_fastp_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/atac-seq/logs"
RAW_DIR="/mnt/may1nov1/u5023/pal/atac-seq/rawdata"
CLEAN_DIR="/mnt/may1nov1/u5023/pal/atac-seq/cleandata"
REPORT_DIR="/mnt/may1nov1/u5023/pal/atac-seq/report"
TMP_DIR="/mnt/may1nov1/u5023/pal/atac-seq/tmp"

mkdir -p ${LOG_DIR} ${CLEAN_DIR} ${REPORT_DIR} ${TMP_DIR}

# 获取所有样本文件
FASTQ_FILES_R1=($(ls ${RAW_DIR}/*_R1_001.fastq.gz | sort -V))
FASTQ_FILES_R2=($(ls ${RAW_DIR}/*_R2_001.fastq.gz | sort -V))

# 检查文件数量
if [ ${#FASTQ_FILES_R1[@]} -eq 0 ] || [ ${#FASTQ_FILES_R2[@]} -eq 0 ] || [ ${#FASTQ_FILES_R1[@]} -ne ${#FASTQ_FILES_R2[@]} ]; then
    echo "错误: 文件数量问题" >&2
    exit 1
fi

# 激活conda环境
eval "$(${CONDA_EXE} shell.bash hook)"
conda activate ${CONDA_ENV}

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "错误:conda环境激活失败"
    exit 1
fi

# 使用GNU Parallel并行处理
echo "开始并行处理 ${#FASTQ_FILES_R1[@]} 个样本..."

parallel --link -j 2 "
    SAMPLE_NAME=\$(basename {1} _R1_001.fastq.gz);
    SAMPLE_LOG=\"${LOG_DIR}/ATACseq_fastp_\${SLURM_JOB_ID}_${SAMPLE_NAME}.log\";

    echo \"正在处理样本: ${SAMPLE_NAME}\" > \${SAMPLE_LOG};
    echo \"内存使用情况:\" >> \${SAMPLE_LOG};
    free -h >> \${SAMPLE_LOG};

    fastp -i {1} -I {2} \
          -o ${CLEAN_DIR}/\${SAMPLE_NAME}_R1_001.fastq.gz \
          -O ${CLEAN_DIR}/\${SAMPLE_NAME}_R2_001.fastq.gz \
          --detect_adapter_for_pe \
          -h ${TMP_DIR}/\${SAMPLE_NAME}.html \
          -j ${TMP_DIR}/\${SAMPLE_NAME}.json >> \${SAMPLE_LOG} 2>&1;

    fastqc -o ${TMP_DIR} -t 2 ${CLEAN_DIR}/\${SAMPLE_NAME}_R1_001.fastq.gz ${CLEAN_DIR}/\${SAMPLE_NAME}_R2_001.fastq.gz >> \${SAMPLE_LOG} 2>&1
" ::: ${FASTQ_FILES_R1[@]} ::: ${FASTQ_FILES_R2[@]}

# 运行multiqc
echo "运行multiqc汇总结果..."
multiqc ${TMP_DIR} -o ${REPORT_DIR} -n "ATACseq_post_QC_report"

# 整理报告文件
mkdir -p ${REPORT_DIR}/fastp_reports
mv ${TMP_DIR}/*.html ${REPORT_DIR}/fastp_reports/
mv ${TMP_DIR}/*.json ${REPORT_DIR}/fastp_reports/

# 清理临时文件
find ${TMP_DIR} -maxdepth 1 -type f -name "*.html" -o -name "*.json" -o -name "*_fastqc*" -delete

echo "所有样本处理完成！报告保存在 ${REPORT_DIR}"
