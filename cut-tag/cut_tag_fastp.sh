#!/bin/bash
#SBATCH --job-name=CUTTAG_fastp
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/cut-tag
#SBATCH -o /mnt/may1nov1/u5023/pal/cut-tag/logs/CUTTAG_fastp_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/cut-tag/logs/CUTTAG_fastp_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
RAW_DIR="/mnt/may1nov1/u5023/pal/cut-tag/rawdata"
CLEAN_DIR="/mnt/may1nov1/u5023/pal/cut-tag/cleandata"
REPORT_DIR="/mnt/may1nov1/u5023/pal/cut-tag/report"
TMP_DIR="/mnt/may1nov1/u5023/pal/cut-tag/tmp"

# 创建所需目录
mkdir -p ${CLEAN_DIR} ${REPORT_DIR} ${TMP_DIR}

# 激活conda环境
eval "$(${CONDA_EXE} shell.bash hook)"
if ! conda activate ${CONDA_ENV}; then
    echo "错误: conda环境激活失败" >&2
    exit 1
fi

# 获取样本列表
SAMPLES=()
for R1 in ${RAW_DIR}/*_1.fq.gz; do
    SAMPLE=$(basename ${R1} _1.fq.gz)
    R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"
    if [[ -f ${R2} ]]; then
        SAMPLES+=(${SAMPLE})
    else
        echo "警告: 样本 ${SAMPLE} 缺少R2文件" >&2
    fi
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "错误: 未找到有效的样本对" >&2
    exit 1
fi

echo "找到 ${#SAMPLES[@]} 个样本"

# 处理函数
process_sample() {
    local SAMPLE=$1
    local R1="${RAW_DIR}/${SAMPLE}_1.fq.gz"
    local R2="${RAW_DIR}/${SAMPLE}_2.fq.gz"

    echo "开始处理样本: ${SAMPLE}"
    echo "输入文件: ${R1} 和 ${R2}"

    # fastp处理
    fastp -i ${R1} -I ${R2} \
          -o ${CLEAN_DIR}/${SAMPLE}_1.fq.gz \
          -O ${CLEAN_DIR}/${SAMPLE}_2.fq.gz \
          --detect_adapter_for_pe \
          -w 4 \
          --compression 9 \
          -h ${TMP_DIR}/${SAMPLE}.html \
          -j ${TMP_DIR}/${SAMPLE}.json

    # FastQC分析
    fastqc -o ${TMP_DIR} -t 4 ${CLEAN_DIR}/${SAMPLE}_1.fq.gz

    echo "样本 ${SAMPLE} 处理完成"
}

# 并行处理
export -f process_sample
export CONDA_EXE CONDA_ENV RAW_DIR CLEAN_DIR REPORT_DIR TMP_DIR

echo "开始并行处理..."
parallel -j ${SLURM_CPUS_PER_TASK} --delay 1 "process_sample {}" ::: ${SAMPLES[@]}

# 汇总结果
echo "运行multiqc..."
multiqc ${TMP_DIR} -o ${REPORT_DIR} -n "CUTTAG_post_QC_report"

# 整理报告
mkdir -p ${REPORT_DIR}/fastp_reports
mv ${TMP_DIR}/*.html ${REPORT_DIR}/fastp_reports/ 2>/dev/null
mv ${TMP_DIR}/*.json ${REPORT_DIR}/fastp_reports/ 2>/dev/null

# 清理
find ${TMP_DIR} -maxdepth 1 -type f -delete

echo "分析完成! 结果保存在: ${REPORT_DIR}"
