#!/bin/bash
#SBATCH --job-name=ATACseq_peak_calling
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/atac-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_peaks_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_peaks_%j.err

# 加载conda环境
source "/mnt/may1nov1/u5023/miniforge3/etc/profile.d/conda.sh"
conda activate macs2

# 设置目录和参数
BAM_DIR="/mnt/may1nov1/u5023/pal/atac-seq/bam_files"
PEAKS_DIR="/mnt/may1nov1/u5023/pal/atac-seq/peaks"
GENOME_SIZE=438609455
THREADS=8

# 创建输出目录
mkdir -p ${PEAKS_DIR}

# 获取所有BAM文件
BAM_FILES=($(ls ${BAM_DIR}/*.filtered.bam | sort -V))

# 并行处理函数
process_sample() {
    local BAM_FILE=$1
    local SAMPLE_NAME=$(basename ${BAM_FILE} .filtered.bam)
    
    echo "正在处理样本: ${SAMPLE_NAME}"
    
    # MACS2峰值检测
    macs2 callpeak -g ${GENOME_SIZE} -t ${BAM_FILE} -f BAMPE \
        --nomodel --extsize 200 --shift -100 \
        --pvalue 0.01 -B -n ${SAMPLE_NAME} \
        --outdir ${PEAKS_DIR}
}

# 导出函数和变量以便并行使用
export -f process_sample
export PEAKS_DIR GENOME_SIZE

# 使用GNU parallel并行运行
printf "%s\n" "${BAM_FILES[@]}" | parallel -j ${THREADS} process_sample

echo "所有样本处理完成！峰值文件保存在 ${PEAKS_DIR}"
