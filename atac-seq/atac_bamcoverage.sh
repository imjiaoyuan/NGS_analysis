#!/bin/bash
#SBATCH --job-name=ATACseq_bw_conversion
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/atac-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bw_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bw_%j.err

# 加载conda环境
CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/atac-seq/logs"
BAM_DIR="/mnt/may1nov1/u5023/pal/atac-seq/bam_files"
BW_DIR="/mnt/may1nov1/u5023/pal/atac-seq/bw_files"
EFFECTIVE_GENOME_SIZE=352383277

# 创建输出目录
mkdir -p ${LOG_DIR} ${BW_DIR}

# 激活conda环境
eval "$(${CONDA_EXE} shell.bash hook)"
conda activate ${CONDA_ENV}

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "conda激活失败" >> ${LOG_DIR}/ATACseq_bw_${SLURM_JOB_ID}.log
    exit 1
fi

# 获取所有BAM文件
BAM_FILES=($(ls ${BAM_DIR}/*.filtered.bam 2>> ${LOG_DIR}/ATACseq_bw_${SLURM_JOB_ID}.log | sort -V))

# 检查文件数量
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "错误: 未找到BAM文件" >&2
    exit 1
fi

# 处理每个样本
for BAM_FILE in "${BAM_FILES[@]}"; do
    SAMPLE_NAME=$(basename ${BAM_FILE} .filtered.bam)
    SAMPLE_LOG="${LOG_DIR}/ATACseq_bw_${SLURM_JOB_ID}_${SAMPLE_NAME}.log"
    
    echo "正在处理样本: ${SAMPLE_NAME}" > ${SAMPLE_LOG}
    echo "内存使用情况:" >> ${SAMPLE_LOG}
    free -h >> ${SAMPLE_LOG}

    # 生成10bp分辨率的bigWig文件
    echo "正在生成10bp分辨率的bigWig文件..." >> ${SAMPLE_LOG}
    bamCoverage --bam ${BAM_FILE} -o ${BW_DIR}/${SAMPLE_NAME}.last-10.bw \
        --binSize 10 --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
        --normalizeUsing RPKM --smoothLength 50 -p 8 >> ${SAMPLE_LOG} 2>&1
    
    # 生成50bp分辨率的bigWig文件
    echo "正在生成50bp分辨率的bigWig文件..." >> ${SAMPLE_LOG}
    bamCoverage --bam ${BAM_FILE} -o ${BW_DIR}/${SAMPLE_NAME}.last-50.bw \
        --binSize 50 --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
        --normalizeUsing RPKM --smoothLength 150 -p 8 >> ${SAMPLE_LOG} 2>&1
    
    echo "${SAMPLE_NAME} 处理完成" >> ${SAMPLE_LOG}
done

echo "所有样本处理完成！bigWig文件保存在 ${BW_DIR}"
