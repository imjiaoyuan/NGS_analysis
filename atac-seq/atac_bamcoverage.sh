#!/bin/bash
#SBATCH --job-name=ATACseq_bw_conversion
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/atac-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bw_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bw_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="deeptools"
LOG_FILE="/mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bw_${SLURM_JOB_ID}.log"
BAM_DIR="/mnt/may1nov1/u5023/pal/atac-seq/bam_files"
BW_DIR="/mnt/may1nov1/u5023/pal/atac-seq/bw_files"
EFFECTIVE_GENOME_SIZE=438609455

mkdir -p ${BW_DIR} $(dirname ${LOG_FILE})

eval "$(${CONDA_EXE} shell.bash hook)"
conda activate ${CONDA_ENV}

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "conda激活失败" >> ${LOG_FILE}
    exit 1
fi

BAM_FILES=($(ls ${BAM_DIR}/*.filtered.bam 2>> ${LOG_FILE} | sort -V))

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "错误: 未找到BAM文件" >> ${LOG_FILE}
    exit 1
fi

for BAM_FILE in "${BAM_FILES[@]}"; do
    SAMPLE_NAME=$(basename ${BAM_FILE} .filtered.bam)
    
    echo "正在处理样本: ${SAMPLE_NAME}" >> ${LOG_FILE}
    echo "正在生成10bp分辨率的bigWig文件..." >> ${LOG_FILE}
    bamCoverage --bam ${BAM_FILE} -o ${BW_DIR}/${SAMPLE_NAME}.last-10.bw \
        --binSize 10 --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
        --normalizeUsing RPKM --smoothLength 50 -p 8 >> ${LOG_FILE} 2>&1

    echo "正在生成50bp分辨率的bigWig文件..." >> ${LOG_FILE}
    bamCoverage --bam ${BAM_FILE} -o ${BW_DIR}/${SAMPLE_NAME}.last-50.bw \
        --binSize 50 --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE} \
        --normalizeUsing RPKM --smoothLength 150 -p 8 >> ${LOG_FILE} 2>&1

    echo "${SAMPLE_NAME} 处理完成" >> ${LOG_FILE}
done

echo "所有样本处理完成！bigWig文件保存在 ${BW_DIR}" >> ${LOG_FILE}
