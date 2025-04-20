#!/bin/bash
#SBATCH --job-name=ATACseq_bwa
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/atac-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bwa_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/atac-seq/logs/ATACseq_bwa_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/atac-seq/logs"
RAW_DIR="/mnt/may1nov1/u5023/pal/atac-seq/cleandata"
BAM_DIR="/mnt/may1nov1/u5023/pal/atac-seq/bam_files"
REPORT_DIR="/mnt/may1nov1/u5023/pal/atac-seq/report"
GENOME_FASTA="/mnt/may1nov1/u5023/pal/Palv.hap2.chr.genome.fa"
BWA_INDEX="/mnt/may1nov1/u5023/pal/atac-seq/bwa_index/Palv.hap2.chr.genome"

mkdir -p ${LOG_DIR} ${BAM_DIR} ${REPORT_DIR} ${BWA_INDEX%/*}

eval "$(${CONDA_EXE} shell.bash hook)"
conda activate ${CONDA_ENV}

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "conda激活失败" >> ${LOG_DIR}/ATACseq_bwa_${SLURM_JOB_ID}.log
    exit 1
fi

# Create BWA index if it doesn't exist
if [ ! -f "${BWA_INDEX}.bwt" ]; then
    bwa index -p ${BWA_INDEX} ${GENOME_FASTA} >> ${LOG_DIR}/ATACseq_bwa_${SLURM_JOB_ID}.log 2>&1
fi

# Get all sample files
FASTQ_FILES_R1=($(ls ${RAW_DIR}/*_R1_001.fastq.gz 2>> ${LOG_DIR}/ATACseq_bwa_${SLURM_JOB_ID}.log | sort -V))
FASTQ_FILES_R2=($(ls ${RAW_DIR}/*_R2_001.fastq.gz 2>> ${LOG_DIR}/ATACseq_bwa_${SLURM_JOB_ID}.log | sort -V))

# Check file counts
if [ ${#FASTQ_FILES_R1[@]} -eq 0 ] || [ ${#FASTQ_FILES_R2[@]} -eq 0 ] || [ ${#FASTQ_FILES_R1[@]} -ne ${#FASTQ_FILES_R2[@]} ]; then
    echo "错误: 文件数量问题" >&2
    exit 1
fi

# Process samples in parallel using GNU Parallel
echo "开始并行处理 ${#FASTQ_FILES_R1[@]} 个样本..."

parallel --link -j 1 "  # Reduced to 1 job per node since each job uses 8 threads
    SAMPLE_NAME=\$(basename {1} _R1_001.fastq.gz);
    SAMPLE_LOG=\"${LOG_DIR}/ATACseq_bwa_\${SLURM_JOB_ID}_${SAMPLE_NAME}.log\";

    echo \"正在处理样本: ${SAMPLE_NAME}\" > \${SAMPLE_LOG};
    echo \"内存使用情况:\" >> \${SAMPLE_LOG};
    free -h >> \${SAMPLE_LOG};

    # Alignment and sorting
    bwa mem -t 8 -M -k 16 ${BWA_INDEX} {1} {2} | \
    samtools sort -@ 8 -o ${BAM_DIR}/\${SAMPLE_NAME}.bam - >> \${SAMPLE_LOG} 2>&1;

    # Flagstat report
    samtools flagstat ${BAM_DIR}/\${SAMPLE_NAME}.bam > ${REPORT_DIR}/\${SAMPLE_NAME}.flagstat.txt 2>> \${SAMPLE_LOG};

    # Filter out unwanted chromosomes
    samtools view -h ${BAM_DIR}/\${SAMPLE_NAME}.bam | \
    grep -v -e \"MW376760.1\" -e \"MZ675536.1\" | \
    samtools view -b -o ${BAM_DIR}/\${SAMPLE_NAME}.filtered.bam 2>> \${SAMPLE_LOG};

    # Index filtered BAM
    samtools index ${BAM_DIR}/\${SAMPLE_NAME}.filtered.bam >> \${SAMPLE_LOG} 2>&1
" ::: ${FASTQ_FILES_R1[@]} ::: ${FASTQ_FILES_R2[@]}

echo "所有样本处理完成！结果保存在 ${BAM_DIR}"
