#!/bin/bash
#SBATCH --job-name=ATACseq_hisat2
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/atac-seq

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/atac-seq/logs"
RAW_DIR="/mnt/may1nov1/u5023/pal/atac-seq/rawdata"
CLEAN_DIR="/mnt/may1nov1/u5023/pal/atac-seq/cleandata"
SAM_DIR="/mnt/may1nov1/u5023/pal/atac-seq/sam_files"
REPORT_DIR="/mnt/may1nov1/u5023/pal/atac-seq/report"
TMP_DIR="/mnt/may1nov1/u5023/pal/atac-seq/tmp"
GENOME_FASTA="/mnt/may1nov1/u5023/pal/Palv.hap2.chr.genome.fa"
INDEX_PREFIX="/mnt/may1nov1/u5023/pal/atac-seq/Palv.hap2.chr.genome"

mkdir -p ${LOG_DIR} ${CLEAN_DIR} ${SAM_DIR} ${REPORT_DIR} ${TMP_DIR}

eval "$(${CONDA_EXE} shell.bash hook)" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log 2>&1
conda activate ${CONDA_ENV} >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log 2>&1

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "错误:conda环境激活失败" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
fi

if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "索引不存在，正在创建HISAT2索引..." >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
    hisat2-build -p 4 ${GENOME_FASTA} ${INDEX_PREFIX} >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log 2>&1
    echo "索引创建完成" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
else
    echo "索引已存在，跳过创建步骤" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
fi

FASTQ_FILES_R1=($(ls ${RAW_DIR}/*_R1_001.fastq.gz 2>> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log | sort -V))
FASTQ_FILES_R2=($(ls ${RAW_DIR}/*_R2_001.fastq.gz 2>> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log | sort -V))

if [ ${#FASTQ_FILES_R1[@]} -eq 0 ]; then
    echo "错误:未找到R1端fastq文件" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
elif [ ${#FASTQ_FILES_R2[@]} -eq 0 ]; then
    echo "错误:未找到R2端fastq文件" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
elif [ ${#FASTQ_FILES_R1[@]} -ne ${#FASTQ_FILES_R2[@]} ]; then
    echo "错误:R1/R2文件数量不匹配" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
    exit 1
fi

echo "找到 ${#FASTQ_FILES_R1[@]} 对双端fastq文件" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log

for ((i=0; i<${#FASTQ_FILES_R1[@]}; i++)); do
    R1=${FASTQ_FILES_R1[$i]}
    R2=${FASTQ_FILES_R2[$i]}
    SAMPLE_NAME=$(basename ${R1} _R1_001.fastq.gz)
    
    echo "正在处理样本: ${SAMPLE_NAME}" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
    
    hisat2 -x ${INDEX_PREFIX} \
           -1 ${R1} \
           -2 ${R2} \
           -no-mixed \
           --no-discordant \
           -S ${SAM_DIR}/${SAMPLE_NAME}.sam \
           --threads 4 >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log 2>&1
done

echo "处理完成" >> ${LOG_DIR}/ATACseq_hisat2_${SLURM_JOB_ID}.log
