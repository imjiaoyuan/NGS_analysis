#!/bin/bash
#SBATCH --job-name=RNAseq_samtools
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/rna-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/rna-seq/logs/RNAseq_samtools_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/rna-seq/logs/RNAseq_samtools_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/rna-seq/logs"
BAM_DIR="/mnt/may1nov1/u5023/pal/rna-seq/bam_files"
TMP_DIR="/mnt/may1nov1/u5023/pal/rna-seq/tmp"
GENOME_FASTA="/mnt/may1nov1/u5023/pal/Palv.hap2.chr.genome.fa"
GTF_FILE="/mnt/may1nov1/u5023/pal/Palv_new.gff"

mkdir -p ${LOG_DIR} ${BAM_DIR} ${TMP_DIR}

# 激活conda环境
eval "$(${CONDA_EXE} shell.bash hook)"
if ! conda activate ${CONDA_ENV}; then
    echo "错误:conda环境激活失败" >&2
    exit 1
fi

# 检查GTF文件是否存在且有效
if [ ! -f "${GTF_FILE}" ]; then
    echo "错误:GTF文件不存在: ${GTF_FILE}" >&2
    exit 1
fi

SAM_FILES=($(ls /mnt/may1nov1/u5023/pal/rna-seq/sam_files/*.sam 2>/dev/null | sort -V))

if [ ${#SAM_FILES[@]} -eq 0 ]; then
    echo "错误:未找到SAM文件" >&2
    exit 1
fi

echo "找到 ${#SAM_FILES[@]} 个SAM文件"

# 定义处理函数
process_sample() {
    local SAM_FILE=$1
    local SAMPLE_NAME=$(basename ${SAM_FILE} .sam)
    local BAM_FILE="${BAM_DIR}/${SAMPLE_NAME}.sorted.bam"

    echo "正在处理样本: ${SAMPLE_NAME}"
    
    # 转换SAM为BAM并排序
    samtools sort -@ 4 -o ${BAM_FILE} ${SAM_FILE}
    
    # 创建索引
    samtools index -@ 4 ${BAM_FILE}
    
    echo "样本 ${SAMPLE_NAME} 处理完成"
    echo ${BAM_FILE}  # 返回BAM文件路径用于后续featureCounts
}

# 导出变量和函数
export -f process_sample
export CONDA_EXE CONDA_ENV LOG_DIR BAM_DIR TMP_DIR

# 并行处理SAM文件
echo "开始并行处理SAM文件..."
BAM_FILES=($(parallel -j ${SLURM_CPUS_PER_TASK} --delay 1 "process_sample {}" ::: ${SAM_FILES[@]}))

# 对所有BAM文件一起进行定量，使用gene_id作为基因标识符
echo "正在对所有BAM文件进行定量..."
featureCounts -T ${SLURM_CPUS_PER_TASK} \
    -t exon \
    -g gene_id \
    -a ${GTF_FILE} \
    -o ${BAM_DIR}/all_samples_counts.txt \
    -p --countReadPairs \
    "${BAM_FILES[@]}"

echo "所有样本处理完成"
