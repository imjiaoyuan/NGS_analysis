#!/bin/bash
#SBATCH --job-name=CUT_TAG_bowtie2
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/cut-tag

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/cut-tag/logs"
RAW_DIR="/mnt/may1nov1/u5023/pal/cut-tag/cleandata"
BAM_DIR="/mnt/may1nov1/u5023/pal/cut-tag/bam_files"
REPORT_DIR="/mnt/may1nov1/u5023/pal/cut-tag/report"
GENOME_FASTA="/mnt/may1nov1/u5023/pal/Palv.hap2.chr.genome.fa"
BOWTIE2_INDEX="/mnt/may1nov1/u5023/pal/cut-tag/bowtie2_index/Palv.hap2.chr.genome"

mkdir -p ${LOG_DIR} ${BAM_DIR} ${REPORT_DIR} ${BOWTIE2_INDEX%/*}

eval "$(${CONDA_EXE} shell.bash hook)"
conda activate ${CONDA_ENV}

if [[ -z "${CONDA_DEFAULT_ENV}" ]]; then
    echo "conda激活失败" >> ${LOG_DIR}/CUT_TAG_bowtie2_${SLURM_JOB_ID}.log
    exit 1
fi

if [ ! -f "${BOWTIE2_INDEX}.1.bt2" ]; then
    bowtie2-build --threads 8 ${GENOME_FASTA} ${BOWTIE2_INDEX} >> ${LOG_DIR}/CUT_TAG_bowtie2_${SLURM_JOB_ID}.log 2>&1
fi

FASTQ_FILES_R1=($(ls ${RAW_DIR}/*_1.fq.gz 2>> ${LOG_DIR}/CUT_TAG_bowtie2_${SLURM_JOB_ID}.log | sort -V))
FASTQ_FILES_R2=($(ls ${RAW_DIR}/*_2.fq.gz 2>> ${LOG_DIR}/CUT_TAG_bowtie2_${SLURM_JOB_ID}.log | sort -V))

for ((i=0; i<${#FASTQ_FILES_R1[@]}; i++)); do
    R1=${FASTQ_FILES_R1[$i]}
    R2=${FASTQ_FILES_R2[$i]}
    SAMPLE_NAME=$(basename ${R1} _1.fq.gz | sed 's/-/_/g')

    bowtie2 -x ${BOWTIE2_INDEX} \
            -1 ${R1} \
            -2 ${R2} \
            --very-sensitive \
            --no-mixed \
            --no-discordant \
            --phred33 \
            -I 10 \
            -X 700 \
            --threads 8 | \
    samtools sort -@ 8 -o ${BAM_DIR}/${SAMPLE_NAME}.bam -
    
    samtools flagstat ${BAM_DIR}/${SAMPLE_NAME}.bam > ${REPORT_DIR}/${SAMPLE_NAME}.flagstat.txt
    
    samtools view -h ${BAM_DIR}/${SAMPLE_NAME}.bam | \
    grep -v -e "MT" -e "chrM" | \
    samtools view -b -o ${BAM_DIR}/${SAMPLE_NAME}.filtered.bam
    
    samtools index ${BAM_DIR}/${SAMPLE_NAME}.filtered.bam
done
