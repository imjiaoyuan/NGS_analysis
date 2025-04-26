#!/bin/bash
#SBATCH --job-name=RNAseq_featureCounts
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -D /mnt/may1nov1/u5023/pal/rna-seq
#SBATCH -o /mnt/may1nov1/u5023/pal/rna-seq/logs/RNAseq_featureCounts_%j.out
#SBATCH -e /mnt/may1nov1/u5023/pal/rna-seq/logs/RNAseq_featureCounts_%j.err

CONDA_EXE="/mnt/may1nov1/u5023/miniforge3/bin/conda"
CONDA_ENV="omics"
LOG_DIR="/mnt/may1nov1/u5023/pal/rna-seq/logs"
BAM_DIR="/mnt/may1nov1/u5023/pal/rna-seq/bam_files"
GTF_FILE="/mnt/may1nov1/u5023/pal/Palv_new.gff"

mkdir -p ${LOG_DIR}

# Activate conda environment
eval "$(${CONDA_EXE} shell.bash hook)"
if ! conda activate ${CONDA_ENV}; then
    echo "Error: Failed to activate conda environment" >&2
    exit 1
fi
samtools sort -@4 -o /mnt/may1nov1/u5023/pal/rna-seq/bam_files/palroot1.sorted.bam  /mnt/may1nov1/u5023/pal/rna-seq/sam_files/palroot1.sam && samtools index -@4 /mnt/may1nov1/u5023/pal/rna-seq/bam_files/palroot1.sorted.bam
# Check GTF file exists
if [ ! -f "${GTF_FILE}" ]; then
    echo "Error: GTF file not found: ${GTF_FILE}" >&2
    exit 1
fi

# Get sorted BAM files
BAM_FILES=($(ls ${BAM_DIR}/*.sorted.bam 2>/dev/null | sort -V))

if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "Error: No BAM files found in ${BAM_DIR}" >&2
    exit 1
fi

echo "Found ${#BAM_FILES[@]} BAM files for quantification"

# Run featureCounts on all BAM files
echo "Running featureCounts quantification..."
featureCounts -T ${SLURM_CPUS_PER_TASK} \
    -t exon \
    -g gene_id \
    -a ${GTF_FILE} \
    -o ./all_samples_counts.txt \
    -p --countReadPairs \
    "${BAM_FILES[@]}"

if [ $? -eq 0 ]; then
    echo "Quantification completed successfully"
else
    echo "Error: featureCounts failed" >&2
    exit 1
fi
