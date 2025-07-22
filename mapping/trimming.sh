#!/bin/bash
# Trim paired-end reads to 50 bp using Trimmomatic (CROP:50)

#SBATCH --job-name=trim50
#SBATCH --output=trim50.out
#SBATCH --error=trim50.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G

echo "Start trimming"
date

module load trimmomatic

trimmomatic PE \
  R1_raw.fastq.gz R2_raw.fastq.gz \
  R1_trimmed.fastq.gz R1_unpaired.fastq.gz \
  R2_trimmed.fastq.gz R2_unpaired.fastq.gz \
  CROP:50

echo "Trimming finished"
date

