#!/bin/bash
#SBATCH --job-name=hicpro_mapping
#SBATCH --cpus-per-task=14
#SBATCH --mem=200GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time=4-00:00:00


echo start
date

module load bowtie2
module load python3/3.10.8
module load conda
conda activate hicpro-3.1.0

HiC-Pro -i 409493 \
        -o 409493_hicpro_output \
        -c config-hicpro_eupsc.txt

echo finish
date

