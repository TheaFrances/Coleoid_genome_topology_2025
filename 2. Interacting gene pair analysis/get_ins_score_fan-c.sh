#!/bin/bash
#SBATCH --job-name=fan-c_ins_score
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 06:00:00


echo start
date

module load fanc

fanc insulation -o bigwig -w 350kb 500kb 700kb 1mb 2mb 3mb -i  409493_intrachrom.allValidPairs.hic@100kb
fanc insulation -o bed -w 350kb 500kb 700kb 1mb 2mb 3mb -i  409493_intrachrom.allValidPairs.hic@100kb

echo finish
date
