#!/bin/bash
#SBATCH --job-name=eupsc_homer
#SBATCH --cpus-per-task=28
#SBATCH --mem=35GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 2-00:00:00



echo start
date
module load homer
findMotifsGenome.pl eupsc_loops_50k+100k_500bp_windows.bed Lachesis_assembly.fasta /path/to/out_dir -size given
echo finish
date