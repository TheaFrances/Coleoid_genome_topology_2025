#!/bin/bash
#SBATCH --job-name=intersect_interactions_repeats
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 7-00:00:00

echo start
date

module load bedtools

bedtools intersect -a 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist_sorted.bed -b repeats_eupsc_sorted.bed -sorted -wa -wb | python3 bedtools_intersect_to_repeat_score_pipe.py > overlaps_eup_100k_repeats_gene_pairs.txt 

echo finish
date

