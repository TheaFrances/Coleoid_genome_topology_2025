#!/bin/bash
#SBATCH --job-name=get_all_repeats_norm
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 7-00:00:00

echo start
date

module load R

Rscript get_all_norm_repeats_per_category_summed_eupsc.R -i overlaps_eup_100k_repeats_gene_pairs.txt -d 409493_100000_EUPvs212489_50000_OBI_genom_dist_interact_threshold_10eupsc_10octbi_with_sof_eupsc_intergenic_genom_dist.txt -o all_repeats_norm_eupsc100k

echo finish
date

