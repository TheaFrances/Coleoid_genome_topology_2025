#!/bin/bash
#SBATCH --job-name=bedtools_atac
#SBATCH --cpus-per-task=24
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 1-00:00:00

ml bedtools

bedtools coverage -sorted \
  -g Lachesis_assembly_chrsize_rm_scaf.sorted.txt \
  -a eupsc_loops_50k+100k_3columns_sorted.bed \
  -b 97309_st29_normalised_peaks_3col_sorted_rm_scaf.bed \
  > cov_eupsc_loop_anchor.tsv

  bedtools coverage -sorted \
  -g Lachesis_assembly_chrsize_rm_scaf.sorted.txt \
  -a eupsc_loops_50k+100k_non_loop_anchor_regions.bed \
  -b 97309_st29_normalised_peaks_3col_sorted_rm_scaf.bed \
  > cov_eupsc_loop_anchor.tsv

  ```bash
bedtools coverage -sorted \
  -g Lachesis_assembly_chrsize_rm_scaf.sorted.txt \
  -a eupsc_loops_50k+100k_interloop_space_sorted.bed \
  -b 97309_st29_normalised_peaks_3col_sorted_rm_scaf.bed \
  > cov_eupsc_29cat_anchor.tsv


bedtools coverage -sorted \
  -g Lachesis_assembly_chrsize_rm_scaf.sorted.txt \
  -a eupsc_loops_50k+100k_non_interloop_space.bed \
  -b 97309_st29_normalised_peaks_3col_sorted_rm_scaf.bed \
  > cov_eupsc_29cat_anchor.tsv
