#!/bin/bash
#SBATCH --job-name=index_map_eupsc_prots_to_sepof
#SBATCH --cpus-per-task=24
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 1-00:00:00


miniprot -t24 --gff GCA_964300435.1_xcSepOffi3.1_genomic.fna eup_prot.fa > sanger_sepof_eup_prot.gff

