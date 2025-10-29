#!/bin/bash
#SBATCH --job-name=interproscan_octbi_ncbi
#SBATCH --cpus-per-task=12
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 4-00:00:00


echo start
date
module load java
module load python3
module load perl
my_interproscan/interproscan-5.62-94.0/interproscan.sh -i obi_ncbi.prot -cpu 12 --goterms -b interproscan_octbi
echo finish
date