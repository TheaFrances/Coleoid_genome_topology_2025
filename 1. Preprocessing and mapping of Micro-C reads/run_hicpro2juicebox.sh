#!/bin/bash
#SBATCH --job-name=hicpro2juicebox
#SBATCH --cpus-per-task=1
#SBATCH --mem=200GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time=03-00:00:00


echo "Start"
date

# Run hicpro2juicebox to convert to .hic with KR normalisation
./hicpro2juicebox.sh \
  -i hic_results/data/409493/409493_intrachrom_rm_scaffolds.allValidPairs \
  -g annotation/Lachesis_assembly_chrsize.txt \
  -j juicer_tools_1.22.01.jar

echo "Finish"
date

