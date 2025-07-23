#!/bin/bash
#SBATCH --job-name=eupsc_repeats_out2gff
#SBATCH --cpus-per-task=1
#SBATCH --mem=300GB
#SBATCH --mail-type=ALL
#SBATCH --output=log/%x-%j.out
#SBATCH --error=log/%x-%j.err
#SBATCH --partition=himem
#SBATCH --time 7-00:00:00

echo start
date

#Script is in the util folder of Repeatmasker install on lisc, must be ran from there
perl /lisc/app/repeatmasker/4.1.7-p1-3.12.9-5.40.1/util/rmOutToGFF3.pl Lachesis_assembly.fasta.out > eupsc_repeats.gff

echo finish
date

