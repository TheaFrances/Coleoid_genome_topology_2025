# -*- coding: utf-8 -*- 
# Check if loops in one species are located on same or across different chromosomes in another species.
# ==============================================================================
# Imports=======================================================================
# ==============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
import re
# ==============================================================================
# Command line options==========================================================
# ==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument("orthos", type=str,
                    help="File of one-to-one orthologs with species 1 in the first column and species 2 in the second column")

parser.add_argument("bed_spp2", type=str,
                    help="bed file of species you want to check loop chromosome status for")

parser.add_argument("loop_genes", type=str,
                    help="loop tsv with genes (.genes file)")
 
if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()
# ==============================================================================
# Main code=====================================================================
# ==============================================================================

def main():

    print("IMPORTANT: Scaffolds must be removed from bed file before running this script. Also make sure you have removed duplicate gene interactions")
    print("NOTE: Ortholog file must have species 1 in the first column and species 2 in the second column. Species 2 should be the species with the bed file inputted and the one you want to check chromosome status in.")
    species_name1 = args.loop_genes.split("/")[-1].split("_")[0]
    species_name2 = args.bed_spp2.split("/")[-1].split(".")[0]
    print("Checking ", species_name2, "chromosome status for loops in", species_name1)

    spp1_ortho_dict = {} #Dictionary of species 1 and other species orthologs
    with open(args.orthos, "r") as orthos:
        for line in orthos:
            col = line.rstrip().split("\t")
            spp1_g = col[0]
            spp2_g = col[1] 
            spp1_ortho_dict[spp1_g] = spp2_g

    spp2_gene_chroms = defaultdict(list)
    with open(args.bed_spp2, "r") as bed2:
        for line in bed2:
            col = line.rstrip().split("\t")
            chrom = col[0]
            gene = col[3]
            spp2_gene_chroms[gene] = chrom

    counter_conserved = 0
    counter_rearranged = 0
    outname = args.loop_genes.split(".tsv.genes_rm_dups")[0] + ".genes_rm_dups." + species_name2 + "_chrom_status.tsv"
    with open(outname, "w") as outfile:
        with open(args.loop_genes, "r") as loop_file:
            outfile.write("chromosome\tbin1_start (bp)\tgenes_bin1\tbin2_start (bp)\tgenes_bin2\tfdr\t")
            outfile.write(species_name2 + " chromosomes included in the loop start\t" + species_name2 + " number of different chromosomes in the loop start\t" + species_name2 + " chromosomes included in the loop end\t" + species_name2 + " number of different chromosomes included in the loop end\t" + species_name2 + " total number of different chromosomes included in the loops start and end\n")
            
            for line in loop_file:
                line = line.rstrip()
                if not line.startswith("chr"):
                    outfile.write(line + "\t")
                    col = line.rstrip().split("\t")
                    
                    genes1 = [gene.strip() for gene in col[2].split(",")]  # This makes genes into a list as well as stripping whitespace
                    genes2 = [gene.strip() for gene in col[4].split(",")]

                    # Process genes1 for chrom1
                    chrom1_results = []
                    chrom1_set = set()  # To keep track of unique chromosomes (excluding 'no_ortho')
                    for gene1 in genes1:
                        if not gene1 in spp1_ortho_dict:
                            chrom1_results.append("no_ortho")
                        else:
                            gene1 = spp1_ortho_dict[gene1]
                            if not gene1 in spp2_gene_chroms:
                                chrom1_results.append("no_ortho_on_chrom")
                            if gene1 in spp2_gene_chroms:
                                chrom1 = spp2_gene_chroms[gene1]
                                chrom1_results.append(chrom1)
                                chrom1_set.add(chrom1)
                    outfile.write(", ".join(chrom1_results) + "\t")  # Write the list of chromosomes (and no_ortho) joined by commas
                    if len(chrom1_set) == 0:
                        outfile.write("NA\t")
                    else:
                        outfile.write(str(len(chrom1_set)) + "\t")  # Write the number of unique chromosomes, excluding 'no_ortho'

                    # Process genes2 for chrom2
                    chrom2_results = []
                    chrom2_set = set()  # To keep track of unique chromosomes (excluding 'no_ortho')
                    for gene2 in genes2:
                        if not gene2 in spp1_ortho_dict:
                            chrom2_results.append("no_ortho")
                        else:
                            gene2 = spp1_ortho_dict[gene2]
                            if not gene2 in spp2_gene_chroms:
                                chrom2_results.append("no_ortho_on_chrom")
                            if gene2 in spp2_gene_chroms:
                                chrom2 = spp2_gene_chroms[gene2]
                                chrom2_results.append(chrom2)
                                chrom2_set.add(chrom2)
                    outfile.write(", ".join(chrom2_results) + "\t")  # Write the list of chromosomes (and no_ortho) joined by commas
                    if len(chrom2_set) == 0:
                        outfile.write("NA\t")
                    else:
                        outfile.write(str(len(chrom2_set)) + "\t")  # Write the number of unique chromosomes, excluding 'no_ortho'
                    filtered_chrom1_results = [chrom for chrom in chrom1_results if chrom != "no_ortho" and chrom != "no_ortho_on_chrom"] #Get chrom1_results without 'no_ortho' and 'no_ortho_on_chrom'
                    filtered_chrom2_results = [chrom for chrom in chrom2_results if chrom != "no_ortho" and chrom != "no_ortho_on_chrom"] #Get chrom2_results without 'no_ortho' and 'no_ortho_on_chrom'
                    # Calculate the total unique chromosomes between loop start and end
                    total_unique_chromosomes = chrom1_set | chrom2_set  #Add sets together to get total unique chromosomes. You don't need to convert to set again because they're already 2 sets and will stay that way.
                    if len(total_unique_chromosomes) == 0:
                        outfile.write("NA_no_genes\n")
                    elif len(filtered_chrom1_results) == 1 and len(filtered_chrom2_results) == 0 or len(filtered_chrom1_results) == 0 and len(filtered_chrom2_results) == 1:
                        outfile.write("NA_only_one_gene\n") # We can't determine if the loop is conserved or rearranged if there is only one orthologous gene in the loop, so write NA_not_enough_genes
                    elif len(filtered_chrom1_results) == 1 and len(filtered_chrom2_results) == 1 and gene1 == gene2: #Just in case the there is only one gene in the loop and it is the same gene in the start and end
                        outfile.write("NA_only_one_gene\n") # We can't determine if the loop is conserved or rearranged if there is only one orthologous gene in the loop, so write NA_not_enough_genes
                    elif len(total_unique_chromosomes) == 1:
                        counter_conserved += 1
                        outfile.write(str(len(total_unique_chromosomes)) + "\n")
                    else:
                        outfile.write(str(len(total_unique_chromosomes)) + "\n")
                        counter_rearranged += 1

    print("Number of loops with genes on the same chromosome in " + species_name2 + " = ", counter_conserved)
    print("Number of loops with genes on different chromosomes in " + species_name2 + " = ", counter_rearranged)
    print("Percentage of loops with genes on different chromosomes in " + species_name2 + " = ", (counter_rearranged / (counter_conserved + counter_rearranged)) * 100, "%")
    print("Output file:", outname)
                    

if __name__ == "__main__":
    main()
