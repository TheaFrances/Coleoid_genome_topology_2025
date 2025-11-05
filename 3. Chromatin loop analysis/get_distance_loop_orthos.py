#-*- coding: utf-8 -*- 
# Add column of average intergenic distance of species orthologs to loop file (and include loop size in the output too, to plot later).
# ==============================================================================
# Imports======================================================================
# ==============================================================================
import argparse
import sys
from collections import defaultdict
import os
import random
import collections
import pandas as pd
# ==============================================================================
# Command line options==========================================================
# ==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument("orthos", type=str,
                    help="File of one-to-one orthologs with species 1 in the first column and species 2 in the second column")

parser.add_argument("bed_spp2", type=str,
                    help="bed file of species you want to check loop chromosome status for")

parser.add_argument("loop_sizes", type=str, 
                    help="loop tsv file with sizes (.loopsize file)")

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
    print("IMPORTANT: Scaffolds must be removed from bed file before running this script.")
    print("NOTE: Ortholog file must have species 1 in the first column and species 2 in the second column. Species 2 should be the species with the bed file inputted and the one you want to check distances in.")
    
    species_name = args.bed_spp2.split("/")[-1].split(".")[0]
    print("Checking loop ortholog distances for species:", species_name)

    spp1_ortho_dict = {}  # Dictionary of species 1 and other species orthologs
    with open(args.orthos, "r") as orthos:
        for line in orthos:
            col = line.rstrip().split("\t")
            spp1_g = col[0]
            spp2_g = col[1]
            spp1_ortho_dict[spp1_g] = spp2_g

    spp2_gene_data = defaultdict(list)
    with open(args.bed_spp2, "r") as bed2:
        for line in bed2:
            col = line.rstrip().split("\t")
            chrom = col[0]
            start_pos = int(col[1])  # Gene start position
            end_pos = int(col[2])    # Gene end position
            gene = col[3]
            spp2_gene_data[gene] = (chrom, start_pos, end_pos)

    # Store loop sizes from the loopsize file
    loop_sizes = {}
    with open(args.loop_sizes, "r") as loopsize_file:
        for line in loopsize_file:
            if not line.startswith("BIN"):
                col = line.rstrip().split("\t")
                chromosome = col[0]
                bin1_start = int(col[1])
                bin2_start = int(col[4])
                fdr = float(col[6])
                loop_size = int(col[-1])  #Loop size from the last column
                loop_sizes[(chromosome, bin1_start, bin2_start, fdr)] = loop_size
               

    counter_dist = 0
    outname = args.loop_genes.split(".genes")[0] + ".genes_rm_dups_" + species_name + "_loop_size_dist.tsv"
    with open(outname, "w") as outfile:
        with open(args.loop_genes, "r") as loop_file:
            # Write the header
            outfile.write("chromosome\tbin1_start_(bp)\tgenes_bin1\tbin2_start_(bp)\tgenes_bin2\tfdr\t")
            outfile.write(species_name + "_total_number_of_different_chromosomes_included_in_the_loops_start_and_end\tloop_size\t" + species_name + "_average_intergenic_distance\n")
            for line in loop_file:
                line = line.rstrip()
                if not line.startswith("chromosome"):
                    col = line.split("\t")
                    
                    chromosome = col[0]
                    bin1_start = int(col[1])  #Start of the first bin
                    genes_bin1 = col[2]
                    bin2_start = int(col[3])  #Start of the second bin
                    genes_bin2 = col[4]
                    fdr = float(col[5])  #Convert FDR to float
                    no_diff_chroms = (col[10])  #Total number of different chromosomes included in the loops start and end

                    genes1 = [gene.strip() for gene in genes_bin1.split(",")]  #Genes in bin 1
                    genes2 = [gene.strip() for gene in genes_bin2.split(",")]  #Genes in bin 2

                    distances = []
                    #Process genes1 and genes2 for orthologous gene pairs
                    for gene1 in genes1:
                        if gene1 in spp1_ortho_dict:
                            gene1_ortho = spp1_ortho_dict[gene1]
                            if gene1_ortho in spp2_gene_data:
                                chrom1, start_pos1, end_pos1 = spp2_gene_data[gene1_ortho]
                                
                                for gene2 in genes2:
                                    if gene2 in spp1_ortho_dict:
                                        gene2_ortho = spp1_ortho_dict[gene2]
                                        if gene2_ortho in spp2_gene_data:
                                            chrom2, start_pos2, end_pos2 = spp2_gene_data[gene2_ortho]

                                            #Only calculate the intergenic distance if on the same chromosome
                                            if chrom1 == chrom2:
                                                if start_pos1 < start_pos2:
                                                    distance = start_pos2 - end_pos1
                                                else:
                                                    distance = start_pos1 - end_pos2
                                                distances.append(distance)

                    # Calculate the average intergenic distance
                    if distances: # If distances is not empty
                        avg_distance = sum(distances) / len(distances)
                        counter_dist += 1 
                    else: # If there are no orthologos or just one ortholog, or orthologs in bin one and bin 2 are on different chromosomes, write NA 
                        avg_distance = "NA" 

                    # Extract the loop size from the loop_sizes dictionary, matching chromosome, bin start, end and the FDR to be extra safe
                    loop_size = loop_sizes.get((chromosome, bin1_start, bin2_start, fdr), "NA") #NA is the default value if the key is not found in the dictionary

                    # Write to outfile
                    outfile.write(f"{chromosome}\t{bin1_start}\t{genes_bin1}\t{bin2_start}\t{genes_bin2}\t{fdr}\t{no_diff_chroms}\t{loop_size}\t{avg_distance}\n")
    
    print("Total number of loops with orthologous intergenic distances printed = ", counter_dist)
    print("Output written to: ", outname)

if __name__ == "__main__":
    
    main()
