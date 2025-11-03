# -*- coding: utf-8 -*- 
# Check if genes in conserved loop files are on the same or different chromosomes other species e.g. O. bimaculoides from conserved loop file output from the script check_conserved_loops.py.
# If genes are on the same chromosome, calculate the minimum distance between the genes in the other species.
# ==============================================================================
# Imports=======================================================================
# ==============================================================================
import argparse
import sys
from collections import defaultdict
# ==============================================================================
# Command line options==========================================================
# ==============================================================================
parser = argparse.ArgumentParser()

parser.add_argument("ortholog_file", type=str, help="File of species 1 orthologs (first column) and species you want to check chromosome status of (second column)")

parser.add_argument("chrom_file", type=str, help="Bed file of chromosome information for species you want to check chromosome status of")

parser.add_argument("conserved_loops_file", type=str, help="File of conserved loops to check chromosome status of")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

# ==============================================================================
# Main code=====================================================================
# ==============================================================================

def calculate_min_distance(genes1, genes2, gene_start_end_dict):
    distances = []  # List to store all calculated distances
    for gene1 in genes1:  # Iterate through first set of genes
        if gene1 in gene_start_end_dict:  # Check if gene has position data
            start1, end1 = map(int, gene_start_end_dict[gene1])  # Convert positions to integers 
            for gene2 in genes2:  # Iterate through second set of genes
                if gene2 in gene_start_end_dict:  # Check if gene has position data
                    start2, end2 = map(int, gene_start_end_dict[gene2])  # Convert positions to integers
                    distance = min(abs(start1 - end2), abs(start2 - end1))  # Calculate minimum distance between genes
                    distances.append(distance)  # Store calculated distance
    if distances:  # If any distances were calculated
        min_distance = min(distances)  # Find the minimum distance
        return min_distance  # Return the minimum distance
    else:
        return "NA_needs_checking"



def main():

     print("IMPORTANT: Scaffolds must be removed from bed file before running this script.")

    chrom_check_species = args.chrom_file.split("/")[-1].split(".")[0]
    print("Checking chromosome status of conserved loops in:", chrom_check_species)

    # Dictionary to store ortholog data
    ortholog_dict = {}
    with open(args.ortholog_file, "r") as infile:
        for line in infile:
            line = line.strip().split("\t")
            gene_spp1 = line[0]
            gene_spp_check = line[1]
            ortholog_dict[gene_spp1] = gene_spp_check

    # Dictionary to store chromosome data
    chrom_dict = defaultdict(list)
    gene_start_end_dict = {}
    with open(args.chrom_file, "r") as infile:
        for line in infile:
            line = line.strip().split("\t")
            chrom = line[0]
            gene_start = line[1]
            gene_end = line[2]
            gene = line[3]
            chrom_dict[gene].append(chrom)
            gene_start_end_dict[gene] = [gene_start, gene_end]

    counter_loop = 0
    counter_same_chrom = 0
    counter_diff_chrom = 0

    # Process conserved loops file
    output = args.conserved_loops_file.split(".tsv")[0] + "_" + chrom_check_species + "_chrom_status_dist.tsv"
    with open(args.conserved_loops_file, "r") as consv_loop_file:
        with open(output, "w") as outfile:
            for line in consv_loop_file:
                line = line.strip()
                if line.endswith("bin2"): # Start of header varies, so use endswith
                    outfile.write(line + "\t" + chrom_check_species +"_chrom_status\t" + chrom_check_species + "_min_distance\n") #Write header
                else:
                    counter_loop += 1
                    species_to_check_orthos_bin1 = 0 
                    species_to_check_orthos_bin2 = 0 
                    line = line.split("\t")
                    spp1_genes_bin1 = line[1].split(", ")
                    spp1_genes_bin2 = line[2].split(", ")
                    chrom_count = set()
                    for gene in spp1_genes_bin1:
                        if gene in ortholog_dict:
                            species_to_check_orthos_bin1 += 1
                            gene_check = ortholog_dict[gene]
                            if gene_check in chrom_dict:
                                chrom_to_check = chrom_dict[gene_check]
                                chrom_count.update(chrom_to_check) # Add doesn't work here because it's a list
                    for gene in spp1_genes_bin2:
                        if gene in ortholog_dict:
                            species_to_check_orthos_bin2 += 1
                            gene_check = ortholog_dict[gene]
                            if gene_check in chrom_dict:
                                chrom_to_check = chrom_dict[gene_check]
                                chrom_count.update(chrom_to_check)
                    if species_to_check_orthos_bin1 >= 1 and species_to_check_orthos_bin2 >= 1: #If there are orthologs in both bins for the species you're checking chromosome status of
                        if len(chrom_count) == 1:
                            # Calculate minimum distance
                            ortho_genes1 = [ortholog_dict[ortho] for ortho in spp1_genes_bin1 if ortho in ortholog_dict and ortholog_dict[ortho] in gene_start_end_dict] #If gene is in ortholog_dict and ortholog_dict[gene] is in gene_start_end_dict, then add ortholog_dict[gene] to genes1
                            ortho_genes2 = [ortholog_dict[ortho] for ortho in spp1_genes_bin2 if ortho in ortholog_dict and ortholog_dict[ortho] in gene_start_end_dict]
                            if len(ortho_genes1) >=1 and len(ortho_genes2) >=1: # Only calculate distance if the orthologs have POSITION DATA!
                                min_distance = calculate_min_distance(ortho_genes1, ortho_genes2, gene_start_end_dict)
                                chrom_status = "same_chrom"
                                counter_same_chrom += 1
                                outfile.write("\t".join(line) + "\t" + chrom_status + "\t" + str(min_distance) + "\n")
                            else:
                                outfile.write("\t".join(line) + "\t" + chrom_status + "\tNA_no_position_data\n")    
                        elif len(chrom_count) > 1:
                            chrom_status = "diff_chrom"
                            counter_diff_chrom += 1
                            outfile.write("\t".join(line) + "\t" + chrom_status + "\tNA_diff_chrom\n")
                    elif (species_to_check_orthos_bin1 + species_to_check_orthos_bin2) == 1: #If there are orthologs in only one bin for the species you're checking chromosome status of
                        chrom_status = "NA"
                        outfile.write("\t".join(line) + "\t" + chrom_status + "\tNA_orthos_in_only_one_bin_" + chrom_check_species + "\n")
                    elif species_to_check_orthos_bin1 == 0 and species_to_check_orthos_bin2 == 0:
                        chrom_status = "NA"
                        outfile.write("\t".join(line) + "\t" + chrom_status + "\tNA_no_orthos_in_" + chrom_check_species + "\n")

    print("Number of conserved loops in input file = ", counter_loop)
    print("Number of conserved loops with genes on the same chromosome = ", counter_same_chrom)
    print("Number of conserved loops with genes on different chromosomes = ", counter_diff_chrom)
    print("Output written to: ", output)

if __name__ == "__main__":
    main()




                    



